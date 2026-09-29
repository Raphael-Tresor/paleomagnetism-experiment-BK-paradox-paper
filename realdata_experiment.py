"""Palaeomagnetic great-circle metric comparison.

Real-data test of the metric-selection rule of the paper's Section 5.2: when
conditioning on a measure-zero set, take the metric that whitens the measurement
noise.

The null set here is physics, not a modelling convenience. A rock carrying two
overlapping magnetization components traces, under progressive demagnetization, a path
whose component vectors span a plane through the origin; the normalized resultant is
confined to the intersection of that plane with S^2, i.e. a great circle. The
palaeomagnetic literature records these as best-fit planes (method code DE-BFP) and
conditions on them routinely.

Five analysts condition the same inputs on the same great circle and differ only in the
metric defining the conditioning limit:

    ambient   round (geodesic) metric on S^2   -- isotropic; also the Hausdorff
                                                  construction of Bungert & Wacker
    whitened  Mahalanobis Sigma^{-1} from the  -- what Section 5.2 predicts is correct
              site's anisotropic scatter
    naive     flat in (dec, inc) coordinates   -- the silent default, drops cos(inc)

    kent      Mahalanobis from Kent's (1982)   -- EXPLORATORY, added 2026-09-18. Same
              FB5 moment estimator, per site      rule as `whitened`, but the anisotropy
                                                  is the field's own FB5 estimator
                                                  rather than a raw orientation tensor.
    pooled    Mahalanobis from the scatter     -- EXPLORATORY, added 2026-09-18. Same
              pooled across ALL sites             rule again, but the metric is a
                                                  property of the measurement process,
                                                  so it is estimated once from every
                                                  site's residuals instead of from the
                                                  3-6 specimens of one site.

`whitened`, `ambient` and `naive` are the primary analysts and alone decide the verdict.
`kent` and `pooled` were added after the primary analysis was run; they are alternative
estimators of the SAME metric-selection rule, addressing the fact that a 2x2 scatter
matrix from 3-6 points is badly determined. They are exploratory.

Data provenance and licence in data/SOURCE.md.

Usage
-----
    python realdata_experiment.py            # full run: results.json + figure
    python realdata_experiment.py --quick    # smoke run, coarse grid, fewer bootstraps
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

SEED = 20260904

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
# The paper includes figures by relative path from the repo root, matching the
# convention of the SO(3) experiment.
REPO_ROOT = HERE.parent.parent
OUT_PDF = REPO_ROOT / "paleomag_great_circle.pdf"
OUT_JSON = HERE / "results.json"

# Primary analysts: the only ones that decide the verdict. Order matters for reporting; do not reorder.
ANALYSTS_PREREGISTERED = ("whitened", "ambient", "naive")

# Exploratory analysts added 2026-09-18, after the primary result was reported. Both
# instantiate the SAME selection rule as `whitened` with a better-conditioned estimate
# of the anisotropy.
ANALYSTS_EXPLORATORY = ("kent", "pooled")

ANALYSTS = ANALYSTS_PREREGISTERED + ANALYSTS_EXPLORATORY


# --------------------------------------------------------------------------- #
# MagIC file reading
# --------------------------------------------------------------------------- #
def read_magic(path: Path) -> list[dict]:
    """Read a MagIC tab-delimited table.

    Line 1 is `tab<TAB><table name>`, line 2 the column headers, data from line 3.
    """
    lines = path.read_text().splitlines()
    header = lines[1].split("\t")
    rows = []
    for line in lines[2:]:
        if not line.strip():
            continue
        values = line.split("\t")
        values += [""] * (len(header) - len(values))
        rows.append(dict(zip(header, values)))
    return rows


def direction_to_vector(dec, inc) -> np.ndarray:
    """(declination, inclination) in degrees -> unit vector in R^3."""
    d, i = np.radians(float(dec)), np.radians(float(inc))
    return np.array([np.cos(i) * np.cos(d), np.cos(i) * np.sin(d), np.sin(i)])


def vector_to_direction(v: np.ndarray) -> tuple[float, float]:
    """Unit vector -> (declination, inclination) in degrees, dec in [0, 360)."""
    v = v / np.linalg.norm(v)
    dec = np.degrees(np.arctan2(v[1], v[0])) % 360.0
    inc = np.degrees(np.arcsin(np.clip(v[2], -1.0, 1.0)))
    return dec, inc


def load_sites() -> dict[str, dict]:
    """Per-site specimen directions, split into lines and great-circle poles.

    Specimen names are prefixed by their site (`sv01a1` -> `sv01`), which is how the
    MagIC specimen table links to sites here.
    """
    specimens = read_magic(DATA / "specimens.txt")
    site_rows = read_magic(DATA / "sites.txt")

    sites: dict[str, dict] = {}
    for row in specimens:
        codes = row.get("method_codes", "")
        name = row.get("specimen", "")
        if not (row.get("dir_dec") and row.get("dir_inc")) or len(name) < 4:
            continue
        try:
            v = direction_to_vector(row["dir_dec"], row["dir_inc"])
        except ValueError:
            continue
        site = name[:4]
        entry = sites.setdefault(site, {"lines": [], "poles": [], "mad": []})
        if "DE-BFL" in codes:
            entry["lines"].append(v)
        elif "DE-BFP" in codes:
            entry["poles"].append(v)

    # Attach the published site mean and Fisher statistics where present.
    for row in site_rows:
        site = row.get("site", "")
        if site in sites and row.get("dir_dec") and row.get("dir_inc"):
            try:
                sites[site]["published_mean"] = direction_to_vector(
                    row["dir_dec"], row["dir_inc"]
                )
                sites[site]["alpha95"] = float(row["dir_alpha95"] or "nan")
                sites[site]["k"] = float(row["dir_k"] or "nan")
            except ValueError:
                pass

    for entry in sites.values():
        entry["lines"] = np.array(entry["lines"]) if entry["lines"] else np.empty((0, 3))
        entry["poles"] = np.array(entry["poles"]) if entry["poles"] else np.empty((0, 3))
    return sites


# --------------------------------------------------------------------------- #
# Directional statistics
# --------------------------------------------------------------------------- #
def fisher_mean(V: np.ndarray) -> tuple[np.ndarray, float]:
    """Fisher mean direction and concentration kappa for unit vectors V."""
    resultant = V.sum(axis=0)
    R = np.linalg.norm(resultant)
    n = len(V)
    mean = resultant / R
    kappa = (n - 1) / (n - R) if n > R else np.inf
    return mean, kappa


def orientation_tensor_axes(V: np.ndarray, mean: np.ndarray):
    """Principal axes and eigenvalues of the scatter of V about `mean`.

    Returns (e1, e2, tau1, tau2): the two orthonormal directions spanning the tangent
    plane at `mean`, ordered by decreasing scatter, with their eigenvalues. This is the
    empirical anisotropy the whitened metric uses.
    """
    # Project onto the tangent plane at `mean`.
    P = np.eye(3) - np.outer(mean, mean)
    residuals = V @ P.T
    T = (residuals.T @ residuals) / max(len(V), 1)
    w, U = np.linalg.eigh(T)
    order = np.argsort(w)[::-1]
    w, U = w[order], U[:, order]
    # The first two eigenvectors span the tangent plane; the third is ~`mean`.
    return U[:, 0], U[:, 1], float(w[0]), float(w[1])


def kent_axes(V: np.ndarray, mean: np.ndarray):
    """Anisotropy axes from Kent's (1982) FB5 moment estimator.

    EXPLORATORY (added 2026-09-18). Same role as `orientation_tensor_axes`: it returns
    (e1, e2, tau1, tau2) for the whitening metric, but obtains them from the moment
    estimator of the Fisher-Bingham (FB5) distribution rather than from a raw second
    moment. FB5 is the field's own anisotropic model (Kent 1982; Gallo et al. 2017) and
    is what PmagPy's `pmag.dokent()` fits, so this aligns the code with the literature
    the manuscript cites.

    Procedure (Kent 1982, section 4; see also Kasarapu 2015, eqs. 6-9):
      1. rotate the sample mean onto a pole, giving the dispersion matrix in that frame;
      2. diagonalise its lower 2x2 block to get the major/minor axis angle psi;
      3. form r1 = |xbar|, r2 = l1 - l2, and take the limiting moment estimates
             kappa ~ (2 - 2 r1 - r2)^-1 + (2 - 2 r1 + r2)^-1
             beta  ~ ((2 - 2 r1 - r2)^-1 - (2 - 2 r1 + r2)^-1) / 2.

    The FB5 concentrations along the two tangent axes are kappa - 2*beta and
    kappa + 2*beta. A Fisher-like law with concentration c has tangential variance ~1/c,
    so the whitening metric consumes tau_j = 1/c_j. Returning variances (not
    concentrations) keeps the contract of `orientation_tensor_axes` exactly, so
    `posterior_on_circle` is untouched.

    Falls back to the orientation tensor when the moment estimate is degenerate, which
    happens for near-coincident training directions.
    """
    n = len(V)
    resultant = V.sum(axis=0)
    norm = np.linalg.norm(resultant)
    if n < 3 or norm < 1e-12:
        return orientation_tensor_axes(V, mean)
    r1 = float(norm / n)

    # Orthonormal frame with `mean` as pole; the other two axes span the tangent plane.
    h1 = mean / np.linalg.norm(mean)
    seed = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(seed, h1)) > 0.9:
        seed = np.array([1.0, 0.0, 0.0])
    h2 = seed - np.dot(seed, h1) * h1
    h2 /= np.linalg.norm(h2)
    h3 = np.cross(h1, h2)

    # Dispersion matrix in that frame; its lower 2x2 block carries the anisotropy.
    Y = np.column_stack([V @ h1, V @ h2, V @ h3])
    S = (Y.T @ Y) / n
    B_lower = S[1:, 1:]
    w, U = np.linalg.eigh(B_lower)
    order = np.argsort(w)[::-1]
    w, U = w[order], U[:, order]
    r2 = float(w[0] - w[1])

    # Kent's limiting moment estimates for the concentration and ovalness.
    d_minus, d_plus = 2.0 - 2.0 * r1 - r2, 2.0 - 2.0 * r1 + r2
    if d_minus <= 1e-9 or d_plus <= 1e-9:
        return orientation_tensor_axes(V, mean)
    kappa = 1.0 / d_minus + 1.0 / d_plus
    beta = 0.5 * (1.0 / d_minus - 1.0 / d_plus)

    # Concentrations along the major/minor tangent axes; FB5 requires 2*beta < kappa.
    c_major, c_minor = kappa - 2.0 * beta, kappa + 2.0 * beta
    if not (np.isfinite(c_major) and np.isfinite(c_minor)) or min(c_major, c_minor) <= 1e-9:
        return orientation_tensor_axes(V, mean)

    # Back to R^3: the eigenvectors of the lower block live in the (h2, h3) plane.
    e_major = U[0, 0] * h2 + U[1, 0] * h3
    e_minor = U[0, 1] * h2 + U[1, 1] * h3
    # Largest variance along the LEAST concentrated axis, matching the (tau1 >= tau2)
    # ordering that `orientation_tensor_axes` guarantees.
    return e_major, e_minor, 1.0 / c_major, 1.0 / c_minor


def pooled_axes(sites: dict):
    """One anisotropy shared by every site, pooled over all training residuals.

    EXPLORATORY (added 2026-09-18). The whitening rule says the metric is the geometry
    of the MEASUREMENT. If the anisotropy comes from the instrument and the field
    protocol (compass and inclinometer precision, demagnetization fitting) it is a
    property of the process, common to every site, and estimating it per site from 3-6
    specimens fits five noisy parameters where one well-determined object would do.

    Each site's residuals are expressed in that site's own tangent frame, so the pooled
    tensor is computed in a COMMON frame rather than by averaging unrelated 3-d vectors.
    The result is one (e1, e2, tau1, tau2) reused at every site and every fold, which
    removes the per-fold estimation noise entirely.

    Caveat, and the reason this is a second analyst rather than a replacement: part of
    the within-site scatter is genuine within-flow variation (magnetic fabric, cooling
    history), which is per-flow physics and does NOT pool. This analyst is the right one
    exactly to the extent that instrumental error dominates.
    """
    coords = []
    for entry in sites.values():
        lines = entry["lines"]
        if len(lines) < 2:
            continue
        mean, _ = fisher_mean(lines)
        h1 = mean / np.linalg.norm(mean)
        seed = np.array([0.0, 0.0, 1.0])
        if abs(np.dot(seed, h1)) > 0.9:
            seed = np.array([1.0, 0.0, 0.0])
        h2 = seed - np.dot(seed, h1) * h1
        h2 /= np.linalg.norm(h2)
        h3 = np.cross(h1, h2)
        for v in lines:
            coords.append([float(v @ h2), float(v @ h3)])
    if len(coords) < 3:
        return None

    C = np.asarray(coords)
    T = (C.T @ C) / len(C)
    w, U = np.linalg.eigh(T)
    order = np.argsort(w)[::-1]
    w, U = w[order], U[:, order]
    # Returned in the per-site tangent frame; `run_loo` maps it into R^3 at each site.
    return (float(U[0, 0]), float(U[1, 0])), (float(U[0, 1]), float(U[1, 1])), float(
        w[0]
    ), float(w[1])


def pooled_axes_at_site(pooled, mean: np.ndarray):
    """Express the pooled anisotropy in R^3 at a site with the given mean direction."""
    (u11, u21), (u12, u22), tau1, tau2 = pooled
    h1 = mean / np.linalg.norm(mean)
    seed = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(seed, h1)) > 0.9:
        seed = np.array([1.0, 0.0, 0.0])
    h2 = seed - np.dot(seed, h1) * h1
    h2 /= np.linalg.norm(h2)
    h3 = np.cross(h1, h2)
    return u11 * h2 + u21 * h3, u12 * h2 + u22 * h3, tau1, tau2


def great_circle_points(pole: np.ndarray, n_grid: int) -> tuple[np.ndarray, np.ndarray]:
    """Uniform arclength parametrisation of the great circle with the given pole.

    Returns (points, phis): unit vectors on {v : <v, pole> = 0} and their angles.
    """
    pole = pole / np.linalg.norm(pole)
    # Any vector not parallel to the pole gives an orthonormal frame of the circle.
    seed = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(seed, pole)) > 0.9:
        seed = np.array([0.0, 0.0, 1.0])
    u = seed - np.dot(seed, pole) * pole
    u /= np.linalg.norm(u)
    w = np.cross(pole, u)
    phis = np.linspace(0.0, 2.0 * np.pi, n_grid, endpoint=False)
    pts = np.cos(phis)[:, None] * u + np.sin(phis)[:, None] * w
    return pts, phis


# --------------------------------------------------------------------------- #
# The three conditional posteriors on a great circle
# --------------------------------------------------------------------------- #
def posterior_on_circle(
    analyst: str,
    pole: np.ndarray,
    mean: np.ndarray,
    kappa: float,
    axes: tuple,
    n_grid: int,
):
    """Posterior density on the great circle {v : <v, pole> = 0} for one analyst.

    All three share the same likelihood -- a Fisher density about the site's estimated
    direction -- and differ only in the geometry used to take the conditioning limit,
    which enters as the reference measure on the circle:

      ambient   uniform arclength (the round metric's Hausdorff measure on the curve)
      whitened  reweighted by the anisotropic scatter, i.e. the Mahalanobis metric
                induced by the site's own error ellipse
      naive     uniform in the (dec, inc) chart, i.e. the arclength element of the flat
                chart metric, which differs from the round one by 1/cos(inc)
      kent      as `whitened`, with the axes from Kent's FB5 moment estimator
      pooled    as `whitened`, with one anisotropy pooled across all sites

    The three anisotropic analysts share one formula and differ ONLY in the (e1, e2,
    tau1, tau2) passed in `axes`; that is the whole point of the comparison, so keep it
    that way.

    Returns (points, phis, density) with density normalised over the circle.
    """
    pts, phis = great_circle_points(pole, n_grid)

    # Shared likelihood: Fisher(kappa) about the site mean.
    log_like = kappa * (pts @ mean)

    if analyst == "ambient":
        log_ref = np.zeros(len(pts))
    elif analyst in ("whitened", "kent", "pooled"):
        # The whitening metric stretches directions of small error and shrinks
        # directions of large error. Along the circle the induced arclength element
        # is sqrt(sum_j (<t, e_j>^2 / tau_j)) for the tangent t; the reference density
        # is proportional to that element.
        e1, e2, tau1, tau2 = axes
        tangents = np.gradient(pts, axis=0)
        tangents /= np.linalg.norm(tangents, axis=1, keepdims=True)
        c1 = tangents @ e1
        c2 = tangents @ e2
        scale = np.sqrt(c1**2 / max(tau1, 1e-12) + c2**2 / max(tau2, 1e-12))
        log_ref = np.log(np.maximum(scale, 1e-300))
    elif analyst == "naive":
        # Flat in (dec, inc): the chart's arclength element omits the cos(inc) factor
        # that the round metric carries, so the reference density picks up 1/cos(inc).
        inc = np.arcsin(np.clip(pts[:, 2], -1.0, 1.0))
        log_ref = -np.log(np.maximum(np.cos(inc), 1e-12))
    else:
        raise ValueError(analyst)

    log_dens = log_like + log_ref
    log_dens -= log_dens.max()
    dens = np.exp(log_dens)
    dens /= dens.sum() * (2 * np.pi / len(pts))  # normalise w.r.t. arclength
    return pts, phis, dens


def nlpd_of_target(pts, dens, target: np.ndarray) -> float:
    """Negative log predictive density of a held-out direction.

    The target does not lie exactly on the circle, so it is projected onto the circle
    and scored there. This is the same treatment for every analyst, so the comparison
    is fair; the projection discards the (shared) off-circle component.
    """
    idx = int(np.argmax(pts @ target))
    return float(-np.log(max(dens[idx], 1e-300)))


def angular_error(pts, dens, target: np.ndarray) -> float:
    """Geodesic angle (degrees) between the target and the posterior circular mean."""
    weights = dens / dens.sum()
    centroid = (weights[:, None] * pts).sum(axis=0)
    norm = np.linalg.norm(centroid)
    if norm < 1e-12:
        return float("nan")
    centroid /= norm
    return float(np.degrees(np.arccos(np.clip(abs(centroid @ target), -1.0, 1.0))))


def credible_arc_covers(pts, dens, target: np.ndarray, level: float) -> bool:
    """Does the highest-density credible set at `level` contain the target?"""
    weights = dens / dens.sum()
    order = np.argsort(weights)[::-1]
    cumulative = np.cumsum(weights[order])
    keep = order[: int(np.searchsorted(cumulative, level) + 1)]
    return int(np.argmax(pts @ target)) in set(keep.tolist())


# --------------------------------------------------------------------------- #
# Leave-one-out evaluation
# --------------------------------------------------------------------------- #
def run_loo(sites: dict, n_grid: int, level: float = 0.90, pooled=None):
    """Leave-one-line-out scoring over every site with >=1 plane and >=2 lines.

    For each held-out line, each great circle at that site yields one scored case.

    `pooled` is the shared anisotropy from `pooled_axes`, or None to skip that analyst.
    Note it is estimated ONCE from every site's lines, including targets held out later;
    that is a mild optimism for the `pooled` analyst alone. The per-site estimators see
    only `train`.
    """
    records = []
    for site, entry in sorted(sites.items()):
        lines, poles = entry["lines"], entry["poles"]
        if len(poles) == 0 or len(lines) < 2:
            continue
        for held in range(len(lines)):
            target = lines[held]
            train = np.delete(lines, held, axis=0)
            mean, kappa = fisher_mean(train)
            if not np.isfinite(kappa):
                continue
            # One set of axes per analyst; everything else about the posterior is shared.
            axes_for = {
                "whitened": orientation_tensor_axes(train, mean),
                "kent": kent_axes(train, mean),
            }
            axes_for["ambient"] = axes_for["naive"] = axes_for["whitened"]  # unused
            if pooled is not None:
                axes_for["pooled"] = pooled_axes_at_site(pooled, mean)
            for pole_index, pole in enumerate(poles):
                row = {"site": site, "held_out": held, "pole": pole_index}
                for analyst in ANALYSTS:
                    if analyst not in axes_for:
                        continue
                    pts, _, dens = posterior_on_circle(
                        analyst, pole, mean, kappa, axes_for[analyst], n_grid
                    )
                    row[f"nlpd_{analyst}"] = nlpd_of_target(pts, dens, target)
                    row[f"angerr_{analyst}"] = angular_error(pts, dens, target)
                    row[f"cover_{analyst}"] = credible_arc_covers(
                        pts, dens, target, level
                    )
                records.append(row)
    return records


def bootstrap_gap(records, a: str, b: str, n_boot: int, rng) -> dict:
    """Bootstrap CI for mean NLPD(a) - NLPD(b), resampling SITES (not cases).

    Resampling at site level respects the clustering of cases within a site.
    """
    by_site: dict[str, list] = {}
    for r in records:
        by_site.setdefault(r["site"], []).append(r)
    site_names = sorted(by_site)

    def gap(rows):
        return float(
            np.mean([r[f"nlpd_{a}"] for r in rows])
            - np.mean([r[f"nlpd_{b}"] for r in rows])
        )

    observed = gap(records)
    draws = []
    for _ in range(n_boot):
        picked = rng.choice(len(site_names), size=len(site_names), replace=True)
        rows = [r for i in picked for r in by_site[site_names[i]]]
        if rows:
            draws.append(gap(rows))
    draws = np.array(draws)
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return {
        "comparison": f"{a} - {b}",
        "observed_gap": observed,
        "ci95": [float(lo), float(hi)],
        "excludes_zero": bool(lo > 0 or hi < 0),
        "favours": a if observed < 0 else b,
        "n_boot": int(len(draws)),
    }


# --------------------------------------------------------------------------- #
# Figure
# --------------------------------------------------------------------------- #
def make_figure(sites, records, summary, out_pdf: Path, n_grid: int, pooled=None):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes_row = plt.subplots(1, 3, figsize=(13.5, 4.2))
    colours = {
        "whitened": "#1b6ca8",
        "ambient": "#e08214",
        "naive": "#7a7a7a",
        "kent": "#4a9d5f",
        "pooled": "#9467bd",
    }
    # Exploratory analysts are drawn dashed so the figure never presents them as part of
    # the primary comparison.
    styles = {a: ("--" if a in ANALYSTS_EXPLORATORY else "-") for a in ANALYSTS}
    shown = [a for a in ANALYSTS if a in summary["mean_nlpd"]]

    # (a) the posteriors on one representative great circle.
    ax = axes_row[0]
    site = "sv01" if "sv01" in sites else sorted(sites)[0]
    entry = sites[site]
    if len(entry["lines"]) >= 2:
        mean, kappa = fisher_mean(entry["lines"])
        anis = {
            "whitened": orientation_tensor_axes(entry["lines"], mean),
            "kent": kent_axes(entry["lines"], mean),
        }
    else:
        mean = entry.get("published_mean", entry["poles"][0])
        kappa = 20.0
        fallback = (np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), 1.0, 1.0)
        anis = {"whitened": fallback, "kent": fallback}
    anis["ambient"] = anis["naive"] = anis["whitened"]
    if pooled is not None:
        anis["pooled"] = pooled_axes_at_site(pooled, mean)
    pole = entry["poles"][0]
    for analyst in shown:
        if analyst not in anis:
            continue
        pts, phis, dens = posterior_on_circle(
            analyst, pole, mean, kappa, anis[analyst], n_grid
        )
        ax.plot(
            np.degrees(phis),
            dens,
            color=colours[analyst],
            label=analyst,
            lw=1.8,
            ls=styles[analyst],
        )
    ax.set_xlabel("position along great circle (deg)")
    ax.set_ylabel("posterior density")
    ax.set_title(f"(a) posteriors on one great circle, site {site}")
    ax.legend(frameon=False, fontsize=8)

    # Exploratory bars are hatched, matching the dashed curves in panel (a).
    hatches = ["//" if a in ANALYSTS_EXPLORATORY else "" for a in shown]

    # (b) mean NLPD per analyst, with bootstrap CI on the primary gap.
    ax = axes_row[1]
    means = [summary["mean_nlpd"][a] for a in shown]
    bars = ax.bar(range(len(shown)), means, color=[colours[a] for a in shown])
    for bar, hatch in zip(bars, hatches):
        bar.set_hatch(hatch)
    ax.set_xticks(range(len(shown)))
    ax.set_xticklabels(shown, fontsize=8)
    ax.set_ylabel("mean NLPD (lower is better)")
    ax.set_title("(b) held-out predictive score")
    for i, m in enumerate(means):
        ax.text(i, m, f"{m:.3f}", ha="center", va="bottom", fontsize=7)

    # (c) coverage of the nominal 90% credible arc.
    ax = axes_row[2]
    cover = [100 * summary["coverage"][a] for a in shown]
    bars = ax.bar(range(len(shown)), cover, color=[colours[a] for a in shown])
    for bar, hatch in zip(bars, hatches):
        bar.set_hatch(hatch)
    ax.axhline(90, ls="--", c="k", lw=1, label="nominal 90%")
    ax.set_xticks(range(len(shown)))
    ax.set_xticklabels(shown, fontsize=8)
    ax.set_ylabel("empirical coverage (%)")
    ax.set_title("(c) calibration of the 90% credible arc")
    ax.legend(frameon=False, fontsize=8)
    fig.text(
        0.5,
        -0.02,
        "hatched / dashed = exploratory (added after the primary analysis)",
        ha="center",
        fontsize=7,
        style="italic",
    )

    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true", help="coarse grid, few bootstraps")
    args = parser.parse_args()

    n_grid = 360 if args.quick else 2880
    n_boot = 200 if args.quick else 2000
    rng = np.random.default_rng(SEED)

    sites = load_sites()
    usable = {
        s: e for s, e in sites.items() if len(e["poles"]) > 0 and len(e["lines"]) >= 2
    }
    print(f"sites total {len(sites)}, usable (>=1 plane, >=2 lines) {len(usable)}")

    pooled = pooled_axes(usable)
    if pooled is not None:
        print(
            f"pooled anisotropy from all sites: tau1/tau2 = "
            f"{pooled[2] / max(pooled[3], 1e-12):.2f}"
        )

    records = run_loo(sites, n_grid=n_grid, pooled=pooled)
    print(f"leave-one-out cases: {len(records)}")

    summary = {
        "mean_nlpd": {
            a: float(np.mean([r[f"nlpd_{a}"] for r in records])) for a in ANALYSTS
        },
        "median_angular_error_deg": {
            a: float(np.median([r[f"angerr_{a}"] for r in records])) for a in ANALYSTS
        },
        "coverage": {
            a: float(np.mean([r[f"cover_{a}"] for r in records])) for a in ANALYSTS
        },
    }

    primary = bootstrap_gap(records, "whitened", "ambient", n_boot=n_boot, rng=rng)
    secondary = bootstrap_gap(records, "whitened", "naive", n_boot=n_boot, rng=rng)

    # Exploratory: do the better-conditioned estimators of the SAME rule beat ambient?
    exploratory = {
        f"{a}_vs_ambient": bootstrap_gap(records, a, "ambient", n_boot=n_boot, rng=rng)
        for a in ANALYSTS_EXPLORATORY
    }

    # The verdict is decided among the PRIMARY analysts only. Letting an exploratory
    # analyst win the primary comparison would convert a post-hoc addition into a
    # confirmatory result.
    prereg_nlpd = {a: summary["mean_nlpd"][a] for a in ANALYSTS_PREREGISTERED}
    best = min(prereg_nlpd, key=prereg_nlpd.get)
    h1_supported = bool(best == "whitened" and primary["excludes_zero"])
    best_overall = min(summary["mean_nlpd"], key=summary["mean_nlpd"].get)

    print("\n--- mean NLPD (lower is better) ---")
    for a in ANALYSTS_PREREGISTERED:
        print(f"  {a:9s} {summary['mean_nlpd'][a]:8.4f}")
    for a in ANALYSTS_EXPLORATORY:
        print(f"  {a:9s} {summary['mean_nlpd'][a]:8.4f}   [exploratory]")
    print("\n--- median angular error (deg) ---")
    for a in ANALYSTS:
        print(f"  {a:9s} {summary['median_angular_error_deg'][a]:8.2f}")
    print("\n--- coverage of nominal 90% arc ---")
    for a in ANALYSTS:
        print(f"  {a:9s} {100*summary['coverage'][a]:7.1f}%")
    print(
        f"\nprimary comparison whitened - ambient: gap {primary['observed_gap']:+.4f}"
        f"  95% CI [{primary['ci95'][0]:+.4f}, {primary['ci95'][1]:+.4f}]"
        f"  excludes zero: {primary['excludes_zero']}"
    )
    print(f"best primary analyst by NLPD: {best}")
    print(
        f"\nPRIMARY VERDICT: H1 {'SUPPORTED' if h1_supported else 'NOT supported'}"
    )

    print("\n--- EXPLORATORY: alternative estimators of the same rule vs ambient ---")
    print("    (no confirmatory weight)")
    for name, comp in exploratory.items():
        print(
            f"  {name:18s} gap {comp['observed_gap']:+.4f}"
            f"  95% CI [{comp['ci95'][0]:+.4f}, {comp['ci95'][1]:+.4f}]"
            f"  excludes zero: {comp['excludes_zero']}"
        )
    if best_overall != best:
        print(f"  (lowest NLPD overall is `{best_overall}`, an exploratory analyst)")

    results = {
        "seed": SEED,
        "n_grid": n_grid,
        "quick": bool(args.quick),
        "n_sites_total": len(sites),
        "n_sites_usable": len(usable),
        "n_loo_cases": len(records),
        "analysts_preregistered": list(ANALYSTS_PREREGISTERED),
        "analysts_exploratory": list(ANALYSTS_EXPLORATORY),
        "summary": summary,
        "primary_comparison": primary,
        "secondary_comparison": secondary,
        "exploratory_comparisons": exploratory,
        "pooled_anisotropy_ratio": (
            float(pooled[2] / max(pooled[3], 1e-12)) if pooled else None
        ),
        "best_analyst_by_nlpd": best,
        "best_analyst_by_nlpd_overall": best_overall,
        "h1_supported": h1_supported,
        "records": records,
    }
    OUT_JSON.write_text(json.dumps(results, indent=2, default=float))
    make_figure(sites, records, summary, OUT_PDF, n_grid=n_grid, pooled=pooled)
    print(f"\nwrote {OUT_JSON}")
    print(f"wrote {OUT_PDF}")


if __name__ == "__main__":
    main()
