# CLAUDE.md — palaeomagnetic great-circle metric experiment

> **Maintenance rule (read first).** This file documents the code in this directory.
> **Whenever you change the code, update this file in the same edit** so it never drifts
> from reality. **Keep it under 300 lines.** If it grows past that, cut prose, not facts.

> **This directory never writes `.tex` files.** The manuscript is edited separately, via
> the `latex-writting-style` skill.

## Purpose

Real-data test of the metric-selection rule of the paper's §5.2 ("take the metric that
whitens the measurement noise"), on a null set that palaeomagnetists genuinely condition
on: a **remagnetization great circle**.

Answers **AE point 2** (real-data application) and **AE report comment 2** (the SO(3)
example's geometry is "known by construction"). Here the metric is read off the
measurement, not assumed.

## Files

| File | Role |
|---|---|
| `PREREGISTRATION.md` | Analysis frozen **before** any score was computed. Read before touching the scoring. |
| `EXPERIMENT_REPORT.md` | Findings, goal/recommendation reasoning, and results (§3d). |
| `DESIGN.md` | Dataset-independent design rationale. |
| `realdata_experiment.py` | The experiment. Run it to reproduce everything. |
| `check_anisotropy.py` | Gate 1: is the error anisotropic at all? |
| `check_null.py` | Gate 1 control: is that anisotropy a small-`n` artifact? |
| `check_geometry.py` | Gate 2: is `DE-BFP`'s (dec,inc) the plane's pole? |
| `diagnose.py` | **Exploratory** post-hoc diagnosis of the null result. Not pre-registered. |
| `sensitivity_outlier.py` | **Exploratory** leave-one-site-out influence analysis. Shows the null is genuine. |
| `naive_check.py` | **Exploratory** site-level test of all three analyst pairings. |
| `data/` | Frozen input + `SOURCE.md` (provenance, CC BY 4.0 licence, citations). |
| `results.json` | Generated. Every computed quantity. Do not hand-edit. |
| `gate_check_results.txt`, `diagnostics.txt` | Generated logs. |
| `../../paleomag_great_circle.pdf` | Generated figure, written to the **repo root**. |

The figure goes to the repo root because the paper includes figures by relative path
(same convention as the SO(3) experiment).

## How to run

```bash
python3 -m venv .venv && ./.venv/bin/pip install numpy scipy pandas matplotlib
./.venv/bin/python realdata_experiment.py           # full: ~10 s
./.venv/bin/python realdata_experiment.py --quick   # smoke: coarse grid, 200 bootstraps
```

`.venv/` is git-ignored. `SEED = 20260904`; reruns are bit-stable.

## Code map (`realdata_experiment.py`)

- `read_magic` — MagIC tab-delimited reader (line 1 `tab<TAB>table`, line 2 headers).
- `direction_to_vector` / `vector_to_direction` — (dec, inc) in degrees ↔ unit vector.
- `load_sites` — groups specimens by site prefix (`sv01a1` → `sv01`), splitting
  `DE-BFL` (lines) from `DE-BFP` (great-circle poles).
- `fisher_mean` — mean direction and concentration κ.
- `orientation_tensor_axes` — principal axes/eigenvalues of tangent-plane scatter; **this
  is the anisotropy the whitened metric consumes.**
- `kent_axes` — **exploratory (2026-09-18).** Same output contract, from Kent's (1982)
  FB5 moment estimator instead of a raw second moment. Returns *variances* `1/c_j` from
  the FB5 concentrations `kappa -/+ 2*beta`, so `posterior_on_circle` is unchanged.
  Falls back to the orientation tensor when the moment estimate is degenerate.
- `pooled_axes` / `pooled_axes_at_site` — **exploratory (2026-09-18).** One anisotropy
  estimated from every site's residuals in their own tangent frames, then mapped back
  into R^3 at each site. Removes per-fold estimation noise.
- `great_circle_points` — uniform-arclength parametrisation of `{v : ⟨v,pole⟩ = 0}`.
- `posterior_on_circle(analyst, …)` — **the heart of the experiment.** All analysts share
  one Fisher likelihood and differ *only* in the reference measure on the circle:
  `ambient` uniform arclength; `whitened`/`kent`/`pooled` reweighted by `Σ⁻¹` along the
  tangent (one formula, three different `axes`); `naive` picks up `1/cos(inc)` from
  treating the chart as flat.
- `nlpd_of_target`, `angular_error`, `credible_arc_covers` — the three scores.
- `run_loo` — leave-one-line-out over sites with ≥1 plane and ≥2 lines. Builds one
  `axes_for[analyst]` map per fold; per-site estimators see only `train`.
- `bootstrap_gap` — CI on the NLPD gap, **resampling sites, not cases** (clustering).
- `make_figure` — 3 panels: posteriors on one circle; mean NLPD; coverage. Exploratory
  analysts are drawn dashed/hatched.

## Invariants the result depends on

1. **`DE-BFP` (dec,inc) is the POLE to the plane**, not a direction on it. Verified:
   planes sit 87.9° from their site mean, lines 5.0° (`check_geometry.py`). If this were
   wrong the null set would be wrong and every number meaningless.
2. **The analysts differ ONLY in `log_ref`** inside `posterior_on_circle`.
   Likelihood, grid, and normalisation are shared. That single term *is* the experiment —
   keep it that way. The three anisotropic analysts additionally share one formula and
   differ only in the `axes` handed to them.
3. **The bootstrap resamples sites**, not leave-one-out cases. Cases within a site are
   dependent; resampling cases would understate the CI and manufacture significance.
4. **Scores are computed on the circle**, with the target projected onto it. The discarded
   off-circle component is identical across analysts, so the comparison stays fair.
5. **The verdict is decided among the PRE-REGISTERED analysts only** (`ANALYSTS_PREREGISTERED`).
   `kent` and `pooled` were added afterwards and must never decide `h1_supported`, or a
   post-hoc addition becomes a confirmatory result.

## Results (full run, 2026-09-04)

11 sites, 76 leave-one-out cases.

| analyst | mean NLPD | median ang. err | coverage (nominal 90%) |
|---|---|---|---|
| `whitened` | **6.5674** | 6.29° | 78.9% |
| `ambient` | 6.5853 | 6.30° | 80.3% |
| `naive` | 6.6092 | 6.64° | 77.6% |

`whitened - ambient` gap **−0.0179**, 95% CI **[−0.0711, +0.0049]** → **includes zero**.

**Pre-registered verdict: H1 NOT supported (H0 not rejected).** The ordering matches the
theory's prediction and is stable across grid resolutions, but the margin is not
distinguishable from zero.

### Exploratory estimators of the same rule (added 2026-09-18)

Both re-run above and **leave every pre-registered number bit-identical** — they add
analysts, they do not touch the primary path.

| analyst | mean NLPD | median ang. err | coverage | vs `ambient` | sites won |
|---|---|---|---|---|---|
| `kent` | 6.5692 | 6.28° | 78.9% | −0.0160, CI [−0.0661, +0.0056] | 6/11 |
| `pooled` | 6.5952 | 6.57° | 78.9% | **+0.0100**, CI [−0.0008, +0.0332] | 3/11 |

**Neither changes the verdict; both sharpen the diagnosis.**

- **`kent` ≈ `whitened`** (6.5692 vs 6.5674). Identical principal axes (0.0° apart at
  every site) but shrunken eigenvalue ratios where `n` is smallest: `sv56` 283 → 4.7,
  `sv61` 161 → 18. So the raw orientation tensor's extreme anisotropies **were** small-`n`
  artifacts, and correcting them changes nothing — the metric choice is immaterial here
  for reasons that survive a better estimator.
- **`pooled` is WORSE than `ambient`** (+0.0100, 3/11 sites), the only analyst to lose to
  it. Its across-site gap SD is **0.0606** vs `whitened`'s 0.1189, i.e. exactly half the
  variance — so pooling did remove the estimation noise, and the result still got worse.
  Reading: the within-site anisotropy is **not** a property of the measurement process
  common to all sites; pooling averages away genuine per-flow geometry (magnetic fabric,
  cooling history). Pooled `τ₁/τ₂ = 2.75` against per-site ratios of 1.2–283.

⚠️ `pooled_axes` is estimated once from **all** sites' lines, including directions later
held out. That is a mild optimism **for `pooled` alone** — and it still loses, so the
conclusion is conservative. Per-site estimators see only `train`.

**The null is genuine, not power-limited** (`sensitivity_outlier.py`): the gap is −0.0179
with all 11 sites but **−0.0017 without `sv61` alone**, whose 3 cases carry nearly the
whole difference. No single-site removal makes the CI exclude zero. At site level
`whitened` wins 6/11, sign p = 0.50.

⚠️ **Two figures withdrawn.** An earlier per-case sign test (50/76, p = 0.0040) ignored
within-site clustering — do not quote it. The pooled angular-error gap (6.64° vs 6.29°)
holds in only 4/11 sites — do not use it as evidence.

**The one supported comparison** (`naive_check.py`, site level): `ambient` beats `naive`
in **9/11 sites, sign p = 0.033**, i.e. the round metric outperforms treating the
(dec, inc) chart as flat. Its CI still includes zero ([−0.0923, +0.0008]), so it is a
consistent sign, not a confirmed effect. `whitened` vs `naive` is 6/11, p = 0.50.

Reading: **the two admissible metrics are indistinguishable here; the flat chart is
worse.** Both are what the theory predicts — these sites are strongly concentrated
(κ ≈ 80–1128), and the paper's SO(3) section already notes the metric matters only where
curvature does.

**Decision (2026-09-04):** report as-is and keep SO(3) as the primary demonstration.
A low-κ dataset search and a κ-stratified secondary analysis were both considered and
**declined**; see `EXPERIMENT_REPORT.md` §3d.

**Scored sites: 11, from 23 with great circles.** 23 → 14 (need ≥2 lines to hold one out)
→ 11 (finite Fisher κ). The three *pure* null-set sites `sv01`/`sv10`/`sv62` have zero
lines and are therefore **excluded** — a real limitation of leave-one-out here (§3g).

## If your change moves these numbers

Either it is a bug, or this table, `EXPERIMENT_REPORT.md` §3d, and any manuscript prose
must be updated together. Do **not** silently re-tune the scoring: the analysis is
pre-registered, and changes belong in `PREREGISTRATION.md` §9 as dated deviations.

## When you change the code — checklist

- [ ] Update the **Code map** / **Invariants** / **Results** above.
- [ ] Re-run the full (not `--quick`) experiment to refresh `results.json` + figure.
- [ ] Log any analysis change in `PREREGISTRATION.md` §9 with its reason.
- [ ] Keep this file **under 300 lines**.
