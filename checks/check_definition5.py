"""Independent check of the posteriors on a great circle against Definition 5.

For a metric g' on S^2 given by a field of 3x3 matrices S_x (acting on tangent vectors),
Definition 5 conditions nu_F on A = {<n,x> = 0} as the limit a -> infinity of
    exp(-a d'(x,A)^2) nu_F(dx).
This script computes that pre-posterior by brute force on a band around A: d'(x,A) is
obtained by direct minimisation of sqrt((x-y)^T S_x (x-y)) over points y of A (the metric
is evaluated at the off-circle point x, small-band approximation of the geodesic distance),
and the mass is accumulated along A. The result is compared with the closed-form
posterior used in the experiment. It does NOT use the tangent-vector formula.

Metrics checked:
  coded        S_x = Sigma^+ / <m,x>^2         (what `kent` computes in realdata_experiment)
  transported  S_x = Kent's ellipse transported from m to x along great circles
  round        S_x = I                          (ambient, sanity check)
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import realdata_experiment as rx  # noqa: E402


def transport(v, m, x):
    """Parallel transport of the tangent vector v at m to x along the great circle m -> x."""
    c = float(np.clip(m @ x, -1, 1))
    axis = np.cross(m, x)
    s = np.linalg.norm(axis)
    if s < 1e-12:
        return v
    k = axis / s
    ang = np.arctan2(s, c)
    return v * np.cos(ang) + np.cross(k, v) * np.sin(ang) + k * (k @ v) * (1 - np.cos(ang))


def S_field(kind, axes, mean):
    e1, e2, t1, t2 = axes
    if kind == "round":
        return lambda x: np.eye(3) - np.outer(x, x)
    if kind == "coded":
        S = np.outer(e1, e1) / t1 + np.outer(e2, e2) / t2
        return lambda x: S / (mean @ x) ** 2
    if kind == "transported":
        def f(x):
            E1, E2 = transport(e1, mean, x), transport(e2, mean, x)
            return np.outer(E1, E1) / t1 + np.outer(E2, E2) / t2
        return f
    raise ValueError(kind)


def closed_form(kind, pts, pole, mean, kappa, axes):
    """Posterior density w.r.t. round arclength, as derived from Definition 5."""
    e1, e2, t1, t2 = axes
    n = pole / np.linalg.norm(pole)
    t = np.cross(n, pts)
    t /= np.linalg.norm(t, axis=1, keepdims=True)
    if kind == "round":
        f = np.ones(len(pts))
    elif kind == "coded":
        f = np.sqrt((t @ e1) ** 2 / t1 + (t @ e2) ** 2 / t2)
    else:
        E1 = np.array([transport(e1, mean, x) for x in pts])
        E2 = np.array([transport(e2, mean, x) for x in pts])
        f = np.sqrt(np.sum(t * E1, 1) ** 2 / t1 + np.sum(t * E2, 1) ** 2 / t2)
    ld = kappa * (pts @ mean) + np.log(f)
    d = np.exp(ld - ld.max())
    return d / (d.sum() * 2 * np.pi / len(pts))


def brute_force(kind, pole, mean, kappa, axes, n_phi=720, n_s=161, half_width=0.02, a=4e5):
    pts, _ = rx.great_circle_points(pole, n_phi)
    n = pole / np.linalg.norm(pole)
    S_of = S_field(kind, axes, mean)
    mass = np.zeros(n_phi)
    s_grid = np.linspace(-half_width, half_width, n_s)
    window = 6
    for i, c in enumerate(pts):
        for s in s_grid:
            x = np.cos(s) * c + np.sin(s) * n  # round-orthogonal offset from A
            S = S_of(x)
            idx = np.arange(i - window, i + window + 1) % n_phi
            diff = x[None, :] - pts[idx]
            d2 = np.einsum("ij,jk,ik->i", diff, S, diff)
            j = idx[int(np.argmin(d2))]
            w = np.exp(-a * d2.min() + kappa * (x @ mean)) * np.cos(s)  # round area element
            mass[j] += w
    return pts, mass / (mass.sum() * 2 * np.pi / n_phi)


def main():
    sites = rx.load_sites()
    print("site  metric        L1 distance (brute force vs closed form)")
    for site in ("sv61", "sv56"):
        e = sites[site]
        train = e["lines"][1:]  # one representative train set
        mean, kappa = rx.fisher_mean(train)
        axes = rx.kent_axes(train, mean)
        pole = e["poles"][0]
        for kind in ("round", "coded", "transported"):
            pts, bf = brute_force(kind, pole, mean, kappa, axes)
            cf = closed_form(kind, pts, pole, mean, kappa, axes)
            l1 = float(np.abs(bf - cf).sum() * 2 * np.pi / len(pts))
            print(f"{site}  {kind:12s}  {l1:.3e}")
        # how far apart the two estimated posteriors are, for scale
        cf_c = closed_form("coded", pts, pole, mean, kappa, axes)
        cf_t = closed_form("transported", pts, pole, mean, kappa, axes)
        print(f"{site}  coded vs transported closed forms: L1 {np.abs(cf_c-cf_t).sum()*2*np.pi/len(pts):.3e}")


if __name__ == "__main__":
    main()
