import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
"""Estimated metric with Kent's axes parallel-transported from the mean (scratch)."""
import numpy as np
import realdata_experiment as rx

def transport(v, m, x):
    """Parallel transport of tangent vector v at m to x along the great circle m->x."""
    c = float(np.clip(m @ x, -1, 1)); axis = np.cross(m, x); s = np.linalg.norm(axis)
    if s < 1e-12: return v
    k = axis / s; ang = np.arctan2(s, c)
    return v*np.cos(ang) + np.cross(k, v)*np.sin(ang) + k*(k@v)*(1-np.cos(ang))

def post_transport(pole, mean, kappa, axes, n_grid):
    pts, phis = rx.great_circle_points(pole, n_grid)
    e1, e2, t1, t2 = axes
    n = pole/np.linalg.norm(pole); t = np.cross(n, pts); t /= np.linalg.norm(t, axis=1, keepdims=True)
    E1 = np.array([transport(e1, mean, x) for x in pts]); E2 = np.array([transport(e2, mean, x) for x in pts])
    scale = np.sqrt(np.sum(t*E1, 1)**2/max(t1,1e-12) + np.sum(t*E2, 1)**2/max(t2,1e-12))
    ld = kappa*(pts@mean) + np.log(np.maximum(scale, 1e-300)); ld -= ld.max()
    d = np.exp(ld); d /= d.sum()*(2*np.pi/len(pts)); return pts, phis, d

sites = rx.load_sites(); per = {}
for site, e in sorted(sites.items()):
    L, P = e["lines"], e["poles"]
    if len(P) == 0 or len(L) < 2: continue
    for h in range(len(L)):
        tr = np.delete(L, h, 0); mean, k = rx.fisher_mean(tr)
        if not np.isfinite(k): continue
        ax = rx.kent_axes(tr, mean)
        for pole in P:
            pa, _, da = rx.posterior_on_circle("ambient", pole, mean, k, ax, 1440)
            pk, _, dk = rx.posterior_on_circle("kent", pole, mean, k, ax, 1440)
            pt, _, dt = post_transport(pole, mean, k, ax, 1440)
            tg = L[h]; base = rx.nlpd_of_target(pa, da, tg)
            per.setdefault(site, []).append((rx.nlpd_of_target(pk, dk, tg)-base, rx.nlpd_of_target(pt, dt, tg)-base))
scat = {"sv56", "sv61"}
for lab, i in [("kent as coded   - ambient", 0), ("kent transported - ambient", 1)]:
    for name, S_ in [("scattered", sorted(scat)), ("concentrated", sorted(set(per)-scat))]:
        v = np.array([np.mean([p[i] for p in per[s]]) for s in S_])
        print(f"  {lab:28s} {name:12s} {v.mean():+.3f} [{v.min():+.3f},{v.max():+.3f}] {int((v<0).sum())}/{len(v)}")
