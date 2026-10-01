import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
"""Check of the Mahalanobis area factor and its effect on Table 1 (scratch, not committed)."""
import numpy as np
import realdata_experiment as rx
rng = np.random.default_rng(0)

# 1) closed form: sqrt(det of g'_x on T_x in a round-orthonormal basis) = |<m,x>| / sqrt(t1 t2)
m = rng.normal(size=3); m /= np.linalg.norm(m)
e1 = np.cross(m, rng.normal(size=3)); e1 /= np.linalg.norm(e1); e2 = np.cross(m, e1)
t1, t2 = 0.3, 0.05
S = np.outer(e1, e1)/t1 + np.outer(e2, e2)/t2
err = 0
for _ in range(1000):
    x = rng.normal(size=3); x /= np.linalg.norm(x)
    a = np.cross(x, rng.normal(size=3)); a /= np.linalg.norm(a); b = np.cross(x, a)
    B = np.stack([a, b], 1)
    lhs = np.sqrt(abs(np.linalg.det(B.T @ S @ B)))
    rhs = abs(m @ x) / np.sqrt(t1*t2)
    err = max(err, abs(lhs-rhs)/max(rhs, 1e-12))
print("max rel err of closed form sqrt(det G_x) = |<m,x>|/sqrt(t1 t2):", err)

# 2) effect on Table 1: kent with the missing area factor 1/|<m,x>|
def post_kent_area(pole, mean, kappa, axes, n_grid):
    pts, phis, _ = rx.posterior_on_circle("kent", pole, mean, kappa, axes, n_grid)
    e1, e2, tau1, tau2 = axes
    n = pole/np.linalg.norm(pole); t = np.cross(n, pts); t /= np.linalg.norm(t, axis=1, keepdims=True)
    scale = np.sqrt((t@e1)**2/max(tau1,1e-12) + (t@e2)**2/max(tau2,1e-12))
    ld = kappa*(pts@mean) + np.log(np.maximum(scale,1e-300)) - np.log(np.maximum(np.abs(pts@mean),1e-300))
    ld -= ld.max(); d = np.exp(ld); d /= d.sum()*(2*np.pi/len(pts))
    return pts, phis, d

sites = rx.load_sites(); per = {}
worst_mass = 0
for site, e in sorted(sites.items()):
    L, P = e["lines"], e["poles"]
    if len(P) == 0 or len(L) < 2: continue
    for h in range(len(L)):
        tr = np.delete(L, h, 0); mean, k = rx.fisher_mean(tr)
        if not np.isfinite(k): continue
        ax = rx.kent_axes(tr, mean)
        for pole in P:
            pa, _, da = rx.posterior_on_circle("ambient", pole, mean, k, ax, 2880)
            pk, _, dk = rx.posterior_on_circle("kent", pole, mean, k, ax, 2880)
            pk2, _, dk2 = post_kent_area(pole, mean, k, ax, 2880)
            far = np.abs(pk2@mean) < 0.2  # within ~78-102 deg of the mean
            worst_mass = max(worst_mass, float(dk2[far].sum()*2*np.pi/len(pk2)))
            tgt = L[h]
            per.setdefault(site, []).append((rx.nlpd_of_target(pk, dk, tgt) - rx.nlpd_of_target(pa, da, tgt),
                                             rx.nlpd_of_target(pk2, dk2, tgt) - rx.nlpd_of_target(pa, da, tgt)))
scat = {"sv56", "sv61"}
for lab, i in [("kent (code) - ambient", 0), ("kent with area factor - ambient", 1)]:
    for name, S_ in [("scattered", sorted(scat)), ("concentrated", sorted(set(per)-scat))]:
        v = np.array([np.mean([p[i] for p in per[s]]) for s in S_])
        print(f"  {lab:32s} {name:12s} {v.mean():+.3f} [{v.min():+.3f},{v.max():+.3f}] {int((v<0).sum())}/{len(v)}")
print("max posterior mass within 12 deg of the degenerate circle <m,x>=0:", worst_mass)
