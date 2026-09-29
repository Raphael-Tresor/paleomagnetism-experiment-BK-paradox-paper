"""How solid is the 'naive is worse' finding? Site-level, clustering-respecting.

Not pre-registered as primary (that was whitened vs ambient), so this is exploratory --
but unlike the withdrawn per-case sign test, this uses the CORRECT unit of analysis.
"""
import json, numpy as np
from scipy.stats import binomtest

r = json.load(open("results.json")); rec = r["records"]
by_site = {}
for x in rec: by_site.setdefault(x["site"], []).append(x)
sites = sorted(by_site)

def gap(rows, a, b):
    return float(np.mean([x[f"nlpd_{a}"] for x in rows]) - np.mean([x[f"nlpd_{b}"] for x in rows]))

def boot(names, a, b, n_boot=4000, seed=20260904):
    rng = np.random.default_rng(seed); d=[]
    for _ in range(n_boot):
        pick = rng.choice(len(names), size=len(names), replace=True)
        rows = [x for i in pick for x in by_site[names[i]]]
        if rows: d.append(gap(rows,a,b))
    d=np.array(d); return float(np.percentile(d,2.5)), float(np.percentile(d,97.5))

for a,b in [("ambient","naive"), ("whitened","naive"), ("whitened","ambient")]:
    g = gap(rec,a,b); lo,hi = boot(sites,a,b)
    sm = np.array([gap(by_site[s],a,b) for s in sites])
    wins = int((sm<0).sum()); n=len(sm)
    bt = binomtest(wins,n,0.5,alternative="greater")
    excl = "YES" if (lo>0 or hi<0) else "no"
    print(f"{a:9s} vs {b:9s}: gap {g:+.4f}  95% CI [{lo:+.4f},{hi:+.4f}]  excludes 0: {excl}")
    print(f"{'':22s} site-level: {a} better in {wins}/{n} sites, sign p={bt.pvalue:.4f}, median {np.median(sm):+.4f}")

    # influence: does any single site drive it?
    worst=None
    for s in sites:
        others=[x for x in sites if x!=s]
        l2,h2 = boot(others,a,b)
        if not (l2>0 or h2<0):
            worst = s if worst is None else worst
    if excl=="YES":
        print(f"{'':22s} robustness: {'survives every single-site removal' if worst is None else f'FAILS when {worst} removed'}")
    print()

# Angular error, site level
print("=== median angular error by analyst (deg) ===")
for a in ("whitened","ambient","naive"):
    v=np.array([x[f"angerr_{a}"] for x in rec])
    print(f"  {a:9s} median {np.median(v):.3f}  mean {v.mean():.3f}")
print()
sm = np.array([np.median([x["angerr_naive"] for x in by_site[s]]) - np.median([x["angerr_ambient"] for x in by_site[s]]) for s in sites])
print(f"per-site (naive - ambient) median angular error: naive worse in {(sm>0).sum()}/{len(sm)} sites")
print(f"  values: {np.round(sm,3)}")
