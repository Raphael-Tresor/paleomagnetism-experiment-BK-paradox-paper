"""POST-HOC SENSITIVITY CHECK -- NOT pre-registered, NOT confirmatory.

Question: the primary result's mean is dominated by one site (-0.410). What happens if
it is removed?

This is leave-one-site-out influence analysis. It CANNOT establish an effect: removing
the observation that drives a result, after seeing that it does, is circular. It can only
answer "is the primary conclusion robust or fragile?"
"""
import json, numpy as np
from itertools import combinations

r = json.load(open("results.json"))
rec = r["records"]
ANAL = ("whitened", "ambient", "naive")

by_site = {}
for x in rec:
    by_site.setdefault(x["site"], []).append(x)
sites = sorted(by_site)

def gap(rows, a="whitened", b="ambient"):
    return float(np.mean([x[f"nlpd_{a}"] for x in rows]) - np.mean([x[f"nlpd_{b}"] for x in rows]))

def boot_ci(rows_by_site, names, n_boot=4000, seed=20260904):
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(n_boot):
        pick = rng.choice(len(names), size=len(names), replace=True)
        rows = [x for i in pick for x in rows_by_site[names[i]]]
        if rows: draws.append(gap(rows))
    d = np.array(draws)
    return float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))

print("=== per-site influence on the primary gap (whitened - ambient) ===")
print(f"{'site':6} {'cases':>5} {'site mean diff':>15} {'gap w/o site':>13} {'CI w/o site':>26}")
print("-"*72)
full_gap = gap(rec)
lo, hi = boot_ci(by_site, sites)
rows_out=[]
for s in sites:
    others = [n for n in sites if n != s]
    rows = [x for n in others for x in by_site[n]]
    g = gap(rows)
    l, h = boot_ci(by_site, others)
    sm = gap(by_site[s])
    rows_out.append((s, len(by_site[s]), sm, g, l, h))
    print(f"{s:6} {len(by_site[s]):5d} {sm:+15.4f} {g:+13.4f}   [{l:+.4f}, {h:+.4f}]{'  <-- CI excludes 0' if (l>0 or h<0) else ''}")

print(f"\nFULL SAMPLE ({len(sites)} sites, {len(rec)} cases): gap {full_gap:+.4f}  95% CI [{lo:+.4f}, {hi:+.4f}]")

# The dominant site
dom = min(rows_out, key=lambda t: t[2])
print(f"\nMost influential site: {dom[0]} (site mean diff {dom[2]:+.4f}, {dom[1]} cases)")
print(f"  gap without it: {dom[3]:+.4f}  95% CI [{dom[4]:+.4f}, {dom[5]:+.4f}]")
if dom[4] > 0 or dom[5] < 0:
    print("  -> CI still excludes zero WITHOUT the outlier")
else:
    print("  -> CI includes zero without it: the effect does NOT survive removal")

# Direction counts, robust view
site_means = np.array([t[2] for t in rows_out])
print(f"\nsites favouring whitened (diff<0): {(site_means<0).sum()}/{len(site_means)}")
print(f"sites favouring ambient  (diff>0): {(site_means>0).sum()}/{len(site_means)}")
print(f"median site-level diff: {np.median(site_means):+.5f}   (mean {site_means.mean():+.5f})")

# Sign test on SITES (not cases) - the clustering-respecting version
from scipy.stats import binomtest
w = int((site_means<0).sum()); n = int(len(site_means))
bt = binomtest(w, n, 0.5, alternative="greater")
print(f"site-level sign test (whitened wins): {w}/{n}, p = {bt.pvalue:.4f}")

# Trimmed mean: drop most extreme site each side
trimmed = np.sort(site_means)[1:-1]
print(f"\ntrimmed (drop 1 each end) mean of site diffs: {trimmed.mean():+.5f}")
print(f"  -> on the NLPD scale this is a {abs(trimmed.mean()):.5f} nat difference")
