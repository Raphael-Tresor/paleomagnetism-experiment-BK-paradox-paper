"""Why is the whitened-vs-ambient effect small? Low power, or no effect?"""
import json, numpy as np
r = json.load(open("results.json"))
rec = r["records"]

print(f"cases {len(rec)}, sites {len(set(x['site'] for x in rec))}")

# Per-case NLPD differences
d = np.array([x["nlpd_whitened"] - x["nlpd_ambient"] for x in rec])
print(f"\nper-case NLPD diff (whitened - ambient):")
print(f"  mean {d.mean():+.4f}  sd {d.std():.4f}  median {np.median(d):+.4f}")
print(f"  whitened better in {(d<0).sum()}/{len(d)} cases ({100*(d<0).mean():.0f}%)")
print(f"  |diff| median {np.median(np.abs(d)):.4f}, max {np.abs(d).max():.4f}")

# Sign test - does whitened win more often than chance?
from scipy.stats import binomtest, wilcoxon
wins = int((d<0).sum()); n = int((d!=0).sum())
bt = binomtest(wins, n, 0.5, alternative="greater")
print(f"\nsign test (whitened wins more often than chance): p = {bt.pvalue:.4f}")
try:
    w = wilcoxon(d, alternative="less")
    print(f"Wilcoxon signed-rank (diff < 0):                p = {w.pvalue:.4g}")
except Exception as e:
    print("wilcoxon failed", e)

# How big is the total spread of NLPD? Context for whether 0.018 is meaningful.
allnl = np.array([x["nlpd_whitened"] for x in rec])
print(f"\nNLPD scale: mean {allnl.mean():.3f}, sd across cases {allnl.std():.3f}")
print(f"effect size (mean diff / sd of diff) = {d.mean()/d.std():.3f}")
print(f"effect as fraction of NLPD sd across cases = {abs(d.mean())/allnl.std():.4f}")

# Site-level clustering: how much does the site bootstrap inflate the CI?
by_site={}
for x in rec: by_site.setdefault(x["site"],[]).append(x["nlpd_whitened"]-x["nlpd_ambient"])
site_means=np.array([np.mean(v) for v in by_site.values()])
print(f"\nsite-level mean diffs: n={len(site_means)}")
print(f"  {np.round(site_means,4)}")
print(f"  mean {site_means.mean():+.4f}  sd {site_means.std():.4f}")
print(f"  sites where whitened better: {(site_means<0).sum()}/{len(site_means)}")
