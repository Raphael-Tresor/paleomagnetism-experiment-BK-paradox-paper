"""Is the observed elongation real, or an artifact of n=4..7 per site?

Simulate genuinely isotropic (Fisher) scatter at the SAME n and kappa as each real
site, and compare the elongation distribution. If real >> simulated, anisotropy is real.
"""
import numpy as np
rng = np.random.default_rng(20260904)

def fisher_sample(n, kappa, rng):
    """Sample n unit vectors from Fisher(kappa) about the north pole."""
    u = rng.random(n)
    w = 1 + np.log(u + (1-u)*np.exp(-2*kappa))/kappa
    th = rng.random(n)*2*np.pi
    s = np.sqrt(np.clip(1-w**2, 0, None))
    return np.column_stack([s*np.cos(th), s*np.sin(th), w])

def elong(V):
    T = (V.T @ V)/len(V)
    w = np.sort(np.linalg.eigvalsh(T))[::-1]
    return 100*(np.sqrt(w[1])-np.sqrt(w[2]))/np.sqrt(w[1]) if w[1] > 0 else 0.0

# (n, k) from the real sites
real = [(4,127.0),(5,79.7),(5,285.8),(5,122.5),(7,293.9),(7,151.1),(5,9.5),(6,174.9),
        (4,184.0),(5,703.8),(4,584.2),(7,109.4),(4,103.2),(4,214.2),(5,394.0),(5,384.1),
        (4,2.2),(4,199.4),(5,30.4),(4,409.3),(4,3.0),(5,19.5),(4,569.7),(4,1128.1),
        (4,489.0),(6,172.4),(4,78.9)]
real_elong = [86.1,79.0,54.9,52.8,34.8,34.5,69.0,81.3,84.0,64.1,64.3,80.8,9.6,19.9,
              51.9,29.3,94.2,72.7,88.0,67.2,93.1,73.9,63.1,69.5,40.7,88.3,79.7]

REPS = 4000
sim_medians = []
for _ in range(REPS):
    es = [elong(fisher_sample(n, k, rng)) for n, k in real]
    sim_medians.append(np.median(es))
sim_medians = np.array(sim_medians)
obs = np.median(real_elong)

print(f"observed median elongation      : {obs:.1f}%")
print(f"isotropic-null median elongation: {sim_medians.mean():.1f}%  "
      f"(95% range {np.percentile(sim_medians,2.5):.1f}-{np.percentile(sim_medians,97.5):.1f}%)")
p = (sim_medians >= obs).mean()
print(f"p-value (one-sided)             : {p:.4f}")
print()
if p > 0.05:
    print("VERDICT: elongation is CONSISTENT with isotropic scatter at these small n.")
    print("  -> Apparent anisotropy is a small-sample artifact. The per-site scatter")
    print("     CANNOT by itself justify an anisotropic metric.")
else:
    print("VERDICT: elongation EXCEEDS the isotropic null -> real anisotropy.")
