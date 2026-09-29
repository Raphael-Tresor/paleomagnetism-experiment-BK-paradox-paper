"""Gate check: is the great-circle geometry / error anisotropy big enough to matter?

If the two candidate metrics (ambient/isotropic vs whitened/anisotropic) give
indistinguishable answers, there is no experiment. Check on the real numbers first.
"""
import numpy as np, csv, sys

def read_magic(path):
    with open(path) as f:
        lines = [l.rstrip("\n") for l in f]
    hdr = lines[1].split("\t")
    rows = []
    for l in lines[2:]:
        if not l.strip():
            continue
        v = l.split("\t")
        v += [""] * (len(hdr) - len(v))
        rows.append(dict(zip(hdr, v)))
    return rows

def dv(dec, inc):
    d, i = np.radians(float(dec)), np.radians(float(inc))
    return np.array([np.cos(i)*np.cos(d), np.cos(i)*np.sin(d), np.sin(i)])

spec = read_magic("data/specimens.txt")
sites = read_magic("data/sites.txt")

lines_by_site, planes_by_site = {}, {}
for r in spec:
    mc = r.get("method_codes", "")
    site = r.get("site") or (r.get("specimen","")[:4])
    if not r.get("dir_dec") or not r.get("dir_inc"):
        continue
    try:
        v = dv(r["dir_dec"], r["dir_inc"])
    except ValueError:
        continue
    if "DE-BFL" in mc:
        lines_by_site.setdefault(site, []).append(v)
    elif "DE-BFP" in mc:
        planes_by_site.setdefault(site, []).append(v)

print(f"specimens: {len(spec)}  sites: {len(sites)}")
print(f"sites with lines : {len(lines_by_site)}")
print(f"sites with planes: {len(planes_by_site)}")

# Site-level Fisher stats and shape of the scatter (isotropy test)
print("\n--- per-site line scatter: eigenvalue structure of the orientation tensor ---")
print(f"{'site':6} {'n':>3} {'k_fisher':>9} {'a95':>7} {'tau2/tau3':>10} {'elong%':>8}")
rows_out = []
for site, vs in sorted(lines_by_site.items()):
    V = np.array(vs)
    n = len(V)
    if n < 4:
        continue
    R = np.linalg.norm(V.sum(axis=0))
    mean = V.sum(axis=0)/R
    k = (n-1)/(n-R) if n > R else np.inf
    a95 = np.degrees(np.arccos(1 - (n-R)/R * ((1/0.05)**(1/(n-1)) - 1))) if n > R else 0.0
    # orientation tensor of residual scatter -> is it circular?
    T = (V.T @ V)/n
    w = np.sort(np.linalg.eigvalsh(T))[::-1]   # tau1 >= tau2 >= tau3
    ratio = w[1]/w[2] if w[2] > 0 else np.inf
    elong = 100*(np.sqrt(w[1]) - np.sqrt(w[2]))/np.sqrt(w[1]) if w[1] > 0 else 0
    print(f"{site:6} {n:3d} {k:9.1f} {a95:7.2f} {ratio:10.2f} {elong:8.1f}")
    rows_out.append((site, n, k, a95, ratio, elong))

if rows_out:
    ratios = np.array([r[4] for r in rows_out])
    elongs = np.array([r[5] for r in rows_out])
    print(f"\ntau2/tau3  median {np.median(ratios):.2f}  max {ratios.max():.2f}")
    print(f"elongation median {np.median(elongs):.1f}%  max {elongs.max():.1f}%")
    print("\nGATE: tau2/tau3 == 1 means circular scatter (isotropic).")
    print("      Values well above 1 mean the error is anisotropic -> the two metrics differ.")
