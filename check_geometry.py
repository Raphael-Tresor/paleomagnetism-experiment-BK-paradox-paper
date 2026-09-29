"""Is a DE-BFP record's (dec,inc) the POLE to the plane, or a direction ON the circle?

Test: for sites that have BOTH lines and planes, the site mean comes from the lines.
If (dec,inc) of a plane is a pole, then pole.mean should be ~0 (mean lies ON the circle,
i.e. perpendicular to the pole). If it were a direction on the circle, pole.mean ~ 1.
"""
import numpy as np

def read_magic(path):
    lines = [l.rstrip("\n") for l in open(path)]
    hdr = lines[1].split("\t"); out=[]
    for l in lines[2:]:
        if not l.strip(): continue
        v=l.split("\t"); v += [""]*(len(hdr)-len(v))
        out.append(dict(zip(hdr,v)))
    return out

def dv(dec,inc):
    d,i = np.radians(float(dec)), np.radians(float(inc))
    return np.array([np.cos(i)*np.cos(d), np.cos(i)*np.sin(d), np.sin(i)])

spec  = read_magic("data/specimens.txt")
sites = read_magic("data/sites.txt")

site_mean = {}
for r in sites:
    if r.get("dir_dec") and r.get("dir_inc"):
        try: site_mean[r["site"]] = dv(r["dir_dec"], r["dir_inc"])
        except ValueError: pass

def site_of(spec_name):
    return spec_name[:4]

dots_plane, dots_line = [], []
for r in spec:
    mc = r.get("method_codes",""); sp = r.get("specimen","")
    if not r.get("dir_dec") or not r.get("dir_inc"): continue
    st = site_of(sp)
    if st not in site_mean: continue
    try: v = dv(r["dir_dec"], r["dir_inc"])
    except ValueError: continue
    d = abs(float(np.dot(v, site_mean[st])))
    if "DE-BFP" in mc: dots_plane.append(d)
    elif "DE-BFL" in mc: dots_line.append(d)

dp, dl = np.array(dots_plane), np.array(dots_line)
print(f"DE-BFP (planes): n={len(dp)}  |cos angle to site mean|  median={np.median(dp):.3f}  mean={dp.mean():.3f}")
print(f"DE-BFL (lines) : n={len(dl)}  |cos angle to site mean|  median={np.median(dl):.3f}  mean={dl.mean():.3f}")
print()
print(f"planes: median angle to site mean = {np.degrees(np.arccos(np.median(dp))):.1f} deg")
print(f"lines : median angle to site mean = {np.degrees(np.arccos(np.median(dl))):.1f} deg")
print()
if np.median(dp) < 0.35:
    print("VERDICT: DE-BFP (dec,inc) is the POLE to the plane (~90 deg from the mean).")
    print("  -> the great circle = { v on S^2 : v . pole = 0 }.  This is the null set.")
else:
    print("VERDICT: DE-BFP (dec,inc) appears to lie ON the circle, not be its pole.")
