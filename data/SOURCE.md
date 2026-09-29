# Data source

## PmagPy `lnp_magic` example dataset — San Francisco Volcanics, Arizona

Basalt lava flows. Palaeomagnetic directions from progressive demagnetization, with
some specimen interpretations recorded as **great circles** (best-fit planes) rather
than lines.

### Download (verified working, no authentication)

```bash
curl -O https://raw.githubusercontent.com/PmagPy/PmagPy/master/data_files/lnp_magic/specimens.txt
curl -O https://raw.githubusercontent.com/PmagPy/PmagPy/master/data_files/lnp_magic/sites.txt
```

Accessed 2026-09-04. Format: MagIC tab-delimited (line 1 = `tab<TAB>table_name`,
line 2 = column headers, data from line 3).

### Verified contents

`specimens.txt` — 503 lines. Method-code tallies (verified by tabulation):

| Code | Count | Meaning |
|---|---|---|
| `DE-BFL` | 154 | Best-fit **line** (a direction) |
| `DE-BFP` | **42** | Best-fit **plane** = great circle |

`sites.txt` — 84 lines, 35 sites. **16 sites use ≥1 great circle**; of those,
**3 are estimated from great circles only**:

| Site | lines | planes | α95 | k |
|---|---|---|---|---|
| `sv01` | 0 | 5 | 6.6 | 286 |
| `sv10` | 0 | 5 | 13.9 | 65.1 |
| `sv62` | 0 | 3 | 38.3 | 308.4 |

Those three are **pure null-set inference cases**: the direction is known only to lie
on the intersection of a plane with `S²`.

Relevant columns — specimens: `dir_dec`, `dir_inc` (degrees), `dir_mad_free`,
`method_codes`. Sites: `dir_dec`, `dir_inc`, `dir_alpha95`, `dir_k`,
`dir_n_specimens_lines`, `dir_n_specimens_planes`, `dir_tilt_correction`.

### Which reference is which — read this before chasing papers

Three different things, often confused:

| Reference | What it is | When you need it |
|---|---|---|
| **Tauxe et al. (2003)**, "Paleomagnetism of the southwestern U.S.A. recorded by 0–5 Ma igneous rocks" | **The data paper.** The field study that collected these rocks. | To understand *where the samples came from* and the geological meaning of the great circles. **Cite this for the data.** |
| **Tauxe et al. (2016)**, "PmagPy: Software package…" | **The software paper.** Says nothing about this dataset. | Only if the analysis *uses PmagPy code*. Our `realdata_experiment.py` uses only numpy, so this is **not currently required** — it becomes required if `pmag.dokent()` is added. |
| MagIC data model / method codes | The **file-format spec**: what `DE-BFP`, `dir_dec`, `dir_alpha95` mean. | To understand *the file you are reading*. This is the practical one. |

**To understand the data itself, none of the papers is the fastest route.** The column
semantics are in the MagIC data model, and the one fact the experiment depends on — that a
`DE-BFP` record's `(dec, inc)` is the **pole** to the plane, not a direction on it — was
verified empirically by `../check_geometry.py` (planes sit 87.9° from their site mean,
lines 5.0°). That check is more reliable than any prose description.

⚠️ **The 2003 paper will not reproduce our numbers.** The files state "Recalculated from
original measurements; supercedes published results", so the values here differ from that
paper's published tables. Cite it for provenance, not for numerical agreement.

### Where the MagIC data model actually lives

There is no single "data model paper". It is a set of machine-readable files in the
PmagPy repository, under `pmagpy/data_model/`. **Downloaded to `datamodel/` here**
(accessed 2026-09-07):

```bash
curl -O https://raw.githubusercontent.com/PmagPy/PmagPy/master/pmagpy/data_model/all_codes.txt
curl -O https://raw.githubusercontent.com/PmagPy/PmagPy/master/pmagpy/data_model/data_model.json
```

- `all_codes.txt` (55 kB) — the **method-code vocabulary** (`DE-BFP`, `LP-DIR-AF`, …).
  JSON despite the `.txt` extension, under a top-level `"definition"` key.
- `data_model.json` (342 kB) — the **column definitions** per table (`specimens`,
  `sites`, …). ⚠️ **Has a UTF-8 BOM**: open with `encoding="utf-8-sig"` or `json.load`
  fails.

Other files in that directory: `method_codes.json`, `controlled_vocabularies*.json`
(several dated snapshots), `MagIC-data-model.txt` (169 kB). The web UI at
`earthref.org/MagIC/data-models` is a JS single-page app and does not serve these to
scripted fetches — use the raw GitHub URLs above.

### Definitions relevant to this experiment (verbatim from those files)

Method codes, from `all_codes.txt`:

| code | official definition |
|---|---|
| `DE-BFL` | "Best fit line" |
| `DE-BFP` | "Best fit plane" |
| **`DE-BFP-G`** | **"Best fit plane: Great circle"** |
| `DE-BFP-S` | "Best fit plane: Small circle" |
| `DE-FM` | "Fisher mean" |
| **`DE-FM-LP`** | **"Fisher mean: Line and planes"** |
| `LP-DIR-AF` | "Directional data: Step-wise alternating field demagnetization" |
| `SO-SUN` | "Sun compass" |
| `SO-CMD-NORTH` | "Correction applied for magnetic declination: True north" |
| `DA-DIR-GEO` | "Direction correction: Adjusted for sample orientation" |

Note the vocabulary distinguishes **great circle** from **small circle** as separate
codes, and has a dedicated code for the Fisher mean *combining lines and planes* — i.e.
the field formalises exactly the conditioning problem the paper studies. Our records carry
the generic `DE-BFP`, not the `-G`/`-S` refinement.

Column definitions, from `data_model.json`:

| column | table | official definition | unit |
|---|---|---|---|
| `dir_dec` | specimens | "Specimen direction in coordinates specified by tilt correction, Declination" | degrees |
| `dir_inc` | specimens | "…Inclination" | degrees |
| `dir_mad_free` | specimens | "Maximum Angular Deviation (MAD) of the free-floating directional PCA fits to the paleomagnetic vector" | degrees |
| `method_codes` | both | "Colon-delimited list of method codes" | — |
| `dir_alpha95` | sites | "Site direction …, Fisher circle" | degrees |
| `dir_k` | sites | "Site direction …, Fisher's dispersion parameter Kappa" | — |
| `dir_n_specimens_lines` | sites | "Number of specimens included in directional calculations based on best-fit lines" | — |
| `dir_n_specimens_planes` | sites | "…based on best-fit planes" | — |
| `dir_tilt_correction` | sites | "Percentage tilt correction applied to the data" | % (0 = geographic, 100 = tilt-corrected) |

⚠️ **The data model does not state whether a `DE-BFP` record's `(dir_dec, dir_inc)` is
the plane's pole or a direction lying on it.** The definition ("Specimen direction…") is
written for the line case and is silent for planes. That is why
`../check_geometry.py` establishes it empirically instead: planes sit **87.9°** from
their site mean, lines **5.0°**, so the stored value is the **pole**. This is the single
most load-bearing fact in the experiment and it rests on that measurement, not on
documentation.

### Source publication — **RESOLVED**

> **Tauxe, L., Constable, C., Johnson, C. L., Koppers, A. A. P., Miller, W. R. &
> Staudigel, H. (2003).** Paleomagnetism of the southwestern U.S.A. recorded by 0–5 Ma
> igneous rocks. *Geochemistry, Geophysics, Geosystems* **4**(4), 8802.
> doi:[10.1029/2002GC000343](https://doi.org/10.1029/2002GC000343)

The `citations` field reads only "This study", which is a legacy MagIC key. It resolves
via four independent checks:

1. MagIC site records for location "San Francisco Volcanics" carry
   `citations = 10.1029/2002GC000343` and `description = "Tauxe et al., 2003"`.
2. A legacy file inside this very directory in the upstream repo,
   `data_files/lnp_magic/zmab0001193tmp02.txt` (line 908), contains the old
   `er_citations` table mapping the literal key `This study` to that DOI.
3. Site count: the upstream file has 47 distinct `sv##` sites; the paper states samples
   were collected from "47 lava flows", of which 35 yield site-mean directions.
4. Crossref confirms all bibliographic fields.

⚠️ **The data are recalculated, not the paper's published tables** — `sites.txt` says
"Recalculated from original measurements; supercedes published results." Cite Tauxe et
al. 2003 for provenance, but **the numbers will not match the paper's tables**; phrase
accordingly.

**MagIC contribution:** ID **19534** (current, v5; supersedes 14854 and 18867),
contribution DOI form `10.7288/V4/MAGIC/19534`, contributor @ltauxe (Scripps
Paleomagnetic Lab). Related: contribution 20079 = Tauxe, Heslop & Gilder (2024), *JGR
Solid Earth* **129**(8), doi:10.1029/2024JB029502 (the PSV10-24 compilation re-using
these data).

### Licence — **RESOLVED: CC BY 4.0, redistribution permitted with attribution**

Authoritative statement at <https://earthref.org/information/disclaimer.htm>, under a
"License" heading, verbatim:

> "All materials posted at EarthRef.org, except where otherwise noted, are licensed
> under the Creative Commons Attribution 4.0 International License."

Link target verified in the raw HTML as `creativecommons.org/licenses/by/4.0/` — so CC BY
4.0, **not** CC0 and **not** CC-BY-NC. Per-contribution JSON-LD also emits
`"license":"https://creativecommons.org/licenses/by/4.0/"`, so the licence attaches to the
data records, not merely to website prose.

⇒ **Shipping these files as supplementary material is permitted, with attribution.**

Two caveats: the "except where otherwise noted" carve-out means the specific contribution
page should be eyeballed once; and `https://earthref.org/terms-of-use` does **not** exist
(returns a File Not Found body under HTTP 200) — the licence lives only on
`disclaimer.htm`.

### Software citation (required if PmagPy is used)

> **Tauxe, L., Shaar, R., Jonestrask, L., Swanson-Hysell, N. L., Minnett, R., Koppers,
> A. A. P., Constable, C. G., Jarboe, N., Gaastra, K. & Fairchild, L. (2016).** PmagPy:
> Software package for paleomagnetic data analysis and a bridge to the Magnetics
> Information Consortium (MagIC) Database. *Geochemistry, Geophysics, Geosystems*
> **17**(6), 2450–2463. doi:[10.1002/2016GC006307](https://doi.org/10.1002/2016GC006307)

PmagPy code licence: 3-clause BSD. (A Zenodo concept DOI 10.5281/zenodo.10478349 exists
but is not linked from the README; the paper above is what the README asks for.)

### Why this dataset

The null set is **physics, not modelling convenience**. When a rock carries two
overlapping magnetization components, progressive demagnetization traces a path whose
two component vectors span a plane through the origin; the normalized resultant is
confined to the intersection of that plane with `S²` — **exactly a great circle**. The
primary direction is known only to lie on that curve.

And the metric is **genuinely disputed in the published literature**, not invented by
us:

| Metric | Advocated by |
|---|---|
| Rotationally symmetric (Fisher, α95 cone) | McFadden & McElhinny (1988); Halls (1976) |
| Elliptical / anisotropic (Kent FB5) | Gallo, Cristallini & Tomezzoli (2017); Tauxe et al. (1991) |

Gallo et al. (2017), §2:

> "In both, Halls (1978) and Mcfadden and McElhinny (1988) methods, Fisher
> distributions are assumed and errors with rotational symmetry around the mean
> direction calculated. As the bias lie along the parallelism of great circles, an
> elongated confidence region should be expected and the assumption of rotational
> symmetry would seem to be insufficient."

Corroborating: Schmidt (1985), "Bias in converging great circle methods," *EPSL*;
Deenen et al. (2011), *GJI* **186**:509–520, arguing Fisher gives the correct mean but
an incorrect circular error estimate.

The discipline's own controlled vocabulary distinguishes the cases the paradox turns
on: `DE-BFP-G` (great circle), `DE-BFP-S` (small circle), `DE-FM-LP` (Fisher mean from
lines and planes).

⚠️ Gallo et al. cite "Halls 1978", but the converging-circles least-squares paper is
**Halls (1976)** (full citation below). Halls (1978) is a different paper.

### Method-paper citations — all verified via Crossref

The two competing treatments of palaeomagnetic great circles:

> **McFadden, P. L. & McElhinny, M. W. (1988).** The combined analysis of remagnetization
> circles and direct observations in palaeomagnetism. *Earth and Planetary Science
> Letters* **87**(1–2), 161–172. doi:10.1016/0012-821X(88)90072-6
> — *the isotropic/Fisher method. This is what PmagPy implements:* `pmag.dolnp()`
> docstring says "Returns fisher mean, a95 for data using the method of McFadden and
> McElhinny 1988 for lines and planes."

> **Halls, H. C. (1976).** A least-squares method to find a remanence direction from
> converging remagnetization circles. *Geophysical Journal of the Royal Astronomical
> Society* **45**(2), 297–304. doi:10.1111/j.1365-246X.1976.tb00327.x

> **Kent, J. T. (1982).** The Fisher–Bingham distribution on the sphere. *Journal of the
> Royal Statistical Society Series B* **44**(1), 71–80.
> doi:10.1111/j.2517-6161.1982.tb01189.x
> — *the anisotropic (FB5) alternative;* `pmag.dokent()` implements it.

> **Deenen, M. H. L.; Langereis, C. G.; van Hinsbergen, D. J. J.; Biggin, A. J. (2011).**
> *Geophysical Journal International* **186**(2), 509–520.
> doi:10.1111/j.1365-246X.2011.05050.x
> ⚠️ **Has an Erratum**: 2014, *GJI* **197**(1), 643, doi:10.1093/gji/ggu021. Crossref
> carries no auto-link, so cite it manually.

> **Gallo, L. C.; Cristallini, E. O.; Tomezzoli, R. N. (2017).** Bootstrapped intersecting
> remagnetization great circles and the subsequent empirical confidence region. *Latinmag
> Letters* **7**, Special Issue, PM07, 1–5.
> ⚠️ **No DOI exists** (Latinmag registers none); conference proceedings, not in
> Crossref/Scopus. Cite the URL:
> <https://www.geofisica.unam.mx/LatinmagLetters/LL17-01-SP/PM/PM07.pdf>

Also relevant: the published exchange **Bailey (1990)** *EPSL* 101:125–126 and
**McFadden (1990)** *EPSL* 101:127–128 debating the MM88 method — useful evidence that the
error-geometry question is genuinely contested in the field.

### Related PmagPy files (also verified present)

- `foldtest/foldtest_example.dat` — n=80: dec, inc, bedding dip-direction, dip.
  A *second* null-set story: conditioning on the bedding-restored plane.
- `gokent/gokent_example.txt` — n=20, for Kent (FB5) fitting.

### MagIC database (broader source, live API, no key)

```bash
curl "https://api.earthref.org/v1/MagIC/search/contributions?query=DE-BFP"  # 91 contributions
curl "https://api.earthref.org/v1/MagIC/data?id=16645" -o contrib.txt
```

7,017 contributions total; 129 hits for "fold test".

### Error models available in PmagPy

- `pmag.fisher_mean()` → `k`, `alpha95`, `csd` — **isotropic**
- `pmag.dokent()` → `Zeta`/`Eta` (major/minor ellipse semi-axes), **each with its own
  dec/inc orientation** — explicitly anisotropic

That contrast is exactly the ambient-vs-whitened metric distinction the experiment needs.
