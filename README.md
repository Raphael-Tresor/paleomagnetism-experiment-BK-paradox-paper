# Palaeomagnetic great-circle experiment

Code and data for the real-data example (Section 5.3 and Supplement Section E) of

> R. Tresor and M. Lukashchuk, *A resolution of the Borel–Kolmogorov paradox via the
> Maximum Entropy Principle*.

A rock carrying two overlapping magnetization components records, under progressive
demagnetization, a direction confined to a great circle of the sphere. The site's
palaeofield direction is therefore conditioned on a set of prior probability zero. The
experiment conditions on these great circles with the MaxEnt posterior of the paper
under several metrics on the sphere, and scores each metric by how well the resulting
posterior predicts held-out measurements (leave-one-out, negative log predictive
density, NLPD).

## Contents

| File | Role |
|---|---|
| `realdata_experiment.py` | Main experiment: loads the data, runs the leave-one-out comparison of the metrics, writes `results.json` and a figure. |
| `check_anisotropy.py` | Per-site scatter of the line directions: Fisher concentration, orientation-tensor eigenvalue ratio and elongation. |
| `check_null.py` | Tests whether the observed elongation is real or a small-sample artifact, against simulated isotropic (Fisher) scatter at the same sample size and concentration (Supplement Section E). |
| `results.json` | Output of a full run of `realdata_experiment.py`, including every leave-one-out record. |
| `data/specimens.txt`, `data/sites.txt` | The data (MagIC format). |
| `data/SOURCE.md` | Provenance, download commands and citation of the data. |
| `data/datamodel/` | MagIC data-model reference files (method codes, column definitions); documentation only, not read by the code. |

## Metrics compared

All analysts share the same prior, data and conditioning great circle, and differ only in
the metric on the sphere. Names in the code and in the paper:

| Code | Paper | Metric |
|---|---|---|
| `naive` | *naive* | flat metric in (declination, inclination) coordinates, i.e. the non-geometric formula (1.2) |
| `ambient` | *ambient* | round (geodesic) metric of the sphere, isotropic |
| `kent` | *estimated* | anisotropic metric from the Fisher–Bingham (FB5) moment estimator of Kent (1982), fitted per site |
| `whitened` | — | anisotropic metric from the raw orientation tensor, per site |
| `pooled` | — | one anisotropic metric pooled over all sites |

## Running

Requires Python 3.9 or later, `numpy` and `matplotlib`.

```bash
pip install numpy matplotlib
python realdata_experiment.py            # full run (about as in the paper)
python realdata_experiment.py --quick    # smoke run: coarse grid, fewer bootstraps
python check_anisotropy.py
python check_null.py
```

The random seed is fixed (`SEED = 20260904`), so a full run reproduces `results.json`.
A full run uses 76 leave-one-out predictions at 11 sites.

## Data

San Francisco Volcanics (Arizona) directions of Tauxe et al. (2003), distributed as the
`lnp_magic` example dataset of PmagPy. Details and download commands in
`data/SOURCE.md`. Please cite

> L. Tauxe et al. (2003). Paleomagnetism of the southwestern U.S.A. recorded by 0–5 Ma
> igneous rocks. *Geochemistry, Geophysics, Geosystems* 4(4), 8802. doi:10.1029/2002GC000343

when using the data.
