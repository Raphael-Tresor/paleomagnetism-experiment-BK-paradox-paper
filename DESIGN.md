# Real-data experiment — design note

Addresses **AE point 2** and **AE report comment 2** (`resubmission_checklist.md` §2.2).
Dataset choice pending; this note fixes the *design*, which is dataset-independent.

## The objection we must answer

> "since the data-generating geometry is **known by construction**, the example does not
> fully address the more difficult practical question of **how the metric should be
> selected when the true generative mechanism is unknown**." — AE report, comment 2

The existing SO(3) experiment (`../../../BK-paper-Bernoulli/experiments/so3_borel_kolmogorov/`)
is *not* weak: it already runs both generative regimes and shows the calibrated analyst is
the one whose metric matches the truth, i.e. **neither metric is universally correct**.
Its limitation is only that *we* set the geometry. So the new experiment must keep that
design and remove our access to the ground truth.

## What the paper already commits to

§5.2 (`subsec:methodology-metric-choice`, l. 1138) states the selection principle:

> "The metric is the geometry of the measurement. When a noise model is available, take
> the metric that **whitens** it (the inverse-covariance, or Mahalanobis, metric with
> tensor `Σ⁻¹`)… for a Gaussian location model it coincides with the Fisher information
> of the observation. When no noise model distinguishes directions, the **ambient
> embedding metric** is the natural default."

**This is the hypothesis the real-data experiment must test.** It also pins the dataset
requirement: the whitening metric differs from the naive geodesic one **only when the
measurement error is anisotropic**. An isotropic-error dataset cannot separate the two
analysts, so the experiment would have nothing to show.

## Dataset requirements (in priority order)

1. **Anisotropic measurement error, documented per observation.** Without this the two
   candidate metrics coincide. On the sphere: Kent/FB5-type elliptical confidence regions
   rather than Fisher circular ones, or a reported error ellipse / covariance per site.
2. **A scientifically motivated null-set** that a practitioner actually wants to impose —
   not one invented for the paper. Candidates: a great circle (remagnetization circle in
   palaeomagnetism), a fixed declination, a bedding/fold plane.
3. **A held-out ground truth** so calibration can be scored without us choosing the
   geometry. This is the crux — see below.
4. Public, citable, licence-clear, and small enough to ship as supplementary material.

## The three analysts

Mirrors the existing two-analyst design, plus the competitor:

| Analyst | Metric | Status |
|---|---|---|
| **Ambient/geodesic** | round metric on `S²` | the isotropic default; our fallback rule |
| **Whitened/Mahalanobis** | `Σ⁻¹` from the reported per-observation error | **what §5.2 predicts is correct** |
| **Naive chart** | flat in (declination, inclination) coordinates | the silent default we warn against |

The Hausdorff/canonical construction of Bungert & Wacker is the ambient-metric analyst on
a Euclidean/embedded space, so this comparison also answers the AE's "compares with
alternative ways of handling the conditioning problem" — see
`../../bkp_alternative_structures.md` C1.

## How to score without knowing the geometry — the key design problem

Three options, best first:

- **(A) Held-out prediction.** Condition on the null-set using each metric, then predict a
  quantity measured *independently* (a later measurement, a different instrument, a
  withheld subset of sites). Score by calibration of the predictive interval / log score.
  **This is the honest version**: no ground-truth metric is ever assumed, and the data
  arbitrate. Requires a dataset with a genuine held-out target.
- **(B) Cross-validation within the dataset.** Leave-one-site-out: condition on the
  constraint using the other sites, predict the held-out site's direction. Weaker than
  (A) — the same measurement process generates train and test — but needs only one
  dataset and is standard practice.
- **(C) Sensitivity/consequence analysis.** Report how much a *decision* changes across
  metrics (e.g. the credible region for a palaeopole, or whether a fold test passes).
  No scoring, but it demonstrates the choice is material. **Weakest — use only as a
  supplement**, since the AE asked for a comparison of inferences, and (C) shows they
  differ without showing which is right.

Aim for **(A)**, fall back to **(B)**, and include **(C)** as the interpretive payoff
either way.

## Deliverables

- `realdata_experiment.py` — one script, seeded, `--quick` smoke mode, writes
  `results.json` + a figure PDF to the repo root (mirroring the SO(3) convention where
  the paper `\includegraphics` by relative path).
- `data/` — the raw download plus a `SOURCE.md` recording URL, access date, licence,
  citation, and the exact preprocessing applied.
- `CLAUDE.md` — code map, invariants, expected numbers, maintenance rule. Copy the
  structure of the SO(3) experiment's `CLAUDE.md`.
- A new manuscript subsection in §5, placed after the SO(3) illustration.

## Environment

`.venv/` here (git-ignore it): numpy 2.0.2, scipy 1.13.1, pandas 2.3.3, matplotlib.
Python 3.9.6, matching the SO(3) experiment's pins. No R available on this machine, so
prefer datasets reachable without R — or extract once to CSV and commit the CSV.

## Honest risks

- **The anisotropy may be too small to matter.** If reported error ellipses are nearly
  circular, the whitened and ambient analysts coincide and the experiment shows nothing.
  **Check this first, on the real numbers, before building anything.**
- **The null-set may be contrived.** If no practitioner conditions on the set we pick, the
  referee's "materially different conclusions" test is not met.
- **§5.2's prediction may fail.** If the whitened metric is *not* the calibrated one, that
  is a finding about the paper's own selection rule and must be reported, not buried. It
  would force a weakening of §5.2 rather than a stronger claim — better to know now.
