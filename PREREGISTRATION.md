# Pre-registration — palaeomagnetic great-circle metric comparison

**Written before any score was computed.** Frozen 2026-09-04. Every choice below was
fixed in advance; deviations must be recorded in §9 with reasons, not silently applied.

Rationale: the experiment tests a rule the paper already commits to (§5.2). If the
scoring rule or the split were chosen after seeing the scores, a favourable result would
be uninterpretable. See `EXPERIMENT_REPORT.md` §3c.

---

## 1. Hypothesis under test

The paper's §5.2 (`subsec:methodology-metric-choice`) states:

> "The metric is the geometry of the measurement. When a noise model is available, take
> the metric that **whitens** it (the inverse-covariance, or Mahalanobis, metric with
> tensor `Σ⁻¹`)… When no noise model distinguishes directions, the ambient embedding
> metric is the natural default."

**H1 (pre-specified, directional):** on real palaeomagnetic great-circle data, the
posterior obtained by conditioning under the **whitened (anisotropic) metric** predicts
held-out directions better than the posterior obtained under the **ambient (isotropic)
metric**, which in turn beats the **naive chart** posterior.

**H0:** the three analysts' predictive scores are indistinguishable.

We commit to reporting whichever of the three outcomes in `EXPERIMENT_REPORT.md` §3c
occurs, including H1 failing.

---

## 2. Data — frozen

`data/specimens.txt`, `data/sites.txt`, downloaded 2026-09-04 from PmagPy
`data_files/lnp_magic` (San Francisco Volcanics, Arizona). Provenance, licence (CC BY
4.0) and citations in `data/SOURCE.md`. **No further data will be added.**

### Verified geometry (`check_geometry.py`)

For a `DE-BFP` record, the stored `(dir_dec, dir_inc)` is the **pole to the best-fit
plane**, confirmed empirically:

| record type | n | median angle to its site mean |
|---|---|---|
| `DE-BFP` (planes) | 36 | **87.9°** |
| `DE-BFL` (lines) | 142 | **5.0°** |

So the conditioning null set for a plane with unit pole `p` is

  `C(p) = { v ∈ S² : ⟨v, p⟩ = 0 }`,

a great circle — a 1-dimensional submanifold of `S²`, measure zero. The site mean lies on
it to within ~2°. **This is a real instance of the paper's conditioning problem, not an
analogy.**

### Counts (verified twice, via `method_codes` and `direction_type`)

- 154 `DE-BFL` lines, **42** `DE-BFP` great circles
- 16 of 35 directional sites use ≥1 great circle
- 3 sites are great-circle-only: `sv01` (5 planes), `sv10` (5), `sv62` (3)

---

## 3. The three analysts

All three condition the same prior on the same null set and differ **only** in the metric
used to define the conditioning limit. Mirrors the two-analyst design of the existing
SO(3) experiment, plus the competitor construction.

| Label | Metric on `S²` | Provenance |
|---|---|---|
| `ambient` | round (geodesic) metric | our isotropic fallback; **also the Bungert–Wacker Hausdorff construction** on an embedded submanifold |
| `whitened` | Mahalanobis `Σ⁻¹` from the site's estimated anisotropic scatter (Kent/FB5 `Zeta`,`Eta` axes) | **what §5.2 predicts is correct** |
| `naive` | flat in `(dec, inc)` coordinates — i.e. drop the `cos(inc)` Jacobian | the silent default the paper warns against |

`ambient` corresponds to the classical Fisher/α95 treatment (McFadden & McElhinny 1988);
`whitened` to the elliptical treatment (Kent 1982; Gallo et al. 2017). **Both are
published positions**, so we are adjudicating an existing dispute, not inventing one.

---

## 4. Scoring — route (B), leave-one-out

Route (A) (independent held-out measurement) is unavailable: the 3 great-circle-only
sites are too few (n=3) to score, and no second instrument's measurements accompany this
dataset. **We therefore pre-commit to route (B).** This is stated now, not after seeing
results.

### Task

For each of the 16 sites having ≥1 great circle and ≥1 line:

1. Hold out **one line direction** `v*` (the target).
2. From the remaining specimens of that site, form the prior/likelihood inputs.
3. For each great circle `C(p)` at that site, compute the posterior on `C(p)` under each
   of the three metrics.
4. Score how well each posterior predicts `v*`.

Iterate over all held-out lines at all qualifying sites (leave-one-out).

### Primary metric — pre-specified

**Mean negative log predictive density** of the held-out direction under each analyst's
posterior, pushed forward to `S²` via the site's likelihood. Lower is better.
Ties broken by the secondary metrics; the primary decides.

### Secondary metrics (reported, not decisive)

- **Angular error**: geodesic angle between `v*` and each posterior's circular mean.
- **Calibration**: empirical coverage of each analyst's nominal 90% credible arc on
  `C(p)`. Target 90%; mirrors the SO(3) experiment's coverage design.

### Success criterion — pre-specified

H1 is supported iff **`whitened` attains the lowest mean NLPD**, and the gap to
`ambient` exceeds its **bootstrap 95% CI** (2000 site-level resamples, resampling sites
not specimens, to respect clustering).

If the CI includes zero → **H0 not rejected**, reported as "the metric choice is
immaterial for this inference," with the observed effect size and CI given.

### Multiplicity

One primary comparison (`whitened` vs `ambient`). `naive` is reported for context and
carries no inferential weight. No other comparison will be promoted to primary.

---

## 5. Implementation constraints

- Single script `realdata_experiment.py`, seeded (`SEED = 20260904`), `--quick` smoke mode.
- Writes `results.json` (all quantities) and a figure PDF.
- **No `.tex` file is written or edited by the experiment.**
- Pins: numpy 2.0.2, scipy 1.13.1, pandas 2.3.3, Python 3.9.6 — matching the SO(3)
  experiment.
- Reruns must be bit-stable.

---

## 6. Gate checks already passed (before pre-registration)

| Check | Script | Result |
|---|---|---|
| Error anisotropy is real, not a small-`n` artifact | `check_null.py` | median elongation 69% vs 52% isotropic null, **p = 0.0008** |
| `DE-BFP` is a pole ⇒ null set is a great circle | `check_geometry.py` | planes 87.9° from site mean, lines 5.0° |

Both were run before this file was frozen and are reported regardless of outcome.

---

## 7. What would falsify H1

Stated explicitly so the result cannot be reinterpreted after the fact:

- `ambient` or `naive` attaining a lower mean NLPD than `whitened`.
- `whitened` winning by a margin whose bootstrap CI includes zero (⇒ H0).
- The anisotropy failing to reproduce at the per-specimen level (§8), which would
  undermine the whitened metric's inputs even if it scored well.

Any of these gets reported as-is, and §5.2 of the manuscript weakened accordingly.

---

## 8. Known limitations, acknowledged in advance

1. **Route (B), not (A).** Train and test share a measurement process, so this is
   weaker than genuine out-of-sample validation. Will be stated in the paper.
2. **Anisotropy currently established at site-mean level.** The whitened metric wants
   per-specimen error geometry. To confirm via `pmag.dokent()` vs `pmag.fisher_mean()`;
   if it fails there, that is a limitation to report (§7).
3. **Small `n`.** 16 usable sites, 4–7 specimens each. The bootstrap CI is the guard;
   low power is a real possibility and would show up as outcome 2 (H0).
4. **Data are recalculated**, not Tauxe et al. (2003)'s published tables
   ("supercedes published results"), so numbers will not match that paper.
5. **Fisher/Kent are the field's models, not necessarily the truth.** We test which
   *metric* predicts better, not which distribution is correct.

---

## 9. Deviations from this pre-registration

**2026-09-04 — none affecting the primary analysis.** The primary comparison ran exactly
as specified: mean NLPD, site-level bootstrap, 2000 resamples, `whitened` vs `ambient`.
Its verdict stands as pre-registered: **H1 not supported** (CI includes zero).

**Post-hoc analyses added after seeing the result** (`diagnose.py`,
`diagnostics.txt`). These were **not** pre-registered and carry **no confirmatory
weight**; they are reported to characterise *why* the primary test came out as it did,
and must be labelled exploratory in any write-up:

- per-case sign test of whether `whitened` wins more often than chance
  — **subsequently withdrawn**: it ignores within-site clustering. The correct site-level
  version gives 6/11, p = 0.50.
- Wilcoxon signed-rank on the per-case differences (same clustering flaw)
- effect-size and site-level decomposition of the NLPD gap
- **site-level comparison of all three analyst pairings** (`naive_check.py`):
  `ambient` beats `naive` in 9/11 sites (sign p = 0.033) — the only comparison with any
  support; `whitened` vs `naive` is 6/11 (p = 0.50); the pooled angular-error gap does
  not reproduce per site (naive worse in only 4/11).
- **leave-one-site-out influence analysis** (`sensitivity_outlier.py`): the gap is
  −0.0179 with all 11 sites and −0.0017 without `sv61` alone, whose 3 cases carry almost
  the whole difference. No site exclusion makes the CI exclude zero. This *supports* the
  pre-registered null rather than overturning it.

Reporting them is not optional: omitting them after having looked would be selective
reporting. But they cannot convert a null primary result into a positive one.

**2026-09-18 — two exploratory analysts added** (`kent`, `pooled`), in
`realdata_experiment.py`. **Not pre-registered, no confirmatory weight.**

*Why.* §3 of this file fixes the *rule* ("take the metric that whitens the noise") but not
the *estimator* of `Σ`. The implementation used a raw orientation tensor of 3–6 training
directions, which is (a) not the Kent/FB5 estimator this file names in §3 as the
anisotropic position being adjudicated, and (b) badly conditioned at that sample size.
The null result therefore bounded one estimator, not the rule. These two analysts vary the
estimator while holding the rule fixed:

- `kent` — Kent's (1982) FB5 moment estimator, i.e. what `pmag.dokent()` fits. Closes the
  gap between what §3 cites and what the code computed.
- `pooled` — one anisotropy estimated from all sites at once, on the reading that the
  measurement geometry is a property of the instrument and protocol rather than of each
  flow.

*Safeguards.* `ANALYSTS_PREREGISTERED` still decides `h1_supported`; the primary
comparison, the bootstrap and the verdict are untouched, and every pre-registered number
reproduces bit-identically. The figure marks the additions dashed/hatched.

*Known limitation of `pooled`.* Its metric is fitted once on all sites' line directions,
including those later held out, so it enjoys a mild optimism the others do not. It loses
to `ambient` anyway, which makes that conclusion conservative rather than fragile.

*Outcome.* Neither changes the verdict. `kent` matches `whitened` (6.5692 vs 6.5674)
despite shrinking the extreme small-`n` anisotropies, and `pooled` is the only analyst
that loses to `ambient` (+0.0100, 3/11 sites) even though it halves the across-site
variance. Read together: the null survives a better-conditioned estimator, and the
anisotropy is per-flow rather than instrumental.
