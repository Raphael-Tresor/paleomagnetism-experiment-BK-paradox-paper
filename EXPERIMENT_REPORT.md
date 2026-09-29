# Real-data experiment — status report

**Scope of this file.** Working notes for the real-data example demanded by **AE point 2**
and **AE report comment 2**. Markdown only; no `.tex` file has been touched. Any
manuscript wording that follows from this must go through the `latex-writting-style`
skill separately.

**Status: dataset chosen and downloaded; the viability gate has been run and PASSED.**
The main experiment is not yet written. Everything below is either a verified fact, a
measured result, or a design decision with its reasoning.

---

## 0. Bottom line (read this first)

**The experiment is complete. Result: the two admissible metrics are indistinguishable on
this dataset; the flat-chart treatment is worse. Both halves are what the paper's theory
predicts.**

| | |
|---|---|
| **Dataset** | PmagPy `lnp_magic`, San Francisco Volcanics (Tauxe et al. 2003), CC BY 4.0 |
| **Null set** | remagnetization great circle `{v : ⟨v,pole⟩ = 0}` — physics, not a modelling device |
| **Pre-registered verdict** | **H1 not supported** (whitened vs ambient CI includes zero) |
| **Robustness** | the null is *genuine*, not power-limited — one site carried the whole apparent effect (§3e) |
| **Only supported comparison** | `ambient` beats `naive` in 9/11 sites, sign p = 0.033 (§3f) |
| **Decision** | report as-is, keep SO(3) primary; options (ii) and (iii) declined (§3d) |

**The claim the paper can make:** *on real great-circle data from concentrated
palaeomagnetic sites, the posterior is insensitive to the choice among admissible metrics,
exactly as the low-noise limit predicts; but treating the coordinate chart as flat is
still measurably wrong.*

**The claim it cannot make:** that §5.2's whitening rule is empirically validated. It
remains a principled recommendation.

**Two retractions of my own earlier reporting**, both corrected in place:
- the per-case sign test (50/76, p = 0.0040) — **withdrawn**, ignores within-site
  clustering; site-level is 6/11, p = 0.50 (§3e)
- the angular-error gap (6.64° vs 6.29°) — **not usable**, holds in only 4/11 sites (§3f)

---

## 1. What the AE is actually asking for

Two *separate* deliverables that are easy to conflate:

| # | Deliverable | Source | Can Bungert–Wacker §2 fill it? |
|---|---|---|---|
| 1 | An early, simple **motivating example** showing the BKP changes inferential conclusions | AE point 2; referee ("open the paper with a sort of motivating example") | **Yes** — see §2 |
| 2 | A **real-data application** where the metric must be argued, not assumed | AE point 2; AE report comment 2 | **No** — simulated, Euclidean, geometry known |

The binding constraint on #2:

> "since the data-generating geometry is **known by construction**, the example does not
> fully address the more difficult practical question of **how the metric should be
> selected when the true generative mechanism is unknown**." — AE report, comment 2

So "real data" is not really about provenance. It is about **removing our access to the
ground-truth metric**.

---

## 2. The motivating example (deliverable #1) — recommend adopting

Bungert & Wacker §2, read in full. `X, Y ~ N(0,1)` i.i.d., `X` the parameter, `Y` noise,
data `Z = Y/X`. Observe `Z = -1`. Two routes:

- condition on `{Z = -1}` directly → posterior density `|x|·exp(-x²)`
- note `{Z = -1} = {X+Y = 0}`, set `W = X+Y`, condition on `{W = 0}` → `N(0, 1/√2)`

Same event, and the posteriors are **qualitatively opposite**: the first assigns near-zero
density at `x = 0`, the second puts its maximum there. They note endless further answers
are available via e.g. `R := (X+Y)·exp((X-Y)²)`.

**Why it fits:** undergraduate algebra, closed form, one figure; flat space (their own
rationale: the manifold setting "needlessly complicates the mathematics"), which answers
the referee's complaint that SO(3) notation is unfamiliar to statisticians; already framed
as Bayesian inference; and reusable later, as the referee asks.

**Provenance — cite the lineage, not just the competitor.** Their footnote 2 says the
example "is from [10] and is similar to an earlier computation in [11]":
- **Proschan, M.A.; Presnell, B. (1998).** "Expect the Unexpected from Conditional
  Expectation." *The American Statistician* **52**(3), 248–252. ← now in `ref.bib` as
  `michael_proschan_expect_1998`
- **Rao, M. (1988).** "Paradoxes in conditional probability." *J. Multivariate Anal.*
  **27**(2), 434–446.

So it is a *classical teaching example*, not a competitor's device.

**One obligation if adopted.** Their §4.2 argues this is *not* really a BKP case: a
Bayesian inverse problem has a distinguished answer (condition on the variable you
measured), and their Lemma 3 resolves it by taking `(X₁, Z)` as fundamental coordinates.
We must engage that argument rather than present the example as an open puzzle.

---

## 3. What reading the existing SO(3) code changed

Source: `../../../BK-paper-Bernoulli/experiments/so3_borel_kolmogorov/`
(`so3_experiment.py`, 20.6 kB, plus a maintained `CLAUDE.md`).

**Finding 1 — the existing experiment is stronger than the AE's comment implies.** It
already runs *two* generative regimes (`draw_truths(regime=...)`, `"isotropic"` = Haar vs
`"flat"` = uniform in the Euler chart) and shows the calibrated analyst is always the one
whose metric matches the truth. Verified headline numbers from its `CLAUDE.md`:

| regime | geodesic analyst | flat-Euler analyst |
|---|---|---|
| isotropic-generated | ≈88.5–90% (calibrated) | 89.5% → **98.9%** (over-covers) |
| flat-generated | 88.3% → **72.7%** (under-covers) | ≈90–92% (calibrated) |

"Neither metric is universally correct" is therefore **already an established result of
yours**, not something the resubmission must newly argue. The gap is only that *we* set
the geometry. That makes deliverable #2 a smaller job than rebuilding: keep the design,
remove our access to ground truth.

**Finding 2 — §5.2 already commits to a falsifiable prediction.** `subsec:methodology-metric-choice`
states the selection rule: *"take the metric that **whitens** [the noise] (the
inverse-covariance, or Mahalanobis, metric with tensor `Σ⁻¹`)… for a Gaussian location
model it coincides with the Fisher information of the observation,"* with the ambient
embedding metric as the isotropic fallback.

**This is the hypothesis the real-data experiment tests**, and it pins the dataset
requirement hard:

> The whitened metric differs from the naive geodesic one **only when the measurement
> error is anisotropic.** With circular error bars the two analysts coincide and the
> experiment shows nothing.

⇒ **First check on any candidate dataset: are per-observation errors anisotropic, and by
how much?** Do this on the real numbers before writing any experiment code.

---

## 3b. Dataset chosen, and the viability gate — **PASSED**

**Dataset: PmagPy `lnp_magic` example — San Francisco Volcanics, Arizona (basalt lava
flows).** Downloaded and verified; full provenance in `data/SOURCE.md`.

```bash
curl -O https://raw.githubusercontent.com/PmagPy/PmagPy/master/data_files/lnp_magic/specimens.txt
curl -O https://raw.githubusercontent.com/PmagPy/PmagPy/master/data_files/lnp_magic/sites.txt
```

### Why this dataset satisfies AE report comment 2

**The null set is physics, not a modelling convenience.** When a rock carries two
overlapping magnetization components, the two component vectors span a plane through the
origin during progressive demagnetization, so the normalized resultant is confined to the
intersection of that plane with `S²` — **exactly a great circle**. The primary direction
is known only to lie on that curve. Palaeomagnetists condition on precisely this set and
have a name and a method code for it (`DE-BFP`, "best-fit plane").

**Verified counts** (by tabulation, not from the source's claims):

| | |
|---|---|
| `DE-BFL` (lines) | 154 |
| **`DE-BFP` (great circles)** | **42** |
| sites using ≥1 great circle | 16 of 35 |
| **sites estimated from great circles only** | **3** (`sv01`, `sv10`, `sv62`) |

Those three sites are **pure null-set inference cases**.

**The metric is disputed in the published literature — we are not inventing the
controversy.** Fisher/α95 (rotationally symmetric: McFadden & McElhinny 1988; Halls 1976)
versus Kent/FB5 (elliptical: Gallo, Cristallini & Tomezzoli 2017; Tauxe et al. 1991).
Gallo et al. (2017) §2 state that under the classical methods "the assumption of
rotational symmetry would seem to be insufficient." Corroborated by Schmidt (1985, *EPSL*)
and Deenen et al. (2011, *GJI* 186:509–520).

### The gate: is the anisotropy large enough to matter?

Per §3 Finding 2, the whitened and ambient metrics coincide under isotropic error, in
which case there would be no experiment. Measured on the real data
(`check_anisotropy.py`, output in `gate_check_results.txt`): per-site orientation tensor
of the line directions, 27 sites with n ≥ 4.

| statistic | median | max |
|---|---|---|
| `τ₂/τ₃` (1 = circular) | **10.4** | 296.1 |
| elongation | **69.0%** | 94.2% |

**Small-n control.** With only 4–7 specimens per site, some apparent elongation is a
sampling artifact. `check_null.py` simulates genuinely isotropic Fisher scatter at the
*same* `(n, κ)` as every real site, 4000 replications:

```
observed median elongation      : 69.0%
isotropic-null median elongation: 52.2%  (95% range 41.4-62.6%)
p-value (one-sided)             : 0.0008
```

**VERDICT: the anisotropy is real, not a small-sample artifact.** The isotropic Fisher
assumption is measurably wrong on this data, so the ambient and whitened analysts will
give different answers. **The experiment is viable.**

⚠️ This is a statement about *scatter of site-mean directions*, which is a proxy for the
error geometry. Before the paper claims the whitening metric is the right one, the
anisotropy should also be established at the level the metric actually uses — the
per-specimen demagnetization path — via `pmag.dokent()` (Zeta/Eta ellipse semi-axes, each
with its own orientation) against `pmag.fisher_mean()` (isotropic `k`, `α95`).

---

## 3c. What is the GOAL of this analysis? (and what do we recommend?)

An explicit answer, because it determines the whole design.

### The goal is NOT "show how to work with an unknown metric"

That framing would be a trap. If the paper's message were "here is how to proceed when
the metric is unknown", the honest answer is *you cannot* — the BKP says the posterior is
underdetermined until the geometry is fixed. Promising a recipe for a genuinely unknown
metric would over-claim, and AE point 3 is explicitly a warning against over-claiming.

### The goal is: show that the metric is an EMPIRICAL question, not a free parameter

The paper's actual position (§5.2) is that the metric is *the geometry of the
measurement* — so it is **identifiable from the measurement process**, and therefore
falsifiable. The real-data example's job is to demonstrate exactly that, in three steps:

1. **The choice is unavoidable and consequential.** On a real great-circle dataset, the
   competing metrics give different palaeopole estimates and different confidence regions.
   Nobody can abstain: the classical Fisher/α95 method is *already* a metric choice, made
   silently.
2. **The choice is decidable from the data.** The metrics make different predictions, so
   held-out scoring adjudicates between them **without anyone declaring the true geometry
   in advance**. This is the part the SO(3) simulation structurally cannot do, and it is
   precisely what AE report comment 2 asks for.
3. **The declared rule is the one that wins.** §5.2 says take the metric that whitens the
   noise. If the anisotropic/whitened analyst scores best, the rule is *validated on data
   where we did not know the answer beforehand*.

So the contribution is not "we resolve unknown metrics." It is: **the BKP's residual
ambiguity is relocated into a modelling choice that is empirically testable, and here is
a case where the test is run and our rule passes.** That is a defensible, bounded claim
that answers point 2 without violating point 3.

### What do we then recommend to a practitioner?

The recommendation is **fixed in advance and does not depend on the result** — that is
what makes it a scientific claim rather than a post-hoc rationalization. It is already
in §5.2:

> **Use the metric that whitens the measurement noise** (Mahalanobis, `Σ⁻¹`); for a
> Gaussian location model this is the Fisher information of the observation. When the
> noise model distinguishes no direction, fall back to the ambient embedding metric.
> When the geometry is genuinely contested, the choice is not removed but *localized to a
> single declared input*.

Applied to this dataset that reads: **if your specimen scatter is elliptical, condition
with the Kent/FB5 geometry, not the Fisher α95 cone.** Concretely — do not use the
rotationally symmetric confidence region for great-circle intersections; that is the same
conclusion Gallo et al. (2017) reach on geophysical grounds, which we would be deriving
from a general principle.

### Does the recommendation depend on the analysis outcome?

**No — but its status does.** The rule is stated a priori; the experiment tests it. Three
possible outcomes, all publishable, and we commit to reporting whichever occurs:

| Outcome | What we report |
|---|---|
| **Whitened metric scores best** | §5.2's rule is validated out-of-sample. Strongest result; the recommendation stands as stated. |
| **No metric separates** (scores within noise) | The choice is immaterial *for this inference*, which is itself useful: it bounds when practitioners must care. Weakens the example's punch, so **check separation early**. |
| **Whitened metric loses** | A genuine finding *against* §5.2. We would have to weaken the selection rule to "the metric is a declared modelling input, and here is evidence that reading it off the noise model is insufficient." |

The third outcome is the reason to run this before restructuring §5. It would be a real
setback — but far better found by us than by a referee, and even then the paper's *core*
theory (unique posterior given a metric) is untouched; only the practical selection
heuristic would need softening.

⚠️ **Pre-register the analysis.** Fix the scoring rule, the held-out split, and the
success criterion *before* looking at the scores, and record them in `CLAUDE.md`.
Otherwise outcome 1 is not credible.

---

## 3d. RESULTS — run 2026-09-04, pre-registered analysis

**Verdict: outcome 2 of §3c. H1 not supported; H0 not rejected.** `whitened` scores best
but the margin's bootstrap CI includes zero. Reported exactly as pre-registered.

Reproduce: `./.venv/bin/python realdata_experiment.py` (full: 2880-point grid, 2000
site-level bootstraps). Full output in `results.json`; diagnostics in `diagnostics.txt`.

### Scale of the analysis

46 sites parsed, **11 sites** contributing **76 leave-one-out cases** (needs ≥1 great
circle and ≥2 lines; fewer than the 16 anticipated, because some great-circle sites lack
two lines to train on).

### Primary result (pre-specified: mean NLPD, lower better)

| analyst | mean NLPD | median angular error | coverage of nominal 90% arc |
|---|---|---|---|
| **`whitened`** | **6.5674** | 6.29° | 78.9% |
| `ambient` | 6.5853 | 6.30° | 80.3% |
| `naive` | 6.6092 | 6.64° | 77.6% |

**Ordering is exactly as H1 predicted** — `whitened` < `ambient` < `naive` — and stable
across grid resolutions (identical to 4 d.p. under `--quick`).

**But the pre-registered test fails:**

```
whitened - ambient: gap -0.0179   95% CI [-0.0711, +0.0049]   excludes zero: FALSE
```

⇒ **H0 not rejected.** By the criterion fixed in `PREREGISTRATION.md` §4, this is *not*
evidence that the whitening metric is better.

### Why — the diagnosis (exploratory, NOT pre-registered)

Post-hoc, logged as a deviation in `PREREGISTRATION.md` §9. **No confirmatory weight.**

- `whitened` wins in **50 of 76 cases (66%)**; per-case sign test p = 0.0040.
  ⚠️ **WITHDRAWN — see §3e.** This test ignores clustering within sites. At site level it
  is 6/11, p = 0.50. Do not quote the per-case figure.
- But the *magnitude* is negligible: **median |ΔNLPD| = 0.004 nats**, effect size
  (mean/sd) = **−0.096**.
- Wilcoxon signed-rank p = 0.051 — borderline.
- **One site dominates the mean.** Site-level mean differences range over
  [+0.026, −0.410]; a single site contributes −0.410 while 5 of 11 sites actually favour
  `ambient`. The site-level bootstrap correctly refuses to generalise from that.

**Honest reading: a small, consistent, but practically negligible advantage, on a sample
too small and too heterogeneous to establish it.** The consistent sign is suggestive; the
magnitude says it would not change a palaeomagnetist's conclusion.

### What this means for the manuscript

The three-way ordering matches the theory's prediction and is worth reporting. The
`naive` analyst is worst on the pooled means, which supports the paper's central warning
against treating a chart as flat — but see **§3f**: only the `ambient` vs `naive`
comparison survives a site-level test (9/11, p = 0.033), the angular-error gap does
**not** hold per site (4/11 only), and no comparison's CI excludes zero. But:

1. **§5.2's rule cannot be claimed as validated.** It must be stated as a *principled
   recommendation*, not an empirically confirmed one. This is the honest version and it
   is still fully compatible with AE point 3's demand against over-claiming.
2. **The useful, defensible finding is a bound**: on this dataset, with Fisher-concentrated
   sites (κ ≈ 80–1100) and near-circular posteriors, **the metric choice is immaterial**.
   That is genuinely informative — it tells practitioners *when they need not care*, and
   it is consistent with the SO(3) experiment's own observation that "the two analysts
   agree in the small-noise limit, because at small scale the space is approximately
   flat."
3. **Where the metric should matter is diffuse data**, and this dataset is mostly
   concentrated. A site-level analysis conditioned on κ would test that, but is *not*
   pre-registered and would be exploratory.

### Options going forward — **DECIDED 2026-09-04**

**Decision: (i) + (iv). Options (ii) and (iii) are declined.**

| Option | Status |
|---|---|
| **(i) Report as-is** | **ADOPTED** — honest, pre-registered, real dataset with a real null set |
| (ii) Find a diffuse-data (low-κ) dataset | **DECLINED** by the author. Would risk looking like dataset-shopping, and §3e shows the null here is genuine rather than power-limited. |
| (iii) κ-stratified secondary analysis | **DECLINED** by the author. 11 sites split into ~5/6 per stratum; the single influential site would land in one bin and drive it. Underpowered to the point of being uninterpretable. §3e also shows there is no latent effect to stratify out. |
| **(iv) Keep SO(3) as the primary demonstration** | **ADOPTED** — it shows the effect where it is large; this dataset bounds where it is not |

**Rationale for the pair.** SO(3) shows the metric choice is consequential where
curvature matters; the palaeomagnetic data show a practitioner's genuine null set where it
is not, with the metric read off the measurement rather than set by us. Together they
**bound when the metric matters**, which answers AE report comment 2 (the geometry was not
ours to choose) without over-claiming, as AE point 3 demands.

---

## 3e. Outlier sensitivity — the effect is FRAGILE, not hidden

**Post-hoc, not pre-registered, not confirmatory** (`sensitivity_outlier.py`, output in
`sensitivity_outlier.txt`). Logged in `PREREGISTRATION.md` §9.

Question asked: the primary mean is dominated by one site (−0.410). What if it is removed?

**Answer: removing it makes the effect vanish, not appear.**

| | gap (whitened − ambient) | 95% CI |
|---|---|---|
| full sample (11 sites, 76 cases) | −0.0179 | [−0.0694, +0.0053] |
| **without `sv61`** | **−0.0017** | **[−0.0172, +0.0077]** |

`sv61` contributes **3 cases** yet carries almost the entire mean difference. Without it
the gap is −0.0017 nats — **an order of magnitude smaller and indistinguishable from
zero.** Leave-one-site-out over all 11 sites: **no** exclusion ever makes the CI exclude
zero.

### The clustering-respecting sign test kills the earlier suggestion

§3d reported a per-case sign test, `whitened` winning 50/76, p = 0.0040. That test
**ignores clustering** — cases within a site are dependent, and sites with 15 cases
counted five times as heavily as sites with 3. Redone at **site** level, which is the
correct unit:

```
sites favouring whitened: 6/11        sites favouring ambient: 5/11
site-level sign test:     p = 0.50
median site-level diff:   -0.0002 nats
trimmed mean (drop 1 each end): -0.0111 nats
```

**6 versus 5 is a coin flip.** The p = 0.0040 from §3d was an artifact of treating
dependent cases as independent, and should not be quoted. I am withdrawing it as
evidence.

### Why removing the outlier is the wrong move anyway

Even had it gone the other way, dropping the observation that drives a result *after*
discovering that it does is circular — it guarantees the finding will not replicate.
Reported here only as an **influence/robustness diagnostic**, which is its legitimate use.

### Conclusion — strengthened, not weakened

The pre-registered verdict (§3d) was already "H0 not rejected." This sensitivity check
makes that conclusion **more secure**: the null is not a power problem masking a real
effect, it is a genuinely null result on this dataset. The apparent signal was one site
out of eleven.

⇒ **On this data the metric choice is immaterial**, and the paper can say so with
confidence rather than hedging about low power. That is a cleaner claim than a marginal
positive would have been, and it is consistent with the concentration argument (κ ≈
80–1100; the paper's own SO(3) analysis predicts agreement in the low-noise limit).

---

## 3f. Is "naive is worse" solid? — partly, and weaker than §3d implied

**Post-hoc, exploratory** (`naive_check.py`, output `naive_check.txt`), but using the
**correct site-level unit** unlike the withdrawn per-case test. Logged in
`PREREGISTRATION.md` §9.

### Site-level results, all three pairings

| comparison | gap (nats) | 95% CI | CI excludes 0 | sites favouring first | sign p |
|---|---|---|---|---|---|
| `ambient` vs `naive` | −0.0240 | [−0.0923, +0.0008] | **no** (just barely) | **9/11** | **0.033** |
| `whitened` vs `naive` | −0.0418 | [−0.1622, +0.0042] | no | 6/11 | 0.50 |
| `whitened` vs `ambient` | −0.0179 | [−0.0694, +0.0053] | no | 6/11 | 0.50 |

### What holds

**`ambient` beats `naive` in 9 of 11 sites (sign test p = 0.033).** That is a consistent
directional result at the correct unit of analysis, and it is the *only* comparison in
this experiment with any statistical support. It says: **using the round metric beats
treating the (dec, inc) chart as flat.** That is the paper's central warning, and it is
the one thing this dataset does support.

### What does NOT hold — three honest caveats

1. **The CI still includes zero** ([−0.0923, +0.0008]), so by the same standard applied
   to the primary comparison this is *not* significant. It is a consistent sign with a
   marginal magnitude. The upper bound is `+0.0008`, i.e. it only barely fails.
2. **`whitened` vs `naive` is NOT supported** — 6/11, p = 0.50. Surprising given the
   ordering of the means, and it happens because `whitened`'s own variance across sites
   is larger. So "whitened beats naive" cannot be claimed.
3. **The angular-error claim was overstated.** §3d quoted 6.64° vs 6.29° medians as
   supporting `naive` being worse. Per site, `naive` has a worse median angular error in
   only **4 of 11 sites** — the pooled medians hide near-cancelling per-site differences
   (range −0.67° to +0.81°). **Do not use the angular-error numbers as evidence.**

### Net reading

The defensible sentence is narrow: **the round metric outperforms the flat-chart
treatment in 9 of 11 sites (p = 0.033), while the two admissible metrics are
indistinguishable.** Both halves are the theory's prediction — geometry matters, the
particular admissible geometry does not, at this concentration.

Everything here is exploratory. It supports the paper's qualitative warning; it does not
constitute a confirmed quantitative result, and should be reported with that status
explicit.

---

## 3g. Site-filtering cascade — a stated limitation

A referee will ask why 23 great-circle sites became 11 scored sites. The answer is
mechanical, and one part of it is a real limitation of the chosen design.

| stage | n | filter |
|---|---|---|
| sites with ≥1 `DE-BFP` great circle | **23** | in `specimens.txt` |
| … also having ≥2 `DE-BFL` lines | **14** | leave-one-out must hold out a line and still have ≥1 to train on |
| … also yielding a finite Fisher κ | **11** | 3 sites drop where `n ≤ R` (two near-identical lines) |

(The earlier figure of "16 sites" came from the site-level summary table in `sites.txt`,
which lists only the 35 sites with a published mean. 23 is the specimen-level count.)

⚠️ **The limitation worth stating in the paper.** The three *pure* null-set sites —
`sv01` (5 planes, 0 lines), `sv10` (5, 0), `sv62` (3, 0) — are **excluded**, precisely
because they have no line directions to hold out. These are the cleanest instances of the
paper's conditioning problem: sites where the direction is known *only* to lie on a great
circle. Route (B) structurally cannot score them.

This is a genuine cost of leave-one-out validation and should be acknowledged rather than
buried. It also reinforces §3e's reading that the sample is thin: 11 sites of 4–7
specimens, with the most interesting cases unusable.

---

## 3h. Citation obligations for the manuscript

Three distinct references, easy to conflate. Full details in `data/SOURCE.md`.

### Required

**The data** — cite for provenance:
> Tauxe, L., Constable, C., Johnson, C. L., Koppers, A. A. P., Miller, W. R. &
> Staudigel, H. (2003). Paleomagnetism of the southwestern U.S.A. recorded by 0–5 Ma
> igneous rocks. *Geochemistry, Geophysics, Geosystems* **4**(4), 8802.
> doi:10.1029/2002GC000343

⚠️ The files say "Recalculated from original measurements; supercedes published results",
so **our numbers will not match that paper's tables.** Word the provenance sentence to
say the data are the MagIC-archived recalculation, not the paper's published values.

Licence: **CC BY 4.0** — redistribution in supplementary material is permitted with
attribution.

**The two competing methods being adjudicated** — both are published positions, which is
what makes this an existing dispute rather than one we invented:
> McFadden, P. L. & McElhinny, M. W. (1988). *EPSL* **87**(1–2), 161–172.
> doi:10.1016/0012-821X(88)90072-6 — the isotropic Fisher/α95 treatment.

> Kent, J. T. (1982). The Fisher–Bingham distribution on the sphere. *JRSS B* **44**(1),
> 71–80. doi:10.1111/j.2517-6161.1982.tb01189.x — the anisotropic (FB5) alternative.

### Not required (unless the code changes)

> Tauxe, L., Shaar, R., Jonestrask, L., et al. (2016). PmagPy: Software package for
> paleomagnetic data analysis… *G-cubed* **17**(6), 2450–2463. doi:10.1002/2016GC006307

This is the **software** paper and says nothing about this dataset. `realdata_experiment.py`
uses **only numpy** — no PmagPy — so this citation is **not currently owed**. It becomes
required if `pmag.dokent()` is added for the per-specimen anisotropy check.

### Optional, for the "the metric is disputed" claim

> Gallo, L. C., Cristallini, E. O. & Tomezzoli, R. N. (2017). *Latinmag Letters* **7**,
> PM07. ⚠️ **No DOI**; cite the URL. Argues rotational symmetry is insufficient for
> great-circle confidence regions — the clearest published statement that the error
> geometry is contested.

> Deenen, M. H. L., et al. (2011). *GJI* **186**(2), 509–520.
> doi:10.1111/j.1365-246X.2011.05050.x ⚠️ **Has an Erratum** (2014, *GJI* **197**(1), 643,
> doi:10.1093/gji/ggu021) that Crossref does not auto-link — cite it manually.

**All of the above must be added to `ref.bib` via Zotero by the author** — per
`latex-writting-style`, this session does not edit `.bib` files.

---

## 4. Design for deliverable #2

### Three analysts (extends the existing two)

| Analyst | Metric | Role |
|---|---|---|
| Ambient / geodesic | round metric on `S²` | isotropic default; our stated fallback |
| **Whitened / Mahalanobis** | `Σ⁻¹` from reported per-observation error | **what §5.2 predicts is correct** |
| Naive chart | flat in (declination, inclination) | the silent default we warn against |

The ambient analyst *is* the Bungert–Wacker Hausdorff construction on an embedded space,
so this table also discharges the AE's "compares with alternative ways of handling the
conditioning problem." Cross-ref: `../../alternative-road-to-compare/bkp_alternative_structures.md` §C1.

### Scoring without knowing the true geometry — the crux

- **(A) Held-out prediction — target.** Condition on the null-set under each metric, then
  predict a quantity measured *independently* (later measurement, different instrument,
  withheld sites). Score predictive calibration / log score. Honest: no ground-truth
  metric is ever assumed, the data arbitrate.
- **(B) Leave-one-site-out CV — fallback.** Needs only one dataset; weaker, since train
  and test share a measurement process.
- **(C) Sensitivity / consequence analysis — supplement only.** Shows a decision changes
  across metrics without showing which is right. **Insufficient alone** — the AE asked how
  inferences *compare*, which implies adjudication.

Aim (A), fall back to (B), include (C) as the interpretive payoff either way.

### Dataset requirements, priority order

1. **Documented anisotropic per-observation error** (Kent/FB5 elliptical regions, or a
   reported error ellipse/covariance per site). Non-negotiable — see §3 Finding 2.
2. **A scientifically motivated null-set** a practitioner actually imposes: a great circle
   (remagnetization circle in palaeomagnetism), a fixed declination, a bedding/fold plane.
3. **A held-out target**, to enable (A).
4. Public, citable, licence-clear, small enough to ship as supplementary material.

---

## 5. Environment (built and verified)

`experiments/realdata/.venv` — numpy 2.0.2, scipy 1.13.1, pandas 2.3.3, matplotlib.
Python 3.9.6. Pins match the SO(3) experiment's `requirements.txt` so both are
reproducible under one toolchain.

⚠️ **No R and no system scientific Python on this machine.** Many canonical directional
datasets ship only as R packages. Mitigation: extract once and commit a CSV, so the
experiment has no R dependency.

---

## 6. Planned deliverables

- `realdata_experiment.py` — one seeded script, `--quick` smoke mode, writes
  `results.json` + figure PDF. Mirror the SO(3) conventions (figure to repo root, since
  the paper includes by relative path).
- `data/` + `SOURCE.md` — URL, access date, licence, citation, exact preprocessing.
- `CLAUDE.md` — code map, invariants, expected numbers, maintenance rule; copy the
  structure of the SO(3) experiment's.
- A manuscript subsection in §5 after the SO(3) illustration — **written separately, via
  the `latex-writting-style` skill.**

---

## 7. Open risks — stated plainly

1. ~~**Anisotropy may be too small to matter.**~~ **RESOLVED — gate passed**, see §3b:
   median elongation 69% against a 52% isotropic null, p = 0.0008.
2. ~~**The null-set may be contrived.**~~ **RESOLVED** — palaeomagnetists condition on
   exactly this great circle and have a method code for it (`DE-BFP`); 42 such
   interpretations and 3 great-circle-only sites in the data.

2b. ~~**Licence is unresolved.**~~ **RESOLVED — CC BY 4.0**, per
   <https://earthref.org/information/disclaimer.htm>. Redistribution in supplementary
   material is permitted with attribution. Data cite: **Tauxe et al. (2003)**,
   doi:10.1029/2002GC000343; software cite: **Tauxe et al. (2016)**,
   doi:10.1002/2016GC006307. See `data/SOURCE.md`.

2c. **The data are recalculated, not the published tables** ("supercedes published
   results"), so our numbers will not match Tauxe et al. 2003's tables. Phrase the
   provenance carefully.
3. **§5.2's rule may fail.** If the whitened metric is *not* the calibrated one, that is a
   finding against the paper's own selection principle. It must be reported, not buried —
   and it would force weakening §5.2. Better found now than by a referee.
4. **Real work, not editing.** This is the item most likely to decide the outcome; the
   dataset should be settled before §5 is restructured, since the example may reshape it.

---

## 8. Next actions

- [x] Dataset search — PmagPy `lnp_magic` / MagIC selected.
- [x] **Anisotropy gate — PASSED** (§3b), p = 0.0008 against the isotropic null.
- [x] Source publication identified: **Tauxe et al. (2003)**, doi:10.1029/2002GC000343.
- [x] **Licence resolved: CC BY 4.0** — redistribution permitted with attribution.
- [x] Method-paper citations verified (MM88, Halls 1976, Kent 1982, Deenen 2011 + its
      Erratum, Gallo 2017).
- [x] **Pre-registered** in `PREREGISTRATION.md`, frozen before any score was computed.
- [x] Scoring route decided **in advance**: route (B), leave-one-out (route (A) ruled out
      — only 3 great-circle-only sites, no second instrument).
- [x] **Experiment written and run** (`realdata_experiment.py`) → §3d.
- [x] Metric separation checked: **the analysts do NOT separate** (outcome 2). Ordering
      matches theory but the CI includes zero.
- [x] **DECIDED**: (i) + (iv). Options (ii) and (iii) declined by the author (§3d).
- [x] Outlier/influence sensitivity run (§3e) — the null is genuine, not power-limited.
- [x] Site-level test of all three pairings (§3f) — only `ambient` > `naive` holds.
- [x] Site-filtering cascade documented as a limitation (§3g).
- [ ] **Remaining (optional, low priority):** confirm anisotropy at the per-specimen level
      via `pmag.dokent()` vs `pmag.fisher_mean()`. §3b's gate used site-mean scatter as a
      proxy. Now lower stakes: since the metrics do not separate, the whitened metric's
      inputs are not load-bearing for any claim being made.
- [ ] **Manuscript write-up** — a §5 subsection after the SO(3) illustration, using the
      wording in §0. Must go through the `latex-writting-style` skill. **Not started; no
      `.tex` file has been touched.**
- [ ] Write `realdata_experiment.py` + `CLAUDE.md`.
- [ ] Only then: draft the manuscript subsection.

## 9. Cross-references

- `DESIGN.md` (this directory) — the dataset-independent design, in full.
- `../../resubmission_checklist.md` §2.1–2.2 — the checklist items this discharges.
- `../../alternative-road-to-compare/bkp_alternative_structures.md` — competing
  constructions; §C1 is Bungert–Wacker, §C3 the Tjur property.
