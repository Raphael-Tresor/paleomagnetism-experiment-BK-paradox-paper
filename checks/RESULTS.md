# Checks of the posteriors against Definition 5 (2026-10-01)

Run from the repository root with `python checks/<script>.py`.

## Posterior on a great circle A under Definition 5

For a metric g' and the shared Fisher prior nu_F (density p w.r.t. the uniform measure),
Definition 5 gives the induced-volume conditional

    nu_F(dx | A)  ∝  (d nu_F / d vol_{g'})(x) dH^1_{g'}(x)
                  =  p(x) (d vol_g / d vol_{g'})(x) (ds_{g'} / ds_g)(x) ds_g(x),

with g the round metric. The experiment writes every posterior against ds_g so that the
predictive densities of the models are comparable.

| analyst | metric g' | posterior w.r.t. ds_g |
|---|---|---|
| ambient | g | p |
| naive (fixed in cd979a8) | flat (inc, dec) chart metric | p sqrt(1 - sin(inc)^2 inc'^2) |
| kent, as coded | Sigma^+ / <m,x>^2 | p sqrt(t^T Sigma^+ t) |
| kent, transported | Kent ellipse transported from m to x | p sqrt(t^T Sigma_x^{-1} t) |

Sigma^+ = e1 e1^T / tau1 + e2 e2^T / tau2 is the rank-two pseudo-inverse built from Kent's
axes at the mean m. With the constant metric Sigma^+ itself, sqrt(det g'_x) =
|<m,x>| / sqrt(tau1 tau2) and the posterior p sqrt(t^T Sigma^+ t) / |<m,x>| is improper.

## `check_definition5.py`: brute-force pre-posterior

Computes exp(-a d'(x,A)^2) nu_F(dx) on a band around A for a = 4e5, with d' obtained by
direct minimisation using the metric at the off-circle point, and bins the mass along A.
Round metric: L1 error 3e-7. For the two anisotropic metrics the brute-force mass matches
each metric's own closed form bin by bin (e.g. sv56, 30-60 deg from the mean: coded
0.585 vs 0.577, transported 0.547 vs 0.544), and not the other one; the residual L1
(1e-1) is discretisation of strongly varying metrics.

## `check_estimated.py`, `check_transport.py`: effect on Table 1 (lower block)

| kent - ambient | scattered sites (sv56, sv61) | concentrated sites |
|---|---|---|
| as coded | -0.236 [-0.407, -0.065] 2/2 | +0.003, 4/9 |
| transported | -0.002 [-0.014, +0.009] 1/2 | +0.000, 5/9 |
| constant Sigma^+ (improper) | +0.132, 1/2 | +0.001, 5/9 |

The gain of the estimated metric at the scattered sites depends on how Kent's tangent
ellipse is extended away from the mean, not on the data.

## Degenerate Kent estimate at sv61

sv61 has 3 directions, so every training set has 2; Kent's second variance is then
tau2 ~ 1e-18 (a degenerate ellipse), clamped to 1e-12 in `posterior_on_circle`.
