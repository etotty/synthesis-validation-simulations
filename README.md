# Synthetic data + validation simulations — v0.1.1

A reproducible simulation foundation for studying synthetic data plus constrained validation as research-design infrastructure in economics and social science. This release establishes the confidential-data DGP, five candidate analyses, population counterparts, and a fixed-design oracle MSE benchmark. The repository is the canonical source for code, parameters, seeds, and generated pilot results.

**v0.1.1 restores cross-group heterogeneity in the outcome slope on X.** Analysis 2 remains correctly specified within the target subgroup. Analysis 5 estimates both group treatment effects but imposes a common X slope, gaining precision while introducing misspecification when gamma is nonzero. The committed, code-generated [pilot summary](paper/pilot_results.md) contains the calibration grid, final rankings, and pairwise 2-versus-5 diagnostics. See [simulation notes](paper/simulation_notes.md) for the scientific rationale and calibration ledger.

## DGP and target

Independently across individuals, draw independent `G ~ Bernoulli(p_G)` and `X ~ N(0,1)`, followed by

```
logit Pr(D=1 | X,G) = alpha0 + alpha1*X + alpha2*G + alpha3*X*G
Y = beta0 + beta1*X + beta2*X^2 + betaG*G + gamma*G*X + tau0*D + Delta*D*G + u
u ~ N(0, sigma^2), independent of (X,D,G)
theta = tau0 + Delta
```

All parameters are in [R/config.R](R/config.R):

| Parameter | Value |
|---|---:|
| n | 1500 |
| p_G | 0.12 |
| (alpha0, alpha1, alpha2, alpha3) | (-0.3, 0.4, 0.9, 1.0) |
| (beta0, beta1, beta2, betaG) | (0, 0.7, 0.5, 0.3) |
| gamma | 0.10 |
| tau0 | 0.25 |
| Delta | 0.25 |
| sigma | 1 |
| theta | 0.50 |

The X slope is 0.7 for G=0 and 0.8 for G=1. The treatment effect remains constant within G: there is no `D*X` effect. None of the candidate regressions includes G×X. The treatment-effect parameters restore the original pre-workaround calibration; the prior v0.1 pair is retained as a comparison case.

## Candidate analyses

All regressions include an intercept. The same subgroup effect is the substantive target even for the misspecified candidates.

| Analysis | Sample and regressors | Reported estimate |
|---|---|---|
| 1 | G=1; D, X | D coefficient |
| 2 | G=1; D, X, X² | D coefficient |
| 3 | G=1 and true propensity in [0.10,0.90]; D, X, X² | D coefficient |
| 4 | Full sample; D, X, X², G | D coefficient |
| 5 | Full sample; D, D×G, X, X², G | D coefficient + D×G coefficient |

Analysis 3 temporarily uses the **known DGP propensity**, not an estimated propensity or empirical range intersection. Endpoints are included; the rule uses no outcomes. Analysis 6 is absent.

Each fit saves `T`, HC3 `SE`, `t = T/SE` (the zero-null statistic), and a conventional homoskedastic `SE_ols` for reference. HC3 is used because analyses 1, 4, and 5 may be misspecified under random sampling. Neither estimated SE is substituted for the oracle variance. No coverage claim for the common target is made for biased candidates.

## Decomposition and oracle loss

`theta_a` is the population OLS coefficient/contrast for candidate a, calculated by deterministic numerical integration over X and exact summation over G,D. It is a probability limit, not necessarily the finite-sample expectation. `b_a = theta_a - theta` is systematic specification distortion. Every draw saves `epsilon = T - theta_a` exactly; epsilon need not have zero finite-sample mean.

For realized design `Z=(X,D,G)`, write an estimator as `T = h'Y`. The code builds h using QR decomposition and uses the known conditional mean m to calculate

```
E[T | Z] = h'm
b_a(Z) = h'm - theta
V_a(Z) = sigma^2 * h'h
L(a,Z) = b_a(Z)^2 + V_a(Z)
```

These are exact conditional moments for the specified DGP, including misspecified linear regressions and the analysis-5 contrast. They condition on treatment assignment as well as covariates. They do not include uncertainty from estimating propensity scores or sigma. Failed/rank-deficient fits are explicitly recorded, and no oracle winner is assigned when the candidate menu is incomplete.

## Reproduce

Use **R 4.3.3** for the reference run. `renv.lock` pins **renv 1.0.3**; the simulation and tests otherwise use only base/recommended-with-R components (`stats`, `graphics`, `grDevices`, `utils`, `tools`). R itself and system libraries are not installed by renv. A fresh clone needs internet access for renv's bootstrap and restoration.

From the repository root:

```sh
Rscript -e 'renv::restore(prompt = FALSE)'
Rscript scripts/run_all.R
```

The committed `.Rprofile` and `renv/activate.R` bootstrap the pinned renv version. Do not use `--vanilla` for normal project execution because it bypasses that activation.

Individual stages:

```sh
Rscript tests/run_tests.R
Rscript scripts/00_calibrate_dgp.R
Rscript scripts/01_population_estimands.R
Rscript scripts/02_pilot_montecarlo.R
```

The pilot script first runs a 20,000-noise-draw fixed-design oracle coherence gate, then computes population estimands and runs 2,000 independent confidential samples. Calibration uses 300 separate design seeds per documented case: gamma in {0, 0.05, 0.10, 0.15, 0.20}, first at the v0.1 treatment effects and then at the restored original effects. The same seeds are reused across cases. Seeds, RNG kinds, and run sizes are centralized in `R/config.R`: stage base + replication number makes each replication reproducible independently of execution order. The scripts create and overwrite their generated outputs. `output/` is ignored by Git. The compact, code-generated `paper/pilot_results.md` is committed for review alongside source, tests, configuration, and documentation. No manually entered pilot results are needed.

To compare complete independent runs (e.g., original checkout and a fresh clone):

```sh
Rscript scripts/check_reproduction.R /path/to/first/checkout /path/to/second/checkout
```

This compares byte hashes of all generated CSVs, PNGs, the saved configuration, and the Markdown result table. Exact bytes were checked in the same R/system environment; other BLAS, R, or graphics-library versions may introduce rounding/rendering differences. The statistical tests use numerical tolerances.

## Structure and outputs

```
R/             config, DGP, analysis menu, oracle loss, diagnostics,
               numerical population estimands, pairwise diagnostics, report generation,
               simulation helpers
scripts/       calibration, population calculation, pilot, run-all,
               independent-run comparison
tests/         DGP, menu/contrasts/failures, oracle moments,
               decomposition and process/order reproducibility
renv/          committed activation and settings; local library ignored
output/results/pilot_draws.csv       one row per replication and analysis
output/results/pairwise_designs.csv  one row per design: 2-versus-5 losses/features
output/results/run_config.R          exact parameters and seed configuration
output/results/session_info.txt      runtime/package provenance
output/tables/                      population estimands, calibration,
                                    pilot summary, epsilon covariance,
                                    oracle check, pairwise summary/quantiles/region comparisons
output/figures/                     estimates/loss, design diagnostics, pairwise_tradeoff.png
paper/simulation_notes.md           substantive decisions and limitations
paper/pilot_results.md              generated compact pilot table
```

`pilot_draws.csv` includes estimates, standard errors, population and conditional biases, conditional variance/loss, oracle winners, total/subgroup/treated/control counts before and after trimming, treatment shares, propensity quantiles and overlap summaries, maximum leverage, and a column-scaled condition number. The summary includes Monte Carlo SEs for means, losses, realized MSE, and oracle frequencies. The epsilon covariance uses joint complete replication vectors; failures and the number of eligible oracle designs are reported explicitly.

## Pairwise interpretation

`Delta L25 = L2-L5`: negative favors subgroup analysis 2; positive favors pooled analysis 5. The exact decomposition is `(V2-V5) - (b5(Z)^2-b2(Z)^2)`. Absolute gaps at most 1e-10 are numerical ties; strict sign frequencies are reported separately. The summary reports signed, absolute, and relative gaps, Monte Carlo uncertainty, and design-feature comparisons by region. Relative gap is `abs(L2-L5)/min(L2,L5)`; small-gap threshold counts are descriptive, not selection rules. Overall oracle frequencies still compare all five analyses. Analysis 3 remains weakly dominated by analysis 2.

The generator checks that the saved calibration configuration matches the pilot configuration and refreshes calibration if needed. Raw CSVs, figures, and configuration records stay local; the generated Markdown summary is tracked. Calibration is selected before inspecting disjoint pilot seeds, and no winner-frequency target is imposed.

## Outside v0.1.1

This release does **not** implement `f_Z`, `rho_epsilon`, `lambda`, PAP vs synthetic vs confidential regimes, Proposition 1 or Proposition 2 experiments, validation-budget or adaptive-validation extensions, or synthesis distortion `delta_a`. Those require review and later versions. Capital `Delta` here is treatment-effect heterogeneity, not synthesis distortion.
