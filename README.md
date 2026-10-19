# Synthetic data + validation simulations — v0.2

A reproducible simulation foundation for studying synthetic data plus constrained validation as research-design infrastructure in economics and social science. This release establishes the confidential-data DGP, five candidate analyses, population counterparts, and a fixed-design oracle MSE benchmark. The repository is the canonical source for code, parameters, seeds, and generated pilot results.

**v0.1.1 restores cross-group heterogeneity in the outcome slope on X.** Analysis 2 remains correctly specified within the target subgroup. Analysis 5 estimates both group treatment effects but imposes a common X slope, gaining precision while introducing misspecification when gamma is nonzero. The committed, code-generated [pilot summary](paper/pilot_results.md) contains the calibration grid, final rankings, and pairwise 2-versus-5 diagnostics. See [simulation notes](paper/simulation_notes.md) for the scientific rationale and calibration ledger.

**v0.1.2 diagnoses the attainable value of design adaptation.** The code-generated [design-value report](paper/design_value_diagnostic.md) compares nine nearby sample-size/subgroup-prevalence cells using perfect pairwise selection against the best fixed choice. It reports absolute and proportional oracle gains with Monte Carlo intervals, loss gaps, bias/variance components, counts, and failures. The canonical v0.1.1 parameters below are unchanged; alternatives require review.

**v0.2 adds a researcher-observable decision layer.** The researcher sees only `Y,X,D,G` and chooses between analyses 2 and 5 using estimated precision plus an estimated pooling-restriction penalty. Oracle quantities are reserved for evaluation. Read the [information model](paper/researcher_information_model.md) and the code-generated [benchmark results](paper/researcher_information_results.md). The same-sample rule is evaluated separately for selected design quality, actual MSE, systematic distortion and treatment-noise selection. No DGP parameters are recalibrated.

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
# Additional v0.1.2 diagnostic (first-pass grid and independent confirmation):
Rscript scripts/03_calibrate_design_value.R
# v0.2 observed-information benchmark, canonical and sensitivity cases:
Rscript scripts/04_researcher_information.R
```

The pilot script first runs a 20,000-noise-draw fixed-design oracle coherence gate, then computes population estimands and runs 2,000 independent confidential samples. Calibration uses 300 separate design seeds per documented case: gamma in {0, 0.05, 0.10, 0.15, 0.20}, first at the v0.1 treatment effects and then at the restored original effects. The same seeds are reused across cases. Seeds, RNG kinds, and run sizes are centralized in `R/config.R`: stage base + replication number makes each replication reproducible independently of execution order. The scripts create and overwrite their generated outputs. `output/` is ignored by Git. The compact, code-generated `paper/pilot_results.md` and `paper/design_value_diagnostic.md` are committed for review alongside source, tests, configuration, and documentation. No manually entered pilot results are needed.

To compare complete independent runs (e.g., original checkout and a fresh clone):

```sh
Rscript scripts/check_reproduction.R /path/to/first/checkout /path/to/second/checkout
# Compare only v0.1.2 outputs after running script 03 in both clones:
Rscript scripts/check_reproduction.R /path/to/first/checkout /path/to/second/checkout design-value
# Compare only v0.2 outputs after running script 04 in both clones:
Rscript scripts/check_reproduction.R /path/to/first/checkout /path/to/second/checkout researcher
```

This compares byte hashes of all generated CSVs, PNGs, the saved configuration, and the Markdown result table. Exact bytes were checked in the same R/system environment; other BLAS, R, or graphics-library versions may introduce rounding/rendering differences. The statistical tests use numerical tolerances.

The existing `run_all.R` retains its pilot scope. Run `03_calibrate_design_value.R` separately for v0.1.2: 3,000 replications in each of the nine cells, followed by 10,000 per documented confirmation cell. Optional `grid` or `confirmation` arguments run one stage; `report` regenerates the combined report and figure from saved stage tables. Cached stages must match their saved configuration before they can be combined. Both stages and all seeds are defined in `R/config.R`. The diagnostic uses shared uniforms and normal draws, with prefixes for smaller n and common thresholds for subgroup membership and treatment. Replicates are independent; cells within a stage are deliberately coupled. The old pilot simulator and seeds are unchanged. Confirmation uses disjoint seeds and its estimates are reported separately, not pooled with the exploratory grid.

The diagnostic's pointwise 95% Monte Carlo intervals use simultaneous intervals for the gains relative to each fixed analysis, then take their minimum. This handles uncertainty about the best fixed analysis near a crossover. They quantify simulation error, not sampling inference about an empirical population. Reported mean relative gaps are distinct from the proportional adaptation gain.

For v0.2, run script 04 separately; `run_all.R` preserves its historical pilot scope. Script 04 uses 10,000 replications in each of the canonical (1500,.12) and sensitivity (2000,.08) cases. Optional `canonical` or `sensitivity` arguments run one case; `report` regenerates the report/figure from compatible saved summaries. The source config records the cases, rule settings and seed bases 700000 (shared design/outcome draws) and 800000 (independent-noise evaluation). No true DGP quantities enter the observed-data API.

The baseline variance proxy is the augmented model's estimated residual variance times the candidate's squared regression-weight norm. The restriction penalty is `[gamma_hat * h5'(GX)]^2`, exactly equal to the squared restricted-versus-augmented treatment-effect difference. The augmented regression is a diagnostic, not a sixth candidate. Analysis 2 receives zero restriction penalty as a maintained reference-model assumption. A separate HC3 comparison retains the same penalty. Ties within 1e-10 favor 2. See the conceptual document for assumptions and unresolved choices.

Because these estimates use Y, mean assigned oracle loss for the selected analysis need not equal its same-sample MSE. The report includes both, plus a simulation-only independent-noise evaluation holding design and choice fixed. It also separates `epsilon` into a realized-design offset and pure outcome noise. This diagnoses selection rather than treating structural diagnostics as automatically independent of focal results.

## Structure and outputs

```
R/             config, DGP, analysis menu, oracle loss, diagnostics,
               numerical population estimands, pairwise diagnostics, report generation,
               simulation helpers, design-value diagnostic and report
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
paper/design_value_diagnostic.md    generated v0.1.2 grid and confirmation report
output/tables/design_value_grid.csv combined diagnostic estimates and MC uncertainty
output/figures/design_value_grid.png oracle gains with MC intervals
output/results/design_value_*.csv   reproducible paired draws for each cell/stage
R/researcher_information.R          observed regressions and structural evidence
R/perceived_loss.R                  feasible loss/choice; no oracle dependencies
R/researcher_evaluation.R            simulation-only evaluation and summaries
R/researcher_report.R                generated report and figure
paper/researcher_information_model.md   information boundary and assumptions
paper/researcher_information_results.md generated v0.2 benchmark
output/results/researcher_*_draws.csv    paired choices, estimates and evaluation
output/results/researcher_*_menu.csv     full five-analysis diagnostics
output/figures/researcher_information.png design loss versus adaptive MSE
```

`pilot_draws.csv` includes estimates, standard errors, population and conditional biases, conditional variance/loss, oracle winners, total/subgroup/treated/control counts before and after trimming, treatment shares, propensity quantiles and overlap summaries, maximum leverage, and a column-scaled condition number. The summary includes Monte Carlo SEs for means, losses, realized MSE, and oracle frequencies. The epsilon covariance uses joint complete replication vectors; failures and the number of eligible oracle designs are reported explicitly.

## Pairwise interpretation

`Delta L25 = L2-L5`: negative favors subgroup analysis 2; positive favors pooled analysis 5. The exact decomposition is `(V2-V5) - (b5(Z)^2-b2(Z)^2)`. Absolute gaps at most 1e-10 are numerical ties; strict sign frequencies are reported separately. The summary reports signed, absolute, and relative gaps, Monte Carlo uncertainty, and design-feature comparisons by region. Relative gap is `abs(L2-L5)/min(L2,L5)`; small-gap threshold counts are descriptive, not selection rules. Overall oracle frequencies still compare all five analyses. Analysis 3 remains weakly dominated by analysis 2.

The generator checks that the saved calibration configuration matches the pilot configuration and refreshes calibration if needed. Raw CSVs, figures, and configuration records stay local; the generated Markdown summary is tracked. Calibration is selected before inspecting disjoint pilot seeds, and no winner-frequency target is imposed.

## Outside v0.2

This release does **not** implement `f_Z`, `rho_epsilon`, `lambda`, PAP vs synthetic vs confidential regimes, Proposition 1 or Proposition 2 experiments, validation-budget or adaptive-validation extensions, or synthesis distortion `delta_a`. Those require review and later versions. Capital `Delta` here is treatment-effect heterogeneity, not synthesis distortion.
