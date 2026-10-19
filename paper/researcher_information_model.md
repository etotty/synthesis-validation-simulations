# v0.2: researcher-observable design information

## Scope and information boundary

The canonical DGP and the five-analysis menu are unchanged. This version makes a feasible decision between analyses 2 and 5 using an analysis dataset containing exactly `(Y,X,D,G)`. The sensitivity case changes only n and p_G to the already reviewed (2000,0.08). No calibration is performed.

The researcher observes the realized sample, subgroup and treatment counts, regression matrices, coefficient estimates, residuals, leverage, estimated standard errors and functions of those quantities. The researcher adopts the specified candidate models and an augmented pooling-restriction diagnostic. Knowing which regressions to consider is a maintained research-design assumption; it is not knowledge of their true coefficients or correctness.

The researcher does not receive true gamma or sigma, the true propensity, the DGP conditional mean m, the target's true value theta, population theta_a, conditional bias, conditional variance or oracle loss. The public interface `researcher_choice(data)` accepts exactly four observed columns and rejects extras. The evaluator explicitly projects its simulation data to those columns before calling it. `R/researcher_information.R` and `R/perceived_loss.R` are loaded into an isolated base-R closure environment with no configuration, DGP, population-estimand or oracle functions. Tests enforce that boundary. The precision and tie settings are decision-rule conventions, not DGP parameters.

The full five-analysis menu remains in the simulation evaluator, including the historical known-propensity trimming diagnostic. That true propensity and its overlap summaries are **evaluation-only**. The feasible rule does not use analysis 3 or any truth-based overlap quantity. The augmented model is a diagnostic, not a sixth selectable analysis.

## Precision baseline

For candidate a, observed design W_a and treatment contrast c_a give weights

```
h_a = W_a (W_a' W_a)^(-1) c_a.
```

The implementation uses QR solves instead of forming the inverse. Analysis 2's weights are extended to the full sample with zeros outside G=1. Analysis 5's contrast selects D plus D*G. Thus `T_a=h_a'Y`.

The augmented diagnostic regression is

```
Y ~ 1 + D + D*G + X + X^2 + G + G*X.
```

It supplies `gamma_hat`, the coefficient on G*X, and `s_aug^2 = RSS_aug/(n-rank(W_aug))`. The baseline precision proxies are

```
Vhat_2 = s_aug^2 * h_2'h_2
Vhat_5 = s_aug^2 * h_5'h_5.
```

This uses an ordinary estimated residual variance, never known sigma squared. A common estimate makes the candidates' precision comparison depend on their information matrices and avoids treating the omitted slope's approximation residual in the restricted model as pure noise. It relies on a maintained homoskedastic error approximation and the adequacy of the augmented mean model. The researcher cannot know those assumptions are true. Their validity under the simulation is an evaluation fact, not an input to the rule.

For continuity with the existing reported standard errors, a diagnostic alternative replaces both variance proxies with the candidates' HC3 squared SEs, retaining the same bias estimate and tie rule. It is not selected based on benchmark performance. HC3 can respond to misspecification and sampling noise differently from the fixed-design oracle variance. The same-sample baseline is not replaced by splitting, cross-fitting, shrinkage, or a significance-screening rule.

## Estimated consequence of the pooling restriction

Write W for the restricted analysis-5 matrix, q=GX, and c for the subgroup contrast. Restricted OLS and the augmented normal equations imply, for **any observed outcome vector**,

```
beta_restricted = beta_augmented,W + (W'W)^(-1) W'q * gamma_hat
T_5 - T_augmented = gamma_hat * c'(W'W)^(-1)W'q
                  = gamma_hat * h_5'q.
```

Define `c5_GX=h_5'q` and `bhat_5=gamma_hat*c5_GX`. The loading uses only observed covariates and treatment assignment. This identity is checked against independent regressions in tests and throughout the benchmark.

Under the maintained simulation DGP the remaining mean terms lie in W, so the evaluator can additionally verify `b_5(Z)=gamma*c5_GX`. The researcher cannot substitute that identity's unknown gamma. `gamma_hat` is noisy evidence about a cross-group slope restriction; `bhat_5` is its estimated treatment-effect consequence relative to the augmented reference model. It is not knowledge of true bias relative to theta.

For G=1 the candidate-2 regression already includes X and X squared. The decision rule assigns it no estimated restriction penalty because it serves as the flexible subgroup reference. This is a maintained asymmetry of the decision problem, not researcher knowledge that bias is zero. An omitted subgroup nonlinearity, unobserved confounding, heterogeneous treatment effects within G, or another invalid restriction could undermine it; none is added in this version.

## Perceived losses and choice

```
perceived_loss_2 = Vhat_2
perceived_loss_5 = Vhat_5 + bhat_5^2
```

Choose 2 if `perceived_loss_2 <= perceived_loss_5 + 1e-10`; otherwise choose 5. The absolute tolerance favors the flexible reference model in numerical ties. Failed/rank-deficient observed fits or invalid precision estimates produce an explicit failure, never an artificial zero loss. The benchmark stops if a required fit fails rather than silently exclude that replication.

The researcher choice does not consult oracle loss. Oracle evaluation uses exact pairwise minima, with exact ties assigned to 2. All five candidates are still computed and saved for diagnostics.

The plug-in square has an important interpretation. Under the augmented model,

```
E[bhat_5^2 | Z] = b_5(Z)^2 + c5_GX^2 Var(gamma_hat | Z).
```

Uncertainty about the restriction therefore adds a positive term to the expected perceived penalty. Version 0.2 diagnoses this behavior rather than silently subtracting the estimation variance, truncating a corrected penalty, imposing a prior, or using true gamma. Whether the decision objective should use a plug-in, uncertainty-adjusted, or posterior penalty remains open for review.

## Same-sample outcome noise and selection

Both gamma_hat and the residual variance depend on Y. Direct invariance to adding a D or D*G mean term does not make the decision independent of treatment-result noise. Let `eta_a=T_a-E[T_a|Z]`, let `r_q` be the residual from regressing GX on W5, and let `h_gamma=r_q/(r_q'r_q)`. Under the simulation's independent homoskedastic errors,

```
Cov(eta_5, gamma_hat-gamma | Z) = sigma^2 h_5'h_gamma = 0
Cov(eta_2, gamma_hat-gamma | Z) = -c5_GX * Var(gamma_hat | Z).
```

The first identity follows from r_q being orthogonal to W5. For the second, the subgroup weights satisfy `h_2'q=0` and reproduce the analysis-5 contrast on W5. Both identities are computationally checked. With positive gamma and c5_GX, choosing pooling when gamma_hat is close to zero can exclude subgroup estimates with positive treatment noise, producing downward selected noise. The sign and size are empirical diagnostics, not a built-in selection target. A square-based rule need not select more positive or more significant treatment effects.

The benchmark reports correlations of the choice indicator with epsilon_2, epsilon_5, eta_2, eta_5 and their difference; gamma estimation variability and error correlations; estimated-versus-oracle bias/loss-gap tracking; and selected distortion with Monte Carlo SEs.

Crucially, the existing decomposition is retained:

```
epsilon_a = T_a-theta_a
          = [E[T_a|Z]-theta_a] + eta_a.
```

Mean selected epsilon is the requested realization-specific component, but it includes design-dependent displacement from the population counterpart. It should not automatically be called pure outcome-noise selection. The report separately gives systematic distortion `E[theta_choice-theta]`, mean epsilon, mean design offset, mean eta, and their sum. The paired mean `eta_choice-eta_2` isolates the change in noise contribution relative to always choosing 2 and improves Monte Carlo precision.

## Assigned design loss is not adaptive-estimator MSE

For fixed a, oracle `L_a(Z)` is its exact conditional MSE. Once the choice depends on Y,

```
D_feasible = E[L_choice(Z)]
```

is an evaluation of which design was selected; it is generally **not** `E[(T_choice-theta)^2]`. Selection changes the distribution of the treatment errors within each chosen region. The benchmark reports both quantities and their paired difference. Calling assigned loss the MSE of the adaptive same-sample estimator would conceal precisely the result-selection channel under investigation.

A simulation-only independent-noise evaluation keeps Z and the already-made choice fixed and computes each selected treatment estimate with a new independent u vector. It does not give the researcher new information or change any choice. Its expected MSE equals the assigned design loss, and its mean selected eta is zero. This is a diagnostic for same-sample selection, not synthetic data, sample splitting, cross-fitting, validation or a new information regime.

## Performance and uncertainty

Benchmarks are always 2, always 5, exact oracle design selection, feasible baseline and the HC3 diagnostic. Report their mean assigned losses, actual MSE, distortions and choice frequencies. Agreement/confusion matrices are descriptive, not substitutes for loss comparisons.

```
D_fixed  = min(E[L2], E[L5])
D_oracle = E[min(L2,L5)]
raw_gain = D_fixed-D_feasible
fraction_captured = raw_gain/(D_fixed-D_oracle).
```

The fraction is not clipped. Zero denominator gives NA. Paired Monte Carlo errors are reported for raw loss differences and a delta-method ratio MCSE is supplied, with a flag if the denominator is within three MCSEs of zero. The best fixed branch is determined by mean loss, never by per-design minima. At a near tie its branch-local MCSE would require caution. The canonical adaptation ceiling is small, so raw loss and gain are the main evidence.

Negative fractions mean feasible selection is worse than the best fixed rule. Although the scalar calculation allows values above one, coherent evaluation on the same replications satisfies `L_choice >= min(L2,L5)` pointwise; an above-one fraction in this benchmark would indicate inconsistent evaluation rather than ordinary sampling noise. No result is clipped to conceal this check.

## Reproduction and scope

Configuration, cases and seeds live in `R/config.R`. Each case uses 10,000 replications. Common uniforms/normals of maximum length 2000 are drawn using seed 700000+r; prefixes and common thresholds preserve both DGP marginals while coupling the cases. Independent-noise evaluation uses seed 800000+r. The original pilot and v0.1.2 streams are unchanged. Run `Rscript tests/run_tests.R`, then `Rscript scripts/04_researcher_information.R` from a restored checkout. Optional `canonical` or `sensitivity` arguments run one case; `report` regenerates the combined Markdown/figure from compatible saved summaries.

The generated results report and this conceptual document are tracked. Detailed CSVs, figures, configuration records and session provenance remain local ignored outputs. No synthetic datasets, fidelity f_Z, result leakage rho_epsilon, sensitivity lambda, PAP/synthetic/confidential regimes, validation, synthesis distortion, new substantive heterogeneity or automated model search are implemented.

## Benchmark assessment

The [generated benchmark](researcher_information_results.md) finds that the straightforward feasible rule has higher assigned design loss than always choosing analysis 2 in both cases. The estimated slope interaction is approximately centered on its true evaluation value but is noisy relative to its magnitude. Estimated bias and perceived loss differences track their oracle counterparts poorly. The squared plug-in penalty is substantially inflated by estimation variance; this does not prevent over-selection of analysis 5 in draws where the estimated interaction happens to be small. These facts concern different features of the penalty distribution and are not contradictory.

Same-sample treatment-noise selection is substantial and downward in both cases. Its paired Monte Carlo error is small relative to the shift. Positive systematic specification distortion and negative selected noise partly cancel, so the small total bias must not be interpreted as absence of either mechanism. Gamma-estimation noise is negatively correlated with subgroup treatment noise, consistent with the conditional covariance calculation above. When a downward gamma-estimation error makes the estimated restriction look small, pooling can replace a subgroup estimate carrying positive treatment noise. Retaining the subgroup model when the estimated penalty is larger tends to retain negative treatment noise instead. This explains the direction without attributing an intentional result preference to the researcher. The HC3 diagnostic gives similar conclusions and does not resolve this issue.

Actual same-sample MSE is reported with its own Monte Carlo error. The canonical difference from assigned loss is modest relative to that error, and the sensitivity difference is less precisely resolved. The clear evidence for a selected mean-noise shift should not be conflated with precise evidence for the magnitude of an MSE change. The independent-noise check is compatible with the assigned-loss benchmark within Monte Carlo uncertainty and removes the selected mean-noise shift.

The observable layer is ready for inspection as a reproducible decision model. It is not yet a demonstrated improvement over a fixed analysis or a clean channel carrying only structural information. Review the treatment of noisy restriction evidence and its same-sample coupling before introducing synthetic fidelity. No baseline rule, penalty correction, splitting procedure or DGP change is adopted in response to these findings.

## Issues requiring review before the next stage

1. Whether the plug-in squared estimated restriction effect is the intended perception rule, or whether uncertainty should be treated differently. These are different decision models, not harmless numerical fixes.
2. Whether a shared augmented residual variance or HC3 better represents the intended researcher's precision reasoning, and which homoskedasticity/model-adequacy assumptions should be maintained.
3. Whether the flexible/reference role of analysis 2 is justified for the eventual applied problem. Its zero assigned penalty is an assumption, not observed truth.
4. Whether same-sample structural diagnostics create enough treatment-noise selection that a separate mechanism or diagnostic protocol is needed. Sample splitting/cross-fitting is not adopted without review.
5. Which evaluation quantity is central: selected design quality, adaptive-estimator MSE, systematic distortion or pure noise selection. They can move in different directions.
6. How an eventual observable information layer should preserve uncertainty and the coupling between structural evidence and focal results. This version does not claim those can be independently tuned.
