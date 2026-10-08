# v0.1 simulation decisions and interpretation

## Scope and process

The user-specified applied-micro DGP and five-analysis menu are retained. G and X are independently drawn, and u is independent Gaussian noise with constant variance. No within-group treatment heterogeneity, synthesis, information regimes, adaptive validation, or sixth analysis has been added. The pilot is a foundation check, not an experiment supporting the paper's propositions. The source PDFs are background materials; no unverified claims from them are needed for this implementation.

## Structural result that overrides the informal calibration clues

Analysis 5 includes every term in the conditional mean. Consequently, for every full-rank realized design, its expected subgroup contrast is theta, and its population counterpart is theta. In particular, the earlier illustrative value theta_5 approximately 0.599 cannot arise as this model's population counterpart with theta=0.50. It would require a different model/menu/target or a calculation error; without the earlier code its provenance cannot be determined.

There is also a finite-design dominance result. Extend the subgroup analysis-2 estimator's outcome weights to the full sample by placing zeros on G=0 observations. These weights satisfy the unbiasedness restrictions for the analysis-5 subgroup contrast. Full-sample OLS is the best linear unbiased estimator of this contrast under independent homoskedastic errors. Therefore

```
V_5(Z) <= V_2(Z), and L(5,Z) <= L(2,Z).
```

Equivalently, analysis 5 can be reparameterized as group-specific intercepts and treatment coefficients with common X and X² slopes. Group-0 observations add information about those common slopes. This is ordinary full-sample OLS under correct cross-group restrictions, not hierarchical shrinkage of subgroup treatment effects.

The same argument extends the trimmed estimator's weights with zeros to the untrimmed subgroup and gives `V_2(Z) <= V_3(Z)`. Analyses 2 and 3 have zero conditional bias because effects are constant within G and the mean model is correct. Thus `L(5,Z) <= L(2,Z) <= L(3,Z)`. Trimming improves propensity spread but cannot improve the specified oracle MSE. These inequalities are tested.

Changing n, subgroup prevalence, assignment coefficients, noise variance, or the existing outcome coefficients cannot reverse this ordering under the stated assumptions. We do not distort the oracle calculation, omit interaction terms from the mean, or add hidden heterogeneity to manufacture analysis-2 wins.

## Calibration ledger

Fixed throughout: n=1500, p_G=0.12, alpha=(-0.3,0.4,0.9,1.0), beta=(0,0.7,0.5,0.3), sigma=1, theta=0.50, and trimming [0.10,0.90]. The subgroup assignment logit is 0.6+1.4X, versus -0.3+0.4X for G=0, producing poorer subgroup overlap. The X² coefficient introduces nonlinear confounding. No coefficient was selected to reproduce the informal example numbers.

1. **Initial:** tau0=0.25, Delta=0.25. In the 300-design calibration, analysis 5 won every design. It was also clear analytically that the intended 2-versus-5 variation was impossible.
2. **Intermediate:** tau0=0.32, Delta=0.18. Holding the target fixed, reducing heterogeneity reduces analysis 4's systematic bias while preserving an intentionally misspecified constant-effect candidate. Analysis 5 still won every calibration design.
3. **Final:** tau0=0.33, Delta=0.17. One further small reduction tests whether analysis 4's precision can offset its remaining bias in some designs. This is the only parameter dimension calibrated. The exact case summaries are generated in `output/tables/calibration_summary.csv`, using common design/noise seeds across cases. Final pilot seeds are disjoint from calibration seeds.

The oracle-choice variation, if present, is a **4-versus-5** tradeoff. This is an explicit departure from the desired qualitative 2-versus-5 tradeoff, not evidence that the original objective was achieved. Analysis 4 wins about 14% of the separate calibration designs, so this final calibration is retained. The observed pilot frequencies are documented in the generated results and release verification below.

## Temporary trimming choice

Analysis 3 retains G=1 observations with **true known propensity** in the inclusive interval [0.10,0.90]. This is a transparent oracle overlap benchmark. It is not fitted logistic regression, a common empirical treated/control range intersection, or a claim about what researchers can observe in future regimes. Selection depends only on the design and is identical for realized outcomes and the oracle calculation. A fixed rule gives an unambiguous population counterpart. All removed treated/control counts are saved. Treatment constancy within G preserves the substantive target after trimming.

The cutoff is a documented temporary choice, not an optimized threshold. Discuss whether future work should use estimated propensity scores, how researchers would obtain them, and whether this candidate is scientifically useful given its dominance under v0.1.

## Population counterpart and realization-specific residual

For each candidate, form population moments `A = E[WW' I(selected)]` and `q = E[W m I(selected)]`. Solve `A beta_a = q` and apply the same coefficient contrast used in the sample. Normalizing both moments by selection probability is unnecessary. The algorithm sums exactly over the two binary variables and uses base-R adaptive numerical integration over the normal covariate. The trimmed interval is found by inverting its known linear propensity logit; no discontinuity is left inside that integral.

Default bounds are [-10,10], with absolute and relative integration tolerances 1e-10. A second calculation expands bounds to [-12,12] and tightens tolerances to 1e-12. The saved refinement differences measure numerical stability; reported moment-integration error estimates are not formal bounds on coefficient error. The normal-tail contribution is negligible at this scale, and exact identities theta_2=theta_3=theta_5=theta independently check the result. This approach avoids large artificial populations and population Monte Carlo error.

`theta_a` is the population projection/probability limit. Misspecified finite-sample OLS need not be centered exactly on that limit. Consequently `epsilon_a=T_a-theta_a` is a decomposition by construction, not an assertion that E[epsilon_a]=0 at n=1500. It contains both conditional noise and design-dependent displacement from the population counterpart:

```
epsilon_a = (E[T_a|Z] - theta_a) + (T_a - E[T_a|Z]).
```

The direction of analysis 1's bias also differs from the informal example. By the omitted-variable formula its bias is beta2 times the coefficient on D in the subgroup population regression of X² on D and X. The present assignment coefficients make that coefficient negative. The sign is a consequence of this calibration, not a change in the estimator or target. Analysis 4's population effect is tau0 plus Delta times the D coefficient from projecting D×G on D, X, X², and G; it need not equal a simple subgroup-prevalence-weighted average.

The covariance is calculated from the joint five-estimate vectors on the same confidential samples. The pilot reports mean epsilon as a diagnostic rather than recentering it to zero.

## Oracle formulas and standard errors

Let W denote the selected regression matrix, c the contrast, and `h=W(W'W)^(-1)c`. In code h is formed with a QR solve, avoiding the explicit crossproduct inverse. Then `T=h'Y`, `E[T|Z]=h'm`, `V= sigma² h'h`, and `L=(h'm-theta)²+V`. Known m is recomputed from the DGP parameters, not estimated from noisy Y. The formula remains valid for misspecified analyses: misspecification affects the conditional mean, not Var(u|Z). Gaussianity allows particularly simple repeated-noise Monte Carlo error checks, but the conditional variance formula only requires independent, constant-variance errors.

The reported SE is HC3, and t tests a zero effect. Conventional OLS SE is also saved. HC3 is a standard random-design choice for the misspecified candidates; it need not estimate the oracle's conditional variance. Fixed-design approximation error appears in fitted residuals but is deterministic conditional on Z. Substituting residual-based SE² into the oracle would therefore change the estimand and can incorrectly alter rankings.

The fixed-design check generates 20,000 independent outcome-noise vectors and compares all five empirical means, variances, and MSEs with their formulas. Checks use six analytic Monte Carlo standard errors to avoid brittle simulation assertions. Independent formula-based lm fits and direct normal-equation calculations also check QR weights and contrasts.

## Diagnostics, failures, and conventions

- Counts and shares are saved for the original design and each analysis sample, with subgroup treated/control counts after selection.
- Overlap diagnostics include subgroup propensity quantiles, mean p(1-p) in each group, fraction outside trim bounds, and the width of the intersection of empirical treated/control propensity ranges. The range intersection is descriptive and never used for trimming. There is no bounded population support violation: logistic treatment probabilities are strictly between zero and one.
- Maximum OLS leverage and the 2-norm condition number of the unit-column-norm regression matrix are saved. Columns are not centered. This is a practical conditioning flag, not a design-free scientific threshold or invariant to all recodings. Population-moment condition numbers are unscaled and are separately labeled.
- Rank deficiency or insufficient rows yields an explicit failure object and NA estimates/losses. Oracle frequencies exclude the entire incomplete design, rather than favor candidates that happened to fit. Failures and denominators are reported. A future calibration with failures needs review; exclusion is not an automatic scientific solution.
- HC3 is undefined at unit leverage; it is returned as NA with a separate `se_status` flag, while a valid coefficient and oracle loss are retained. The summary counts SE failures separately.
- Oracle ties within 1e-12 select the lowest index and record the tie count. This has no substantive effect when minima are distinct.
- CSVs preserve joint draws and generated tables include Monte Carlo uncertainty. PNG figures suppress boxplot outlier points for readability; tails remain in raw draws and diagnostic quantiles.

## Before v0.2

The substantive choice is whether to accept the honest 4-versus-5 selection problem in this baseline or explicitly revise assumptions to create the intended 2-versus-5 tradeoff. Possibilities for discussion include cross-group differences in outcome control slopes or other stated misspecification of the pooled model; none is implemented here. Noise or assignment calibration alone cannot repair the structural dominance. Likewise, meaningful oracle gains from trimming require revisiting its rationale or assumptions.

Separately confirm the known-propensity trimming convention, HC3 reporting choice, population-counterpart interpretation, and whether the observed dominance of three menu items leaves a sufficiently rich design problem for the later experiments. Future information regimes must not inadvertently expose the oracle diagnostics as researcher-observable information. No regime experiments should be run until these issues are reviewed.

## Release verification

All four automated test suites passed for the final calibration, covering fixed-seed DGP reproduction, all five fits, independent lm contrasts, trimming boundaries, rank/empty-sample handling, unit-leverage SE handling, analytic and repeated-noise oracle moments, nonnegative loss, exact epsilon decomposition, population integration refinement, structural dominance, execution-order invariance, and fresh-process reproducibility.

The final 2,000-design pilot completed without fit or SE failures. A fresh local Git clone downloaded and installed renv 1.0.3, restored the locked environment, removed the existing generated outputs, and reran the complete pipeline under R 4.3.3. The independent-run comparison verified identical hashes for all 11 generated tables, raw draws, figures, saved configuration, and the Markdown results table. Session provenance is excluded from hash comparison because it describes each runtime. This tests fresh-clone restoration on the same operating system, not cross-platform byte identity.

The population integration refinement changed no reported estimand at a practically relevant precision. Fixed-design Monte Carlo variances agreed with their analytic counterparts to within 0.84%, with all mean/variance/MSE checks passing their Monte Carlo tolerances. See generated CSVs for exact diagnostics and uncertainty. The high subgroup leverage tail is retained in the data rather than discarded or silently winsorized. Version 0.1 stops here for review.

## Repository output policy

Generated simulation data, detailed tables, figures, and runtime records under `output/` are local outputs excluded from version control. The compact, code-generated `paper/pilot_results.md` is committed so results can be reviewed directly on GitHub. Run `Rscript scripts/run_all.R` from a restored clone to recreate them. No result transport, CSV splitting, or data reconstruction step is required. This storage-only change leaves the tested simulation code and parameters unchanged.
