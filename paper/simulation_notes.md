# v0.1.1 simulation decisions and interpretation

## Scope and process

The applied-micro DGP now includes the explicitly requested gamma*G*X outcome term; the five-analysis menu is retained. G and X are independently drawn, and u is independent Gaussian noise with constant variance. No within-group treatment heterogeneity, synthesis, information regimes, adaptive validation, or sixth analysis has been added. The pilot is a foundation check, not an experiment supporting the paper's propositions. The source PDFs are background materials; no unverified claims from them are needed for this implementation.

## Why the revision restores the intended design problem

The original v0.1 prompt omitted cross-group outcome-slope heterogeneity that, as clarified by the user in the v0.1.1 request, had been present in earlier exploratory work. With that omission, analysis 5 contained every term in the conditional mean. Its population counterpart and fixed-design expectation were exactly theta. The earlier example theta_5 near 0.599 was therefore incompatible with the written v0.1 DGP, but a nonzero gamma*G*X term supplies the missing source of misspecification.

The scientific choice is between estimating the outcome relationship in a small target subgroup and pooling observations under cross-group restrictions that improve precision but may be imperfect. Version 0.1.1 restores that distinction by allowing X's outcome slope to differ by group. This is not D*X treatment heterogeneity: the potential-outcome difference remains tau0+Delta*G, with common target theta=tau0+Delta=0.50 for every candidate. No candidate regression is changed to include G*X.

For G=1, gamma*G*X is simply gamma*X. Analyses 2 and 3 already include X and X², so both remain exactly conditionally unbiased. Analysis 1 still omits X²; adding a component proportional to its included X regressor does not change its treatment estimate or bias at fixed theta/design/noise. Analysis 5 restricts the X slopes to be common across groups. Its conditional bias is now exactly `gamma * h5'(G*X)`, because all remaining mean terms are included in its regression.

## Original dominance result: the gamma=0 special case

There is also a finite-design dominance result. Extend the subgroup analysis-2 estimator's outcome weights to the full sample by placing zeros on G=0 observations. These weights satisfy the unbiasedness restrictions for the analysis-5 subgroup contrast. Full-sample OLS is the best linear unbiased estimator of this contrast under independent homoskedastic errors. Therefore

```
V_5(Z) <= V_2(Z), and L(5,Z) <= L(2,Z).
```

Equivalently, analysis 5 can be reparameterized as group-specific intercepts and treatment coefficients with common X and X² slopes. Group-0 observations add information about those common slopes. This is ordinary full-sample OLS under correct cross-group restrictions, not hierarchical shrinkage of subgroup treatment effects.

The same argument extends the trimmed estimator's weights with zeros to the untrimmed subgroup and gives `V_2(Z) <= V_3(Z)`. Analyses 2 and 3 have zero conditional bias because effects are constant within G and the mean model is correct. When gamma=0, `L(5,Z) <= L(2,Z) <= L(3,Z)`; at nonzero gamma only `L(2,Z) <= L(3,Z)` is guaranteed. Trimming improves propensity spread but cannot improve the specified oracle MSE. The gamma=0 dominance and the general trimming inequality are tested separately.

At gamma=0, changing n, subgroup prevalence, assignment coefficients, noise variance, or the other outcome coefficients cannot reverse this ordering. With nonzero gamma, the variance ordering V5<=V2 still holds because the design matrices and noise variance are unchanged. The loss comparison instead becomes

```
L2 - L5 = (V2 - V5) - gamma^2 * [h5'(G*X)]^2.
```

Both the precision gain and the omitted-slope distortion vary across realized designs. This is the intended mechanism for switching; the conditional-loss definition is unchanged.

## Historical v0.1 calibration ledger

Fixed throughout: n=1500, p_G=0.12, alpha=(-0.3,0.4,0.9,1.0), beta=(0,0.7,0.5,0.3), sigma=1, theta=0.50, and trimming [0.10,0.90]. The subgroup assignment logit is 0.6+1.4X, versus -0.3+0.4X for G=0, producing poorer subgroup overlap. The X² coefficient introduces nonlinear confounding. No coefficient was selected to reproduce the informal example numbers.

1. **Initial:** tau0=0.25, Delta=0.25. In the 300-design calibration, analysis 5 won every design. It was also clear analytically that the intended 2-versus-5 variation was impossible.
2. **Intermediate:** tau0=0.32, Delta=0.18. Holding the target fixed, reducing heterogeneity reduces analysis 4's systematic bias while preserving an intentionally misspecified constant-effect candidate. Analysis 5 still won every calibration design.
3. **Final:** tau0=0.33, Delta=0.17. One further small reduction tests whether analysis 4's precision can offset its remaining bias in some designs. This is the only parameter dimension calibrated. These historical cases remain documented in the v0.1 source history; v0.1.1 uses the new gamma grid described below. Final pilot seeds are disjoint from calibration seeds.

The historical v0.1 oracle-choice variation was a **4-versus-5** tradeoff. This is an explicit departure from the desired qualitative 2-versus-5 tradeoff, not evidence that the original objective was achieved. Analysis 4 wins about 14% of the separate calibration designs, so v0.1 retained that calibration. Its final pilot shares were 16.25% for analysis 4 and 83.75% for analysis 5.

## v0.1.1 calibration ledger and selection

The initial grid is gamma in {0, 0.05, 0.10, 0.15, 0.20}, at the v0.1 values tau0=0.33 and Delta=0.17, with all other parameters fixed. The same 300 seeds are used in every case. At gamma=0.10, the pairwise 2-versus-5 crossover is present, but analysis 4 wins all calibration designs. This is not a failure of the pairwise mechanism: v0.1 had deliberately reduced Delta to make the constant-effect candidate competitive when the slope heterogeneity was absent. Positive gamma also partly offsets analysis 4's negative treatment-effect distortion.

The only subsequent adjustment is to restore the original, documented pre-workaround treatment-effect values tau0=0.25 and Delta=0.25, preserving theta=0.50. The same gamma grid and seeds are run again. This removes the compensating v0.1 calibration adjustment now that the omitted substantive feature has been restored; no search over treatment-effect parameters is performed. The v0.1 parameter case remains in every generated calibration table, so its stronger performance for analysis 4 remains visible. Analysis 4 is never excluded or penalized by the selection rule.

Restoring tau0 and Delta does not change the pairwise 2-versus-5 loss comparison: both analyses allow the correct subgroup treatment coefficient, and analysis 5 includes D and D*G. At fixed theta, gamma, and realized design, their conditional biases and variances are identical across the two effect cases. The adjustment only changes the constant-effect candidate's treatment misspecification (and related observed estimates).

Final selection is gamma=0.10 at the restored effects: the grid confirms meaningful overall shares for both 2 and 5. Gamma=0.10 is the starting value motivated by the restored substantive feature, not a frequency target. It changes the X slope from 0.70 to 0.80. No additional gamma values or other parameter changes are considered. The narrower values below and above it establish whether the tradeoff is local and whether the change in bias behaves coherently. The pilot uses 2,000 seeds disjoint from calibration; parameters will not be retuned to pilot frequencies. Exact case results and the final suitability assessment are reported in the code-generated `pilot_results.md` and release verification below.

## Pairwise diagnostics

Each pilot design saves `Delta L25=L2-L5`, its absolute magnitude, the variance gain `V2-V5`, and the squared-bias penalty `b5(Z)^2-b2(Z)^2`. Pairwise numerical ties use the absolute tolerance 1e-10. We report both strict sign frequencies and tolerance-based disjoint 2-wins/5-wins/tie shares. Signed, absolute, and relative gap quantiles expose near-ties; the relative gap divides by the smaller loss. Fractions below absolute gaps of 1e-6, 1e-4, and 1e-3 are descriptive, with no declaration that an arbitrary threshold defines a meaningful research benefit.

Region comparison tables show counts, treatment balance, propensity spread/overlap, oracle variances and biases, leverage, and scaled condition numbers. Region means and medians are descriptive associations, not causal predictors. Pairwise comparisons require only valid fits for 2 and 5; overall winner frequencies still require the full five-analysis menu. The denominator and invalid comparisons are reported separately. The figure shows the loss-difference distribution and its association with subgroup sample size, without excluding the tails.

## Temporary trimming choice

Analysis 3 retains G=1 observations with **true known propensity** in the inclusive interval [0.10,0.90]. This is a transparent oracle overlap benchmark. It is not fitted logistic regression, a common empirical treated/control range intersection, or a claim about what researchers can observe in future regimes. Selection depends only on the design and is identical for realized outcomes and the oracle calculation. A fixed rule gives an unambiguous population counterpart. All removed treated/control counts are saved. Treatment constancy within G preserves the substantive target after trimming.

The cutoff is a documented temporary choice, not an optimized threshold. Discuss whether future work should use estimated propensity scores, how researchers would obtain them, and whether this candidate is scientifically useful given its dominance under both versions.

## Population counterpart and realization-specific residual

For each candidate, form population moments `A = E[WW' I(selected)]` and `q = E[W m I(selected)]`. Solve `A beta_a = q` and apply the same coefficient contrast used in the sample. Normalizing both moments by selection probability is unnecessary. The algorithm sums exactly over the two binary variables and uses base-R adaptive numerical integration over the normal covariate. The trimmed interval is found by inverting its known linear propensity logit; no discontinuity is left inside that integral.

Default bounds are [-10,10], with absolute and relative integration tolerances 1e-10. A second calculation expands bounds to [-12,12] and tightens tolerances to 1e-12. The saved refinement differences measure numerical stability; reported moment-integration error estimates are not formal bounds on coefficient error. The normal-tail contribution is negligible at this scale, and exact identities theta_2=theta_3=theta independently check the result; theta_5=theta is tested only at gamma=0. This approach avoids large artificial populations and population Monte Carlo error.

`theta_a` is the population projection/probability limit. Misspecified finite-sample OLS need not be centered exactly on that limit. Consequently `epsilon_a=T_a-theta_a` is a decomposition by construction, not an assertion that E[epsilon_a]=0 at n=1500. It contains both conditional noise and design-dependent displacement from the population counterpart:

```
epsilon_a = (E[T_a|Z] - theta_a) + (T_a - E[T_a|Z]).
```

The direction of analysis 1's bias also differs from the informal example. By the omitted-variable formula its bias is beta2 times the coefficient on D in the subgroup population regression of X² on D and X. The present assignment coefficients make that coefficient negative. The sign is a consequence of this calibration, not a change in the estimator or target. Analysis 4's population effect is tau0 plus Delta times the D coefficient from projecting D×G on D, X, X², and G, plus gamma times the D coefficient from projecting G×X on those regressors; it need not equal a simple subgroup-prevalence-weighted average.

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

## Before later information-regime experiments

The present task ends after calibration and reproduction checks. Review the final 2-versus-5 shares, sizes of the loss gaps, and the explicit restoration of the pre-workaround treatment effects. The known-propensity trimming convention and usefulness of a dominated analysis 3 remain open design questions. Future information regimes must not expose these oracle quantities as automatically researcher-observable. No information-regime experiments are implemented here.

## Release verification

The revised tests cover fixed-seed reproduction; the gamma term in both m and Y under common draws; a constant-within-group treatment contrast; independent lm contrasts; trimming boundaries and failure handling; exact fixed-design moments; population integration refinement; gamma=0 dominance; nonzero population and conditional bias for analysis 5; nonnegative loss and exact epsilon decomposition; pairwise alignment and tie boundaries; and fresh-process/order reproducibility. The fixed-design repeated-noise check remains 20,000 draws and must pass before each pilot.

The code-generated pilot report records failures, numerical refinement, and coherence discrepancies. The final pilot supports the intended subgroup-versus-pooling mechanism: both analyses win with gaps far above floating-point error in almost all designs, while analysis 3 remains dominated. Smaller subgroups and fewer subgroup controls characterize the pooling-win region; its variance gain is larger even though its conditional bias is slightly higher. Propensity summaries do not all point toward poorer overlap in that region, so sample size and overlap should not be conflated. These comparisons are descriptive.

The baseline is suitable for a reviewed next stage, with a quantitative limitation: always choosing analysis 2 is already close to the pairwise oracle in average loss. The generated report calculates the gain available from perfect pairwise design selection, so later experiments can be judged against an appropriately small benchmark. No additional calibration or regime experiments are performed to enlarge it. A fresh local Git clone restored the pinned renv environment and ran the complete test, 10-case calibration, and pilot pipeline under R 4.3.3. All four test suites passed; no pilot fit or SE failures occurred. The independent comparison found identical hashes for all 18 generated tables, raw draws, figures, configuration records, and the Markdown report. Session provenance is excluded from byte comparisons. This establishes reproduction in the same R/system environment, not cross-platform byte identity.

## Repository output policy

Generated simulation data, detailed tables, figures, and runtime records under `output/` are local outputs excluded from version control. The compact, code-generated `paper/pilot_results.md` is committed so results can be reviewed directly on GitHub. Run `Rscript scripts/run_all.R` from a restored clone to recreate them. No result transport, CSV splitting, or data reconstruction step is required. The report is generated by code; numerical pilot tables are not manually copied into tracked prose.
