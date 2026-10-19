# Simulation decisions and interpretation through v0.2

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


## v0.1.2: attainable value of design adaptation

Winner frequency is insufficient for calibration. A candidate can win many designs by very small margins, or win rarely with larger margins; neither frequency alone nor the unconditional absolute gap measures the gain over a good fixed choice. The relevant pairwise ceiling is `G25 = min(E[L2],E[L5]) - E[min(L2,L5)]`, and `g25 = G25/min(E[L2],E[L5])`. Equivalently, compute the expected nonnegative excess loss relative to the oracle for each fixed candidate and take the smaller expectation. Equivalently, `G25 = (E[|L2-L5|] - |E[L2-L5]|)/2`: variation in the loss difference must overcome the average advantage of a fixed choice. This comparison does not assume researchers observe the oracle or specify an information regime.

The requested nine-cell grid crosses n={1000,1500,2000} with p_G={0.08,0.12,0.16}. Each first-pass cell uses 3,000 independent designs. All outcome and assignment coefficients, gamma=0.10, tau0=0.25, Delta=0.25, sigma=1, the five regressions, trimming and the exact fixed-design loss remain unchanged. These sample sizes and prevalence values characterize nearby study structures, not new outcome mechanisms. No extra grid values are added. The canonical configuration remains n=1500 and p_G=0.12 regardless of the diagnostic results.

A separate diagnostic draw generator uses a fixed maximum length of 2000. For replicate r it draws independent arrays of U_G, X, U_D and u. Cells reuse those arrays; smaller n uses a prefix, G=1[U_G<p_G], and D=1[U_D<ps(X,G)]. This supplies common random numbers across both dimensions while preserving each cell's DGP. It does not modify the v0.1.1 sampler or reproduce that sampler's exact seed-to-design mapping. Differences between the old pilot and the diagnostic at the same parameters are Monte Carlo variation. Seed bases 500000 (grid) and 600000 (confirmation), plus r, are centralized in config. Final confirmation estimates are kept separate from the first pass.

Monte Carlo intervals account for the minimum over two expected losses. Write Ak=Lk-min(L2,L5). For G, simultaneous Bonferroni normal intervals cover E[A2] and E[A5]; taking minima of their lower/upper endpoints covers their minimum. For g, apply the same construction to E[Ak]/E[Lk] with paired delta-method influence values `(Ak-rk*Lk)/E[Lk]`. The intervals are pointwise asymptotic 95% Monte Carlo intervals, not simultaneous across the nine cells. They avoid relying on a differentiable minimum when the best fixed choice is uncertain. Branch-local standard errors are saved as supplementary summaries; the interval construction, not those local errors alone, handles a mean-loss crossover. Independent confirmation reduces the risk of interpreting a first-pass high estimate as a stable result. No extra resampling randomness is introduced.

All five analyses are still fitted. Pairwise gains require valid fits for both 2 and 5; the diagnostic stops rather than silently discard failed pairs. Failures for all five, overall winners, trimming dominance checks, detailed pairwise quantiles, and design-region summaries are generated locally. The report records count distributions and the bias/variance decomposition so any apparent gain can be evaluated alongside sample plausibility and numerical stability. CSVs and figures remain ignored; `paper/design_value_diagnostic.md` is generated and committed. The old pilot report is retained as the historical v0.1.1 result.

### Confirmation-cell rationale

The independent confirmation uses 10,000 designs each at the canonical (1500,0.12), (2000,0.08), and (1000,0.16). This choice was recorded after seeing the three p_G=0.08 first-pass cells, before seeing the completed p_G=0.16 grid. At p_G=0.08, increasing n moved the comparison from clear pooling dominance toward a crossover. The (2000,0.08) cell therefore characterizes that boundary. The (1000,0.16) cell is a scientific counterpart with the same expected subgroup size of 160 but a much smaller outside group. This distinguishes subgroup size from the information available for pooling; it is not a search for 50/50 shares or a mechanical selection of the highest grid estimate. These are candidate study structures for review, not changes to the canonical baseline. All confirmation cells use the same independent seed stream.

Zero or very rare observed reversals need separate care. At zero observed gain the local MCSE is also zero, which does not bound the magnitude of unobserved losses. The combined report therefore leaves the upper gain endpoint unavailable (NA), rather than presenting [0,0] as an informative interval. The CSV supplies the exact one-sided 95% binomial upper bound on reversal probability when none is observed; it is not a bound on expected gain. Normal Monte Carlo intervals are also fragile with only one or a few reversals. The larger confirmation runs focus on the baseline and plausible crossover cells, where both winner regions have substantial support, rather than trying to precisely estimate negligible tail gains everywhere.

### Assessment for review

The evaluated grid supports the mechanism but supplies only modest attainable pairwise gains. Neither candidate study structure turns it into a large-value adaptation problem: (2000,0.08) is near the 1% boundary, and (1000,0.16) is below 1%, compared with roughly 0.3% at the canonical baseline. Neither the 3–5% nor above-5% range is supported. The exact first-pass and independent confirmation estimates and intervals are generated in `design_value_diagnostic.md`; classification near 1% should be interpreted with those intervals rather than as a threshold crossing.

The two alternatives are scientifically interpretable sensitivity candidates. Both expect 160 subgroup observations, compared with 180 in the baseline, but differ markedly in outside-group information. Their increased switching reflects the same variance-versus-squared-bias mechanism, with gaps far above the 1e-10 tie tolerance. Their subgroup treated/control counts are ordinary for these sample sizes, and the diagnostics show no fit/rank failure problem. However, their best fixed losses are higher than the baseline's: greater attainable adaptation gain is not the same as better absolute estimation quality. No empirical study context has been supplied that would make either sample structure substantively preferable.

Recommendation: retain v0.1.1 as canonical and keep these two cells available for sensitivity discussion. The improvements in the pairwise ceiling are real but small in absolute and proportional terms, so they do not justify replacing the baseline simply to amplify a later result. This conclusion applies to the evaluated grid; it is not a global impossibility claim for every nearby calibration. No extra grid values, gamma changes, heterogeneity, analysis-menu edits, or loss redefinitions were pursued. Any later design-information experiment must be interpreted against this limited pairwise ceiling, and this release does not implement one.

The full-menu context also matters. In the sparsest cell (1000,0.08), analysis 4 wins most designs even though analysis 5 wins every pairwise 2-versus-5 comparison. Its bias can be outweighed by precision when the subgroup is small. Its occasional wins in other cells are preserved and reported, not tuned away. Analysis 3 remains weakly dominated under the maintained assumptions. All nine first-pass cells show zero fit, rank and SE failures. The confirmation failure counts and trimming checks are included in the generated report.

### v0.1.2 verification

All five test suites passed under R 4.3.3 with renv 1.0.3, including all existing tests and the 20,000-draw fixed-design oracle coherence test. New tests check hand-computable gains, dominance/equal-loss cases, nonnegativity, the oracle-minimum inequality, ratio symmetry, invalid-input rejection, nested common random numbers, unchanged canonical parameters, and reproduction of the supplied-data path in a separate R process. The diagnostic fitted all five analyses in 57,000 cell-replications (27,000 grid plus 30,000 confirmation); no fit, rank or SE failures or trimming-dominance violations occurred.

A fresh local Git clone restored the pinned environment, passed all five suites, and independently generated all nine grid cells and all three confirmation cells. `scripts/check_reproduction.R ... design-value` found identical hashes for all 41 diagnostic CSVs, configuration records, the figure and the generated Markdown report. Session information is recorded but excluded from byte comparisons. The figure was visually checked. This establishes reproduction in the same R/system environment; it is not a claim of cross-platform byte identity. A direct comparison with the v0.1.1 configuration also verified that `config$dgp` is exactly unchanged. No later simulation stage has been implemented.


## v0.2: researcher-observable design information

The canonical DGP, five-analysis menu and oracle machinery are retained exactly. The only additional case is the previously evaluated (2000,0.08) sensitivity, with no recalibration. The information boundary, OVB algebra, decision assumptions, selection diagnostics and unresolved issues are documented in `researcher_information_model.md`. Results are code-generated in `researcher_information_results.md`.

The new researcher interface accepts exactly Y,X,D,G. It has no configuration or oracle dependency and rejects simulation-only columns. Analysis 5's estimated omitted-slope consequence is gamma_hat*h5'(GX), which equals its treatment contrast minus the augmented model's contrast for any Y. The augmented model is a structural diagnostic and never a sixth menu choice. Analysis 2 is assigned no estimated restriction penalty as the flexible reference, a maintained decision assumption rather than knowledge of true zero bias.

The baseline precision estimate uses a common augmented-model residual variance times each estimator's squared weights. This transparent ordinary variance estimator avoids contaminating the shared noise estimate with the restricted model's omitted-slope residual. It maintains homoskedasticity and adequacy of the augmented diagnostic model; these assumptions are documented rather than supplied as truth. HC3 squared SEs are compared diagnostically with the same estimated bias penalty. No result-driven choice between precision estimators, penalty de-biasing, shrinkage, significance gate, sample splitting or cross-fitting is adopted.

The plug-in squared penalty contains estimation variance. Its conditional expected inflation over oracle squared bias is c5_GX^2 Var(gamma_hat|Z). The benchmark reports this quantity in evaluation only. The decision rule uses the unadjusted estimated square exactly as requested. Numerical ties within 1e-10 prefer analysis 2; invalid observed fits cause an explicit failure and stop the benchmark.

Outcome-noise selection is a separate issue from noisy design classification. Conditional on Z, gamma_hat noise is orthogonal to analysis 5's treatment noise but has covariance -c5_GX Var(gamma_hat|Z) with analysis 2's treatment noise. Consequently a rule based on estimated structural evidence can select treatment noise without directly rewarding a treatment estimate. Evaluation retains epsilon=T-theta_a, but also separates it into conditional-mean displacement from theta_a and pure outcome noise. Assigned oracle loss E[L_choice(Z)] is reported separately from actual same-sample MSE; these are generally unequal for a Y-dependent choice. A fresh independent error vector, with design and choice fixed, evaluates this distinction without changing the researcher rule or introducing an information regime.

The benchmark uses 10,000 common-random-number replications per case, with disjoint new seed bases 700000 and 800000, preserving all earlier seed streams. Detailed full-menu and paired outputs remain ignored under the existing policy; only generated Markdown summaries are committed. The code does not introduce synthetic datasets or any later fidelity, leakage, analyst-sensitivity, regime or validation stage.

### v0.2 assessment for review

The feasible same-sample rule worsens assigned design loss relative to always analysis 2 in both the canonical and sensitivity cases. The normalized fraction captured is negative; its large magnitude in the canonical case reflects the small attainable denominator and should be read alongside the raw losses. Estimates of gamma are approximately centered but noisy, and perceived loss differences correlate weakly with oracle differences. A noisy squared restriction estimate can have an inflated mean and still induce frequent pooling choices when its realized magnitude is small. The generated report gives the variance-inflation decomposition and confusion matrices rather than using agreement alone as a performance criterion.

The principal conceptual finding is a substantial downward pure treatment-noise selection shift under the feasible rule. It partially offsets positive systematic specification distortion, leaving a deceptively small net bias. This is present with both the baseline precision proxy and the HC3 diagnostic. The rule never explicitly rewards a favorable treatment result, yet noisy structural evidence and focal treatment errors share outcomes. Epsilon is therefore decomposed into design offset and pure noise before interpreting this mechanism. Actual MSE and assigned design loss remain separate outcomes, with Monte Carlo uncertainty reported for their difference.

The implementation is suitable for reviewing the observable decision layer, but these findings do not justify treating it as a successful design optimizer or as independent of realized treatment noise. Decisions about the plug-in squared penalty, reference-model assumption, precision estimator and same-sample coupling remain substantive choices for review. No correction, splitting, recalibration or later regime is introduced to improve the result. The canonical configuration remains exactly unchanged.

### v0.2 verification

All six test suites passed in the working checkout and a fresh local Git clone under R 4.3.3 and the restored renv 1.0.3 environment. The existing tests, including the repeated-noise oracle check, remain intact. The added suite checks the isolated observable closure and rejection of truth columns, deterministic augmented coefficients, the exact omitted-variable contrast identity, both feasible variance calculations against independent regressions, nonnegative perceived loss and deterministic ties, explicit rank failures, hand-computed unclipped gain fractions, score covariance identities, decomposition identities and fixed-seed/separate-process reproduction.

Both checkouts independently ran 10,000 replications in each case, retaining all five candidate fits. The final `scripts/check_reproduction.R ... researcher` comparison found identical hashes for all 21 generated CSVs, configuration records, the figure and the Markdown report. An initially incomplete working-checkout raw menu file was detected by this comparison; the canonical case was rerun from source and the entire comparison then passed. Final row-integrity checks confirmed 10,000 paired rows and 50,000 menu rows per case. No coefficient-fit, rank or feasible-precision failures occurred. The figure was visually checked. Session provenance is recorded but excluded from hashes; byte identity is established for this R/system environment, not every platform.

A direct comparison with the pre-v0.2 configuration verified that `config$dgp` is exactly unchanged. Earlier pilot and design-value reports remain historical records; only the new benchmark is regenerated for this version. All simulation data and figures remain ignored, with source and the code-generated Markdown report committed. Work stops at the observable-information layer for review.
