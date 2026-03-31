---
name: movement-science-stats
description: |
  Statistical quality checklist for movement science and neuroscience research.
  Auto-triggers when analyzing data, interpreting results, running statistics,
  writing results sections, or reviewing analysis code. Based on Makin & Orban
  de Xivry (2019, eLife).
triggers:
  - statistical analysis
  - kinematics data
  - results interpretation
  - correlation
  - t-test
  - ANOVA
  - p-value
  - effect size
  - sample size
  - pre post
  - intervention
  - group comparison
  - writing results section
  - reviewing analysis
version: 1.0
source: "Makin TR & Orban de Xivry JJ (2019). Ten common statistical mistakes
         to watch out for when writing or reviewing a manuscript. eLife 8:e48175.
         https://pmc.ncbi.nlm.nih.gov/articles/PMC6785265/"
---

# Movement Science Statistics Checklist

Before finalising any analysis or results section, systematically check for
these 10 common mistakes.

---

## MISTAKE 1 — Missing or inadequate control condition

**The problem**: Pre-post changes can arise from habituation, learning, or
passage of time — not your intervention. A control group that lacks a sham
intervention, is underpowered, or has a different baseline is inadequate.

**Check**:
- Is there a control group/condition for every time-based comparison?
- Is the control group adequately powered (not just a token small group)?
- Were control and experimental groups sampled simultaneously with randomised allocation?
- Are experimenters blinded to expected outcomes?

**Fix**: If no proper control exists, frame conclusions as tentative.

---

## MISTAKE 2 — Comparing two effects without directly testing the difference

**The problem**: Concluding that effect A > effect B because A is significant
(p≤0.05) and B is not, WITHOUT directly comparing A and B. This is wrong even
if the p-values look very different. Extremely common in movement science and
neuroscience.

**Check**:
- Does any conclusion involve "significant in group A but not in group B"?
- Is there a direct statistical test of the *difference* between effects?

**Fix**: Use a direct comparison test (unpaired t-test, ANOVA interaction term,
Monte Carlo simulation for correlations). Never infer group differences from
separate significance tests. Point reviewers to Nieuwenhuis et al. (2011,
Nature Neuroscience) if needed.

---

## MISTAKE 3 — Inflating the unit of analysis

**The problem**: Using number of observations (trials, timepoints, limbs) as
degrees of freedom instead of number of independent subjects. Artificially
inflates df, lowers the significance threshold, makes spurious results look
real.

**Common in kinematics**: Using individual movement cycles, individual trials,
or bilateral limb data as independent observations.

**Check**:
- Are degrees of freedom based on number of subjects, not number of observations?
- Are within-subject repeated measures treated as independent?
- If bilateral data: are both limbs treated as independent subjects?

**Fix**: Use linear mixed models (LMM) with subject as random effect, or
aggregate to one value per subject before testing. Never use number of
trials/observations as df for between-subject inferences.

---

## MISTAKE 4 — Spurious correlations

**The problem**: Outliers or subgroup clustering can create or destroy
correlations that don't reflect true relationships. A single extreme data point
can dramatically inflate Pearson's R.

**Check**:
- Is a scatterplot shown for every reported correlation?
- Are outliers present, and if removed/retained, is the justification given?
- Are subgroups (e.g., patients vs. controls) pooled in a correlation without accounting for group?
- Are assumptions of linearity and approximate normality of residuals met?

**Fix**: Use robust correlation methods (bootstrapping, Winsorizing, skipped
correlations). If there are several groups or conditions, z-score per groupxcondition and then perform the analys.
Always plot data. Justify any data exclusions a priori.

---

## MISTAKE 5 — Small sample sizes

**The problem**: Small N only detects large effects. False positives with small
N have artificially large effect sizes — making spurious results look compelling.
Underpowered studies also miss real effects.

**Check**:
- Is an a priori power analysis reported, based on an independent effect size estimate?
- Is the sample size justified for the claimed effect size?
- Are extraordinary claims based on very small N flagged as preliminary?

**Fix**: Report power analysis. If N is inherently limited (clinical populations,
non-human primates), provide replications within and between cases. 

---

## MISTAKE 6 — Circular analysis ("double dipping")

**The problem**: Using the same data to both SELECT features of interest AND
then TEST hypotheses about those features. Recruits noise to inflate effects.
Common in neuroimaging, EEG, and kinematic ROI analyses.

**Examples in movement science**:
- Selecting a time window of interest based on where the effect looks biggest,
  then testing significance in that window
- Sub-grouping participants based on their response magnitude, then testing
  group differences in that same response
- Correlating post-intervention values with pre-post difference scores
  (both depend on the same post measure)

**Check**:
- Were selection criteria (time windows, ROIs, subgroups) defined before seeing the data?
- Is there any overlap between the data used for selection and the data used for testing?

**Fix**: Pre-register analysis criteria. Use independent data splits (separate
participants or trials) for selection and testing. Use cross-validation or
bootstrapping approaches.

---

## MISTAKE 7 — p-hacking / flexibility of analysis

**The problem**: Trying multiple analysis pipelines, outcome variables,
covariates, or exclusion criteria until p≤0.05 — without disclosing this.
Inflates false positive rate dramatically.

**Check**:
- Are all measured variables and analysis decisions reported?
- Is the analysis pipeline consistent with previous publications from this lab?
- Are pre-planned vs. exploratory analyses clearly distinguished?
- Were any outcome variables or preprocessing steps changed after data collection?

**Fix**: Distinguish confirmatory from exploratory analyses explicitly. Report
all tested variables. Pre-register where possible. Exploratory analyses are
valid but must be labelled as such and cannot support strong conclusions alone.

---

## MISTAKE 8 — Failing to correct for multiple comparisons

**The problem**: Running many tests increases false positive probability.
In a 2×3×3 design, the chance of at least one spurious significant effect is
~30% even with no real effects.

**Particularly relevant for**: Kinematic time-series analysis, EEG/EMG
multi-channel analysis, multiple joint/segment comparisons, multiple outcome
variables.

**Check**:
- How many independent statistical tests were run?
- Is a correction applied (Bonferroni, FDR, cluster-based permutation for time-series)?
- If no correction, is there explicit justification?

**Fix**: Apply appropriate corrections. For kinematic time-series data,
consider SPM (Statistical Parametric Mapping) or cluster-based permutation
tests rather than running t-tests at each timepoint independently.

---

## MISTAKE 9 — Over-interpreting non-significant results

**The problem**: A non-significant p-value does NOT mean no effect exists.
It could mean: (a) truly no effect, (b) insufficient power, or (c) poor
experimental design. Treating p>0.05 as "proof of no difference" is invalid.

**Check**:
- Does any conclusion state or imply that a non-significant result means
  "no effect" or "no difference"?
- Is effect size reported alongside p-values for all results?
- For null results claimed as meaningful: is Bayesian evidence or equivalence
  testing used?

**Fix**: Always report effect sizes AND confidence intervals alongside p-values.
For meaningful null results, use Bayesian statistics (e.g., Bayes Factor) or
equivalence tests to distinguish "no effect" from "underpowered". Never write
"there was no significant difference between groups" when you mean "p>0.05".

---

## MISTAKE 10 — Confusing correlation with causation

**The problem**: A significant correlation between two variables does not mean
one causes the other. Could reflect reverse causation, a third variable,
or coincidence.

**Common in movement science**: Correlating kinematic features with clinical
scores, neural measures with behavioural outcomes, training amount with
performance change.

**Check**:
- Is causal language used (affects, drives, leads to, causes) for correlational data?
- Is the manipulation truly experimental (randomised, controlled) or observational?

**Fix**: Use only associative language for correlational data ("was associated
with", "predicted", "co-varied with"). Reserve causal language for randomised
controlled experimental manipulations. Consider mediation/moderation analyses
or competing model testing to strengthen causal inference where appropriate.

---

## Quick Reference Checklist

Before submitting or finalising analysis, verify:

- [ ] Control condition adequate and properly powered
- [ ] Group/condition differences tested directly (not inferred from separate tests)
- [ ] Unit of analysis = number of subjects, not observations
- [ ] Scatterplots shown for all correlations; outliers addressed
- [ ] Sample size justified with power analysis
- [ ] Analysis criteria defined independently of results (no double dipping)
- [ ] All analysis decisions reported; pre-planned vs. exploratory distinguished
- [ ] Multiple comparisons corrected where needed
- [ ] Effect sizes AND confidence intervals reported alongside all p-values
- [ ] Causal language only used for true experimental manipulations

---

## Notes for Kinematics Data Specifically

- **Time-series data**: Never run independent t-tests at each timepoint.
  Use SPM1D, cluster-based permutation tests, or functional ANOVA.
- **Bilateral data**: Treat as within-subject repeated measure, not independent.
- **Trial-level data**: Aggregate to participant means before between-subject tests,
  or use LMM with participant as random effect.
- **Pre-post designs**: Always use ANCOVA (controlling for baseline) rather than
  simple post-only comparison; ANCOVA has more power in randomised studies.
- **MATLAB users**: Functions like `corrcoef`, `ttest2`, `anova1` do not protect
  against any of these mistakes — responsibility is entirely yours.
- **Python users**: `scipy.stats` and `pingouin` are your friends;
  `pingouin` has built-in effect sizes, power analysis, and robust correlations.

---

## Reference

Makin TR & Orban de Xivry JJ (2019). Ten common statistical mistakes to watch
out for when writing or reviewing a manuscript. *eLife* 8:e48175.
https://doi.org/10.7554/eLife.48175