## Child growthcleanr permissiveness specs
***NOTE*** This is for future implementation.
Date created: 2026-03-25
Date last updated: 2026-04-17

The goal of this document is to define modifiable parameters that will adjust the permissiveness of the child growth cleaner algorithm. This is based on a framework first developed for adults. See the reference below. Overall, users will be able to select from one of several levels of permissiveness, which will adjust the settings for multiple parameters to predefined values. In addition to setting a level of permissiveness, users can override any individual parameter. Although the main permissiveness option will set adult and child permissiveness to the same level, users could also override that to set either adult or child to a different level.


Reference document: /Users/Shared/ClaudeCode/adult-growthcleanr/repo-adult-gc/gc-adult-permissiveness-spec-2026-03-24.md


## Overall Framework

| Level | Intended Use | Weight | Height | CF Handling* |
|-------|-------------|--------|--------|-----------|
| `loosest` | Retain all plausible values for children with a very wide range of medical conditions | Allows for extreme weight changes, such as with bariatric surgery or rapid fluid status change | Allows for a large degree height measurement error and unusual height patterns | Only excludes CFs with very unlikely degrees of z-score change from original (still more likely to exclude CFs than not) |
| `looser` | Retain most plausible values for children with a broad range of medical conditions | Allows for slightly less extreme weight changes than loosest limits | Allows for a large degree of height measurement error and many unusual height patterns | Excludes CFs with unlikely degrees of z-score change from original |
| `tighter` | Exclude some plausible values that likely reflect illness with large effect on body size | Excludes weights and weight changes that are likely due to conditions such as severe cachexia, very rapid weight change from bariatric surgery, rapid change in fluid status | Less permissive limits for height error and atypical height | Excludes CFs with relatively unlikely degrees of z-score change from original |
| `tightest` | Exclude extreme measurements even if plausible | Exclude larger weight changes even if plausible | Less permissive for height error or unusual height patterns | Excludes CFs with relatively unlikely degrees of z-score change from original |

*CF Note: All levels have some degree of special CF handling; users can override this and completely turn off special handling of CFs.

---

## Step 6: Carried Forward Detection and Rescue

CF rescue determines which carried-forward values are re-included based on how much the z-score changed between the original and repeated value. The `include.carryforward` parameter can disable CF detection entirely (all CFs kept as Include).

### CF rescue delta-Z thresholds

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `cf_deltaZ_single` (Code 4: single CF rescue threshold) | | 0.05 | | |
| `cf_deltaZ_single_wholehalfimp` (Code 5: single CF + imperial rounding, upper threshold) | | 0.1 | | |
| `cf_deltaZ_multi_teen` (Code 6: 2+ CFs, adolescent rescue threshold) | | 0.05 | | |
| `cf_deltaZ_multi_teen_wholehalfimp` (Code 7: 2+ CFs, adolescent + imperial, upper threshold) | | 0.1 | | |

### CF adolescent age thresholds

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `cf_teen_age_female` (years) | | 16 | | |
| `cf_teen_age_male` (years) | | 17 | | |

### CF global switch

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `include.carryforward` (TRUE = skip CF detection entirely) | | FALSE | | |

---

## Step 7: Biologically Implausible Values (BIV)

### Absolute BIV limits

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| **Weight** | | | | |
| `biv_wt_min_birth` (kg) | | 0.2 | | |
| `biv_wt_min_after_birth` (kg) | | 1 | | |
| `biv_wt_max_birth` (kg) | | 10.5 | | |
| `biv_wt_max_under2y` (kg) | | 35 | | |
| `biv_wt_max_all` (kg) | | 600 | | |
| **Height** | | | | |
| `biv_ht_min` (cm) | | 18 | | |
| `biv_ht_max` (cm) | | 244 | | |
| `biv_ht_max_birth` (cm) | | 65 | | |
| **Head Circumference** | | | | |
| `biv_hc_min` (cm) | | 13 | | |
| `biv_hc_max` (cm) | | 75 | | |
| `biv_hc_max_birth` (cm) | | 50 | | |

### Standardized BIV limits (z-score thresholds, unrecentered)

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `biv_wt_z_low_under1y` | | -25 | | |
| `biv_wt_z_low_1y_plus` | | -15 | | |
| `biv_wt_z_high` | | 22 | | |
| `biv_ht_z_low_under1y` | | -25 | | |
| `biv_ht_z_low_1y_plus` | | -15 | | |
| `biv_ht_z_high` | | 8 | | |
| `biv_hc_z_low` | | -15 | | |
| `biv_hc_z_high` | | 15 | | |

---

## Step 9: Evil Twins (OTL Detection)

Evil twins identifies adjacent extreme value pairs based on the difference in recentered z-scores between consecutive measurements.

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `otl_threshold` (z-score difference threshold for OTL) | | 5 | | |

---

## Step 11: Extreme EWMA (EWMA1)

### EWMA exponent (age-dependent)

Controls the weight decay for the exponentially weighted moving average. More negative = faster decay (less influence from distant measurements).

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `ewma.exp` (age ≤ 1 year) | | -1.5 | | |
| `ewma.exp` (age ≥ 3 years) | | -3.5 | | |
| `ewma.exp` (1–3 years) | | linear interpolation | | |
| `ewma_window` (max observations on each side) | | 15 | | |

### EWMA1 pre-filter

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `ewma1_prefilter_tbc` (min \|tbc.sd\| to process group) | | 3.5 | | |

### EWMA1 exclusion thresholds

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `ewma1_dewma_all` (dewma.all threshold) | | 3.5 | | |
| `ewma1_dewma_dir` (dewma.before / dewma.after threshold) | | 3 | | |
| `ewma1_tbc` (tbc.sd threshold) | | 3.5 | | |
| `ewma1_cdewma_all` (c.dewma.all threshold) | | 3.5 | | |

---

## Step 13: Final SDE Resolution

SDE resolution uses EWMA to choose among same-day duplicate values. The EWMA parameters above (`ewma.exp`, `ewma_window`) also apply here.

No additional permissiveness-specific parameters identified for this step beyond the shared EWMA parameters.

---

## Step 15: Moderate EWMA (EWMA2)

### EWMA2 pre-filter

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `ewma2_prefilter_tbc_range` (min tbc.sd range to process group) | | 1 | | |

### EWMA2 additional criteria (addcrit) thresholds

These control the directional confirmation required for EWMA2 exclusions. Both dewma.before and dewma.after must exceed this threshold, AND the tbc differences to adjacent measurements must also exceed this threshold.

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `ewma2_addcrit` (dewma.before/after and tbc_diff threshold) | | 1 | | |

### EWMA2 middle exclusion thresholds

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `ewma2_middle_dewma` (dewma.all threshold for middle values) | | 1 | | |
| `ewma2_middle_cdewma` (c.dewma.all threshold for middle values) | | 1 | | |

### EWMA2 birth weight thresholds

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `ewma2_birth_wt_dewma_near` (next meas < 1 year) | | 3 | | |
| `ewma2_birth_wt_dewma_far` (next meas ≥ 1 year) | | 4 | | |

### EWMA2 first measurement thresholds

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `ewma2_first_dewma_near` (next meas < 1 year) | | 2 | | |
| `ewma2_first_dewma_far` (next meas ≥ 1 year) | | 3 | | |

### EWMA2 last measurement thresholds

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `ewma2_last_dewma` (gap < 2 years, \|prior tbc\| < 2) | | 2 | | |
| `ewma2_last_tbc_prev_threshold` (threshold for "high" variant) | | 2 | | |
| `ewma2_last_ext_dewma` (gap ≥ 2 years, \|prior tbc\| < 2) | | 3 | | |
| `ewma2_last_ext_dop_diff` (tbc.sd − tbc_dop threshold for ext) | | 4 | | |
| `ewma2_last_ext_high_dewma_add` (additive: dewma > 1 + \|prior tbc\|) | | 1 | | |

---

## Step 16: Birth HT/HC EWMA2

Same structure as Step 15 EWMA2 but applied specifically to birth height and head circumference (agedays = 0). Uses same thresholds as EWMA2 birth weight rules above (3 for near, 4 for far).

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `ewma2_birth_hthc_dewma_near` | | 3 | | |
| `ewma2_birth_hthc_dewma_far` | | 4 | | |

---

## Step 17: Height/HC Velocity (mindiff/maxdiff)

### Height velocity scaling parameters

These control how Tanner/WHO reference velocities are scaled to determine the minimum allowed decrease (mindiff) and maximum allowed increase (maxdiff) in height between consecutive measurements.

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| **mindiff (height decrease tolerance)** | | | | |
| `mindiff_min_vel_mult` (multiplier on min velocity) | | 0.5 | | |
| `mindiff_constant` (subtracted constant, cm) | | 3 | | |
| `mindiff_short_gap_exp` (exponent for gap < 1 year) | | 2 | | |
| **maxdiff (height increase tolerance)** | | | | |
| `maxdiff_max_vel_mult` (multiplier on max velocity) | | 2 | | |
| `maxdiff_constant` (added constant, cm) | | 5.5 | | |
| `maxdiff_short_gap_exp` (exponent for gap < 1 year) | | 0.33 | | |
| `maxdiff_long_gap_exp` (exponent for gap > 1 year) | | 1.5 | | |

### Height velocity floors (max.ht.vel minimum values)

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `ht_vel_floor_base` (cm) | | 2.54 | | |
| `ht_vel_floor_2mo` (gap > 2 months, cm) | | 5.08 | | |
| `ht_vel_floor_6mo` (gap > 6 months, cm) | | 10.16 | | |
| `ht_vel_floor_1yr` (gap > 1 year, cm) | | 20.32 | | |

### WHO-based scaling for height (gap < 9 months)

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `who_ht_mindiff_mult` (multiplier on WHO mindiff) | | 0.5 | | |
| `who_ht_mindiff_constant` (subtracted, cm) | | 3 | | |
| `who_ht_maxdiff_mult` (multiplier on WHO maxdiff) | | 2 | | |
| `who_ht_maxdiff_constant` (added, cm) | | 3 | | |

### Birth adjustments for height

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `ht_birth_tolerance` (added to mindiff/maxdiff at birth, cm) | | 1.5 | | |

### Default fallback for mindiff (when no velocity reference available)

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `mindiff_default_ht` (cm) | | -3 | | |

### HC velocity parameters

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `who_hc_mindiff_mult` (multiplier on WHO mindiff) | | 0.5 | | |
| `who_hc_mindiff_constant` (subtracted, cm) | | 1.5 | | |
| `who_hc_maxdiff_mult` (multiplier on WHO maxdiff) | | 2 | | |
| `who_hc_maxdiff_constant` (added, cm) | | 1.5 | | |
| `hc_birth_tolerance` (added to mindiff/maxdiff at birth, cm) | | 0.5 | | |
| `mindiff_default_hc` (cm) | | -1.5 | | |

---

## Step 19: Pairs and Singles

### Two-measurement exclusion thresholds

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `pairs_tbc_diff_long` (z-score diff threshold, gap ≥ 1 year) | | 4 | | |
| `pairs_tbc_diff_short` (z-score diff threshold, gap < 1 year) | | 2.5 | | |

### Single-measurement exclusion thresholds

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `single_tbc_threshold` (\|tbc.sd\| threshold when DOP available) | | 3 | | |
| `single_comp_diff_threshold` (comp_diff threshold when DOP available) | | 5 | | |
| `single_tbc_threshold_nodop` (\|tbc.sd\| threshold when no DOP) | | 5 | | |

### Pairs/singles mode

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `lt3.exclude.mode` | | "default" | | |

---

## Step 21: Error Load

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `error.load.threshold` (error ratio to trigger exclusion of all remaining) | | 0.4 | | |
| `error.load.mincount` (minimum error count before evaluating) | | 2 | | |

---

## Other Configurable Parameters

These are declared parameters that may not directly affect exclusion rates but are listed for completeness.

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `recover.unit.error` (attempt unit error correction) | | FALSE | | |
| `biv.z.wt.low.young` (Step 7 lower WT cutoff, age<1y) | | -25 | | |
| `biv.z.wt.low.old` (Step 7 lower WT cutoff, age>=1y) | | -15 | | |
| `biv.z.wt.high` (Step 7 upper WT cutoff, all ages) | | 22 | | |
| `biv.z.ht.low.young` (Step 7 lower HT cutoff, age<1y) | | -25 | | |
| `biv.z.ht.low.old` (Step 7 lower HT cutoff, age>=1y) | | -15 | | |
| `biv.z.ht.high` (Step 7 upper HT cutoff, all ages) | | 8 | | |
| `biv.z.hc.low` (Step 7 lower HC cutoff, all ages) | | -15 | | |
| `biv.z.hc.high` (Step 7 upper HC cutoff, all ages) | | 15 | | |
| `height.tolerance.cm` (declared but not used in child algorithm) | | 2.5 | | |
| `adult_cutpoint` (age in years dividing child/adult) | | 20 | | |

---

## Summary Table

**Current default is looser**

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| **Step 6: Carried Forward** | | | | |
| CF single rescue delta-Z | | 0.05 | | |
| CF single+imperial rescue delta-Z | | 0.1 | | |
| CF multi-teen rescue delta-Z | | 0.05 | | |
| CF multi-teen+imperial rescue delta-Z | | 0.1 | | |

| **Step 7: BIV** | | | | |
| Absolute BIV limits | | current defaults | | |
| Standardized BIV z-score limits | | current defaults | | |

| **Step 9: Evil Twins** | | | | |
| OTL z-diff threshold | 6 | 5 | 4 | 3 |

| **Step 11: EWMA1 (Extreme)** | | | | |
| dewma.all / tbc.sd threshold | | 3.5 | | |
| dewma directional threshold | | 3 | | |

| **Step 15: EWMA2 (Moderate)** | | | | |
| addcrit threshold | | 1 | | |
| middle dewma | | 1 | | |
| birth WT dewma (near/far) | | 3 / 4 | | |
| first dewma (near/far) | | 2 / 3 | | |
| last dewma | | 2 | | |
| last-ext dewma | | 3 | | |
| last-ext DOP diff | | 4 | | |

| **Step 17: Height/HC Velocity** | | | | |
| mindiff constant (HT) in cm | -5 | -3  | - 3  | - 1  |
| maxdiff constant (HT) | | +5.5  | | |
| mindiff constant (HC) | | -1.5 | | |
| maxdiff constant (HC) | | +1.5  | | |
| HT birth tolerance | | ±1.5  | | |
| HC birth tolerance | | ±0.5 | | |

| **Step 19: Pairs/Singles** | | | | |
| Pair z-diff (≥1 year) | | 4 | | |
| Pair z-diff (<1 year) | | 2.5 | | |
| Single \|tbc\| (with DOP) | | 3 | | |
| Single comp_diff | | 5 | | |
| Single \|tbc\| (no DOP) | | 5 | | |

| **Step 21: Error Load** | | | | |
| error.load.threshold | 0.41 | 0.41 | 0.29 | 0.29 |
| error.load.mincount | | 2 | | |
NOTE: changed default; values above should be > (so for 0.41, 4/10 errors would not trigger)
---

## User Overrides

The `permissiveness` parameter sets defaults for all sub-parameters. Individual sub-parameters can be overridden by passing them explicitly to `cleangrowth()`. Explicit values always take precedence over the permissiveness preset.

---

## Notes

1. **`height.tolerance.cm`** is declared as a parameter in `cleangrowth()` but is not actually used in the child algorithm. The child algorithm uses Tanner/WHO velocity-based mindiff/maxdiff instead. This parameter may be relevant for the adult algorithm only.

2. **EWMA exponent** is age-dependent in the child algorithm (unlike the adult algorithm where it's a fixed parameter). The age-dependent interpolation between -1.5 (≤1 year) and -3.5 (≥3 years) could itself be parameterized.
Thinking through this: Exponent changes how much closer values dominate. Larger magnitude (lower) values result in closer values dominating more. 

3. **Steps not listed above** (Early 13: SDE-Identicals, Step 5: Temporary SDE) do not have obvious permissiveness-relevant thresholds — they use median-based selection rather than threshold comparisons. However, the EWMA parameters they rely on (for SDE resolution in Step 13) are covered above.
