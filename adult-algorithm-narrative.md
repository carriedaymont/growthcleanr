# Adult Growthcleanr Algorithm Narrative

## Key concepts

***EWMA:*** Exponentially Weighted Moving Average — a weighted average of a subject's other weight measurements, where measurements closer in time receive exponentially more weight. The exponent is e=−5, making the nearest neighbor almost entirely dominant; this is closer to a near-neighbor comparison than a traditional moving average. Three variants are computed for each observation: `ewma_all` (all other values), `ewma_before` (excluding the immediately preceding value), and `ewma_after` (excluding the immediately following value). The before/after variants are used as confirmation — a value is only excluded when the overall and both directional deviations agree. Edge values (first/last in subject) have only one neighbor, so before/after falls back to `ewma_all`. Window parameter `ewma_window` (default 15) limits how many observations on each side contribute.

***Ordered steps:*** Removing more extreme problems first, getting more refined at the end

***Iterations:*** Many exclusion steps use a dynamic loop: one value is excluded per round (the worst candidate by a defined scoring rule), then the EWMA and thresholds are recomputed from scratch on the remaining values, and the process repeats. The loop runs until no candidates remain or fewer than 3 values are left. This single-worst-out-per-round approach prevents a cluster of outliers from masking each other — removing the most deviant value first often exposes the next, which would have appeared plausible if evaluated simultaneously.

***Exclusion rescue:*** The adult algorithm does not rescue exclusions the way the child algorithm rescues carried-forward values. The one exception is the moderate EWMA step (11Wb): when 4 or more consecutive candidates are found in a single pass, they are not individually excluded — instead, the entire subject is flagged for error load escalation. This treats a long run of EWMA failures as evidence of overall data unreliability rather than a sequence of isolated errors.

***Flexibility:*** Defaults to `looser` (moderate permissiveness). Can be made more permissive (`loosest`) or more restrictive (`tighter`, `tightest`).

---

## Rounding Tolerance

All threshold comparisons in the algorithm apply a **0.12 unit tolerance** (cm for height, kg for weight). This accommodates rounding that may occur during data extraction, transformation, and loading — specifically, rounding to the nearest 0.1 cm or 0.1 kg, plus a small margin for floating-point precision.

The tolerance makes all threshold comparisons slightly more permissive: a value must exceed a limit by more than 0.12 to be excluded. For example, if the BIV minimum height is 50 cm, a recorded value of 49.9 cm is not excluded (it is within 0.12 of the threshold). When comparing two measured values (e.g., height difference in a band check), the tolerance is applied once (0.12), not doubled.

**Where tolerance is applied:**
- BIV limits (Step 1): value vs fixed min/max thresholds
- Weight cap EWMA checks (Step 4W): `dewma_all` vs 50 kg; `dewma_before`/`dewma_after` vs 45 kg (0.9 factor applied to 50 kg base); adjacency diffs vs 50 kg
- Evil twins caps (Step 9Wa): weight difference vs dynamic cap
- Extreme EWMA (Step 9Wb): dewma vs dynamic threshold
- Height distinct (Steps 10Ha/10Hb): band checks, gain/loss direction and rescue
- Weight 2D ordered and non-ordered (Steps 11Wa/11Wa2): weight difference vs wtallow
- Moderate EWMA (Step 11Wb): dewma vs wtallow (standard and alternate pathways)
- 1D limits (Step 13): value vs fixed min/max thresholds

**Where tolerance is NOT applied:**
- BMI thresholds (derived from ht+wt; 0.12 in either produces <0.1 BMI unit change)
- Perclimit ratio checks (conversion factors cancel in ratios)
- Error load ratio (Step 14; count-based, not measurement-based)
- Exact equality checks (RV identification, SDE identical detection — same conversion applied uniformly, so identical inputs produce identical outputs)

---

## Permissiveness Framework

The algorithm's aggressiveness is controlled by a `permissiveness` parameter with four levels. Each level sets defaults for multiple sub-parameters. Individual sub-parameters can be overridden; explicit values always take precedence over the preset.

| Level | Intended Use |
|-------|-------------|
| `loosest` | Retain all plausible values. Allows extreme weight changes (bariatric surgery, rapid fluid status change), significant height measurement error, height gain in young adults, and height loss. |
| `looser` (default) | Retain most plausible values. Slightly less extreme weight changes allowed. Height gain allowed but not height loss. |
| `tighter` | Exclude plausible values likely reflecting illness with large effect on body size. Less permissive height limits. Height gain allowed, not loss. Repeated values treated as less informative. |
| `tightest` | Exclude extreme measurements even if plausible. Tightest weight and height limits. No height gain or loss allowance. Repeated values treated as less informative. |

### Parameters controlled by permissiveness

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| Overall (BIV) HT limits | 50–244 cm | 120–230 cm | 142–213 cm | 147–208 cm |
| Overall (BIV) WT limits | 20–500 kg | 30–270 kg | 36–159 kg | 39–136 kg |
| Overall (BIV) BMI limits | 5–300 | 12–65 | 16–45 | 18–40 |
| Single (1D) limits | wider BMI-dependent split (see spec) | same as BIV | same as BIV | same as BIV |
| `wtallow_formula` | PW-H (piecewise) | PW-H (piecewise) | PW-L (piecewise-lower) | allofus15 |
| UW scaling | UW-based (see wtallow-formulas.md) | UW-based | UW-based | cap limited by PW-L |
| ET caps | wtallow cap + 20 | wtallow cap + 20 | wtallow cap + 20 | allofus15-cap-12m |
| `perclimit_low` (wt ≤45 kg) | 0.5 | 0.5 | 0.7 | 0.7 |
| `perclimit_mid` (45<wt≤80 kg) | 0.4 | 0.4 | 0.4 | 0.4 |
| `perclimit_high` (wt >80 kg) | 0 (disabled) | 0 (disabled) | 0.4 | 0.4 |
| `error_load_threshold` | 0.41 | 0.41 | 0.29 | 0.29 |
| `mod_ewma_f` | 0.75 | 0.75 | 0.60 | 0.60 |
| `ht_band` | 3" | 3" | 2" | 2" |
| `allow_ht_loss` | TRUE | FALSE | FALSE | FALSE |
| `allow_ht_gain` | TRUE | TRUE | TRUE | FALSE |
| `repval_handling` | independent | independent | linked | linked |

### Implementation

The `permissiveness` parameter defaults to `"looser"`. All sub-parameters default to `NULL` in the function signature. At the start of `cleanadult()`, `resolve_permissiveness()` fills NULLs from the preset table and passes through any non-NULL (user-specified) values. This means:
- `cleanadult(df)` uses all looser defaults
- `cleanadult(df, permissiveness = "tighter")` uses all tighter defaults
- `cleanadult(df, permissiveness = "tighter", allow_ht_loss = TRUE)` uses tighter defaults but overrides height loss to be allowed

Parameters NOT controlled by permissiveness: `scale_max_lbs` (default Inf), `ewma_window` (default 15), `quietly` (default FALSE).

---

## List of steps

***For the most part, HT and WT are evaluated separately, although there are some shared steps. See separate tables for HT and WT steps below.***

***Step numbering is not sequential. Steps 5–8 and 12W are not present in the algorithm; numbering is preserved for consistency with prior documentation.***

Height steps (include steps shared by HT and WT):

| Step # | Step title | Brief description |
|--------|-----------|-------------------|
| 1H | BIV | Exclude biologically implausible heights |
| 3H | Temp SDE | Temporarily flag same-day height duplicates |
| 9H | Final SDE Height | Resolve same-day height duplicates: exclude identical, choose keeper by category |
| 10H | Height Distinct | Evaluate 2D pairs (band, loss/gain rescue, frequency) and 3+D windows (w2, loss/gain groups) |
| 11H | Mean Height | Calculate mean height per subject (group-aware for loss/gain groups) |

Weight steps (include steps shared by HT and WT):

| Step # | Step title | Brief description |
|--------|-----------|-------------------|
| 1W | BIV | Exclude biologically implausible weights |
| 2W | Repeated Values | Identify groups of identical weight values; mark first vs. subsequent |
| 3W | Temp SDE | Temporarily flag same-day weight duplicates |
| 4W | Scale Max | Exclude weights at physical scale maximum unless EWMA confirms |
| 9Wa | Evil Twins | Exclude one member of adjacent pairs with implausibly large weight differences |
| 9Wb | Extreme EWMA | Exclude weight outliers with large EWMA deviations (interval-specific caps) |
| 10W | Final SDE Weight | Resolve same-day weight duplicates: exclude identical, choose keeper by category |
| 11Wa | 2D Ordered Weight | Exclude 2D ordered weight pairs outside wtallow or perclimit |
| 11Wa2 | 2D Non-Ordered Weight | Exclude 2D interleaved weight pairs using wtallow + dominance rules |
| 11Wb | Moderate EWMA | Exclude weight outliers using 7-step EWMA flow with trajectory rescue and error load |
| 13 | Single Distinct (1D) | Exclude 1D height/weight values outside BMI-based or no-BMI limits |
| 14 | Error Load | Exclude all remaining values when error ratio exceeds threshold |

---

## Step 1: Biologically Implausible Values (BIV)

***Exclude values outside absolute biological limits.***

| | |
|---|---|
| Scope | Height and Weight (individually), then same-day ht+wt pairs (BMI check) |
| Prior Step | None (first step) |
| Next Step | 2W: Repeated Values (Weight); 3H: Temp SDE (Height) |
| Distinct/RVs | N/A — runs before RV identification |
| Exclusion Code | `Exclude-A-HT-BIV`, `Exclude-A-WT-BIV` |
| Configurable parameter names | `overall_ht_min`, `overall_ht_max`, `overall_wt_min`, `overall_wt_max`, `overall_bmi_min`, `overall_bmi_max` |

***Overview:***

The first step removes values that are biologically impossible for an adult human. These absolute outer limits are deliberately very wide — they are meant to catch only values that could not possibly be correct (e.g., data entry errors, unit errors that produce extreme values). More refined exclusion of implausible-but-not-impossible values happens in later steps.

***Key terms and variable names:***

- `biv_df` — data frame holding the BIV thresholds, indexed by "height" and "weight" with columns "low" and "high"
- `remove_biv()` — support function that returns a logical vector identifying values outside BIV limits

***Configurable parameter defaults and options:***

| Parameter | Default (loosest) | Unit | Notes |
|-----------|---------|------|-------|
| `overall_ht_min` | 50 | cm | |
| `overall_ht_max` | 244 | cm | |
| `overall_wt_min` | 20 | kg | Not pounds — specified because scale_max_lbs uses pounds |
| `overall_wt_max` | 500 | kg | Not pounds — specified because scale_max_lbs uses pounds |
| `overall_bmi_min` | 5 | kg/m² | Only meaningfully different from `single_bmi_min` at loosest (5 vs 10); same at other levels |
| `overall_bmi_max` | 300 | kg/m² | Only meaningfully different from `single_bmi_max` at loosest (300 vs 65); same at other levels |

These thresholds use strict inequalities with 0.12 rounding tolerance (see Rounding Tolerance section): a value must be more than 0.12 below the minimum or above the maximum to be excluded. The `permissiveness` parameter provides preset less-permissive defaults that tighten these (and other) limits; explicit parameter values override permissiveness presets.

***Logic and implementation:***

1. For each subject, evaluate all height values against the height BIV limits
2. Exclude any height outside [`overall_ht_min`, `overall_ht_max`] (with 0.12 tolerance)
3. For each subject, evaluate all weight values against the weight BIV limits
4. Exclude any weight outside [`overall_wt_min`, `overall_wt_max`] (with 0.12 tolerance)
5. **BMI BIV check (Step 1B):** Runs only when `overall_bmi_min < single_bmi_min` or `overall_bmi_max > single_bmi_max` — i.e., only when the BIV BMI limits are wider than the Step 13 single-value limits. In practice this means loosest only (overall BMI [5,300] vs single BMI [10,250]). At other levels the two sets of limits are equal, so skipping the check avoids applying thresholds earlier than the algorithm intends.
6. For each ageday where both a surviving height and a surviving weight exist, pair the first height and first weight (by id order) and compute BMI = weight_kg / (height_cm / 100)²
7. Exclude that ht and wt pair if BMI < `overall_bmi_min` or BMI > `overall_bmi_max`. No rounding tolerance is added to the BMI comparison: tolerance is already embedded in the measurement-level checks above, and BMI is a derived value
8. Excluded values are removed from the working dataframes so they do not participate in subsequent steps
9. There is no need to check for same-day duplicates before this step because all values outside these limits are excluded regardless

***Rationale for selected decisions:***

- The weight maximum (500 kg) is in kilograms, not pounds. This is noted explicitly because the algorithm also handles pound-to-kg conversion, and confusion between units is a common source of data error.
- These limits are intentionally very broad. For example, a height of 50 cm (~20 inches) would be implausible for an adult, but unusually small adults and data that might represent seated heights or other measurement contexts are handled by later, more nuanced steps.
- BIV values count toward the error load denominator in Step 14 (Too Many Errors). A subject with many BIV values plus a few other exclusions could trigger error load on their remaining included values.

---

## Step 2W: Repeated Value Markers (Weight only)

***Identify groups of identical weight values and distinguish first occurrence from subsequent repeats.***

| | |
|---|---|
| Scope | Weight only |
| Prior Step | 1W: BIV |
| Next Step | 3W: Temp SDE |
| Distinct/RVs | Creates RV markers (this is where RVs are defined) |
| Exclusion Code | None — this step marks but does not exclude |
| Configurable parameter names | `repval_handling` (affects later EWMA steps, not this step) |

***Overview:***

This step identifies repeated values (RVs) — weight measurements that appear more than once for a subject. The first occurrence is distinguished from subsequent identical values. RVs are not excluded at this point; instead they are flagged so that later steps (particularly the EWMA steps) can handle them differently depending on the `repval_handling` mode.

Height does not use RV tracking. While height does identify same-day duplicates (Step 3H), it does not distinguish first vs. subsequent occurrences of the same height value across different days.

***Key terms and variable names:***

- **Repeated value (RV):** A weight measurement value that appears more than once for a subject (exact match; 78.1 and 78.101 are not identical)
- **First RV:** The earliest occurrence (by age, with internal_id as tiebreaker) of a repeated value
- `is_rv` — logical flag: TRUE for non-first occurrences of a repeated value
- `is_first_rv` — logical flag: TRUE for the first occurrence of a repeated value
- Unique values have both `is_rv = FALSE` and `is_first_rv = FALSE`
- `identify_rv()` — support function that sets RV flags based on current included values
- `redo_identify_rv()` — support function that re-identifies RVs specifically after temp SDE marks a first_rv as extraneous

***Configurable parameter defaults and options:***

| Parameter | Default | Options | Notes |
|-----------|---------|---------|-------|
| `repval_handling` | `"independent"` | `"independent"`, `"linked"` | Does not affect this step directly; controls how RVs are handled in EWMA steps (9Wb, 11Wb) |

- **independent:** RVs are treated as full participants in a single EWMA pass. Best when RVs are likely real measurements (or approximations).
- **linked:** Two-pass structure (firstRV/allRV) with RV propagation. Best when RVs are likely carried-forward values, not independent measurements.

***Logic and implementation:***

1. For each included weight value (in `meas_m`), count how many times it appears for the subject
2. Values appearing only once: `is_first_rv = FALSE`, `is_rv = FALSE`
3. Values appearing more than once:
   - The first occurrence (earliest age, internal_id tiebreaker): `is_first_rv = TRUE`
   - All subsequent occurrences: `is_rv = TRUE`
4. RV identification depends on the dataframe being pre-sorted by age and internal_id
5. After later exclusion steps (weight cap, evil twins), RVs are re-identified on the remaining values:
   - `identify_rv()` is called fresh on the remaining rows
   - `temp_sde()` is re-run
   - `redo_identify_rv()` cleans up if SDE marked a first_rv as extraneous

***Rationale for selected decisions:***

- RVs are weight-only because repeated identical heights across visits are common and expected (adults don't grow), while repeated identical weights are more likely to represent carried-forward or default values.
- The first occurrence is kept as the "real" measurement because it is the earliest recorded instance; subsequent identical values may be copies.
- Re-identification after exclusions is necessary because if the first occurrence of a value is excluded, the next occurrence should become the new first_rv.

---

## Step 3: Temp Same-Day Extraneous (SDE)

***Temporarily flag same-day values that deviate most from patient median.***

| | |
|---|---|
| Scope | Height (3H) and Weight (3W), evaluated separately |
| Prior Step | 1H: BIV (Height); 2W: Repeated Values (Weight) |
| Next Step | 4W: Weight Cap (Weight); 9H: Final SDE Height Resolution (Height) |
| Distinct/RVs | WT median uses non-RV values only; HT median uses all values |
| Exclusion Code | None — values are flagged as `extraneous`, not excluded yet |
| Configurable parameter names | None |

***Overview:***

When a subject has multiple measurements on the same day, this step temporarily identifies which value to keep and which to flag as extraneous. For each same-day group, the value closest to the patient's overall median survives; the others are flagged `extraneous = TRUE`. This is a temporary designation — flagged values remain in the working dataframe but are skipped by most subsequent steps. Final SDE resolution (where extraneous values are actually excluded or rescued) happens later (Step 9H for height, Step 10W for weight).

***Key terms and variable names:***

- **Same-day extraneous (SDE):** A measurement flagged because another measurement of the same type on the same day is a better fit for the patient
- **Extraneous:** The temporary flag; `extraneous = TRUE` means this value is tentatively set aside
- `temp_sde()` — support function that computes the median, diffs, and sets `extraneous` flags
- `redo_identify_rv()` — called after 3W to update RV identification if a first_rv was flagged extraneous

***Configurable parameter defaults and options:***

None. This step has no configurable parameters.

***Logic and implementation:***

1. Identify days with multiple measurements for the subject (same `age_days`)
2. Calculate the patient's overall median:
   - **Weight:** Median of all included non-RV weight values
   - **Height:** Median of all included height values
3. For each value on a same-day group, compute `diff = abs(measurement - median)`
4. Within each same-day group, keep the value with the smallest diff (highest internal_id tiebreaker); mark all others as `extraneous = TRUE`
5. **Weight only:** After flagging, call `redo_identify_rv()` to update RV identification in case a first_rv was flagged extraneous
6. Store extraneous flags in the subject-level tracking vector (`h_extraneous` or `w_extraneous`)

***Rationale for selected decisions:***

- The median is calculated from all included values (not just non-same-day values) because the goal is to identify which same-day value best represents the patient's overall pattern.
- Weight uses non-RV values for the median because repeated values may be carried-forward and could bias the median. Height uses all values because the RV concept does not apply to height.
- This is a temporary step because the "correct" same-day value may change after later exclusion steps alter the patient's value distribution. Final resolution happens in Steps 9H and 10W.

---

## Step 4W: Scale Max (Weight Cap)

***Exclude weights at a physical scale maximum unless EWMA and adjacency checks confirm the value fits the patient's pattern.***

| | |
|---|---|
| Scope | Weight only |
| Prior Step | 3W: Temp SDE |
| Next Step | 9Wa: Evil Twins |
| Distinct/RVs | EWMA and adjacency use firstRV (non-RV) values; RVs of excluded cap values are also excluded |
| Exclusion Codes | `Exclude-A-WT-Scale-Max`, `Exclude-A-WT-Scale-Max-Identical`, `Exclude-A-WT-Scale-Max-RV-Propagated` |
| Configurable parameter names | `scale_max_lbs` |

***Overview:***

Many clinical scales have a physical upper limit (commonly 400 lbs). When a weight exactly equals this limit, it is likely a truncated reading rather than a true weight. This step identifies weights near the scale maximum and excludes them unless they are consistent with the patient's other weights (as determined by EWMA and adjacency checks). This step is optional — it only runs when `scale_max_lbs` is set to a finite value.

There are two cases:
1. **All distinct weights at cap:** If a subject's only distinct non-RV weight is the cap value, all of the subject's weights are excluded (code: `Scale-Max-Identical`).
2. **Mixed cap and non-cap values:** Individual cap values are evaluated using EWMA delta and adjacency difference criteria. Values must be >50 kg away from expected by both measures to be excluded (code: `Scale-Max`). RVs of excluded values are also excluded (code: `Scale-Max-RV-propagated` in linked mode).

***Key terms and variable names:***

- **Scale max / weight cap:** The physical upper limit of a clinical scale, specified in pounds
- **EWMA (Exponentially Weighted Moving Average):** A weighted average of a subject's other measurements, where closer-in-time measurements get exponentially more weight
- **dewma:** Delta from EWMA = measurement - EWMA (positive means above expected, negative means below)
- **Before/after EWMA:** EWMA recalculated excluding the immediately adjacent measurement (to avoid the value being evaluated from dominating its own reference)
- **Adjacency diff:** Difference between a value and its chronologically neighboring values
- `adult_ewma_cache_init()` — support function computing EWMA with before/after variants using matrix operations
- `check_between()` — support function checking if values fall within a range

***Configurable parameter defaults and options:***

| Parameter | Default | Unit | Notes |
|-----------|---------|------|-------|
| `scale_max_lbs` | `Inf` | pounds | `Inf` disables this step entirely. Common value: 400 (lbs). |

The cap detection range is: `round(scale_max_lbs / 2.2046226, 1) ± 0.1` kg, using inclusive bounds.

***Logic and implementation:***

1. Skip if `scale_max_lbs == Inf` or fewer than 2 included non-extraneous values
2. Convert `scale_max_lbs` to kg; define detection range as `round(kg, 0.1) ± 0.1`
3. Identify non-RV values at the cap and count distinct non-RV weights
4. **If only 1 distinct non-RV weight and it's the cap:** Exclude all of the subject's weights with `Scale-Max-Identical`
5. **If mixed cap and non-cap, and >1 distinct non-RV weight:**
   a. Calculate EWMA on non-RV values (exponent = -5), including `ewma_before` (EWMA excluding the immediately preceding observation) and `ewma_after` (EWMA excluding the immediately following observation). At edges, before/after EWMA equals the full EWMA.
   b. Calculate adjacency diffs (difference from previous and next non-RV values)
   c. For positive outliers (dewma50p): `dewma_all > 50 + 0.12` AND `dewma_before > 45 + 0.12` AND `dewma_after > 45 + 0.12` AND `dewma_all` is not NA
   d. For negative outliers (dewma50m): `dewma_all < −(50 + 0.12)` AND `dewma_before < −(45 + 0.12)` AND `dewma_after < −(45 + 0.12)` AND `dewma_all` is not NA
   e. For positive adjacency (d50p): `d_prev > 50 + 0.12` AND `d_next > 50 + 0.12`. Missing neighbors (first/last in series): treated as +Inf (confirms exclusion).
   f. For negative adjacency (d50m): `d_prev < −(50 + 0.12)` AND `d_next < −(50 + 0.12)`. Missing neighbors: treated as −Inf (confirms exclusion).
   g. Exclude cap values where (dewma50p AND d50p) OR (dewma50m AND d50m)
6. Propagate exclusion to RVs of excluded cap values
7. After exclusions: re-run `identify_rv()`, `temp_sde()`, `redo_identify_rv()` on remaining values

***Rationale for selected decisions:***

- The 50 kg threshold with 0.9 factor (45 kg) on before/after EWMA provides a conservative exclusion criterion — values near the cap are excluded only when clearly inconsistent with the patient's pattern. All comparisons include the standard 0.12 kg rounding tolerance.
- Missing adjacency neighbors (first or last value in series) are treated as confirming exclusion for both positive and negative directions. A cap value at the edge of a patient's data with no neighbor on one side should still be excludable if the other neighbor confirms.
- The step does not exclude heights, even in the "all weights at cap" case. Height exclusion is handled by other steps.
- The 0.9 factor on before/after EWMA allows slight variation when excluding a value changes the EWMA reference.

---

## Step 9Wa: Evil Twins

***Exclude one member of adjacent weight pairs with implausibly large differences.***

| | |
|---|---|
| Scope | Weight only |
| Prior Step | 4W: Scale Max |
| Next Step | 9H: Final SDE Height Resolution; 9Wb: Extreme EWMA |
| Distinct/RVs | Both Inc and RV values participate in OOB detection, median, and pairs guard |
| Exclusion Code | `Exclude-A-Evil-Twins` |
| Configurable parameter names | `wtallow_formula` (shared with EWMA steps) |

***Overview:***

Evil twins are adjacent weight values separated by an implausibly large difference — so large that they appear to belong to different people. This step identifies such pairs and excludes one member (the one farthest from the patient's median, with a plausibility guardrail for extreme values). It runs before the Extreme EWMA step (9Wb) to prevent mutually-anchoring outlier pairs from distorting EWMA calculations.

The step iterates, each time re-evaluating the remaining values, until no OOB pairs remain or fewer than 3 values are left.

***Key terms and variable names:***

- **Evil twin pair:** Two adjacent weight values (sorted by age) whose absolute difference exceeds an interval-specific threshold
- **OOB (out of bounds):** Flag marking both members of an evil twin pair
- **Upper weight (UW):** The larger of the two weights in a pair, used for UW-based ET cap scaling. See `wtallow-formulas.md` for full specification.
- **Plausibility guardrail:** Values <38 kg or >180 kg that are in OOB pairs get maximum deviation score (99999), ensuring they are excluded first
- **Pairs guard:** Requires ≥3 included values to exclude; with only 2 values, it's ambiguous which is wrong
- `evil_twins()` — support function implementing the iterative detection and exclusion logic. Signature: `evil_twins(w_subj_df, wtallow_formula)`
- `compute_et_limit()` — computes the ET cap for a given interval and formula. Signature: `compute_et_limit(interval_months, formula, uw)`. The ET cap is derived from the wtallow formula's 6m and 12m caps + 20.

***Configurable parameter defaults and options:***

The ET cap is derived from the wtallow formula and UW (upper weight of the pair). For PW-H and PW-L formulas, ET cap = wtallow cap at 6m or 12m + 20. For allofus15, ET cap = allofus15-cap-12m (flat). For custom CSV, ET caps are fixed at 70/100.

| Formula | ET cap ≤6m (UW=120) | ET cap >6m (UW=120) | UW scaling |
|---------|---------------------|---------------------|------------|
| PW-H (loosest/looser) | 70 | 100 | Yes — caps scale with UW; see `wtallow-formulas.md` |
| PW-L (tighter) | 53.33 | 73.33 | Yes for UW<120 (scales down); No for UW>120 |
| allofus15 (tightest) | allofus15-cap-12m | allofus15-cap-12m | Cap limited by PW-L effective cap at 12m |
| Custom CSV | 70 | 100 | No |

All thresholds include the standard 0.12 kg rounding tolerance (see Rounding Tolerance section). See `wtallow-formulas.md` for full formula details including UW adjustment for high and low UW.

***Logic and implementation:***

1. **Pre-check:** Skip if <3 non-extraneous values, or if max-min range ≤ minimum possible ET cap + 0.12 (computed via `compute_et_limit(1, formula, uw = min(weights))`, no pair could exceed any interval-specific cap)
2. **Repeat until no OOB pairs found or <3 values remain:**
   a. Remove previously excluded values; break if <3 remain
   b. Sort by age (internal_id tiebreaker)
   c. Compute adjacent weight differences and age intervals (in months, using 30.4375 days/month)
   d. Compute dynamic ET cap per pair: `compute_et_limit(interval, formula, uw = uw_pair)` + 0.12 rounding tolerance. `uw_pair` is the larger of the two weights in the pair. ET caps scale with UW (see wtallow-formulas.md).
   e. Flag OOB pairs: both members of any pair where diff > cap
   f. If no OOB pairs, break (done)
   g. Compute subject median from all remaining values (Inc + RV)
   h. Compute absolute deviation from median for each value
   i. **Plausibility guardrail:** Set deviation to 99999 for OOB values <38 kg or >180 kg
   j. **Pairs guard:** Break if <3 values remain
   k. Among OOB values, sort by: most deviant from median → highest weight → lowest id
   l. Exclude the top-ranked OOB value
3. After all exclusions: re-run `identify_rv()`, `temp_sde()`, `redo_identify_rv()` on remaining values

***Rationale for selected decisions:***

- Two-tier interval thresholds (≤6m / >6m) reflect that larger weight changes are more plausible over longer time periods. UW-based scaling accommodates larger absolute changes for heavier patients and applies tighter limits for lighter patients (see wtallow-formulas.md).
- The plausibility guardrail ensures that biologically extreme values (<38 or >180 kg) are excluded first when they are part of an OOB pair, regardless of which value is closer to the median.
- The pairs guard prevents excluding values when there are only 2 included weights — with only a pair, there is no reference to determine which is correct.
- Both Inc and RV values participate fully in all aspects (OOB detection, median, pairs guard) because in the algorithm's default independent mode, RVs are treated as real measurements.
- The step runs before Extreme EWMA because mutually-anchoring outlier pairs (e.g., two very high values) would pull the EWMA toward them, making them appear more plausible than they are.
- Dynamic iteration (rather than a fixed number of rounds) ensures all OOB pairs are resolved before passing to EWMA steps, regardless of how many cascading exclusions are needed.

---

## Step 9H: Final SDE Height Resolution

***Resolve same-day height duplicates: exclude identical values, then choose keepers using category-based median selection.***

| | |
|---|---|
| Scope | Height only |
| Prior Step | 3H: Temp SDE (temporary flagging) |
| Next Step | 10H: Height Distinct Pairs |
| Distinct/RVs | Height does not use RV tracking |
| Exclusion Codes | `Exclude-A-HT-Identical`, `Exclude-A-HT-Extraneous` |
| Configurable parameter names | None |

***Overview:***

This step does the final resolution of same-day height duplicates that were temporarily flagged in Step 3H. It has two parts:

- **Part A (Identical):** When all same-day values are identical, exclude all but one (keep lowest id). If removing identical values eliminates all duplicates on a day, that day is no longer treated as SDE.
- **Part B (Non-identical):** For remaining days with multiple different values, categorize the subject and choose which value to keep based on category-specific median comparisons.

The keeper retains its original measurement value — no averaging or mean replacement is performed.

***Key terms and variable names:***

- **SDE day:** A day with more than one remaining included measurement after identical removal
- **daymed:** Median of all values within a single day
- **nonsdemed:** Subject-level median of values on non-SDE days only
- **medmed:** Median of day-medians (one median per day, balancing weighting across days with different numbers of values)
- **sderatio:** Ratio of SDE days to total distinct days

***Configurable parameter defaults and options:***

None. This step has no configurable parameters.

***Logic and implementation:***

**Part A — Identical values:**
1. All values (including those temporarily marked extraneous by Step 3H) remain in `h_subj_df`; extraneous flags identify which days have duplicates
2. On each day with extraneous values, identify measurements that appear more than once (exact match on meas_m)
3. Keep the first occurrence (lowest id); exclude the rest with `Exclude-A-HT-Identical`
4. Re-run `temp_sde()` to update extraneous flags on remaining values

**Part B — Non-identical SDEs:**
1. Identify remaining days with SDEs
2. Calculate three reference medians:
   - `daymed`: median within each day
   - `nonsdemed`: median of values on non-SDE days (NA if no non-SDE days exist)
   - `medmed`: median of per-day medians
3. Categorize subject:
   - **Category 1 (one day):** `daystot == 1`
   - **Category 2 (high non-SDE):** `daystot >= 4` AND `sderatio < 0.5`
   - **Category 3 (low non-SDE):** `daystot == 2 or 3` OR `sderatio >= 0.5`
4. For each SDE day, choose keeper by sorting:
   - **Category 1:** closest to day-median, tiebreaker highest id
   - **Category 2:** closest to non-SDE median, then day-median, then highest id
   - **Category 3:** closest to median-of-medians, then non-SDE median, then day-median, then highest id
5. Keep the first value after sorting; exclude the rest with `Exclude-A-HT-Extraneous`

***Rationale for selected decisions:***

- **Tiebreaker — lowest id for identical, highest id for non-identical:** For identical values there is no meaningful difference, so keep the first-entered. For non-identical values, the later measurement may represent a more careful re-measurement.
- **Three categories:** The choice of reference median depends on how much non-SDE data exists. With many non-SDE days (Category 2), the non-SDE median is the best reference. With few non-SDE days (Category 3), the median-of-medians balances across days. With only one day (Category 1), only the day-median is available.
- **No averaging:** The keeper retains its original value. Mean height calculation is a separate step that happens later (Step 11H).

---

## Step 9Wb: Extreme EWMA

***Exclude weight outliers whose EWMA deviation exceeds interval-specific caps.***

| | |
|---|---|
| Scope | Weight only |
| Prior Step | 9Wa: Evil Twins |
| Next Step | 10W: Final SDE Weight Resolution |
| Distinct/RVs | Independent: all values participate. Linked: firstRV pass (non-RV only), then propagate, then allRV pass. |
| Exclusion Codes | `Exclude-A-WT-Traj-Ext-N` (independent), `Exclude-A-WT-Traj-Extreme-firstRV-N` (linked), `Exclude-A-WT-Traj-Extreme-allRV-N` (linked) |
| Configurable parameter names | `wtallow_formula`, `repval_handling` |

***Overview:***

This step uses Exponentially Weighted Moving Averages to identify weight values that deviate far from the patient's expected trajectory. The EWMA gives exponentially more weight to measurements closest in time — with the exponent e=-5, the nearest neighbor almost completely dominates, making this a near-neighbor comparison rather than a traditional moving average.

A value is excluded when its EWMA deviation exceeds an interval-specific ET cap AND both directional deviations (before/after) confirm. The step iterates, removing one value per round (the worst outlier), until no candidates remain or fewer than 3 values are left.

***Key terms and variable names:***

- **EWMA (Exponentially Weighted Moving Average):** Weighted average of other measurements, with weights proportional to |age difference|^e. With e=-5, nearest neighbor dominates.
- **dewma_all:** EWMA deviation using all other values
- **dewma_before / dewma_after:** EWMA deviation excluding the immediately prior / next measurement. Used for the 90% confirmation rule.
- **min_gap_months:** Minimum gap (in months) to either chronological neighbor. Determines which tier (≤6m or >6m) applies.
- **Upper weight (UW):** `max(wt, ewma)` for each observation. Used for UW-based ET cap scaling.
- `remove_ewma_wt()` — support function implementing the iterative exclusion logic (shared with Moderate EWMA, Step 11Wb). Signature: `remove_ewma_wt(subj_df, wtallow_formula, ...)`
- `compute_et_limit()` — computes the ET cap per observation: `compute_et_limit(min_gap_months, formula, uw)`. See `wtallow-formulas.md` for full specification.
- `adult_ewma_cache_init()` — EWMA computation with position-based window (default 15 on each side)

***Configurable parameter defaults and options:***

| Parameter | loosest/looser | tighter | tightest | Unit | Notes |
|-----------|------|------|------|------|------|
| `wtallow_formula` | `"piecewise"` | `"piecewise-lower"` | `"allofus15"` | — | Determines ET cap structure (ET cap = wtallow cap + 20 for PW-H/PW-L; see wtallow-formulas.md) |
| `repval_handling` | `"independent"` | `"linked"` | `"linked"` | — | `"independent"` or `"linked"` (controls single vs two-pass) |
| `ewma_window` | 15 | 15 | 15 | observations | Max observations on each side for EWMA (not permissiveness-controlled) |

The threshold is computed dynamically per observation: `compute_et_limit(min_gap_months, formula, uw)` where `uw = max(wt, ewma)`. ET caps scale with UW — for higher UW they increase (PW-H only, capped at UW=180), for lower UW they decrease (PW-H and PW-L, with 2/3×UW ceiling). See `wtallow-formulas.md` for complete formula details.

***Logic and implementation:***

**Guards:**
1. Skip if no weight values remain
2. Skip if <3 non-extraneous values
3. **Range gate:** Skip if max-min of included values ≤ minimum possible ET cap + 0.12 (computed via `compute_et_limit(1, formula, uw = min(weights))`). If the entire range is within the smallest possible cap plus rounding tolerance, no value can exceed any cap.

**Independent mode (single pass):**
1. Call `remove_ewma_wt()` with all non-extraneous values (including RVs)
2. Each round:
   a. Compute minimum neighbor gap in months for each value. Missing gaps (edge values) → Inf (Stata missing-as-infinity convention).
   b. Compute per-observation threshold: `compute_et_limit(min_gap_months, formula, uw)` where `uw = max(wt, ewma_all)`
   c. Compute EWMA (exponent e=-5, anchor a=0, window = ewma_window observations on each side)
   d. **Positive outlier:** dewma_all > threshold + 0.12 AND dewma_before > 0.9×threshold + 0.12 AND dewma_after > 0.9×threshold + 0.12
   e. **Negative outlier:** dewma_all < -(threshold + 0.12) AND dewma_before < -(0.9×threshold + 0.12) AND dewma_after < -(0.9×threshold + 0.12)
   f. Missing directional dewma (edge values) → Inf, confirming exclusion
   g. Among candidates, select the one with largest |dewma_all|, then earliest age, then lowest id
   h. Exclude with code `"Exclude-A-WT-Traj-Ext-N"`
   i. Rebuild EWMA for next round (window boundaries shift after removal)
3. Iterate until no candidates or <3 values remain (no artificial round limit)
4. After exclusions: re-identify RVs on remaining values

**Linked mode (two-pass):**
1. **firstRV pass:** Run `remove_ewma_wt()` on non-RV values only, with label `"Exclude-A-WT-Traj-Extreme-firstRV"`. Apply same range gate (cap_short + 0.12) to firstRV subset.
2. Propagate firstRV exclusions to their RV copies via `propagate_to_rv()`
3. Re-identify RVs on remaining values
4. **allRV pass:** Run `remove_ewma_wt()` on all remaining non-extraneous values (RVs now participate), with label `"Exclude-A-WT-Traj-Extreme-allRV"`. Apply range gate (minimum ET cap + 0.12) to remaining subset.
5. Re-identify RVs after allRV exclusions

***Rationale for selected decisions:***

- **No round limit:** Unlike the original Stata implementation (capped at 3 rounds), R iterates until convergence. In practice, extreme EWMA rarely exceeds 3 rounds, but removing the artificial cap ensures all genuine outliers are caught.
- **90% rule:** The directional confirmation (before/after must exceed 90% of threshold) prevents excluding values that are outliers in one direction but consistent with the trend in the other. With e=-5, the 90% rule is effectively unreachable (nearest neighbor dominates so heavily that dewma_before/after are nearly equal to dewma_all), but it is retained at zero cost for consistency with the algorithm's theoretical framework.
- **Missing-as-infinity for edge values:** When a value has no prior or next neighbor, the missing directional dewma is treated as confirming exclusion. An extreme outlier at the edge of a patient's data should still be excludable based on the neighbors that do exist.
- **Two-tier structure:** Both Evil Twins and Extreme EWMA use the same ET cap structure derived from the wtallow formula (ET cap = wtallow cap + 20 for PW-H/PW-L). The distinction between ET and Extreme EWMA is the detection method (adjacent pairs vs EWMA deviation), not the cap structure.
- **Sort tiebreaker:** The R implementation includes `internal_id` as a final tiebreaker in all passes for deterministic sort order.

---

## Step 10H: Height Distinct Values

***Evaluate subjects with 2 or 3+ distinct height values: exclude pairs outside the height band (with loss/gain and frequency rescue), or values outside the best height window.***

| | |
|---|---|
| Scope | Height only |
| Prior Step | 9H: Final SDE Height Resolution |
| Next Step | 10W: Final SDE Weight Resolution |
| Distinct/RVs | Height does not use RV tracking |
| Exclusion Codes | `Exclude-A-HT-Ord-Pair-All`, `Exclude-A-HT-Ord-Pair`, `Exclude-A-HT-Window-All`, `Exclude-A-HT-Window` |
| Configurable parameter names | `ht_band`, `allow_ht_loss`, `allow_ht_gain` |

***Overview:***

This step evaluates subjects whose remaining included heights have 2 or more distinct values. It routes to one of two substeps:

- **10Ha (2 distinct):** Check if the two values differ by more than `ht_band` inches. If so, attempt loss rescue (height decreased within plausible limits), gain rescue (young adult height gain within allowance), and frequency rescue (one value appears ≥4/3 as often as the other).
- **10Hb (3+ distinct):** Find the best height window (w2) — a range of `ht_band` inches that captures the most values with sufficient ratio to outliers. Values outside the best window are excluded, unless rescued by loss-group or gain-group analysis.

Subjects with only 1 distinct height skip this step entirely.

***Key terms and variable names:***

- **ht_band:** Allowed height range in inches (default 3). Used with cm conversion and 0.12 cm rounding tolerance (see Rounding Tolerance section): effective threshold = `ht_band * 2.54 + 0.12` cm.
- **pairhtloss / pairhtgain:** Flags for whether a 2D pair qualifies for loss or gain rescue
- **htallow:** Height growth allowance based on age and velocity formula: `velocity * (ln(age2 - 16.9) - ln(age1 - 16.9))`
- **w2 window:** A height range of `ht_band * 2.54 + 0.12` cm starting from each distinct height
- **w2/o2 ratio:** Number of values inside vs. outside a window; must be ≥ 3/2 for the window to be viable
- **Loss group / gain group:** Sequential grouping of heights into clusters within 2 inches (5.08 + 0.12 cm) of each other
- `ht_allow()` — support function computing height growth allowance
- `ht_change_groups()` — support function forming sequential height groups
- `ht_3d_growth_compare()` — support function validating gain patterns across groups
- `check_between()` — support function for inclusive range checks

***Configurable parameter defaults and options:***

| Parameter | loosest | looser | tighter | tightest | Unit | Notes |
|-----------|------|------|------|------|------|------|
| `ht_band` | 3 | 3 | 2 | 2 | inches | Effective threshold in cm: `ht_band * 2.54 + 0.12` |
| `allow_ht_loss` | `TRUE` | `FALSE` | `FALSE` | `FALSE` | — | When FALSE, disables loss rescue in both 2D and 3D |
| `allow_ht_gain` | `TRUE` | `TRUE` | `TRUE` | `FALSE` | — | When FALSE, disables gain rescue in both 2D and 3D |

***Logic and implementation:***

**Step 10Ha: 2 Distinct Heights**

1. Extract the two distinct heights (ht_1, ht_2) in age order (both in cm)
2. **Band check:** `abs(ht_1 - ht_2) > (ht_band * 2.54 + 0.12)` — if FALSE, values are within band; skip (all remain included)
3. If band exceeded:
   a. **Gain rescue (pairhtgain):** Only evaluated if `ageyears1 < 25` (young adult). Checks:
      - Height increased: `(ht_2 - ht_1) > -0.12` (in cm, with rounding tolerance)
      - Within allowance: `(ht_2 - ht_1) ≤ (htallow + 2) * 2.54 + 0.12` (htallow in inches, converted to cm with tolerance; 2-inch buffer)
      - Temporal ordering: all ht_2 measurements occur after all ht_1 measurements
      - htallow coefficient by age gap: <2yr → 15.5, 2–3yr (inclusive) → 13, >3yr → 12
      - Age capped at 25 for htallow formula
   b. **Loss rescue (pairhtloss):** Only evaluated if `allow_ht_loss = TRUE`. Checks:
      - Height decreased: `ht_1 > ht_2 - 0.12` (in cm, with rounding tolerance)
      - Within loss limit: `(ht_1 - ht_2) ≤ 5 * 2.54 + 0.12` for age at ht_2 < 50, or `≤ 7 * 2.54 + 0.12` for age at ht_2 ≥ 50
      - Temporal ordering: all ht_2 measurements occur after all ht_1 measurements
   c. If neither rescue applies: mark all values for exclusion
   d. **Frequency rescue:** If `keepht1` (ht_1 count ≥ ht_2 count × 4/3), rescue ht_1 values. Same for keepht2.
4. **Exclusion code:**
   - `Exclude-A-HT-Ord-Pair` if frequency rescue saved some values
   - `Exclude-A-HT-Ord-Pair-All` if no frequency rescue (all excluded)

**Step 10Hb: 3+ Distinct Heights**

***Part 1 — Window (w2) evaluation:***
1. Define w2 window for each distinct height: `[ht, ht + ht_band * 2.54 + 0.12]`
2. Count observations and distinct values inside each window (w2) and outside (o2)
3. Require w2/o2 ratio ≥ 3/2 for a window to be viable
4. If multiple viable windows: select by highest score (`tot_w2 + 0.5 * numdistinct_w2`); ties broken by largest mean absolute distance from w2 mean to o2 values
5. If one viable window: use it
6. If no viable window: exclude all values (`Exclude-A-HT-Window-All`)
7. Exclude values outside the best window (`Exclude-A-HT-Window`)

***Part 2 — Loss group rescue (if `allow_ht_loss = TRUE`):***
1. Only evaluated if any values are still marked for exclusion after w2 evaluation
2. Form sequential height groups by iterating through values in age order:
   - Current group range within 2 inches (5.08 + 0.12 cm): extend group
   - New value > 2 inches below group minimum: start new loss group
   - New value > 2 inches above group maximum: abort (wrong direction for loss)
   - Maximum 3 groups; abort if exceeded
3. Validate loss pattern (all comparisons in cm with 0.12 tolerance):
   - G2-G1 direction: must decrease by at least 0.12 cm
   - G2-G1 magnitude: loss ≤ 5 inches (age < 50) or ≤ 7 inches (age ≥ 50)
   - G3-G2: same direction and magnitude checks
   - G3-G1 cumulative: ≤ 6/8/9 inches depending on age pattern
4. If all checks pass: rescue all values (clear all exclusion criteria)

***Part 3 — Gain group rescue:***
1. Only evaluated if any values are still marked for exclusion AND min age < 25 AND not already rescued by loss
2. Form sequential height groups (same as loss, but maximum 6 groups; abort if height decreases)
3. Validate gain pattern using `ht_3d_growth_compare()`:
   - Adjacent check (compare = "before"): each group vs. previous group
   - Cumulative check (compare = "first"): each group vs. group 1
   - Direction: gain must not drop more than 0.12 cm
   - Magnitude: gain must not exceed htallow + 0.12 cm
   - htallow coefficient by age gap: <1yr → 20, 1-3yr → 15, >3yr → 12
4. If all checks pass: rescue all values

***Rationale for selected decisions:***

- **0.12 rounding tolerance:** See Rounding Tolerance section. Applied to all height and weight threshold comparisons throughout the algorithm.
- **Loss rescue guarded by `allow_ht_loss`:** This is a configurable parameter because height loss (e.g., from kyphosis, osteoporosis) is real but uncommon. Setting `allow_ht_loss = FALSE` makes the algorithm more conservative — any large height difference results in exclusion.
- **Gain rescue limited to age < 25:** Adults over 25 do not grow, so height gain rescue only applies to young adults who might still be growing.
- **Two exclusion codes per substep:** Distinguishing "all excluded" from "some rescued by frequency" (2D) or "no viable window" from "outside window" (3D) provides more granular information about why values were excluded.
- **Frequency rescue (4/3 ratio):** When one height value appears substantially more often, it is likely the correct measurement. The 4/3 ratio threshold prevents marginal frequency differences from driving rescue.
- **Loss groups use means, not individual values:** Group means are compared because individual measurements may have rounding noise; averaging within a group produces a more stable reference.

---

## Step 11H: Mean Height

***Calculate the mean height for each subject, respecting loss/gain group boundaries from Step 10H.***

| | |
|---|---|
| Scope | Height only |
| Prior Step | 10H: Height Distinct Values |
| Next Step | 11Wa: 2D Ordered Weight Pairs |
| Distinct/RVs | N/A — height does not use RV tracking |
| Exclusion Codes | None — this step calculates a value, not exclusions |
| Configurable parameter names | None |

***Overview:***

This step computes a representative mean height for each subject, used later for BMI calculations. The mean is calculated from `meas_m` (metric measurement values) of remaining included heights. How the mean is computed depends on the outcome of Step 10H:

- **Pair loss or gain (2D):** Each observation retains its individual original measurement as its "mean" (no averaging across the two distinct values, since they represent different physical states)
- **3D with validated loss groups:** Mean within each loss group (observations in the same group share one mean)
- **3D with validated gain groups:** Mean within each gain group
- **All other subjects:** Simple arithmetic mean of all included measurements

***Key terms and variable names:***

- `mean_ht` — output column: the calculated mean height for each observation
- `meas_m` — the metric measurement value (used for mean calculation)
- `loss_groups` / `gain_groups` — group assignments from Step 10H
- `pairhtloss` / `pairhtgain` — flags from Step 10Ha indicating pair rescue type

***Logic and implementation:***

1. Initialize `meanht` as NA for all height observations
2. If no included heights remain, leave as NA
3. Otherwise, compute mean height based on Step 10H outcome:
   - **pairhtloss or pairhtgain:** Use each observation's `meas_m` directly (no averaging)
   - **Loss groups validated:** For each group in `glist_loss`, compute `mean(meas_m)` for observations in that group; all observations in the group receive the same mean
   - **Gain groups validated:** Same logic using `glist_gain`
   - **Default:** Compute `mean(meas_m)` across all included heights; all observations receive the same value
4. Store in the `mean_ht` output column

***Rationale for selected decisions:***

- **Uses `meas_m` (metric measurement):** Since SDE resolution retains a keeper's original value without averaging (Step 9H), `meas_m` equals the raw recorded measurement for all surviving heights. The mean is based on these raw values, not any intermediate calculation.
- **Pair loss/gain keeps individual values:** When a subject has two distinct heights representing different physical states (e.g., height loss from aging), averaging them would produce a value that doesn't represent either state. Each observation keeps its own measurement.
- **Group-aware averaging:** Loss/gain groups represent different time periods with different true heights. Averaging within groups but not across them preserves the temporal signal.

---

## Step 10W: Final SDE Weight Resolution

***Resolve same-day weight duplicates: exclude identical values, then choose keepers using category-based median selection.***

| | |
|---|---|
| Scope | Weight only |
| Prior Step | 9Wb: Extreme EWMA |
| Next Step | 11H: Mean Height; 11Wa: Distinct Ordered Pairs (Weight) |
| Distinct/RVs | Non-SDE median uses firstRV (non-RV) values only; re-identifies RVs after exclusions |
| Exclusion Codes | `Exclude-A-WT-Identical`, `Exclude-A-WT-Extraneous` |
| Configurable parameter names | None |

***Overview:***

Same logic as Step 9H (Final SDE Height Resolution), applied to weight. Resolves same-day weight duplicates that were temporarily flagged in Step 3W. Two parts:

- **Part A (Identical):** Exclude all but one when same-day values are identical (keep lowest id).
- **Part B (Non-identical):** Categorize subject and choose keeper using category-specific median comparisons.

The keeper retains its original measurement value — no averaging is performed. After exclusions, RVs are re-identified since a first_rv may have been excluded.

***Key terms and variable names:***

Same as Step 9H, with one weight-specific difference:
- **nonsdemed (weight):** Uses firstRV (non-RV) values on non-SDE days only, unlike height which uses all values. This prevents carried-forward RV values from biasing the reference median.

***Configurable parameter defaults and options:***

None.

***Logic and implementation:***

Identical to Step 9H with these weight-specific differences:
1. After identical removal, `temp_sde()` is re-run with `ptype = "weight"` (non-RV median) and `redo_identify_rv()` updates RV flags
2. Non-SDE median (`nonsdemed`) is calculated from non-RV values on non-SDE days only
3. After all SDE exclusions, RVs are re-identified via `identify_rv()`

See Step 9H for full details on categories, sorting, and tiebreakers.

***Rationale for selected decisions:***

- Using firstRV (non-RV) values for the non-SDE median prevents carried-forward values from dominating the reference. Height does not need this because repeated identical heights are expected (adults don't grow).
- All other rationale is the same as Step 9H.

---

## Step 11Wa: 2D Ordered Weight Pairs

***Exclude subjects with exactly 2 distinct weight values in time order when the weight difference exceeds the allowed amount or the weight ratio is too extreme.***

| | |
|---|---|
| Scope | Weight only |
| Prior Step | 10W: Final SDE Weight Resolution |
| Next Step | 11Wa2: 2D Non-Ordered Weight Pairs |
| Distinct/RVs | Routing uses firstRV (non-RV) numdistinct and ordering; allRV override for non-ordered detection; evaluation uses firstRV values |
| Exclusion Codes | `Exclude-A-WT-2D-Ordered` |
| Configurable parameter names | `wtallow_formula` |

***Overview:***

This step handles subjects whose remaining included weights have exactly 2 distinct values that are time-ordered (all instances of value 1 occur before all instances of value 2). It applies two independent tests: (1) whether the absolute weight difference exceeds a time-dependent allowance (`wtallow`), and (2) whether the ratio of the smaller weight to the larger weight falls below a threshold (`perclimit`). If either test fails, all weight observations for the subject are excluded.

This step also performs the routing that determines whether subjects go to 2D ordered (this step), 2D non-ordered (Step 11Wa2), or 3+D/EWMA (Step 11Wb). Subjects with 1 distinct value go to Step 13 (Distinct Single).

***Key terms and variable names:***

- **numdistinct** — number of distinct weight values among firstRV (non-RV) included observations
- **pair_distinct** — TRUE if firstRV has exactly 2 distinct values and they are time-ordered (max age of value 1 < min age of value 2)
- **is_2d_nonord** — TRUE if allRV has exactly 2 distinct values but they are NOT time-ordered (overrides pair_distinct)
- **wt_first / wt_last** — first and last weight values (by age) among firstRV observations
- **wt_diff** — `wt_last - wt_first`
- **ageyears_diff** — time gap between the two value groups: `min age of value 2 - max age of value 1` (in years)
- **wtallow (wta)** — maximum allowed weight change for the given time interval, computed by `compute_wtallow(months, formula, uw)`. UW is the upper weight of the pair (`max(wt_first, wt_last)`).
- **wt_perc** — ratio of smaller to larger weight: `min(wt_first, wt_last) / max(wt_first, wt_last)`
- **perclimit** — threshold for wt_perc; subject-level, max of `compute_perc_limit()` across observations. Weight-dependent and permissiveness-dependent (see permissiveness table above).

***Configurable parameter defaults and options:***

- **`wtallow_formula`** (default `"piecewise"`): Selects the formula for computing weight allowance from the time interval between measurements. All formulas are UW-adjusted (scaled by the upper weight of the pair). Options:

  - **`"piecewise"` (PW-H)** — Three-segment curve. Default for loosest/looser. At UW=120 (base):
    - 0–1 month: log-curved from ~10.85 kg to 20 kg (`10 + 10 × log(1 + 5×months) / log(6)`)
    - 1–6 months: linear from 20 kg to 50 kg
    - 6–12 months: linear from 50 kg to 80 kg
    - >12 months: flat at 80 kg
    - UW > 120: caps increase by 0.25×(UW−120), capped at UW=180; 1m value = 25
    - UW < 120: caps scale as `(base_cap − 20) × (UW/120) + 20`; ceiling of UW × 2/3

  - **`"piecewise-lower"` (PW-L)** — Same shape as PW-H but with lower caps. Default for tighter. At UW=120 (base):
    - 0–1 month: same log curve as PW-H
    - 1–6 months: linear from 20 kg to 33.33 kg
    - 6–12 months: linear from 33.33 kg to 53.33 kg
    - >12 months: flat at 53.33 kg
    - UW > 120: NOT scaled up (same as base formula)
    - UW < 120: caps scale as `(base_cap − 20) × (UW/120) + 20`; ceiling of UW × 2/3

  - **`"allofus15"`** — Step function for short intervals, linear ramp for long intervals. Default for tightest.
    - 0–2 days: 5 kg; 3–7 days: 10 kg; 8 days–<6 months: 15 kg
    - 6–12 months: linear from 15 kg to allofus15-cap-12m
    - >12 months: flat at allofus15-cap-12m
    - allofus15-cap-12m = min(40, effective PW-L at 12m for that UW)

  - Custom CSV path with columns `months` and `wtallow` (linearly interpolated, no UW adjustment).

See `wtallow-formulas.md` for the complete specification of all formulas, UW adjustments, ET caps, and ceilings.

***Logic and implementation:***

**Routing (numdistinct and ordering):**

1. Filter to firstRV (non-RV) included weight observations (`w_nonrv`)
2. Count distinct values among firstRV observations
3. If firstRV has exactly 2 distinct values, check ordering: `max(ages of value 1) < min(ages of value 2)` (strict less-than)
4. **allRV override:** If allRV (all included values, including RVs) also has exactly 2 distinct values but they are NOT time-ordered, override the firstRV result — route to 2D non-ordered (Step 11Wa2) instead
5. Routing outcome:
   - 1 distinct value → Step 13 (Distinct Single)
   - 2 distinct, ordered → this step (11Wa)
   - 2 distinct, not ordered → Step 11Wa2
   - 3+ distinct → Step 11Wb (Moderate EWMA)

**Evaluation (if 2D ordered):**

1. Compute `wt_first` (earliest weight) and `wt_last` (latest weight) from firstRV observations
2. Compute `ageyears_diff` = `min age of value 2 - max age of value 1` (gap between the two groups)
3. Compute `wtallow` using `compute_wtallow(ageyears_diff * 12, formula, uw = max(wt_first, wt_last))`
4. Test 1 (wtallow): `abs(wt_diff) > wtallow + 0.12` (with rounding tolerance)
5. Compute `wt_perc` = `min(wt_first, wt_last) / max(wt_first, wt_last)`
6. Determine `perclimit`: subject-level, max of `compute_perc_limit()` across all observations. At loosest/looser: ≤45 kg → 0.5, 45–80 kg → 0.4, >80 kg → 0 (disabled). At tighter/tightest: ≤45 kg → 0.7, >45 kg → 0.4.
7. Test 2 (perclimit): `wt_perc < perclimit` (strict less-than)
8. If either test fails, exclude ALL weight observations for the subject (including RVs)

***Rationale for selected decisions:***

- **firstRV for routing, allRV for non-ordered override:** firstRV gives the correct numdistinct count (RVs are duplicates of existing values). However, firstRV with 2 distinct values always appears ordered (each value has one observation). The allRV check catches cases where interleaved RVs reveal that the values aren't truly time-separated.
- **Subject-level perclimit:** Uses the maximum `compute_perc_limit()` across all observations. At loosest/looser, perclimit is disabled (0) for weights >80 kg, so the subject-level value is driven by any lighter observations.
- **Two independent exclusion criteria:** wtallow catches absolute differences that are too large for the time interval. perclimit catches extreme ratios (e.g., 40 kg vs 100 kg) even when the absolute difference might be within wtallow for a long interval.
- **Parameterized wtallow:** Different datasets may warrant different sensitivity. The piecewise formula (loosest/looser default) is more permissive at short intervals and stricter at long intervals; tighter/tightest use lower-cap variants.

---

## Step 11Wa2: 2D Non-Ordered Weight Pairs

***Exclude subjects with exactly 2 distinct weight values that are interleaved in time, using a four-rule decision tree based on wtallow, prior exclusions, and value dominance.***

| | |
|---|---|
| Scope | Weight only |
| Prior Step | 11Wa: 2D Ordered Weight Pairs |
| Next Step | 11Wb: Moderate EWMA |
| Distinct/RVs | Uses allRV (all included values) for both detection and evaluation |
| Exclusion Codes | `Exclude-A-WT-2D-Non-Ordered` |
| Configurable parameter names | `wtallow_formula` |

***Overview:***

This step handles subjects whose remaining included weights have exactly 2 distinct values that are NOT time-ordered — meaning the values are interleaved (e.g., 70, 110, 70, 110). Unlike 2D ordered pairs (Step 11Wa), the values cannot be separated into "before" and "after" groups, so a different evaluation strategy is used.

The step applies a four-rule decision tree: (1) if all adjacent different-value pairs are within wtallow, keep everything; (2) if any pair exceeds wtallow and the subject has prior non-SDE exclusions, exclude all; (3) if one value is dominant (>65% of observations), exclude only the minority; (4) if neither value dominates, exclude all.

***Key terms and variable names:***

- **any_outside** — TRUE if any adjacent pair of different-value observations has a weight difference exceeding wtallow for that pair's time gap
- **prior_nonSDE** — TRUE if the subject has any prior weight exclusions that are not SDE-related (not `Include`, `temp extraneous`, `Exclude-A-HT-Identical`, `Exclude-A-WT-Identical`, `Exclude-A-HT-Extraneous`, `Exclude-A-WT-Extraneous`, empty string, or NA)
- **dominant_pct** — proportion of observations belonging to the more common value: `max(count_v1, count_v2) / total`
- **minority_val** — the less common of the two distinct values

***Configurable parameter defaults and options:***

Same as Step 11Wa — uses `compute_wtallow()` with `wtallow_formula`. UW for each adjacent pair is the larger of the two weights in the pair (`uw_pair`).

***Logic and implementation:***

1. Sort observations by age (tiebreaker: id)
2. For each adjacent pair where the weight values differ:
   - Compute time gap in months: `abs(age_days[i+1] - age_days[i]) / 30.4375`
   - Compute wtallow for that gap using `compute_wtallow()`
   - If `abs(wt[i+1] - wt[i]) > wtallow + 0.12` → set `any_outside = TRUE` and stop checking
   - Skip adjacent pairs with identical weight values (diff = 0, never outside wtallow)
3. Apply decision rules:
   - **Rule 1:** `any_outside = FALSE` → all observations stay included (no action)
   - **Rule 2:** `any_outside = TRUE` AND subject has prior non-SDE exclusions → exclude ALL observations
   - **Rule 3:** `any_outside = TRUE`, no prior non-SDE exclusions, AND one value is dominant (`dominant_pct > 0.65`, strict) → exclude only the minority value observations
   - **Rule 4:** `any_outside = TRUE`, no prior non-SDE exclusions, AND neither value dominates (`dominant_pct <= 0.65`) → exclude ALL observations

***Prior non-SDE exclusion check:***

Checks the full subject weight history (`w_subj_keep`), not just currently included observations. Uses exact code matching (not substring search) to identify SDE codes. Codes that do NOT count as "prior non-SDE": `Include`, `temp extraneous`, `Exclude-A-HT-Identical`, `Exclude-A-WT-Identical`, `Exclude-A-HT-Extraneous`, `Exclude-A-WT-Extraneous`, empty string, NA. Everything else — including `Exclude-A-WT-Scale-Max-Identical` — counts as a prior non-SDE exclusion.

***Rationale for selected decisions:***

- **Interleaved values need different logic than ordered:** Ordered pairs have a clear "before vs. after" that can be evaluated with a single wtallow check. Interleaved values suggest either real oscillation (keep) or data quality issues (exclude), so the decision tree weighs multiple factors.
- **Prior exclusions as a signal:** If the subject already has non-SDE exclusions, there's existing evidence of data quality issues, so the threshold for excluding remaining values is lower.
- **Dominance rule (65%):** When one value clearly predominates, the minority value is more likely to be erroneous. The 65% threshold prevents excluding values that are roughly equally represented, where the "error" is ambiguous.
- **All-or-nothing for balanced splits:** When neither value dominates and the difference exceeds wtallow, both values are suspect and all are excluded.

---

## Step 11Wb: Moderate EWMA

***Exclude weight outliers using a 7-step EWMA flow with trajectory rescue, alternate pathway, percentage criterion, and error load detection.***

| | |
|---|---|
| Scope | Weight only |
| Prior Step | 11Wa2: 2D Non-Ordered Weight Pairs |
| Next Step | 13: Distinct Single (1D) |
| Distinct/RVs | Independent: single pass, RVs as full participants. Linked: firstRV pass (non-RV only) with RV propagation, then allRV pass. |
| Exclusion Codes | `Exclude-A-WT-Traj-Moderate-N`, `Exclude-A-WT-Traj-Moderate-allRV-N`, `Exclude-A-WT-Traj-Moderate-Error-Load-N`, `Exclude-A-WT-Traj-Moderate-Error-Load-RV-N` |
| Configurable parameter names | `wtallow_formula`, `ewma_window`, `repval_handling`, `mod_ewma_f`, `perclimit_low`, `perclimit_mid`, `perclimit_high` |

***Overview:***

This step handles subjects with 3+ distinct weight values. It uses EWMA deviation to identify outliers, but with a more nuanced 7-step flow compared to Extreme EWMA (Step 9Wb). The key differences from Extreme EWMA are:

- Uses **wtallow** (time-interval-dependent) as the threshold instead of fixed caps
- Includes **trajectory rescue** — values that follow the interpolation/extrapolation of surrounding measurements are protected from exclusion
- Has an **alternate pathway** for values adjacent to unreliable neighbors (≤14 days, >wtallow apart)
- Includes a **percentage criterion** (wt/EWMA < perclimit) as an independent exclusion pathway
- Detects **error load** (4+ consecutive candidates) to flag subjects with widespread data quality issues
- Uses **graduated scoring** to prioritize which candidate to exclude when multiple exist

Each round excludes at most one non-error-load value (the highest-scoring candidate). Error load values are all excluded immediately. The step iterates until no candidates remain, fewer than 3 values are left, or the maximum round limit (100) is reached.

***Key terms and variable names:***

- **dewma_all / dewma_bef / dewma_aft** — EWMA deviation using all, all-except-predecessor, all-except-successor
- **wtallow (wta)** — time-interval-dependent weight allowance, computed from `compute_wtallow(months, formula, uw)` using the minimum of the before and after age gaps. UW = `max(wt, ewma)` for the observation.
- **wta_base** — wtallow at UW=120 (no UW adjustment), used as the denominator in prioritization scoring. This ensures scoring is not distorted by UW-dependent wtallow.
- **minagediff** — minimum of before and after age gaps (years), used to determine wtallow
- **perclimit** — observation-level ratio threshold via `compute_perc_limit(meas, permissiveness)`. Weight-dependent and permissiveness-dependent (see permissiveness table). Note: observation-level here, unlike Step 11Wa which uses subject-level max.
- **exc_stand** — standard pathway candidates (Step 1)
- **exc_pair** — alternate pathway candidates (Step 3)
- **exc_wt_i** — accumulated exclusion candidates across Steps 2, 4, 5
- **error_load** — observations in 4+ consecutive exc_wt_i runs
- **ewma_window** — number of observations on each side used for EWMA computation (default 15, so up to 30 total neighbors)
- **trajectory rescue** — three methods (interpolation, prior extrapolation, next extrapolation) that can protect a value from standard-pathway exclusion

***Configurable parameter defaults and options:***

| Parameter | Default | Unit | Notes |
|-----------|---------|------|-------|
| `wtallow_formula` | `"piecewise"` | — | Formula for wtallow (see Step 11Wa and wtallow-formulas.md) |
| `mod_ewma_f` | 0.75 | — | Directional factor: loosest/looser = 0.75; tighter/tightest = 0.60 |
| `perclimit_low` | 0.5 | — | % criterion for wt ≤45 kg (loosest/looser); 0.7 (tighter/tightest) |
| `perclimit_mid` | 0.4 | — | % criterion for 45<wt≤80 kg (all levels) |
| `perclimit_high` | 0.0 | — | % criterion for wt >80 kg; 0 = disabled (loosest/looser); 0.4 (tighter/tightest) |
| `ewma_window` | 15 | observations | Max observations on each side for EWMA. All EWMAs use this window. |
| `repval_handling` | `"independent"` | — | Controls single vs two-pass EWMA |
| `max_rounds` | 100 | — | Maximum EWMA iteration rounds |

***Logic and implementation:***

**Guards:**
1. Only subjects with 3+ distinct weight values enter this step (others went to 11Wa, 11Wa2, or 13)
2. Must have ≥ 3 included observations
3. **Trigger check:** Skip if weight range ≤ min wtallow AND (perclimit disabled OR min percent ratio ≥ max perclimit) — impossible to trigger any exclusion

**EWMA computation:**
- Exponent e = -5, anchor a = 0 (same as Extreme EWMA)
- Position-based window: only the nearest `ewma_window` (default 15) observations on each side contribute
- Rebuilt from scratch each round (window boundaries shift when observations are removed)

**7-step exclusion flow (each round):**

**Step 1 — Standard pathway + trajectory rescue:**
- Standard criteria: `|dewma_all| > wtallow + 0.12` AND both directional dewma exceed `0.75 × wtallow + 0.12` in the same direction
- Trajectory rescue checks three methods (all use ±5 kg error margin):
  1. **Interpolation:** Value falls between p1 and n1 (± 5 kg). Uses non-strict inequalities (≥, ≤).
  2. **Extrapolation prior:** Value falls within the range defined by p2 and the linear extrapolation from p2→p1. Rounded to nearest 0.2. Distance guard: extrapolation distance must be ≤ 2× source interval.
  3. **Extrapolation next:** Same as prior, using n1→n2 backwards. Same distance guard.
- A value is rescued if ANY of the three methods succeeds. Must fail ALL three to remain a candidate.
- `exc_stand = standard_criteria AND fails_all_trajectory`

**Step 2 — Standard run detection + pair/trio prioritization:**
- Detect consecutive runs of `exc_stand` within each subject
- **4+ consecutive** → error load (all excluded immediately in Step 7)
- **Isolated (run length 1)** → directly becomes exc_wt_i
- **Pairs/trios (run length 2-3)** → score each member and pick the highest:
  - First of run: `|dewma_aft / wtallow_aft_unrel|` (wtallow computed as if before-neighbor were removed)
  - Last of run: `|dewma_bef / wtallow_bef_unrel|`
  - Middle of trio: `|dewma_bef/wtallow_bef_unrel + dewma_aft/wtallow_aft_unrel|`
  - Highest score → exc_wt_i

**Step 3 — Alternate pathway:**
- Catches values adjacent to unreliable neighbors (≤14 days AND >wtallow apart)
- If prior is unreliable: only check aft dewma (`|dewma_all| > wtallow + 0.12` AND `|dewma_aft| > 0.75×wtallow + 0.12`, same direction)
- If next is unreliable: only check bef dewma (same logic)
- Must not already be exc_wt_i
- Result: `exc_pair`

**Step 4 — Alternate run detection + pair/trio prioritization:**
- Same run detection as Step 2, applied to exc_pair candidates
- Isolated → directly exc_wt_i
- Pairs/trios → scored with same approach, accounting for committed neighbors (exc_wt_i from prior steps that are adjacent to the pair/trio)
- Effective position adjusts when a committed neighbor is adjacent (first becomes middle, etc.)

**Step 5 — Percentage criterion:**
- Independent of EWMA deviation direction
- `wt/ewma_all < perclimit AND wt/ewma_bef < perclimit AND wt/ewma_aft < perclimit`
- Must not already be exc_wt_i
- Uses observation-level perclimit via `compute_perc_limit(meas, permissiveness)`

**Step 6 — Consecutive check (all pathways combined):**
- After Steps 2, 4, 5 accumulate exc_wt_i candidates, re-check for 4+ consecutive
- Newly identified 4+ runs → error load

**Step 7 — Final prioritization:**
- **Error load values:** All excluded immediately with code `"...-Error-Load round N"`
- **Non-error-load candidates:** Score each and exclude only the single highest-scoring:
  - Edge (first/last in subject): `pmax(0, |dewma_all| - wta) / wta_base` where `wta` is UW-adjusted wtallow and `wta_base` is wtallow at UW=120 (no adjustment)
  - Interior: graduated multiplier × `(excess_bef + excess_aft)` where each `excess = pmax(0, |dewma| - wta) / wta_base`
    - `min_excess = min(excess_bef, excess_aft)`
    - `multiplier = min(1.0, 0.6 + 0.4 × min_excess)`
    - Full score when both directions independently exceed wtallow; 0.6× discount when only one direction is bad
  - Scoring uses excess above UW-adjusted wtallow divided by base wtallow, so a 10 kg excess at any UW always scores the same
  - Tiebreakers: smallest wtallow → closest to median age → earliest position → lowest id

**Mode handling:**

**Independent mode:** Single pass with all values (including RVs) as full participants.

**Linked mode:**
1. **firstRV pass:** Run `remove_mod_ewma_wt()` on non-RV values only. Label: `"Exclude-A-WT-Traj-Moderate"`.
   - **Error load RV escalation:** If any error-load value has RV copies, escalate the entire patient — all remaining Include values get `"Exclude-A-WT-Traj-Moderate-Error-Load-RV-N"`.
   - **Non-error-load propagation:** Propagate exclusions to RV copies via `propagate_to_rv()` (appends `-RV-propagated` suffix).
   - Re-identify RVs on remaining values.
2. **allRV pass:** Run `remove_mod_ewma_wt()` on all remaining values (RVs now participate). Label: `"Exclude-A-WT-Traj-Moderate-allRV"`.

***Rationale for selected decisions:***

- **wtallow instead of fixed caps:** Moderate EWMA handles subjects with 3+ distinct values over varying time spans. Fixed caps (appropriate for Extreme EWMA's large outliers) would be too coarse; wtallow scales with the time interval to catch subtler errors.
- **Trajectory rescue:** Values that follow the linear trend of surrounding measurements are likely real, even if their EWMA deviation is large (e.g., a genuine rapid weight change). Rescue prevents false positives. Must fail all three methods to ensure the value genuinely doesn't fit any trajectory.
- **Alternate pathway:** When a neighbor is ≤14 days away and >wtallow apart, that neighbor is unreliable and distorts the directional dewma on that side. The alternate pathway drops the unreliable side's confirmation requirement.
- **Graduated multiplier for interior scoring:** Interior values have both before and after context, while edge values have only one direction. The graduated multiplier bridges these: when both directions independently confirm (min_ratio ≥ 1), the score is full; when only one direction is bad, the score is discounted (0.6× floor), preventing interior values from being unfairly ranked against edge values.
- **Observation-level perclimit (not subject-level):** Unlike Step 11Wa where the max perclimit across observations is used, here each observation uses its own weight and the current permissiveness level to determine its threshold via `compute_perc_limit()`. At loosest/looser, observations >80 kg have perclimit=0 (disabled).
- **Error load escalation in linked mode:** If a patient has both repeated values (RVs) and 4+ consecutive flagged values, the data quality is too poor to trust any remaining values.
- **EWMA window (default 15):** Limits each EWMA to the 15 nearest observations on each side. For subjects with many observations, distant measurements (which contribute negligibly due to the e=-5 exponent) are excluded from the computation, improving both efficiency and preventing numerical noise from very distant observations.
- **Max rounds (default 100):** Allows the algorithm to handle subjects with many observations without being artificially limited. In practice, convergence typically occurs within a few rounds.

---

## Step 13: Single Distinct Value (1D) Evaluation

***Exclude height or weight values when a subject has only one distinct value for that parameter, using BMI-dependent or no-BMI limits.***

| | |
|---|---|
| Scope | Height and Weight |
| Prior Step | 11Wb: Moderate EWMA |
| Next Step | 14: Error Load |
| Distinct/RVs | Weight uses firstRV for numdistinct. Remaining RVs are treated as Include. |
| Exclusion Codes | `Exclude-A-HT-Single`, `Exclude-A-WT-Single` |
| Configurable parameter names | `single_ht_min_bmi`, `single_ht_max_bmi`, `single_wt_min_bmi`, `single_wt_max_bmi`, `single_ht_min_nobmi`, `single_ht_max_nobmi`, `single_wt_min_nobmi`, `single_wt_max_nobmi`, `single_bmi_min`, `single_bmi_max` |

***Overview:***

After all prior steps, some subjects may have only one distinct value remaining for height or weight. These cannot be evaluated by EWMA or pair comparisons, so this step applies simple range limits. The limits depend on whether BMI can be computed (i.e., the subject has both a height and a weight on the same day):

- **BMI available and within range:** Use loose limits (wide range of acceptable values)
- **BMI available but extreme:** Exclude all 1D values for the parameter (BMI outside normal range suggests data quality issue)
- **No BMI available:** Use tighter limits (can't validate with BMI, so be more conservative)

***Key terms and variable names:***

- **1D** — subject has exactly 1 distinct included value for a parameter (height or weight)
- **BMI** — computed as `weight_kg / (height_cm / 100)²` from same-day height/weight pairs
- **bmi_extreme** — TRUE if BMI is outside `[single_bmi_min, single_bmi_max]` (strict `<` / `>`)
- **has_bmi** — TRUE if the subject has at least one day with both an included height and an included weight

***Configurable parameter defaults and options:***

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| `single_ht_min_bmi` (cm) | 60 | 120 | 142 | 147 |
| `single_ht_max_bmi` (cm) | 245 | 230 | 213 | 208 |
| `single_wt_min_bmi` (kg) | 12 | 30 | 36 | 39 |
| `single_wt_max_bmi` (kg) | 350 | 270 | 159 | 136 |
| `single_ht_min_nobmi` (cm) | 122 | 120 | 142 | 147 |
| `single_ht_max_nobmi` (cm) | 245 | 230 | 213 | 208 |
| `single_wt_min_nobmi` (kg) | 30 | 30 | 36 | 39 |
| `single_wt_max_nobmi` (kg) | 350 | 270 | 159 | 136 |
| `single_bmi_min` (kg/m²) | 10 | 12 | 16 | 18 |
| `single_bmi_max` (kg/m²) | 250 | 65 | 45 | 40 |

**Note:** For looser/tighter/tightest, the BMI-available and no-BMI limits are identical — BMI availability does not change which values are excluded at these levels. The BMI/no-BMI split only matters at loosest (wider limits when BMI is available: 60 cm vs 122 cm for ht_min; 12 kg vs 30 kg for wt_min). `single_bmi_min`/`single_bmi_max` match the overall BIV BMI limits at each level.

***Logic and implementation:***

The step runs as a two-pass loop over each subject (both passes use the same logic, but pass 2 catches cases where a pass 1 exclusion removes a same-day partner, causing BMI to become unavailable):

**Pass 1:**
1. Determine which parameters are 1D: exactly 1 distinct included value and at least 1 observation
2. If neither height nor weight is 1D, skip this subject
3. Check BMI availability: find days where both an included height and an included weight exist
4. **If BMI available:**
   - Compute BMI from same-day height/weight pair
   - If BMI is extreme (`bmi < single_bmi_min` or `bmi > single_bmi_max`, strict): exclude ALL 1D height observations AND all 1D weight observations
   - If BMI is within range: check each 1D parameter against the BMI-available limits (`single_ht_min_bmi`/`single_ht_max_bmi`, `single_wt_min_bmi`/`single_wt_max_bmi`). Exclude if outside limits (strict `<` / `>`).
5. **If no BMI available:** check each 1D parameter against the no-BMI limits (`single_ht_min_nobmi`/`single_ht_max_nobmi`, `single_wt_min_nobmi`/`single_wt_max_nobmi`). Exclude if outside limits (strict `<` / `>`).

**Pass 2:**
- Re-evaluate after pass 1 exclusions. If pass 1 excluded a height (or weight), the partner weight (or height) may have lost its BMI availability. In that case, re-evaluate using the tighter no-BMI limits.
- If BMI is still available, nothing new is checked (already handled in pass 1).

***Rationale for selected decisions:***

- **BMI-dependent limits:** When both height and weight are available on the same day, BMI provides a cross-validation that allows looser individual limits. Without BMI, tighter limits compensate for the lack of cross-validation.
- **BMI extreme → exclude all 1D:** If the combination of height and weight produces an implausible BMI, at least one must be wrong. Since both are 1D (no other values to compare against), there's no way to determine which, so both are excluded.
- **Two-pass structure:** A pass 1 BMI-extreme exclusion of height could leave weight without BMI availability. Pass 2 catches this by re-evaluating with tighter no-BMI limits.
- **Strict inequalities:** Values exactly at the limit boundaries are kept, consistent with other steps.

---

## Step 14: Error Load

***Exclude all remaining included values for a parameter when the proportion of errors is too high, indicating widespread data quality issues.***

| | |
|---|---|
| Scope | Height and Weight (evaluated separately) |
| Prior Step | 13: Single Distinct (1D) |
| Next Step | None (final step) |
| Distinct/RVs | Weight uses firstRV exclusion codes for error counting |
| Exclusion Codes | `Exclude-A-HT-Too-Many-Errors`, `Exclude-A-WT-Too-Many-Errors` |
| Configurable parameter names | `error_load_threshold` |

***Overview:***

This is the final step. It evaluates whether a subject has too many errors relative to the total number of meaningful observations. If the ratio of errors to (errors + included values) exceeds the threshold, all remaining included values for that parameter are excluded. Height and weight are evaluated independently — a subject can have error load for weight but not height, or vice versa.

***Key terms and variable names:***

- **error ratio** — `n_errors / (n_errors + n_included)` for one parameter within one subject
- **n_errors** — count of observations with error exclusion codes (see below)
- **n_included** — count of observations still marked "Include"
- **denominator** — `n_errors + n_included` (must be ≥ 3 to evaluate)
- **SDE exclusions** — Same-day identical and extraneous exclusions are NOT counted as errors and are excluded from the denominator entirely

***Configurable parameter defaults and options:***

| Parameter | Default | Description |
|-----------|---------|-------------|
| `error_load_threshold` | 0.41 | Error ratio above this triggers exclusion (strict `>`) |

***Which exclusion codes count as errors:***

**Height errors:**
- `Exclude-A-HT-BIV`
- `Exclude-A-HT-Single` (1D)
- `Exclude-A-HT-Ord-Pair`
- `Exclude-A-HT-Ord-Pair-All`
- `Exclude-A-HT-Window`
- `Exclude-A-HT-Window-All`

**Weight errors:**
- `Exclude-A-WT-BIV`
- `Exclude-A-WT-Single` (1D)
- `Exclude-A-WT-2D-Ordered`
- `Exclude-A-WT-2D-Non-Ordered`
- `Exclude-A-WT-Scale-Max`
- `Exclude-A-WT-Scale-Max-Identical`
- `Exclude-A-WT-Scale-Max-RV-Propagated`
- `Exclude-A-Evil-Twins`
- All `Exclude-A-WT-Traj-*` codes (Extreme and Moderate, all rounds and variants)

**NOT counted as errors:**
- Same-day identical / extraneous (SDE) — excluded from denominator entirely
- EWMA-RV-propagated codes — these represent the same underlying error as the source exclusion; counting both would inflate the ratio

***Logic and implementation:***

1. For each subject, process height and weight separately
2. Remove SDE observations from consideration (not errors, not in denominator)
3. Count errors (observations matching error codes above)
4. Count remaining included observations
5. Compute denominator = errors + included. If denominator < 3, skip (too few observations to evaluate)
6. Compute ratio = errors / denominator
7. If ratio > 0.41 (strict `>`): exclude all remaining included values for that parameter

***Rationale for selected decisions:***

- **Threshold 0.41 (not 0.40):** The 0.01 buffer (effectively 0.4 + 0.01) prevents borderline cases from triggering. A ratio of exactly 0.40 (e.g., 2 errors out of 5) is not considered excessive.
- **Minimum denominator of 3:** With fewer than 3 meaningful observations, the ratio is too unstable to be informative (e.g., 1 error out of 2 = 0.50, but this single remaining value may be perfectly valid).
- **SDEs excluded from denominator:** Same-day duplicates are a data collection artifact, not measurement errors. Including them would dilute the ratio.
- **EWMA-RV-propagated not counted:** In linked mode, when a firstRV value is excluded by EWMA, its RV copies are propagated. These represent the same measurement error — counting both the source and propagated copies would double-count the error.
- **Scale-Max-RV-propagated IS counted:** Unlike EWMA-RV-propagated, scale-max RV copies represent independent observations at the scale ceiling. Each contributes evidence of a data quality pattern.
- **Height and weight evaluated independently:** A subject with many height errors but good weights should only lose the heights, not the weights.

