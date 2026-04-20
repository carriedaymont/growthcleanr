# growthcleanr Parameters and Thresholds Reference

Reference Last Updated: 2026-04-18. Initial draft by: Claude (Opus), with Carrie Daymont.

Code Last Updated: 2026-04-18

Authoritative index of every user-configurable parameter, hardcoded numeric threshold, and reference-data-driven value used by the growthcleanr algorithm. Once established, this reference is the single source of truth — algorithm reference documents (`child-algorithm-reference.md`, `wrapper-narrative-2026-04-17.md`, `adult-algorithm-narrative.md`) and CLAUDE.md files should cite this reference rather than duplicating defaults.

When this reference disagrees with code, reference documents, or CLAUDE.md, raise for discussion — do not auto-reconcile to any one source. Code is a strong signal but not automatically correct; a mismatch can mean code, a reference, or the master table should change.

---

## How to Read This Reference

**Three categories:**

1. **User-configurable parameters** — exposed as arguments to `cleangrowth()` (and internally to `cleanchild()` / `cleanadult()`). A caller can change these. Organized by algorithm scope.
2. **Hardcoded thresholds** — numeric values baked into the code; a caller cannot change them without editing the package. Organized by step.
3. **Reference-data-driven values** — loaded from external files in `inst/extdata/` or defined in large lookup tables. Change the file, change the value. Some are user-overridable (e.g., `sd.recenter`, `ref.data.path`).

**Sections per category:**

- Parameter/threshold name (as written in code)
- Value (default for configurable; literal for hardcoded; "see file X" for ref-data)
- Scope (step, phase, or algorithm)
- Code location (file + function + sub-section)
- Brief role

**Conventions:**

- `<` / `>` = strict inequalities unless noted
- Param values: `"WEIGHTKG"`, `"HEIGHTCM"`, `"HEADCM"` (after input conversion)
- Sex: 0 = male, 1 = female
- All adult threshold comparisons include the 0.12 cm/kg adult rounding tolerance unless explicitly noted otherwise

---

## 1. User-configurable parameters

### 1A. Wrapper-level (exposed in `cleangrowth()` signature)

| Parameter | Default | Used in | Role |
|-----------|---------|---------|------|
| `id` | `NULL` (required) | `cleangrowth` preprocessing; output | Unique row identifier; preserved untouched into output. Must be provided. |
| `subjid` | (required, vector) | `cleangrowth` preprocessing | Subject identifier vector |
| `param` | (required, vector) | `cleangrowth` preprocessing | Measurement-type vector (`WEIGHTKG`, `HEIGHTCM`, `HEADCM`, `HEIGHTIN`, `WEIGHTLBS`, `LENGTHCM`) |
| `agedays` | (required, vector) | `cleangrowth` preprocessing | Age in days at measurement |
| `sex` | (required, vector) | `cleangrowth` preprocessing | Sex indicator (0/`"m"`/`"M"` = male; 1/`"f"`/`"F"` = female) |
| `measurement` | (required, vector) | `cleangrowth` preprocessing | Recorded value in units specified by `param` |
| `adult_cutpoint` | `20` | Preprocessing (pediatric/adult split) | Age (years) dividing pediatric/adult processing. Clamped to `[18, 20]`. |
| `sd.recenter` | `NA` | Recentering (Phase 11) | Optional user-supplied recentering medians. `NA` uses built-in `rcfile-2023-08-15_format.csv.gz`. |
| `ref.data.path` | `""` | Preprocessing; Child Step 17 | Override directory for reference files. `""` uses package `inst/extdata/`. |
| `batch_size` | `2000` | Outer-wrapper batching | Subjects per processing batch |
| `parallel` | `FALSE` | Inner batching | Run batches in parallel (plyr::ddply dispatch) |
| `num.batches` | `NA` | Inner batching | Number of parallel batches; auto-set when `parallel = TRUE` |
| `log.path` | `NA` | Parallel logging | Log directory for parallel mode |
| `ref_tables` | `NULL` | All reads | Pre-loaded reference closures from `gc_preload_refs()`; skips per-call disk reads (~0.93 sec/call) |
| `cached_results` | `NULL` | Partial-run mode | data.table from prior `cleangrowth()` call. Auto-detects changed subjects if `changed_subjids` is NULL. |
| `changed_subjids` | `NULL` | Partial-run mode | Optional explicit vector of subject IDs to re-process |
| `tri_exclude` | `FALSE` | Output assembly | Include `tri_exclude` column in output |
| `quietly` | `TRUE` | All steps | Suppress progress messages |

### 1B. Child algorithm (passed through `cleangrowth()` into `cleanchild()`)

| Parameter | Default | Used in | Role |
|-----------|---------|---------|------|
| `biv.z.wt.low.young` | `-25` | Child Step 7 | Lower unrecentered CSD z cutoff for WEIGHTKG at `ageyears < 1` |
| `biv.z.wt.low.old` | `-15` | Child Step 7 | Lower unrecentered CSD z cutoff for WEIGHTKG at `ageyears >= 1` |
| `biv.z.wt.high` | `22` | Child Step 7 | Upper unrecentered CSD z cutoff for WEIGHTKG (all ages) |
| `biv.z.ht.low.young` | `-25` | Child Step 7 | Lower unrecentered CSD z cutoff for HEIGHTCM at `ageyears < 1` |
| `biv.z.ht.low.old` | `-15` | Child Step 7 | Lower unrecentered CSD z cutoff for HEIGHTCM at `ageyears >= 1` |
| `biv.z.ht.high` | `8` | Child Step 7 | Upper unrecentered CSD z cutoff for HEIGHTCM (all ages) |
| `biv.z.hc.low` | `-15` | Child Step 7 | Lower unrecentered CSD z cutoff for HEADCM (all ages) |
| `biv.z.hc.high` | `15` | Child Step 7 | Upper unrecentered CSD z cutoff for HEADCM (all ages) |
| `cf_rescue` | `"standard"` | Child Step 6 | CF rescue mode: `"standard"` (lookup), `"none"` (all CFs excluded), `"all"` (all rescued). Validated via `match.arg()`. |
| `cf_detail` | `FALSE` | Child Step 6; output | If TRUE, add diagnostic columns `cf_status` and `cf_deltaZ` to output |
| `ewma_window` | `15` | Child Steps 11, 13, 15/16, 17; Adult Steps 9Wb, 11Wb | Max observations on each side that contribute to EWMA weighting. Single parameter in `cleangrowth()` affecting both child and adult. |
| `error.load.mincount` | `2` | Child Step 21 | Minimum "real" errors (not SDE / CF / Missing / Not-Cleaned) before error-load is evaluated |
| `error.load.threshold` | `0.5` | Child Step 21 | Error ratio above which all remaining Include rows for the subject-param are excluded. Strict `>`. |
| `length.adjust` | `FALSE` | Preprocessing (before z-scores) | If TRUE, subtracts 0.7 cm from LENGTHCM measurements with `agedays > 2 × 365.25` before z-score calculation. Converts post-infancy recumbent-labeled lengths to the standing equivalent assumed by WHO/CDC references at those ages. No effect for `agedays ≤ 730`. `child_clean.R` → `cleangrowth()`, inside the outer batch loop before `param == 'LENGTHCM'` relabeling. |

### 1C. Adult algorithm (wrapper-exposed)

Three adult parameters are exposed directly through `cleangrowth()`; the rest are set by the permissiveness preset (Section 1D).

| Parameter | Default | Used in | Role |
|-----------|---------|---------|------|
| `adult_permissiveness` | `"looser"` | `cleanadult` preset resolver | Preset level: `"loosest"`, `"looser"`, `"tighter"`, `"tightest"`. Sets defaults for all adult sub-parameters. |
| `adult_scale_max_lbs` | `Inf` | Adult Step 4W | Physical scale upper limit in pounds. Weights at this cap (± rounding) get `Exclude-A-Scale-Max`. |
| `ewma_window` | `15` | Adult Steps 9Wb, 11Wb | See Section 1B — shared parameter with child; single value in `cleangrowth()` affects both algorithms. |

### 1D. Adult permissiveness presets (set by `adult_permissiveness`)

All adult sub-parameters below are set by the permissiveness preset and can be individually overridden by passing them directly to `cleanadult()`. `cleangrowth()` currently exposes only `adult_permissiveness` and `adult_scale_max_lbs`, not individual adult sub-parameters.

| Parameter | loosest | looser | tighter | tightest | Role |
|-----------|---------|--------|---------|----------|------|
| `overall_ht_min` (cm) | 50 | 120 | 142 | 147 | BIV lower HT limit (Adult Step 1) |
| `overall_ht_max` (cm) | 244 | 230 | 213 | 208 | BIV upper HT limit (Adult Step 1) |
| `overall_wt_min` (kg) | 20 | 30 | 36 | 39 | BIV lower WT limit (Adult Step 1) |
| `overall_wt_max` (kg) | 500 | 270 | 159 | 136 | BIV upper WT limit (Adult Step 1) |
| `overall_bmi_min` | 5 | 12 | 16 | 18 | BIV lower BMI limit (Adult Step 1) |
| `overall_bmi_max` | 300 | 65 | 45 | 40 | BIV upper BMI limit (Adult Step 1) |
| `single_ht_min_bmi` (cm) | 60 | 120 | 142 | 147 | 1D HT lower limit when BMI available (Adult Step 13) |
| `single_ht_max_bmi` (cm) | 245 | 230 | 213 | 208 | 1D HT upper limit when BMI available |
| `single_wt_min_bmi` (kg) | 12 | 30 | 36 | 39 | 1D WT lower limit when BMI available |
| `single_wt_max_bmi` (kg) | 350 | 270 | 159 | 136 | 1D WT upper limit when BMI available |
| `single_ht_min_nobmi` (cm) | 122 | 120 | 142 | 147 | 1D HT lower limit when BMI not available |
| `single_ht_max_nobmi` (cm) | 245 | 230 | 213 | 208 | 1D HT upper limit when BMI not available |
| `single_wt_min_nobmi` (kg) | 30 | 30 | 36 | 39 | 1D WT lower limit when BMI not available |
| `single_wt_max_nobmi` (kg) | 350 | 270 | 159 | 136 | 1D WT upper limit when BMI not available |
| `single_bmi_min` | 10 | 12 | 16 | 18 | 1D BMI lower limit |
| `single_bmi_max` | 250 | 65 | 45 | 40 | 1D BMI upper limit |
| `wtallow_formula` | PW-H | PW-H | PW-L | allofus15 | Weight-allowance formula family (see `wtallow-formulas.md`). Affects Adult Steps 9Wa, 9Wb, 11Wa, 11Wa2, 11Wb. |
| `perclimit_low` (wt ≤45 kg) | 0.5 | 0.5 | 0.7 | 0.7 | Percentage limit (Adult Steps 11Wa, 11Wb) |
| `perclimit_mid` (wt 45–80 kg) | 0.4 | 0.4 | 0.4 | 0.4 | Percentage limit |
| `perclimit_high` (wt >80 kg) | 0 (disabled) | 0 (disabled) | 0.4 | 0.4 | Percentage limit |
| `mod_ewma_f` | 0.75 | 0.75 | 0.60 | 0.60 | Moderate EWMA directional factor (Adult Step 11Wb) |
| `error_load_threshold` | 0.41 | 0.41 | 0.29 | 0.29 | Error ratio above which all remaining Includes excluded (Adult Step 14). Strict `>`. |
| `ht_band` (inches) | 3 | 3 | 2 | 2 | Height band tolerance (Adult Step 10H). Converted to cm via `× 2.54 + 0.12`. |
| `allow_ht_loss` | TRUE | FALSE | FALSE | FALSE | Whether height loss over time is allowed (Adult Step 10H) |
| `allow_ht_gain` | TRUE | TRUE | TRUE | FALSE | Whether height gain over time is allowed (Adult Step 10H) |
| `repval_handling` | independent | independent | linked | linked | Repeated-value handling mode (Adult Steps 2W, 9Wb, 11Wb) |
| `max_rounds` | 100 | 100 | 100 | 100 | Maximum iterations for Adult Step 11Wb moderate EWMA loop |

### 1E. Removed parameters

No currently supported deprecated parameters. `include.carryforward` and `weight_cap` were removed from the `cleangrowth()` signature. Callers using pre-v3.0.0 code should migrate to `cf_rescue` (replaces `include.carryforward`) and `adult_scale_max_lbs` (replaces `weight_cap`).

---

## 2. Hardcoded thresholds (child)

### 2A. Preprocessing, Z-score calculation (wrapper phases, but scoped to child path)

See Section 4 for preprocessing-level values that apply before dispatch.

### 2B. Early Child Step 13 (SDE-Identicals)

| Threshold | Value | Location | Role |
|-----------|-------|----------|------|
| Age-dependent tiebreaker (birth) | `agedays == 0` → keep lowest `internal_id` | `cleanchild`, Early Child Step 13 inline block | Prefer earliest measurement at birth (before postnatal fluid shifts) |
| Age-dependent tiebreaker (other) | keep highest `internal_id` | `cleanchild`, Early Child Step 13 | Prefer later measurement (likely more careful) |

### 2C. Child Step 5 (Temporary SDE)

| Threshold | Value | Location | Role |
|-----------|-------|----------|------|
| Fallback SP median is NA | use DOP median | `identify_temp_sde()` in `child_clean.R` | Use cross-parameter median when subject-parameter median unavailable |
| Fallback both medians NA | treat 0 as median (`absdmedian.spz = |tbc.sd|`) | `identify_temp_sde()` | Graceful degradation when no reference point |
| Fallback DOP median NA | `absdmedian.dopz = Inf` (sorts last) | `identify_temp_sde()` | Deprioritize rows with no DOP reference |

### 2D. Child Step 6 (Carried Forwards)

| Threshold | Value | Location | Role |
|-----------|-------|----------|------|
| CF detection | Exact equality on `v.orig`; prior unique ageday must have exactly one measurement | `cleanchild` Child Step 6 Phase 1 | Row is CF when `!is.na(prior_single_val) & v.orig == prior_single_val` |
| `wholehalfimp` WEIGHTKG | `abs((v.orig * 2.20462262) %% 1) < 0.01` | `cleanchild` Child Step 6 Phase 3a | Whole pounds detection |
| `wholehalfimp` HEIGHTCM / HEADCM | `abs((v.orig / 2.54) %% 0.5) < 0.01` | `cleanchild` Child Step 6 Phase 3a | Whole or half inches detection |
| Lookup rescue condition | threshold is non-NA, > 0, and `absdiff < threshold` | `cleanchild` Child Step 6 Phase 4; `.cf_get_thresholds()` | Standard-mode rescue |

CF rescue lookup thresholds — reference-data (Section 4).

### 2E. Child Step 7 (BIV) — absolute thresholds

All strict `<` or `>`. Standardized z cutoffs (configurable) are in Section 1B.

| Threshold | Value | Scope | Role |
|-----------|-------|-------|------|
| WEIGHTKG low, year 1 | `v < 0.2` kg | `agedays <= 365` | Published minimum viable birth weight |
| WEIGHTKG low, after year 1 | `v < 1` kg | `agedays > 365` | Minimum plausible outpatient weight |
| WEIGHTKG high at birth | `v > 10.5` kg | `agedays == 0` | Maximum plausible birth weight |
| WEIGHTKG high, under 2y | `v > 35` kg | `ageyears < 2` | Highest plausible weight under 2y |
| WEIGHTKG high, any age | `v > 600` kg | All ages | Published maximum human weight |
| HEIGHTCM low | `v < 18` cm | All ages | Fenton z = -6 at 22 0/7 weeks |
| HEIGHTCM high | `v > 244` cm | All ages | Published maximum |
| HEIGHTCM high at birth | `v > 65` cm | `agedays == 0` | Fenton z = 6 at 40 0/7 weeks |
| HEADCM low | `v < 13` cm | All ages | Fenton z = -6 at 22 0/7 weeks |
| HEADCM high | `v > 75` cm | All ages | Published maximum |
| HEADCM high at birth | `v > 50` cm | `agedays == 0` | Birth maximum |
| Standardized-BIV age boundary | `ageyears < 1` vs `ageyears >= 1` | Child Step 7 standardized block | Splits WT/HT lower cutoffs into "young" and "old" bands |
| Guard pattern | `!grepl("^Exclude-C-BIV$", exclude)` | Child Step 7 standardized block | Skip rows just assigned `Exclude-C-BIV` by absolute block |

### 2F. Child Step 9 (Evil Twins)

| Threshold | Value | Location | Role |
|-----------|-------|----------|------|
| OTL threshold | both `|tbc.sd_next - tbc.sd_current| > 5` AND `|ctbc.sd_next - ctbc.sd_current| > 5` (strict, dual) | `calc_otl_evil_twins()` in `child_clean.R` | Adjacent extreme pair |
| Pre-filter count | `sp_count_9 > 2L` | `cleanchild` Child Step 9a | Requires 3+ total rows per subject-param |
| Minimum valid rows | `< 2L` breaks group loop | `cleanchild` Child Step 9b | Stop iteration when too few rows remain |
| Tiebreaker | lowest `internal_id` | `cleanchild` Child Step 9b | Deterministic (not age-dependent) |

### 2G. Child Step 11 (EWMA1 — Extreme EWMA)

| Threshold | Value | Location | Role |
|-----------|-------|----------|------|
| Count pre-filter | skip subject-params with `<= 2` Include values | `cleanchild` Child Step 11a | EWMA needs 3+ points |
| Extreme pre-filter | skip if no Include has `|tbc.sd| > 3.5` | `cleanchild` Child Step 11a | Cheap bypass for subject-params that can't meet criterion |
| EWMA exponent, gap ≤ 1y | `-1.5` | `cleanchild` Child Step 11b (age-gap-dependent exponent) | Less negative → nearby neighbors carry more weight |
| EWMA exponent, gap ≥ 3y | `-3.5` | `cleanchild` Child Step 11b | More negative → only closest values carry weight |
| EWMA exponent interpolation | `-1.5 - (ageyears - 1)` for gap between 1 and 3y | `cleanchild` Child Step 11b | Linear interpolation |
| Positive outlier `dewma.all` | `> 3.5` | `cleanchild` Child Step 11b pot_excl | Strict |
| Positive outlier `dewma.before` | `> 3` | Same | Strict |
| Positive outlier `dewma.after` | `> 3` | Same | Strict |
| Positive outlier `tbc.sd` | `> 3.5` | Same | Strict |
| Positive outlier corrected | `c.dewma.all > 3.5` OR `ctbc.sd` is NA | Same | NA escape auto-satisfies for non-potcorr |
| Negative outlier mirrors | `dewma.all < -3.5`, `dewma.before < -3`, `dewma.after < -3`, `tbc.sd < -3.5`, `c.dewma.all < -3.5` | Same | Mirrors positive side |
| Endpoint rule | First and last Include values never excluded | `cleanchild` Child Step 11b | EWMA unreliable at endpoints |
| Per-iteration exclusion cap | At most 1 per subject-param | `cleanchild` Child Step 11b | Prevents cascade |
| Worst-value sort key | `abs(tbc.sd + dewma.all)` descending; tiebreaker lowest `internal_id` | `cleanchild` Child Step 11b | Combines z-score magnitude and EWMA deviation |

### 2H. Child Step 13 (Final SDE Resolution)

| Threshold | Value | Location | Role |
|-----------|-------|----------|------|
| One-day SDE-All-Extreme | `absdmedian.spz > 2` (strict) | `cleanchild` Child Step 13 Phase B2 | If closest to median is > 2 SD away, all same-day values excluded |
| EWMA SDE-All-Extreme | `min_absdewma > 1` (strict) | `cleanchild` Child Step 13 Phase B3 | If closest EWMA-wise is > 1 SD away, all same-day values excluded |
| Age-dependent tiebreaker | Same as Early Child Step 13 | `identify_temp_sde()` | Birth: lowest; other ages: highest |

### 2I. Child Steps 15/16 (EWMA2 — Moderate EWMA)

Perturbation magnitudes:

| Parameter | Value | Scope |
|-----------|-------|-------|
| WEIGHTKG `p_plus` | `1.05 × v` | Child Step 15 pre-calc |
| WEIGHTKG `p_minus` | `0.95 × v` | Child Step 15 pre-calc |
| HEIGHTCM `p_plus` | `v + 1` cm | Child Step 15 pre-calc |
| HEIGHTCM `p_minus` | `v - 1` cm | Child Step 15 pre-calc |
| HEADCM `p_plus` | `v + 1` cm | Child Step 15 pre-calc |
| HEADCM `p_minus` | `v - 1` cm | Child Step 15 pre-calc |

`addcrit` thresholds: all `> 1` (high) or `< -1` (low), strict.

Child Step 15 exclusion rules (all thresholds strict, hardcoded):

| Rule | Position | Gap (days) | `dewma.all` | `c.dewma.all` | Extra condition |
|------|----------|------------|-------------|---------------|-----------------|
| middle | Not first/last | — | `> 1` / `< -1` | `> 1` / `< -1` | addcrit |
| birth-WT | `agedays == 0` | next `< 365.25` | `> 3` / `< -3` | `> 3` / `< -3` | addcrit |
| birth-WT-ext | `agedays == 0` | next `>= 365.25` | `> 4` / `< -4` | `> 4` / `< -4` | addcrit |
| first | `first_meas` | next `< 365.25` | `> 2` / `< -2` | `> 2` / `< -2` | addcrit |
| first-ext | `first_meas` | next `>= 365.25` | `> 3` / `< -3` | `> 3` / `< -3` | addcrit |
| last | Last row | prev `< 2*365.25` | `> 2` / `< -2` | `> 2` / `< -2` | `abs(tbc_prev) < 2` |
| last-high | Last row | prev `< 2*365.25` | `> abs(tbc_prev)` / `< -abs(tbc_prev)` | **`> 3` / `< -3`** | `abs(tbc_prev) >= 2` |
| last-ext | Last row | prev `>= 2*365.25` | `> 3` / `< -3` | `> 3` / `< -3` | `abs(tbc_prev) < 2`; DOP diff |
| last-ext-high | Last row | prev `>= 2*365.25` | `> 1 + abs(tbc_prev)` / `< -1 - abs(tbc_prev)` | **`> 3` / `< -3`** | `abs(tbc_prev) >= 2`; DOP diff |

DOP diff thresholds (last-ext, last-ext-high): `tbc.sd - tbc_dop > 4` or `< -4`.

Child Step 16 (Birth HT/HC — separate step):

| Rule | Gap to next | `dewma.all` | `c.dewma.all` |
|------|-------------|-------------|---------------|
| birth-HT-HC | `< 365.25` | `> 3` / `< -3` | `> 3` / `< -3` |
| birth-HT-HC-ext | `>= 365.25` | `> 4` / `< -4` | `> 4` / `< -4` |

### 2J. Child Step 17 (Height/HC Velocity)

Velocity reference tables are reference-data (Section 4). Tolerance formulas below are hardcoded.

**HEIGHTCM:**

| Item | Value | Role |
|------|-------|------|
| Tanner age floor | midpoint age `≥ 30 months` | Only apply Tanner above this |
| max.ht.vel floor (always) | `2.54 cm` (1 inch) | Minimum allowable increase |
| max.ht.vel floor (gap > 2 months) | `5.08 cm` | |
| max.ht.vel floor (gap > 6 months) | `10.16 cm` | |
| max.ht.vel floor (gap > 1 year) | `20.32 cm` | |
| mindiff (Tanner, gap < 1y) | `0.5 × min.ht.vel × (gap_years)^2 - 3` | |
| mindiff (Tanner, gap > 1y) | `0.5 × min.ht.vel - 3` | |
| maxdiff (Tanner, gap < 1y) | `2 × max.ht.vel × (gap_years)^0.33 + 5.5` | |
| maxdiff (Tanner, gap > 1y) | `2 × max.ht.vel × (gap_years)^1.5 + 5.5` | |
| WHO age cap | `round(agedays / 30.4375)` capped at 24 | |
| WHO gap→increment map | `<20d, 20–45, 46–75, 76–106, 107–152, 153–198, ≥200` → 1/1/2/3/4/6/6 months | |
| WHO mindiff (<9m) | `who_mindiff × 0.5 - 3` | |
| WHO maxdiff (<9m) | `who_maxdiff × 2 + 3` | |
| WHO fallback (≥9m, no Tanner) | Same transform | |
| Default mindiff (no ref) | `-3` cm | |
| Birth adjustment (`agedays == 0`) | `mindiff -= 1.5`, `maxdiff += 1.5` | |

**HEADCM:**

| Item | Value | Role |
|------|-------|------|
| WHO age groups | Same as HT gap→increment map | |
| mindiff transform | `who_mindiff × 0.5 - 1.5` | Tighter than HT |
| maxdiff transform | `who_maxdiff × 2 + 1.5` | Tighter than HT |
| Default mindiff (no ref) | `-1.5` cm | |
| Birth adjustment (`agedays == 0`) | `mindiff -= 0.5`, `maxdiff += 0.5` | Tighter than HT |

**Violation condition:** `diff_prev < mindiff_prior` OR `diff_next < mindiff` OR `diff_prev > maxdiff_prior` OR `diff_next > maxdiff`. Strict inequalities.

**Resolution with 3+ measurements:** EWMA-based tiebreaker (`bef.g.aftm1`, `aft.g.aftm1`), internal_id tiebreaker.
**Resolution with 2 measurements:** Higher `abs(tbc.sd)`, internal_id tiebreaker.

### 2K. Child Step 19 (Pairs and Singles)

| Threshold | Value | Role |
|-----------|-------|------|
| Pair, gap ≥ 365.25 | `|diff_tbc| > 4` AND `|diff_ctbc| > 4` (or NA) | Large-gap pair |
| Pair, gap < 365.25 | `|diff_tbc| > 2.5` AND `|diff_ctbc| > 2.5` (or NA) | Small-gap pair |
| Single with DOP | `|tbc.sd| > 3` AND `comp_diff > 5` | Single-measurement exclusion with cross-check |
| Single without DOP | `|tbc.sd| > 5` | Very extreme isolated measurement |
| Gap boundary | `365.25` days (inclusive ≥ vs strict <) | |

### 2L. Child Step 21 (Error Load)

Configurable thresholds in Section 1B. Hardcoded component:

| Item | Value |
|------|-------|
| `non_error_codes` list | `Exclude-C-Identical`, `Exclude-C-Extraneous`, `Exclude-C-CF`, `Exclude-Missing`, `Exclude-Not-Cleaned` |

Note: `Exclude-C-Temp-Same-Day` is NOT in `non_error_codes`. By the time Child Step 21 runs, Child Step 13 has resolved every temp SDE; any `Exclude-C-Temp-Same-Day` row at this point would be a bug and would correctly count as an error.

---

## 3. Hardcoded thresholds (adult)

All adult threshold comparisons apply the 0.12 cm/kg rounding tolerance unless explicitly noted otherwise.

### 3A. Global

| Threshold | Value | Location | Role |
|-----------|-------|----------|------|
| Adult rounding tolerance | `0.12` cm/kg | All adult comparisons | Added to threshold comparisons (not BMI) |
| Adult EWMA exponent | `-5` | All adult EWMA steps | More negative than child's `-1.5` |
| Adult EWMA anchor | `0` | All adult EWMA steps | |
| BMI formula | `weight_kg / (height_cm / 100)²` | Multiple | |

### 3B. Adult Step 1 (BIV)

Configurable HT/WT/BMI limits in Section 1D. Guard:

| Item | Value | Role |
|------|-------|------|
| BMI BIV guard | Runs only when `overall_bmi_min < single_bmi_min` OR `overall_bmi_max > single_bmi_max` | Skip BMI-BIV when 1D limits are already tighter |
| BMI BIV rounding tolerance | Not applied | Unlike other comparisons, BMI BIV uses exact values |

### 3C. Adult Step 4W (Scale Max)

| Threshold | Value | Role |
|-----------|-------|------|
| Scale cap conversion | `round(adult_scale_max_lbs / 2.2046226, 1) ± 0.1` kg | Convert lbs to kg with detection tolerance |
| `dewma50p` threshold | `dewma_all > 50 + 0.12` | Large positive deviation |
| `dewma50p` before/after | `> 45 + 0.12` (0.9 × 50) | Directional check |
| `dewma50m` mirror | `< -(50 + 0.12)`; before/after `< -(45 + 0.12)` | |
| d50p adjacency | `d_prev > 50 + 0.12` AND `d_next > 50 + 0.12` | |
| Missing adjacency neighbors | Treated as +Inf (positive) / -Inf (negative) | Confirms exclusion at edges |
| Minimum included | 2+ | |

### 3D. Adult Step 9Wa (Evil Twins)

| Threshold | Value | Role |
|-----------|-------|------|
| Interval boundary | `≤6m` vs `>6m` | Two ET-cap tiers |
| Minimum values | 3+ | Pairs guard |
| Low OOB score penalty | 99999 for values `<38 kg` | Guardrail |
| High OOB score penalty | 99999 for values `>180 kg` | Guardrail |
| Range gate | `max - min ≤ min ET cap + 0.12` | Evaluation gate |

### 3E. Adult Step 9H (Final SDE Height)

| Threshold | Value | Role |
|-----------|-------|------|
| Category 2 | `daystot >= 4` AND `sderatio < 0.5` | |
| Category 3 | `daystot == 2` or `3` OR `sderatio >= 0.5` | |
| Identical tiebreaker | lowest `internal_id` | |
| Non-identical tiebreaker | highest `internal_id` | Differs from Category 2 identical |

### 3F. Adult Step 9Wb (Extreme EWMA)

| Threshold | Value | Role |
|-----------|-------|------|
| Positive outlier (all) | `dewma_all > threshold + 0.12` | Main criterion |
| Positive outlier (before/after) | `dewma_before > 0.9 × threshold + 0.12` AND `dewma_after > 0.9 × threshold + 0.12` | Directional |
| Negative outlier | Mirrors of positive with `<` | |
| Directional factor | `0.9` | Before/after relaxation |
| Missing directional | Treated as Inf (confirms exclusion) | |
| Minimum non-extraneous | 3+ | |
| Threshold source | `compute_et_limit(min_gap_months, formula, uw)` where `uw = max(wt, ewma_all)` | |
| Round limit | None (R implementation) | |

### 3G. Adult Step 10H (Height Distinct)

| Threshold | Value | Role |
|-----------|-------|------|
| Effective band | `ht_band × 2.54 + 0.12` cm | Converts inches to cm |
| Gain rescue age boundary | `ageyears1 < 25` | Below this age, allow height gain |
| htallow formula | `velocity × (ln(age2 - 16.9) - ln(age1 - 16.9))` | |
| htallow coefficient (gap <2yr) | 15.5 | 2D gain |
| htallow coefficient (gap 2–3yr inclusive) | 13 | 2D gain |
| htallow coefficient (gap >3yr) | 12 | 2D gain |
| Age cap for htallow | 25 | |
| Gain rescue tolerance | `(htallow + 2) × 2.54 + 0.12` cm | 2-inch buffer |
| Loss rescue tolerance (age <50) | `5 × 2.54 + 0.12` cm | 5 inches |
| Loss rescue tolerance (age ≥50) | `7 × 2.54 + 0.12` cm | 7 inches |
| Frequency rescue ratio | `count ≥ count × 4/3` | |
| w2 window step | `ht_band × 2.54 + 0.12` cm | |
| w2/o2 ratio minimum | `≥ 3/2` | |
| w2 selection score | `tot_w2 + 0.5 × numdistinct_w2` | |
| Loss group band | `2 inches = 5.08 + 0.12 cm` | |
| Loss group max count | 3 | |
| Gain group max count | 6 | |
| Loss group cumulative limit | `6/8/9 inches` depending on age pattern | |
| Gain group htallow (<1yr) | 20 | |
| Gain group htallow (1–3yr) | 15 | |
| Gain group htallow (>3yr) | 12 | |
| Loss magnitude 3D (age <50) | `≤ 5 inches` | |
| Loss magnitude 3D (age ≥50) | `≤ 7 inches` | |

### 3H. Adult Step 11Wa (2D Ordered Weight)

| Threshold | Value | Role |
|-----------|-------|------|
| Rounding tolerance | `0.12` kg | |
| wtallow source | `compute_wtallow(months, formula, uw)` (see `wtallow-formulas.md`) | |

### 3I. Adult Step 11Wa2 (2D Non-Ordered Weight)

| Threshold | Value | Role |
|-----------|-------|------|
| Rounding tolerance | `0.12` kg | |
| Dominance threshold | `dominant_pct > 0.65` (strict) | |
| Interval conversion | `abs(age_days[i+1] - age_days[i]) / 30.4375` | Days to months |
| Prior non-SDE exclusion set | `Include`, temp extraneous, `Exclude-A-Identical`, `Exclude-A-Extraneous`, empty, NA | |

### 3J. Adult Step 11Wb (Moderate EWMA)

| Threshold | Value | Role |
|-----------|-------|------|
| Minimum distinct values | 3+ | |
| Minimum included | 3+ | |
| Standard pathway main | `|dewma_all| > wtallow + 0.12` | |
| Standard pathway directional | Both directional `> mod_ewma_f × wtallow + 0.12` (same direction) | |
| Trajectory rescue error margin | `±5 kg` | |
| Trajectory rescue rounding | Nearest 0.2 kg | |
| Extrapolation distance | `≤ 2× source interval` | |
| Run-length thresholds | 4+ consecutive → error load; 1 → isolated; 2–3 → pair/trio | |
| Unreliable neighbor threshold | `≤ 14 days` AND `> wtallow apart` | |
| Alternate pathway criterion | `|dewma_all| > wtallow + 0.12` AND `|dewma_aft| > 0.75 × wtallow + 0.12` (or `|dewma_bef|`) | |
| Alternate pathway factor | `0.75` | Relaxed directional factor |
| Percentage criterion | `wt/ewma_all < perclimit` AND `wt/ewma_bef < perclimit` AND `wt/ewma_aft < perclimit` | All three must pass |
| Graduated scoring multiplier | `min(1.0, 0.6 + 0.4 × min_excess)` | |
| Edge scoring formula | `pmax(0, |dewma_all| - wta) / wta_base` | |
| `wta_base` | wtallow at `UW = 120` (no adjustment) | |

### 3K. Adult Step 13 (Single Distinct 1D)

All limits configurable (Section 1D). No additional hardcoded thresholds besides BMI formula.

### 3L. Adult Step 14 (Error Load)

| Threshold | Value | Role |
|-----------|-------|------|
| Minimum denominator | `≥ 3` | |
| SDE excluded from counts | `Exclude-A-Identical`, `Exclude-A-Extraneous` | |
| RV-propagated codes | Excluded from error count | |

---

## 4. Reference-data-driven values

### 4A. Growth reference tables (`inst/extdata/`)

| File | Loaded by | Consumed by | User-overridable? |
|------|-----------|-------------|-------------------|
| `growthfile_who.csv.gz` | `read_anthro()`, `gc_preload_refs()` | Z-score calculation (wrapper preprocessing) | Via `ref.data.path` |
| `growthfile_cdc_ext_infants.csv.gz` | `read_anthro()`, `gc_preload_refs()` | Z-score calculation | Via `ref.data.path` |
| `fent_foraga.csv.gz` | `cleangrowth()` Phase 10 (GA correction) | Fenton weight-to-GA derivation for potcorr | Via `ref.data.path` |
| `fenton2025_ms_lookup_smoothed.csv` | `cleangrowth()` Phase 10 | Fenton M/S lookup for corrected z-scores | Via `ref.data.path` |
| `rcfile-2023-08-15_format.csv.gz` | Recentering (Phase 11) | `sd.median` lookup for recentering | Via `sd.recenter` parameter |
| Tanner height velocity table | `cleanchild` Child Step 17 | HT velocity tolerance | Via `ref.data.path` |
| WHO HT velocity table | `cleanchild` Child Step 17 | HT velocity (<9m) | Via `ref.data.path` |
| WHO HC velocity table | `cleanchild` Child Step 17 | HC velocity | Via `ref.data.path` |

### 4B. Lookup tables defined in code

| Lookup | Location | Consumed by | Overridable? |
|--------|----------|-------------|--------------|
| CF rescue thresholds (HT/WT, other/imperial × 8 age bins × 5 interval bins) | `.cf_rescue_lookup()` in `child_clean.R` | Child Step 6 | Not user-overridable (edit code) |
| CF rescue age bins | 0–3 mo, 3–6 mo, 6–12 mo, 1–2 y, 2–5 y, 5–10 y, 10–15 y, 15–20 y | Child Step 6 | |
| CF rescue interval bins | <1 wk, 1 wk – 1 mo, 1–6 mo, 6 mo – 1 yr, >1 yr | Child Step 6 | |
| CF rescue threshold levels | 0.05 / 0.20 / 0.40 / NR (no rescue) / `--` (impossible cell) | Child Step 6 | |
| HEADCM CF rescue thresholds | **Placeholder: reuses HEIGHTCM matrices** (see Known Issues D18) | Child Step 6 | |
| DOP mapping | WT↔HT, HC→HT | `get_dop()` in `child_clean.R` | Not overridable |
| Adult permissiveness presets | Section 1D values | `resolve_permissiveness()` in `adult_support.R` | Override individual params directly to `cleanadult()` |
| Adult wtallow formulas (PW-H, PW-L, allofus15) | `compute_wtallow()` / `compute_et_limit()` in `adult_support.R`; full spec in `wtallow-formulas.md` | Adult Steps 9Wa, 9Wb, 11Wa, 11Wa2, 11Wb | Custom CSV also supported |

Full CF rescue threshold tables: see [`cf-rescue-thresholds.md`](cf-rescue-thresholds.md).
Full wtallow formula spec: see [`wtallow-formulas.md`](wtallow-formulas.md).

### 4C. Z-score / age-blending constants

| Item | Value | Location |
|------|-------|----------|
| WHO-only range (HT/WT) | `< 2 years` | Z-score infrastructure |
| WHO/CDC blend range (HT/WT) | `2–5 years`; formula `(cdc × (age-2) + who × (5-age)) / 3` | Z-score infrastructure |
| CDC-only range (HT/WT) | `> 5 years` | Z-score infrastructure |
| HEADCM age policy | WHO only at all ages through 5y | Z-score infrastructure |
| HEADCM age cutoff | `agedays > 3 × 365.25` → `Exclude-Not-Cleaned`; `agedays >= 5 × 365.25` → `Exclude-Not-Cleaned` | Preprocessing + Phase 11 |
| Imperial → metric (HEIGHTIN) | `× 2.54` | Preprocessing Phase 7 |
| Imperial → metric (WEIGHTLBS) | `÷ 2.2046226` | Preprocessing Phase 7 |
| Sex encoding | `0`/`"m"`/`"M"` → 0L; `1`/`"f"`/`"F"` → 1L | Preprocessing Phase 3 |
| `measurement == 0 → NaN` | Hardcoded conversion | Preprocessing Phase 3 |

### 4D. Gestational Age Correction (Child Step 2b / Phase 10)

| Item | Value | Role |
|------|-------|------|
| potcorr_wt condition | First weight for subject AND `sd.orig < -2` AND `age < 10 months` | Flag potentially preterm subject |
| `intwt` formula | `trunc(v × 100) × 10` (grams) | Integer weight for Fenton merge |
| `intwt` floor | Values in `[100, 500)` → `500` | Clamp low values to Fenton table minimum |
| `pmagedays` | `agedays + fengadays` | Postmenstrual age in days |
| `cagedays` | `pmagedays - 280` | Corrected age (days) |
| Fenton term cutoff | `fengadays < 259` | Below this → considered preterm for correction |
| Fenton weight scaling | `v_fenton = v × 1000` for WT; `v` otherwise | Fenton table is in grams |
| Fenton CSD z-score | `(v_fenton − M) / (S_upper × M)` if `v_fenton >= M`, else `/ (S_lower × M)` | Split-normal CSD |
| WHO supine/standing offset | `+0.8 cm` | Age > 2 but corrected age ≤ 2 |
| CDC supine/standing offset | `+0.7 cm` | Age > 2 but corrected age ≤ 2 |
| Fenton precedence ≤ 2y | `sd.c := unmod_zscore` when `ageyears_2b ≤ 2` | |
| Smoothing range | 2–4 years; formula `(sd.orig × (4 - age) + sd.c × (age - 2)) / 2` | |
| Above-4y policy | `sd.corr = sd.orig` | |
| uncorr revert rule | Revert correction if `sd.corr_abssumdiff > sd.orig_abssumdiff` | Correction worsens trajectory |
| uncorr window | First 4 weights AND `ageyears_2b < 2` | Where revert check runs |

---

## 5. Discrepancies (resolved 2026-04-18)

All 11 discrepancies flagged in the initial build have been resolved. Full record below; resolutions also applied to child-algorithm-reference, gc CLAUDE.md, and wrapper-narrative.

### 5.1. `include.carryforward` — removed from documentation

- **Resolution:** Removed from child-algorithm-reference.md summary and Child Step 6 per-step tables; removed from gc CLAUDE.md child Configurable Parameters table. The parameter is no longer in the `cleangrowth()` signature (post-D27), so it is no longer referenced as a supported parameter anywhere in user-facing documentation. Section 1E notes the removal for migration reference.

### 5.2. `weight_cap` — breadcrumb dropped

- **Resolution:** Removed "(formerly `weight_cap`)" parenthetical from gc CLAUDE.md adult Configurable Parameters table. Section 1E notes the removal.

### 5.3. `ewma_window` scope — all steps listed

- **Resolution:** Verified in code that `ewma_window` is used in Child Steps 11, 13, 15, 16, and 17 (5 sites total: inline ewma() calls and ewma_cache_init() calls). Updated the "Used in" column to "Child Steps 11, 13, 15/16, 17" in: parameters-reference.md Section 1B; child-algorithm-reference.md summary and Steps 11, 15/16 per-step tables; gc CLAUDE.md child Configurable Parameters table.

### 5.4. `sd.recenter` schema — documented in wrapper reference

- **Resolution:** Added a "Required schema for a custom `sd.recenter` data.table" sub-section to wrapper-narrative-2026-04-17.md in the Z-Score Infrastructure / Recentering block. Required columns: `param`, `sex`, `agedays`, `sd.median`. gc CLAUDE.md updated to mention the column list inline.

### 5.5. HEADCM CF rescue thresholds — no action

- **Status:** N/A. Already tracked as deferred item D18. Master table Section 4B describes the HT-matrix placeholder; no documentation change needed.

### 5.6. Child Step 21 error-load `non_error_codes` list — verified, references corrected

- **Code actual list:** `Exclude-C-Identical`, `Exclude-C-Extraneous`, `Exclude-C-CF`, `Exclude-Missing`, `Exclude-Not-Cleaned`.
- **Child reference previously said:** `Exclude-C-Identical`, `Exclude-C-Extraneous`, `Exclude-C-CF`, `Exclude-C-Temp-Same-Day`, `Exclude-Missing` — incorrect on two counts (missing `Exclude-Not-Cleaned`; added `Exclude-C-Temp-Same-Day` that is not in the code list).
- **Resolution:** Updated parameters-reference.md Section 2L and child-algorithm-reference.md Child Step 21 to match code. Added a note explaining why `Exclude-C-Temp-Same-Day` is intentionally not in the list (Child Step 13 resolves all temp SDEs before this step; a remaining one would indicate a bug and should count as an error).

### 5.7. Adult `adult_scale_max_lbs` kg conversion — no action

- **Status:** Informational only. No change needed. Flagged to make this hardcoded conversion discoverable.

### 5.8. Adult `mod_ewma_f` — CLAUDE.md cross-references master table

- **Resolution:** Added a cross-reference note beneath the adult permissiveness preset table in gc CLAUDE.md pointing to `parameters-reference.md` for the role and code location of each preset column (`mod_ewma_f`, `perclimit_*`, `error_load_threshold`, etc.). The numeric preset table remains in CLAUDE.md for quick reference.

### 5.9. `max_rounds` — added to CLAUDE.md preset table

- **Resolution:** Added `max_rounds` row (100 across all levels) to the gc CLAUDE.md Adult Permissiveness Presets table.

### 5.10. `wfl_ht_range` — no action (pipeline-level, correctly excluded)

- **Status:** Noted in this record. `wfl_ht_range` is a pipeline-level parameter used by `compute_outcomes()` in the Quality-Pipeline (Qual-AD) error-impact package, not by growthcleanr. Correctly excluded from this master table.

### 5.11. Adult sections — point to `wtallow-formulas.md`

- **Resolution:** Section 4B already points to `wtallow-formulas.md` for the full formula spec. Sections 3G–3J are left as-is (high-level threshold summary pointing out the formula families and key coefficients; `wtallow-formulas.md` is the authoritative source for every formula constant and the derivation logic).

---

## 6. Change log

| Date | Change |
|------|--------|
| 2026-04-18 | Initial version. Built from a code scan, three references scans, and two CLAUDE.md scans. |
| 2026-04-18 | All 11 initial discrepancies resolved. Updates to child-algorithm-reference, gc CLAUDE.md, and wrapper-narrative applied accordingly. Full record in Section 5. |
| 2026-04-20 | `ewma_window` scope extended to adult (Steps 9Wb, 11Wb); now passed through `cleangrowth()` to `cleanadult()`. Section 1B updated; Section 1C row added. |
