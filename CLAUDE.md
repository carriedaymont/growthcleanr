# CLAUDE.md — gc-github-latest (growthcleanr)

**Last updated:** 2026-04-12 (internal_id rename complete in both algorithms; error load bug fix)

## Overview

growthcleanr (gc) is an R package that identifies and flags
implausible anthropometric measurements (height, weight, head
circumference) in electronic health record data. It does not
remove any values — each row gets an `exclude` code:
`"Include"` or an exclusion code naming the step and reason.

Repository: https://github.com/carriedaymont/growthcleanr
Branch: `efficiency-updates` (v3.0.0, not yet merged to main)
CRAN: v2.2.0 (outdated — do not use for development)

This CLAUDE.md is for **Claude Code agents working on gc
code**. For pipeline/Qual-AD context, see
`__Pipeline/CLAUDE.md`.

**CRITICAL:** The only valid local working copy of growthcleanr
is `__Pipeline/gc-github-latest/`. Stale copies exist in
`/Users/Shared/` — do NOT use them. Never edit files outside
`__Pipeline/` without Carrie's explicit written permission in
chat.

---

## Algorithms

| Algorithm | Status | Default path | Notes |
|-----------|--------|--------------|-------|
| Child | Active, primary | `child_clean.R` | Default pediatric path in v3.0.0 |
| Adult | **Closed pending validation** | `adult_clean.R` + `adult_support.R` | Permissiveness framework, 4 exclusion levels. Do not modify without checking with Carrie first. |
| Legacy pediatric | Deprecated | `pediatric_clean_legacy.R` | `use_legacy_algorithm = TRUE` |
| `adjustcarryforward()` | Deprecated | removed in v3.0.0 | CF adjustment utility; present on Main/CRAN, removed from `efficiency-updates` |

The child algorithm replaces the legacy pediatric algorithm.
The legacy pediatric algorithm is maintained for backward
compatibility (`use_legacy_algorithm = TRUE`).
`prelim_infants` is deprecated — it triggers a warning and
maps to `use_legacy_algorithm` (verified; tested in
test-child-regression.R and test-cleangrowth.R). The
user-facing `prelim_infants` parameter has no other effect;
internal `read_anthro(prelim_infants = TRUE)` calls are
hardcoded and independent of the user-facing parameter.

This CLAUDE.md covers both the child and adult algorithms.

---

## Code Structure

### Key files

| File | What it contains |
|------|------------------|
| `R/child_clean.R` | `cleangrowth()` entry point + exports (top of file); `gc_preload_refs()` (pre-loads reference closures for repeated calls); `cleanbatch_child()` main algorithm + support functions (`valid()`, `temporary_extraneous_infants()`, `calc_otl_evil_twins()`, `calc_and_recenter_z_scores()`, `ewma()`, `ewma_cache_init()`/`ewma_cache_update()`, `get_dop()`, `read_anthro()`) |
| `R/pediatric_clean_legacy.R` | Legacy pediatric algorithm (`use_legacy_algorithm = TRUE`); deprecated, retained for backward compatibility |
| ~~`R/growth.R`~~ | Removed in v3.0.0 — `cleangrowth()` moved into `child_clean.R` |
| `R/utils.R` | Shared utilities (does NOT contain `valid()` — see note below) |
| `R/adult_clean.R` | `cleanadult()` main algorithm — permissiveness framework, 14 steps |
| `R/adult_support.R` | Adult support functions — permissiveness presets, EWMA, BIV, height groups, evil twins, error load, etc. |
| `inst/extdata/` | Reference tables (growth charts, recentering file, velocity tables) |
| `tests/testthat/` | Test suite (see Testing section) |

### Preprocessing → main algorithm split

`cleangrowth()` (in `child_clean.R`) handles:
- Input validation, data.table construction
- Imperial → metric conversion (HEIGHTIN, WEIGHTLBS)
- LENGTHCM → HEIGHTCM reclassification (param only,
  no measurement adjustment — known limitation)
- Age cutpoint split (pediatric vs. adult)
- Z-score calculation (CSD method, WHO/CDC blending)
- GA correction (Step 2b) for potcorr subjects
- Recentering → tbc.sd and ctbc.sd
- Missing value identification
- Batching and dispatch to `cleanbatch_child()`

`cleanbatch_child()` processes one batch through all
cleaning steps (see Complete List of Steps below).

### Batching

Two layered batching systems:

1. **Outer wrapper** (memory management): Divides subjects
   into batches of 2,000. Batch size is hard-coded — should
   be parameterized (future work).

2. **Inner batching** (parallelism): Subdivides for parallel
   processing. With `parallel = FALSE` (current default),
   this is a no-op.

---

## Input/Output Format

### Required input (individual vectors, not a dataframe)

| Parameter | Type | Description |
|-----------|------|-------------|
| `subjid` | any | Subject identifier |
| `param` | character | `"WEIGHTKG"`, `"WEIGHTLBS"`, `"HEIGHTCM"`, `"HEIGHTIN"`, `"LENGTHCM"`, or `"HEADCM"` |
| `agedays` | numeric | Age in days at measurement |
| `sex` | any | `0`/`"m"`/`"M"` = male; `1`/`"f"`/`"F"` = female |
| `measurement` | numeric | Recorded value in units specified by `param` |
| `id` | any (required) | Unique row identifier; preserved in output. Can be numeric, character, UUID, etc. |

### Input handling

- `internal_id` (sequential `1:N`) is created at entry for all
  internal sorting and tiebreaking. The user's `id` is never
  used for algorithm logic — only for output.

- `measurement == 0` → replaced with `NaN` (treated as missing)
- Imperial converted to metric before processing
- LENGTHCM relabeled to HEIGHTCM (no measurement adjustment)
- Data split by `adult_cutpoint` (default 20 years)

### Sorting (critical for deterministic results)

```
setkey(data.df, subjid, param, agedays, id)
```

The `id` field breaks ties for same-day duplicates. Without
it in the sort key, SDE resolution order is undefined.

### Output

`cleangrowth()` returns a data.table with:

| Column | Description |
|--------|-------------|
| `id` | User-provided row identifier (preserved as-is, any type) |
| `exclude` | Exclusion code or `"Include"` |
| `param` | Growth parameter |
| `cf_rescued` | CF rescue reason (empty if not rescued) |
| `sd.orig_who` | WHO CSD z-score |
| `sd.orig_cdc` | CDC CSD z-score |
| `sd.orig` | Blended CSD z-score |
| `tbc.sd` | Recentered blended z-score |
| `ctbc.sd` | Recentered corrected z-score |
| `final_tbc` | ctbc.sd for potcorr, tbc.sd for others |

**Adult-specific output columns** (NA for child rows):

| Column | Description |
|--------|-------------|
| `mean_ht` | Subject mean included height (used in adult 2D WT step) |
| `bin_result` | Binary `"Include"`/`"Exclude"` (adult rows only) |

**Child-specific output columns** (NA for adult rows):
`cf_rescued`, `sd.orig_who`, `sd.orig_cdc`, `sd.orig`,
`tbc.sd`, `ctbc.sd`, `final_tbc`.

**Breaking change from v2.2.0:** v3.0.0 returns a data.table
with multiple columns (not a character vector of exclusion
codes like CRAN v2.2.0). Any code consuming `cleangrowth()`
output must expect a data.table, not a vector.

---

## Complete List of Exclusion Codes (child algorithm)

| Code | Step | Param | Description |
|------|------|-------|-------------|
| `Include` | — | All | Value passes all checks |
| `Missing` | Init | All | NA, NaN, agedays < 0; also HC ≥ 5y |
| `Not cleaned` | Init | HEADCM | HC with agedays > 3 × 365.25 |
| `Unit-Error-High` | 3 | HT, WT | Unit error detected (too high); corrected |
| `Unit-Error-Low` | 3 | HT, WT | Unit error detected (too low); corrected |
| `Swapped-Measurements` | 4 | HT, WT | HT/WT swapped; corrected |
| `Exclude-Carried-Forward` | 6 | All | Identical to prior-day value (not rescued) |
| `Exclude-Absolute-BIV` | 7 | All | Outside absolute biological limits |
| `Exclude-Standardized-BIV` | 7 | All | Z-score beyond age-dependent cutoffs |
| `Exclude-Evil-Twins` | 9 | All | Adjacent extreme value pair/group |
| `Exclude-EWMA1-Extreme` | 11 | All | Extreme EWMA outlier |
| `Exclude-SDE-Identical` | 13 | All | Same-day duplicate, identical value |
| `Exclude-SDE-All-Exclude` | 13 | All | All same-day values excluded |
| `Exclude-SDE-All-Extreme` | 13 | All | All same-day values extreme |
| `Exclude-SDE-EWMA` | 13 | All | SDE resolved by EWMA comparison |
| `Exclude-SDE-EWMA-All-Extreme` | 13 | All | All SDE values extreme by EWMA |
| `Exclude-SDE-One-Day` | 13 | All | All measurements on single day |
| `Exclude-Temporary-Extraneous-Same-Day` | 5 | All | Temp SDE (may persist) |
| `Exclude-EWMA2-middle` | 15 | All | Moderate EWMA: middle value |
| `Exclude-EWMA2-birth-WT` | 16 | WT | Birth weight EWMA outlier |
| `Exclude-EWMA2-birth-WT-ext` | 16 | WT | Birth weight EWMA (extended) |
| `Exclude-EWMA2-first` | 15 | All | Moderate EWMA: first value |
| `Exclude-EWMA2-first-ext` | 15 | All | Moderate EWMA: first (extended) |
| `Exclude-EWMA2-last` | 15 | All | Moderate EWMA: last value |
| `Exclude-EWMA2-last-high` | 15 | All | Moderate EWMA: last (high) |
| `Exclude-EWMA2-last-ext` | 15 | All | Moderate EWMA: last (extended) |
| `Exclude-EWMA2-last-ext-high` | 15 | All | Moderate EWMA: last (ext, high) |
| `Exclude-EWMA2-birth-HT-HC` | 16 | HT, HC | Birth HT/HC EWMA outlier |
| `Exclude-EWMA2-birth-HT-HC-ext` | 16 | HT, HC | Birth HT/HC EWMA (extended) |
| `Exclude-Min-diff` | 17 | HT, HC | Height/HC decreased > allowed |
| `Exclude-Max-diff` | 17 | HT, HC | Height/HC increased > allowed |
| `Exclude-2-meas->1-year` | 19 | All | 2 meas >1 yr apart, one excluded |
| `Exclude-2-meas-<1-year` | 19 | All | 2 meas <1 yr apart, one excluded |
| `Exclude-1-meas` | 19 | All | Single measurement excluded |
| `Exclude-Error-load` | 21 | All | Error ratio exceeds threshold |

The `exclude.levels` list includes many legacy codes from
the original pediatric algorithm for backward compatibility.
All exclusion codes will be renamed at the end of development.

### CF rescue reason codes (cf_rescued column)

| Code | Meaning |
|------|---------|
| `""` | Not CF, or CF not rescued |
| `1 CF deltaZ <0.05` | Single CF, z-score diff < 0.05 |
| `1 CF deltaZ <0.1 wholehalfimp` | Single CF, diff < 0.1, imperial |
| `Teen 2 plus CF deltaZ <0.05` | ≥2 CFs, teen, diff < 0.05 |
| `Teen 2 plus CF deltaZ <0.1 wholehalfimp` | ≥2 CFs, teen, diff < 0.1, imperial |

CF rescue thresholds may be tuned in future updates.

---

## Complete List of Steps (child algorithm)

Step numbers are not consecutive — they align with the Stata
implementation. Steps 3, 4, 8, 10, 12, 14, 18, 20 either
do not exist or are handled within other steps.

| Step | Name | Brief description |
|------|------|-------------------|
| (preprocessing) | Z-score calculation | WHO/CDC CSD z-scores, blending, rounding |
| (preprocessing) | Step 2b: GA correction | Fenton-based correction for potcorr subjects |
| (preprocessing) | Recentering | Subtract population medians → tbc.sd/ctbc.sd |
| Early 13 | SDE-Identicals | Remove same-day identical values before CF detection |
| 5 | Temporary SDE | Temporarily flag same-day duplicates |
| 6 | Carried Forwards | Identify and optionally rescue carried-forward values |
| 7 | BIV | Exclude biologically implausible values |
| 9 | Evil Twins | Exclude adjacent extreme values |
| 11 | EWMA1 (Extreme) | Exclude extreme EWMA outliers |
| 13 | Final SDE | Resolve remaining same-day duplicates |
| 15 | EWMA2 (Moderate) | Exclude moderate EWMA outliers |
| 16 | Birth HT/HC | EWMA2 variant for birth values |
| 17 | Height/HC Velocity | Exclude values exceeding velocity limits |
| 19 | Pairs and Singles | Evaluate subjects with 1–2 remaining measurements |
| 21 | Error Load | Exclude all if error ratio too high |
| 22 | Output | Assemble return columns |

---

## Adult Algorithm

### Overview

The adult algorithm (`cleanadult()`) uses a permissiveness
framework with 4 preset levels (`loosest`, `looser`, `tighter`,
`tightest`) that control all thresholds simultaneously. Default
is `looser`. Individual parameters can override presets.
Unlike the child algorithm, the adult algorithm does not compute
z-scores — it works directly with raw measurements. It uses
BMI (computed internally) for some threshold decisions.

### cleanadult() Interface

Called internally by `cleangrowth()`. Input: data.table with
`id`, `internal_id`, `subjid`, `sex`, `agedays`, `param`,
`measurement`. `internal_id` is a sequential numeric created
by `cleangrowth()` for deterministic sorting and tiebreaking.
Accepts HEIGHTCM/HEIGHTIN and WEIGHTKG/WEIGHTLBS (converts
internally). Sex is not used by the algorithm but required by
the package interface.

**Output:** Returns input df plus `result` (exclusion code),
`mean_ht`, and optionally `bin_result` (default ON),
`extraneous`, `loss_groups`, `gain_groups`. Internal columns
(`meas_m`, `ageyears`, `age_days`) are dropped. The user's
original `id` is preserved: `cleanadult()` swaps `internal_id`
into the `id` column at entry for internal processing and
restores the user's `id` at exit.

### Adult Algorithm Steps

1H/1W BIV → 2W RV markers → 3H/3W Temp SDE → 4W Weight Cap →
9Wa Evil Twins → 9H HT SDE → 9Wb Extreme EWMA →
10H HT Distinct → 10W WT SDE → 11H Mean HT →
11Wa 2D Ord WT → 11Wa2 2D Non-Ord WT → 11Wb Moderate EWMA →
13 Distinct 1D → 14 Error Load

**There is no Step 12W.**

| Step | Name | Brief description |
|------|------|-------------------|
| 1H/1W | BIV | Biologically implausible value exclusion (HT, WT, BMI) |
| 2W | RV markers | Mark repeated values for linked mode |
| 3H/3W | Temp SDE | Temporarily flag same-day duplicates |
| 4W | Weight Cap | Exclude weights at physical scale maximum |
| 9Wa | Evil Twins | Adjacent pair with implausible weight difference |
| 9H | HT SDE | Same-day height resolution (identical + extraneous) |
| 9Wb | Extreme EWMA | Extreme EWMA weight outliers |
| 10H | HT Distinct | 2D height pairs and 3+D height windows |
| 10W | WT SDE | Same-day weight resolution |
| 11H | Mean HT | Compute subject mean height (used by 2D WT) |
| 11Wa | 2D Ord WT | 2D ordered weight pairs (wtallow/perclimit) |
| 11Wa2 | 2D Non-Ord WT | 2D non-ordered weight pairs |
| 11Wb | Moderate EWMA | Moderate EWMA weight outliers + error load escalation |
| 13 | Distinct 1D | Single-measurement exclusion (HT, WT, BMI limits) |
| 14 | Error Load | Exclude all if error ratio exceeds threshold |

### Adult Exclusion Codes

#### Non-SDE Codes

| Code | Step | Description |
|------|------|-------------|
| `Include` | — | Not excluded |
| `Exclude-A-HT-BIV` | 1H | Biologically implausible height |
| `Exclude-A-WT-BIV` | 1W | Biologically implausible weight |
| `Exclude-A-WT-Scale-Max` | 4W | Weight at scale maximum |
| `Exclude-A-WT-Scale-Max-Identical` | 4W | All weights identical at scale max |
| `Exclude-A-WT-Scale-Max-RV-Propagated` | 4W | RV copy of scale-max exclusion (linked mode) |
| `Exclude-A-Evil-Twins` | 9Wa | Adjacent pair with implausible weight difference |
| `Exclude-A-WT-Traj-Ext-N` | 9Wb | Extreme EWMA outlier (independent mode) |
| `Exclude-A-WT-Traj-Extreme-firstRV-N` | 9Wb | Extreme EWMA outlier (linked firstRV pass) |
| `Exclude-A-WT-Traj-Extreme-allRV-N` | 9Wb | Extreme EWMA outlier (linked allRV pass) |
| `Exclude-A-HT-Ord-Pair` | 10Ha | 2D height pair outside band (one excluded) |
| `Exclude-A-HT-Ord-Pair-All` | 10Ha | 2D height pair outside band (all excluded) |
| `Exclude-A-HT-Window` | 10Hb | 3+D height outside window (one excluded) |
| `Exclude-A-HT-Window-All` | 10Hb | 3+D height outside window (all excluded) |
| `Exclude-A-WT-2D-Ordered` | 11Wa | 2D ordered weight pair outside wtallow/perclimit |
| `Exclude-A-WT-2D-Non-Ordered` | 11Wa2 | 2D non-ordered weight pair |
| `Exclude-A-WT-Traj-Moderate-N` | 11Wb | Moderate EWMA outlier (independent or firstRV) |
| `Exclude-A-WT-Traj-Moderate-allRV-N` | 11Wb | Moderate EWMA outlier (linked allRV pass) |
| `Exclude-A-WT-Traj-Moderate-Error-Load-N` | 11Wb | 4+ consecutive moderate EWMA candidates |
| `Exclude-A-WT-Traj-Moderate-Error-Load-RV-N` | 11Wb | Error load escalation to entire patient (linked) |
| `Exclude-A-HT-Single` | 13 | 1D height outside limits |
| `Exclude-A-WT-Single` | 13 | 1D weight outside limits |
| `Exclude-A-HT-Too-Many-Errors` | 14 | Error ratio > threshold |
| `Exclude-A-WT-Too-Many-Errors` | 14 | Error ratio > threshold |

**Note:** N in trajectory (Traj) codes = the iteration of
the EWMA exclusion loop in which the value was excluded.

#### SDE Codes

| Code | Steps | Description |
|------|-------|-------------|
| `Exclude-A-HT-Identical` | 9H | Same-day identical heights (keep one) |
| `Exclude-A-HT-Extraneous` | 9H | Same-day non-identical height (SDE loser) |
| `Exclude-A-WT-Identical` | 10W | Same-day identical weights (keep one) |
| `Exclude-A-WT-Extraneous` | 10W | Same-day non-identical weight (SDE loser) |

### Adult Permissiveness Presets

Default: `"looser"`

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| BIV HT (cm) | 50–244 | 120–230 | 142–213 | 147–208 |
| BIV WT (kg) | 20–500 | 30–270 | 36–159 | 39–136 |
| BIV BMI | 5–300 | 12–65 | 16–45 | 18–40 |
| 1D limits | split (BMI/no-BMI) | same as BIV | same as BIV | same as BIV |
| wtallow formula | PW-H (piecewise) | PW-H (piecewise) | PW-L (piecewise-lower) | allofus15 |
| UW scaling | UW-based (see wtallow-formulas.md) | UW-based | UW-based | cap limited by PW-L |
| ET caps | wtallow cap + 20 | wtallow cap + 20 | wtallow cap + 20 | allofus15-cap-12m |
| perclimit ≤45 kg | 0.5 | 0.5 | 0.7 | 0.7 |
| perclimit 45–80 kg | 0.4 | 0.4 | 0.4 | 0.4 |
| perclimit >80 kg | 0 (disabled) | 0 (disabled) | 0.4 | 0.4 |
| error_load_threshold | 0.41 | 0.41 | 0.29 | 0.29 |
| mod_ewma_f | 0.75 | 0.75 | 0.60 | 0.60 |
| ht_band | 3" | 3" | 2" | 2" |
| allow_ht_loss | TRUE | FALSE | FALSE | FALSE |
| allow_ht_gain | TRUE | TRUE | TRUE | FALSE |
| repval_handling | independent | independent | linked | linked |

### Key Differences: Adult vs. Child

| Aspect | Child | Adult |
|--------|-------|-------|
| Z-scores | CSD z-scores (WHO/CDC blend) | No z-scores — raw measurements |
| Data structure | Single `data.df`, never removes rows | Copies rows into shrinking dataframes |
| SDE tiebreaking at birth | Keep lowest `id` (pre-postnatal shift) | N/A (no births) — always keeps highest `id` |
| Permissiveness levels | Not yet implemented (planned) | 4 levels: loosest/looser/tighter/tightest |
| Parameters | param = HT, WT, HC | param = HT, WT only (no HC) |
| Rounding tolerance | None (removed) | 0.12 cm/kg on all threshold comparisons |
| Head circumference | Supported (WHO only, ≤3y cleaned) | Not applicable |
| `perclimit` scope | N/A | 11Wa: subject-level max wt; 11Wb: observation-level |
| Sort determinism | `setkey(data.df, subjid, param, agedays, internal_id)` | All sorts include `internal_id` (as character) as final tiebreaker |
| Missing-as-infinity | N/A | `ifelse(is.na(...), Inf, ...)` for edge EWMA values |

---

## Configurable Parameters (child)

| Parameter | Default | Used in | Description |
|-----------|---------|---------|-------------|
| `recover.unit.error` | FALSE | Step 3 | Attempt to identify and correct unit errors |
| `sd.extreme` | 25 | Step 7 | SD-score cutoff for extreme exclusion |
| `z.extreme` | 25 | Step 7 | Z-score cutoff for extreme exclusion |
| `height.tolerance.cm` | 2.5 | Step 17 | Maximum height decrease tolerated |
| `error.load.mincount` | 2 | Step 21 | Min exclusions before evaluating error load |
| `error.load.threshold` | 0.5 | Step 21 | Error ratio above this excludes all |
| `sd.recenter` | NA | Recentering | Child always uses fixed file |
| `include.carryforward` | FALSE | Step 6 | If TRUE, skip CF detection entirely |
| `ewma.exp` | -1.5 | EWMA steps | Exponent for EWMA weighting |
| `lt3.exclude.mode` | "default" | Step 19 | Mode for pairs/singles exclusion |
| `ewma_window` | 15 | EWMA steps | Max observations on each side for EWMA |
| `adult_cutpoint` | 20 | Preprocessing | Age (years) dividing pediatric/adult |
| `use_legacy_algorithm` | FALSE | Dispatch | If TRUE, use legacy pediatric algorithm |
| `quietly` | TRUE | All | Suppress progress messages |
| `ref_tables` | NULL | All reads | Pre-loaded closures from `gc_preload_refs()`; skips disk reads |
| `cached_results` | NULL | Partial run | data.table from prior `cleangrowth()` call; paired with `changed_subjids` |
| `changed_subjids` | NULL | Partial run | Vector of subject IDs to re-run; unchanged subjects use `cached_results` |

### gc_preload_refs()

Exported function that pre-loads all three `read_anthro()` closures once
and returns them as a named list. Avoids ~0.93 sec of disk reads per
`cleangrowth()` call (replaces 3 `read_anthro()` calls, ~0.31 sec each).
Benchmarked on 200 subjects: standard 4.46 sec → preloaded 3.40 sec
(1.06 sec saved/call, ~15 hours over 50K reps). Load time: 0.46 sec (once).

Combined with `changed_subjids` partial run (50K reps, 200 subjects):

| % subjects changed | Time/rep | vs standard | 50K saving |
|---|---|---|---|
| 100% (full run, refs only) | 3.40 sec | 1.3× | 15 hrs |
| 50% changed | 1.93 sec | 2.3× | 35 hrs |
| 25% changed | 1.16 sec | 3.8× | 46 hrs |
| 10% changed | 0.58 sec | 7.7× | 54 hrs |
| 5% changed | 0.40 sec | 11.1× | 56 hrs |

~0.2–0.3 sec fixed overhead per partial call (batching setup + merge).

```r
refs <- gc_preload_refs()
result <- cleangrowth(..., ref_tables = refs)
```

Returns: `list(mtz_cdc_prelim, mtz_who_prelim, mtz_cdc)`.

### Partial run (changed_subjids)

When measurements for only a subset of subjects have changed between
runs, pass `cached_results` (full prior results) and `changed_subjids`
(vector of subject IDs to re-run). Only those subjects are processed;
the rest are taken from `cached_results`. Subjects are independent in
all GC operations (by-group on subjid/param, fixed reference-based
recentering), so this is safe.

```r
# First run
res_baseline <- cleangrowth(..., ref_tables = refs)

# Later run — only subjects 3, 7, 12 changed
res_updated <- cleangrowth(...,
  ref_tables = refs,
  cached_results = res_baseline,
  changed_subjids = c(3, 7, 12))
```

Notes:
- Both `cached_results` and `changed_subjids` must be provided together
- If `changed_subjids` contains subject IDs not present in the input
  data, they are silently ignored
- Output row order matches input order (sorted by `id`)
- Primarily designed for Eric's secure-environment use case (Qual-AD)

### Configurable Parameters (adult, via cleangrowth)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `adult_permissiveness` | `"looser"` | Sets defaults for all adult sub-parameters |
| `adult_scale_max_lbs` | `Inf` | Physical scale upper limit in lbs (formerly `weight_cap`) |

All adult sub-parameters (BIV limits, 1D limits, wtallow
formula, etc.) can be passed individually to `cleanadult()`
to override the preset. See `adult_clean.R` roxygen for the
full list. `cleangrowth()` currently exposes only
`adult_permissiveness` and `adult_scale_max_lbs`.

---

## Key Concepts for Code Work

### The `valid()` function

**Two separate definitions exist:** `child_clean.R` (line
5044) and `pediatric_support_legacy.R` (line 10). They serve
different algorithm paths. **Do not create additional copies.**
A previous bug was caused by duplicate `valid()` definitions
where R's alphabetical file load order caused the wrong copy
to take precedence (see Fixed issues below).

Critical gatekeeper — returns a logical vector of which rows
are eligible for each step. The `include.*` flags control
whether temporarily excluded categories participate:

- `valid(df)` — only non-excluded, non-Missing, non-"Not cleaned"
- `valid(df, include.temporary.extraneous = TRUE)` — also temp SDEs
- `valid(df, include.extraneous = TRUE)` — also permanent SDEs
- `valid(df, include.carryforward = TRUE)` — also CFs

Works by checking the text of the `exclude` column, not
factor ordering.

### CSD z-scores (not LMS)

All z-scores use the Conditional Standard Deviation method:
```
If measurement < M:  sd = (measurement - M) / csd_neg
If measurement >= M: sd = (measurement - M) / csd_pos
```
This is intentionally more sensitive to extreme high values
than standard LMS z-scores. Intermediate z-score rounding
has been removed (previously used Stata-style rounding for
cross-platform validation).

### WHO/CDC age blending

| Age range | HT/WT formula | HC |
|-----------|---------------|-----|
| < 2 years | WHO only | WHO only |
| 2–5 years | `(cdc × (age-2) + who × (5-age)) / 3` | WHO only |
| > 5 years | CDC only | WHO only (through 5y) |

HC is WHO-only at all ages. HC > 3y is "Not cleaned"; HC ≥ 5y
is "Missing" (no reference data).

### Age-dependent id tiebreaking

At birth (agedays == 0): keep lowest id (earliest measurement,
before postnatal fluid shifts). All other ages: keep highest
id (later measurement, likely more careful). This differs from
the adult algorithm, which always keeps the highest id.

### Designated Other Parameter (DOP)

Weight's DOP is height; height's DOP is weight; HC's DOP is
height. Used in Steps 5, 6, 11, 13, 15, 19 for cross-parameter
plausibility checks.

### data.table reference semantics

The child algorithm keeps all rows in a single `data.df`
data.table — rows are never physically removed. Exclusions
work by setting the `exclude` column and using `valid()` to
filter. This differs from the adult algorithm, which copies
rows into separate shrinking dataframes.

`:=` modifies `data.df` in place. Be careful not to
accidentally modify copies or joined tables.

---

## Testing

### testthat suite (child)

4 test files in `tests/testthat/`:

| File | Tests | Assertions | Coverage |
|------|-------|------------|----------|
| `test-cleangrowth.R` | 4 | 54 | Legacy pediatric, adult, edge cases |
| `test-child-regression.R` | 7 | 49 | Structural invariants, frozen counts, spot checks, cross-sample stability, Missing handling, HC, prelim_infants compat |
| `test-child-parameters.R` | 10 | 24 | All configurable parameters |
| `test-child-edge-cases.R` | 12 | 23 | Single subject, sparse data, all-NA, mixed NA, SDE-Identical, negative agedays, HEADCM >3yr, extreme values, density mix, CF, deterministic |

Total child (as of 2026-03-24): 33 tests, 150 assertions.
Existing `test-cdc.R` and `test-utils.R` not modified.

Run with:
```bash
NOT_CRAN=true Rscript -e \
  'testthat::test_file("tests/testthat/test-child-regression.R")'
```

### testthat suite (adult)

| File | Tests | Assertions | Coverage |
|------|-------|------------|----------|
| `test-adult-clean.R` | 198 | 198 | All 14 steps, all 4 permissiveness levels, edge cases |

Regression test: `tests/test_harness.R` runs `cleanadult()`
against `inst/testdata/adult-gc-test-ALL-PHASES.csv` — 1508
rows at all 4 permissiveness levels.

```bash
# Unit tests
NOT_CRAN=true Rscript -e \
  'testthat::test_file("tests/testthat/test-adult-clean.R")'

# Regression test (any permissiveness level)
Rscript tests/test_harness.R [loosest|looser|tighter|tightest]
```

**Test CSV columns:** `id`, `subjid`, `param`, `agedays`,
`sex`, `measurement`, plus `expected_loosest`,
`expected_looser`, `expected_tighter`, `expected_tightest`,
`test_description`, `test_category`, `precision_sensitive`.

**Deferred test gaps:** Error load with -5 exponent, weight
scaling at permissiveness levels.

### Runtime benchmarks

| Dataset | Size | Sequential | Parallel 2 | Parallel 4 |
|---------|------|-----------|------------|------------|
| syngrowth (full) | 77,721 rows, 2,697 subj | ~114 sec | ~65 sec (1.7×) | ~30 sec (3.9×) |
| syn.csv (1,036 subj subset) | 28,434 rows | ~18 sec | — | — |
| 500-subject subsample | ~14K rows | ~11 sec | ~8 sec | ~6 sec |

Parallel produces identical results to sequential (verified 77,721 rows, 0 mismatches).
Requires installed package — see Known Issues. Run in background from Claude Code
(`run_in_background: true`).

---

## Known Issues

### Open (child)

- [x] **`stop()` → `warning()` in Step 5:** Already resolved.
  All duplicate-Include safety checks (Step 5 line ~2795,
  Step 13 lines ~3502–3504 and ~3514) use `warning()`.
  CLAUDE.md item was never ticked.
- [x] **Parallel processing fixed (2026-04-04):** Three bugs
  were causing failures:
  1. `internal_id` not passed to `temporary_extraneous_infants()`
     at 7 call sites — fixed by adding to column subset
  2. `internal_id` dropped inside the function at line 2373
     (small copy didn't include it) — fixed in keyby select
  3. `ewma_cache_init`/`ewma_cache_update` missing from child
     `var_for_par` — added. `.paropts` updated to load
     `"growthcleanr"` on workers.
  **Requirement:** Package must be installed (`devtools::install_local()`)
  before using `parallel = TRUE`. Will not work with
  `devtools::load_all()` only (workers need the installed
  package to access extdata via `system.file()`).
- [ ] **Batch size hard-coded at 2,000:** Should be a
  `batch_size` parameter on `cleangrowth()`.
- [ ] **Batch-invariant operations inside loop:**
  `exclude.levels` definition, Tanner/WHO velocity reads,
  and `read_anthro()` calls are repeated per batch. Move
  before the loop.
- [ ] **Exclusion code rename pending:** All codes will be
  renamed at end of development. Legacy codes in
  `exclude.levels` should be cleaned up at the same time.
- [ ] **CF rescue thresholds:** Current thresholds (0.05/0.1
  with wholehalfimp, teen age cutoffs) may be tuned in
  future updates.
- [ ] **LENGTHCM → HEIGHTCM:** No measurement adjustment for
  the ~0.5–0.7 cm supine/standing difference. Known
  limitation; future update planned.
- [ ] **Unused variables in Step 19:** `abs_tbd.sd`,
  `abs_ctbd.sd`, `med_dop`, `med_cdop` are computed but
  never used. Leftover from refactoring. Harmless but
  should be cleaned up.

### Open (adult)

- [ ] **Deferred test gaps:** Error load with -5 exponent,
  UW scaling edge cases (very low/high UW).
- [ ] **Performance:** `setkey(df, subjid)` optimization
  deferred.
- [ ] **Integration tests through `cleangrowth()`:** Need
  to add tests that exercise the full `cleangrowth()` →
  `cleanadult()` path, not just direct `cleanadult()` calls.

### Fixed (recent)

- [x] **`internal_id` for all internal sorting/tiebreaking
  (2026-04-12):** Both child and adult algorithms use
  `internal_id` (sequential integer `1:N`, created by
  `cleangrowth()`) for all internal sorting, tiebreaking,
  and ordering. The user's original `id` (any type) is
  preserved untouched and returned in output.
  **Child:** `internal_id` added to all 7 `.SDcols` and
  column-select locations that were missing it (lines 3344,
  3521, 4132, 4276, 4375, 4907, 5027 in `child_clean.R`).
  **Adult:** `cleanadult()` converts `internal_id` to
  character at entry (`as.character()`) so all named vector
  indexing (which uses character keys in R) works correctly.
  The user's `id` is saved, `internal_id` is swapped into
  the `id` column for internal processing, and the user's
  `id` is restored at exit. All 28 `$id` references in
  `adult_support.R` renamed to `$internal_id`.
  `cleangrowth()` no longer overwrites `data.adult$id`
  with `line`.
- [x] **Child error load bug fix (2026-04-12):** Step 21
  used hardcoded `.4` threshold instead of the configurable
  `error.load.threshold` parameter (default 0.5). Fixed to
  use the parameter. This changes behavior: 2 rows in
  syngrowth that were excluded at 0.4 are now included at
  the correct 0.5 threshold.
- [x] **Outer batching wrapper defeated:** Lines 440–441
  overwrote batch filter with all subjects. Fixed.
- [x] **Missing data bug:** NA measurements coded as
  SDE-EWMA instead of Missing. Root cause: duplicate
  `valid()` definitions in two files; R's alphabetical
  load order caused the unfixed copy to overwrite the
  fixed one. Fixed in v3.0.0.
- [x] **BIV threshold for preterm weights:** Previous
  threshold excluded legitimate preterm weights (0.7–1.0
  kg). Changed to <0.2 kg for agedays ≤ 365 and <1 kg
  for agedays > 365.
- [x] **CF rescue re-inclusion:** Rescued CFs are now set
  back to `exclude = "Include"` so they participate in
  downstream steps.
- [x] **`calc_and_recenter_z_scores()` blending bug:**
  CDC-only cutoff was `>= 4 years` instead of `> 5 years`.
  Confirmed bug — confused WHO/CDC blending window with
  corrected/uncorrected smoothing window.

---

## Next Priorities

1. **Child algorithm work** — current focus. Adult is closed
   pending clinician validation.
2. **Design extreme/clinical test patients from literature:**
   Create synthetic patients based on real clinical scenarios
   with literature-sourced growth values (e.g., hydrocephalus
   with shunt placement, craniosynostosis, failure to thrive,
   severe obesity, Turner syndrome, growth hormone deficiency).
   Goal: verify gc handles physiologically real but extreme
   trajectories correctly — not just random perturbations.

### Adult status (closed 2026-04-09)

**Do not modify adult code or tests without checking with
Carrie first.** Adult algorithm is closed pending clinician
validation. Walk-through complete (2026-04-09): 198/198 unit
tests, 1508/1508 regression tests at all 4 permissiveness
levels. Deferred items in `walkthrough-todo-2026-04-09.md`.
Completed items:
- Adult integration into package (2026-04-03)
- wtallow redesign (2026-04-09)
- Algorithm walk-through with wtallow reconciliation (2026-04-09)

---

## CRAN Preparation Checklist

Identified 2026-04-03. Items marked [x] are done.

**Timeline note (updated 2026-04-04):** CRAN is not the immediate
deadline. Priority order: (1) adult validation, (2) child
validation, (3) CRAN cleanup (can overlap with validation as long
as changes don't affect algorithm performance).

### Critical (R CMD check ERRORs/WARNINGs)

- [ ] **`parallel` not in DESCRIPTION Imports:** Imported in
  NAMESPACE but not declared in DESCRIPTION. Add to Imports.
  (~2 min)
- [ ] **Vignettes excluded via `.Rbuildignore`:** `^vignettes$`
  prevents vignette building. Either fix vignettes to build
  cleanly or remove vignette files. (~1–4 hours depending on
  vignette state)
- [ ] **Package tarball > 5 MB:** `inst/extdata/` alone is
  ~5.3 MB. `test_syngrowth_sas_output_compare.csv.gz` (1.8 MB)
  is the largest — move to a companion data package or remove
  if only needed for development. `stress_test_data.csv`
  (4.4 MB) in `tests/testthat/` also ships. (~1–2 hours to
  audit and reorganize)

### High (R CMD check NOTEs)

- [x] **Missing NULL declarations in `cleanadult()`:** Added
  `result <- mean_ht <- ... <- NULL` for data.table `:=`
  columns. (Fixed 2026-04-03)
- [x] **`cat()`/`print()` in `cleanadult()`:** Changed to
  `message()` for CRAN-preferred output handling. (Fixed
  2026-04-03)
- [ ] **`R/deprec/` and `R/modified_source_code/` not in
  `.Rbuildignore`:** R doesn't load subdirectories, but they
  bloat the tarball. Add `^R/deprec$` and
  `^R/modified_source_code$`. (~2 min)
- [ ] **`cat()`/`print()` in child algorithm:** Same issue as
  adult — `child_clean.R` and legacy files use `cat()`
  throughout. Convert to `message()`. (~30 min)

### Medium (pre-submission polish)

- [ ] **Verify `LICENSE` file exists:** DESCRIPTION says
  `MIT + file LICENSE`. Confirm plain `LICENSE` file (not
  just `LICENSE.md`) is at package root.
- [ ] **Roxygen return value outdated for `cleangrowth()`:**
  Still describes "Vector of exclusion codes" — should describe
  the data.table with all columns including adult-specific
  ones (`mean_ht`, `bin_result`). (~15 min)
- [ ] **`plyr` dependency:** Only used for `ddply` in parallel
  adult dispatch. Long-term: replace with
  `data.table`/`foreach` equivalent. (~2 hours)
- [ ] **Test runtime:** Ensure total test suite completes in
  < 10 min. Long-running tests should use `skip_on_cran()`.
  (~30 min to audit and add skips)

---

## Coding Conventions

- **Sex:** 0 = male, 1 = female
- **Param values:** `"WEIGHTKG"`, `"HEIGHTCM"`, `"HEADCM"`
  (after input conversion)
- **Use `data.table`** throughout — consistent with gc internals
- **Rounding:** Stata-style rounding (`round_stata()`) has
  been removed from the child algorithm z-score pipeline.
  Z-scores are no longer rounded to 0.001 at intermediate
  steps. The function was only needed for Stata comparison.
- **Factor levels for `exclude`:** Assigning an unlisted
  string to a factor silently produces NA — always verify
  new codes exist in `exclude.levels`
- **Sort order:** `setkey(data.df, subjid, param, agedays, id)`
  — assumed by many steps. Re-sort after any modification that
  could change order.
- **`by`-group correctness:** `by = subjid` vs.
  `by = .(subjid, param)` — wrong grouping silently produces
  wrong results. Check carefully.
- **C++/Rcpp ruled out** — CRAN target, Rtools not always
  available on enterprise systems
- **Report weight before height** in summaries and output

---

## Pitfalls (for Claude agents)

Things that have caused bugs before or fail silently:

- **Do not create a second `valid()` definition.** Two copies
  in different files caused the missing data bug (R's
  alphabetical load order made the wrong copy win). Child
  version is in `child_clean.R`; legacy version is in
  `pediatric_support_legacy.R`. Do not add a third.
- **`parallel = TRUE` requires installed package.** Will fail
  with `load_all()` only — workers need `system.file()` access
  to extdata. Install with `devtools::install_local(".")` first.
- **Do not modify sort order without re-sorting.** Many steps
  assume `setkey(data.df, subjid, param, agedays, id)`. If
  you add/modify rows, call `setkey()` again.
- **Do not assign unlisted strings to the `exclude` factor.**
  Unlisted values silently become NA. Always verify new
  exclusion codes exist in `exclude.levels` first.
- **Do not use `by = subjid` when you need
  `by = .(subjid, param)`.** Wrong grouping produces wrong
  results with no error or warning.
- **No intermediate z-score rounding.** `round_stata()` has
  been removed from the child algorithm. Do not reintroduce
  rounding at intermediate z-score steps.
- **Do not have duplicate support files in `R/`.** R loads
  all `.R` files in `R/` alphabetically. If two files define
  the same function, the later one wins silently. This caused
  breakage when both `02a_support.R` (new) and
  `adult_support.R` (old) coexisted — the old file's
  `check_between`/`round_pt` overwrote the new versions.
  **Rule:** Only one adult support file should exist in `R/`.
- **Adult `result` vs child `exclude`:** The adult algorithm
  uses `result` as its exclusion column name. `cleangrowth()`
  maps `res$result` → `exclude` in the combined output.
  Do not rename `result` to `exclude` inside `cleanadult()`.
- **Adult rounding tolerance:** 0.12 cm/kg on all threshold
  comparisons. This is intentional and should not be removed.

---

## GitHub and Distribution

- **GitHub:** `carriedaymont/growthcleanr`, branch
  `efficiency-updates`
- **CRAN:** v2.2.0 (does not include child algorithm or
  v3.0.0 changes)
- **Install from local:**
  ```r
  devtools::install_local("gc-github-latest",
    force = TRUE, upgrade = "never")
  ```
- **Install from GitHub:**
  ```r
  devtools::install_github("carriedaymont/growthcleanr",
    ref = "efficiency-updates",
    force = TRUE, upgrade = "never")
  ```
- v3.0.0 has not been released on GitHub or CRAN yet

---

## Child Narrative Document

The technical narrative for the child algorithm is in
`child-gc-narrative-2026-03-18.md` (in gc-github-latest).
It documents what the current R code does, serving as both
a debugging aid and long-term reference.

### Status

**Front matter: COMPLETE.** Covers Key Concepts, Architecture,
Data Requirements, Z-Score Infrastructure, Working/Output
Dataframe, Variable Glossary, Output Format (with full
exclusion code tables), Configurable Parameters, and Complete
List of Steps.

**Step-by-step documentation: COMPLETE for all steps.**
Steps 5, 6, 7, 9, 11, 13, 15, 16, 17, 19, 21, 22 all
documented with summary tables, overview, logic details,
rationale, and code review checklist findings.

**Open questions and issues:** Tracked in the narrative's
"Open Questions and Issues to Investigate" section. Items
marked `[x]` are resolved; `[ ]` items are still open.

### Structure

Each step section includes:
- Summary table (scope, prior/next step, exclusion codes)
- Overview
- Key terms and variable names
- Logic and implementation
- Rationale for selected decisions
- Code review checklist findings (12-point checklist)

---

## Terms Reference

This is the **canonical** terms reference (primary source).
`__Pipeline/CLAUDE.md` has a convenience copy; update both
if adding terms.

| Term | Full name |
|------|-----------|
| ewma | Exponentially weighted moving average |
| dewma | Deviation from EWMA |
| absdewma | Absolute value of dewma |
| spa | Subject-parameter-ageday combination |
| sde | Same-day extraneous value |
| dop | Designated other parameter |
| tbc | To-be-cleaned (recentered z-score) |
| ctbc | Corrected to-be-cleaned (recentered corrected z-score) |
| potcorr | Potentially correctable (prematurity) |
| biv | Biologically implausible value |
| csd | Conditional standard deviation |
| nnte | No need to evaluate (legacy; always FALSE in current code) |
| otl | Out of line (used in Evil Twins step; formerly "oob"/"out of bounds") |

---

## Running growthcleanr

### Standalone (from gc-github-latest directory)

```r
library(growthcleanr)
library(data.table)

dt <- fread("path/to/data.csv")
result <- cleangrowth(
  subjid = dt$subjid,
  param = dt$param,
  agedays = dt$agedays,
  sex = dt$sex,
  measurement = dt$measurement,
  parallel = FALSE
)
```

### From Pipeline (error-impact)

See `__Pipeline/CLAUDE.md` for pipeline-specific wrapper
instructions. Working directory must be `error-impact/` so
`here::here()` resolves correctly.
