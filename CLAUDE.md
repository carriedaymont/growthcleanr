# CLAUDE.md — gc-github-latest (growthcleanr)

**Last updated:** 2026-04-16 (exclusion code walkthrough: Traj-Ext→Traj-Extreme, missing factor levels, narrative sync)

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
| `R/child_clean.R` | `cleangrowth()` entry point + exports (top of file); `gc_preload_refs()` (pre-loads reference closures for repeated calls); `cleanchild()` main algorithm + support functions (`.child_valid()`, `.child_exc()`, `temporary_extraneous_infants()`, `calc_otl_evil_twins()`, `calc_and_recenter_z_scores()`, `ewma()`, `ewma_cache_init()`/`ewma_cache_update()`, `get_dop()`, `read_anthro()`) |
| `R/pediatric_clean_legacy.R` | Legacy pediatric algorithm (`use_legacy_algorithm = TRUE`); deprecated, retained for backward compatibility |
| ~~`R/growth.R`~~ | Removed in v3.0.0 — `cleangrowth()` moved into `child_clean.R` |
| `R/utils.R` | Shared utilities (does NOT contain `.child_valid()` or `valid()` — see note below) |
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
- Batching and dispatch to `cleanchild()`

`cleanchild()` processes one batch through all
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

Child exclusion codes follow the `Exclude-C-{Reason}` format.
Codes are NOT param-specific — the param for each row is in
the data (the `param` column). The `.child_exc()` helper
function (line ~142 of `child_clean.R`) generates these codes:
`.child_exc(reason)` → `Exclude-C-{reason}`.

| Code | Step | Param | Description |
|------|------|-------|-------------|
| `Include` | — | All | Value passes all checks |
| `Exclude-Missing` | Init | All | NA, NaN, agedays < 0; also HC ≥ 5y |
| `Exclude-Not-Cleaned` | Init | HEADCM | HC with agedays > 3 × 365.25 |
| `Exclude-C-CF` | 6 | All | Carried forward: identical to prior-day value (not rescued) |
| `Exclude-C-BIV` | 7 | All | Biologically implausible value (absolute or standardized) |
| `Exclude-C-Evil-Twins` | 9 | All | Adjacent extreme value pair/group |
| `Exclude-C-Traj-Extreme` | 11 | All | Extreme EWMA outlier (EWMA1) |
| `Exclude-C-Identical` | Early 13, 13 | All | Same-day duplicate, identical value |
| `Exclude-C-Extraneous` | 13 | All | Same-day extraneous value (SDE loser by EWMA or other criteria) |
| `Exclude-C-Traj` | 15, 16 | All | Moderate EWMA outlier (EWMA2): covers middle, first, last, birth WT, birth HT/HC, and all sub-variants |
| `Exclude-C-Abs-Diff` | 17 | HT, HC | Height/HC velocity violation (increase > max or decrease > min allowed) |
| `Exclude-C-Pair` | 19 | All | Subject with 2 measurements, one excluded |
| `Exclude-C-Single` | 19 | All | Subject with 1 measurement, excluded |
| `Exclude-C-Too-Many-Errors` | 21 | All | Error ratio exceeds threshold |
| `Exclude-C-Temp-Same-Day` | 5 | All | **Internal only** — temporary SDE flag, must be resolved before output |

Codes are shared across all params — use the `param` column
in the data to determine which parameter a code applies to.

**`Exclude-Missing` and `Exclude-Not-Cleaned`** are shared codes
assigned in `cleangrowth()` before dispatch.

**`exclude.levels`** has been updated to include all child and
adult codes (was previously missing adult codes, causing silent
NAs in factor output).

**Codes NOT implemented (Steps 3/4):** Unit error recovery
(`recover.unit.error`) and measurement swapping are legacy
features scheduled for redesign. No `Unit-Error-*` or
`Swapped-Measurements` codes are assigned by the child
algorithm.

**Simplified vs. detailed codes:** Unlike the adult algorithm
which has distinct codes for each sub-step (e.g., separate
moderate EWMA codes for firstRV vs allRV), the child algorithm
uses simplified codes. For example, all EWMA2 exclusions use
`Traj` regardless of whether the value was a middle, first,
last, or birth measurement. All SDE resolutions use either
`Identical` or `Extraneous`. The specific sub-step logic is
in the code comments but not reflected in the exclusion code.

**Legacy codes:** Exclusion codes in `pediatric_clean_legacy.R`
were NOT renamed — the legacy algorithm retains the old code
names for backward compatibility.

### CF rescue

**Parameter:** `cf_rescue = c("standard", "none", "all")`
- `"standard"` (default): age/interval/param-specific lookup thresholds
- `"none"`: no rescue (all CFs excluded)
- `"all"`: all CFs rescued (no CF exclusions)

**cf_rescued column values:**

| Code | Meaning |
|------|---------|
| `""` | Not CF, or CF not rescued |
| `"Rescued"` | CF rescued by lookup threshold (standard mode) |
| `"Rescued-All"` | CF rescued because cf_rescue="all" |

**Optional cf_detail columns** (enabled by `cf_detail = TRUE`):
- `cf_status`: NA (not CF candidate), "CF-NR" (excluded), or "CF-Resc" (rescued)
- `cf_deltaZ`: absolute z-score difference from originator

CF rescue thresholds are defined in `cf-rescue-thresholds.md` —
lookup tables by age bin × interval bin × param × rounding type.
Three levels (0.05 / 0.20 / 0.40) plus NR (no rescue).
Derived from 100K synthetic tracker population; full methodology
in `__Pipeline/CF-exploration/cf-threshold-schemes.md`.

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
| `Exclude-A-BIV` | 1H/1W | Biologically implausible value (height or weight) |
| `Exclude-A-Scale-Max` | 4W | Weight at scale maximum |
| `Exclude-A-Scale-Max-Identical` | 4W | All weights identical at scale max |
| `Exclude-A-Scale-Max-RV-Propagated` | 4W | RV copy of scale-max exclusion (linked mode) |
| `Exclude-A-Evil-Twins` | 9Wa | Adjacent pair with implausible weight difference |
| `Exclude-A-Traj-Extreme` | 9Wb | Extreme EWMA outlier (independent mode) |
| `Exclude-A-Traj-Extreme-firstRV` | 9Wb | Extreme EWMA outlier (linked firstRV pass) |
| `Exclude-A-Traj-Extreme-firstRV-RV-Propagated` | 9Wb | RV copy of firstRV extreme exclusion (linked mode) |
| `Exclude-A-Traj-Extreme-allRV` | 9Wb | Extreme EWMA outlier (linked allRV pass) |
| `Exclude-A-Ord-Pair` | 10Ha | 2D height pair outside band (one excluded) |
| `Exclude-A-Ord-Pair-All` | 10Ha | 2D height pair outside band (all excluded) |
| `Exclude-A-Window` | 10Hb | 3+D height outside window (one excluded) |
| `Exclude-A-Window-All` | 10Hb | 3+D height outside window (all excluded) |
| `Exclude-A-2D-Ordered` | 11Wa | 2D ordered weight pair outside wtallow/perclimit |
| `Exclude-A-2D-Non-Ordered` | 11Wa2 | 2D non-ordered weight pair |
| `Exclude-A-Traj-Moderate` | 11Wb | Moderate EWMA outlier (independent or firstRV) |
| `Exclude-A-Traj-Moderate-RV-Propagated` | 11Wb | RV copy of firstRV moderate exclusion (linked mode) |
| `Exclude-A-Traj-Moderate-allRV` | 11Wb | Moderate EWMA outlier (linked allRV pass) |
| `Exclude-A-Traj-Moderate-Error-Load` | 11Wb | 4+ consecutive moderate EWMA candidates |
| `Exclude-A-Traj-Moderate-Error-Load-RV` | 11Wb | Error load escalation to entire patient (linked) |
| `Exclude-A-Single` | 13 | 1D measurement outside limits (height or weight) |
| `Exclude-A-Too-Many-Errors` | 14 | Error ratio > threshold |

Adult codes are NOT param-specific — the param for each row
is in the data (the `param` column).

**Note:** The `-N` round suffix (iteration number) was removed
from adult trajectory codes in the 2026-04-14 exclusion code
rename. Trajectory codes no longer include the loop iteration.

#### SDE Codes

| Code | Steps | Description |
|------|-------|-------------|
| `Exclude-A-Identical` | 9H/10W | Same-day identical values (keep one) |
| `Exclude-A-Extraneous` | 9H/10W | Same-day non-identical value (SDE loser) |

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
| `error.load.mincount` | 2 | Step 21 | Min exclusions before evaluating error load |
| `error.load.threshold` | 0.5 | Step 21 | Error ratio above this excludes all |
| `sd.recenter` | NA | Recentering | Child always uses fixed file |
| `include.carryforward` | FALSE | Step 6 | **Deprecated** — use `cf_rescue` instead. If TRUE, skip CF detection entirely |
| `cf_rescue` | `"standard"` | Step 6 | CF rescue mode: `"standard"` (age/interval/param-specific lookup), `"none"` (all CFs excluded), `"all"` (all CFs rescued) |
| `cf_detail` | FALSE | Step 6 | If TRUE, add `cf_status` and `cf_deltaZ` columns to output |
| `ewma.exp` | -1.5 | EWMA steps | **Legacy only.** Exponent for EWMA weighting. Child algorithm uses age-dependent exponents internally |
| `ewma_window` | 15 | EWMA steps | Max observations on each side for EWMA |
| `adult_cutpoint` | 20 | Preprocessing | Age (years) dividing pediatric/adult |
| `lt3.exclude.mode` | `"default"` | Legacy | Passed to `cleanlegacy()` — re-added after accidental removal |
| `use_legacy_algorithm` | FALSE | Dispatch | If TRUE, use legacy pediatric algorithm |
| `quietly` | TRUE | All | Suppress progress messages |
| `ref_tables` | NULL | All reads | Pre-loaded closures from `gc_preload_refs()`; skips disk reads |
| `cached_results` | NULL | Partial run | data.table from prior `cleangrowth()` call. Auto-detects changed subjects if `changed_subjids` is NULL; uses explicit list if both provided |
| `changed_subjids` | NULL | Partial run | Optional vector of subject IDs to re-run. If NULL and `cached_results` provided, changed subjects are auto-detected |
| `batch_size` | 2000 | Batching | Number of subjects per processing batch (was hardcoded, now configurable) |

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

### Partial run (cached_results)

When running GC repeatedly on data that differs for only a subset of
subjects (e.g., baseline vs. error-injected data, or updated data
extracts), pass `cached_results` from a prior full run. Only subjects
with changed data are re-processed; unchanged subjects receive their
cached results. Subjects are independent in all GC operations
(by-group on subjid/param, fixed reference-based recentering), so
this is safe.

**Two modes:**

1. **Auto-detect** (recommended): Pass `cached_results` only.
   `cleangrowth()` compares the incoming input data against the
   cached results to automatically identify subjects with any
   row-level differences (added/removed/modified measurements).
   Only those subjects are re-processed.

2. **Explicit**: Pass both `cached_results` and `changed_subjids`.
   Only subjects in `changed_subjids` are re-processed. No
   auto-detection is performed.

```r
# First run (full baseline)
refs <- gc_preload_refs()
res_baseline <- cleangrowth(..., ref_tables = refs)

# Auto-detect mode: pass modified data + cached baseline
# GC figures out which subjects differ
res_error <- cleangrowth(
  subjid = d_err$subjid, param = d_err$param,
  agedays = d_err$agedays, sex = d_err$sex,
  measurement = d_err$measurement, id = d_err$id,
  ref_tables = refs,
  cached_results = res_baseline)

# Explicit mode: caller provides the list of changed subjects
res_error <- cleangrowth(...,
  ref_tables = refs,
  cached_results = res_baseline,
  changed_subjids = c(3, 7, 12))
```

**Auto-detection details:**
- Compares input subjid, param, agedays, sex, and measurement
  (with the same 0→NaN transformation GC applies internally)
- Detects added subjects (in input but not cache), removed
  subjects (in cache but not input), and modified subjects
  (any row-level difference in the comparison columns)
- Removed subjects simply won't appear in output (they're not
  in the input data)
- With `quietly = FALSE`, prints a summary:
  `Auto-detected N changed subjects (X added, Y modified, Z removed, W unchanged)`

**Performance (50 subjects, 3 modified):**
- Full run: ~1.05 sec
- Auto-detect partial: ~1.0 sec (comparison overhead dominates
  at small scale; savings increase with larger datasets and
  fewer changed subjects)
- Explicit partial: ~0.29 sec (no comparison overhead)
- No changes detected: ~0.11 sec (returns cache immediately)

Notes:
- If `changed_subjids` contains subject IDs not present in the
  input data, they are silently ignored
- Output row order matches input order (sorted by `internal_id`)
- Designed for two use cases: (1) error-injection pipeline
  workflows where most subjects are unmodified, and (2) Eric's
  secure-environment workflow with updated data extracts

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

### The `.child_valid()` function (formerly `valid()`)

**Renamed from `valid()` to `.child_valid()` (2026-04-14)** to
avoid collision with the legacy `valid()` in
`pediatric_support_legacy.R`. R's alphabetical file loading
caused the legacy version to overwrite the child version when
both were named `valid()`. The legacy `valid()` remains
unchanged. **Do not create additional copies of either function.**

Critical gatekeeper — returns a logical vector of which rows
are eligible for each step. The `include.*` flags control
whether temporarily excluded categories participate:

- `.child_valid(df)` — only non-excluded, non-Exclude-Missing, non-Exclude-Not-Cleaned
- `.child_valid(df, include.temporary.extraneous = TRUE)` — also temp SDEs
- `.child_valid(df, include.extraneous = TRUE)` — also permanent SDEs
- `.child_valid(df, include.carryforward = TRUE)` — also CFs

Works by checking the text of the `exclude` column, not
factor ordering. Uses the new `Exclude-Missing` and
`Exclude-Not-Cleaned` code names.

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

### testthat suite

Full details in `testing-reference.md` (test-by-test listing,
coverage gaps, run instructions).

6 actively maintained test files in `tests/testthat/`:

| File | Tests | Assertions | Scope |
|------|-------|------------|-------|
| `test-cleangrowth.R` | 14 | 92 | `cleangrowth()` API: legacy pediatric, adult integration (7 tests), child-adult spanning (3 tests), deprecation |
| `test-child-regression.R` | 9 | 54 | Frozen counts, spot checks, cross-sample stability, Missing, HC, preload_refs, changed_subjids |
| `test-child-algorithms.R` | 24 | ~40 | CF rescue (modes + threshold cells), Evil Twins/OTL (1/2/3 unit errors + collateral), Error Load (threshold/mincount/constructed), HC boundary, cf_detail, parallel, GA correction (potcorr + near-potcorr), Birth EWMA2 (extreme/normal WT, extreme HT), LENGTHCM identity |
| `test-child-parameters.R` | 9 | 21 | use_legacy_algorithm, sd.recenter, include.carryforward, sd.extreme, ewma_window, error.load params, recover.unit.error, imperial units, LENGTHCM |
| `test-child-edge-cases.R` | 12 | 24 | Single subject, sparse data, all-NA, mixed NA, SDE-Identical, negative agedays, HEADCM >3yr, extreme values, density mix, CF, deterministic |
| `test-adult-clean.R` | 198 | 198 | All 14 adult steps, 4 permissiveness levels, edge cases (via `cleanadult()` directly) |
| **Total** | **266** | **~429** | |

Additional test files not actively maintained:
`test-cdc.R` (not modified in v3.0.0); `test-utils.R` (10
failures — old API without `id` parameter, tests utility
functions `splitinput`, `recode_sex`, `longwide`, `simple_bmi`).

**Not tested (by design):** Steps 3 (unit error recovery) and
4 (swapped measurements) are legacy features scheduled for
redesign. Do not write tests for the current implementation —
they would need to be rewritten anyway and would slow down
the redesign.

Run with:
```bash
NOT_CRAN=true Rscript -e \
  'testthat::test_file("tests/testthat/test-child-regression.R")'
```

### Adult regression harness (separate from testthat)

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

**IMPORTANT: Reinstall before testing.** Tests load
`library(growthcleanr)`, which uses the *installed* package,
not source files. After any code changes, always reinstall
before running tests:
```r
devtools::install_local("gc-github-latest",
  force = TRUE, upgrade = "never")
```
Without this, tests run against stale installed code and
results are misleading — tests may pass or fail for the
wrong reasons. (`devtools::load_all()` loads source code
in the current R session but does not affect fresh `Rscript`
processes used by `testthat::test_file()`.)

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
- [x] **Batch size parameterized (2026-04-15):** `batch_size`
  parameter added to `cleangrowth()` (default 2000).
- [x] **Batch-invariant operations moved before loop
  (2026-04-13):** `exclude.levels` definition and
  Tanner/WHO velocity reads moved before the outer loop.
  `read_anthro()` calls already optimized via `ref_tables`.
- [x] **Exclusion code rename (2026-04-14):** Child codes
  renamed to `Exclude-C-{Reason}` format; adult `-N` round
  suffix removed from trajectory codes; `Missing` →
  `Exclude-Missing`; `Not cleaned` → `Exclude-Not-Cleaned`.
  `exclude.levels` updated to include all adult codes. Legacy
  codes in `pediatric_clean_legacy.R` were NOT renamed.
  Param indicators (-WT-, -HT-, -HC-) removed from all codes
  (2026-04-16) — codes are not param-specific; param is in
  the data.
- [x] **CF rescue thresholds (2026-04-14):** Replaced fixed
  thresholds (0.05/0.10 with wholehalfimp, teen age cutoffs)
  with age/interval/param/rounding-specific lookup via
  `.cf_rescue_lookup()`. Three levels: 0.05 (slow growth),
  0.20 (moderate), 0.40 (fast growth), plus NR (no rescue).
  No string-length or teen restriction. Three modes via
  `cf_rescue` parameter: "standard" (lookup), "none", "all".
  `include.carryforward` deprecated. Thresholds derived from
  100K synthetic tracker population; methodology in
  `__Pipeline/CF-exploration/cf-threshold-schemes.md`.
- [ ] **LENGTHCM → HEIGHTCM:** No measurement adjustment for
  the ~0.5–0.7 cm supine/standing difference. Known
  limitation; future update planned.
- [x] **Unused variables in Step 19:** `abs_tbd.sd`,
  `abs_ctbd.sd`, `med_dop`, `med_cdop` — confirmed removed
  in current code (verified 2026-04-15).

### Robustness audit (2026-04-14)

Comprehensive audit of code patterns that could warn or
degrade on large/unusual datasets. No blocking issues found.

**Investigated and confirmed safe:**
- **Division by zero in adult slope calculations**
  (`adult_support.R` ~lines 823, 836): `(ap1 - ap2)` and
  `(an1 - an2)` denominators cannot be zero because all
  same-day duplicates are resolved before these steps run
  (only one Include per subject-param-ageday). If zero
  somehow occurs, it indicates an upstream bug and should
  break loudly.
- **Division by zero in adult BMI** (`adult_support.R`
  ~line 1476): `ht_val` comes from included heights that
  already passed BIV (min 50cm at loosest). Cannot be zero.
- **HT+WT same-day pairing assumption:** Both algorithms
  handle unpaired measurements correctly. Adult uses
  `intersect()` for BMI days — skips if no shared days.
  Child DOP is a secondary tiebreaker with explicit NA→Inf
  fallback when the other param doesn't exist for a subject.

**Low priority cleanup (functional but could be cleaner):**
- [ ] **`suppressWarnings(min/max)` + `is.infinite()` pattern**
  (child_clean.R ~lines 4006, 4109, 4117, 4126): Suppresses
  warnings from `min()`/`max()` on empty subsets, then
  converts `-Inf`/`Inf` to `NA`. Works correctly but masks
  any unexpected warnings. Could be rewritten with explicit
  `if (length(...) == 0)` guards to match the style of the
  5 locations fixed in commit `64c186b`.
- [ ] **`which.min`/`which.max` on empty vectors** (child,
  ~lines 4442, 5201): Returns `integer(0)` on empty input.
  Surrounding logic has row-count checks, so safe in
  practice.
- [ ] **Growing vectors with `c()` in loops** (adult,
  ~lines 291-292, 559): O(n²) pattern — `exc_ids <- c(exc_ids, new_id)`.
  Negligible for typical subject sizes but could slow down
  on subjects with thousands of measurements.

### Open (adult)

- [ ] **Deferred test gaps:** Error load with -5 exponent,
  UW scaling edge cases (very low/high UW).
- [ ] **Performance:** `setkey(df, subjid)` optimization
  deferred.
- [x] **Integration tests through `cleangrowth()`:** 7 tests
  added to `test-cleangrowth.R` exercising the full
  `cleangrowth()` → `cleanadult()` path: output columns,
  imperial units, permissiveness passthrough, scale_max_lbs,
  weight_cap deprecation, id preservation, spanning subject
  output structure.

### Fixed (recent)

- [x] **Exclusion code walkthrough (2026-04-16):**
  Full walkthrough focused on exclusion codes. Found and fixed:
  missing `Exclude-A-Traj-Extreme-firstRV-RV-Propagated` factor
  level (latent bug — would silently produce NA if triggered);
  `Exclude-A-Traj-Ext` → `Exclude-A-Traj-Extreme` rename for
  consistency with linked mode codes; silently broken error-load
  test (`test-child-parameters.R` used old code name, always
  matched 0 rows); stale `-N` suffix in adult narrative (4
  locations); child narrative updated to remove param indicators
  (~37 occurrences); 3 stale code comments fixed; added
  `Exclude-A-Traj-Moderate-RV-Propagated` and
  `Exclude-A-Traj-Extreme-firstRV-RV-Propagated` to CLAUDE.md
  canonical table. All tests pass. Details in
  `walkthrough-todo-2026-04-16.md`.
- [x] **Remove param from exclusion codes (2026-04-16):**
  Removed param indicators (-WT-, -HT-, -HC-) from all child
  and adult exclusion codes. Child: `Exclude-C-{WT|HT|HC}-{Reason}`
  → `Exclude-C-{Reason}` (e.g., `Exclude-C-WT-BIV` →
  `Exclude-C-BIV`). Adult: `Exclude-A-{WT|HT}-{Reason}` →
  `Exclude-A-{Reason}` (e.g., `Exclude-A-WT-BIV` →
  `Exclude-A-BIV`). Codes are not param-specific — the param
  for each row is in the data. `.child_exc()` updated to take
  only reason (not param). `exclude.levels` updated.
- [x] **Full walkthrough cleanup (2026-04-15):** Resolved 30+
  deferred items from the full code walkthrough. Key changes:
  all `cat()`/`print()` → `message()` for CRAN compliance
  (50+ calls); ~70 lines of dead/commented-out code removed;
  `batch_size` parameter added (was hardcoded 2000);
  unnecessary `data.adult[, id := line]` removed; batching
  dplyr → data.table; `|` → `||` short-circuit fix; `names(
  table(subjid) > 1)` correctness fix; stale `"valid"` →
  `".child_valid"` in `var_for_par`; NA fallback logic added
  to `calc_and_recenter_z_scores`; `as_matrix_delta()` copied
  to child_clean.R; `ewma.exp` documented as legacy-only;
  conflict warning for `include.carryforward` + `cf_rescue`;
  Stata comments condensed; clarifying comments throughout.
  Full findings in `walkthrough-todo-2026-04-15.md`. Commit
  `0bea2ba`.
- [x] **min/max warnings on empty subsets (2026-04-14):**
  Fixed 5 locations in `child_clean.R` where `min()`/`max()`
  on empty Include subsets produced user-visible `-Inf`
  warnings. All replaced with explicit `if (length(...) == 0)`
  guards. Commit `64c186b`.
- [x] **Child BIV tests (2026-04-14):** Added 3 separate
  focused BIV tests (WT, HT, HC) to `test-child-edge-cases.R`.
  Each uses a minimal dataset with one extreme value and
  verifies the correct BIV exclusion code. Commit `8201e9c`.
- [x] **Test suite overhaul (2026-04-14):** Adult tests
  fixed for `library()` loading (was missing `library(growthcleanr)`
  and calling unexported `cleanadult()` directly). Added 7
  adult integration tests and 3 spanning subject tests to
  `test-cleangrowth.R`. `permissiveness_presets()` exported.
  `testing-reference.md` created. Commit `1797135`.
- [x] **Exclusion code rename (2026-04-14):** All child
  exclusion codes renamed from `Exclude-{Reason}` to
  `Exclude-C-{Reason}`, using `.child_exc(reason)` helper.
  `Missing` → `Exclude-Missing`, `Not cleaned` →
  `Exclude-Not-Cleaned` (shared codes assigned in
  `cleangrowth()`). Adult: `-N` round suffix removed from
  trajectory codes. `exclude.levels` updated to include all
  adult codes (was missing them, causing silent NAs in factor
  output). `valid()` renamed to `.child_valid()` in
  `child_clean.R` to avoid collision with legacy `valid()` in
  `pediatric_support_legacy.R`. `lt3.exclude.mode` parameter
  re-added to `cleangrowth()` with default `"default"` (was
  accidentally removed but `cleanlegacy()` still requires it).
  Legacy codes in `pediatric_clean_legacy.R` NOT renamed.
  Tests: child 200/207 (7 pre-existing CF rescue + parallel
  failures), adult 198/198 unit + 1508/1508 regression at
  all 4 levels, test-utils.R 6 pre-existing failures.
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
  fixed one. Fixed in v3.0.0. Further prevented by
  renaming child version to `.child_valid()` (2026-04-14).
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
- [x] **`cat()`/`print()` in child algorithm (2026-04-15):**
  All 50+ `cat()`/`print()` calls in `child_clean.R`
  converted to `message()`. Legacy files not changed.

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

- **Do not create additional `valid()`/`.child_valid()`
  definitions.** Two copies in different files caused the
  missing data bug (R's alphabetical load order made the wrong
  copy win). Child version is `.child_valid()` in
  `child_clean.R`; legacy version is `valid()` in
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
