# CLAUDE.md — gc-github-latest (growthcleanr)

**Last updated:** 2026-04-17 (Session 4a CF walkthrough: removed run_cf_detection optimization; cf_rescue="all" now rescues every detected CF including shared-SPA ones, with Step 13 resolving multi-Include SPAs; HEADCM threshold placeholder moved to Known Issues; pre-session cleanup migrated D5/D6/D9+D10 to Known Issues → Open (wrapper) and fixed D4/D16 in child_clean.R)

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

**Legacy removed (2026-04-16):** The legacy pediatric algorithm
(`pediatric_clean_legacy.R`, `pediatric_support_legacy.R`),
`R/deprec/`, `R/modified_source_code/`, and associated
parameters (`use_legacy_algorithm`, `prelim_infants`,
`lt3.exclude.mode`, `ewma.exp`, `recover.unit.error`) have all
been removed. Prior CRAN versions (v2.2.0) retain the legacy
code for backward compatibility.

This CLAUDE.md covers both the child and adult algorithms.

---

## Code Structure

### Key files

| File | What it contains |
|------|------------------|
| `R/child_clean.R` | `cleangrowth()` entry point + exports (top of file); `gc_preload_refs()` (pre-loads reference closures for repeated calls); `cleanchild()` main algorithm + support functions (`.child_valid()`, `.child_exc()`, `identify_temp_sde()`, `calc_otl_evil_twins()`, `calc_and_recenter_z_scores()`, `ewma()`, `as_matrix_delta()` (internal), `ewma_cache_init()`/`ewma_cache_update()`, `get_dop()`, `read_anthro()`) |
| `R/utils.R` | Shared utilities (does NOT contain `.child_valid()`) |
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
   into batches (default 2,000, configurable via `batch_size`).

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

**Simplified vs. detailed codes:** Unlike the adult algorithm
which has distinct codes for each sub-step (e.g., separate
moderate EWMA codes for firstRV vs allRV), the child algorithm
uses simplified codes. For example, all EWMA2 exclusions use
`Traj` regardless of whether the value was a middle, first,
last, or birth measurement. All SDE resolutions use either
`Identical` or `Extraneous`. The specific sub-step logic is
in the code comments but not reflected in the exclusion code.

### CF rescue

**Parameter:** `cf_rescue = c("standard", "none", "all")`
- `"standard"` (default): age/interval/param-specific lookup thresholds
- `"none"`: no rescue (all CFs excluded)
- `"all"`: every detected CF rescued, including CFs on an SPA that
  also has another Include. Step 13 final-SDE resolution handles any
  resulting multi-Include SPAs. Use when the caller wants to treat
  CFs as plausible and let downstream SDE logic pick among candidates.

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
`measurement`. `internal_id` is a sequential integer created
by `cleangrowth()` for deterministic sorting and tiebreaking.
Accepts HEIGHTCM/HEIGHTIN and WEIGHTKG/WEIGHTLBS (converts
internally). Sex is not used by the algorithm but required by
the package interface.

**Output:** Returns input df plus `result` (exclusion code),
`mean_ht`, and optionally `bin_result` (default ON),
`extraneous`, `loss_groups`, `gain_groups`. Internal columns
(`meas_m`, `ageyears`, `age_days`) are dropped. The user's
original `id` is preserved untouched. `internal_id` (integer)
is used for all internal sorting and named vector indexing
(via `as.character()` for names, `as.integer()` for join keys).

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
| Sort determinism | `setkey(data.df, subjid, param, agedays, internal_id)` | All sorts include `internal_id` (integer) as final tiebreaker; `as.character()` used for named vector keys |
| Missing-as-infinity | N/A | `ifelse(is.na(...), Inf, ...)` for edge EWMA values |

---

## Configurable Parameters (child)

| Parameter | Default | Used in | Description |
|-----------|---------|---------|-------------|
| `sd.extreme` | 25 | Step 7 | SD-score cutoff for extreme exclusion |
| `z.extreme` | 25 | Step 7 | Z-score cutoff for extreme exclusion |
| `error.load.mincount` | 2 | Step 21 | Min exclusions before evaluating error load |
| `error.load.threshold` | 0.5 | Step 21 | Error ratio above this excludes all |
| `sd.recenter` | NA | Recentering | Default uses built-in rcfile; pass a data.table for custom recentering |
| `include.carryforward` | FALSE | Step 6 | **Deprecated** — use `cf_rescue` instead. If TRUE, skip CF detection entirely |
| `cf_rescue` | `"standard"` | Step 6 | CF rescue mode: `"standard"` (age/interval/param-specific lookup), `"none"` (all CFs excluded), `"all"` (every detected CF rescued; Step 13 resolves any multi-Include SPAs) |
| `cf_detail` | FALSE | Step 6 | If TRUE, add `cf_status` and `cf_deltaZ` columns to output |
| `debug` | FALSE | Step 11, output | If TRUE, emit 6 `ewma1_it1.*` columns capturing EWMA values from the first iteration of Step 11 (EWMA1 Extreme), useful for diagnosing why a row was flagged. Off by default to keep standard output lean. |
| `ewma_window` | 15 | EWMA steps | Max observations on each side for EWMA |
| `adult_cutpoint` | 20 | Preprocessing | Age (years) dividing pediatric/adult |
| `quietly` | TRUE | All | Suppress progress messages |
| `ref_tables` | NULL | All reads | Pre-loaded closures from `gc_preload_refs()`; skips disk reads |
| `cached_results` | NULL | Partial run | data.table from prior `cleangrowth()` call. Auto-detects changed subjects if `changed_subjids` is NULL; uses explicit list if both provided |
| `changed_subjids` | NULL | Partial run | Optional vector of subject IDs to re-run. If NULL and `cached_results` provided, changed subjects are auto-detected |
| `batch_size` | 2000 | Batching | Number of subjects per processing batch |

### gc_preload_refs() and partial runs (cached_results)

Canonical documentation for both features lives in
`child-gc-narrative-2026-04-13.md` under "Partial runs and
preloaded references" (in the Batching and Dispatch section).
Quick reference:

- `gc_preload_refs()` returns a list of CDC/WHO reference
  closures (`list(mtz_cdc_prelim, mtz_who_prelim)`) that can be
  passed to `cleangrowth()` via `ref_tables = refs` to skip
  per-call disk reads (~0.9 sec/call saved).
- `cached_results` + optional `changed_subjids` let you
  re-process only subjects whose input rows changed vs. a prior
  run. Auto-detect mode (omit `changed_subjids`) compares
  subjid/param/agedays/sex/measurement against the cache.

Observed speedups (200 subjects, 50 000 calls):

| % subjects changed | Time/rep | vs standard | 50K saving |
|---|---|---|---|
| 100% (full run, refs only) | 3.40 s | 1.3× | 15 hrs |
| 50% changed | 1.93 s | 2.3× | 35 hrs |
| 25% changed | 1.16 s | 3.8× | 46 hrs |
| 10% changed | 0.58 s | 7.7× | 54 hrs |
| 5% changed | 0.40 s | 11.1× | 56 hrs |

Designed for two workflows: (1) error-injection pipelines where
most subjects are unchanged per rep, and (2) secure-environment
runs where only a subset of subjects has new data.

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

### The `.child_valid()` function

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
| `test-cleangrowth.R` | 13 | ~80 | `cleangrowth()` API: adult integration (7 tests), child-adult spanning (3 tests) |
| `test-child-regression.R` | 8 | ~48 | Frozen counts, spot checks, cross-sample stability, Missing, HC, preload_refs, changed_subjids |
| `test-child-algorithms.R` | 24 | ~40 | CF rescue (modes + threshold cells), Evil Twins/OTL (1/2/3 unit errors + collateral), Error Load (threshold/mincount/constructed), HC boundary, cf_detail, parallel, GA correction (potcorr + near-potcorr), Birth EWMA2 (extreme/normal WT, extreme HT), LENGTHCM identity |
| `test-child-parameters.R` | 7 | ~15 | include.carryforward, sd.extreme, ewma_window, error.load params, imperial units, LENGTHCM |
| `test-child-edge-cases.R` | 12 | 24 | Single subject, sparse data, all-NA, mixed NA, SDE-Identical, negative agedays, HEADCM >3yr, extreme values, density mix, CF, deterministic |
| `test-adult-clean.R` | 198 | 198 | All 14 adult steps, 4 permissiveness levels, edge cases (via `cleanadult()` directly) |
| **Total** | **262** | **~405** | |

Additional test files not actively maintained:
`test-cdc.R` (not modified in v3.0.0); `test-utils.R` (10
failures — old API without `id` parameter, tests utility
functions `splitinput`, `recode_sex`, `longwide`, `simple_bmi`).

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

- [ ] **LENGTHCM → HEIGHTCM:** No measurement adjustment for
  the ~0.5–0.7 cm supine/standing difference. Known
  limitation; future update planned.
- [ ] **HEADCM CF rescue thresholds are HEIGHTCM-derived
  (placeholder).** `.cf_rescue_lookup()` reuses the HEIGHTCM
  `other` and `imperial` matrices for HEADCM, with an inline
  `# placeholder: use HT thresholds` comment. The 100K
  synthetic tracker population used to derive the current
  thresholds (see `__Pipeline/CF-exploration/cf-threshold-schemes.md`)
  did not include HC. HC-specific thresholds would require a
  dedicated derivation study on a synthetic HC population.
  Decision deferred; HT thresholds are a reasonable starting
  point given similar reference-standard structure.
- **`parallel = TRUE` requires installed package.** Will fail
  with `load_all()` only — workers need `system.file()` access
  to extdata. Install with `devtools::install_local(".")` first.

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

### Open (wrapper)

Deferred from walkthrough sessions 1–3 (2026-04-16):

- [ ] **`cleangrowth()` `@return` roxygen incomplete.** Roxygen
  lists only part of the returned data.table's columns. The
  actual return also includes `internal_id, subjid, agedays,
  sex, line, bin_exclude, tri_exclude`, conditional columns
  (`cf_status`/`cf_deltaZ` under `cf_detail=TRUE`, `ewma1_it1.*`
  under `debug=TRUE`), and internal working columns (`v`,
  `v_adult`, `measurement`, `age_years`, `newbatch`, `sd.corr`,
  `sd.orig_uncorr`, `z.orig`, `uncorr`, `potcorr`). Before
  updating the roxygen, decide which internal columns should
  be dropped from output vs. kept as part of the public
  contract, then align `@return` to the final schema.
- [ ] **Partial-run output row ordering.** `cleangrowth()` in
  partial-run mode rbinds changed-subject output to cached
  output and sorts by `internal_id`. Because `internal_id` is
  reassigned 1..K over only the changed-subject input and
  collides with 1..N in `cached_results`, the resulting
  whole-dataset row order is deterministic but semantically
  meaningless. Per-row data is correct. Acceptable as-is;
  if revisited, sort by user `id` (contract-guaranteed unique)
  or renumber `internal_id` after rbind.
- [ ] **Velocity reference tables not in `gc_preload_refs()`.**
  `tanner_ht_vel` is loaded once per `cleangrowth()` call;
  WHO velocity tables (`who_ht_vel_for_age`,
  `who_hc_vel_for_age`) are loaded once per batch inside
  `cleanchild()`. Both are tiny (< 1 kB each gzipped) and the
  `fread()` calls are fast, so correctness is unaffected.
  Extension: add the three tables to `gc_preload_refs()` and
  plumb optional parameters through `cleangrowth()` /
  `cleanchild()`; keep existing `fread()` fallbacks.

### Open (adult)

- [ ] **Deferred test gaps:** Error load with -5 exponent,
  UW scaling edge cases (very low/high UW).
- [ ] **Performance:** `setkey(df, subjid)` optimization
  deferred.
- [ ] **Deferred test gaps:** Error load with -5 exponent,
  UW scaling edge cases (very low/high UW).
- [ ] **Performance:** `setkey(df, subjid)` optimization
  deferred.

### Fixed (recent — condensed 2026-04-16)

See git history for full details on each fix.

- **Pre-walkthrough cleanup (2026-04-16):** Added `debug = FALSE`
  parameter to `cleangrowth()` and `cleanchild()`; gates 6
  `ewma1_it1.*` columns (created in Step 11 only when `debug=TRUE`,
  emitted in output only when `debug=TRUE`). Resolved deferred
  walkthrough items: Step 17 dead inner `if (count_exclude >= 1)`
  branch collapsed; adult dead `keeper_id` removed (2 locations);
  adult dead `"temp extraneous"` filter removed in `eval_2d_nonord()`;
  Step 15 Birth WT block now uses `.child_exc(param, "Traj")` for
  consistency with surrounding code (was 4 hardcoded literals).
- **Doc/dead-code cleanup (2026-04-16):** Deleted 6 orphaned
  man pages for removed extdata files. Removed dead `nnte`
  column assignment and unused NULL declarations (`sum_val`,
  `no_dup_val`, `no_outliers`, `no_bigdiff`, `nottoofar`,
  `nnte_full`). Fixed 3 stale comments (legacy file refs).
  Cleaned `_pkgdown.yml` (removed deleted functions/datasets).
  Updated vignettes (removed `prelim-infants-algorithm.Rmd`,
  fixed stale params in `configuration.Rmd`, `usage.Rmd`,
  `large-data-sets.Rmd`). Fixed `ewma()` roxygen (added
  missing `@param` tags, separated from `as_matrix_delta`).
  `as_matrix_delta` now internal (removed from exports).
  Regenerated all man pages via `roxygenise()`.
- **Legacy removal (2026-04-16):** Deleted
  `pediatric_clean_legacy.R`, `pediatric_support_legacy.R`,
  `R/deprec/`, `R/modified_source_code/`, 5 legacy extdata
  files. Removed parameters: `use_legacy_algorithm`,
  `prelim_infants`, `lt3.exclude.mode`, `ewma.exp`,
  `recover.unit.error`. Simplified `read_anthro()` (removed
  individual WHO txt path), `gc_preload_refs()` (2 closures,
  not 3), recentering (removed NHANES/derive string paths),
  dispatch (direct `cleanchild()` only), `var_for_par`
  (removed legacy functions).
- **Exclusion code work (2026-04-14–16):** Codes renamed to
  `Exclude-C-{Reason}` / `Exclude-A-{Reason}` format; param
  indicators removed; `Missing` → `Exclude-Missing`; adult
  `-N` round suffix removed; `.child_exc()` helper;
  `exclude.levels` updated with all adult codes.
- **CF rescue thresholds (2026-04-14):** Age/interval/param-
  specific lookup tables replace fixed thresholds.
- **internal_id (2026-04-12–16):** Sequential integer for all
  internal sorting; user `id` preserved untouched.
- **Walkthrough cleanups (2026-04-15–16):** 30+ items: `cat()`
  → `message()`, dead code removal, `batch_size` parameter,
  correctness fixes.
- **Earlier fixes:** BIV preterm threshold, error load bug
  (hardcoded .4 → parameter), missing data bug (duplicate
  `valid()` definitions), CF rescue re-inclusion, z-score
  blending bug, parallel processing fixes, outer batching
  wrapper fix.

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
  prevents vignette building. Stale parameter references cleaned
  up (2026-04-16): removed `prelim-infants-algorithm.Rmd`,
  updated `configuration.Rmd`, `usage.Rmd`, `large-data-sets.Rmd`.
  Still need to verify vignettes build cleanly or remove them.
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
- [x] **`R/deprec/` and `R/modified_source_code/` removed
  (2026-04-16):** Entire directories deleted.
- [x] **`cat()`/`print()` in child algorithm (2026-04-15):**
  All `cat()`/`print()` calls converted to `message()`.

### Medium (pre-submission polish)

- [ ] **Verify `LICENSE` file exists:** DESCRIPTION says
  `MIT + file LICENSE`. Confirm plain `LICENSE` file (not
  just `LICENSE.md`) is at package root.
- [x] **Orphaned man pages cleaned up (2026-04-16):** Deleted
  `bmianthro.Rd`, `lenanthro.Rd`, `weianthro.Rd`,
  `nhanes-reference-medians.Rd`, `growth_cdc_ext.Rd`,
  `fentlms_forz.Rd`. Regenerated `cleangrowth.Rd`, `ewma.Rd`,
  `fenton2025_ms_lookup_smoothed.Rd` via `roxygenise()`.
  `as_matrix_delta` removed from exports (internal helper).
- [x] **Roxygen return value for `cleangrowth()`:** Updated to
  describe data.table with all columns. (Fixed prior to
  2026-04-16; confirmed in generated `cleangrowth.Rd`.)
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

- **Do not create additional `.child_valid()` definitions.**
  `.child_valid()` is defined in `child_clean.R`. Do not
  duplicate it in other files — R loads all `.R` files
  alphabetically and the later definition wins silently.
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
