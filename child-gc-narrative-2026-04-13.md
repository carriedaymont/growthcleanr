# Child Growthcleanr Algorithm Narrative

Technical specification for `child_clean.R` (v3.0.0) and
supporting functions.

Date of Narrative: 2026-03-18 (updated 2026-04-13). Initial
draft by: Claude (Opus), with Carrie Daymont.

Date of Most Recent Code: 2026-04-13

This document describes what the current R code does. It is
intended as a debugging aid and a long-term technical reference.

**Note:** Line number references have been removed (2026-04-13)
because they become stale with each code edit. Use function
names and step numbers to locate code.

---

## How to Read This Document

This narrative follows the same general structure as the adult
growthcleanr narrative review. Each algorithm step has:

- A summary table (scope, prior/next step, exclusion codes)
- Overview
- Key terms and variable names
- Logic and implementation
- Rationale for selected decisions

### Code review checklist (applied to each step)

1. **Comments and narrative match code:** No stale comments
 describing removed logic; narrative accurately reflects
 what the code does
2. **Exact boundaries:** `<` vs `<=`, `>` vs `>=`, and
 exact numeric tolerances (e.g., 0.4 vs 0.41) verified
 against spec/Stata
3. **Uninitialized or unnecessary variables:** No variables
 used before assignment; no leftover variables from prior
 refactors
4. **Efficiency:** Caching opportunities, unnecessary
 merges, redundant computations, operations that could
 be vectorized
5. **Edge cases:** 0, 1, or 2 valid rows; empty groups;
 all-NA values — does code handle gracefully or crash?
6. **`valid()` call correctness:** Which `include.*` flags
 are passed? Consistent with step's purpose?
7. **`by`-group correctness:** `by = subjid` vs
 `by = .(subjid, param)` etc. — wrong grouping silently
 produces wrong results
8. **Sort order assumptions:** Does code assume a particular
 sort? Is that guaranteed by a prior `setkey`/`setorder`?
9. **data.table reference semantics:** `:=` on `data.df`
 (intended) vs. accidentally modifying a copy or joined
 table?
10. **Factor level issues:** Exclusion codes assigned must
 exist in `exclude.levels`; assigning an unlisted string
 to a factor silently produces NA
11. **Parameter scope:** Does the step handle all 3 params
 (HT, WT, HC) or only some? Is HC exclusion intentional
 and documented?
12. **Interaction with later steps:** Variables stored for
 downstream use — correctly scoped in `data.df`?

Before the step-by-step walkthrough, this document covers:

- Key Concepts (how the algorithm works at a high level)
- Architecture (preprocessing vs. main algorithm)
- Data Requirements and Input Format
- Z-Score Infrastructure
- Working Dataframe and Output
- Variable Glossary
- Output Format
- Configurable Parameters
- Complete List of Steps

***Changes from pediatric:*** There are multiple changes from the previous pediatric algorithm and the current child algorithm. These differences are discussed in a separate document. This document focuses on what the algorithm does now.

---

## Key Concepts

***Data retention:*** Although we refer to values as being excluded, growthcleanr does not exclude any values from the dataset returned to users. It adds a column to user' input recommending inclusion or exclusion. It does exclude values from further consideration in cleaning. For example, a value of 690kg will be excluded early in the algorithm and is ignored thereafter. As shorthand, we refer to this as excluding a value. 

***Three growth parameters:*** The child algorithm cleans weight (WEIGHTKG), height/length (HEIGHTCM/LENGTHCM), and head circumference (HEADCM). LENGTHCM is relabeled to HEIGHTCM — only the `param` value changes; the measurement value is **not** altered. No adjustment is made for the ~0.5–0.7 cm difference between supine length and standing height. After relabeling, z-scores are calculated identically to HEIGHTCM. Code comments and 725–727 acknowledge this as a known limitation to be addressed in a future update. HEADCM is only cleaned through age 3 years (marked "Exclude-Not-Cleaned" for agedays > 3 × 365.25); z-score data for HC is available only through age 5 years (marked "Exclude-Missing" for agedays ≥ 5 × 365.25).

***CSD z-scores, not measurements:*** The child algorithm bases most decisions on z-scores, which standardize measurements by age and sex. The algorithm uses the Conditional Standard Deviation (CSD) method for calculating z-scores (also referred to as SD-scores), not the standard LMS method. [ADD CDC REFERENCE] The LMS accounts for significant skewness in the distribution of weights, which is often useful, but in cleaning it reduces the ability of z-scores to distinguish large absolute differences in weight for heavier children, reducing cleaning accuracy. Several variations of z-scores (SD-scores) are used. These are described further below and are listed in the "### Z-score variants" section below. The two primary variations are sd.orig (the unaltered SD score) and tbc.sd. tbc.sd is calculated by subtracting a recentering value (described further below) from sd.orig.

***Dual z-scores for prematurity correction:*** One z-score variant is corrected z-scores. For subjects whose data are consistent with prematurity or SGA, a corrected z-score is calculated using values derived from the Fenton 2025 percentile reference. [ADD Fenton 2025 REF ] Please note that we do not use the Fenton LMS values that would be used in clinical care. Instead, we extracted approximate standard deviations from the plotted percentiles. Corrected z-scores are used as an extra check to prevent exclusion of values that are plausible for a baby with very low birth weight. 

***Recentering:*** After calculating z-scores, the algorithm subtracts a population-level median z-score for each parameter, sex, and age to create "to-be-cleaned" SD scores (tbc.sd). Because median z-score changes with age, using an unrecentered z-score can lead to bias in exclusions when comparing z-scores across long time periods. The recentering medians come from a fixed reference file (rcfile-2023-08-15_format.csv.gz) that was derived from data from outpatient records from the Children's Hospital of Philadelphia.

***Designated Other Parameter (DOP):*** Many steps compare a measurement against a different growth parameter for the same
subject. The designated other parameter (DOP) for weight is height; the DOP for height is weight; the DOP for HC is height. DOP comparisons are used in temporary SDE resolution (Step 5), CF rescue, EWMA steps, final SDE resolution (Step 13), and pairs/singles evaluation (Step 19). In Step 13, DOP medians are computed for one-day SDE resolution (cross-parameter comparison) and DOP is used as a secondary sort key for selecting which SDE to keep.

***Ordered steps:*** Like the adult algorithm, the child algorithm processes steps in a specific order, removing more
extreme problems first and getting more refined later. Same-day extraneous values (described below) are handled in repeated steps because of bidirectional impact of same-day selection and other exclusions. Some errors cause over- or under-detection of other errors, so they must be excluded first. The order of steps is central to the algorithm.

***Iterations within steps:*** Several steps (Evil Twins, EWMA1, EWMA2) are iterative. If more than one value is a candidate for exclusion, they select the most extreme candidate and then re-evaluate for further exclusions. This continues until no candidates remain or too few values are left to evaluate. This is a central part of the algorithm and reduces the likelihood of erroneous values leading to exclusion of true values.

***Temporary vs. permanent exclusions:*** Same-day extraneous (SDE) values are multiple values for the same parameter on the same day. growthcleanr selects the most plausible SDE for inclusion in the final dataset. are first flagged temporarily (Step 5) and later
resolved permanently (Step 13). Carried-forward values may be
"rescued" and returned to Include status based on z-score
similarity criteria. These temporary mechanisms allow the
algorithm to refine decisions as the dataset gets cleaner.

***No separate parameter pipelines:*** Unlike the adult
algorithm, which processes height and weight in completely
separate working dataframes, the child algorithm keeps all
parameters in a single data.table (data.df) with parameter
tracked in a `param` column. Steps that operate on a single
parameter use `by = .(subjid, param)` groupings or explicit
filtering. This design supports the cross-parameter
comparisons (DOP) that are central to the child algorithm.

***CF rescue:*** The child algorithm has a detailed
carried-forward rescue system (Step 6) that attempts to
distinguish genuinely repeated measurements from data entry
artifacts. Rescue criteria are based on z-score similarity
between the CF value and its originator, whether the
measurement is in whole/half imperial units, the subject's
age, and the length of the CF string. Rescued CFs are returned
to "Include" status with the rescue reason stored in a
separate `cf_rescued` column.

***EWMA (Exponentially Weighted Moving Average):*** The
central mechanism for detecting implausible values. For each
measurement, the EWMA computes a weighted average of the
subject's other measurements for the same parameter, giving
exponentially more weight to temporally closer values. The
deviation from EWMA (dewma) indicates how far a value is from
what we would predict. Default exponent is -1.5 (child) vs.
-5 (adult). With -1.5, the weighting is less dominated by the
nearest neighbor, spreading influence across a wider time
range — appropriate for pediatric data where growth trajectories
change more rapidly.

***The `valid()` function:*** A critical gatekeeper that
determines which rows participate in each step. It returns a
logical vector indicating which rows are currently "valid"
(eligible for processing). Its `include.*` flags control
whether temporarily excluded categories (temp SDEs, permanent
SDEs, carried-forward values) are included:

- `valid(df)` — only rows with `exclude` that does not start
 with "Exclude"
- `valid(df, include.temporary.extraneous = TRUE)` — also
 includes "Exclude-C-Temp-Same-Day"
- `valid(df, include.extraneous = TRUE)` — also includes
 "Exclude-C-{WT|HT|HC}-Extraneous"
- `valid(df, include.carryforward = TRUE)` — also includes
 "Exclude-C-{WT|HT|HC}-CF"

The valid function works by checking the text of the exclude
column rather than using factor ordering. It uses `grepl("^Exclude", ...)`
to identify excluded rows, then checks for specific patterns to
re-include temporarily excluded categories. Non-Exclude codes like
"Include", "Unit-Error-High", "Unit-Error-Low", "Swapped-Measurements"
pass through.

***Age-dependent id tiebreaking:*** When resolving same-day
duplicates, the child algorithm uses an age-dependent rule:
at birth (agedays == 0), keep the lowest id (earliest
measurement, before postnatal fluid shifts and interventions);
at all other ages, keep the highest id (later measurement,
which may represent a more careful re-measurement). This
differs from the adult algorithm, which always keeps the
highest id.

---

## Architecture: Preprocessing vs. Main Algorithm

The child algorithm is split across two functions:

### `cleangrowth()` — preprocessing (child_clean.R)

Handles all data setup before the main cleaning loop:

1. Input validation and data.table construction
2. Imperial-to-metric conversion (HEIGHTIN → HEIGHTCM,
 WEIGHTLBS → WEIGHTKG) — done by modifying `v` and
 changing `param`
3. LENGTHCM → HEIGHTCM reclassification
4. Age cutpoint split (pediatric vs. adult data)
5. Z-score calculation:
 a. CSD z-scores from WHO and CDC references
 b. Age-dependent blending of WHO/CDC z-scores
6. Gestational age correction (Step 2b) for potcorr subjects
7. SD-score recentering to create tbc.sd and ctbc.sd
8. Missing value identification
9. Batching and dispatch to `cleanchild()`
 (see "Batching and Dispatch" section below for details,
 including a known bug where the outer wrapper is
 currently defeated)
10. Adult data dispatch to `cleanadult()`
11. Result merging and output assembly

### `cleanchild()` — main algorithm (child_clean.R)

Processes one batch of pediatric data through all cleaning
steps:

- Early Step 13: SDE-Identicals
- Step 5: Temporary SDEs
- Step 6: Carried Forwards (with rescue)
- Step 7: BIV (Absolute and Standardized)
- Step 9: Evil Twins
- Step 11: Extreme EWMA (EWMA1)
- Step 13: Final SDE Resolution
- Step 15: Moderate EWMA (EWMA2)
- Step 16: Birth HT/HC (EWMA2 variant for birth values)
- Step 17: Height/HC Velocity Checks
- Step 19: Pairs and Singles Evaluation
- Step 21: Error Load
- Step 22: Output Preparation

Support functions are defined in `child_clean.R` (after
`cleanchild`) and `utils.R`:

- `valid()` — row eligibility filter
- `temporary_extraneous_infants()` — temp SDE resolution
- `calc_otl_evil_twins()` — evil twins OTL calculation
- `calc_and_recenter_z_scores()` — z-score recalculation
 for CF rescue
- `ewma()` — EWMA computation with matrix operations
- `ewma_cache_init()` / `ewma_cache_update()` — incremental
 EWMA cache for iterative steps
- `get_dop()` — designated other parameter lookup
- `read_anthro()` — growth reference table loader

---

## Batching and Dispatch

### Purpose

The outer batching wrapper (added by Chris Palmisano)
divides subjects into groups for memory management. Without
batching, preprocessing (z-score calculation, recentering,
reference table merges) for the full dataset must fit in
memory at once. For large datasets (tens of thousands or
millions of subjects), this can cause memory pressure. The
wrapper processes one batch at a time, then combines results.

### Subject-level batching

Batching occurs at the **subject level**, not the row level.
All of a given subject's observations — across all parameters
(WEIGHTKG, HEIGHTCM, HEADCM) and all ages — are assigned to
the same batch. This is essential because the cleaning
algorithm uses each subject's full longitudinal trajectory
when evaluating individual measurements (EWMA, velocity
checks, carried-forward detection, etc.).

### Interaction with the child/adult age split

After batching by subject, `cleangrowth()` splits each
batch's data by `adult_cutpoint` (default 20 years, range
18–20) into pediatric (`data.all`) and adult (`data.adult`)
subsets. A subject whose observations span the cutpoint will
have some rows sent to `cleanchild()` and others to
`cleanadult()`. This is intentional: the child and adult
algorithms are independent pipelines with no interaction
between them. A subject's child measurements are cleaned
using only their other child measurements, and their adult
measurements are cleaned using only their other adult
measurements.

Both sets of rows for such a subject are guaranteed to be in
the same batch because batching is done on `subjid` before
the age split. The result assembly at the end of each batch
iteration recombines the child and adult results for the
batch's subjects, and `rbindlist()` combines across batches.

### Two batching systems in the code

The code contains two layered batching systems with different
purposes:

**System 1 — Outer wrapper (Chris's addition, for memory
management):**
- In `child_clean.R`
- Creates batches of up to 2,000 subjects using
 `dplyr::mutate(batch = (row_number() - 1) %/% 2000 + 1)`
- Iterates with a `for` loop; each iteration filters
 `data.all.ages` to the current batch's `subjid` values
- Collects results in `results_list`, then combines with
 `rbindlist()` at the end

**System 2 — Inner batching (for parallelism):**
- In `child_clean.R`
- Assigns random batch numbers to subjects in the
 already-filtered `data.all` using
 `sample(num.batches, ...)`
- When `parallel = FALSE`, `num.batches = 1`, making this a
 no-op: all subjects in the outer batch go to a single
 `cleanchild()` call
- When `parallel = TRUE`, dispatches via
 `ddply(..., .parallel = TRUE)` to process batches in
 parallel across CPU cores

The inner batching subdivides whatever data the outer wrapper
provides into parallel work units. With `parallel = FALSE`,
it has no effect. The adult algorithm has its own inner
parallel batching system, structured
identically.

### Outer wrapper bug (FIXED)

The outer batching wrapper previously had a bug where
`data.all` and `data.adult` were overwritten with ALL
subjects (not just the current batch's subjects), defeating
the memory management purpose of the wrapper. Fixed by:

1. Batch filter now correctly limits `data.all` and
 `data.adult` to the current batch's subjects
2. A batch-level reference (`data.batch`) is used for
 result assembly instead of `data.all.ages`

Additionally, several batch-invariant operations have been
moved before the loop:
- `exclude.levels` definition
- Tanner/WHO velocity reference file reads

These were not bugs, but caused redundant I/O and
computation per batch.

### Batch size

The current batch size of 2,000 subjects is hard-coded.
Chris tested several sizes to find a good default. This
should be parameterized (e.g., a `batch_size` argument to
`cleangrowth()` with default 2000) for two reasons:

- Efficiency improvements since the original testing may
 shift the optimal batch size
- In HPC environments with large memory, much larger batches
 (or no batching) may be preferable for datasets with
 millions of subjects, avoiding per-batch overhead

### Parallel processing — known limitation

The inner parallel batching system (`parallel = TRUE`) is
currently disabled because it produces incorrect results on
macOS. The root cause has not been fully investigated but is
most likely the well-known incompatibility between
`data.table`'s internal threading and R's forked parallel
processes (`doParallel`/`foreach` with `type = "FORK"`).
macOS has additional restrictions on `fork()` with
multithreaded processes that can silently corrupt shared
state.

Potential fixes to investigate:
- Use `type = "PSOCK"` socket clusters instead of forked
 processes
- Replace `ddply`/`foreach` with `data.table`'s own
 `setDTthreads()` parallelism
- Use `future`/`furrr` framework with `plan(multisession)`

This is a future priority — parallel processing would
significantly reduce runtime on large datasets, particularly
on macOS where the outer batching wrapper's sequential
processing is the only option.

---

## Data Requirements and Input Format

### Required input columns

`cleangrowth()` accepts individual vectors (not a dataframe):

| Parameter | Type | Description |
|-----------|------|-------------|
| `subjid` | any | Subject identifier |
| `param` | character | `"WEIGHTKG"`, `"WEIGHTLBS"`, `"HEIGHTCM"`, `"HEIGHTIN"`, `"LENGTHCM"`, or `"HEADCM"` |
| `agedays` | numeric | Age in days at measurement |
| `sex` | any | `0`, `"m"`, or `"M"` for male; `1`, `"f"`, or `"F"` for female |
| `measurement` | numeric | Recorded value in units specified by `param` |
| `id` | numeric or character (optional) | Unique row identifier; if NULL, generated as `1:N` |

### Input handling

- `measurement == 0` is replaced with `NaN` (treated as
 missing) in the internal variable `v`
- Imperial measurements are converted to metric before
 processing: HEIGHTIN multiplied by 2.54 and relabeled
 HEIGHTCM; WEIGHTLBS divided by 2.2046226 and relabeled
 WEIGHTKG
- LENGTHCM is relabeled to HEIGHTCM; only the
 `param` value changes — the measurement value is not
 adjusted. No supine-to-standing correction is applied.
 Z-scores are calculated identically to HEIGHTCM after
 relabeling. This is a known limitation (see code comments, 725–727).
- `sex` is recoded to integer: 0 = male, 1 = female
- `agedays` is cast to integer
- Data are split by `adult_cutpoint` (default 20 years):
 rows with `agedays < cutpoint * 365.25` go to the child
 algorithm; rows at or above go to the adult algorithm

### Sorting

The primary sort order is:

```
setkey(data.df, subjid, param, agedays, id)
```

This sort is critical for deterministic results. The `id`
field breaks ties when multiple measurements share the same
`subjid`, `param`, and `agedays` (same-day duplicates). If
`id` is not included in the sort key, SDE resolution order
is undefined and results may differ between sequential and
parallel runs.

### Internal identifiers

- `line` — row number in original input order; used to
 restore output order
- `index` — sequential row number after sorting; used as
 internal key within `cleanchild()`
- `id` — user-provided or auto-generated; used as tiebreaker
 in sorting and SDE resolution
- `v` — working measurement value (metric, with 0 → NaN)
- `v.orig` — copy of `v` before any transformations (unit
 error recovery, swapped parameter correction)

---

## Z-Score Infrastructure

The child algorithm computes and uses multiple z-score
variants. Understanding which is used where is essential for
debugging.

### Calculation method: CSD (Conditional Standard Deviation)

All z-scores in the cleaning algorithm use the CSD method:

```
If measurement < M: sd = (measurement - M) / csd_neg
If measurement >= M: sd = (measurement - M) / csd_pos
```

Where:
- `M` = age/sex-specific median from reference table
- `csd_pos` = half the distance from M to the value at
 z = +2, precomputed in reference tables
- `csd_neg` = half the distance from M to the value at
 z = -2, precomputed in reference tables

This is NOT the standard LMS z-score used in clinical
practice. CSD z-scores are intentionally more sensitive to
extreme high values (especially for weight), which helps
the algorithm detect implausible measurements.

### Reference tables

| Table file | Source | Used for |
|------------|--------|----------|
| `growthfile_who.csv.gz` | WHO 2006 | HT/WT age < 5y; all HEADCM |
| `growthfile_cdc_ext_infants.csv.gz` | CDC 2000 extended | HT/WT age ≥ 2y |
| `fentlms_foraga.csv.gz` | Fenton | Weight → est. gestational age |
| `fentlms_forz.csv.gz` | Fenton | GA-corrected z-scores |
| `rcfile-2023-08-15_format.csv.gz` | Derived | Recentering medians |

### Z-score variants

| Variable | Source | Description |
|----------|--------|-------------|
| `sd.orig_who` | WHO CSD | WHO-only CSD z-score |
| `sd.orig_cdc` | CDC CSD | CDC-only CSD z-score |
| `sd.orig` | Blended | Age-blended WHO/CDC CSD z-score (see below) |
| `sd.orig_uncorr` | Copy | Copy of `sd.orig` before GA correction; used in CF rescue (Step 6) and BIV (Step 7) |
| `sd.corr` | GA-corrected | Fenton-corrected z-score for potcorr subjects; equals `sd.orig` for others |
| `tbc.sd` | Recentered | `sd.orig - sd.median`; primary score used by most algorithm steps |
| `ctbc.sd` | Recentered corrected | `sd.corr - sd.median`; used alongside `tbc.sd` in Evil Twins and EWMA steps |
| `final_tbc` | Output | `ctbc.sd` for potcorr subjects, `tbc.sd` for others |

Intermediate z-score rounding has been removed.
Previously all z-score variables were rounded to 0.001
using Stata-style rounding (half away from zero) for
cross-platform validation; this is no longer needed.

### Age blending — main z-score (`sd.orig`)

No intermediate rounding is applied — `sd.orig_who` and
`sd.orig_cdc` are blended at full precision.

| Age range | Parameter | Formula |
|-----------|-----------|---------|
| < 2 years | any | `sd.orig = sd.orig_who` |
| 2–5 years | HT, WT | `(sd.orig_cdc × (age-2) + sd.orig_who × (5-age)) / 3` |
| > 5 years | HT, WT | `sd.orig = sd.orig_cdc` |
| any age | HEADCM | `sd.orig = sd.orig_who` |

Note: The blending window is 2–5 years with divisor 3.
At exactly 2 years, the formula yields 100% WHO. At exactly
5 years, it yields 100% CDC.

***Code-verified detail:*** The code assigns WHO for
`ageyears < 2` and CDC for `ageyears > 5` (strict
inequality). For ages exactly 2.0 or 5.0, the smooth_val
condition (`>= 2 & <= 5`) applies, so they go through
the blending formula. This means ages exactly at the
boundaries use the blending formula, which gives 100%
WHO at 2.0 and 100% CDC at 5.0.

### Age blending — corrected z-score (`sd.corr`)

For potcorr subjects, `sd.corr` is calculated in two
stages:

**Stage 1 — Compute `sd.c` (the corrected z-score before
smoothing into original):**

First, corrected WHO and CDC z-scores (`sd.c_who`,
`sd.c_cdc`) are calculated using corrected age (`cagedays`)
with an adjustment for standing vs. supine position when
chronological age > 730 days but corrected age ≤ 730 days
(+0.8 cm for WHO, +0.7 cm for CDC; the relevant code section).

Then `sd.c` is assigned using WHO/CDC blending with the
same 2–5 year window as `sd.orig`:
- HEADCM: WHO at all ages
- HT/WT age ≤ 2: WHO
- HT/WT age 2–5: weighted blend `(who × (5-age) + cdc × (age-2)) / 3`
- HT/WT age ≥ 5: CDC

For potcorr subjects age ≤ 2 with corrected WHO available
and post-menstrual age ≥ 350 days, `sd.corr` is set to
`sd.c_who`.

Then Fenton is preferred over corrected WHO for ages ≤ 2:
if Fenton (`unmod_zscore`) is available, it overwrites
`sd.c`; corrected WHO is used only as fallback.

**Stage 2 — Smoothing corrected into original:**

| Age range | Formula |
|-----------|---------|
| ≤ 2 years (potcorr) | `sd.corr = sd.c` (fully corrected) |
| 2–4 years | `(sd.orig × (4 - age) + sd.c × (age - 2)) / 2` |
| > 4 years | `sd.corr = sd.orig` (no correction) |

Note: This smoothing window is 2–4 years with divisor 2,
which is narrower than the main WHO/CDC blending window
(2–5 years, divisor 3).

### Age blending — `calc_and_recenter_z_scores()` (CF rescue)

This helper function, used during carry-forward rescue in
Step 6, recalculates z-scores for modified measurement
values. It should use the same WHO/CDC blending formula
as the main z-score calculation (2–5 year window).

***What the code currently does:***

| Age range | Parameter | Formula |
|-----------|-----------|---------|
| < 2 years | any | WHO only |
| 2–5 years | HT, WT | `(cdc × (age-2) + who × (5-age)) / 3` |
| > 5 years | HT, WT | CDC only |
| any age | HEADCM | WHO only |

This now matches the main z-score blending formula.

***BUG FIXED (2026-03-20):*** The CDC-only cutoff was
previously `>= 4 years`, which conflated the WHO/CDC
blending window (2–5 years) with the corrected/uncorrected
smoothing window (2–4 years). Fixed to `> 5` to match the
main z-score calculation.

### Recentering

After z-score calculation, the algorithm subtracts age-,
sex-, and parameter-specific median z-scores to produce
tbc.sd:

```
tbc.sd = sd.orig - sd.median
ctbc.sd = sd.corr - sd.median
```

The child algorithm always uses a fixed reference file
(`rcfile-2023-08-15_format.csv.gz`) for recentering medians,
regardless of dataset size. The legacy algorithm uses NHANES
medians for small datasets (< 5000 observations) or
derives medians from the input for larger datasets.

The recentering median file is indexed by param, sex, and
agedays. It is merged into the data by a rolling join on
these keys.

---

## Working Dataframe and Output

### How exclusions work

The algorithm never removes rows from the output. Every
input row receives an `exclude` code — either `"Include"`
or an exclusion code indicating why and at which step the
value was excluded.

Internally, most steps operate on the full `data.df`
data.table but filter to currently valid rows using the
`valid()` function. When the narrative says a value is
"excluded," this means:

1. The value's `exclude` column is set to the appropriate
 exclusion code
2. Subsequent calls to `valid()` will exclude this row
 from processing
3. The value remains in `data.df` with its exclusion code

There is also a softer mechanism: the **temporary SDE**
flag. Values marked `"Exclude-C-Temp-Same-Day"`
are excluded by default `valid()` calls but can be included
by passing `include.temporary.extraneous = TRUE`. This
allows steps to optionally include or exclude temp SDEs.
Permanent SDE resolution (Step 13) replaces temp SDE codes
with final codes.

### Key difference from adult algorithm

The adult algorithm copies rows into separate working
dataframes (`h_subj_df`, `w_subj_df`) that shrink as values
are excluded. The child algorithm keeps all rows in a single
`data.df` and uses the `exclude` column plus `valid()` to
control which rows participate. Rows are never physically
removed from `data.df` during processing (except for a
temporary removal of SDE-Identicals during CF string
detection, which are re-added after).

---

## Variable Glossary

### Input and identification

| Variable | Type | Created | Description |
|----------|------|---------|-------------|
| `line` | integer | Preprocessing | Original row number; used to restore output order |
| `id` | numeric/char | Input or generated | User-provided unique row id; tiebreaker for sorting |
| `index` | integer | cleanchild | Sequential row number after sorting; internal key |
| `subjid` | factor | Input | Subject identifier |
| `param` | character | Input | `"WEIGHTKG"`, `"HEIGHTCM"`, or `"HEADCM"` (after conversion) |
| `agedays` | integer | Input | Age in days |
| `sex` | integer | Input | 0 = male, 1 = female |
| `v` | numeric | Preprocessing | Working measurement value (metric; 0 → NaN) |
| `v.orig` | numeric | Step 5 | Copy of `v` before transformations |

### Z-score variables

| Variable | Type | Created | Description |
|----------|------|---------|-------------|
| `sd.orig_who` | numeric | Preprocessing | WHO-only CSD z-score |
| `sd.orig_cdc` | numeric | Preprocessing | CDC-only CSD z-score |
| `sd.orig` | numeric | Preprocessing | Blended WHO/CDC CSD z-score |
| `sd.orig_uncorr` | numeric | Preprocessing | Copy of sd.orig before GA correction |
| `sd.corr` | numeric | Step 2b | GA-corrected CSD z-score |
| `sd.median` | numeric | Recentering | Population median z-score for recentering |
| `tbc.sd` | numeric | Recentering | Recentered blended z-score (`sd.orig - sd.median`) |
| `ctbc.sd` | numeric | Recentering | Recentered corrected z-score (`sd.corr - sd.median`) |
| `final_tbc` | numeric | Output | ctbc.sd for potcorr subjects, tbc.sd for others |

### Prematurity correction variables

| Variable | Type | Created | Description |
|----------|------|---------|-------------|
| `potcorr` | logical | Step 2b | TRUE if subject may need GA correction (subject-level) |
| `potcorr_wt` | logical | Step 2b | TRUE for the qualifying first weight row (row-level) |
| `uncorr` | integer | Step 2b | 1 if correction was reverted (made z-scores more extreme) |
| `fengadays` | numeric | Step 2b | Estimated gestational age in days (from Fenton) |
| `pmagedays` | numeric | Step 2b | Post-menstrual age in days |
| `cagedays` | numeric | Step 2b | Corrected age in days |

### Exclusion and status variables

| Variable | Type | Created | Description |
|----------|------|---------|-------------|
| `exclude` | factor | Preprocessing | Exclusion code or "Include" |
| `cf_rescued` | character | Step 6 | CF rescue reason; empty string for non-rescued rows |
| `nnte` | logical | Preprocessing | Always FALSE in current code (legacy optimization removed) |

### Carried-forward variables (Step 6)

| Variable | Type | Created | Description |
|----------|------|---------|-------------|
| `cf_binary` | logical | Step 6 | TRUE if row is a carried-forward value |
| `originator` | logical | Step 6 | TRUE for Include values whose next value is CF |
| `cf_string_num` | integer | Step 6 | Sequential string number per subject-param |
| `originator_z` | numeric | Step 6 | sd.orig_uncorr of the originator |
| `cf_string_length` | integer | Step 6 | Number of CFs in the string |
| `wholehalfimp` | logical | Step 6 | TRUE if measurement is in whole/half imperial units |

### EWMA variables (internal, not in output)

| Variable | Type | Description |
|----------|------|-------------|
| `ewma.all` | numeric | EWMA excluding only the current observation |
| `ewma.before` | numeric | EWMA excluding current + prior observation |
| `ewma.after` | numeric | EWMA excluding current + next observation |
| `dewma.all` | numeric | tbc.sd - ewma.all |
| `dewma.before` | numeric | tbc.sd - ewma.before |
| `dewma.after` | numeric | tbc.sd - ewma.after |

---

## Output Format

`cleangrowth()` returns a data.table with all original input
columns plus additional columns. The columns returned from
`cleanchild()` are:

| Column | Type | Description |
|--------|------|-------------|
| `id` | original | User-provided row identifier |
| `line` | integer | Original row order |
| `exclude` | factor→character | Exclusion code or "Include" |
| `param` | character | Growth parameter |
| `cf_rescued` | character | CF rescue reason (empty if not rescued) |
| `sd.orig_who` | numeric | WHO CSD z-score (if available) |
| `sd.orig_cdc` | numeric | CDC CSD z-score (if available) |
| `sd.orig` | numeric | Blended CSD z-score |
| `tbc.sd` | numeric | Recentered blended z-score |
| `ctbc.sd` | numeric | Recentered corrected z-score |
| `final_tbc` | numeric | ctbc.sd for potcorr, tbc.sd for others |

After merging with adult results, `cleangrowth()` also
includes checkpoint diagnostic columns (sd.corr, potcorr,
uncorr, sd.orig_uncorr) when the child algorithm was used.

### Complete list of exclusion codes (child algorithm)

| Code | Step | Param | Description |
|------|------|-------|-------------|
| `Include` | — | All | Value passes all checks |
| `Exclude-Missing` | Init | All | Measurement is NA, NaN, or agedays < 0; also HC ≥ 5y |
| `Exclude-Not-Cleaned` | Init | HEADCM | HC with agedays > 3 × 365.25 |
| `Unit-Error-High` | 3 | HT, WT | Unit error detected (value too high); corrected in-place |
| `Unit-Error-Low` | 3 | HT, WT | Unit error detected (value too low); corrected in-place |
| `Swapped-Measurements` | 4 | HT, WT | Height and weight appear swapped; corrected in-place |
| `Exclude-C-Temp-Same-Day` | 5 | All | Temp SDE (may persist if not resolved by Step 13) |
| `Exclude-C-{WT\|HT\|HC}-CF` | 6 | All | Value identical to prior-day value (not rescued) |
| `Exclude-C-{WT\|HT\|HC}-BIV` | 7 | All | Outside absolute or standardized biological limits |
| `Exclude-C-{WT\|HT\|HC}-Evil-Twins` | 9 | All | Adjacent extreme value pair/group |
| `Exclude-C-{WT\|HT\|HC}-Traj-Extreme` | 11 | All | Extreme EWMA outlier |
| `Exclude-C-{WT\|HT\|HC}-Identical` | 13 | All | Same-day duplicate with identical value |
| `Exclude-C-{WT\|HT\|HC}-Extraneous` | 13 | All | Same-day extraneous (resolved by EWMA, all extreme, or single day) |
| `Exclude-C-{WT\|HT\|HC}-Traj` | 15/16 | All | Moderate EWMA outlier (all sub-rules collapsed) |
| `Exclude-C-{WT\|HT\|HC}-Abs-Diff` | 17 | HT, HC | Height/HC velocity exceeded allowed limits |
| `Exclude-C-{WT\|HT\|HC}-Pair` | 19 | All | 2 measurements, one excluded |
| `Exclude-C-{WT\|HT\|HC}-Single` | 19 | All | Single measurement excluded |
| `Exclude-C-{WT\|HT\|HC}-Too-Many-Errors` | 21 | All | Error ratio exceeds threshold |

Note: Exclusion codes are not parameter-specific — the param
is in the data row (e.g., `Exclude-C-CF` applies to WT, HT, and HC).
The `exclude.levels` list also includes legacy codes from the
original pediatric algorithm for backward compatibility with
datasets cleaned by older versions.

### CF rescue reason codes (cf_rescued column)

| Code | Meaning |
|------|---------|
| `""` (empty) | Not a carried-forward value, or CF not rescued |
| `Rescued` | Single CF, z-score difference < 0.05 |
| `Rescued-Imperial` | Single CF, z-score diff < 0.1, whole/half imperial |
| `Rescued-Adol` | ≥2 CFs, teen age, z-score diff < 0.05 |
| `Rescued-Adol-Imperial` | ≥2 CFs, teen age, z-score diff < 0.1, imperial |

Rescued CFs have their `exclude` column set back to
`"Include"` — the rescue reason is only preserved in the
`cf_rescued` column.

---

## Configurable Parameters

| Parameter | Default | Used in | Description |
|-----------|---------|---------|-------------|
| `recover.unit.error` | FALSE | Step 3 | Attempt to identify and correct unit errors |
| `sd.extreme` | 25 | Step 7 | SD-score cutoff for extreme exclusion |
| `z.extreme` | 25 | Step 7 | Z-score cutoff for extreme exclusion |
| `error.load.mincount` | 2 | Step 21 | Min exclusions before evaluating error load |
| `error.load.threshold` | 0.5 | Step 21 | Error ratio above this excludes all |
| `sd.recenter` | NA | Recentering | Recentering method (child always uses fixed file) |
| `include.carryforward` | FALSE | Step 6 | If TRUE, skip CF detection entirely |
| `ewma.exp` | -1.5 | EWMA steps | Exponent for EWMA weighting |
| `ewma_window` | 15 | EWMA steps | Max observations on each side for EWMA |
| `adult_cutpoint` | 20 | Preprocessing | Age (years) dividing pediatric/adult |
| `use_legacy_algorithm` | FALSE | Dispatch | If TRUE, use legacy pediatric algorithm |
| `quietly` | TRUE | All | Suppress progress messages |

---

## Complete List of Steps

| Step | Name | Brief description |
|------|------|-------------------|
| (preprocessing) | Z-score calculation | WHO/CDC CSD z-scores, blending |
| (preprocessing) | Step 2b: GA correction | Fenton-based correction for potcorr subjects |
| (preprocessing) | Recentering | Subtract population medians to create tbc.sd/ctbc.sd |
| Early 13 | SDE-Identicals | Remove same-day identical values before CF detection |
| 5 | Temporary SDE | Temporarily flag same-day duplicates |
| 6 | Carried Forwards | Identify and optionally rescue carried-forward values |
| 7 | BIV | Exclude biologically implausible values (absolute + standardized) |
| 9 | Evil Twins | Exclude adjacent extreme values |
| 11 | EWMA1 (Extreme) | Exclude extreme EWMA outliers |
| 13 | Final SDE | Resolve remaining same-day duplicates |
| 15 | EWMA2 (Moderate) | Exclude moderate EWMA outliers |
| 16 | Birth HT/HC | EWMA2 variant for birth height and HC values |
| 17 | Height/HC Velocity | Exclude values exceeding velocity limits |
| 19 | Pairs and Singles | Evaluate subjects with 1–2 remaining measurements |
| 21 | Error Load | Exclude all if error ratio too high |
| 22 | Output | Assemble return columns |

Step numbers are not consecutive because they are aligned
with the parallel Stata implementation. Steps 3 (unit error),
4 (swaps), 8, 10, 12, 14, 18, and 20 either do not exist
or are handled within other steps in the child algorithm.

---

## Open Questions and Issues to Investigate

- [x] **BUG FIXED: Z-score blending in
 `calc_and_recenter_z_scores()`:** CDC-only cutoff was
 `>= 4 years`; fixed to `> 5` to match main z-score
 blending. See "Age blending — calc_and_recenter_z_scores()"
 section above.
- [x] **HC age limits:** "Exclude-Not-Cleaned" is applied at
 `agedays > 3 × 365.25`, but "Exclude-Missing" is
 applied at `agedays >= 5 × 365.25`.
 **Resolved:** "Exclude-Not-Cleaned" at > 3y fires first during
 initialization, so HC > 3y is never cleaned. The
 >= 5y Missing assignment is functionally redundant for
 HC but catches the z-score reference gap (WHO HC
 reference data only available through age 5). Not a bug;
 consider adding a comment in code for clarity.
- [x] **BUG FIXED: Batching wrapper defeats its own purpose.**
 Fixed by removing the overwrite and creating a batch-level
 reference for result assembly.
- [x] **Debug output fixed:** Step 5 `stop()` on duplicate
 Include values converted to `warning()`. Also fixed in
 Step 13 (two locations).
- [ ] **Exclusion code rename and cleanup:** All exclusion
 code names will be renamed at the end of development.
 The exclude.levels list also includes many legacy codes
 from the original pediatric algorithm that the child
 algorithm does not produce — clean up at the same time.
- [x] **Parallel processing fixed (2026-04-04).** Three
 bugs resolved (see CLAUDE.md). Requires installed package
 via `devtools::install_local()`.
- [ ] **Batch size should be parameterized.** Currently
 hard-coded at 2,000 subjects. Add a `batch_size` argument
 to `cleangrowth()`. See "Batching and Dispatch — Batch
 size" section.
- [x] **Batch-invariant operations moved outside the loop
 (2026-04-13).** `exclude.levels` definition and Tanner/WHO
 velocity reference reads moved before the loop.
 `read_anthro()` calls already optimized via `ref_tables`.

---

## Step 13: SDE Resolution (including Early Step 13)

| | |
|---|---|
| **Scope** | All parameters |
| **Exclusion codes** | `Exclude-C-{WT|HT|HC}-Identical`, `Exclude-C-{WT|HT|HC}-Extraneous` |
| **Code location** | See code |

### Overview

Same-day extraneous (SDE) resolution runs in two phases:

1. **Early Step 13** (before Steps 5/6): Removes same-day
 identical values so they don't inflate CF string
 detection
2. **Main Step 13** (after Step 11): Permanently resolves
 all remaining same-day duplicates using EWMA-based
 selection

### Early Step 13: SDE-Identicals

| | |
|---|---|
| **Prior step** | Preprocessing |
| **Next step** | Step 5 (Temporary SDEs) |

Runs on `data.df` before any other cleaning step.

**Logic:**
1. Sort by `subjid, param, agedays, id`
2. For each `(subjid, param, agedays, v)` group — same
 subject, parameter, day, AND measurement value — count
 Include rows
3. If count > 1, select one to keep via **age-dependent id
 rule** (birth: lowest id; other ages: highest id)
4. Others → `Exclude-C-{WT|HT|HC}-Identical`

Handles partial identicals: [10.5, 10.5, 11.2] excludes
one 10.5 but leaves 11.2 untouched.

### Main Step 13

| | |
|---|---|
| **Prior step** | Step 11 (EWMA1) |
| **Next step** | Step 15 (EWMA2) |

After Steps 5–11 have removed BIVs, Evil Twins, extreme
EWMA outliers, and carried forwards, Step 13 permanently
resolves all remaining same-day duplicates. It works on a
subset (`data.sde`) containing only subjects with same-day
measurements, then merges results back.

**Phase A: Setup**

1. Capture temp SDE ids for DOP median exclusion
 (`temp_sde_ids_step13`)
2. Safety checks: warn if Include duplicates exist (these
 would indicate a bug in prior steps)
3. Re-run `temporary_extraneous_infants()` with
 `exclude_from_dop_ids` — in Step 13, temp SDEs are
 excluded from DOP median calculation (unlike Step 5
 where all values contribute)
4. Pre-filter to subjects with same-day measurements

**Phase B1: SDE-Identicals**

Runs again on `data.sde` to catch identicals that emerged
after Steps 5–11 removed intervening values. Two checks:

1. **Whole-day identical**: All values on
 a day are the same → keep one, exclude rest
2. **Partial identical**: Duplicate values
 mixed with different values → mark duplicates of each
 value, keeping one per value group

Both use the age-dependent id rule. Only Include rows are
candidates (`exclude == "Include"` guard).

**Phase B2: One-Day SDE + SDE-All-Extreme**

For subject-params with data on only ONE unique day
(`one_day_sde_flag`):

1. **SDE-All-Extreme**: If the closest
 value to the day's median is still `> 2` SD away, ALL
 values on that day are marked
 `Exclude-C-{WT|HT|HC}-Extraneous`. Threshold: strict `> 2`.

2. **DOP medians**: For one-day SDEs,
 cross-parameter medians provide a secondary sort key.
 DOP mapping: WT↔HT, HC→HT. Uses only fully included
 values (`!was_temp_sde`).

3. **One-Day selection**: Sort by:
 - `absdiff_rel_to_median` ascending
 - `absdiff_dop_for_sort` ascending (Inf if NA)
 - Age-dependent `internal_id` tiebreaker
 Keep first; others → `Exclude-C-{WT|HT|HC}-Extraneous`

**Phase B3: SDE-EWMA resolution**

For subject-params with data on multiple days:

1. Compute EWMA from fully included non-temp-SDE values
 only
2. Assign EWMA values to temp SDEs by copying the
 `ewma.all` from their non-SDE same-day counterpart
3. `spa_ewma = max(ewma.all)` per `(subjid, param,
 agedays)` — used as the reference point
4. `absdewma = |tbc.sd - spa_ewma|` for each row
5. **SDE-All-Extreme**: If
 `min_absdewma > 1` for a group with 2+ eligible values,
 all → `Exclude-C-{WT|HT|HC}-Extraneous`. Threshold: strict
 `> 1`.
6. **SDE-EWMA selection**: Sort eligible
 values by `absdewma` ascending, age-dependent
 `internal_id` tiebreaker. Keep first → Include; others
 → `Exclude-C-{WT|HT|HC}-Extraneous`

**Phase B4: Merge back**

SDE results are merged back to `data.df` by `id`. Only
rows that were Include or Temp-SDE are overwritten —
permanent exclusions from prior steps (BIV, Evil Twins,
etc.) are preserved. Extra columns from SDE processing
are dropped.

### Checklist findings

1. **Debug `stop()` blocks converted to `warning()`:**
 Two duplicate-Include checks (pre and post temp SDE
 marking) were hard stops with CSV file writing. Converted
 to warnings matching the Step 5 pattern.
2. **Removed debug save points:** Commented-out saveRDS
 blocks for Step 13 input and output.
3. **Removed unnecessary `_rounded` aliases:**
 `absdiff_median_rounded` and `absdiff_dop_rounded` were
 exact copies. Replaced references with originals.
4. **DOP mapping consistent:** HC → HT in all three DOP
 median calculations.
5. **Boundaries:** SDE-All-Extreme `> 2` (one-day) and
 `> 1` (EWMA-based) — both strict. Age-dependent id
 tiebreaker uses `internal_id` (not `id`) in Step 13.
6. **`valid()` calls:** Temp SDEs included when building
 `data.sde` subset. EWMA computed on Include-only
 (excludes temp SDEs). DOP medians exclude temp SDEs
 via `!was_temp_sde`.
7. **Merge safety:** Only Include/Temp-SDE rows overwritten
 on merge back — permanent exclusions preserved.
8. **Parameter scope:** All 3 params handled. DOP has
 param-specific mapping.
9. **Factor levels:** All 5 SDE exclusion codes exist in
 `exclude.levels`.
10. **Comment:** "h. no need to calculate for
 R." — cryptic reference to a Stata step. Harmless.

---

## Step 5: Temporary SDE Resolution

| | |
|---|---|
| **Scope** | All parameters |
| **Operates on** | Valid rows (via `valid()`, includes temp SDEs) |
| **Prior step** | Early Step 13 (SDE-Identicals) |
| **Next step** | Step 6 (Carried Forwards) |
| **Exclusion code** | `Exclude-C-Temp-Same-Day` |
| **Code location** | See code |

### Overview

After identical same-day values have been removed, remaining
same-day duplicates (SDEs with **different** values) need
temporary resolution. Step 5 selects the most plausible value
for each same-day group and temporarily excludes the others.
These temporary exclusions are revisited in Step 13 (final
SDE resolution) after the dataset has been further cleaned.

The function `temporary_extraneous_infants()` is also reused
in Step 13 with different parameters (see Step 13).

### Key terms

| Term | Meaning |
|---|---|
| SP median | Subject-parameter median: median `tbc.sd` for a given subject and parameter, across all valid values |
| DOP median | Designated-other-parameter median: the SP median of the other parameter for the same subject (WT↔HT, HC→HT) |
| `absdmedian.spz` | Absolute distance from `tbc.sd` to SP median |
| `absdmedian.dopz` | Absolute distance from `tbc.sd` to DOP median |

### Logic

**Input:** A data.table with columns `id, internal_id,
subjid, param, agedays, tbc.sd, exclude`.

**Step 5a — Identify SDE days:**
For each `(subjid, param, agedays)` group, if more than one
valid row exists, all rows in the group are flagged as
`extraneous.this.day = TRUE`.

**Step 5b — Compute SP medians:**
For each `(subjid, param)`, compute `median.spz` = median
of `tbc.sd` across ALL valid values (not just non-SDE days).
This matches the Stata implementation.

**Step 5c — Compute DOP medians:**
For each parameter, look up the SP median of its designated
other parameter for the same subject:
- WEIGHTKG → uses HEIGHTCM median
- HEIGHTCM → uses WEIGHTKG median
- HEADCM → uses HEIGHTCM median

In Step 5, all valid values contribute to DOP medians.
In Step 13 (when this function is reused), temp SDEs are
excluded from DOP median calculation via the
`exclude_from_dop_ids` parameter.

**Step 5d — Compute distances:**
For each SDE-day row:
- `absdmedian.spz = |tbc.sd - median.spz|`
- `absdmedian.dopz = |tbc.sd - median.dopz|`

Fallback logic when medians are unavailable:
1. If SP median is NA but DOP median exists: use DOP median
 as the primary distance (`absdmedian.spz = |tbc.sd - median.dopz|`)
2. If both medians are NA: use `|tbc.sd|` (treat 0 as
 median)
3. If DOP median is NA: set `absdmedian.dopz = Inf` (sorts
 last)

**Step 5e — Select value to keep:**
Within each `(subjid, param, agedays)` group, sort by:
1. `absdmedian.spz` ascending (closest to SP median first)
2. `absdmedian.dopz` ascending (closest to DOP median as
 tiebreaker)
3. Age-dependent id tiebreaker:
 - At birth (`agedays == 0`): `id` ascending (keep lowest)
 - At other ages: `id` descending (keep highest)

The first value after sorting is kept; all others are marked
`extraneous = TRUE`.

**Step 5f — Return result:**
The function returns a logical vector in the original row
order (not the `keyby`-sorted order), indicating which rows
should be excluded. The caller sets these rows to
`Exclude-C-Temp-Same-Day`.

### Rationale

- **SP median as primary:** The value closest to the
 subject's own trajectory (as captured by their median
 z-score) is most likely correct.
- **DOP median as tiebreaker:** When two values are equally
 close to the SP median, the one whose z-score is also
 closer to the other parameter's median is preferred. This
 cross-parameter check adds robustness.
- **Temporary resolution:** This is intentionally
 conservative — the dataset still contains errors that
 could distort medians. Step 13 re-resolves SDEs after
 BIV, Evil Twins, and EWMA1 have removed more extreme
 values.

### Safety check (after Step 5, the code)

After Step 5, a check verifies that no `(subjid, param,
agedays)` group has more than one `Include` row. If
duplicates are found, a warning is issued. This would
indicate a bug in `temporary_extraneous_infants()` —
the check has never fired in testing.

### Also creates

- `v.orig`: A copy of the measurement values
 before any transformations. Created just before Step 5
 runs, used in CF detection (exact value comparison) and
 wholehalfimp calculation. Dropped after Step 6.
- `subj.dup`: List of subject IDs that have
 at least one temporary SDE. Used for efficiency in
 later steps that only need to process SDE subjects.

---

## Step 6: Carried Forwards

| | |
|---|---|
| **Scope** | All parameters |
| **Operates on** | Valid rows including temp SDEs, excluding single-measurement subject-params |
| **Prior step** | Step 5 (Temporary SDEs) |
| **Next step** | Step 7 (BIV) |
| **Exclusion code** | `Exclude-C-{WT|HT|HC}-CF` |
| **Rescue codes** | `Rescued`, `Rescued-Imperial`, `Rescued-Adol`, `Rescued-Adol-Imperial` (in `cf_rescued` column) |
| **Code location** | See code |
| **Controlled by** | `include.carryforward` parameter (if TRUE, skip entirely) |

### Overview

A carried-forward (CF) value is a measurement identical to
the measurement on the immediately prior day for the same
subject and parameter. CFs typically represent data entry
artifacts (copying the previous value) rather than true
measurements. Step 6 identifies CFs, organizes them into
strings (consecutive sequences), evaluates each for possible
rescue based on z-score similarity, and either excludes or
re-includes them.

The step has four phases:
1. CF detection
2. Temp SDE re-evaluation (after CFs removed)
3. CF string construction and rescue evaluation
4. SDE-Identical row restoration

### Phase 1: CF detection

**Pre-filtering:** Only subject-params with at least 2
valid measurements AND at least one duplicate value
(`uniqueN(v.orig) < .N`) are processed. This skip
optimization is logged.

**CF logic:**
1. Count values per `(subjid, param, agedays)` — days with
 multiple values (SDEs) are noted
2. For each unique ageday per subject-param, find the prior
 ageday (via `shift(agedays, type = "lag")`)
3. Get the measurement from the prior ageday — but ONLY if
 that day had exactly 1 value. If the prior day had
 multiple values (SDEs), `single_val` is set to NA,
 preventing CF matching.
4. A row is CF if `prior_single_val` is not NA AND
 `v.orig == prior_single_val` (exact equality, no
 tolerance)

**Key constraint:** CFs are only detected against single-
value prior days. This prevents false CF detection when a
prior day has SDEs with different values, one of which
happens to match.

### Phase 2: Temp SDE re-evaluation

After CFs are marked, the dataset's remaining valid
measurements may have changed enough to affect which SDE
value is most plausible. If any temp SDEs exist:
1. Reset all `Exclude-C-Temp-Same-Day` back
 to `Include`
2. Re-run `temporary_extraneous_infants()` (the same
 function from Step 5) on the updated dataset

This ensures SDE selection accounts for CF removals.

### Phase 3: CF string construction and rescue

**3a. wholehalfimp calculation:**

Determines whether each measurement is in whole or half
imperial units (row-level flag, not subject-level):
- WEIGHTKG: `v.orig * 2.20462262 mod 1 < 0.01` (whole
 pounds)
- HEIGHTCM: `v.orig / 2.54 mod 1 < 0.01` OR
 `v.orig / 2.54 mod 0.5 < 0.01` (whole or half inches)
- HEADCM: same as HEIGHTCM

Note: For HEIGHTCM/HEADCM, the `mod 1` check is redundant
since `mod 0.5` is a superset. Not a bug.

**3b. SDE-Identical temporary removal:**

SDE-Identical rows are temporarily removed from `data.df`
before string construction. They would break the
positional CF string logic (rle, originator detection).
They are re-added after rescue processing.

**3c. Positional string detection:**

1. `ageday_has_include`: For each
 `(subjid, param, agedays)`, TRUE if any row has
 `exclude == "Include"`. CFs on days that also have
 Includes are **not eligible for rescue**.

2. `cf_binary`: TRUE for CF rows.

3. `originator`: An Include row whose NEXT
 row (in subject-param order) is a CF. Any Include can
 be an originator regardless of `ageday_has_include`.

4. `cf_string_num`: Sequential string
 number per subject-param. Assigned to originators first,
 then propagated forward through consecutive CFs (via a
 loop, the relevant code section). Propagation only reaches CFs
 where `ageday_has_include` is FALSE — CFs on days with
 Includes don't get string numbers and are not eligible
 for rescue.

5. `originator_z`: The originator's
 `sd.orig_uncorr` z-score. Propagated alongside
 `cf_string_num`.

6. `seq_win`: Position within the
 string — 0 for the originator, 1/2/3... for CFs.

7. `absdiff`: `|sd.orig_uncorr - originator_z|`
 — how far the CF's z-score is from the originator's.

**3d. Rescue evaluation:**

Rescue is applied per `(subjid, param, cs)` group. Only
CFs with valid `seq_win` AND `ageday_has_include == FALSE`
are evaluated. Rescue criteria depend on string length:

**Single CF (`max(seq_win) == 1`):**
- `absdiff < 0.05` → rescue code
 `Rescued`
- `absdiff >= 0.05 & absdiff < 0.1 & wholehalfimp` →
 rescue code `Rescued-Imperial`

**Multi-CF string (`max(seq_win) > 1`):**
Only adolescents are eligible:
- Female (`sex == 1`) age `>= 16` years, OR
 Male (`sex == 0`) age `>= 17` years
- Same `absdiff` thresholds as single CF, with rescue
 codes `Rescued-Adol` and
 `Rescued-Adol-Imperial`

**3e. Re-inclusion:**

Rescued CFs have their rescue reason stored in
`cf_rescued`, then `exclude` is set back to `"Include"`.
This means rescued CFs participate in all downstream
algorithm steps.

### Phase 4: Cleanup

1. SDE-Identical rows are added back to `data.df` with
 `rbind(..., fill = TRUE)` and re-sorted
2. `v.orig` and `wholehalfimp` are dropped (no longer
 needed)

### Checklist findings

1. **Stale comments:** Large blocks of commented-out old
 dplyr/map_lgl CF logic remain. Should be removed. A
 comment says "non NNTE values" but nnte was removed —
 comment is stale.
2. **Unnecessary variable:** `cf_string_length`
 is computed and immediately deleted. The
 rescue logic uses `max(seq_win)` instead. Should be
 removed.
3. **Uncleaned variables:** `seq_win`, `cs`, `absdiff`,
 and `ageday_has_include` persist in `data.df` after
 Step 6 but are never used again. Should be added to
 `cols_to_drop_6`.
4. **Redundant wholehalfimp check:** For HEIGHTCM/HEADCM,
 the `mod 1` check is a subset of the `mod 0.5` check.
 Not a bug but could be simplified.
5. **Boundary: `absdiff < 0.05`** — strict less-than.
 A CF with exactly `absdiff == 0.05` is NOT rescued by
 the first criterion but IS eligible for the second
 (`>= 0.05 & < 0.1 & wholehalfimp`).
6. **Boundary: teen age `>= 16` / `>= 17`** — inclusive.
 Exactly 16-year-old females and 17-year-old males are
 eligible for multi-CF rescue.
7. **Edge case:** If `include.carryforward = TRUE`, the
 entire step is skipped, but `sde_identical_rows` is
 still initialized and the SDE-Identical
 restoration still runs (it's outside the
 `if` block). This is correct — SDE-Identicals were
 removed before Step 6 regardless.
8. **`valid()` call:** Uses `include.temporary.extraneous
 = TRUE`, which is correct — temp SDEs should participate
 in CF detection (they might be CFs themselves).
9. **Parameter scope:** All 3 parameters handled.
 wholehalfimp has HC-specific logic.
10. **Factor levels:** All 4 rescue codes exist in
 `exclude.levels`. `Exclude-C-{WT|HT|HC}-CF` also exists.
11. **Sort order:** `setorder(cf_subset, subjid, param,
 agedays, id)` ensures deterministic CF
 detection order.

---

## Step 7: BIV (Biologically Implausible Values)

| | |
|---|---|
| **Scope** | All parameters |
| **Operates on** | Valid rows including temp SDEs |
| **Prior step** | Step 6 (Carried Forwards) |
| **Next step** | Step 9 (Evil Twins) |
| **Exclusion codes** | `Exclude-C-{WT|HT|HC}-BIV` |
| **Code location** | See code |

### Overview

BIV identifies measurements that are biologically impossible
or extremely unlikely, using two approaches applied in
sequence:

1. **Absolute BIV** — fixed measurement limits (e.g., weight
 cannot be negative, height cannot exceed 244 cm)
2. **Standardized BIV** — z-score cutoffs using unrecentered
 CSD z-scores (`sd.orig_uncorr`)

After both, temp SDEs are re-evaluated (Step 5 rerun).

### Absolute BIV thresholds

**Weight (WEIGHTKG):**

| Condition | Threshold | Rationale |
|---|---|---|
| First year (`agedays <= 365`) | `v < 0.2` kg | Published minimum viable birth weight; allows preterm |
| After first year (`agedays > 365`) | `v < 1` kg | Minimum plausible outpatient weight |
| Birth | `v > 10.5` kg | Maximum plausible birth weight |
| Age < 2 years | `v > 35` kg | Highest plausible weight under 2y |
| Any age | `v > 600` kg | Published maximum human weight |

***BUG FIXED (2026-03-20):*** The code previously applied
`v < 0.2` only at `agedays == 0`, leaving ages 1–364 days
with no minimum weight check. Fixed to `agedays <= 365`
to match the intended behavior documented in CLAUDE.md.

**Height (HEIGHTCM):**

| Condition | Threshold |
|---|---|
| Any age | `v < 18` cm |
| Any age | `v > 244` cm |
| Birth | `v > 65` cm |

**Head circumference (HEADCM):**

| Condition | Threshold |
|---|---|
| Any age | `v < 13` cm |
| Any age | `v > 75` cm |
| Birth | `v > 50` cm |

### Standardized BIV thresholds

Uses **unrecentered** z-scores (`sd.orig_uncorr`). Rows
already marked `Exclude-C-{WT|HT|HC}-BIV` are skipped (the
`exclude != abs_biv` guard) to preserve the more specific
code.

**Weight:**

| Condition | Threshold |
|---|---|
| Age < 1 year | `sd.orig_uncorr < -25` |
| Age >= 1 year | `sd.orig_uncorr < -15` |
| Any age | `sd.orig_uncorr > 22` |

**Height:**

| Condition | Threshold |
|---|---|
| Age < 1 year | `sd.orig_uncorr < -25` |
| Age >= 1 year | `sd.orig_uncorr < -15` |
| Any age | `sd.orig_uncorr > 8` |

**Head circumference:**

| Condition | Threshold |
|---|---|
| Any age | `sd.orig_uncorr < -15` |
| Any age | `sd.orig_uncorr > 15` |

### Temp SDE re-evaluation

After BIV exclusions, temp SDEs are reset to Include and
`temporary_extraneous_infants()` is rerun. This is the same
pattern as after CF detection — removing extreme values
changes which SDE value is most plausible.

### Variables created and dropped

- `ageyears`: Created for age threshold checks,
 dropped after Step 7 completes.

### Checklist findings

1. **BUG FIXED: Minimum weight gap for agedays 1–364.**
 Code used `agedays == 0` for the `v < 0.2` threshold;
 fixed to `agedays <= 365`. Also changed `v < 1` from
 `agedays >= 365` to `agedays > 365` to avoid overlap
 at exactly day 365.
2. **Cleaned up:** Removed commented-out duplicate code
 and stale nnte comment.
4. **Boundaries — all strict inequalities:** `< 0.2`,
 `< 1`, `> 10.5`, etc. A weight of exactly 0.2 kg at
 birth is NOT excluded. A z-score of exactly -25 is NOT
 excluded. These are strict `<` and `>`, not `<=`/`>=`.
5. **Boundaries — age thresholds:** `ageyears < 1` and
 `ageyears >= 1` for Standardized-BIV use strict `<`
 and `>=`. A child at exactly 365.25 days (1.0 years)
 gets the `>= 1` thresholds.
6. **`valid()` call:** Correctly includes temp SDEs.
 `valid_set` is computed once and not refreshed between
 Absolute and Standardized BIV, but the `!= abs_biv`
 guard prevents double-exclusion.
7. **v == 0 handler removed (2026-04-13):** Was dead code —
 `v` is already NaN for zero measurements after
 preprocessing. Removed.
8. **Parameter scope:** All 3 params handled with
 param-specific thresholds.
9. **Factor levels:** `Exclude-C-{WT|HT|HC}-BIV` exists in
 `exclude.levels` (absolute and standardized BIV collapsed
 into one code).
10. **Efficiency:** `valid_set` could be refreshed after
 Absolute-BIV to avoid checking already-excluded rows
 in Standardized-BIV. Minor — only matters for very
 large datasets.

---

## Step 9: Evil Twins

| | |
|---|---|
| **Scope** | All parameters |
| **Operates on** | Valid rows, excluding temp SDEs, requiring 3+ total rows per subject-param |
| **Prior step** | Step 7 (BIV) |
| **Next step** | Step 11 (EWMA1 — Extreme EWMA) |
| **Exclusion code** | `Exclude-C-{WT|HT|HC}-Evil-Twins` |
| **Code location** | See code |

### Overview

Evil Twins identifies adjacent measurements that are both
extreme — cases where two or more consecutive values are
far from the subject's trajectory. The original pediatric
algorithm often missed these because each value's EWMA was
distorted by its extreme neighbor(s). This step catches
them before EWMA processing.

The step is **iterative**: it finds the single worst OTL
value, excludes it, recalculates OTL on the remaining
values, and repeats until no OTL values remain.

### OTL calculation (`calc_otl_evil_twins`, the relevant code section)

A measurement is "over the limit" (OTL) if it differs from
at least one adjacent measurement (in time order within
the same subject-param) by more than 5 in **both** `tbc.sd`
and `ctbc.sd`:

```
OTL = (|tbc_next - tbc_current| > 5 AND
 |ctbc_next - ctbc_current| > 5)
 OR (|tbc_prev - tbc_current| > 5 AND
 |ctbc_prev - ctbc_current| > 5)
```

**Key details:**
- Threshold is **strict `> 5`** — a difference of exactly
 5.0 is NOT OTL
- Both `tbc.sd` AND `ctbc.sd` must exceed 5 for the
 **same** pair (not mixed across neighbors). This dual
 requirement prevents false positives when correction
 moves the z-score close to the boundary.
- Adjacent means the next/prior row within the same
 `(subjid, param)` in agedays order
- Cross-subject/param boundaries are handled by padding
 with `Inf` (always exceeds threshold, but `same_sp_next/
 prev` flags prevent false matches)
- Edge case: < 2 rows returns all `otl = FALSE`

### Step 9 main logic

**9a. Setup:**
- Sort by `(subjid, param, agedays, id)` — MUST happen
 before `valid_set` computation (comment explains why)
- `valid_set`: excludes temp SDEs
 (`include.temporary.extraneous = FALSE`), requires
 `sp_count_9 > 2` (total row count, loose pre-filter)
- Initial `calc_otl_evil_twins` on all valid rows to check
 if any OTL exists (cheap early exit)

**9b. Per-group processing:**
Only subject-params with at least one OTL value are
processed. For each:

1. Extract valid rows for this `(subjid, param)`
2. Skip if < 2 valid rows
3. **While loop** — iterate until no OTL remains:
 a. Compute median of `tbc.sd` for Include rows
 b. Compute `med_diff = |tbc.sd - median|` for each
 Include row
 c. Find the worst OTL value by sort priority:
 - Highest `med_diff` (furthest from median)
 - Highest `|tbc.sd|` (most extreme overall)
 - Lowest `id` (deterministic tiebreaker)
 d. Exclude that one value (`Exclude-C-{WT|HT|HC}-Evil-Twins`)
 e. Recalculate OTL on remaining Include values
 f. Break if < 2 Include values remain

**9c. Apply exclusions:**
All exclusions are collected as `line` values and applied
back to `data.df` in bulk.

**9d. Cleanup and temp SDE rerun:**
- `sp_count_9` dropped
- Temp SDEs reset to Include and re-evaluated

### Checklist findings

1. **Unnecessary `_rounded` aliases in
 `calc_otl_evil_twins`:** The code creates
 `tbc_next_diff_rounded` etc. as exact copies of the
 unrounded values. The comment says "Handles
 floating-point noise" but no rounding is applied.
 These are unnecessary variables — should use the
 originals directly.
2. **Stale comment:** the code mentions "non NNTE values."
 NNTE was removed.
3. **`sp_count_9` counts all rows, not just valid:** This
 is a loose pre-filter. A subject-param with 3 total
 rows (2 excluded + 1 valid) passes the count check
 but is caught by `length(grp_idx) < 2L`.
 Not a bug, just slightly wasteful.
4. **Boundaries:** OTL threshold is strict `> 5`. Tiebreaker
 uses lowest `id` (not age-dependent like SDE steps).
5. **`valid()` call:** Correctly excludes temp SDEs
 (`include.temporary.extraneous = FALSE`). Evil Twins
 operates on the "clean" view without SDE candidates.
6. **Sort order:** Explicit sort before
 `valid_set` computation. Comment explains this is
 critical — sorting after would misalign the boolean
 vector.
7. **Both `tbc.sd` and `ctbc.sd` used:** The dual
 threshold ensures corrected subjects (potcorr) aren't
 falsely flagged when correction moves one z-score
 close to the neighbor but not the other.
8. **Parameter scope:** All 3 params handled uniformly.
9. **Factor levels:** `Exclude-C-{WT|HT|HC}-Evil-Twins` exists in
 `exclude.levels`.
10. **Efficiency:** Good — pre-filter checks for any OTL
 globally before entering per-group loop. Only
 subject-params with OTL are processed. The while loop
 typically runs 1–3 iterations per group.

---

## Step 11: EWMA1 (Extreme EWMA)

| | |
|---|---|
| **Scope** | All parameters |
| **Operates on** | Include values only (no temp SDEs), requiring 3+ per subject-param |
| **Prior step** | Step 9 (Evil Twins) |
| **Next step** | Step 13 (Final SDE Resolution) |
| **Exclusion code** | `Exclude-C-{WT|HT|HC}-Traj-Extreme` |
| **Code location** | See code |

### Overview

EWMA1 identifies extreme outliers using Exponentially
Weighted Moving Averages. For each measurement, the EWMA
predicts what the value "should" be based on the subject's
other measurements for the same parameter. Values that
deviate extremely from their EWMA prediction are excluded.

The step is **iterative at two levels**:
1. Per subject-param: exclude at most ONE worst value
2. Globally: repeat across all subject-params until no
 new exclusions occur

This prevents one extreme value from distorting the EWMA
of nearby values, causing a cascade of false exclusions.

### Pre-filtering

Two pre-filters narrow the work:
1. **Count filter:** Subject-params with <= 2 Include
 values are skipped (EWMA needs 3+ points)
2. **Extreme filter:** Subject-params where no Include
 value has `|tbc.sd| > 3.5` are skipped (impossible to
 meet exclusion criteria)

### EWMA calculation

For each subject-param group with `> 2` Include values:

1. **Exponent calculation**: Age-gap-
 dependent exponent using linear interpolation:
 - Gap <= 1 year: exponent = -1.5
 - Gap >= 3 years: exponent = -3.5
 - Between: linear interpolation
 (`-1.5 - (ageyears - 1)`)
 - `ageyears` is based on the maximum of the gap
 before and after each measurement

2. **EWMA computation**: Calls
 `ewma()` function with windowed matrix operations.
 Three variants computed:
 - `ewma.all`: excluding only the current observation
 - `ewma.before`: excluding current + prior
 - `ewma.after`: excluding current + next
 - Same three computed for `ctbc.sd` (corrected),
 prefixed with `c.`. If `ctbc.sd == tbc.sd` for all
 rows (vast majority — non-potcorr subjects), the
 ctbc columns are copied rather than recomputed.

3. **Deviation (dewma)**:
 - `dewma.all = tbc.sd - ewma.all`
 - `dewma.before = tbc.sd - ewma.before`
 - `dewma.after = tbc.sd - ewma.after`
 - `c.dewma.all = ctbc.sd - c.ewma.all`

### Exclusion criteria

A value is a potential exclusion if ALL of:

**Positive outlier:**
- `dewma.all > 3.5`
- `dewma.before > 3`
- `dewma.after > 3`
- `tbc.sd > 3.5`
- `c.dewma.all > 3.5` OR `ctbc.sd` is NA

**OR negative outlier:**
- `dewma.all < -3.5`
- `dewma.before < -3`
- `dewma.after < -3`
- `tbc.sd < -3.5`
- `c.dewma.all < -3.5` OR `ctbc.sd` is NA

The `c.dewma.all` check uses `ctbc.sd` (not `c.dewma.all`)
for the NA test. When `ctbc.sd` is NA (non-potcorr subjects
or missing corrected z-score), the corrected dewma check
passes automatically. When `ctbc.sd` is available,
`c.dewma.all = ctbc.sd - c.ewma.all` must also exceed the
threshold.

**Pre-filter:** Before per-group processing, subject-params
with `max(|tbc.sd|) <= 3.5` are skipped entirely (if the
most extreme z-score doesn't exceed 3.5, no value can meet
the `tbc.sd > 3.5` criterion).

**Additional constraints:**
- First and last Include values are never excluded
- At most ONE value excluded per subject-param per
 iteration

### Worst-value selection

When multiple values meet criteria, the worst is selected
by: `order(pot_excl, abs(tbc.sd + dewma.all), -id,
decreasing = TRUE)`. The sort key
`abs(tbc.sd + dewma.all)` combines the z-score magnitude
with its deviation — values that are both extreme and far
from EWMA prediction score highest. Lowest `id` breaks
ties.

### Global iteration loop

After each per-group pass:
1. Identify subject-params with NEW exclusions (comparing
 before/after)
2. Recalculate temp SDEs only for subjects that had both
 new exclusions AND existing SDEs (targeted, not global)
3. Next iteration processes only subject-params with new
 exclusions
4. Loop terminates when no new exclusions occur

After the loop, a final global temp SDE recalculation runs.

### Debug columns preserved

Iteration 1 EWMA values are saved to permanent columns
(e.g., `ewma1_it1.ewma_all`, `ewma1_it1.dewma_all`) for
diagnostic output. These persist through to the final
result.

### Checklist findings

1. **Removed `nnte` from `.SDcols`:** `nnte` was passed
 into the closure but never used (filter was removed).
 Fixed.
2. **Removed stale comments:** "Removed nnte filter" (3
 instances), "Round to 3 decimals" (no rounding done).
3. **Removed debug exit block:** Commented-out
 `saveRDS`/`return` after Step 11.
4. **Boundaries — all strict:** `> 3.5`, `> 3`, `< -3.5`,
 `< -3` for dewma; `> 3.5` / `< -3.5` for `tbc.sd`
 and `c.dewma.all`.
5. **`valid()` call:** Correctly excludes temp SDEs.
 Only Include values participate in EWMA calculation
 and exclusion candidates.
6. **ctbc optimization:** Skips ctbc EWMA computation
 when `ctbc.sd == tbc.sd` (copies instead). This is
 correct for non-potcorr subjects.
7. **Parameter scope:** All 3 params handled uniformly.
8. **Factor levels:** `Exclude-C-{WT|HT|HC}-Traj-Extreme` exists in
 `exclude.levels`.
9. **Efficiency:** Good — two levels of pre-filtering
 (count + extreme), targeted SDE recalc, and
 subject-params drop out of the loop as they stop
 producing exclusions.
10. **`ewma_window` parameter:** Passed through to
 `ewma()` — controls max observations on each side
 (default 15).

---

## Steps 15/16: EWMA2 (Moderate EWMA)

| | |
|---|---|
| **Scope** | Step 15: all params (birth HT/HC excluded); Step 16: birth HT/HC only |
| **Operates on** | Include values only (no temp SDEs), 3+ per subject-param |
| **Prior step** | Step 13 (Final SDE Resolution) |
| **Next step** | Step 17 (Height/HC Velocity) |
| **Exclusion codes** | `Exclude-C-{WT|HT|HC}-Traj` (all moderate EWMA sub-rules collapsed into one code) |
| **Code location** | See code |

### Overview

EWMA2 uses the same EWMA mechanism as EWMA1 but with
**lower thresholds** and **position-specific rules**. Where
EWMA1 catches extreme outliers (thresholds ~3.5), EWMA2
catches moderate outliers (thresholds ~1–4 depending on
position and context).

Step 15 handles all measurements except birth HT/HC. Step
16 handles birth HT/HC separately because birth values
have different clinical expectations and need higher
thresholds.

### Pre-calculation (once, before iteration loop)

**p_plus / p_minus:** For each valid
measurement, compute a ±5% (weight) or ±1 cm (height/HC)
perturbation:
- WEIGHTKG: `p_plus = 1.05*v`, `p_minus = 0.95*v`
- HEIGHTCM: `p_plus = v+1`, `p_minus = v-1`
- HEADCM: `p_plus = v+1`, `p_minus = v-1`

These are converted to z-scores (`tbc.p_plus`, `tbc.p_minus`)
using `calc_and_recenter_z_scores()`. The ±5% rule tests
whether the measurement would still be an outlier if
perturbed slightly — a robustness check that reduces false
positives for values near exclusion boundaries.

**first_meas:** Marks the first
non-birth Include measurement per subject-param. Used for
the "first" exclusion rules. HT/HC exclude birth from
Step 15 entirely (birth handled in Step 16), so for HT/HC
`first_meas` is the first non-birth value.

### Additional criteria (addcrit)

All EWMA2 rules require `addcrit` — a supplementary check
that the measurement differs from BOTH neighbors even
after ±5%/±1cm perturbation. This prevents excluding
values that are close to a neighbor — even if EWMA says
the value is outlying, if it's within 1 SD of a neighbor
(even after perturbation), it's likely real.

**Neighbor difference variables** (computed within each
subject-param group, in agedays order):
- `tbc_diff_next = tbc.sd - tbc.sd[next]`
- `tbc_diff_prior = tbc.sd - tbc.sd[prior]`
- `tbc_diff_plus_next = tbc.p_plus - tbc.sd[next]`
- `tbc_diff_plus_prior = tbc.p_plus - tbc.sd[prior]`
- `tbc_diff_minus_next = tbc.p_minus - tbc.sd[next]`
- `tbc_diff_minus_prior = tbc.p_minus - tbc.sd[prior]`

Where `tbc.p_plus` and `tbc.p_minus` are the recentered
z-scores of `p_plus` and `p_minus` respectively, computed
via `calc_and_recenter_z_scores()` using the same WHO/CDC
blending and recentering as the original z-scores.

**addcrithigh:**
```
dewma.before > 1 & dewma.after > 1 &
 ((tbc_diff_next > 1 & tbc_diff_plus_next > 1 &
 tbc_diff_minus_next > 1) | is.na(tbc_diff_next)) &
 ((tbc_diff_prior > 1 & tbc_diff_plus_prior > 1 &
 tbc_diff_minus_prior > 1) | is.na(tbc_diff_prior))
```

**addcritlow:** Mirror with `< -1` for all comparisons.

NA handling: endpoints have NA for their missing neighbor
(first row has NA prior, last row has NA next), so
`addcrit` passes on that side automatically.

### Step 15 exclusion rules

All rules require:
1. The appropriate `addcrit` (addcrithigh or addcritlow)
2. `dewma.all` threshold (uncorrected)
3. `c.dewma.all` threshold (corrected) — when `c.dewma.all`
 is available, it must also exceed the threshold. If
 `c.dewma.all` is NA (non-potcorr subjects), the check
 passes automatically:
 `(c.dewma.all > threshold | is.na(c.dewma.all))`

| Rule | Position | Gap | dewma.all | c.dewma.all | Extra condition |
|---|---|---|---|---|---|
| middle | Not first/last | — | `> 1` / `< -1` | `> 1` / `< -1` | — |
| birth-WT | `agedays == 0` | next `< 365.25` | `> 3` / `< -3` | `> 3` / `< -3` | — |
| birth-WT-ext | `agedays == 0` | next `>= 365.25` | `> 4` / `< -4` | `> 4` / `< -4` | — |
| first | `first_meas` | next `< 365.25` | `> 2` / `< -2` | `> 2` / `< -2` | — |
| first-ext | `first_meas` | next `>= 365.25` | `> 3` / `< -3` | `> 3` / `< -3` | — |
| last | Last row | prev `< 2*365.25` | `> 2` / `< -2` | `> 2` / `< -2` | `abs(tbc_prev) < 2` |
| last-high | Last row | prev `< 2*365.25` | `> abs(tbc_prev)` / `< -abs(tbc_prev)` | **`> 3` / `< -3`** | `abs(tbc_prev) >= 2` |
| last-ext | Last row | prev `>= 2*365.25` | `> 3` / `< -3` | `> 3` / `< -3` | `abs(tbc_prev) < 2`, DOP diff |
| last-ext-high | Last row | prev `>= 2*365.25` | `> 1+abs(tbc_prev)` / `< -1-abs(tbc_prev)` | **`> 3` / `< -3`** | `abs(tbc_prev) >= 2`, DOP diff |

**Note:** For last-high and last-ext-high, `dewma.all`
scales with `abs(tbc_prev)` but `c.dewma.all` is hardcoded
at 3. This means the corrected check is fixed regardless
of how extreme the prior value is.

**Gap definitions:**
- birth-WT, first: gap = `agedays[next] - agedays` (to
 next measurement)
- last: gap = `agedays[last] - agedays[last-1]` (from
 prior measurement)

**first_meas definition:**
- WEIGHTKG: first Include row per subject-param where
 `agedays > 0` (birth WT is handled by birth-WT rules)
- HEIGHTCM/HEADCM: first Include row among non-birth rows
 per subject-param (birth HT/HC excluded from Step 15
 entirely — handled by Step 16)

**DOP diff (last-ext and last-ext-high):** Cross-parameter
check. For positive outliers: `tbc.sd - tbc_dop > 4` (or
`tbc_dop` is NA). For negative: `tbc.sd - tbc_dop < -4`
(or `tbc_dop` is NA). Where `tbc_dop` is the same-day DOP
value if available, otherwise the median `tbc.sd` of all
DOP values. If no DOP data exists, `tbc_dop` is NA and the
check passes.

**Key boundaries:**
- Gap thresholds: `< 365.25` vs `>= 365.25` (strict `<`);
 `< 2*365.25` vs `>= 2*365.25` for last
- All dewma thresholds are strict `>`/`<`
- `abs(tbc_prev) < 2` vs `>= 2` distinguishes "normal
 prior" from "extreme prior" for last-value rules

### Step 16: Birth HT/HC

Separate step for birth height and HC values, using the
same EWMA infrastructure but with birth-specific rules
only. Uses the same addcrit definitions and c.dewma.all
pattern as Step 15.

| Rule | Gap to next | dewma.all | c.dewma.all |
|---|---|---|---|
| birth-HT-HC | `< 365.25` | `> 3` / `< -3` | `> 3` / `< -3` |
| birth-HT-HC-ext | `>= 365.25` | `> 4` / `< -4` | `> 4` / `< -4` |

Only processes HT/HC subject-params with a birth
measurement (agedays == 0) and 3+ Include values.
Structurally identical to Step 15 (iterative, cached EWMA,
one exclusion per iteration).

### Worst-value selection

Both steps use `abs(tbc.sd + dewma.all)` as the sort key
(same as EWMA1), with lowest `id` as tiebreaker.

### Global iteration

Same pattern as EWMA1: subject-params drop out when they
stop producing new exclusions. EWMA caches (`ewma2_caches`,
`ewma2b_caches`) persist across iterations for incremental
updates. DOP snapshot refreshed each iteration.

### Variables created and dropped

- `p_plus`, `p_minus`, `tbc.p_plus`, `tbc.p_minus`,
 `first_meas`: created before Step 15, dropped after
 Step 16
- `sp_key`: created for Step 15/16, dropped after Step 16

### Checklist findings

1. **Removed debug save points:** Three commented-out
 blocks (Step 15 input, output; Step 16 output).
2. **Boundaries verified:** All dewma thresholds are
 strict `>` / `<`. Gap thresholds use strict `<` vs
 `>=` consistently.
3. **`valid()` calls:** Step 15 excludes temp SDEs AND
 birth HT/HC. Step 16 excludes temp SDEs, includes
 only birth HT/HC subjects.
4. **p_plus/p_minus values:** WT uses multiplicative
 (±5%), HT/HC uses additive (±1 cm). Consistent with
 clinical interpretation (5% of weight is proportional;
 1 cm is the measurement precision for length).
5. **first_meas:** Correctly recalculated each iteration
 (exclusions change which row is "first").
6. **DOP lookup:** Uses pre-computed keyed snapshot
 (`dop_snap`) refreshed each iteration. O(log n) lookup
 via `get_dop()`.
7. **Parameter scope:** Step 15 handles all 3 params
 (birth HT/HC excluded). Step 16 handles HT/HC birth
 only (WT birth is in Step 15).
8. **Factor levels:** All 11 EWMA2 exclusion codes exist
 in `exclude.levels`.
9. **Efficiency:** Pre-filter by tbc.sd range > 1 (same
 logic as EWMA1's > 3.5 filter but with lower bound).
 Incremental EWMA caching via `ewma_cache_init/update`.
10. **No temp SDE rerun after Step 15/16.** Unlike Steps
 7/9/11, there is no temp SDE recalculation after EWMA2.
 This is correct — temp SDEs were permanently resolved
 in Step 13.

---

## Step 17: Height/HC Velocity

| | |
|---|---|
| **Scope** | HEIGHTCM and HEADCM only (not WEIGHTKG) |
| **Operates on** | Valid rows excluding temp SDEs, 2+ per subject-param |
| **Prior step** | Steps 15/16 (EWMA2) |
| **Next step** | Step 19 (Pairs and Singles) |
| **Exclusion codes** | `Exclude-C-{HT|HC}-Abs-Diff` |
| **Code location** | See code |

### Overview

Step 17 checks whether the raw measurement change between
consecutive observations is physiologically plausible.
Unlike EWMA steps that work on z-scores, this step works
on raw cm values. Height/HC should not decrease more than
allowed (mindiff) or increase more than allowed (maxdiff).
The allowed change depends on:
- Age (via Tanner velocity tables for HT, WHO velocity for
 HT under 9 months and HC under 24 months)
- Time gap between measurements
- Parameter (HC has tighter tolerances than HT)

Weight is excluded because weight can legitimately decrease.

### Velocity references and threshold computation

**Height** uses two reference systems applied in priority
order. Each pair of adjacent measurements gets a `mindiff`
(minimum allowed change, typically negative = allowed
decrease) and `maxdiff` (maximum allowed increase):

**Tanner (ages 2.5+ years):**
- `tanner.months = 6 + 12 * round(midpoint_age_years)`
 where midpoint is the mean of the two agedays values.
 Only used when midpoint age ≥ 30 months.
- Merge with Tanner reference table → `min.ht.vel` and
 `max.ht.vel` (annual velocity values by sex and age)
- **max.ht.vel floors** (gap-dependent):
 - Always: floor at 2.54 cm (1 inch)
 - If gap > 2 months: floor at 5.08 cm (2 inches)
 - If gap > 6 months: floor at 10.16 cm (4 inches)
 - If gap > 1 year: floor at 20.32 cm (8 inches)
- **mindiff** (gap < 1 year):
 `0.5 * min.ht.vel * (gap_years)^2 - 3`
- **mindiff** (gap > 1 year):
 `0.5 * min.ht.vel - 3`
- **maxdiff** (gap < 1 year):
 `2 * max.ht.vel * (gap_years)^0.33 + 5.5`
- **maxdiff** (gap > 1 year):
 `2 * max.ht.vel * (gap_years)^1.5 + 5.5`
- When gap equals exactly 1 year, the `< 365.25` and
 `> 365.25` branches both miss — mindiff/maxdiff use
 the `< 365.25` formula from the first condition that
 applies.

**WHO (ages 0–24 months, height):**
- `whoagegrp.ht = round(agedays / 30.4375)`, capped at 24
- WHO increment age groups based on gap between adjacent
 measurements:
 - gap < 20 days → 1 month
 - 20 ≤ gap < 46 → 1 month
 - 46 ≤ gap < 76 → 2 months
 - 76 ≤ gap < 107 → 3 months
 - 107 ≤ gap < 153 → 4 months
 - 153 ≤ gap < 199 → 6 months
 - gap ≥ 200 → 6 months
- Merge with WHO velocity reference → `who_mindiff_ht`,
 `who_maxdiff_ht` at appropriate increment intervals
- **Gap scaling:** If actual gap < reference interval
 (e.g., 45 days for a 2-month reference), scale mindiff
 proportionally. If gap > reference interval, scale
 maxdiff proportionally.
- **For gaps < 9 months** (WHO preferred):
 - Transform: `mindiff = who_mindiff * 0.5 - 3`
 - Transform: `maxdiff = who_maxdiff * 2 + 3`
- **For gaps ≥ 9 months** (Tanner preferred, WHO fallback):
 - Only use WHO if Tanner data unavailable (no
 `min.ht.vel`)
 - WHO fallback transform:
 `mindiff = who_mindiff * 0.5 - 3`,
 `maxdiff = who_maxdiff * 2 + 3`
- **Default** when no reference data available:
 `mindiff = -3` (allow up to 3 cm decrease)
- **Birth adjustment (agedays == 0):**
 `mindiff -= 1.5`, `maxdiff += 1.5`

**HC** uses WHO only (no Tanner):
- WHO increment age groups (same gap-based mapping as HT
 but using HC-specific reference tables):
 - gap < 20 days → 1 month
 - 20 ≤ gap < 46 → 1 month
 - 46 ≤ gap < 76 → 2 months
 - 76 ≤ gap < 107 → 3 months
 - 107 ≤ gap < 153 → 4 months
 - 153 ≤ gap < 199 → 6 months
 - gap ≥ 200 → 6 months
- **Gap scaling:** Same proportional scaling as HT
- **Tolerance transform:**
 `mindiff = who_mindiff * 0.5 - 1.5`
 `maxdiff = who_maxdiff * 2 + 1.5`
 (tighter than HT: ±1.5 cm vs ±3 cm)
- **Default** when no WHO data: `mindiff = -1.5`
- **Birth adjustment (agedays == 0):**
 `mindiff -= 0.5`, `maxdiff += 0.5`
 (tighter than HT: ±0.5 cm vs ±1.5 cm)

### Violation detection and resolution

Step 17 uses a pre-filter that computes all thresholds
in one vectorized pass, then only enters the per-group
while loop for subject-params with at least one violation.

**Violation condition:** A pair violates if:
`diff_prev < mindiff_prior` (decreased too much from prior)
OR `diff_next < mindiff` (next decreased too much from this)
OR `diff_prev > maxdiff_prior` (increased too much from prior)
OR `diff_next > maxdiff` (next increased too much from this)

Where `diff_prev = v - v[prior]`, `diff_next = v[next] - v`
(raw cm differences, not z-scores).

**Resolution within each group:**

For **3+ measurements**, EWMA-based tiebreaker:
1. Compute EWMA on raw `tbc.sd` values
2. For each violating pair, flag which member has a worse
 EWMA relationship:
 - `bef.g.aftm1`: TRUE if `|dewma.before|` exceeds the
 prior row's `|dewma.after|`
 - `aft.g.aftm1`: TRUE if `|dewma.after|` exceeds the
 next row's `|dewma.before|`
3. Exclude the flagged member. Among multiple candidates,
 select by highest `|dewma.before|` (codes 1,3) or
 `|dewma.after|` (codes 2,4), with `internal_id`
 tiebreaker.
4. Remove excluded row, recalculate, repeat until no
 violations remain.

For **exactly 2 measurements:**
- Exclude the one with higher `abs(tbc.sd)`, with
 `internal_id` tiebreaker.
- No EWMA computation.

### Checklist findings

1. **Removed debug save point** (Step 17 output).
2. **Removed commented-out old code** (~30 lines of old
 exclusion logic and CSV debug output).
3. **Removed CP markers** (`### CP ADJUST`, `### CP ADD
 SECTION`, `### CP change`).
4. **Stale nnte comment** removed (was
 "Removed nnte_full filter").
5. **Boundaries:** `diff < mindiff` and `diff > maxdiff`
 — strict inequalities. WHO interval boundaries use
 `>=` / `<` consistently.
6. **`valid()` call:** Excludes temp SDEs, excludes weight.
7. **Sort order:** Explicit `order(agedays, id)` inside
 each while-loop iteration.
8. **Parameter scope:** HT and HC only. HC has separate
 velocity tables, tighter tolerances, and no Tanner.
9. **Factor levels:** `Exclude-C-{HT|HC}-Abs-Diff` exists in
 `exclude.levels`.
10. **Unused variables in Step 19 closure:** `abs_tbd.sd`,
 `abs_ctbd.sd`, `med_dop`, `med_cdop` are computed but
 never used (linter warnings). These appear to be
 leftover from prior refactoring. Not harmful but should
 be cleaned up.

---

## Step 19: Pairs and Singles Evaluation

| | |
|---|---|
| **Scope** | All parameters |
| **Operates on** | Include rows with 1–2 remaining measurements per subject-param |
| **Prior step** | Step 17 (Height/HC Velocity) |
| **Next step** | Step 21 (Error Load) |
| **Exclusion codes** | `Exclude-C-{WT|HT|HC}-Pair`, `Exclude-C-{WT|HT|HC}-Single` |
| **Code location** | See code |

### Overview

After all cleaning steps, some subject-params may have only
1–2 Include measurements remaining. Step 19 evaluates
whether these isolated values are plausible by comparing
them against the designated other parameter (DOP).

### DOP snapshot

A snapshot of all Include rows is taken BEFORE by-group
processing. This prevents a cross-parameter ordering issue
where excluding one param's values during processing could
make the DOP lookup fail for the other param.

### Pair evaluation

For 2 remaining measurements:
1. Compute `diff_tbc.sd` and `diff_ctbc.sd` (z-score
 difference between the two)
2. Compute DOP comparison: for each measurement, if a
 same-day DOP value exists, use `|dop_tbc - tbc|`;
 otherwise use `|median(dop_tbc) - tbc|`
3. Select which to exclude: highest `comp_diff` (or
 highest `|tbc.sd|` if both DOP comparisons are NA),
 with lowest `id` as tiebreaker
4. Apply exclusion:
 - `|diff_tbc| > 4 & (|diff_ctbc| > 4 or NA) &
 gap >= 365.25` → `Exclude-C-{WT|HT|HC}-Pair`
 - `|diff_tbc| > 2.5 & (|diff_ctbc| > 2.5 or NA) &
 gap < 365.25` → `Exclude-C-{WT|HT|HC}-Pair`
5. If one excluded, re-evaluate remaining as single

### Single evaluation

For 1 remaining measurement:
- `(|tbc.sd| > 3 & comp_diff > 5)` OR
 `(|tbc.sd| > 5 & no DOP data)`
 → `Exclude-C-{WT|HT|HC}-Single`

### Checklist findings

1. **Boundaries:** Pair: `> 4` / `> 2.5` strict; gap
 `>= 365.25` inclusive. Single: `> 3` / `> 5` strict.
2. **`valid()` call:** Excludes temp SDEs. Counts only
 Include rows for singles/pairs determination (matches
 Stata).
3. **DOP snapshot:** Correct — prevents cross-parameter
 interference.
4. **Both `tbc.sd` and `ctbc.sd` checked** for pairs
 (corrected must also agree when available).
5. **Parameter scope:** All 3 params handled uniformly.
6. **Factor levels:** All 3 codes exist.

---

## Step 21: Error Load

| | |
|---|---|
| **Scope** | All parameters |
| **Operates on** | All rows |
| **Prior step** | Step 19 (Pairs and Singles) |
| **Next step** | Step 22 (Output) |
| **Exclusion code** | `Exclude-C-{WT|HT|HC}-Too-Many-Errors` |
| **Code location** | See code |

### Overview

If a subject-param has a high ratio of excluded values
(excluding SDEs, CFs, and Missing from the denominator),
all remaining Include values are also excluded.

### Logic

1. `non_error_codes`: SDE-Identical, SDE-All-Exclude,
 SDE-All-Extreme, SDE-EWMA, SDE-One-Day,
 Carried-Forward, Missing — excluded from both
 numerator and denominator
2. `n_errors = count of excluded rows NOT in
 non_error_codes`
3. `n_includes = count of Include rows`
4. `err_ratio = n_errors / (n_errors + n_includes)`
5. If `err_ratio > error.load.threshold` (default 0.5)
 AND `n_errors >= error.load.mincount`
 (default 2): all Include rows →
 `Exclude-C-{WT|HT|HC}-Too-Many-Errors`

### Checklist findings

1. **Boundary:** `err_ratio > error.load.threshold` — strict
 `>`. An exactly 50% error rate does NOT trigger.
 **Bug fix (2026-04-12):** Was hardcoded `> .4` instead of
 using the `error.load.threshold` parameter (default 0.5).
 Fixed to use the parameter. This changes behavior: 2 rows
 in syngrowth that were excluded at 0.4 are now included
 at the correct 0.5 threshold.
2. **`error.load.mincount`:** Default 2. Prevents a
 subject-param with 1 error and 1 include (50% ratio)
 from triggering.
3. **Denominator excludes SDEs/CFs:** Correct — these
 are data structure artifacts, not cleaning errors.
4. **Parameter scope:** All 3 params, grouped by
 `(subjid, param)`.
5. **Factor level:** `Exclude-C-{WT|HT|HC}-Too-Many-Errors` exists.
6. **Stale nnte comment** removed.

---

## Step 22: Output Assembly

| | |
|---|---|
| **Code location** | See code |

### Overview

Assembles the return value from `cleanchild()`:

1. **`final_tbc`**: For potcorr
 subjects, returns `ctbc.sd`; for others, `tbc.sd`.
 Falls back to `tbc.sd` if potcorr/ctbc columns don't
 exist.

2. **Return columns**:
 - Essential: `id`, `line`, `exclude`, `param`,
 `cf_rescued`
 - Z-scores (if present): `sd.orig_who`, `sd.orig_cdc`,
 `sd.orig`, `tbc.sd`, `ctbc.sd`
 - EWMA1 iteration 1 diagnostics (if present):
 `ewma1_it1.*`
 - `final_tbc`

The caller (`cleangrowth()`) merges this with the original
input data and checkpoint diagnostics before returning to
the user.
