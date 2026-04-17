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
 against spec
3. **Uninitialized or unnecessary variables:** No variables
 used before assignment; no leftover variables from prior
 refactors
4. **Efficiency:** Caching opportunities, unnecessary
 merges, redundant computations, operations that could
 be vectorized
5. **Edge cases:** 0, 1, or 2 valid rows; empty groups;
 all-NA values — does code handle gracefully or crash?
6. **`.child_valid()` call correctness:** Which `include.*` flags
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

***Three growth parameters:*** The child algorithm cleans weight (WEIGHTKG), height/length (HEIGHTCM/LENGTHCM), and head circumference (HEADCM). LENGTHCM is relabeled to HEIGHTCM — only the `param` value changes; the measurement value is **not** altered. No adjustment is made for the ~0.5–0.7 cm difference between supine length and standing height. After relabeling, z-scores are calculated identically to HEIGHTCM. This is a known limitation to be addressed in a future update. HEADCM is only cleaned through age 3 years (marked "Exclude-Not-Cleaned" for agedays > 3 × 365.25); z-score data for HC is available only through age 5 years (marked "Exclude-Missing" for agedays ≥ 5 × 365.25).

***CSD z-scores, not measurements:*** The child algorithm bases most decisions on z-scores, which standardize measurements by age and sex. The algorithm uses the Conditional Standard Deviation (CSD) method for calculating z-scores (also referred to as SD-scores), not the standard LMS method. [ADD CDC REFERENCE] The LMS accounts for significant skewness in the distribution of weights, which is often useful, but in cleaning it reduces the ability of z-scores to distinguish large absolute differences in weight for heavier children, reducing cleaning accuracy. Several variations of z-scores (SD-scores) are used. These are described further below and are listed in the "### Z-score variants" section below. The two primary variations are sd.orig (the unaltered SD score) and tbc.sd. tbc.sd is calculated by subtracting a recentering value (described further below) from sd.orig.

***Dual z-scores for prematurity correction:*** One z-score variant is corrected z-scores. For subjects whose data are consistent with prematurity or SGA, a corrected z-score is calculated using values derived from the Fenton 2025 percentile reference. [ADD Fenton 2025 REF] We do not use the Fenton LMS parameters. Instead, we extracted approximate medians (M) and standard deviations from the plotted Fenton 2025 percentile curves for weight, length, and head circumference, by sex and gestational age in days. Separate upper and lower standard deviations (S_upper and S_lower) were extracted from the curves above and below the median, respectively. Values were smoothed to remove extraction artifacts. These are stored in `inst/extdata/fenton2025_ms_lookup_smoothed.csv` (long format: one row per sex × param × ga_days, with columns M, S_upper, S_lower). Z-scores are then calculated using the same CSD method used for the main WHO/CDC z-scores: `z = (measurement - M) / (S_upper × M)` when measurement >= M, and `z = (measurement - M) / (S_lower × M)` when measurement < M, where S_upper and S_lower are proportional (CV-like) standard deviations so that the absolute SD equals S × M. Corrected z-scores are used as an extra check to prevent exclusion of values that are plausible for a baby with very low birth weight.

***Recentering:*** After calculating z-scores, the algorithm subtracts a population-level median z-score for each parameter, sex, and age to create "to-be-cleaned" SD scores (tbc.sd). Because median z-score changes with age, using an unrecentered z-score can lead to bias in exclusions when comparing z-scores across long time periods. The recentering medians come from a fixed reference file (rcfile-2023-08-15_format.csv.gz) that was derived from data from outpatient records from the Children's Hospital of Philadelphia.

***Designated Other Parameter (DOP):*** Many steps compare a measurement against a different growth parameter for the same
subject. The designated other parameter (DOP) for weight is height; the DOP for height is weight; the DOP for HC is height. DOP comparisons are used in temporary SDE resolution (Step 5), CF rescue, EWMA steps, final SDE resolution (Step 13), and pairs/singles evaluation (Step 19). In Step 13, DOP medians are computed for one-day SDE resolution (cross-parameter comparison) and DOP is used as a secondary sort key for selecting which SDE to keep.

***Ordered steps:*** Like the adult algorithm, the child algorithm processes steps in a specific order, removing more
extreme problems first and getting more refined later. Same-day extraneous values (described below) are handled in repeated steps because of bidirectional impact of same-day selection and other exclusions. Some errors cause over- or under-detection of other errors, so they must be excluded first. The order of steps is central to the algorithm.

***Iterations within steps:*** Several steps (Evil Twins, EWMA1, EWMA2) are iterative. If more than one value is a candidate for exclusion, they select the most extreme candidate and then re-evaluate for further exclusions. This continues until no candidates remain or too few values are left to evaluate. This is a central part of the algorithm and reduces the likelihood of erroneous values leading to exclusion of true values.

***Temporary vs. permanent exclusions:*** Same-day extraneous (SDE) values are multiple values for the same parameter on the same day. growthcleanr selects the most plausible SDE for inclusion in the final dataset. SDEs are first flagged temporarily (Step 5) and later
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

***CF rescue:*** The child algorithm has a lookup-based
carried-forward rescue system (Step 6) that attempts to
distinguish genuinely repeated measurements from data-entry
artifacts. In the default `cf_rescue = "standard"` mode,
rescue thresholds vary by age bin, interval-from-originator
bin, parameter, and rounding type (whole/half imperial vs.
other); a CF is rescued when its absolute z-score difference
from the originator falls below the lookup cell's threshold.
`cf_rescue = "none"` and `cf_rescue = "all"` provide all-or-
nothing behavior. Rescued CFs are returned to `Include` with
the rescue reason stored in the separate `cf_rescued` column.
See Step 6 and `cf-rescue-thresholds.md` for details.

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

***The `.child_valid()` function:*** A critical gatekeeper that
determines which rows participate in each step. It returns a
logical vector indicating which rows are currently "valid"
(eligible for processing). Its `include.*` flags control
whether temporarily excluded categories (temp SDEs, permanent
SDEs, carried-forward values) are included:

- `.child_valid(df)` — only rows with `exclude` that does not start
 with "Exclude"
- `.child_valid(df, include.temporary.extraneous = TRUE)` — also
 includes "Exclude-C-Temp-Same-Day"
- `.child_valid(df, include.extraneous = TRUE)` — also includes
 "Exclude-C-Extraneous"
- `.child_valid(df, include.carryforward = TRUE)` — also includes
 "Exclude-C-CF"

`.child_valid()` works by checking the text of the `exclude`
column rather than using factor ordering. It uses
`grepl("^Exclude", ...)` to identify excluded rows, then
checks for specific codes to conditionally re-include the
temporary categories (`"Exclude-C-Temp-Same-Day"`,
`"Exclude-C-Extraneous"`, `"Exclude-C-CF"`). `"Include"` is the
only non-Exclude value the algorithm produces, so it passes
through unconditionally.

***Age-dependent internal_id tiebreaking:*** When resolving
same-day duplicates, the child algorithm uses an age-dependent
rule based on `internal_id` (the sequential integer assigned in
id-sorted order at session start): at birth (agedays == 0),
keep the lowest `internal_id` (earliest measurement, before
postnatal fluid shifts and interventions); at all other ages,
keep the highest `internal_id` (later measurement, which may
represent a more careful re-measurement). This differs from the
adult algorithm, which always keeps the highest `internal_id`.
The user's original `id` column is never used as a tiebreaker;
it is preserved untouched into output.

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
9. Batching and dispatch to `cleanchild()` (see "Batching and
 Dispatch" section below)
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

- `.child_valid()` — row eligibility filter
- `identify_temp_sde()` — temp SDE resolution
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

**System 1 — Outer wrapper (for memory management):**
- In `child_clean.R`, inside `cleangrowth()`.
- Creates batches of up to `batch_size` subjects
 (default 2000; configurable via the `batch_size` parameter):
 `patients[, batch := (seq_len(.N) - 1L) %/% batch_size + 1L]`.
- Iterates the batches with a `for` loop. Each iteration
 filters the input to the current batch's `subjid` values
 (both the pediatric `data.all` and the adult `data.adult`
 are limited to those subjects), runs preprocessing and the
 child/adult algorithms on just that batch, and pushes the
 per-batch result onto `results_list`.
- After the loop, `rbindlist(results_list)` combines
 per-batch results and the merge against `data.batch`
 restores user input columns.

**System 2 — Inner batching (for parallelism):**
- Also in `child_clean.R`, inside each outer-batch
 iteration.
- Assigns round-robin batch numbers to subjects in the
 already-filtered `data.all`:
 `(seq_along(subjid.unique) - 1L) %% num.batches + 1L`.
- When `parallel = FALSE`, `num.batches = 1`, making this a
 no-op: all subjects in the outer batch go to a single
 `cleanchild()` call.
- When `parallel = TRUE`, dispatches via
 `ddply(..., .parallel = TRUE)` to process inner batches in
 parallel across CPU cores.

The inner batching subdivides whatever data the outer wrapper
provides into parallel work units. With `parallel = FALSE`,
it has no effect. The adult algorithm has its own inner
parallel batching system, structured identically.

### Parallel processing

`parallel = TRUE` runs the inner batches concurrently. It
requires the growthcleanr package to be installed (not just
`devtools::load_all()`-ed) because the worker processes need
`system.file()` access to the `inst/extdata/` reference
tables. Results are bit-identical to sequential processing
(verified on the 77,721-row `syngrowth` dataset). Typical
observed speedups on 2,697 subjects: ~114 s sequential, ~65 s
at `num.batches = 2` (1.7×), ~30 s at `num.batches = 4`
(3.9×).

Batch-invariant setup (the `exclude.levels` factor and the
Tanner / WHO velocity reference files) is performed once
outside the outer loop so the cost is not repeated per batch.

### Partial runs and preloaded references

Two `cleangrowth()` parameters support workflows that call
the function many times on overlapping data (simulation
loops, error-injection pipelines, secure-environment
workflows with updated data extracts). Both are optional and
independent; either or both can be used.

**`ref_tables` + `gc_preload_refs()`.** `gc_preload_refs()`
reads the CDC and WHO growth-chart reference files once and
returns a named list of closures (`mtz_cdc_prelim` and
`mtz_who_prelim`). Passing that list as the `ref_tables`
argument to subsequent `cleangrowth()` calls skips the
per-call `read_anthro()` disk reads (≈0.9 seconds per call).
The return value is a session-lifetime object; rebuild only
if the reference files themselves change.

```r
refs <- gc_preload_refs()
res  <- cleangrowth(..., ref_tables = refs)
```

**`cached_results` (+ optional `changed_subjids`).** When
rerunning `cleangrowth()` on input that differs from a prior
run for only a subset of subjects, pass the prior result as
`cached_results`. Only subjects whose input rows differ are
re-processed; unchanged subjects are copied from the cache.
Subjects are independent in all cleaning operations
(by-subject grouping, fixed-reference recentering), so
partial runs produce the same output a full run would for
the subjects that changed, and bit-identical output for the
ones that didn't.

Two modes:

- **Auto-detect** (preferred — `cached_results` provided,
 `changed_subjids = NULL`): `cleangrowth()` compares the
 incoming `subjid`, `param`, `agedays`, `sex`, and
 `measurement` against the cache (with the same `0 → NaN`
 transformation applied internally) and re-processes any
 subjects with added, removed, or modified rows.
- **Explicit** (`cached_results` + `changed_subjids`
 provided): only the listed subjects are re-processed; the
 comparison step is skipped.

With `quietly = FALSE`, the wrapper prints a one-line
summary: `Auto-detected N changed subjects (X added, Y
modified, Z removed, W unchanged)`. Output row order matches
input order (sorted by the session-assigned `internal_id`).

Observed runtimes (200 subjects, 50 000 calls, 5 % of
subjects changed per call): ≈0.40 s/call vs ≈3.4 s/call for
a full run with preloaded refs — roughly an 11× speedup and
~56 hours saved over the simulation. See
`gc-github-latest/CLAUDE.md` for the full benchmark table.

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
| `fent_foraga.csv.gz` | Fenton | Weight → est. gestational age |
| `fenton2025_ms_lookup_smoothed.csv` | Fenton 2025 | M, S_upper, S_lower for CSD z-scores |
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
(half away from zero) for cross-platform validation; this
is no longer needed.

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

This matches the main z-score blending formula.

### Recentering

After z-score calculation, the algorithm subtracts age-,
sex-, and parameter-specific median z-scores to produce
tbc.sd:

```
tbc.sd = sd.orig - sd.median
ctbc.sd = sd.corr - sd.median
```

Recentering uses a precomputed reference file
(`inst/extdata/rcfile-2023-08-15_format.csv.gz`) whose
medians were built once via the procedure implemented in
`sd_median()` (year-of-age medians treated as midyear-age,
linearly interpolated by day of age, clamped outside the
covered range). Callers can override by passing a
`sd.recenter` data.table to `cleangrowth()`.

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
`.child_valid()` function. When the narrative says a value is
"excluded," this means:

1. The value's `exclude` column is set to the appropriate
 exclusion code
2. Subsequent calls to `.child_valid()` will exclude this row
 from processing
3. The value remains in `data.df` with its exclusion code

There is also a softer mechanism: the **temporary SDE**
flag. Values marked `"Exclude-C-Temp-Same-Day"`
are excluded by default `.child_valid()` calls but can be included
by passing `include.temporary.extraneous = TRUE`. This
allows steps to optionally include or exclude temp SDEs.
Permanent SDE resolution (Step 13) replaces temp SDE codes
with final codes.

### Key difference from adult algorithm

The adult algorithm copies rows into separate working
dataframes (`h_subj_df`, `w_subj_df`) that shrink as values
are excluded. The child algorithm keeps all rows in a single
`data.df` and uses the `exclude` column plus `.child_valid()` to
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
| `id` | any (required) | Input | User-provided row identifier. Required (errors if `NULL`). Preserved untouched into output; not used as a sort key or tiebreaker by the algorithm. |
| `internal_id` | integer | `cleangrowth()` | Sequential integer `1:N` assigned in id-sorted order at session start. Used throughout as the final sort-tiebreaker and for named-vector indexing where applicable (age-dependent tiebreaking keeps the lowest `internal_id` at `agedays == 0`, the highest at `agedays > 0`). |
| `index` | integer | `cleanchild` | Sequential `1:.N` within each outer-batch dataframe, reassigned at the top of `cleanchild()` so merges and subset operations within the batch have contiguous keys. |
| `subjid` | any | Input | Subject identifier (any type; preserved as-is into output). |
| `param` | character | Input (post-conversion) | `"WEIGHTKG"`, `"HEIGHTCM"`, or `"HEADCM"` after imperial→metric relabeling. |
| `agedays` | integer | Input | Age in days at the measurement. |
| `sex` | integer | Input | 0 = male, 1 = female. |
| `v` | numeric | Preprocessing | Working measurement value (metric; 0 → NaN). |
| `v.orig` | numeric | Step 5 | Copy of `v` captured at the start of Step 5, used in Step 6 CF detection. |

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

All of these are internal Step 6 working columns; all are
dropped before Step 6 exits (via `cols_to_drop_6` or earlier
cleanup).

| Variable | Type | Role |
|----------|------|------|
| `v.orig` | numeric | Copy of the raw measurement (created in Step 5, used in Step 6 detection). |
| `wholehalfimp` | logical | TRUE if the measurement is in whole or half imperial units. Drives the rounding dimension of the rescue lookup. |
| `cf_binary` | logical | TRUE if the row is currently `Exclude-C-CF`. Temp, cleaned mid-step. |
| `originator` | logical | TRUE for an Include whose next row (subject/param/ageday sort order) is a CF. Temp, cleaned mid-step. |
| `originator_seq` | integer | `cumsum(originator)` per subject-param. Temp, cleaned mid-step. |
| `cf_string_num` / `cs` | integer | Sequential string number per subject-param. `cs` is the alias retained through the rescue phase. |
| `originator_z` | numeric | The originator's `sd.orig_uncorr` z-score, propagated forward to its CFs. Temp, cleaned mid-step. |
| `seq_win` | integer | Position within a CF string: `0` for originator, `1, 2, 3, ...` for subsequent CFs. |
| `absdiff` | numeric | \|sd.orig_uncorr − originator_z\| for each CF. |
| `ageday_has_include` | logical | TRUE if any row on the same (subjid, param, agedays) has `exclude == "Include"`. Restricts CF rescue eligibility (not originator assignment). |
| `orig_ageday` | integer | The originator's agedays, propagated forward to its CFs. |
| `cf_interval` | integer | `agedays − orig_ageday` for each CF. Drives the interval dimension of the rescue lookup. |

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
| `Exclude-C-Temp-Same-Day` | 5 | All | Temp SDE (may persist if not resolved by Step 13) |
| `Exclude-C-CF` | 6 | All | Value identical to prior-day value (not rescued) |
| `Exclude-C-BIV` | 7 | All | Outside absolute or standardized biological limits |
| `Exclude-C-Evil-Twins` | 9 | All | Adjacent extreme value pair/group |
| `Exclude-C-Traj-Extreme` | 11 | All | Extreme EWMA outlier |
| `Exclude-C-Identical` | 13 | All | Same-day duplicate with identical value |
| `Exclude-C-Extraneous` | 13 | All | Same-day extraneous (resolved by EWMA, all extreme, or single day) |
| `Exclude-C-Traj` | 15/16 | All | Moderate EWMA outlier (all sub-rules collapsed) |
| `Exclude-C-Abs-Diff` | 17 | HT, HC | Height/HC velocity exceeded allowed limits |
| `Exclude-C-Pair` | 19 | All | 2 measurements, one excluded |
| `Exclude-C-Single` | 19 | All | Single measurement excluded |
| `Exclude-C-Too-Many-Errors` | 21 | All | Error ratio exceeds threshold |

Note: Exclusion codes are not parameter-specific — the param
is in the data row (e.g., `Exclude-C-CF` applies to WT, HT, and
HC). The `exclude.levels` factor defined in `cleangrowth()`
enumerates exactly the codes listed here, plus the adult
codes; values outside this list would silently become `NA` when
assigned to the `exclude` factor, so new codes must be added
to `exclude.levels` before use.

### CF rescue reason codes (cf_rescued column)

| Code | Meaning |
|------|---------|
| `""` (empty) | Not a CF candidate, or CF not rescued (the row stays `Exclude-C-CF`) |
| `Rescued` | CF rescued in standard mode — \|Δz\| below the lookup threshold for the row's (age bin, interval bin, param, rounding type) cell |
| `Rescued-All` | CF rescued because `cf_rescue = "all"` |

Rescued CFs have their `exclude` column set back to
`"Include"`; the rescue reason is preserved only in the
`cf_rescued` column.

---

## Configurable Parameters

The table is audited incrementally by walk-through session.
Rows for parameters used in steps that have not yet been
walked are shown with their current defaults but without a
detailed description; those rows will be filled in and
verified against the code when each step is walked.

| Parameter | Default | Used in | Description |
|-----------|---------|---------|-------------|
| `id` | NULL | Preprocessing | Required row identifier preserved into output. Can be any type (numeric, character, UUID). If `NULL`, `cleangrowth()` errors. |
| `adult_cutpoint` | 20 | Preprocessing | Age (years) dividing the pediatric and adult paths. Subjects with measurements spanning the cutpoint have their pediatric and adult rows cleaned by separate algorithms. |
| `sd.recenter` | NA | Preprocessing (recentering) | If left `NA`, the built-in rcfile is used (`inst/extdata/rcfile-2023-08-15_format.csv.gz`). Pass a data.table to override with custom recentering medians. |
| `ref.data.path` | `""` | Preprocessing | Directory override for reference files (CDC/WHO z-score tables, rcfile, velocity tables, Fenton data). Default empty string means use the installed package's `inst/extdata/`. |
| `batch_size` | 2000 | Batching (outer wrapper) | Number of subjects per outer-wrapper batch. Subjects are never split across batches. |
| `parallel` | FALSE | Batching (inner) | If `TRUE`, inner batches dispatch via `ddply(..., .parallel = TRUE)` for concurrent processing. Requires the package to be installed (not just `load_all()`-ed). |
| `num.batches` | NA | Batching (inner) | Number of inner (parallel) batches. Auto-set when `parallel = TRUE`; ignored when `parallel = FALSE` (acts as 1). |
| `log.path` | NA | Batching (inner) | If non-NA and `parallel = TRUE`, each worker sinks messages to `{log.path}/cleangrowth-{date}-batch-{n}.log`. |
| `ref_tables` | NULL | Performance | Pre-loaded reference closures from `gc_preload_refs()`. Skips repeated disk reads across many `cleangrowth()` calls; see "Partial runs and preloaded references" below. |
| `cached_results` | NULL | Performance | Prior `cleangrowth()` output. Unchanged subjects are copied from cache; changed subjects are re-processed. Either auto-detected or listed explicitly via `changed_subjids`. See "Partial runs and preloaded references" below. |
| `changed_subjids` | NULL | Performance | Optional explicit list of subject IDs to re-process when using `cached_results`. If `NULL`, changed subjects are auto-detected by comparing input to cache. |
| `debug` | FALSE | Diagnostics | If `TRUE`, emits 6 `ewma1_it1.*` columns capturing first-iteration Step 11 EWMA values. Off by default to keep standard output lean. |
| `tri_exclude` | FALSE | Output | If `TRUE`, adds a `tri_exclude` column with a three-level summary (Include / Exclude / N/A). Convenience for downstream filtering. |
| `quietly` | TRUE | All | Suppress progress messages. When `FALSE`, the wrapper prints batch and step progress via `message()`. |
| `include.carryforward` | FALSE | Step 6 | **Deprecated** — mapped to `cf_rescue = "all"` with a warning. Use `cf_rescue` directly. |
| `cf_rescue` | `"standard"` | Step 6 | See Step 6. |
| `cf_detail` | FALSE | Step 6 | See Step 6. |
| `weight_cap` | (missing) | Adult | **Deprecated** — mapped to `adult_scale_max_lbs` with a warning. |
| `adult_permissiveness` | `"looser"` | Adult | Not walked yet; see `adult-algorithm-narrative.md`. |
| `adult_scale_max_lbs` | `Inf` | Adult | Not walked yet. |
| `biv.z.wt.low.young` | -25 | Step 7 | Lower CSD z cutoff for WEIGHTKG at `ageyears < 1`. |
| `biv.z.wt.low.old` | -15 | Step 7 | Lower CSD z cutoff for WEIGHTKG at `ageyears >= 1`. |
| `biv.z.wt.high` | 22 | Step 7 | Upper CSD z cutoff for WEIGHTKG (all ages). |
| `biv.z.ht.low.young` | -25 | Step 7 | Lower CSD z cutoff for HEIGHTCM at `ageyears < 1`. |
| `biv.z.ht.low.old` | -15 | Step 7 | Lower CSD z cutoff for HEIGHTCM at `ageyears >= 1`. |
| `biv.z.ht.high` | 8 | Step 7 | Upper CSD z cutoff for HEIGHTCM (all ages). |
| `biv.z.hc.low` | -15 | Step 7 | Lower CSD z cutoff for HEADCM (all ages). |
| `biv.z.hc.high` | 15 | Step 7 | Upper CSD z cutoff for HEADCM (all ages). |
| `ewma_window` | 15 | EWMA steps | Not walked yet. |
| `error.load.mincount` | 2 | Step 21 | Not walked yet. |
| `error.load.threshold` | 0.5 | Step 21 | Not walked yet. |

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

Step numbers are not consecutive. Steps 3 (unit error),
4 (swaps), 8, 10, 12, 14, 18, and 20 either do not exist
or are handled within other steps in the child algorithm.

---

## Step 2b: Gestational Age Correction (preprocessing)

| | |
|---|---|
| **Scope** | Subjects whose first weight suggests prematurity (potcorr subjects) |
| **Operates on** | All rows of potcorr subjects |
| **Prior step** | Z-score calculation (WHO/CDC) |
| **Next step** | Recentering |
| **Output columns** | `sd.corr`, `uncorr`, `potcorr` (merged back onto all rows) |
| **Exclusion code** | None — computes corrected z-scores only |
| **Code location** | `cleangrowth()` in `child_clean.R` ~lines 737–1005 |
| **Controlled by** | No user parameters; uses built-in Fenton 2025 reference |

### Overview

Very preterm subjects whose anthropometric references should
be preterm growth standards are identified as "potentially
correctable" (`potcorr`) and have their z-scores recomputed
against the Fenton 2025 reference using postmenstrual age. If
the correction happens to worsen the early trajectory
(corrected first-to-fourth-weight deviations sum larger than
uncorrected), the correction is reverted for that subject
(`uncorr := 1`, `sd.corr := sd.orig`). All three columns
(`sd.corr`, `uncorr`, `potcorr`) are written back to
`data.all` for every row; downstream steps read them through
`final_tbc` and related recentered fields.

### Phase 1: potcorr identification

Subjects are flagged by a two-level rule.

- **Sort.** Sort by `(subjid, param, agedays, id_sort)` where
 `id_sort = internal_id` at `agedays == 0` and `-internal_id`
 otherwise. At birth this keeps the earliest-recorded row as
 the "first" weight; at later ages it keeps the most recent
 internal_id — the age-dependent tiebreaker used throughout
 the algorithm.
- **`potcorr_wt`** (row-level): a weight row is
 `potcorr_wt = TRUE` iff (a) it is the first weight row for
 that subject under the sort above (`seq_len(.N) == 1L` with
 `by = subjid`), (b) `sd.orig < -2`, and (c) age is less than
 10 months. NA-safe: `sd.orig` is `NA` for `Exclude-Missing`
 rows so they never flag.
- **`potcorr`** (subject-level): `any(potcorr_wt, na.rm = TRUE)`
 per subject. Propagated to all rows of that subject.

### Phase 2: fast path when no potcorr

`has_potcorr <- any(data.all$potcorr, na.rm = TRUE)`. If
FALSE, Step 2b short-circuits: `sd.corr := sd.orig`,
`uncorr := 0`. Fenton reference files are not read. This
typical-case fast path matters because potcorr subjects are
usually < 1 % of a dataset.

### Phase 3: Fenton correction (potcorr subjects only)

All work happens on a copy `pc <- data.all[subjid %in% pc_ids]`.

1. **Integer weight.** `intwt := trunc(v * 100) * 10` (grams
 rounded to 10 g); values in `[250, 500)` are floored to 500
 to hit the Fenton table's minimum.
2. **Fenton merge #1 — weight → GA.** Merge `pc` with
 `fent_foraga` (from
 `inst/extdata/fent_foraga.csv.gz`) on `(sex, intwt)` to get
 `fengadays` (Fenton gestational age in days for that
 weight).
3. **Propagate.** Use the qualifying first-weight row's
 `fengadays` for every row of the subject
 (`fengadays_subj := min(fengadays[potcorr_wt], ...)`, by
 subjid).
4. **Ages.** `pmagedays := agedays + fengadays` (postmenstrual
 age); `cagedays := pmagedays - 280` (corrected age). Both
 computed only when `fengadays < 259` (term cutoff).
5. **Fenton merge #2 — GA → M/S.** Map gc param names
 (`WEIGHTKG`/`HEIGHTCM`/`HEADCM`) to Fenton param names
 (`weight`/`length`/`headcirc`), then merge with
 `fenton2025_ms_lookup_smoothed.csv` on `(sex, fengadays,
 param)` to pull `M`, `S_upper`, `S_lower`.
6. **Fenton CSD z-score.** `v_fenton` is `v * 1000` for weight
 and `v` otherwise. `unmod_zscore := (v_fenton − M) /
 (S_upper * M)` if `v_fenton >= M`, else `/ (S_lower * M)`.
 This uses the same CSD method as WHO/CDC elsewhere in the
 algorithm.
7. **Correction failure.** If Fenton merge failed (no M),
 reset `potcorr_wt` to FALSE for that row and recompute the
 subject-level `potcorr`. A subject whose weight is outside
 Fenton's range silently falls out of the correction.
8. **Initial `sd.corr`.** Rows with `unmod_zscore` non-NA get
 `sd.corr := unmod_zscore`. Rows still `NA` fall back to
 `sd.orig`.
9. **Corrected WHO/CDC z-scores.** Compute `sd.c_who` and
 `sd.c_cdc` via the in-memory `measurement.to.z_who` /
 `measurement.to.z` closures at `cagedays`. For HEIGHTCM
 between ages 2 and 2+ at Fenton-corrected-age-≤2, apply the
 supine/standing offsets (+0.8 cm for WHO, +0.7 cm for CDC)
 before the z-score call.
10. **`sd.c` (corrected blend).** HEADCM uses WHO at all ages.
 Other params: WHO for `ageyears_2b ≤ 2`, smoothed
 WHO→CDC blend for `2 < ageyears_2b < 5`, CDC for
 `ageyears_2b ≥ 5`.
11. **Fenton takes precedence under 2.** For `ageyears_2b ≤ 2`
 with Fenton available, `sd.c := unmod_zscore`.
12. **Final `sd.corr`.** Age-dependent:
 - `ageyears_2b ≤ 2` with potcorr: `sd.c`
 - `2 < ageyears_2b ≤ 4`: smoothed blend
 `(sd.orig * (4 − ageyears_2b) + sd.c * (ageyears_2b − 2)) / 2`
 - otherwise: `sd.orig`

### Phase 4: uncorr check (revert when correction worsens trajectory)

Only weight rows of potcorr subjects participate. A
temporary `tmp` copy is built:

- SDE-Identicals are filtered via the age-dependent
 `keep_id` rule (match Early Step 13 behavior).
- `seq_win := sequence(.N)` per subject; keep rows with
 `seq_win ≤ 4 AND ageyears_2b < 2`.
- Drop subjects that have only one value or no non-NA
 z-scores.
- Compute `sd.corr_abssumdiff := |sum(sd.corr[1] − sd.corr)|`
 and `sd.orig_abssumdiff := |sum(sd.orig[1] − sd.orig)|`
 per subject.
- A subject's correction is reverted if
 `sd.corr_abssumdiff > sd.orig_abssumdiff` (and the first
 weight exists): in `pc`, `sd.corr := sd.orig`,
 `uncorr := 1`.

### Phase 5: merge results back

`data.all[, sd.corr := sd.orig]` and `uncorr := 0L` baseline,
then merge the pc results back by `id`:

```r
data.all[pc_result, on = "id", `:=`(
  sd.corr = i.sd.corr,
  uncorr = i.uncorr,
  potcorr = i.potcorr
)]
```

Temporary columns (`agemonths`, `ageyears_2b`, `potcorr_wt`,
etc.) are dropped before the data flows into recentering.

### Checklist findings (2026-04-17 walk)

1. **Sort determinism.** Age-dependent sort key
 (`id_sort`) is assigned explicitly before potcorr_wt is
 marked.
2. **Fenton reference files.** Loaded only under
 `has_potcorr`; `fent_foraga.csv.gz` and
 `fenton2025_ms_lookup_smoothed.csv` both in
 `inst/extdata/`.
3. **Fast-path correctness.** When no subjects flag, the
 fast path writes `sd.corr := sd.orig` and `uncorr := 0L`
 — every row still has a non-missing `sd.corr` for
 downstream use.
4. **NA safety.** The `potcorr_wt` predicate uses
 `!is.na(sd.orig) & sd.orig < -2`, so `Exclude-Missing`
 rows (where `sd.orig` is `NA`) cannot flag.
5. **Fenton merge failure.** Correctly resets `potcorr_wt`
 to FALSE and recomputes subject-level `potcorr` so the
 subject is excluded from downstream Fenton logic without
 silently using partial results.
6. **Correction revert rule.** The abssumdiff comparison is
 computed only over weight rows with `ageyears_2b < 2` and
 `seq_win ≤ 4`. The check is first-weight-only
 (`is_first`), which matches its purpose (the correction
 is driven by the first qualifying weight).
7. **HEIGHTCM supine/standing offsets** (+0.8 WHO, +0.7
 CDC) are applied only when `agedays > 730 & cagedays ≤
 730` — the band where chronological age is post-2y but
 corrected age is still under 2y.

---

## Step 13: SDE Resolution (including Early Step 13)

| | |
|---|---|
| **Scope** | All parameters |
| **Exclusion codes** | `Exclude-C-Identical`, `Exclude-C-Extraneous` |
| **Code location** | Early Step 13 is inline at the top of `cleanchild()` in `child_clean.R` (~lines 2518–2537). Main Step 13 is later in `cleanchild()` and reuses `identify_temp_sde()` (`child_clean.R` ~lines 2021+). |

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
4. Others → `Exclude-C-Identical`

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
3. Re-run `identify_temp_sde()` with
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
 `Exclude-C-Extraneous`. Threshold: strict `> 2`.

2. **DOP medians**: For one-day SDEs,
 cross-parameter medians provide a secondary sort key.
 DOP mapping: WT↔HT, HC→HT. Uses only fully included
 values (`!was_temp_sde`).

3. **One-Day selection**: Sort by:
 - `absdiff_rel_to_median` ascending
 - `absdiff_dop_for_sort` ascending (Inf if NA)
 - Age-dependent `internal_id` tiebreaker
 Keep first; others → `Exclude-C-Extraneous`

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
 all → `Exclude-C-Extraneous`. Threshold: strict
 `> 1`.
6. **SDE-EWMA selection**: Sort eligible
 values by `absdewma` ascending, age-dependent
 `internal_id` tiebreaker. Keep first → Include; others
 → `Exclude-C-Extraneous`

**Phase B4: Merge back**

SDE results are merged back to `data.df` by `id`. Only
rows that were Include or Temp-SDE are overwritten —
permanent exclusions from prior steps (BIV, Evil Twins,
etc.) are preserved. Extra columns from SDE processing
are dropped.

### Checklist findings

1. **DOP mapping consistent:** HC → HT in all three DOP
 median calculations.
2. **Boundaries:** SDE-All-Extreme `> 2` (one-day) and
 `> 1` (EWMA-based) — both strict. Age-dependent id
 tiebreaker uses `internal_id` (not the user's `id`) in
 Step 13.
3. **`.child_valid()` calls:** Temp SDEs included when
 building the `data.sde` subset. EWMA computed on
 Include-only rows (excludes temp SDEs). DOP medians
 exclude temp SDEs via `!was_temp_sde`.
4. **Merge safety:** Only Include and Temp-SDE rows are
 overwritten on merge-back — permanent exclusions from
 earlier steps are preserved.
5. **Parameter scope:** All three params are handled. DOP
 has a param-specific mapping (see get_dop()).
6. **Factor levels:** All 5 SDE exclusion codes exist in
 `exclude.levels`.

---

## Step 5: Temporary SDE Resolution

| | |
|---|---|
| **Scope** | All parameters |
| **Operates on** | Valid rows (via `.child_valid()`, includes temp SDEs) |
| **Prior step** | Early Step 13 (SDE-Identicals) |
| **Next step** | Step 6 (Carried Forwards) |
| **Exclusion code** | `Exclude-C-Temp-Same-Day` |
| **Code location** | Inline in `cleanchild()` at `child_clean.R` ~lines 2540–2556; dispatches to `identify_temp_sde()` (`child_clean.R` ~lines 2021+) without the Step-13 `exclude_from_dop_ids` argument. |

### Overview

After identical same-day values have been removed, remaining
same-day duplicates (SDEs with **different** values) need
temporary resolution. Step 5 selects the most plausible value
for each same-day group and temporarily excludes the others.
These temporary exclusions are revisited in Step 13 (final
SDE resolution) after the dataset has been further cleaned.

The function `identify_temp_sde()` is also reused
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
3. Age-dependent `internal_id` tiebreaker:
 - At birth (`agedays == 0`): `internal_id` ascending (keep lowest)
 - At other ages: `internal_id` descending (keep highest)

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
indicate a bug in `identify_temp_sde()` —
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
| **Operates on** | Valid rows (including temp SDEs), restricted to subject-parameter combinations with at least one duplicate `v.orig` value |
| **Prior step** | Step 5 (Temporary SDEs) |
| **Next step** | Step 7 (BIV) |
| **Exclusion code** | `Exclude-C-CF` |
| **Rescue codes** | `Rescued`, `Rescued-All` (in `cf_rescued` column) |
| **Code location** | `child_clean.R` ~lines 2565–2880 (CF logic), support helpers `.cf_rescue_lookup()` and `.cf_get_thresholds()` at ~lines 2297–2400 |
| **Controlled by** | `cf_rescue` parameter (`"standard"` default / `"none"` / `"all"`); `cf_detail` parameter for optional diagnostic output columns |

### Overview

A carried-forward (CF) value is a measurement identical to
the measurement on the immediately prior day for the same
subject and parameter. CFs typically represent data-entry
artifacts (copying the previous value forward) rather than
independent re-measurements. Step 6 identifies CFs, organizes
them into "strings" that start with a non-CF originator and
continue through consecutive CFs, and then either leaves them
excluded or rescues them back to `Include` based on the
selected `cf_rescue` mode.

Three rescue modes:

- `cf_rescue = "standard"` (default): use the age ×
 interval × param × rounding deltaZ lookup implemented in
 `.cf_rescue_lookup()` / `.cf_get_thresholds()`. A CF is
 rescued when |Δz| = |sd.orig_uncorr − originator_z| is
 below the lookup cell's threshold. Full threshold tables
 are in `cf-rescue-thresholds.md`.
- `cf_rescue = "none"`: no rescue. Every detected CF stays
 `Exclude-C-CF`.
- `cf_rescue = "all"`: every detected CF is rescued, even
 CFs on a same-day SPA as another non-CF Include. Step 13
 final-SDE resolution handles the resulting multi-Include
 SPAs.

An optional `cf_detail = TRUE` adds two diagnostic columns
to the output:

- `cf_status`: `NA` for non-candidates, `"CF-Resc"` for
 rescued CFs, `"CF-NR"` for detected CFs that were not
 rescued.
- `cf_deltaZ`: |sd.orig_uncorr − originator_z| for every
 CF candidate.

### Phase 1: CF detection

**Pre-filter.** Detection runs only on subject-parameter
combinations that have (a) more than one valid measurement
(`.child_valid(data.df, include.temporary.extraneous = TRUE)`
returns TRUE) and (b) at least one duplicate `v.orig` value
(`uniqueN(v.orig) < .N`). Subject-params without duplicates
skip CF detection entirely; their CF eligibility is logged.

**Detection logic.** For each (subjid, param) that passes the
pre-filter:

1. Count values per (subjid, param, agedays) — days with
 multiple values (SDEs) are identified.
2. For each unique ageday per subject-param, look up the
 prior unique ageday with `shift(agedays, type = "lag")`.
3. Extract the prior ageday's value, but only if that day
 had exactly one measurement. Days with multiple values
 set `single_val` to `NA`, which disables CF matching.
4. A row is CF when `!is.na(prior_single_val) & v.orig ==
 prior_single_val` — exact equality, no tolerance.

Detected CFs have `exclude` set to `Exclude-C-CF` via
`.child_exc(param, "CF")`.

**Why only single-value prior days?** Prevents false CF
detection when a prior day has an SDE group whose values
happen to span the current day's measurement.

### Phase 2: Temp SDE re-evaluation

If any rows are still `Exclude-C-Temp-Same-Day` after Phase
1, their Step-5 selection may be stale now that some
originally Include values on the same day have been
excluded as CFs. Step 6 re-runs temp SDE resolution:

1. Reset every `Exclude-C-Temp-Same-Day` row back to
 `Include`.
2. Re-run `identify_temp_sde()` on the updated data — the
 same function used in Step 5. Because CF-excluded rows
 are not valid, temp SDE selection now accounts for the
 CF removals.

### Phase 3: CF string construction

**3a. `wholehalfimp` flag.** Row-level indicator of whether a
measurement is in whole or half imperial units; drives the
rounding dimension of the rescue lookup.

- `WEIGHTKG`: `abs((v.orig * 2.20462262) %% 1) < 0.01`
 (whole pounds)
- `HEIGHTCM`: `abs((v.orig / 2.54) %% 0.5) < 0.01` (whole
 or half inches)
- `HEADCM`: same formula as `HEIGHTCM`

**3b. SDE-Identical temporary removal.** Rows flagged
`Exclude-C-Identical` in Early Step 13 would break the
positional string logic below (adjacent-row `shift()`
operations, originator detection). They are moved out of
`data.df` into `sde_identical_rows` and re-added in Phase 5
after rescue.

**3c. `ageday_has_include`.** Per (subjid, param, agedays),
TRUE if any row on that day has `exclude == "Include"`. CFs
on days that also have an Include are not eligible for
rescue — rescuing them would create multiple Includes on the
same day. Temp SDEs are treated as non-Include here, so they
do not block CF rescue.

**3d. Positional string detection.** Per (subjid, param):

- `cf_binary`: TRUE for rows with `exclude ==
 "Exclude-C-CF"`.
- `originator`: TRUE for `Include` rows whose immediate next
 row (in subject/param/ageday sort order) is `cf_binary ==
 TRUE`. Any Include can be an originator regardless of its
 own `ageday_has_include` value.
- `cf_string_num` (alias `cs`): sequential string number
 per subject-param; `cumsum(originator)` assigns the
 originator's number, then a loop propagates it forward
 through consecutive CFs. Propagation reaches a CF only
 when `ageday_has_include == FALSE`, so CFs sharing a day
 with an Include never receive a string number and are not
 eligible for rescue.
- `originator_z`: the originator's `sd.orig_uncorr`
 (pre-GA-correction) z-score, propagated forward alongside
 `cf_string_num`.
- `seq_win`: position within the string — `0` for the
 originator, `1, 2, 3, ...` for subsequent CFs.
- `absdiff`: |sd.orig_uncorr − originator_z| for each CF.
- `orig_ageday`: the originator's agedays, propagated
 forward to the string's CFs so that `cf_interval =
 agedays − orig_ageday` (interval from originator, in
 days) can be used as the interval bin for the rescue
 lookup.

### Phase 4: Rescue evaluation

Rescue mode behavior:

**`cf_rescue = "none"`.** No action; every detected CF
stays `Exclude-C-CF`.

**`cf_rescue = "all"`.** Every row with `exclude ==
"Exclude-C-CF"` is re-included: `cf_rescued := "Rescued-All"`,
`exclude := "Include"`. Unlike standard mode, `"all"` does
not apply the `ageday_has_include` check — a CF that shares
an SPA with another Include is still rescued. This preserves
the caller's "ignore CFs" intent: a CF consistent with the
trajectory should not be preferentially excluded just
because another value landed on the same day. Any multi-
Include SPAs produced this way are resolved by Step 13 final
SDE.

**`cf_rescue = "standard"`** (default). Rescue considers
only CFs that are eligible — `seq_win > 0` and
`ageday_has_include == FALSE`. For each eligible CF,
`.cf_get_thresholds()` looks up a deltaZ threshold using:

- **Age bin** (from `agedays`): 0–3 mo, 3–6 mo, 6–12 mo,
 1–2 y, 2–5 y, 5–10 y, 10–15 y, 15–20 y.
- **Interval bin** (from `cf_interval`): <1 wk,
 1 wk – 1 mo, 1–6 mo, 6 mo – 1 yr, >1 yr.
- **Param** (`HEIGHTCM`, `WEIGHTKG`, `HEADCM`).
- **Rounding type** (`wholehalfimp` + age ≥ 2 y →
 `imperial`, otherwise `other`).

The CF is rescued when the lookup cell's threshold is
non-NA, greater than zero, and `absdiff < threshold`. Cells
with threshold `0` encode "no rescue" (NR); cells with
threshold `NA` encode an impossible combination (e.g.
imperial rounding for an age bin that shouldn't plausibly
encounter it). Rescued CFs set `cf_rescued := "Rescued"`
and `exclude := "Include"`.

Full threshold tables and the rationale for each cell are
in `cf-rescue-thresholds.md`; the matrices themselves are
in `.cf_rescue_lookup()` in `child_clean.R`.

**HEADCM placeholder.** `.cf_rescue_lookup()` currently
reuses HEIGHTCM's matrices for HEADCM pending HC-specific
analysis (see the inline `# placeholder: use HT thresholds`
comment in the helper).

### Phase 5: SDE-Identical restoration, cf_detail, cleanup

1. **SDE-Identical restoration.** If `sde_identical_rows`
 has rows, rbind them back into `data.df` (`fill = TRUE`)
 and re-sort by (subjid, param, agedays, internal_id).
2. **cf_detail columns.** When `cf_detail = TRUE`:
 - `cf_status := "CF-Resc"` for rows with non-empty
 `cf_rescued`; `"CF-NR"` for CF candidates that were not
 rescued; `NA` otherwise.
 - `cf_deltaZ := absdiff` for all CF candidates.
3. **Drop temp columns.** `cols_to_drop_6` drops `v.orig`,
 `wholehalfimp`, `seq_win`, `cs`, `absdiff`,
 `ageday_has_include`, `orig_ageday`, `cf_interval` at the
 end of Step 6.

### Checklist findings (2026-04-17 walk)

1. **Sort determinism.** All sorts include `internal_id` as
 the final tiebreaker: `setkey(data.df, subjid, param,
 agedays, internal_id)` at the top of `cleanchild()`,
 `cf_subset[order(subjid, param, agedays, internal_id)]`
 inside the CF block, and
 `data.df[order(subjid, param, agedays, internal_id)]`
 after the SDE-Identical rbind.
2. **Birth tiebreaking.** Step 6 does not perform
 age-dependent tiebreaking. The birth vs non-birth rule
 lives in Early Step 13 and Step 5.
3. **Z-score correctness.** Detection uses `v.orig` (raw
 measurement, exact equality). The originator z-score
 and `absdiff` use `sd.orig_uncorr` (pre-GA-correction),
 which is the correct reference-based raw z-score.
4. **Dead code.** None found. Pre-2026-04-14 artifacts
 (`cf_string_length`, `Rescued-Imperial`,
 `Rescued-Adol{,-Imperial}`, commented-out dplyr/`map_lgl`
 blocks) are all removed. The previously-dead
 `cf_threshold` store on `data.df` was removed in this
 walk (F17). The `run_cf_detection` optimization that
 conditionally skipped the Step 6 body when `cf_rescue =
 "all"` and `cf_detail = FALSE` was removed (F19);
 detection and rescue now always run, producing consistent
 `cf_rescued` labels regardless of `cf_detail`.
5. **Exclusion code names.** `Exclude-C-CF` appears in
 `exclude.levels.peds`. `Rescued` and `Rescued-All` are
 the only values written to `cf_rescued`.
6. **Column names.** `data.df`, `v.orig`, `sd.orig_uncorr`,
 `tbc.sd`, `exclude` all match current conventions.
7. **Configurable parameter defaults.** `cf_rescue =
 "standard"` and `cf_detail = FALSE` match the defaults
 documented in `gc-github-latest/CLAUDE.md`.
8. **Step linkage.** Prior: Step 5 (Temporary SDEs); next:
 Step 7 (BIV).
9. **Grepl vs. exact matching.** Step 6 uses exact
 equality throughout (`exclude == "Include"`, `exclude
 == "Exclude-C-CF"`, etc.); no `grepl()` on exclusion
 codes.
10. **Inline comment accuracy.** Comments reviewed in this
 walk; no stale references remain after F17 (removed
 the stale "Store thresholds on data.df for cf_detail
 output" comment alongside its dead assignment).
11. **Stale content.** Removed the dead `cf_threshold`
 store and updated this section to reflect the
 2026-04-14 lookup-based rescue design.
12. **Threshold reconciliation.** `.cf_rescue_lookup()`
 matrices match `cf-rescue-thresholds.md`; HEADCM reuses
 HEIGHTCM as a known placeholder (flagged above).
13. **`.child_valid()` flags.** Phase 1 uses
 `.child_valid(data.df, include.temporary.extraneous =
 TRUE)` — correct. Temp SDEs are legitimate CF
 candidates themselves.
14. **DOP logic.** Step 6 does not use the designated
 other parameter.

---

## Step 7: BIV (Biologically Implausible Values)

| | |
|---|---|
| **Scope** | All parameters |
| **Operates on** | Valid rows including temp SDEs |
| **Prior step** | Step 6 (Carried Forwards) |
| **Next step** | Step 9 (Evil Twins) |
| **Exclusion codes** | `Exclude-C-BIV` |
| **Code location** | `R/child_clean.R` lines ~2882–2976 |

### Overview

BIV identifies measurements that are biologically impossible
or extremely unlikely, using two blocks applied in sequence:

1. **Absolute BIV** — fixed limits on raw measurements (`v`).
2. **Standardized BIV** — cutoffs on unrecentered CSD
 z-scores (`sd.orig_uncorr`).

After both blocks, the temp-SDE identification is rerun:
all rows currently flagged `Exclude-C-Temp-Same-Day` are
reset to `Include`, then `identify_temp_sde()` is called
against the post-BIV state.

The working `valid_set` is computed once at the top of the
step via `.child_valid(data.df, include.temporary.extraneous
= TRUE)` and reused across both BIV blocks. It contains
`Include` rows and rows flagged `Exclude-C-Temp-Same-Day`;
any other non-`Include` row is not in `valid_set` and so
cannot be assigned `Exclude-C-BIV` here.

### Absolute BIV thresholds

All conditions use strict `<` or `>` (not `<=` / `>=`). Rows
that match get `exclude := "Exclude-C-BIV"`.

**Weight (WEIGHTKG):**

| Condition | Threshold | Rationale |
|---|---|---|
| First year (`agedays <= 365`) | `v < 0.2` kg | Published minimum viable birth weight; allows preterm |
| After first year (`agedays > 365`) | `v < 1` kg | Minimum plausible outpatient weight |
| Birth (`agedays == 0`) | `v > 10.5` kg | Maximum plausible birth weight |
| `ageyears < 2` | `v > 35` kg | Highest plausible weight under 2y |
| Any age | `v > 600` kg | Published maximum human weight |

**Height (HEIGHTCM):**

| Condition | Threshold | Rationale |
|---|---|---|
| Any age | `v < 18` cm | Fenton z=-6 at 22 0/7 weeks |
| Any age | `v > 244` cm | Published maximum |
| Birth (`agedays == 0`) | `v > 65` cm | Fenton z=6 at 40 0/7 weeks |

**Head circumference (HEADCM):**

| Condition | Threshold | Rationale |
|---|---|---|
| Any age | `v < 13` cm | Fenton z=-6 at 22 0/7 weeks |
| Any age | `v > 75` cm | Published maximum |
| Birth (`agedays == 0`) | `v > 50` cm | Birth maximum |

### Standardized BIV thresholds

Uses **unrecentered** z-scores (`sd.orig_uncorr`). Each
cutoff is a user-settable `cleangrowth()` parameter; the
default values (shown in the tables below) preserve the
prior hardcoded behavior. The standardized block adds a
`!grepl(biv_pattern, exclude)` guard where
`biv_pattern <- "^Exclude-C-BIV$"`. The guard skips rows
that the absolute block just assigned `Exclude-C-BIV`,
because `valid_set` was computed before the absolute block
ran and therefore still treats those rows as valid. Rows
originally flagged `Exclude-C-Temp-Same-Day` that were not
hit by absolute BIV can be overwritten to `Exclude-C-BIV`
by the standardized block; the 7d temp-SDE rerun afterward
re-flags any SPA that still has a duplicate Include.

**Weight:**

| Condition | Threshold | Parameter (default) |
|---|---|---|
| `ageyears < 1` | `sd.orig_uncorr < biv.z.wt.low.young` | `biv.z.wt.low.young` (-25) |
| `ageyears >= 1` | `sd.orig_uncorr < biv.z.wt.low.old` | `biv.z.wt.low.old` (-15) |
| Any age | `sd.orig_uncorr > biv.z.wt.high` | `biv.z.wt.high` (22) |

**Height:**

| Condition | Threshold | Parameter (default) |
|---|---|---|
| `ageyears < 1` | `sd.orig_uncorr < biv.z.ht.low.young` | `biv.z.ht.low.young` (-25) |
| `ageyears >= 1` | `sd.orig_uncorr < biv.z.ht.low.old` | `biv.z.ht.low.old` (-15) |
| Any age | `sd.orig_uncorr > biv.z.ht.high` | `biv.z.ht.high` (8) |

The default upper HT cutoff (`biv.z.ht.high = 8`) is
tighter than the default upper WT cutoff
(`biv.z.wt.high = 22`). The inline comment attributes the
tighter HT default to analysis of CHOP data showing the
`±15`/`±25` range was too loose for heights.

**Head circumference:**

| Condition | Threshold | Parameter (default) |
|---|---|---|
| Any age | `sd.orig_uncorr < biv.z.hc.low` | `biv.z.hc.low` (-15) |
| Any age | `sd.orig_uncorr > biv.z.hc.high` | `biv.z.hc.high` (15) |

### Temp SDE re-evaluation (7d)

After both BIV blocks, all rows currently flagged
`Exclude-C-Temp-Same-Day` are reset to `Include`, then
`identify_temp_sde()` is rerun on the full dataset. This
mirrors the Step-5 temp-SDE logic against the post-BIV
state: if the prior temp-SDE keeper on an SPA has just
been excluded as BIV, another value in that SPA becomes
the flagged duplicate.

### Configurable parameters in scope for Step 7

Eight user-settable standardized-BIV cutoffs:
`biv.z.wt.low.young` (-25), `biv.z.wt.low.old` (-15),
`biv.z.wt.high` (22), `biv.z.ht.low.young` (-25),
`biv.z.ht.low.old` (-15), `biv.z.ht.high` (8),
`biv.z.hc.low` (-15), `biv.z.hc.high` (15). Defaults
reproduce the pre-configurable hardcoded thresholds.
Absolute-BIV thresholds (`v < 0.2` kg, `v > 244` cm, etc.)
are not user-configurable.

### Variables created and dropped

- `ageyears`: created at the top of Step 7 for age-
 threshold checks, dropped at the end of the step.
- `biv_pattern`: local character scalar
 (`"^Exclude-C-BIV$"`) used only as the
 standardized-block guard; not stored on `data.df`.

### Checklist findings

1. **Operates-on set.** `.child_valid(data.df,
 include.temporary.extraneous = TRUE)` is the correct
 selector for BIV: Include plus temp-SDE rows.
2. **Boundaries — all strict inequalities.** `< 0.2`,
 `< 1`, `> 10.5`, `> 244`, `< -25`, `> 22`, etc. A
 weight of exactly 0.2 kg is not excluded; a
 standardized z of exactly `-25` is not excluded.
3. **Age boundary at 1 year.** Standardized BIV uses
 `ageyears < 1` vs. `ageyears >= 1`. A child whose
 `ageyears` equals 1 exactly falls in the `>= 1` band.
4. **Valid-set refresh.** `valid_set` is computed once
 and not refreshed between absolute and standardized
 BIV. The standardized block relies on the
 `!grepl(biv_pattern, exclude)` guard to skip rows just
 assigned `Exclude-C-BIV` by the absolute block. Non-
 `Include` rows other than `Exclude-C-Temp-Same-Day`
 are not in `valid_set` and so cannot be overwritten by
 Step 7.
5. **Parameter scope.** All three params (WEIGHTKG,
 HEIGHTCM, HEADCM) have their own absolute and
 standardized thresholds.
6. **Factor levels.** `Exclude-C-BIV` is present in
 `exclude.levels`; absolute and standardized BIV share
 the single code.
7. **Standardized-BIV cutoffs are parameterized.** Eight
 `biv.z.*` cutoffs flow from `cleangrowth()` through
 `cleanchild()` into Step 7. Defaults reproduce prior
 hardcoded thresholds. See the "Configurable parameters
 in scope for Step 7" subsection above.
8. **Efficiency (minor).** `valid_set` could be refreshed
 after the absolute block to skip newly-excluded rows
 in the standardized block. Current design relies on
 the `!grepl(biv_pattern, exclude)` guard; impact is
 limited to very large datasets.

---

## Step 9: Evil Twins

| | |
|---|---|
| **Scope** | All parameters |
| **Operates on** | Valid rows, excluding temp SDEs, requiring 3+ total rows per subject-param |
| **Prior step** | Step 7 (BIV) |
| **Next step** | Step 11 (EWMA1 — Extreme EWMA) |
| **Exclusion code** | `Exclude-C-Evil-Twins` |
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
 d. Exclude that one value (`Exclude-C-Evil-Twins`)
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
5. **`.child_valid()` call:** Correctly excludes temp SDEs
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
9. **Factor levels:** `Exclude-C-Evil-Twins` exists in
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
| **Exclusion code** | `Exclude-C-Traj-Extreme` |
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
5. **`.child_valid()` call:** Correctly excludes temp SDEs.
 Only Include values participate in EWMA calculation
 and exclusion candidates.
6. **ctbc optimization:** Skips ctbc EWMA computation
 when `ctbc.sd == tbc.sd` (copies instead). This is
 correct for non-potcorr subjects.
7. **Parameter scope:** All 3 params handled uniformly.
8. **Factor levels:** `Exclude-C-Traj-Extreme` exists in
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
| **Exclusion codes** | `Exclude-C-Traj` (all moderate EWMA sub-rules collapsed into one code) |
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
3. **`.child_valid()` calls:** Step 15 excludes temp SDEs AND
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
| **Exclusion codes** | `Exclude-C-Abs-Diff` |
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
6. **`.child_valid()` call:** Excludes temp SDEs, excludes weight.
7. **Sort order:** Explicit `order(agedays, id)` inside
 each while-loop iteration.
8. **Parameter scope:** HT and HC only. HC has separate
 velocity tables, tighter tolerances, and no Tanner.
9. **Factor levels:** `Exclude-C-Abs-Diff` exists in
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
| **Exclusion codes** | `Exclude-C-Pair`, `Exclude-C-Single` |
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
 gap >= 365.25` → `Exclude-C-Pair`
 - `|diff_tbc| > 2.5 & (|diff_ctbc| > 2.5 or NA) &
 gap < 365.25` → `Exclude-C-Pair`
5. If one excluded, re-evaluate remaining as single

### Single evaluation

For 1 remaining measurement:
- `(|tbc.sd| > 3 & comp_diff > 5)` OR
 `(|tbc.sd| > 5 & no DOP data)`
 → `Exclude-C-Single`

### Checklist findings

1. **Boundaries:** Pair: `> 4` / `> 2.5` strict; gap
 `>= 365.25` inclusive. Single: `> 3` / `> 5` strict.
2. **`.child_valid()` call:** Excludes temp SDEs. Counts only
 Include rows for singles/pairs determination.
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
| **Exclusion code** | `Exclude-C-Too-Many-Errors` |
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
 `Exclude-C-Too-Many-Errors`

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
5. **Factor level:** `Exclude-C-Too-Many-Errors` exists.
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
