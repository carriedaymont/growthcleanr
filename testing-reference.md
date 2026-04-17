# growthcleanr Testing Reference

**Last updated:** 2026-04-14

## Summary

| File | Tests | Assertions | Scope |
|------|-------|------------|-------|
| `test-cleangrowth.R` | 14 | 92 | `cleangrowth()` API: legacy pediatric, adult integration, child-adult spanning, deprecation |
| `test-child-regression.R` | 9 | 54 | Frozen counts, spot checks, cross-sample stability, Missing, HC, preload_refs, changed_subjids |
| `test-child-algorithms.R` | 24 | ~40 | CF rescue, Evil Twins/OTL, Error Load, HC boundary, cf_detail, parallel, GA correction, Birth EWMA2, LENGTHCM |
| `test-child-parameters.R` | 9 | 21 | Parameter acceptance and behavioral changes |
| `test-child-edge-cases.R` | 12 | 24 | Sparse data, missing data, extreme values, determinism |
| `test-adult-clean.R` | 198 | 198 | All 14 adult steps, 4 permissiveness levels, edge cases (via `cleanadult()` directly) |
| **Total** | **266** | **~429** | |

Additional test files not actively maintained:
- `test-utils.R`: 10 failures (old API without `id` parameter; tests `splitinput`, `recode_sex`, `longwide`, `simple_bmi`)
- `test-cdc.R`: not modified in v3.0.0

---

## Running Tests

```bash
# Single file
NOT_CRAN=true Rscript -e 'testthat::test_file("tests/testthat/test-child-algorithms.R")'

# All files
NOT_CRAN=true Rscript -e 'testthat::test_dir("tests/testthat/")'

# Adult regression harness (separate from testthat)
Rscript tests/test_harness.R [loosest|looser|tighter|tightest]
```

**IMPORTANT:** Tests use `library(growthcleanr)` which loads the *installed* package. After any code changes, reinstall before testing:
```r
devtools::install_local("gc-github-latest", force = TRUE, upgrade = "never")
```

---

## test-cleangrowth.R

Tests for the `cleangrowth()` entry point — the public API that dispatches to child and adult algorithms and assembles combined output.

### Legacy pediatric

| Test | What it verifies |
|------|------------------|
| legacy algorithm works on pediatric synthetic data | NHANES vs derived recentering, spot-check specific records, frozen exclusion counts, Missing handling |

### Adult (syngrowth-based)

| Test | What it verifies |
|------|------------------|
| works as expected on adult synthetic data | Spot-check records, cutpoint=20 vs 18 comparison, frozen counts, Missing handling |
| works without either adult or pediatric data | Adult-only, pediatric-only, and empty data all return correct row counts |

### Child-adult spanning subjects

| Test | What it verifies |
|------|------------------|
| all rows preserved through child/adult split | Real syngrowth subject (30 child + 23 adult rows): no rows lost, no NA codes |
| child rows get child codes, adult rows get adult codes | Code prefix correctness: `Exclude-C-*` for child, `Exclude-A-*` for adult |
| constructed subject processed correctly across boundary | Synthetic subject 18y–22y: all rows returned, correct code prefixes, normal values mostly included |

### Adult integration (cleangrowth → cleanadult)

| Test | What it verifies |
|------|------------------|
| output columns present and correctly typed | `mean_ht` (numeric, populated), `bin_result`, `bin_exclude` present; `exclude` is factor, no NAs |
| imperial units processed correctly | HEIGHTIN/WEIGHTLBS work through full pipeline; param preserved as imperial in output; `mean_ht` in cm |
| adult_permissiveness passes through | 115cm HT → BIV at `looser`, passes BIV at `loosest` |
| adult_scale_max_lbs passes through | 136kg flagged with cap=300lbs, included without cap |
| weight_cap deprecation warns and works | Deprecated parameter triggers warning, still functions |
| character/UUID ids preserved | String ids survive the round trip intact |
| spanning subject output structure | Child rows have NA `mean_ht`; adult rows have populated `mean_ht`; `bin_exclude` non-NA for all |

### Deprecation

| Test | What it verifies |
|------|------------------|
| prelim_infants = TRUE runs child algorithm with deprecation warning | Warning message, correct results, valid factor output |

---

## test-child-regression.R

Frozen-count and stability tests using syngrowth data. These break when algorithm behavior changes intentionally — update expected values when that happens.

| Test | What it verifies |
|------|------------------|
| structural invariants on 100-subject subset | Row counts, column presence, factor type, no NAs, at least 1 Include |
| exclusion counts match expected on 100 subjects | Frozen counts by exclusion code (regression baseline) |
| spot-check individual records | Specific record ids produce specific exclusion codes |
| within-subject exclusions stable at 50 vs 100 subjects | Same subject produces same code regardless of batch context |
| missing measurements return Missing | NA measurements → `Exclude-Missing` |
| HEADCM measurements are processed | HC rows present in output with valid codes |
| prelim_infants = TRUE produces deprecation warning | Backward compat: warns and runs child algorithm |
| gc_preload_refs: identical results to standard run | `ref_tables` parameter produces byte-identical output |
| changed_subjids: partial run produces identical results | Cached + partial rerun matches full run |

---

## test-child-algorithms.R

Tests that construct specific scenarios and verify the algorithm responds correctly. Organized by algorithm step.

### Section 1: CF rescue (Step 6)

| Test | What it verifies |
|------|------------------|
| standard rescues more CFs than none | `cf_rescue="standard"` leaves fewer CF exclusions than `"none"` |
| cf_rescued column populated for rescued CFs | Rescued rows have `cf_rescued` codes, don't retain CF exclusion |
| none mode rescues no CFs | `cf_rescue="none"` → zero rescued rows |
| all mode rescues all CFs | `cf_rescue="all"` → zero CF exclusions |
| cf_detail produces cf_status and cf_deltaZ | `cf_detail=TRUE` adds columns with valid values |
| cf_detail columns absent by default | Columns not present when `cf_detail=FALSE` |

### Section 2: CF rescue threshold-specific (Step 6)

| Test | What it verifies |
|------|------------------|
| infant short-interval CF with deltaZ=0.35 is rescued | 0-3mo, <1wk, HT-Other → threshold 0.40; deltaZ 0.35 rescued |
| NR cell (10-15y, >1y interval, WT) excludes CF | NR cell (threshold 0.00) excludes even zero-deltaZ CF |
| NR cell (6-12mo, >1mo interval) excludes all CFs | Second NR cell verification |

### Section 3: Evil Twins / OTL (Step 9)

| Test | What it verifies |
|------|------------------|
| single unit error (height x2.54) is excluded | Injected unit error on syngrowth subject is caught |
| two consecutive unit errors are both excluded | Both injected errors caught |
| three consecutive unit errors are all excluded | All three caught |
| unit errors don't exclude neighboring clean measurements | Collateral damage ≤ 2 rows |

### Section 4: Error Load (Step 21)

| Test | What it verifies |
|------|------------------|
| threshold parameter controls exclusion behavior | `error.load.threshold=0.4` catches more than 0.5 |
| mincount parameter controls exclusion behavior | `error.load.mincount=1` catches more than 10 |
| subject with many errors triggers Error-load | Corrupted subject → Too-Many-Errors on remaining includes |

### Section 5: Age boundaries

| Test | What it verifies |
|------|------------------|
| HEADCM cleaned under 3yr, not cleaned over 3yr | HC at 1095 days → cleaned; at 1097 → `Exclude-Not-Cleaned` |

### Section 6: Parallel execution

| Test | What it verifies |
|------|------------------|
| parallel=TRUE matches parallel=FALSE | Identical exclusion codes on 100 subjects |

### Section 7: GA correction / potcorr (Step 2b)

| Test | What it verifies |
|------|------------------|
| very low birth weight gets potcorr correction | 1.5kg birth → `ctbc.sd ≠ tbc.sd`; normal birth → `ctbc.sd == tbc.sd` |
| birth weight z >= -2 does not trigger potcorr | 2.6kg birth → `ctbc.sd == tbc.sd` (no correction) |

### Section 8: Birth EWMA2 (Steps 15-16)

| Test | What it verifies |
|------|------------------|
| extreme birth weight is excluded | 8kg at birth → excluded |
| normal birth weight is included | 3.5kg at birth → included |
| extreme birth height is excluded | 70cm at birth → excluded |

### Section 9: LENGTHCM

| Test | What it verifies |
|------|------------------|
| LENGTHCM identical to HEIGHTCM | Same measurements produce identical exclusion codes |

---

## test-child-parameters.R

Tests that each configurable parameter is accepted and behaves as documented.

| Test | What it verifies |
|------|------------------|
| use_legacy_algorithm = TRUE | Runs legacy code, produces different results than child |
| sd.recenter parameter | Accepted, runs without error |
| include.carryforward | TRUE keeps CFs, FALSE excludes them |
| biv.z.wt.high parameter | Changes standardized-BIV (Step 7) exclusion behavior |
| ewma_window parameter | Accepted, doesn't crash |
| error load parameters | Affect error-load exclusions |
| recover.unit.error parameter | Accepted, doesn't crash |
| imperial units (HEIGHTIN, WEIGHTLBS) | Converted and processed correctly |
| LENGTHCM param | Accepted and processed |

---

## test-child-edge-cases.R

Tests with minimal, degenerate, or boundary-condition data.

| Test | What it verifies |
|------|------------------|
| single subject | Algorithm runs on 1 subject |
| single measurement per param | 1 HT + 1 WT → valid output |
| exactly 2 measurements per param | 2 per param → valid output |
| all-NA measurements | All NA → all `Exclude-Missing` |
| mix of NA and valid | NAs get Missing, valid values get real codes |
| same-day identical measurements | Detected and marked |
| negative agedays | → `Exclude-Missing` |
| HEADCM > 3 years | → `Exclude-Not-Cleaned` |
| biologically implausible values | Extreme values excluded |
| mix of data-rich and data-sparse subjects | Both handled correctly |
| carried-forward values | Detected by algorithm |
| deterministic results | Same input → same output (seed-independent) |

---

## test-adult-clean.R

198 unit tests calling `cleanadult()` directly (via `growthcleanr:::`). Organized by algorithm step.

| Section | Tests | What it covers |
|---------|-------|----------------|
| BIV-only datasets | ~15 | HT/WT/BMI BIV boundaries at all 4 permissiveness levels |
| Minimal / edge cases | ~10 | Single obs, HT-only, WT-only, multi-subject, mixed |
| Weight Cap (4W) | ~5 | Scale max detection, Inf default, identical values |
| Same-day (9H, 10W) | ~10 | Identical HT/WT, extraneous resolution |
| Evil Twins (9Wa) | ~8 | Implausible WT differences, dynamic OOB pairs |
| Extreme EWMA (9Wb) | ~5 | Outlier detection with wide spacing |
| Height Distinct (10H) | ~15 | 2D pairs, 3+D windows, ht_band parameter |
| Moderate EWMA (11Wb) | ~5 | Moderate outliers among normal values |
| 1D Evaluation (13) | ~10 | Single HT/WT/BMI limits, custom limits |
| Error Load (14) | ~5 | Threshold-based exclusion |
| Parameter Validation | ~5 | Invalid formulas, missing columns |
| Output Structure | ~10 | Column presence, optional columns, row preservation |
| Repeated Values (2W) | ~5 | Repeated weight handling |
| 2D Ordered WT (11Wa) | ~10 | Perclimit, ratio thresholds |
| 2D Non-Ordered WT (11Wa2) | ~5 | Non-ordered pair logic |
| Weight Cap Influence | ~5 | Remaining value marking at scale max |
| Full Regression | 1 | 1220-row dataset at all 4 permissiveness levels |
| Imperial Units | ~5 | HEIGHTIN/WEIGHTLBS conversion |
| Mean Height | ~5 | Subject mean HT calculation |
| Multi-Subject Isolation | ~5 | Exclusions don't cross subjects |
| Permissiveness Framework | ~50 | All 4 levels, overrides, specific parameter effects |

### Adult regression harness (separate)

`tests/test_harness.R` runs `cleanadult()` against `inst/testdata/adult-gc-test-ALL-PHASES.csv` (1508 rows) at all 4 permissiveness levels. Not a testthat test — run separately:

```bash
Rscript tests/test_harness.R [loosest|looser|tighter|tightest]
```

---

## Coverage Gaps

### Not tested (by design)

- **Steps 3-4 (unit error recovery / swapped measurements):** Legacy features scheduled for redesign. Tests would need rewriting after redesign.

### Known gaps

- **test-utils.R:** 10 failures — old API without `id` parameter. Tests `splitinput`, `recode_sex`, `longwide`, `simple_bmi`. Need updating for v3.0.0 data.table return format.
- **Step 15 (Moderate EWMA) child:** Only birth WT path tested (Steps 15/16). Middle, first, last EWMA2 sub-steps not directly tested with constructed data (covered indirectly via regression counts).
- **Step 17 (Height/HC velocity):** Not directly tested with constructed data (covered via regression counts).
- **Step 19 (Pairs and singles):** Not directly tested with constructed data (covered via regression counts).
- **Step 13 (Final SDE resolution):** Only same-day identical tested; EWMA-based SDE resolution not directly tested.
- **Adult error load with -5 exponent:** Deferred.
- **Adult UW scaling edge cases:** Very low/high upper weight not tested.
- **Adult integration tests through `cleangrowth()`:** Now added (7 tests). Previously only direct `cleanadult()` calls were tested.
