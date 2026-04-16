# Adult Growthcleanr Integration Guide

**Last edited:** 2026-04-03

Reference for integrating adult growthcleanr into the combined growthcleanr R package.

---

## 1. File Replacement

Replace the package's adult files with our versions:

| Package file | Replace with | Notes |
|---|---|---|
| `R/adult_clean.R` | `Carrie-R-Code/adult_clean.R` | Main `cleanadult()` function |
| `R/adult_support.R` | `Carrie-R-Code/02a_support.R` | Rename to `adult_support.R` on copy |

The package's old `adult_support.R` has Wen's original versions of `check_between`, `ewma_dn`, `remove_biv`, etc. Our `02a_support.R` has rewritten versions of all these plus new functions (permissiveness presets, compute_et_limit, compute_perc_limit, etc.).

**Also copy:**

| File | Destination | Notes |
|---|---|---|
| `Carrie-R-Code/test_cleanadult.R` | `tests/testthat/test-adult-clean.R` | Remove `source()` calls at top |
| `adult-gc-test-ALL-PHASES.csv` | `tests/testthat/` or `inst/testdata/` | Regression test dataset |
| `Carrie-R-Code/test_harness.R` | `tests/` | Update paths after move |

---

## 2. cleanadult() Interface

```r
cleanadult(df,
           permissiveness = "looser",        # "loosest"/"looser"/"tighter"/"tightest"
           # BIV limits (NULL = use preset)
           overall_ht_min, overall_ht_max, overall_wt_min, overall_wt_max,
           overall_bmi_min, overall_bmi_max,
           # 1D limits (NULL = use preset)
           single_ht_min_bmi, single_ht_max_bmi, single_wt_min_bmi, single_wt_max_bmi,
           single_ht_min_nobmi, single_ht_max_nobmi, single_wt_min_nobmi, single_wt_max_nobmi,
           single_bmi_min, single_bmi_max,
           # Algorithm parameters (NULL = use preset)
           repval_handling, ht_band, allow_ht_loss, allow_ht_gain,
           wtallow_formula, error_load_threshold,
           ewma_cap_short, ewma_cap_long,
           wt_scale_et, wt_scale_wtallow, wt_scale_threshold,
           mod_ewma_f,
           # Non-permissiveness parameters
           scale_max_lbs = Inf, ewma_window = 15,
           # Output options
           include_bin_result = TRUE, include_extraneous = FALSE,
           include_ht_groups = FALSE,
           quietly = FALSE)
```

**Input:** data.table with columns `id`, `subjid`, `param`, `agedays`, `sex`, `measurement`. Accepts HEIGHTCM/HEIGHTIN and WEIGHTKG/WEIGHTLBS (converts internally).

**Output:** Input df plus `result` (exclusion code), `mean_ht`, `bin_result` (default). Internal columns (`meas_m`, `ageyears`, `age_days`) dropped. Original `id` type preserved.

---

## 3. Add Adult Parameters to cleangrowth()

The existing `cleangrowth()` has only `weight_cap = Inf` for adult control. With the permissiveness framework, most users only need to pass `permissiveness`:

```r
# Add to cleangrowth() signature:
adult_permissiveness = "looser",    # Sets defaults for all adult sub-parameters
adult_scale_max_lbs = Inf,          # Physical scale upper limit in lbs
```

Optionally expose sub-parameters for advanced users who want to override specific settings while keeping a permissiveness preset for everything else. All sub-parameters default to NULL (use preset).

**Rename `weight_cap`:** The existing `weight_cap` parameter should be renamed to `adult_scale_max_lbs` (or `scale_max_lbs`). Add a deprecation shim:

```r
if (!missing(weight_cap)) {
  warning("weight_cap has been renamed to adult_scale_max_lbs.", call. = FALSE)
  adult_scale_max_lbs <- weight_cap
}
```

---

## 4. Parameter Passthrough

Update the adult dispatch section in `cleangrowth()` (currently around line 1491):

```r
# Current:
res <- cleanadult(data.adult, weight_cap = weight_cap)

# New:
res <- cleanadult(data.adult,
                  permissiveness = adult_permissiveness,
                  scale_max_lbs = adult_scale_max_lbs,
                  quietly = quietly)
```

Also update the parallel `ddply` call similarly.

---

## 5. Update var_for_par (Parallel Worker Exports)

The `var_for_par` list (around line 1457) names every function needed by parallel workers. Replace with:

```r
var_for_par <- c(
  # Main function
  "cleanadult",
  # Permissiveness
  "permissiveness_presets", "resolve_permissiveness",
  # Convenience
  "check_between", "round_pt",
  # Dynamic threshold helpers
  "compute_et_limit", "compute_perc_limit", "compute_wtallow",
  # EWMA
  "as.matrix.delta_dn", "ewma_dn",
  "ewma_cache_init", "ewma_cache_update",
  # BIV / RV / SDE
  "remove_biv", "remove_biv_high", "remove_biv_low",
  "identify_rv", "temp_sde", "redo_identify_rv",
  # Height
  "ht_allow", "ht_change_groups", "ht_3d_growth_compare",
  # Moderate EWMA support
  "detect_runs", "compute_trajectory_fails",
  # Evil Twins
  "evil_twins",
  # RV propagation (linked mode)
  "propagate_to_rv",
  # Extreme EWMA
  "remove_ewma_wt",
  # Moderate EWMA
  "remove_mod_ewma_wt",
  # 2D Non-Ordered
  "eval_2d_nonord",
  # 1D
  "eval_1d",
  # Error Load
  "eval_error_load"
)
```

---

## 6. Column Mapping (Already Handled by cleangrowth)

`cleangrowth()` sets up these columns before calling `cleanadult()`:
- `age_years := agedays/365.25` — our code renames to `ageyears`
- `measurement := v_adult` — matches our expectation
- `id := line` — our code uses `id` as primary key; `line` passes through untouched
- `subjid`, `param`, `agedays`, `sex` — present and correct

The join code after `cleanadult()` returns uses `res$line` and `res$result`. Our code preserves `line` (doesn't modify it) and outputs `result`. Compatible.

---

## 7. Removed Functions (vs. Old Package Code)

Functions in the old `adult_support.R` that are NOT in our version (they were pediatric-specific):
- `get_float_rem`, `rem_hundreds`, `rem_unit_errors`
- `get_num_places`, `switch_tens_ones`, `rem_transpositions`

These are only referenced from the old `adult_clean.R` which is also being replaced.

---

## 8. NAMESPACE / roxygen

No new exports needed — `cleanadult()` and all support functions are `@noRd` (internal). The only export is `cleangrowth()`.

Dependencies (all already imported by package):
- `data.table`
- `stats::median`
- No new package dependencies

---

## 9. Algorithm Steps

1H/1W BIV → 2W RV markers → 3H/3W Temp SDE → 4W Weight Cap → 9Wa Evil Twins → 9H HT SDE → 9Wb Extreme EWMA → 10H HT Distinct → 10W WT SDE → 11H Mean HT → 11Wa 2D Ord WT → 11Wa2 2D Non-Ord WT → 11Wb Moderate EWMA → 13 Distinct 1D → 14 Error Load

**There is no Step 12W.**

---

## 10. Exclusion Codes

### Non-SDE Codes

| Code | Step | Description |
|------|------|-------------|
| `Include` | — | Not excluded |
| `Exclude-A-BIV` | 1H/1W | Biologically implausible value (height, weight, or BMI) |
| `Exclude-A-Scale-Max` | 4W | Weight at scale maximum |
| `Exclude-A-Scale-Max-Identical` | 4W | All weights identical at scale max |
| `Exclude-A-Scale-Max-RV-Propagated` | 4W | RV copy of scale-max exclusion (linked mode) |
| `Exclude-A-Evil-Twins` | 9Wa | Adjacent pair with implausible weight difference |
| `Exclude-A-Traj-Ext` | 9Wb | Extreme EWMA outlier (independent mode) |
| `Exclude-A-Traj-Extreme-firstRV` | 9Wb | Extreme EWMA outlier (linked firstRV pass) |
| `Exclude-A-Traj-Extreme-allRV` | 9Wb | Extreme EWMA outlier (linked allRV pass) |
| `Exclude-A-Ord-Pair` | 10Ha | 2D height pair outside band (one excluded) |
| `Exclude-A-Ord-Pair-All` | 10Ha | 2D height pair outside band (all excluded) |
| `Exclude-A-Window` | 10Hb | 3+D height outside window (one excluded) |
| `Exclude-A-Window-All` | 10Hb | 3+D height outside window (all excluded) |
| `Exclude-A-2D-Ordered` | 11Wa | 2D ordered weight pair outside wtallow/perclimit |
| `Exclude-A-2D-Non-Ordered` | 11Wa2 | 2D non-ordered weight pair |
| `Exclude-A-Traj-Moderate` | 11Wb | Moderate EWMA outlier (independent or firstRV) |
| `Exclude-A-Traj-Moderate-allRV` | 11Wb | Moderate EWMA outlier (linked allRV pass) |
| `Exclude-A-Traj-Moderate-Error-Load` | 11Wb | 4+ consecutive moderate EWMA candidates |
| `Exclude-A-Traj-Moderate-Error-Load-RV` | 11Wb | Error load escalation to entire patient (linked) |
| `Exclude-A-Single` | 13 | 1D value outside limits (height, weight, or BMI) |
| `Exclude-A-Too-Many-Errors` | 14 | Error ratio > threshold |

### Same-Day Event (SDE) Codes

| Code | Steps | Description |
|------|-------|-------------|
| `Exclude-A-Identical` | 9H/10W | Same-day identical values (keep one) |
| `Exclude-A-Extraneous` | 9H/10W | Same-day non-identical value (SDE loser) |

Exclusion codes are no longer param-specific — the same code is used regardless of whether the excluded value is height or weight.

---

## 11. Permissiveness Summary

Default: `"looser"`

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| BIV HT (cm) | 50–244 | 120–230 | 142–213 | 147–208 |
| BIV WT (kg) | 20–500 | 30–270 | 36–159 | 39–136 |
| BIV BMI | 5–300 | 12–65 | 16–45 | 18–40 |
| 1D limits | split (BMI/no-BMI) | same as BIV | same as BIV | same as BIV |
| wtallow formula | piecewise | piecewise | piecewise-lower | allofus15 |
| wtallow scaling | +0.50×(maxwt−120)⁺ | +0.50×(maxwt−120)⁺ | none | none |
| EWMA cap <6m | 50+0.70×scale | 50+0.70×scale | 40 | 40 |
| EWMA cap ≥6m | 80+0.70×scale | 80+0.70×scale | 60 | 40 |
| perclimit ≤45 kg | 0.5 | 0.5 | 0.7 | 0.7 |
| perclimit 45–80 kg | 0.4 | 0.4 | 0.4 | 0.4 |
| perclimit >80 kg | 0 (disabled) | 0 (disabled) | 0.4 | 0.4 |
| error_load_threshold | 0.41 | 0.41 | 0.29 | 0.29 |
| mod_ewma_f | 0.75 | 0.75 | 0.60 | 0.60 |
| ht_band | 3" | 3" | 2" | 2" |
| allow_ht_loss | TRUE | FALSE | FALSE | FALSE |
| allow_ht_gain | TRUE | TRUE | TRUE | FALSE |
| repval_handling | independent | independent | linked | linked |

---

## 12. Testing

```bash
# Regression test at any permissiveness level (1483/1483 at all levels)
Rscript test_harness.R [loosest|looser|tighter|tightest]

# Unit tests (210/210)
Rscript -e 'testthat::test_file("test_cleanadult.R")'
```

**Test CSV columns:** `id`, `subjid`, `param`, `agedays`, `sex`, `measurement`, plus legacy columns, plus `expected_loosest`, `expected_looser`, `expected_tighter`, `expected_tightest`, `test_description`, `test_category`, `precision_sensitive`.

**Test file changes for package:**
- Remove `source()` calls (lines 8-9 of test_cleanadult.R) — package loading handles this
- `run_clean()` helper calls `cleanadult()` directly, which is fine for unit testing internals
- Add integration tests through `cleangrowth()` after merge

---

## 13. Quick Compatibility Test

After replacing files and updating the passthrough:

```r
library(growthcleanr)
result <- cleangrowth(
  subjid = c("A", "A", "A"),
  param = c("WEIGHTKG", "WEIGHTKG", "WEIGHTKG"),
  agedays = c(10000, 10365, 10730),
  sex = c(0, 0, 0),
  measurement = c(75, 76, 74),
  id = c("1", "2", "3"),
  quietly = TRUE
)
# All three should be "Include"
stopifnot(all(result$exclude == "Include"))
```

---

## 14. Design Notes

- **Rounding tolerance:** 0.12 cm/kg on all threshold comparisons. See algorithm narrative Rounding Tolerance section.
- **perclimit scope:** 11Wa uses subject-level (max across obs). 11Wb uses observation-level.
- **sex variable:** Not used by adult algorithm but required by package. Arbitrary value OK.
- **Sort determinism:** All sorts include `internal_id` as final tiebreaker.
- **Missing-as-infinity:** R uses `ifelse(is.na(...), Inf, ...)` for edge EWMA values.
- **ID handling:** `internal_id` (sequential integer from `cleangrowth()`) is used for all internal processing. Converted to character inside `cleanadult()` for named vector indexing. User's original `id` is preserved and restored on output.

---

## 15. Remaining Work

- Package integration
- Deferred test gaps: error load with -5 exponent, weight scaling at perm levels
- Performance: `setkey(df, subjid)` optimization (deferred)
