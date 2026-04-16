testthat::skip_on_cran()
library(growthcleanr)
library(data.table)

# =============================================================================
# Stress test: 400 subjects (8 archetypes × 50) with layered errors
#
# Uses a pre-generated fixture (stress_test_data.rds) containing:
# - 8 growth trajectory archetypes (tracker, falter, catchup, sgaAsym,
#   sgaSym, preterm, latePreterm, PGR)
# - 9 error types independently applied to 10% of patients each
#   (CFs, CF chains, SDE similar, SDE extreme, ±10%, ±50%,
#    unit errors, outlier spikes, swapped HT↔WT)
# - ~33K rows total (31K clean + ~1.3K added SDE rows)
#
# Regenerate fixture with:
#   cd "__Pipeline/error-impact"
#   Rscript "../gc-github-latest/scripts/generate_stress_test_data.R"
#
# Primary purposes:
# 1. Catch algorithm changes via frozen exclusion counts
# 2. Performance benchmark on a realistic, error-laden dataset
# =============================================================================

# ---------------------------------------------------------------------------
# Load fixture
# ---------------------------------------------------------------------------
fixture_path <- file.path(
  testthat::test_path(), "stress_test_data.rds"
)

skip_if(
  !file.exists(fixture_path),
  "Stress test fixture not found. Run scripts/generate_stress_test_data.R from error-impact/ first."
)

fixture <- readRDS(fixture_path)
stress_data <- fixture$data

# ---------------------------------------------------------------------------
# Test 1: Structural invariants
# ---------------------------------------------------------------------------
test_that("stress test: structural invariants hold on 400-subject errored dataset", {

  t0 <- proc.time()
  res <- cleangrowth(
    subjid      = stress_data$subjid,
    param       = stress_data$param,
    agedays     = stress_data$agedays,
    sex         = stress_data$sex,
    measurement = stress_data$measurement,
    id          = stress_data$id,
    quietly     = TRUE
  )
  elapsed <- (proc.time() - t0)[["elapsed"]]

  # Report timing
  cat(sprintf(
    "\n=== STRESS BENCHMARK: %d subjects (%d rows) completed in %.1f sec ===\n",
    fixture$n_subjects, nrow(res), elapsed
  ), file = stderr())

  # 1a. Output rows match input
  expect_equal(nrow(res), nrow(stress_data),
               info = "Output row count must match input")

  # 1b. All input IDs present
  expect_true(all(stress_data$id %in% res$id),
              info = "Every input id must appear in output")

  # 1c. No NA exclusion codes
  expect_false(any(is.na(res$exclude)),
               info = "No exclusion codes should be NA")

  # 1d. Valid factor levels
  expect_true(is.factor(res$exclude))

  # 1e. Include + exclusions = total
  counts <- res[, .N, by = exclude]
  expect_equal(sum(counts$N), nrow(stress_data))

  # 1f. Generous performance ceiling (should be ~20-40 sec normally)
  expect_lt(elapsed, 180,
            label = sprintf("Runtime %.1f sec exceeded 180 sec ceiling",
                            elapsed))

  # Save result for subsequent tests
  assign("stress_result", res, envir = parent.env(environment()))
  assign("stress_elapsed", elapsed, envir = parent.env(environment()))
})

# ---------------------------------------------------------------------------
# Test 2: Frozen exclusion counts — update when algorithm changes
# ---------------------------------------------------------------------------
test_that("stress test: exclusion category counts match expected", {

  skip_if(!exists("stress_result"),
          "Skipping: stress_result not available (Test 1 must run first)")

  res <- stress_result
  catcount <- function(cat) {
    n <- res[exclude == cat, .N]
    if (length(n) == 0) return(0L)
    return(n)
  }

  # Print current counts for easy reference when updating
  cat("\n--- Stress test exclusion counts ---\n", file = stderr())
  counts <- res[, .N, by = exclude][order(-N)]
  for (i in seq_len(nrow(counts))) {
    cat(sprintf("  %-40s %5d\n",
                as.character(counts$exclude[i]), counts$N[i]),
        file = stderr())
  }
  cat("---\n", file = stderr())

  # Total rows
  expect_equal(nrow(res), 33101)

  # Freeze top categories — these are the values to update when the
  # algorithm intentionally changes. Run the test once with new code,
  # read the printed counts, and update here.
  #
  # Frozen counts — updated 2026-04-16: exclusion codes no longer param-specific;
  # per-param counts merged. Update these when algorithm intentionally changes.
  expect_equal(catcount("Include"), 29194)
  expect_equal(catcount("Exclude-C-Extraneous"), 929)     # was HT 312 + WT 617
  expect_equal(catcount("Exclude-C-Evil-Twins"), 734)      # was HT 630 + WT 104
  expect_equal(catcount("Exclude-C-Traj"), 613)            # was HT 360 + WT 253
  expect_equal(catcount("Exclude-C-BIV"), 606)             # was HT 545 + WT 61
  expect_equal(catcount("Exclude-C-Too-Many-Errors"), 401)
  expect_equal(catcount("Exclude-C-CF"), 378)              # was HT 209 + WT 169
  expect_equal(catcount("Exclude-C-Abs-Diff"), 156)
  expect_equal(catcount("Exclude-C-Traj-Extreme"), 58)     # was HT 39 + WT 19
  expect_equal(catcount("Exclude-C-Identical"), 31)        # was HT 17 + WT 14
  expect_equal(catcount("Exclude-Missing"), 1)

  # 11 distinct categories present (was 18 with param-specific codes)
  expect_equal(
    length(unique(as.character(res$exclude))), 11
  )
})

# ---------------------------------------------------------------------------
# Test 3: Per-archetype summary
# ---------------------------------------------------------------------------
test_that("stress test: all archetypes processed without crash", {

  skip_if(!exists("stress_result"),
          "Skipping: stress_result not available")

  res <- stress_result

  # Merge archetype info back (coerce to character for join compatibility)
  arch_map <- unique(stress_data[, .(subjid, archetype)])
  arch_map[, subjid := as.character(subjid)]
  res_copy <- copy(res)
  res_copy[, subjid := as.character(subjid)]
  res_arch <- merge(res_copy, arch_map, by = "subjid", all.x = TRUE)

  # Every archetype should have results
  archetypes_present <- unique(res_arch$archetype)
  expect_true(all(fixture$archetypes %in% archetypes_present),
              info = "All 8 archetypes should be in results")

  # Print per-archetype include rates
  cat("\n--- Per-archetype Include rates ---\n", file = stderr())
  arch_summary <- res_arch[, .(
    n = .N,
    n_include = sum(exclude == "Include"),
    pct_include = round(100 * sum(exclude == "Include") / .N, 1)
  ), by = archetype][order(archetype)]

  for (i in seq_len(nrow(arch_summary))) {
    cat(sprintf("  %-15s %5d rows, %5d Include (%5.1f%%)\n",
                arch_summary$archetype[i],
                arch_summary$n[i],
                arch_summary$n_include[i],
                arch_summary$pct_include[i]),
        file = stderr())
  }
  cat("---\n", file = stderr())

  # Each archetype should have some Includes (not all excluded)
  for (arch in fixture$archetypes) {
    arch_inc <- res_arch[archetype == arch & exclude == "Include", .N]
    expect_gt(arch_inc, 0,
              label = sprintf("Archetype '%s' should have some Includes", arch))
  }
})
