# generate_stress_test_data.R
#
# Generates a stress-test dataset for growthcleanr regression testing:
# - 400 subjects (50 per archetype × 8 archetypes)
# - 10 error types, each independently applied to 10% of patients
# - Some patients will accumulate multiple error types
# - Tracks errors at the measurement level (not just patient level)
#
# Output columns (in addition to standard GC input):
#   measurement_clean  — original value before any errors (NA for SDE-added rows)
#   err_cf ... err_rounding — binary flags for each error type applied to this row
#   is_sde_added       — TRUE for rows added by SDE scenarios
#
# Output: tests/testthat/stress_test_data.rds, stress_test_data.csv
#
# Run from the error-impact project directory:
#   cd "__Pipeline/error-impact"
#   Rscript "../gc-github-latest/scripts/generate_stress_test_data.R"
#
# Requires: error-impact source files, growthcleanr installed

library(data.table)

# Source error-impact modules
source(here::here("R", "config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "synthpop.R"))
source(here::here("R", "dataset.R"))
source(here::here("R", "induce_errors.R"))

cat("=== Generating stress test dataset ===\n")

# -------------------------------------------------------------------------
# Step 1: Build clean archetype population (50 per archetype = 400 total)
# -------------------------------------------------------------------------
cat("Building archetype population...\n")

ds <- build_clean_dataset(
  archetype_mix = c(
    tracker     = 50,
    falter      = 50,
    catchup     = 50,
    sgaAsym     = 50,
    sgaSym      = 50,
    preterm     = 50,
    latePreterm = 50,
    PGR         = 50
  ),
  include_hc    = TRUE,
  return_true_z = FALSE,
  seed          = 42L,
  verbose       = TRUE
)

n_subj <- length(unique(ds$subjid))
n_clean_rows <- nrow(ds)
cat(sprintf("  %d subjects, %d rows\n", n_subj, n_clean_rows))

# -------------------------------------------------------------------------
# Step 2: Define error types and apply each to 10% of patients
# -------------------------------------------------------------------------

# Each error type: list(scenario_name, overrides)
# For types needing both HT+WT, we apply two scenarios sequentially
error_type_names <- c("err_cf", "err_cf_chain", "err_sde_similar",
                      "err_sde_extreme", "err_small", "err_large",
                      "err_unit", "err_outlier", "err_swap", "err_rounding")

error_types <- list(
  err_cf          = list(scenarios = "cf_both_5p"),
  err_cf_chain    = list(scenarios = "cf_chain_both_5p"),
  err_sde_similar = list(scenarios = "sde_similar_both"),
  err_sde_extreme = list(scenarios = "sde_extreme_both"),
  err_small       = list(scenarios = c("pct_rand_ht_10", "pct_rand_wt_10")),
  err_large       = list(
    scenarios = c("pct_rand_ht_20", "pct_rand_wt_20"),
    overrides = list(pct_value = 0.50)
  ),
  err_unit        = list(scenarios = "unit_ht_in"),
  err_outlier     = list(scenarios = "outlier_ht_20cm"),
  err_swap        = list(scenarios = "swap_5p"),
  err_rounding    = list(scenarios = "rnd_both")
)

# Save clean measurement and initialize error tracking columns
ds[, measurement_clean := measurement]
ds[, is_sde_added := FALSE]
for (en in error_type_names) ds[, (en) := FALSE]

# Create a stable row key for matching (before SDE rows exist)
ds[, orig_row_key := paste(subjid, param, agedays, id, sep = "_")]

all_specs <- scenario_specs()
all_subjids <- unique(ds$subjid)
prop_affected <- 0.10
n_affected <- round(n_subj * prop_affected)

# Track which patients got which errors (patient-level summary)
error_assignments <- data.table(subjid = all_subjids, key = "subjid")

# Master seed for patient selection (not error generation)
set.seed(99L)

for (err_name in names(error_types)) {
  err_def <- error_types[[err_name]]
  scenarios <- err_def$scenarios
  overrides <- err_def$overrides  # NULL if not specified

  # Independently select 10% of patients
  selected <- sample(all_subjids, n_affected)
  error_assignments[subjid %in% selected, (err_name) := TRUE]

  cat(sprintf("  Applying %-16s to %d patients (%s)...\n",
              err_name, length(selected),
              paste(scenarios, collapse = " + ")))

  # Extract selected patients' rows, save pre-error measurements
  subset_before <- ds[subjid %in% selected]
  n_before <- nrow(subset_before)
  meas_before <- subset_before$measurement

  # Apply each scenario sequentially
  subset_err <- copy(subset_before)
  for (sc_name in scenarios) {
    spec <- all_specs[scenario == sc_name]
    if (nrow(spec) != 1) {
      stop(sprintf("Scenario '%s' not found or ambiguous", sc_name))
    }
    subset_err <- induce_errors(
      subset_err, spec,
      seed = 12345L,
      overrides = overrides
    )
  }

  # Identify which rows were modified vs added
  n_after <- nrow(subset_err)
  n_added <- n_after - n_before

  # Mark error flag on existing rows where measurement changed
  # (first n_before rows correspond to the original rows)
  changed <- subset_err$measurement[seq_len(n_before)] != meas_before
  changed[is.na(changed)] <- FALSE  # NA != x gives NA
  subset_before[changed, (err_name) := TRUE]

  # Update measurement on changed rows
  subset_before[, measurement := subset_err$measurement[seq_len(n_before)]]

  # Handle SDE-added rows (rows beyond n_before)
  if (n_added > 0) {
    added_rows <- subset_err[(n_before + 1):n_after]
    # Build new rows with tracking columns
    new_rows <- added_rows[, .(id, subjid, param, agedays, sex, measurement,
                               archetype)]
    new_rows[, measurement_clean := NA_real_]
    new_rows[, is_sde_added := TRUE]
    for (en in error_type_names) {
      new_rows[, (en) := (en == err_name)]
    }
    new_rows[, orig_row_key := NA_character_]

    # Add visit_type if it exists in ds
    if ("visit_type" %in% names(ds)) {
      new_rows[, visit_type := NA_character_]
    }

    # Ensure column alignment and append
    missing_cols <- setdiff(names(ds), names(new_rows))
    for (mc in missing_cols) {
      new_rows[, (mc) := ds[[mc]][NA_integer_]]
    }
    new_rows <- new_rows[, names(ds), with = FALSE]
    ds <- rbind(ds[!subjid %in% selected], subset_before, new_rows,
                use.names = TRUE)
  } else {
    ds <- rbind(ds[!subjid %in% selected], subset_before,
                use.names = TRUE)
  }
}

# Reassign ids (SDE errors added rows)
setorder(ds, subjid, param, agedays)
ds[, id := seq_len(.N)]
ds[, orig_row_key := NULL]  # no longer needed

cat(sprintf("  Final dataset: %d subjects, %d rows (%d clean + %d SDE-added)\n",
            length(unique(ds$subjid)), nrow(ds),
            sum(!ds$is_sde_added), sum(ds$is_sde_added)))

# Fill NA error assignment flags with FALSE
for (err_name in names(error_types)) {
  error_assignments[is.na(get(err_name)), (err_name) := FALSE]
}

# Summary stats
n_errors_per_patient <- rowSums(
  as.matrix(error_assignments[, names(error_types), with = FALSE])
)
cat(sprintf("\n  Patients with 0 errors: %d\n",
            sum(n_errors_per_patient == 0)))
cat(sprintf("  Patients with 1 error:  %d\n",
            sum(n_errors_per_patient == 1)))
cat(sprintf("  Patients with 2+ errors: %d\n",
            sum(n_errors_per_patient >= 2)))
cat(sprintf("  Max errors on one patient: %d\n",
            max(n_errors_per_patient)))

# Row-level error summary
any_err <- rowSums(as.matrix(ds[, error_type_names, with = FALSE])) > 0
cat(sprintf("\n  Rows with any error: %d / %d (%.1f%%)\n",
            sum(any_err), nrow(ds), 100 * sum(any_err) / nrow(ds)))
cat(sprintf("  SDE-added rows: %d\n", sum(ds$is_sde_added)))
cat(sprintf("  Rows rounded to 0 (will be Missing): %d\n",
            sum(ds$measurement == 0, na.rm = TRUE)))

# -------------------------------------------------------------------------
# Step 3: Save
# -------------------------------------------------------------------------
out_dir <- file.path(dirname(here::here()), "gc-github-latest",
                     "tests", "testthat")

# RDS with full metadata
saveRDS(
  list(
    data              = ds,
    error_assignments = error_assignments,
    n_subjects        = n_subj,
    n_rows            = nrow(ds),
    n_clean_rows      = n_clean_rows,
    archetypes        = c("tracker", "falter", "catchup", "sgaAsym",
                          "sgaSym", "preterm", "latePreterm", "PGR"),
    error_types       = names(error_types),
    error_type_names  = error_type_names,
    prop_affected     = prop_affected,
    generation_date   = Sys.Date(),
    generation_seed   = 42L
  ),
  file.path(out_dir, "stress_test_data.rds")
)

# CSV
data.table::fwrite(ds, file.path(out_dir, "stress_test_data.csv"))

cat(sprintf("\nSaved to: %s\n", out_dir))
cat(sprintf("  RDS: %.1f KB\n",
            file.size(file.path(out_dir, "stress_test_data.rds")) / 1024))
cat(sprintf("  CSV: %.1f KB\n",
            file.size(file.path(out_dir, "stress_test_data.csv")) / 1024))
cat("=== Done ===\n")
