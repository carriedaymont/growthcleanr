# Test harness for adult growthcleanr R code
# Compares R output against expected values from adult-gc-test-ALL-PHASES.csv
# Usage: Rscript test_harness.R [loosest|looser|tighter|tightest]
#   Default: looser

library(data.table)

# When run from within the package, functions are available after loading
# Source directly if running standalone (outside devtools::load_all)
tryCatch({
  # Try package namespace first
  cleanadult
}, error = function(e) {
  source(file.path(here::here(), "R", "adult_support.R"))
  source(file.path(here::here(), "R", "adult_clean.R"))
})

# Parse command-line argument for permissiveness level
args <- commandArgs(trailingOnly = TRUE)
perm_level <- if (length(args) >= 1) args[1] else "looser"
valid_levels <- c("loosest", "looser", "tighter", "tightest")
if (!perm_level %in% valid_levels) {
  stop(paste0("Invalid permissiveness level: '", perm_level,
              "'. Must be one of: ", paste(valid_levels, collapse = ", ")))
}

# Load data
input_file <- file.path(here::here(), "inst", "testdata", "adult-gc-test-ALL-PHASES.csv")
df <- fread(input_file)

# Run algorithm at selected permissiveness level
cat(sprintf("Running algorithm at permissiveness = '%s'...\n", perm_level))
df <- cleanadult(df, permissiveness = perm_level, scale_max_lbs = 400)

# =============================================================================
# VALIDATION: Compare results to expected values
# =============================================================================

cat(sprintf("\n=== VALIDATION (permissiveness = %s) ===\n\n", perm_level))

# Load original test CSV for expected values
test_csv <- fread(input_file)
test_csv[, id := as.character(id)]

# Merge results back
results <- df[, .(id = as.character(id), result)]
merged <- merge(test_csv, results, by = "id", all.x = TRUE)

# Use the permissiveness-specific expected column
expected_col <- paste0("expected_", perm_level)

if (!expected_col %in% names(merged)) {
  stop(paste0("Column '", expected_col, "' not found in CSV. ",
              "Available columns: ", paste(names(merged), collapse = ", ")))
}

# Filter to rows where expected value is non-empty
check_rows <- merged[!is.na(get(expected_col)) & get(expected_col) != ""]
check_rows[, expected := get(expected_col)]
check_rows[, match := result == expected]

pass <- sum(check_rows$match, na.rm = TRUE)
total <- nrow(check_rows)
fail <- total - pass

cat(sprintf("Result: %d/%d passing (%d mismatches)\n", pass, total, fail))

if (fail > 0) {
  cat("\nMismatches:\n")
  mismatches <- check_rows[match == FALSE, .(id, subjid, param, measurement,
                                              expected, got = result,
                                              test_description)]
  print(mismatches, nrows = 100)
}

# Summary by expected code
cat("\n=== Summary by expected code ===\n")
summary <- check_rows[, .(pass = sum(match, na.rm = TRUE),
                           total = .N,
                           fail = sum(!match, na.rm = TRUE)),
                       by = expected]
summary <- summary[order(-fail, expected)]
print(summary)

cat(sprintf("\nTOTAL: %d/%d passing (%.1f%%)\n", pass, total, 100 * pass / total))
