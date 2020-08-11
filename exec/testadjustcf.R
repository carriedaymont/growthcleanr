#!/usr/bin/env Rscript

### CARRY FORWARD ADJUSTMENT SCRIPT
# The goal of this script is to consider the height values that growthcleanr
# excludes as “carried forward” for potential re-inclusion by using a reverse
# absolute height velocity check based on step 15 of the Daymont et al.
# algorithm.

library(argparse, quietly = T)

parser <-
  ArgumentParser(description = 'CLI driver for carry forward adjustments')
parser$add_argument(
  'infile',
  metavar = 'INFILE',
  type = 'character',
  nargs = 1,
  help = 'Input file, cleaned by cleangrowth()'
)
parser$add_argument('--gridlength',
                    default = 9,
                    type = "integer",
                    help = "Number of steps in grid to search")
parser$add_argument('--quietly',
                    default = FALSE,
                    action = 'store_true',
                    help = 'Verbose debugging info')
parser$add_argument('--debug',
                    default = FALSE,
                    action = 'store_true',
                    help = 'Produce extra data files for debugging; use w/--outdir')
parser$add_argument("--maxrecs",
                    default = 0,
                    type = "integer",
                    help = "Limit to specified # recs, default 0 (no limit)")
parser$add_argument("--outdir",
                    default = "output",
                    type = "character",
                    help = "Directory for output files, default 'output' (no trailing slash)")

args <- parser$parse_args()

grid.length <- args$gridlength
outdir <- args$outdir
quietly <- args$quietly

# Load plyr before dplyr intentionally
library(plyr, quietly = T)
library(dplyr, quietly = T)

library(data.table, quietly = T)

library(growthcleanr)


######################## Prepare the data #############
# Read in data
dm <- fread(args$infile)
setkey(dm, subjid, param, agedays)

# limit dataset for debugging/performance reasons if specified
if (args$maxrecs > 0) {
  if (!quietly)
    cat(sprintf(
      "[%s] Taking max %s recs for debugging\n",
      Sys.time(),
      args$maxrecs
    ))
  dm_filt <- dm[subjid %in% unique(dm[, subjid])[1:args$maxrecs], ]
} else {
  dm_filt <- dm
}

# Create outdir if necessary; ignore warnings if it already exists
dir.create(file.path(outdir), showWarnings = FALSE)

dm_filt$n <- as.numeric(rownames(dm_filt))

subjid <- with(dm_filt, subjid)
param <- with(dm_filt, param)
sex <- with(dm_filt, sex)
agedays <- with(dm_filt, agedays)
measurement <- with(dm_filt, measurement)
orig.exclude <- with(dm_filt, exclude)

### Execute a sweep of parameters ###
if (!quietly)
  cat(sprintf(
    "[%s] Running parameter sweep of adjustcarryforward on cleaned dataset\n",
    Sys.time()
  ))

# Define the params to test
v_minfactor <- seq(0, 1, length = grid.length)
v_maxfactor <- seq(0, 4, length = grid.length)
v_banddiff <-  seq(0, 6, length = grid.length)
v_banddiff_plus <- seq(0, 11, length = grid.length)
v_min_ht.exp_under <- seq(0, 4, length = grid.length)
v_min_ht.exp_over <- seq(-1, 1, length = grid.length)
v_max_ht.exp_under <- seq(0, 0.66, length = grid.length)
v_max_ht.exp_over <- seq(0, 3, length = grid.length)

# Execute
for (index in 1:length(v_minfactor)) {
  if (!quietly)
    cat(sprintf(
      "[%s] Calling adjustcarryforward(), run %s\n",
      Sys.time(),
      index
    ))
  out <-
    adjustcarryforward(
      subjid,
      param,
      agedays,
      sex,
      measurement,
      orig.exclude,
      n,
      quietly = quietly,
      minfactor = v_minfactor[index],
      maxfactor = v_maxfactor[index],
      banddiff = v_banddiff[index],
      banddiff_plus = v_banddiff_plus[index],
      min_ht.exp_under = v_min_ht.exp_under[index],
      v_min_ht.exp_over[index],
      max_ht.exp_under = v_max_ht.exp_under[index],
      max_ht.exp_over = v_max_ht.exp_over[index]
    )

  setnames(out, "adjustcarryforward", sprintf("run-%s", index))
  export <-
    merge(as.data.frame(dm_filt), out, by = "n", all = T)

  # Combine into one wide set for simplicity
  if (index == 1) {
    combo <- export
  } else {
    combo <- merge(combo, out, by = "n", all = T) # %>% select(-n)
  }
}

# Combined adjusted set
fwrite(combo %>% select(-n), file.path(outdir, "all-adjusted.csv"), row.names = F)

# Record the sweep parameters for review
grid <-
  data.frame(
    run = 1:grid.length,
    minfactor = v_minfactor,
    maxfactor = v_maxfactor,
    banddiff = v_banddiff,
    banddiff_plus = v_banddiff_plus,
    min_ht.exp_under = v_min_ht.exp_under,
    min_ht.exp_over = v_min_ht.exp_over,
    max_ht.exp_under = v_max_ht.exp_under,
    max_ht.exp_over = v_max_ht.exp_over
  )
fwrite(grid, file.path(outdir, "params.csv"), row.names = F)

# Any additional debugging output
if (args$debug) {
  if (!quietly) {
    cat(sprintf("[%s] Writing debugging output files\n", Sys.time()))
  }
  # Filtered data
  fwrite(dm_filt, file.path(outdir, "dm_filt.csv"), row.names = F)
}
