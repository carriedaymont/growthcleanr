#!/usr/bin/env Rscript

### CARRY FORWARD ADJUSTMENT SCRIPT
# The goal of this script is to consider the height values that growthcleanr
# excludes as “carried forward” for potential re-inclusion by using a reverse
# absolute height velocity check based on step 15 of the Daymont et al.
# algorithm.

# Specify libraries ----

library(argparse, quietly = T)
# Load plyr before dplyr intentionally
library(plyr, quietly = T)
library(dplyr, quietly = T)

library(data.table, quietly = T)

library(growthcleanr)

# Specify functions ----

# Function to make a grid vector for future sweep
# inputs:
# - low: number to start vector
# - high: number to end vector
# - grid.length: length of vector
# - searchtype: type of vector to make ("random","line-grid","full-grid")
# - is_included:
# - fg_val: value to use for non-used param in full-grid search
# outputs
# - numeric vector with specified qualities
make_grid_vect <- function(low,
                           high,
                           grid.length,
                           searchtype,
                           is_included = T,
                           fg_val = (low + high)/2){
  if (searchtype == "random"){ # random parameter sweep
    mid <- (low + high)/2
    # return half random below, midpoint, half random above
    gv <- c(runif(floor(grid.length/2), low, mid),
            mid,
            runif(floor(grid.length/2), mid, high))

    return(gv)
  } else if (searchtype == "line-grid") { # line-grid parameter sweep
    return(seq(low, high, length = grid.length))
  } else if (searchtype == "full-grid"){
    if (is_included){
      # full grid vectors are put together outside the list
      return(seq(low, high, length = grid.length))
    } else {
      if (is.na(fg_val)){
        return(c((low + high)/2))
      } else {
        return(c(fg_val))
      }
    }
  } else {
    stop("Invalid search type")
  }
}

# Function to execute sweep of parameter values for adjustcarryforward
# inputs:
# - grid_df: data frame of parameters to sweep
# outputs:
# - combined input file and result of each run
exec_sweep <- function(grid_df){
  for (index in 1:nrow(grid_df)) {
    if (!quietly)
      cat(sprintf(
        "[%s] Calling adjustcarryforward(), run %s of %s\n",
        Sys.time(),
        index,
        nrow(grid_df)
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
        minfactor = grid_df$minfactor[index],
        maxfactor = grid_df$maxfactor[index],
        banddiff = grid_df$banddiff[index],
        banddiff_plus = grid_df$banddiff_plus[index],
        min_ht.exp_under = grid_df$min_ht.exp_under[index],
        min_ht.exp_over = grid_df$min_ht.exp_over[index],
        max_ht.exp_under = grid_df$max_ht.exp_under[index],
        max_ht.exp_over = grid_df$max_ht.exp_over[index]
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

  return(combo)
}

# Specify parser ----

parser <-
  ArgumentParser(description = 'CLI driver for carry forward adjustments')
parser$add_argument(
  'infile',
  metavar = 'INFILE',
  type = 'character',
  nargs = 1,
  help = 'Input file, cleaned by cleangrowth()'
)
parser$add_argument('--searchtype',
                    default = 'random',
                    type = "character",
                    help = "Type of search to perform: random (default), line-grid, full-grid")
parser$add_argument('--gridlength',
                    default = 9,
                    type = "integer",
                    help = "Number of steps in grid to search")
parser$add_argument('--seed',
                    default = 7,
                    type = "integer",
                    help = "Random seed, used only when performing random search")
parser$add_argument("--param",
                    default = "none",
                    type = "character",
                    help = "CSV to specify which parameters to run full search on, and values to use if not, used only when performing full-grid search")
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
parser$add_argument("--outfile",
                    default = paste0("test_adjustcarryforward_",
                                     format(Sys.time(), "%m-%d-%Y_%H-%M-%S")),
                    type = "character",
                    help = "Output file name, default 'test_adjustcarrforward_DATE_TIME', where DATE is the current system date and time")

# Parse arguments for ease ----

args <- parser$parse_args()

seed <- args$seed
searchtype <- args$searchtype
grid.length <- args$gridlength
outfile <- args$outfile
quietly <- args$quietly
param_choice <- args$param

if (searchtype == "random"){
  # if the line search isn't odd, add one so that it includes the midpoint
  grid.length <- if (grid.length %% 2 == 0){ grid.length + 1 } else { grid.length }
}

# Prepare the data ----

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

dm_filt$n <- as.numeric(rownames(dm_filt))

subjid <- with(dm_filt, subjid)
param <- with(dm_filt, param)
sex <- with(dm_filt, sex)
agedays <- with(dm_filt, agedays)
measurement <- with(dm_filt, measurement)
orig.exclude <- with(dm_filt, clean_value)

# Execute parameter sweep ----

if (!quietly){
  cat(sprintf(
    "[%s] Running parameter sweep of adjustcarryforward on cleaned dataset\n",
    Sys.time()
  ))
}

# Specify seed, if using random
if (searchtype == "random"){
  set.seed(seed)
}

# Have a data frame of ranges
param_n <- c("minfactor",
             "maxfactor",
             "banddiff",
             "banddiff_plus",
             "min_ht.exp_under",
             "min_ht.exp_over",
             "max_ht.exp_under",
             "max_ht.exp_over")

p_range <- data.frame(
  "param" = param_n,
  low = c(0,0,0,0,0,-1,0,0),
  high = c(1,4,6,11,4,1,0.66,3)
)
rownames(p_range) <- p_range$param

# Define the params to test
if (searchtype == "full-grid"){
  # specify parameter choices -- default to choosing all
  if (param_choice != "none"){
    param_df <- read.csv(param_choice)
    colnames(param_df) <- c("param", "include", "value")
  } else {
    param_df <- data.frame(
      "param" = param_n,
      "include" = T,
      "value" = NA
    )
  }

  if (searchtype == "full-grid"){
    cat(paste0(
      "[",Sys.time(),
      "] Note: due to the search type specified (full-grid), there will be ",
      grid.length^sum(param_df$include)," runs\n"
    ))
  }

  # now pass them in one by one
  grid_list <- list()
  for (i in 1:nrow(param_df)){
    grid_list[[param_df$param[i]]] <-
      make_grid_vect(p_range[param_df$param[i], "low"],
                     p_range[param_df$param[i], "high"],
                     grid.length,
                     searchtype,
                     as.logical(param_df$include[i]),
                     as.numeric(param_df$value[i]))
  }

  # make the full grid based on the list
  grid_df <- expand.grid(grid_list)
} else { # for random and line-grid options
  grid_df <- as.data.frame(matrix(NA, nrow = grid.length, ncol = nrow(p_range)))
  colnames(grid_df) <- p_range$param
  for (i in 1:ncol(grid_df)){
    grid_df[,p_range$param[i]] <-
      make_grid_vect(p_range[p_range$param[i], "low"],
                     p_range[p_range$param[i], "high"],
                     grid.length,
                     searchtype)
  }
}


# Execute
combo <- exec_sweep(grid_df)

# Write out results ----

# Combined adjusted set
fwrite(combo %>% select(-n), paste0(outfile, ".csv"), row.names = F)

# Record the sweep parameters for review
grid <- cbind("run" = 1:nrow(grid_df), grid_df)
fwrite(grid, paste0(outfile, "_parameters.csv"), row.names = F)

# Any additional debugging output
if (args$debug) {
  if (!quietly) {
    cat(sprintf("[%s] Writing debugging output files\n", Sys.time()))
  }
  # Filtered data
  fwrite(dm_filt, paste0(outfile, "_debug_filtered_data.csv"), row.names = F)
}
