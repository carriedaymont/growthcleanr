#!/usr/bin/env Rscript

library(argparser)
library(data.table)

library(growthcleanr)

parser <- arg_parser("CLI driver for growthcleanr")

parser <- add_argument(
  parser,
  "infile",
  type = "character",
  nargs = 1,
  help = "input file"
)
parser <- add_argument(
  parser,
  "outfile",
  type = "character",
  nargs = 1,
  help = "output file"
)
parser <- add_argument(
  parser,
  "--sdrecenter",
  type = "character",
  nargs = 1,
  default = "",
  help = "sd.recenter data file"
)
parser <- add_argument(
  parser,
  "--adult_cutpoint",
  type = "numeric",
  default = 20,
  help = "adult cutpoint"
)
parser <- add_argument(
  parser,
  "--weightcap",
  type = "numeric",
  default = Inf,
  help = "weight cap"
)
parser <- add_argument(
  parser,
  "--numbatches",
  type = "numeric",
  nargs = 1,
  default = 1,
  help = "Number of batches"
)
parser <- add_argument(
  parser,
  "--quietly",
  flag = TRUE,
  help = "Disable verbose output"
)
parser <- add_argument(
  parser,
  "--adult_split",
  type = "numeric",
  default = Inf,
  help = "Number of splits to run data on"
)

argv <- parse_args(parser)
print(argv)

log.path <- sprintf("output/log/%s/%s", Sys.Date(), basename(argv$infile))

if (argv$sdrecenter != "") {
  if (argv$sdrecenter == "nhanes") {
    sdrecenter = "nhanes"
  } else {
    sdrecenter <- fread(argv$sdrecenter)
  }
} else {
  sdrecenter <- ""
}

parallel <- argv$numbatches > 1
num.batches <- argv$numbatches

df_in <- fread(argv$infile)
df_in$id_order <- 1:nrow(df_in)

# handle adult split arguments
adult_split <- argv$adult_split

# we'll process the adult data using a split, then do the work
if (adult_split < Inf & length(unique(df_in$subjid)) < adult_split){
  adult_split <- length(unique(df_in$subjid))
}

if (adult_split < Inf & adult_split > 1){
  # will warn if they're not exact multiples
  # split by subject
  subj_split <- suppressWarnings(split(unique(df_in$subjid), 1:adult_split))
  # map indices to subjects
  # subj_split <- data.frame(
  #   entry = rep(seq_along(subj_split), length(subj_split)),
  #   subjid = as.character(unlist(subj_split))
  # )
  subj_split_df <- data.frame()
  for (spl in 1:length(subj_split)){
    subj_split_df <- rbind(
      subj_split_df,
      data.frame(
        "entry" = spl,
        "subjid" = as.character(subj_split[[spl]])
      ))
  }
  df_in$subjid <- as.character(df_in$subjid)

  # add batch to df_in
  df_in <- merge(df_in, subj_split_df, by = "subjid")

  # split based on subject id
  split.list <- suppressWarnings(split(df_in, df_in$entry))

} else {
  # no need for copy, they can refer to the same thing
  split.df <- df_in
  adult_split <- 1 # for Inf, just converting for later
}

# do split based on number of adult splits
df_out <- lapply(1:adult_split, function(x){
  if (!argv$quietly){
    cat(sprintf(
      "[%s] Processing adult data split %g, using %d batch(es)...\n",
      Sys.time(), x, num.batches))
  }

  if (adult_split > 1){
    split.df <- split.list[[x]]
  } # otherwise, split adult has already been created

  # Separate the logs or they'll overwrite each other
  split.log.path <- sprintf("%s/split-%03d", log.path, x)

  split.df$exclude <- cleangrowth(
    split.df$subjid,
    split.df$param,
    split.df$agedays,
    split.df$sex,
    split.df$measurement,
    sd.recenter = sdrecenter,
    weight_cap = argv$weightcap,
    adult_cutpoint = argv$adult_cutpoint,
    log.path = split.log.path,
    num.batches = num.batches,
    parallel = parallel,
    quietly = argv$quietly
  )

  return(split.df)
})

# put all the splits together
df_out <- rbindlist(df_out)
# get it in the original order
df_out <- df_out[order(df_out$id_order),]
# remove indexing variables
idx_var <- if (adult_split == 1){
  c("id_order")
} else {
  c("id_order", "entry")
}
df_out <- df_out[, !idx_var, with = FALSE]

fwrite(df_out, argv$outfile, row.names = FALSE)
