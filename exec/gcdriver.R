#!/usr/bin/env Rscript

library(argparser)
library(data.table)

library(growthcleanr)

parser <- arg_parser("CLI driver for growthcleanr")

parser <- add_argument(parser,
  "infile",
  type = "character",
  nargs = 1,
  help = "input file"
)
parser <- add_argument(parser,
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
parser <- add_argument(parser,
  "--quietly",
  flag = TRUE,
  help = "Disable verbose output"
)

argv <- parse_args(parser)
print(argv)

logfile <- sprintf("output/log/log-%s.txt", argv$infile)

if (argv$sdrecenter != "") {
  sdrecenter <- fread(argv$sdrecenter)
} else {
  sdrecenter <- ""
}

df_in <- fread(argv$infile)
df_out <- df_in[, exclude :=
  cleangrowth(
    subjid,
    param,
    agedays,
    sex,
    measurement,
    sd.recenter = sdrecenter,
    log.path = logfile,
    quietly = argv$quietly
  )]
fwrite(df_out, argv$outfile, row.names = FALSE)
