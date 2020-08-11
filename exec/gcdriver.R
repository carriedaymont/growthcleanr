#!/usr/bin/env Rscript

library(argparse)
library(data.table)

library(growthcleanr)

parser <-
  ArgumentParser(description = 'CLI driver for growthcleanr')
parser$add_argument(
  'infile',
  metavar = 'INFILE',
  type = "character",
  nargs = 1,
  help = 'input file'
)
parser$add_argument(
  'outfile',
  metavar = 'OUTFILE',
  type = 'character',
  nargs = 1,
  help = 'output file'
)
parser$add_argument(
  '--sdrecenter',
  type = 'character',
  nargs = 1,
  default = '',
  help = 'sd.recenter data file'
)
parser$add_argument('--quietly',
                    default = F,
                    action = 'store_true',
                    help = 'Disable verbose output')
args <- parser$parse_args()

logfile <- sprintf('output/log/log-%s.txt', args$infile)

if (args$sdrecenter != '') {
  sdrecenter <- fread(args$sdrecenter)
} else {
  sdrecenter <- ''
}

df_in <- fread(args$infile)
df_out <- df_in[, exclude :=
                  cleangrowth(
                    subjid,
                    param,
                    agedays,
                    sex,
                    measurement,
                    sd.recenter = sdrecenter,
                    log.path = logfile,
                    quietly = args$quietly
                  )]
fwrite(df_out, args$outfile, row.names = FALSE)
