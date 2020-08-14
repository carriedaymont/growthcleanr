# growthcleanr

## [1.2.1] - 2020-08-14

### Added

- New tests in `tests/testthat/test-utils.R` and `tests/testthat/test-cdc.R`
  to support newly added functions

### Changed

- Improved error handling in `longwide()`; fixed missing import in DESCRIPTION

## [1.2] - 2020-07-24

### Added

- New CDC BMI calculation function `ext_bmiz()`, comparable to SAS program
  published at https://www.cdc.gov/nccdphp/dnpao/growthcharts/resources/sas.htm
- Reference data file `inst/extdata/CDCref_d.csv` from CDC for use with
  `ext_bmiz()`
- New function `longwide()` for transforming `cleangrowth()` output for use with
  `ext_bmiz()`
- New function `recode_sex()` for recoding input data column values for `sex` to
  match `cleangrowth()` or `ext_bmiz()` requirements
- New `exec/gcdriver.R` command-line script for CLI execution of `cleangrowth()`
- New `Dockerfile` (and `.dockerignore`) enabling containerized use of
  `growthcleanr`
- Started test suite in `tests`
- New experimental function `adjustcarryforward()` in `R/adjustcarryforward.R`
  and driver script `exec/testadjustcf.R` (see README-adjustcarryforward.md for
  details)

### Changed

- Reorganized code from `R/growth.R` into separate files for clarity and easier
  maintenance (all utility functions not directly used by `cleangrowth()` are
  now in `R/utils.R`)
- Updated README with details and examples for added functions

## [1.1] - 2020-02-07

### Added

- New options to add flexibility:
  - `error.load.mincount` and `error.load.threshold`
  - `lt3.exclude.mode` with default (same as before) and `flag.both` mode for
    handling unmatched pairs
  - `sdmedian.filename` and `sdrecentered.filename`
- New `splitinput()` function
- New example synthetic data set `syngrowth` loads automatically.

### Changed

- Several updates to improve performance, including eliminating use of
  data.table in ewma function.
- Updated README with link to paper, detailed introduction, more installation
  details, examples, notes on handling large datasets, lists of parameters
  and exclusions.

## [1.0.0] - 2018-09-11

### Added

- Initial version posted to GitHub.
