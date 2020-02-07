# growthcleanr

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
