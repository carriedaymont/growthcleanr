
<!-- README.md is generated from README.Rmd. Please edit that file -->

# growthcleanr

<!-- badges: start -->

[![R build
status](https://github.com/carriedaymont/growthcleanr/workflows/R-CMD-check/badge.svg)](https://github.com/carriedaymont/growthcleanr/actions)
<!-- badges: end -->

R package for cleaning data from Electronic Health Record systems,
focused on cleaning height and weight measurements.

<a name="cite"></a> This package implements the [Daymont et
al. algorithm](https://academic.oup.com/jamia/article/24/6/1080/3767271),
as specified in Supplemental File 3 within the [Supplementary
Material](https://academic.oup.com/jamia/article/24/6/1080/3767271#97610899)
published with that paper.

> Carrie Daymont, Michelle E Ross, A Russell Localio, Alexander G Fiks,
> Richard C Wasserman, Robert W Grundmeier, Automated identification of
> implausible values in growth data from pediatric electronic health
> records, Journal of the American Medical Informatics Association,
> Volume 24, Issue 6, November 2017, Pages 1080–1087,
> <https://doi.org/10.1093/jamia/ocx037>

This package also includes an R version of the [SAS macro published by
the
CDC](https://www.cdc.gov/nccdphp/dnpao/growthcharts/resources/sas.htm)
for calculating percentiles and Z-scores of pediatric growth
observations and utilities for working with both functions.

## Summary

The `growthcleanr` package processes data prepared in a specific format
to identify biologically implausible height and weight measurements. It
bases these evaluations on techniques which use patient-specific
longitudinal analysis and variations from published growth trajectory
charts. These techniques are performed in a specific order which refines
and improves results throughout the process.

Results from `growthcleanr` include a flag for each measurement
indicating whether it is to be included or excluded based on
plausibility, with a variety of specific types of exclusions identified
distinctly. These flags can be analyzed further by researchers studying
growth data to determine which measurements to include or exclude in
their own studies. No values are deleted or otherwise removed; each is
only flagged in a new column.

To start running `growthcleanr`, an R installation with a variety of
additional packages is required, as is a growth measurement dataset
prepared for use in `growthcleanr`.

The rest of this documentation includes:

  - [Quickstart](#quickstart) - a brief tour of using growthcleanr,
    including data preparation
  - [Installation](#installation) - including required dependencies and
    OS-specific notes
  - [Usage](#usage) - examples of cleaning data, multiple options,
    example data
  - [Configuration options](#configuration) - changing growthcleanr
    operational settings
  - [Understanding output](#output) - the exclusion types growthcleanr
    identifies
  - [Working with large datasets](#largedata) - notes and suggestions
  - [Computing BMI percentiles and Z-scores](#bmi) - additional
    functions for determining percentiles and Z-scores using the CDC
    method
  - [Related tools](#related) - additional software complementary to
    growthcleanr

## <a name="quickstart"></a>Quickstart

### <a name="setup"></a>R setup

The `growthcleanr` package must be installed, which will in turn install
the following packages:

  - `data.table`
  - `foreach`
  - `doParallel`
  - `labelled`
  - `plyr`

In addition, an R environment with `devtools` installed is also
required. Further installation details and notes can be found under
[Installation](#installation) below.

### Data preparation

Once the `growthcleanr` package and its dependencies are installed, the
minimal data vectors needed are:

1.  `subjid` - a unique identifier for each subject; if you have long
    (64-bit) integer `subjid` values, add the `bit64` R package to your
    install list to ensure R will process these correctly
2.  `param` - a string designating the measurement being ‘`HEIGHTCM`’,
    ‘`LENGTHCM`’, or ‘`WEIGHTKG`’
3.  `agedays` - a numeric value for age in days at time of measurement;
    will be rounded down (floor) to nearest integer if needed (note that
    there are 30.4375 days in one month; see also the additional [notes
    on calculating
    `agemos`](https://www.cdc.gov/nccdphp/dnpao/growthcharts/resources/sas.htm)
    from whole numbers of months
4.  `sex` - a numeric value of `0` (male) or `1` (female)
5.  `measurement` - numeric value of measurement in cm for height or kg
    for weight

**Note**: the dataset should be sorted by `subjid`, `param`, and
`agedays`, as demonstrated in the Example below.

### <a name="example"></a>Example

Starting with a dataset such as the following:

    #>    id     subjid sex agedays    param measurement
    #> 1   1 1575544570   1    2177 HEIGHTCM      104.40
    #> 2   2 1575544570   1    2548 HEIGHTCM      117.70
    #> 3   3 1575544570   1    2548 HEIGHTCM      117.74
    #> 4   4 1575544570   1    2919 HEIGHTCM      117.70
    #> 5   5 1575544570   1    3290 HEIGHTCM      129.20
    #> 6   6 1575544570   1    3661 HEIGHTCM      133.90
    #> 7   7 1575544570   1    3661 HEIGHTCM      133.64
    #> 8   8 1575544570   1    4032 HEIGHTCM      139.60
    #> 9   9 1575544570   1    4032 HEIGHTCM      140.02
    #> 10 10 1575544570   1    4403 HEIGHTCM      139.60
    #> 11 11 1575544570   1    4774 HEIGHTCM      152.90
    #> 12 12 1575544570   1    4774 HEIGHTCM      153.17
    #> 13 13 1575544570   1    5145 HEIGHTCM      152.90
    #> 14 14 1575544570   1    5516 HEIGHTCM      157.80
    #> 15 15 1575544570   1    2177 WEIGHTKG       18.90

Note that in this Example, we see data from two subjects, sorted as
described above, with two pairs of height/weight measurements for the
first patient and four pairs for the second. Most of these values pass
the “eyeball test”, but the third weight measurement for `subjid` 2
looks like it might be a carried forward value, given that their same
ageday height has changed by several centimeters and their third height
value is exactly the same as their second. This is an example of an
implausible value and exclusion type that `growthcleanr` might identify.

The `cleangrowth()` function provided by `growthcleanr` will return a
vector that specifies, for each measurement, whether to:

  - “`Include`” the measurement,
  - Conclude that the measurement is “`Missing`”,
  - Or to “`Exclude`” the measurement (numerous variations of `Exclude`
    exist, as specified below in Understanding output).

For a data.frame object `source_data` containing growth data:

``` r
library(growthcleanr)
library(data.table)

source_data <- as.data.frame(syngrowth)
source_data <- source_data[source_data$subjid %in% unique(source_data[, "subjid"])[1:5], ]

# prepare data as a data.table
data <- as.data.table(source_data)

# set the data.table key for better indexing
setkey(data, subjid, param, agedays)

# generate new exclusion flag field using function
cleaned_data <- data[, clean_value := cleangrowth(subjid, param, agedays, sex, measurement)]

# extract data limited only to values flagged for inclusion:
only_included_data <- cleaned_data[clean_value == 'Include']
```

If our Example dataset above were named `source_data`, examining
`cleaned_data` would show:

    #>     id    subjid sex agedays    param measurement             clean_value
    #>  1: 56 881283195   1     889 HEIGHTCM      89.300                 Include
    #>  2: 57 881283195   1    1071 HEIGHTCM      98.100                 Include
    #>  3: 58 881283195   1    1211 HEIGHTCM     100.600                 Include
    #>  4: 59 881283195   1    1253 HEIGHTCM     100.600 Exclude-Carried-Forward
    #>  5: 60 881283195   1    1253 HEIGHTCM     102.190                 Include
    #>  6: 61 881283195   1    1435 HEIGHTCM     100.600 Exclude-Carried-Forward
    #>  7: 62 881283195   1    1806 HEIGHTCM     112.900                 Include
    #>  8: 63 881283195   1    2177 HEIGHTCM     116.320                 Include
    #>  9: 64 881283195   1    2548 HEIGHTCM     116.320 Exclude-Carried-Forward
    #> 10: 65 881283195   1    2919 HEIGHTCM     134.100                 Include
    #> 11: 66 881283195   1    3290 HEIGHTCM     134.100 Exclude-Carried-Forward
    #> 12: 67 881283195   1    3290 HEIGHTCM     140.350                 Include
    #> 13: 68 881283195   1     889 WEIGHTKG      15.500                 Include
    #> 14: 69 881283195   1    1071 WEIGHTKG      16.700                 Include
    #> 15: 70 881283195   1    1071 WEIGHTKG      16.976       Exclude-Duplicate

In this sample output, we see that `growthcleanr` has indeed flagged the
4 weight measurement for `subjid` 881283195 as an apparent carried
forward value. The measurement remains in the dataset, and it is up to
the researcher to determine which exclusions to accept directly or
investigate further.

To compute percentiles and Z-scores using the included CDC macro,
continue with these next steps:

``` r
# Use the built-in utility function to convert the input observations to wide
# format for BMI calculation
cleaned_data_wide <- long_wide(cleaned_data)

# Compute Z-scores and percentiles
cleaned_data_bmi <- ext_bmiz(cleaned_data_wide)
```

The `long_wide()` function defaults to converting only records included
by `cleangrowth()`, with options to change this default, and
`ext_bmiz()` will return the wide version of the cleaned data with many
metrics added. Both of these functions are [described in detail
below](#bmi).

## <a name="installation"></a>Installation

`growthcleanr` has been developed and tested using R version 3 (most
recently 3.6+). It should work using R on Windows, macOS, or Unix/Linux,
although there are some additional [platform-specific notes](#platform)
you may wish to review.

In an up-to-date R environment such as RStudio, first install the
dependencies:

``` r
# Required packages for growthcleanr operation
install.packages("devtools")
# Optional additional packages referenced in this document
install.packages("argparse")
install.packages("bit64")
```

Note that this might take some time.

You can install the `growthcleanr` package directly from GitHub using
`devtools` in the R console with:

``` r
devtools::install_github("carriedaymont/growthcleanr")
```

### Source-level install for developers

If you want to work with and potentially change the `growthcleanr` code
itself, you can download or clone the `growthcleanr` source code and
then install it from source. To clone the source using `git`:

``` bash
git clone https://github.com/carriedaymont/growthcleanr.git
```

Either way, once you have the `growthcleanr` package source, open an R
session from the `growthcleanr` base directory. Then install
growthcleanr using the R `devtools` package:

``` r
devtools::install(".")
```

You can also install the package from an installation file if one is
obtained.

### <a name="platform"></a>Platform-specific notes

#### Windows

We have observed an issue with using `growthcleanr` with Windows in a
large agency where some R packages are installed by an administrator for
a user who may not have write access permissions on the package folder.
Similar issues may occur with some networked drives. This caused
problems with `growthcleanr`’s parallel processing option. If possible,
install R and its packages in locations hosted on the same local machine
and folder(s) for which the primary user has write permissions. These
steps should help to avoid this problem.

Some users have reported that `growthcleanr` runs more slowly on Windows
compared with Linux or macOS.

#### macOS

`growthcleanr` uses the
[`data.table`](https://rdatatable.gitlab.io/data.table/) package for R
extensively. `data.table` provides a faster version of R’s data frames,
and is used to improve `growthcleanr` performance. Typically
`data.table` installs in a manner that will be able to take advantage of
multiple threads. You will know it worked successfully if, when you load
the `data.table` library in R, you see something like the following:

``` r
library(data.table)
```

That data.table reports “using 2 threads” indicates that installation
has succeeded. If the message reports using only one thread, see the
advice under the [“OpenMP enabled compiler for Mac”
instructions](https://github.com/Rdatatable/data.table/wiki/Installation#openmp-enabled-compiler-for-mac)
to re-install `data.table`.

Users have reported errors running multiple batches with `parallel =
TRUE` from within RStudio. If this happens, the problem may be resolved
by running from RGui or from the command line using `Rscript`. An
example standalone script that may be used for this purpose is
documented below under [Working with large data sets](#largedata).

#### Docker

This package includes a `Dockerfile` that enables easy installation of R
and `growthcleanr` on a machine with [Docker](https://www.docker.com/)
installed. It requires an up-to-date Docker install, and a few
command-line steps, but can save time over installing R and
`growthcleanr`’s dependencies manually.

To install and run `growthcleanr` using Docker, open the PowerShell on
Windows, or open the Terminal on macOS, and enter this `docker` command:

``` bash
docker run -it mitre/growthcleanr:latest R
```

The first time this command is run, it might take a few minutes to
download and extract several necessary components, but this should be
fully automated. If successful, you should see an R prompt, from which
you can use `growthcleanr` as described below.

This R environment is virtualized inside Docker, however, and isolated
from your local machine. Because of this, you will need to map a local
folder on your computer into the Docker environment to work with your
own data. For example, if you are on Windows, and your data is in
`C:\Users\exampleuser\analysis`, specify a mapping using the added `-v`
step below:

``` bash
docker run -it -v C:\Users\exampleusers\analysis:/usr/src/app mitre/growthcleanr:latest R
```

Note that the slashes in file paths reverse direction from the reference
to the folder location on your Windows machine (before the colon) to the
folder location on the Docker container (after the colon); this is
intentional, and accounts for how the two different environments
reference disk locations.

Note also that when mapping a folder on Windows, you may be prompted to
confirm that you indeed want to “Share” the folder. This is a standard
Windows security practice, and it is okay to confirm and proceed.

If you are on macOS, and your data is in `/Users/exampleuser/analysis`,
specify a folder mapping like this:

``` bash
docker run -it -v /Users/exampleuser/analysis:/usr/src/app mitre/growthcleanr:latest R
```

If you mapped a folder, then inside the Docker environment’s R prompt,
when you then issue a command like `list.files()`, you should see a list
of the same files in the R session that you see in that folder on your
desktop. You can now open and read your data files, run `cleangrowth()`
and other analyses, and write result files to that same directory.

Exit the Docker R environment with `quit()` as you normally would. Any
new files you saved will appear in the desktop folder you mapped.

## <a name="usage"></a>Usage

With the `growthcleanr` package installed [as described
above](#installation), it can be loaded as a package:

``` r
library(growthcleanr)
```

For convenience, `growthcleanr` ships with an example synthetic dataset
created using [Synthea](https://synthetichealth.github.io/synthea/),
with support for simulated pediatric growth measurement errors based on
the protocol included as supplementary material to the [Daymont et
al. paper](#cite). This dataset is automatically loaded with the
`growthcleanr` library, and is called `syngrowth`.

``` r
dim(syngrowth)
#> [1] 56703     6
head(syngrowth)
#>   id     subjid sex agedays    param measurement
#> 1  1 1575544570   1    2177 HEIGHTCM      104.40
#> 2  2 1575544570   1    2548 HEIGHTCM      117.70
#> 3  3 1575544570   1    2548 HEIGHTCM      117.74
#> 4  4 1575544570   1    2919 HEIGHTCM      117.70
#> 5  5 1575544570   1    3290 HEIGHTCM      129.20
#> 6  6 1575544570   1    3661 HEIGHTCM      133.90
```

It can be used as in the example above. Note that processing the example
data with `cleangrowth()` will likely take a few minutes to complete.

``` r
data <- as.data.table(syngrowth)
setkey(data, subjid, param, agedays)
cleaned_data <- data[, clean_value := cleangrowth(subjid, param, agedays, sex, measurement)]
cleaned_data[, .N, by = "clean_value"][order(-N)]
#>                   clean_value     N
#>  1:                   Include 38875
#>  2:         Exclude-Duplicate 10546
#>  3:   Exclude-Carried-Forward  6694
#>  4:         Exclude-SD-Cutoff   168
#>  5:            Exclude-EWMA-8   135
#>  6:      Exclude-EWMA-Extreme    95
#>  7:            Exclude-EWMA-9    93
#>  8: Exclude-Min-Height-Change    65
#>  9:      Swapped-Measurements    16
#> 10:   Exclude-Too-Many-Errors     6
#> 11:           Exclude-EWMA-11     5
#> 12:     Exclude-Pair-Delta-18     2
#> 13:           Exclude-EWMA-12     2
#> 14: Exclude-Max-Height-Change     1
```

If you are able to run these steps and see a similar result, you have
the `growthcleanr` package installed correctly. The resulting
`cleaned_data` can be reviewed, subsetted, and compared in more detail
using all the tools R provides.

For complete information about the options that can be set on the
`cleangrowth()` function, see [Configuration](#configuration). Below are
a few additional examples:

This example shows three configuration options in use:

``` r
cleaned_data <- data[, 
  clean_value_both := cleangrowth(
    subjid, param, agedays, sex, measurement,
    lt3.exclude.mode = "flag.both",
    ref.data.path = "inst/extdata",
    quietly = FALSE
  )
]
```

  - `lt3.exclude.mode = "flag.both"` will exclude both measurements for
    a subject if they only have two unexcluded measurements of one
    parameter type with at least one implausible value and no same
    age-day measurements of the other parameter.

  - `ref.data.path = "inst/extdata"` shouldn’t be necessary if you have
    the `growthcleanr` package installed, but if you are running it from
    its source directly you may need to specify its full path.

  - `quietly = FALSE` enables verbose output, marking the progress of
    the algorithm through its many processing steps. This can be very
    helpful while testing.

This example shows built-in options for processing data in batches,
which can speed the process while working with large data sets:

``` r
cleaned_data <- data[, 
  clean_value_both := cleangrowth(
    subjid, param, agedays, sex, measurement,
    parallel = TRUE,
    num.batches = 4,
    log.path = "logs"
  )
]
```

The best `num.batches` value for your environment may vary depending on
the computing resources you have available. If you do not specify
`num.batches`, `growthcleanr` will estimate a batch count based on R
functions for checking the system hardware. If you are working with
large datasets, you may need to experiment with these options to
determine the best settings for your needs. Note that there are
additional notes on [working with large datasets](#largedata) below.

The default value of `log.path` is `"."`, the current working directory.
`growthcleanr` will write batch-specific log files to the `log.path`
directory, and will create the directory if necessary.

Note that if you run `growthcleanr` in parallel with multiple batches
and see warning errors such as the following, they can be ignored.

``` r
#> Warning messages:
#> 1: <anonymous>: ... may be used in an incorrect context: ‘.fun(piece, ...)’
#> 2: <anonymous>: ... may be used in an incorrect context: ‘.fun(piece, ...)’
```

This set of warnings is coming from one of `growthcleanr`’s
dependencies, and does not indicate either failure or improper
execution.

## <a name="configuration"></a>Configuration options

`growthcleanr` offers several options for configuration, with default
values set to address common cases. You may wish to experiment with the
settings to discover which work best for your research and your dataset.

All of the following options may be set as additional parameters in the
call to `cleangrowth()`.

### Algorithm-related options

The following options change the behavior of the growthcleanr algorithm.

  - `recover.unit.error` - default `FALSE`; when `FALSE`, measurements
    identified as unit errors (e.g., apparent height values in inches
    instead of centimeters) will be flagged but not corrected, when
    `TRUE` these values will be corrected and included as valid
    measurements for cleaning.

  - `sd.extreme` - default `25`; a very extreme value check on modified
    (recentered) Z-scores used as a first-pass elimination of clearly
    implausable values, often due to misplaced decimals.

  - `z.extreme` - default `25`; similar usage as `sd.extreme`, for
    absolute Z-scores.

**NOTE**: many different steps in the `growthcleanr` algorithm use
highly refined techniques to identify implausible values that require
more subtlety to detect. These default values for `sd.extreme` and
`z.extreme` are set very high by design to eliminate completely
implausible values early in the process; these extreme values require
none of `growthcleanr`’s additional, more subtle approaches to be
identified as exclusions. If `sd.extreme` and `z.extreme` are configured
to be lower values, measurements eliminated at the early step which
checks against these extreme values will not be further considered using
later techniques.

  - `include.carryforward` - default `FALSE`; if set to `TRUE`,
    `growthcleanr` will skip algorithm step 9, which identifies carried
    forward measurements, and will not flag these values for exclusion.

  - `ewma.exp` - default `-1.5`; the exponent used for weighting
    measurements when calculating exponentially weighted moving average
    (EWMA). This exponent should be negative to weight growth
    measurements closer to the measurement being evaluated more
    strongly. Exponents that are further from zero (e.g., `-3`) will
    increase the relative influence of measurements close in time to the
    measurement being evaluated compared to using the default exponent.

  - `ref.data.path` - defaults to using CDC reference data from year
    2000; supply a file path to use alternate reference data. Note that
    when running from an installed `growthcleanr` package (e.g. having
    called `library(growthcleanr)`), this path does not need to be
    specified. Developers testing the source code directly from the
    source directory will need to specify this as well.

  - `error.load.mincount` - default `2`; minimum count of exclusions on
    parameter for one subject before considering excluding all
    measurements.

  - `error.load.threshold` - default `0.5`; threshold to exceed for
    percentage of excluded measurement count relative to count of
    included other measurement (e.g., if 3 of 5 WTs are excluded, and 5
    corresponding HTs are included, this exceeds the 0.5 threshold and
    the other two WTs will be excluded).

  - `lt3.exclude.mode` - default `default`; determines type of exclusion
    procedure to use for 1 or 2 measurements of one type without
    matching same ageday measurements for the other parameter. Options
    include:
    
      - `default` - standard growthcleanr approach
      - `flag.both` - in case of two measurements with at least one
        beyond thresholds, flag both instead of one (as in default)

  - `sd.recenter` - defaults to NA; data frame or table w/median
    SD-scores per day of life by gender and parameter. Columns must
    include param, sex, agedays, and sd.median (referred to elsewhere as
    “modified Z-score”). By default, median values will be calculated
    using growth data to be cleaned.

### Operational options

The following options change the execution of the program overall, with
no effect on the algorithm itself.

  - `parallel` - default `FALSE`; set to `TRUE` to run `growthcleanr` in
    parallel. Running in parallel will split the input data into batches
    (while keeping all records for each subject together) and process
    each batch on a different processor/core to maximize throughput.
    Recommended for large datasets with more than 100K rows.

  - `num.batches` - specifies the number of batches to run in parallel;
    only applies if `parallel` is set to `TRUE`. Defaults to the number
    of workers returned by the `getDoParWorkers()` function in the
    `foreach` package. Note that processing in parallel may affect
    overall system performance.

  - `sdmedian.filename` - filename for optionally saving sd.median data
    calculated on the input dataset to as CSV. Defaults to `""`, for
    which this data will not be saved. Use for extracting medians for
    parallel processing scenarios other than the built-in parallel
    option. See notes on [large data sets](#largedata) for details.

  - `sdrecentered.filename` - filename to save re-centered data to as
    CSV. Defaults to "", for which this data will not be saved. Useful
    for post-processing and debugging.

  - `quietly` - default `TRUE`; when `TRUE`, displays function messages
    and will output log files when `parallel` is `TRUE`.

  - `log.path` - default `"."`; sets directory for batch log file output
    when processing with `parallel = TRUE`. A new directory will be
    created if necessary.

## <a name="output"></a>Understanding output

As seen [in the example above](#example), `growthcleanr` creates a
vector of flags identifying which measurements should be included or
excluded, and why. This vector can be added to the original dataset for
further review. The flags `growthcleanr` sets are text labels, most of
which are descriptive. The text does not always map directly back to the
precise step in the algorithm.

For convenience, the following table maps which `algorithm step` in the
algorithm, as defined by Supplemental File 3 within the [Supplementary
Material](https://academic.oup.com/jamia/article/24/6/1080/3767271#97610899)
of the Daymont et al. paper, define exclusions in the algorithm by the
exclusion identifier also prescribed in the algorithm (“`algorithm exc
id`”). The `exclusion type` column identifies the text of the exclusion
flag set by `growthcleanr`. The `notes` column at right identifies minor
discrepancies between the algorithm’s step labels and labels used in
comment text in `growthcleanr`.

| algorithm step | algorithm exc id | exclusion type                           | notes                       |
| -------------- | ---------------- | ---------------------------------------- | --------------------------- |
| 2d             | 0                | Include                                  | \-                          |
| 2d             | 1                | Missing                                  | \-                          |
| 5b             | 2                | Exclude-Temporary-Duplicate              | \-                          |
| 7d             | \-               | Swapped-Measurement                      | \-                          |
| 8f             | \-               | Unit-Error-High                          | \-                          |
| 8f             | \-               | Unit-Error-Low                           | \-                          |
| 8f             | \-               | Unit-Error-Possible                      | Not set in R code           |
| 9c             | 3                | Exclude-Carried-Forward                  | \-                          |
| 10c            | 4                | Exclude-SD-Cutoff                        | 10d, 10e                    |
| 11d            | 5                | Exclude-EWMA-Extreme                     | 11e                         |
| 11f.ii         | 6                | Exclude-EWMA-Extreme-Pair                | 11i (R only)                |
| 12d.i          | 7                | Exclude-Duplicate                        | 12diii, 12ei, 12f           |
| 14f.i          | 8                | Exclude-EWMA-8                           | Set in 14h (in R)           |
| 14f.ii         | 9                | Exclude-EWMA-9                           | Set in 14h (in R)           |
| 14f.iii        | 10               | Exclude-EWMA-10                          | Set in 14h (in R)           |
| 14f.iv         | 11               | Exclude-EWMA-11                          | Set in 14h (in R)           |
| 14f.v          | 12               | Exclude-EWMA-12                          | Set in 14h (in R)           |
| 14f.vi         | 13               | Exclude-EWMA-13                          | Set in 14h (in R)           |
| 14f.vii        | 14               | Exclude-EWMA-14                          | Set in 14h (in R)           |
| 15s            | 15               | Exclude-Min-Height-Change                | Set in 15q (in R), 15s, 15t |
| 15s            | 16               | Exclude-Max-Height-Change                | Set in 15q (in R), 15s, 15t |
| 16b            | 17               | Exclude-Pair-Delta-17                    | \-                          |
| 16b            | 18               | Exclude-Pair-Delta-18                    | \-                          |
| \-             | \-               | Exclude-Pair-Delta-19                    | Added `flag.both` mode      |
| 16d            | 19               | Exclude-Single-Outlier                   | \-                          |
| 17a            | 20               | Exclude-Too-Many-Errors                  | \-                          |
| 17a            | 21               | Exclude-Too-Many-Errors-Other\_Parameter | \-                          |

A researcher wanting to compare two datasets – for example, one
including only measurements `growthcleanr` flags for inclusion, and a
second which also includes carried forward values – could extract only
measurements where the flag is set to “`Include`” for the first set and
either “`Include`” or “`Exclude-Carried-Forward`” for the second set.
This flexibility gives each researcher control over how to incorporate
growthcleanr results into their own research.

## <a name="largedata"></a> Working with large data sets

The nature of the `growthcleanr` algorithm is repetitive. It runs many
checks over each individual subject’s measurements, in a specific
sequence, and revises its assessments as it goes. Because of this, it
can take some time to process even a modest dataset. For reference, the
`syngrowth` synthetic data example packaged with `growthcleanr` takes
2-3 minutes to process on a contemporary laptop. `growthcleanr` uses
optimized libraries (e.g. `data.table`) to go as fast as possible, but
there are limits given the repeated passes over data required.

If you have under one million records, using `growthcleanr` may just
necessitate planning to wait a little while the job processes. There are
several strategies to improve performance you might consider:

  - The built-in `parallel` option should take advantage of multiple CPU
    cores when your hardware allows.
  - Running `growthcleanr` on a machine with more CPU cores and more RAM
    can improve this further.
  - Verify that supporting libraries like `data.table` are [installed
    properly](#platform) for your machine
  - Anecdotally, Windows users have sometimes seen slower performance
    when compared with similarly-appointed hardware running Linux or
    macOS.

For very large datasets, on the order of millions of records or more,
cleaning data can take several hours or more. When running this kind of
job, the risk of failures due to external factors such as power outages
increases. The current structure of `growthcleanr` does not checkpoint
results as progress is made, so if a job has to run overnight and fails
in the middle, a researcher would likely need to start over from
scratch.

Because `growthcleanr` operates for the most part on individual subjects
one at a time, however, this issue might be mitigated by splitting the
input data into many small files, then running `growthcleanr` separately
on each file, with results re-combined at the end. The primary benefit
of this strategy would be the saving of the results of each file as it
completes, allowing a researcher whose job fails overnight to resume
processing on only the as yet incomplete input files. A secondary
benefit might be a slight improvement of use of available RAM.

Adopting this approach might require some custom code, and there are few
pitfalls to avoid. The following lays out a rough approach:

  - Extract `sd.median` values for the entire dataset as a whole using
    the `sdmedian.filename` option. This is critical, as `growthcleanr`
    will generate these values itself for an input dataset if they are
    not provided. If a large dataset is split into 1,000 smaller files
    to be cleaned separately, each of those 1,000 `growthcleanr` jobs
    needs to use the same `sd.median` values to recenter or the results
    will be inconsistent. Pass this data in to `cleangrowth()` using the
    parameter `sd.recenter`.
  - Split the data into many small files, but keep each individual
    subject’s measurements together in one file. If, for example, a
    subject has 12 measurements, all 12 should be in one and only one of
    the smaller files to maximize the longitudinal analysis of their
    values. The `splitinput()` function will perform this split safely,
    e.g.:

<!-- end list -->

``` r
count <- split_input(syngrowth, fname = "mydata", fdir = tempdir())
count
#> [1] 5
list.files(tempdir(), pattern = "mydata.*")
#> [1] "mydata.00000.csv" "mydata.00001.csv" "mydata.00002.csv" "mydata.00003.csv"
#> [5] "mydata.00004.csv" "mydata.00005.csv"
```

  - Write a standalone driver script to execute `growthcleanr` on each
    separate file, then write out its results when complete. An example
    is below.
  - Invoke the job using a tool like [GNU
    Parallel](https://www.gnu.org/software/parallel/), which can be run
    on all major platforms. An example invocation is below as well.
  - If a failure occurs, set aside both the completed inputs and their
    results. Then re-run the job on only the remaining input data.
  - When complete, re-combine the data into one file or as otherwise
    appropriate.

A standalone driver script might look like the following (note that it
uses the `argparse` library for processing options; this should make it
easy to read, but requires additional installation steps from what is
recommended above under [Installation](#installation)).

``` r
library(argparse)
library(data.table)
library(growthcleanr)

parser <- ArgumentParser(description = "Run cleangrowth() from script")
parser$add_argument("infile", 
  metavar = "INFILE", type = "character", nargs = 1, help = "input file"
)
parser$add_argument("outfile",
  metavar = "OUTFILE", type = "character", nargs = 1, help = "output file"
)
parser$add_argument("--sdrecenter",
  type = "character", nargs = 1, default = "", help = "sd.recenter data file"
)
parser$add_argument("--quietly",
  default = FALSE, action = "store_true", help = "Disable verbose output"
)
args <- parser$parse_args()

if (args$sdrecenter != "") {
  sdrecenter <- fread(args$sdrecenter)
} else {
  sdrecenter <- ""
}

df_in <- fread(args$infile)
df_out <- df_in[, 
  clean_value := cleangrowth(subjid, param, agedays, sex, measurement,
    sd.recenter = sdrecenter,
    quietly = args$quietly
  )
]
fwrite(df_out, args$outfile, row.names = FALSE)
```

A file like this could be saved as, for example, `gcdriver.R`, which can
then be invoked on a single input file using `Rscript`:

``` bash
Rscript gcdriver.R --quietly --sdrecenter sdmedians.csv mydata.csv mydata-cleaned.csv
```

If you have many small files saved, for example, as `mydata.00001.csv`,
`mydata.00002.csv`, etc., this driver can be invoked using `parallel`:

``` bash
ls mydata.?????.csv | parallel -j2 --eta \
  "Rscript gcdriver.R --quietly --sdrecenter sdmedians.csv {} {}-clean"
```

This lists your input files, and passes the filename list to parallel to
use in invoking the driver script, one file at a time, and to process as
parallel jobs as resources are available. Each run of the job then saves
the cleaned output with `-clean` appended to the input filename, and as
each completes, the next file on the list will be started until all are
complete. A few things to note:

  - The careful use of naming and wildcards when listing input files and
    naming output files can save accidental re-running of data from
    output files
  - The `-j2` option specifies running two jobs at once; `-j0` would use
    as many CPU cores as a machine has available. This might need to
    vary to match specific hardware.
  - The `--eta` option will report on progress.
  - The `--quietly` option on the driver script will make it easier to
    monitor progress with less verbose output coming from
    `growthcleanr`.
  - The `--sdrecenter` option on the driver script should be set to
    ensure each individual file is recentered using the same set for the
    entire input.
  - If multiple cores are available, the job should proceed with speedup
    roughly similar to what can be gained using the parallel batching
    feature built in to `growthcleanr`, with the main difference being
    the saving of intermediate output as each smaller file completes.

If your data to be cleaned is very large, it might help to store it
compressed, for example with `gzip` and its corresponding `.gz` filename
extension. `growthcleanr` can read in `.gz` input, but you might need to
install the `R.utils` package first. `R` will provide a message if this
is required.

## <a name="bmi"></a>Computing BMI percentiles and Z-scores

The CDC provides a [SAS macro for computing BMI percentiles and
Z-scores](https://www.cdc.gov/nccdphp/dnpao/growthcharts/resources/sas.htm);
the function `ext_bmiz()`, included in `growthcleanr`, provides an
equivalent feature. `ext_bmiz()` calculates the sigma (scale parameter
for the half-normal distribution, extended BMI percentile), extended
BMIz, and the CDC LMS z-scores for weight, height, and BMI. Note that
for BMIs ≤ 95th percentile of the CDC growth charts, the extended values
for BMI are equal to the LMS values. The extended values differ only for
children who have a BMI \> 95th percentile.

The function assumes a variable ‘sex’ (coded as 1=boys / 2=girls) and
variables for age in months, weight (kg), height (cm), and BMI
(weight/ht2). Please be careful with age - the units should be months
and use the most accurate information available (e.g., 23.4928 months.
The extended BMIz is the inverse cumulative distribution function (CDF)
of the extended BMI percentile. If the extended percentile is very close
to 100, the `qnorm` function in R produces an infinite value. The occurs
only if the extended BMI percentile is \> 99.99999999999999. This occurs
infrequently, such as a 48-month-old with a BMI \> 39, and it is likely
that these BMIs represent data entry errors. For these cases, extended
BMIz is set to 8.21, a value that is slightly greater than the largest
value that can be calculated.

Because `ext_bmiz()` performs cross-sectional analysis of BMI,
observation data must be in a wide format, i.e. with height and weight
information on the same row. This is distinct from `cleangrowth()`,
which performs longitudinal analysis on all observations for each
subject, presented in a long format with one observation per row. To
facilitate use of both functions, `growthcleanr` includes a utility
function to transform data used with `cleangrowth()` for use with
`ext_bmiz()`. It is optimized to move data directly from the output of
`cleangrowth()` into input for `ext_bmiz()`, but has options to support
independent use as well.

Using the `syngrowth` example dataset, to convert the data after it has
been cleaned by `cleangrowth()` for use with `ext_bmiz()`, use
`long_wide()`:

``` r
# Use the built-in utility function to convert the input observations to wide
# format for BMI calculation
cleaned_data_wide <- long_wide(cleaned_data)

# Compute Z-scores and percentiles
cleaned_data_bmi <- ext_bmiz(cleaned_data_wide)
```

Note that this assumes that `cleaned_data` has the same structure as
before:

``` r
names(cleaned_data)
#> [1] "id"          "subjid"      "sex"         "agedays"     "param"      
#> [6] "measurement" "clean_value"
```

The wide dataset `cleaned_data_wide` will include rows with aligned
height and weight measurements drawn from the observations in
`cleaned_data` marked by `cleangrowth()` for inclusion. As such, it will
be a shorter dataset (fewer rows) based on fewer observations.

``` r
dim(cleaned_data)
#> [1] 56703     7

dim(cleaned_data_wide)
#> [1] 17191    11

head(cleaned_data_wide)
#>     subjid   agey     age      bmi sex   wt wt_id     ht ht_id agedays seq_
#> 1:  775155 2.4339 29.2068 16.51604   1 13.1  1518  89.06  1511     889    1
#> 2:  775155 2.9322 35.1864 17.18042   1 14.7  1519  92.50  1512    1071    2
#> 3:  775155 3.4305 41.1660 16.53260   1 15.3  1520  96.20  1513    1253    3
#> 4:  775155 5.9603 71.5236 16.01739   1 20.2  1523 112.30  1517    2177    4
#> 5: 1340377 3.9288 47.1456 15.64862   2 15.9  7965 100.80  7950    1435    5
#> 6: 1340377 4.8871 58.6452 16.04127   2 18.4  7966 107.10  7951    1785    6
```

In this example, the subject identifiers previously marked as `subjid`
are now in the `id` column; individual identifiers for observations of a
single parameter are not present.

`long_wide()` can be called with name mapping parameters if your input
set uses different column names. For example, if `my_cleaned_data`
specifies age in days as `aged` and parameter type as `type`, specify
each, with quotes:

``` r
my_cleaned_data <- copy(cleaned_data)
setnames(my_cleaned_data, old = c("agedays", "param"), new = c("aged", "type"))
head(my_cleaned_data)
#>      id subjid sex aged     type measurement             clean_value
#> 1: 1510 775155   0  889 HEIGHTCM      84.900       Exclude-Duplicate
#> 2: 1511 775155   0  889 HEIGHTCM      89.060                 Include
#> 3: 1512 775155   0 1071 HEIGHTCM      92.500                 Include
#> 4: 1513 775155   0 1253 HEIGHTCM      96.200                 Include
#> 5: 1514 775155   0 1435 HEIGHTCM      96.200 Exclude-Carried-Forward
#> 6: 1515 775155   0 1435 HEIGHTCM      99.692                 Include

cleaned_data_wide <- long_wide(my_cleaned_data, agedays = "aged", param = "type")
cleaned_data_wide
#>            subjid    agey     agem      bmi sex   wt wt_id     ht ht_id agedays
#>     1:     775155  2.4339  29.2068 16.51604   1 13.1  1518  89.06  1511     889
#>     2:     775155  2.9322  35.1864 17.18042   1 14.7  1519  92.50  1512    1071
#>     3:     775155  3.4305  41.1660 16.53260   1 15.3  1520  96.20  1513    1253
#>     4:     775155  5.9603  71.5236 16.01739   1 20.2  1523 112.30  1517    2177
#>     5:    1340377  3.9288  47.1456 15.64862   2 15.9  7965 100.80  7950    1435
#>    ---                                                                         
#> 17187: 2146098000 10.0233 120.2796 15.88932   2 31.5 56627 140.80 56616    3661
#> 17188: 2146138598  2.4339  29.2068 15.47239   2 14.2 15552  95.80 15545     889
#> 17189: 2146138598  3.4305  41.1660 15.51762   2 16.8 15555 104.05 15548    1253
#> 17190: 2146138598  3.5838  43.0056 15.60233   2 17.3 15556 105.30 15549    1309
#> 17191: 2146138598  3.9288  47.1456 15.57547   2 18.1 15557 107.80 15550    1435
```

By default, `long_wide()` will only transform records flagged by
`cleangrowth()` for inclusion. To include more records, specify each
category assigned by `cleangrowth()` using the `inclusion_types` option.
For example, to include carried forward values along with included
records for the BMI calculation:

``` r
cleaned_data_wide_cf <- long_wide(
  long_df = cleaned_data,
  inclusion_types = c("Include", "Exclude-Carried-Forward")
)
cleaned_data_wide_cf
#>            subjid   agey    agem      bmi sex   wt wt_id      ht ht_id agedays
#>     1:     775155 2.4339 29.2068 16.51604   1 13.1  1518  89.060  1511     889
#>     2:     775155 2.9322 35.1864 17.18042   1 14.7  1519  92.500  1512    1071
#>     3:     775155 3.4305 41.1660 16.53260   1 15.3  1520  96.200  1513    1253
#>     4:     775155 3.9288 47.1456 16.53260   1 15.3  1521  96.200  1514    1435
#>     5:     775155 3.9288 47.1456 15.39469   1 15.3  1521  99.692  1515    1435
#>    ---                                                                        
#> 23132: 2146138598 2.4339 29.2068 15.47239   2 14.2 15552  95.800 15545     889
#> 23133: 2146138598 2.9322 35.1864 14.14337   2 14.2 15553 100.200 15546    1071
#> 23134: 2146138598 3.4305 41.1660 15.51762   2 16.8 15555 104.050 15548    1253
#> 23135: 2146138598 3.5838 43.0056 15.60233   2 17.3 15556 105.300 15549    1309
#> 23136: 2146138598 3.9288 47.1456 15.57547   2 18.1 15557 107.800 15550    1435
```

An additional option, `include_all`, set to `FALSE` by default, will
include all observations for transformation.

See `?long_wide` for full details.

With wide data in hand, output taken directly from `long_wide()` can be
passed to `ext_bmiz()`. Alternatively, you can provide a similarly
formatted data frame directly to `ext_bmiz()`.

Note that `ext_bmiz()` expects the `sex` variable to be coded as a
numeric value of `1` (male) or `2` (female). `long_wide()` will handle
this given the expected input values of `0` (male) or `1` (female).

If you are using different input data for `ext_bmiz()`, ensure the `sex`
variable in your data set is properly coded. To do so, use
`recode_sex()`. For example, if you have data in the PCORnet CDM format:

``` r
recode_sex(
  input_data = cleaned_data,
  sourcecol = "sex",
  sourcem = "M",
  sourcef = "F",
  targetm = 1L,
  targetf = 2L
)
#>           id     subjid sex agedays    param measurement
#>     1:  1510     775155   0     889 HEIGHTCM      84.900
#>     2:  1511     775155   0     889 HEIGHTCM      89.060
#>     3:  1512     775155   0    1071 HEIGHTCM      92.500
#>     4:  1513     775155   0    1253 HEIGHTCM      96.200
#>     5:  1514     775155   0    1435 HEIGHTCM      96.200
#>    ---                                                  
#> 56699: 15553 2146138598   1    1071 WEIGHTKG      14.200
#> 56700: 15554 2146138598   1    1071 WEIGHTKG      13.844
#> 56701: 15555 2146138598   1    1253 WEIGHTKG      16.800
#> 56702: 15556 2146138598   1    1309 WEIGHTKG      17.300
#> 56703: 15557 2146138598   1    1435 WEIGHTKG      18.100
#>                    clean_value sex_recoded
#>     1:       Exclude-Duplicate          NA
#>     2:                 Include          NA
#>     3:                 Include          NA
#>     4:                 Include          NA
#>     5: Exclude-Carried-Forward          NA
#>    ---                                    
#> 56699: Exclude-Carried-Forward          NA
#> 56700:       Exclude-Duplicate          NA
#> 56701:                 Include          NA
#> 56702:                 Include          NA
#> 56703:                 Include          NA
```

With data in wide format and with sex variables properly coded,
`ext_bmiz()` can be called:

``` r
cleaned_data_bmi <- ext_bmiz(cleaned_data_wide)
head(cleaned_data_bmi)
#>     subjid   agey     age      bmi sex   wt wt_id     ht ht_id agedays
#> 1:  775155 2.4339 29.2068 16.51604   1 13.1  1518  89.06  1511     889
#> 2:  775155 2.9322 35.1864 17.18042   1 14.7  1519  92.50  1512    1071
#> 3:  775155 3.4305 41.1660 16.53260   1 15.3  1520  96.20  1513    1253
#> 4:  775155 5.9603 71.5236 16.01739   1 20.2  1523 112.30  1517    2177
#> 5: 1340377 3.9288 47.1456 15.64862   2 15.9  7965 100.80  7950    1435
#> 6: 1340377 4.8871 58.6452 16.04127   2 18.4  7966 107.10  7951    1785
#>         bmiz     bmip        waz       wp         haz       hp      p95
#> 1: 0.1628453 56.46799 -0.1985513 42.13069 -0.36884065 35.61233 18.78685
#> 2: 0.8854337 81.20386  0.2999721 61.79008 -0.51324051 30.38915 18.31932
#> 3: 0.5828371 71.99985  0.1093221 54.35265 -0.50266522 30.75998 18.01381
#> 4: 0.4623661 67.80906 -0.1389032 44.47633 -0.56390741 28.64086 18.36444
#> 5: 0.2541728 60.03190  0.1210300 54.81663  0.12282474 54.88771 18.02950
#> 6: 0.6196201 73.22461  0.2780873 60.95273  0.04849898 51.93407 18.19673
#>         p97   bmip95  mod_bmiz     mod_waz     mod_haz    sigma ext_bmip
#> 1: 19.23180 87.91274 0.1306365 -0.22517132 -0.37732507 1.583547 56.46799
#> 2: 18.71228 93.78308 0.7875695  0.25504032 -0.53302684 1.818131 81.20386
#> 3: 18.38801 91.77740 0.5066643  0.08834330 -0.51529883 2.048196 71.99985
#> 4: 19.04130 87.21960 0.3150702 -0.17522959 -0.56091593 3.146493 67.80906
#> 5: 18.60078 86.79455 0.1773321  0.08808784  0.11920984 2.274792 60.03190
#> 6: 18.90966 88.15470 0.4237897  0.19696572  0.04639902 2.621219 73.22461
#>     ext_bmiz sev_obese obese
#> 1: 0.1628453         0     0
#> 2: 0.8854337         0     0
#> 3: 0.5828371         0     0
#> 4: 0.4623661         0     0
#> 5: 0.2541728         0     0
#> 6: 0.6196201         0     0
```

The output columns include:

| variable                     | description                                                                                                                                                                       |
| ---------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| waz, wp                      | LMS Weight-for-sex/age z-score and percentile                                                                                                                                     |
| haz, hp                      | LMS Height-for-sex/age z-score and percentile                                                                                                                                     |
| bmiz, bmip                   | LMS BMI-for-sex/age z-score and percentile                                                                                                                                        |
| mod\_waz, mod\_haz, mod\_bmi | Modified z-scores for identifying outliers (see the information in the [CDC SAS growth charts program website](https://www.cdc.gov/nccdphp/dnpao/growthcharts/resources/sas.htm)) |
| bmip95                       | BMI expressed as a percentage of the 95th percentile. A value ≥ 120 is widely used as the cut point for severe obesity.                                                           |
| sigma                        | Scale parameter of the half-normal distribution                                                                                                                                   |
| ext\_bmip                    | Extended BMI percentile                                                                                                                                                           |
| ext\_bmiz                    | Extended BMI z-score                                                                                                                                                              |

For convenience, these labels are available on the output of
`ext_bmiz()`, e.g., when viewed in RStudio with
`View(cleaned_data_bmi)`, or on the console:

``` r
labelled::var_label(cleaned_data_bmi)
#> $subjid
#> NULL
#> 
#> $agey
#> NULL
#> 
#> $age
#> [1] "Age (months)"
#> 
#> $bmi
#> [1] "BMI"
#> 
#> $sex
#> NULL
#> 
#> $wt
#> [1] "Weight (kg)"
#> 
#> $wt_id
#> NULL
#> 
#> $ht
#> [1] "Height (cm)"
#> 
#> $ht_id
#> NULL
#> 
#> $agedays
#> NULL
#> 
#> $bmiz
#> [1] "LMS BMI-for-sex/age z-score"
#> 
#> $bmip
#> [1] "LMS BMI-for-sex/age percentile"
#> 
#> $waz
#> [1] "LMS Weight-for-sex/age z-score"
#> 
#> $wp
#> [1] "LMS Weight-for-sex/age percentile"
#> 
#> $haz
#> [1] "LMS Height-for-sex/age z-score"
#> 
#> $hp
#> [1] "LMS Height-for-sex/age percentile"
#> 
#> $p95
#> [1] "95th percentile of BMI in growth charts"
#> 
#> $p97
#> [1] "97th percentile of BMI in growth charts"
#> 
#> $bmip95
#> [1] "BMI as a percentage of the 95th percentile"
#> 
#> $mod_bmiz
#> [1] "Modified BMI-for-age z-score"
#> 
#> $mod_waz
#> [1] "Modified Weight-for-age z-score"
#> 
#> $mod_haz
#> [1] "Modified Height-for-age z-score"
#> 
#> $sigma
#> [1] "Scale parameter for half-normal distribution"
#> 
#> $ext_bmip
#> [1] "Extended BMI percentile"
#> 
#> $ext_bmiz
#> [1] "Extended BMI z-score"
#> 
#> $sev_obese
#> [1] "BMI >= 120% of 95th percentile (0/1)"
#> 
#> $obese
#> [1] "BMI >= 95th percentile (0/1)"
```

Like `long_wide()`, `ext_bmiz()` also includes options for mapping
alternate column names, for age, weight, height, and BMI. The default
column names are the same as the output from `long_wide()` for
convenience. If you have different column names, specify the column
names without quotes. For example, for a dataset using “heightcm” and
“weightkg” instead of “ht” and “wt”:

``` r
my_cleaned_data_wide <- copy(cleaned_data_wide)
setnames(my_cleaned_data_wide, old = c("ht", "wt"), new = c("heightcm", "weightkg"))
head(my_cleaned_data_wide)
#>     subjid   agey     age      bmi sex weightkg wt_id heightcm ht_id agedays
#> 1:  775155 2.4339 29.2068 16.51604   1     13.1  1518    89.06  1511     889
#> 2:  775155 2.9322 35.1864 17.18042   1     14.7  1519    92.50  1512    1071
#> 3:  775155 3.4305 41.1660 16.53260   1     15.3  1520    96.20  1513    1253
#> 4:  775155 5.9603 71.5236 16.01739   1     20.2  1523   112.30  1517    2177
#> 5: 1340377 3.9288 47.1456 15.64862   2     15.9  7965   100.80  7950    1435
#> 6: 1340377 4.8871 58.6452 16.04127   2     18.4  7966   107.10  7951    1785
#>    seq_
#> 1:    1
#> 2:    2
#> 3:    3
#> 4:    4
#> 5:    5
#> 6:    6

ext_bmiz(my_cleaned_data_wide, age = "age", ht = "heightcm", wt = "weightkg")
#>            subjid    agey      age      bmi sex   wt wt_id     ht ht_id agedays
#>     1:     775155  2.4339  29.2068 16.51604   1 13.1  1518  89.06  1511     889
#>     2:     775155  2.9322  35.1864 17.18042   1 14.7  1519  92.50  1512    1071
#>     3:     775155  3.4305  41.1660 16.53260   1 15.3  1520  96.20  1513    1253
#>     4:     775155  5.9603  71.5236 16.01739   1 20.2  1523 112.30  1517    2177
#>     5:    1340377  3.9288  47.1456 15.64862   2 15.9  7965 100.80  7950    1435
#>    ---                                                                         
#> 17187: 2146098000 10.0233 120.2796 15.88932   2 31.5 56627 140.80 56616    3661
#> 17188: 2146138598  2.4339  29.2068 15.47239   2 14.2 15552  95.80 15545     889
#> 17189: 2146138598  3.4305  41.1660 15.51762   2 16.8 15555 104.05 15548    1253
#> 17190: 2146138598  3.5838  43.0056 15.60233   2 17.3 15556 105.30 15549    1309
#> 17191: 2146138598  3.9288  47.1456 15.57547   2 18.1 15557 107.80 15550    1435
#>                bmiz     bmip        waz       wp        haz       hp      p95
#>     1:  0.162845256 56.46799 -0.1985513 42.13069 -0.3688406 35.61233 18.78685
#>     2:  0.885433660 81.20386  0.2999721 61.79008 -0.5132405 30.38915 18.31932
#>     3:  0.582837141 71.99985  0.1093221 54.35265 -0.5026652 30.75998 18.01381
#>     4:  0.462366057 67.80906 -0.1389032 44.47633 -0.5639074 28.64086 18.36444
#>     5:  0.254172837 60.03190  0.1210300 54.81663  0.1228247 54.88771 18.02950
#>    ---                                                                       
#> 17187: -0.458129542 32.34297 -0.2404025 40.50091  0.3969069 65.42820 22.96109
#> 17188: -0.486411301 31.33378  0.8546829 80.36366  1.7235811 95.76083 18.66623
#> 17189:  0.005468304 50.21815  1.0128840 84.44422  1.7023375 95.56539 18.10574
#> 17190:  0.125610747 54.99800  1.0593271 85.52746  1.7220601 95.74707 18.06912
#> 17191:  0.197525686 57.82919  1.0256103 84.74624  1.6956671 95.50255 18.02950
#>             p97   bmip95     mod_bmiz     mod_waz    mod_haz    sigma ext_bmip
#>     1: 19.23180 87.91274  0.130636512 -0.22517132 -0.3773251 1.583547 56.46799
#>     2: 18.71228 93.78308  0.787569526  0.25504032 -0.5330268 1.818131 81.20386
#>     3: 18.38801 91.77740  0.506664311  0.08834330 -0.5152988 2.048196 71.99985
#>     4: 19.04130 87.21960  0.315070178 -0.17522959 -0.5609159 3.146493 67.80906
#>     5: 18.60078 86.79455  0.177332068  0.08808784  0.1192098 2.274792 60.03190
#>    ---                                                                        
#> 17187: 24.57902 69.20105 -0.590757876 -0.31445677  0.3859532 4.443536 32.34297
#> 17188: 19.12367 82.88974 -0.554892094  0.74219573  1.7224761 1.730347 31.33378
#> 17189: 18.62133 85.70552  0.003901076  0.86190237  1.6961075 2.093856 50.21815
#> 17190: 18.60021 86.34799  0.089438272  0.90273897  1.7156579 2.149579 54.99800
#> 17191: 18.60078 86.38881  0.136551089  0.85844390  1.6874890 2.274792 57.82919
#>            ext_bmiz sev_obese obese
#>     1:  0.162845256         0     0
#>     2:  0.885433660         0     0
#>     3:  0.582837141         0     0
#>     4:  0.462366057         0     0
#>     5:  0.254172837         0     0
#>    ---                             
#> 17187: -0.458129542         0     0
#> 17188: -0.486411301         0     0
#> 17189:  0.005468304         0     0
#> 17190:  0.125610747         0     0
#> 17191:  0.197525686         0     0
```

For `ext_bmiz()`, use the most precise age in months available. If an
input dataset only has age in months as integer values, by default
`ext_bmiz()` will automatically convert these to double values and add
`0.5` to account for the distribution of actual ages over the range of
days within a month. This is enabled with the option
`adjust.integer.age`, set to `TRUE` by default. Specify `FALSE` to
disable.

``` r
ext_bmiz(cleaned_data_wide, age = "age", adjust.integer.age = FALSE)
#>            subjid    agey      age      bmi sex   wt wt_id     ht ht_id agedays
#>     1:     775155  2.4339  29.2068 16.51604   1 13.1  1518  89.06  1511     889
#>     2:     775155  2.9322  35.1864 17.18042   1 14.7  1519  92.50  1512    1071
#>     3:     775155  3.4305  41.1660 16.53260   1 15.3  1520  96.20  1513    1253
#>     4:     775155  5.9603  71.5236 16.01739   1 20.2  1523 112.30  1517    2177
#>     5:    1340377  3.9288  47.1456 15.64862   2 15.9  7965 100.80  7950    1435
#>    ---                                                                         
#> 17187: 2146098000 10.0233 120.2796 15.88932   2 31.5 56627 140.80 56616    3661
#> 17188: 2146138598  2.4339  29.2068 15.47239   2 14.2 15552  95.80 15545     889
#> 17189: 2146138598  3.4305  41.1660 15.51762   2 16.8 15555 104.05 15548    1253
#> 17190: 2146138598  3.5838  43.0056 15.60233   2 17.3 15556 105.30 15549    1309
#> 17191: 2146138598  3.9288  47.1456 15.57547   2 18.1 15557 107.80 15550    1435
#>                bmiz     bmip        waz       wp        haz       hp      p95
#>     1:  0.162845256 56.46799 -0.1985513 42.13069 -0.3688406 35.61233 18.78685
#>     2:  0.885433660 81.20386  0.2999721 61.79008 -0.5132405 30.38915 18.31932
#>     3:  0.582837141 71.99985  0.1093221 54.35265 -0.5026652 30.75998 18.01381
#>     4:  0.462366057 67.80906 -0.1389032 44.47633 -0.5639074 28.64086 18.36444
#>     5:  0.254172837 60.03190  0.1210300 54.81663  0.1228247 54.88771 18.02950
#>    ---                                                                       
#> 17187: -0.458129542 32.34297 -0.2404025 40.50091  0.3969069 65.42820 22.96109
#> 17188: -0.486411301 31.33378  0.8546829 80.36366  1.7235811 95.76083 18.66623
#> 17189:  0.005468304 50.21815  1.0128840 84.44422  1.7023375 95.56539 18.10574
#> 17190:  0.125610747 54.99800  1.0593271 85.52746  1.7220601 95.74707 18.06912
#> 17191:  0.197525686 57.82919  1.0256103 84.74624  1.6956671 95.50255 18.02950
#>             p97   bmip95     mod_bmiz     mod_waz    mod_haz    sigma ext_bmip
#>     1: 19.23180 87.91274  0.130636512 -0.22517132 -0.3773251 1.583547 56.46799
#>     2: 18.71228 93.78308  0.787569526  0.25504032 -0.5330268 1.818131 81.20386
#>     3: 18.38801 91.77740  0.506664311  0.08834330 -0.5152988 2.048196 71.99985
#>     4: 19.04130 87.21960  0.315070178 -0.17522959 -0.5609159 3.146493 67.80906
#>     5: 18.60078 86.79455  0.177332068  0.08808784  0.1192098 2.274792 60.03190
#>    ---                                                                        
#> 17187: 24.57902 69.20105 -0.590757876 -0.31445677  0.3859532 4.443536 32.34297
#> 17188: 19.12367 82.88974 -0.554892094  0.74219573  1.7224761 1.730347 31.33378
#> 17189: 18.62133 85.70552  0.003901076  0.86190237  1.6961075 2.093856 50.21815
#> 17190: 18.60021 86.34799  0.089438272  0.90273897  1.7156579 2.149579 54.99800
#> 17191: 18.60078 86.38881  0.136551089  0.85844390  1.6874890 2.274792 57.82919
#>            ext_bmiz sev_obese obese
#>     1:  0.162845256         0     0
#>     2:  0.885433660         0     0
#>     3:  0.582837141         0     0
#>     4:  0.462366057         0     0
#>     5:  0.254172837         0     0
#>    ---                             
#> 17187: -0.458129542         0     0
#> 17188: -0.486411301         0     0
#> 17189:  0.005468304         0     0
#> 17190:  0.125610747         0     0
#> 17191:  0.197525686         0     0
```

`ext_bmiz()` uses reference data provided by the CDC, included in the
`growthcleanr` package as `inst/extdata/CDCref_d.csv`. This file is
automatically loaded and used by default. If you are working with a
different reference dataset or developing the `growthcleanr` package,
specify an alternate path to this file with `ref.data.path`, as for
`cleangrowth()`.

## <a name="related"></a>Related tools

The CDC provides a [SAS Program for the 2000 CDC Growth
Charts](https://www.cdc.gov/nccdphp/dnpao/growthcharts/resources/sas.htm)
which can also be used to identify biologically implausible values using
a different approach, as also implemented for `growthcleanr` in the
function `ext_bmiz()`, described above under [Computing BMI percentiles
and Z-scores](#bmi).

[GrowthViz](https://github.com/mitre/GrowthViz) provides insights into
how `growthcleanr` assesses data, packaged in a Jupyter notebook. It
ships with the same `syngrowth` synthetic example dataset as
`growthcleanr`, with results included.
