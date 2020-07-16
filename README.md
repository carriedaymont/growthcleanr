# growthcleanr

R package for cleaning data from Electronic Health Record systems, focused on
cleaning height and weight measurements.

<a name="cite"></a>
This package implements the
[Daymont et al. algorithm](https://academic.oup.com/jamia/article/24/6/1080/3767271),
as specified in Supplemental File 3 within the
[Supplementary Material](https://academic.oup.com/jamia/article/24/6/1080/3767271#97610899)
published with that paper.

> Carrie Daymont, Michelle E Ross, A Russell Localio, Alexander G Fiks, Richard
> C Wasserman, Robert W Grundmeier, Automated identification of implausible
> values in growth data from pediatric electronic health records, Journal of the
> American Medical Informatics Association, Volume 24, Issue 6, November 2017,
> Pages 1080–1087, https://doi.org/10.1093/jamia/ocx037

## Summary

The `growthcleanr` package processes data prepared in a specific format to
identify biologically implausible height and weight measurements. It bases these
evaluations on techniques which use patient-specific longitudinal analysis and
variations from published growth trajectory charts. These techniques are
performed in a specific order which refines and improves results throughout the
process.

Results from `growthcleanr` include a flag for each measurement indicating
whether it is to be included or excluded based on plausibility, with a variety
of specific types of exclusions identified distinctly. These flags can be
analyzed further by researchers studying growth data to determine which
measurements to include or exclude in their own studies. No values are deleted
or otherwise removed; each is only flagged in a new column.

To start running `growthcleanr`, an R installation with a variety of additional
packages is required, as is a growth measurement dataset prepared for use in
`growthcleanr`.

The rest of this documentation includes:

* [Quickstart](#quickstart) - a brief tour of using growthcleanr, including data
  preparation
* [Installation](#installation) - including required dependencies and
  OS-specific notes
* [Usage](#usage) - examples of cleaning data, multiple options, example data
* [Configuration options](#configuration) - changing growthcleanr operational
  settings
* [Understanding output](#output) - the exclusion types growthcleanr identifies
* [Working with large datasets](#largedata) - notes and suggestions
* [Related tools](#related) - additional software complementary to growthcleanr

## <a name="quickstart"></a>Quickstart

### <a name="setup"></a>R setup

An R environment with the following packages installed is required to run
growthcleanr:

* `devtools`
* `dplyr`
* `data.table`
* `foreach`
* `doParallel`
* `Hmisc`

The `growthcleanr` package itself must also be installed. Further installation
details and notes can be found under [Installation](#installation) below.

### Data preparation

Once the `growthcleanr` package and its dependencies are installed, the minimal
data vectors needed are:

1. `subjid` - a unique identifier for each subject; if you have long (64-bit)
   integer `subjid` values, add the `bit64` R package to your install list to
   ensure R will process these correctly
2. `param` - a string designating the measurement being '`HEIGHTCM`',
   '`LENGTHCM`', or '`WEIGHTKG`'
3. `agedays` - a numeric value for age in days at time of measurement; will be
   rounded down (floor) to nearest integer if needed (note that there are
   30.4375 days in one month; see also the additional
   [notes on calculating `agemos`](https://www.cdc.gov/nccdphp/dnpao/growthcharts/resources/sas.htm)
   from whole numbers of months
3. `sex` - a numeric value of `0` (male) or `1` (female)
4. `measurement` - numeric value of measurement in cm for height or kg for
   weight

**Note**: the dataset should be sorted by `subjid`, `param`, and `agedays`, as
demonstrated in the Example below.

### <a name="example"></a>Example

Starting with a dataset such as the following:

| subjid | param | agedays | sex | measurement |
| - | - | - | - | - |
| 1 | HEIGHTCM | 2790 | 0 | 118.5 |
| 1 | HEIGHTCM | 3677 | 0 | 148.22 |
| 1 | WEIGHTKG | 2790 | 0 | 23.8 |
| 1 | WEIGHTKG | 3677 | 0 | 38.41 |
| 2 | HEIGHTCM | 2112 | 1 | 118.2 |
| 2 | HEIGHTCM | 2410 | 1 | 123.9 |
| 2 | HEIGHTCM | 2708 | 1 | 128.03 |
| 2 | HEIGHTCM | 3029 | 1 | 135.3 |
| 2 | WEIGHTKG | 2112 | 1 | 25.84 |
| 2 | WEIGHTKG | 2410 | 1 | 28.61 |
| 2 | WEIGHTKG | 2708 | 1 | 28.61 |
| 2 | WEIGHTKG | 3029 | 1 | 34.5 |
| ... | ... | ... | ... | ... |

Note that in this Example, we see data from two subjects, sorted as described
above, with two pairs of height/weight measurements for the first patient and
four pairs for the second. Most of these values pass the "eyeball test", but the
third weight measurement for `subjid` 2 looks like it might be a carried forward
value, given that their same ageday height has changed by several centimeters
and their third height value is exactly the same as their second. This is an
example of an implausible value and exclusion type that `growthcleanr` might
identify.

The `cleangrowth()` function provided by `growthcleanr` will return a vector
that specifies, for each measurement, whether to:

* "`Include`" the measurement,
* Conclude that the measurement is "`Missing`",
* Or to "`Exclude`" the measurement (numerous variations of `Exclude` exist, as
  specified below in Understanding output).

For a data.frame object `source_data` containing growth data:

```R
library(growthcleanr)

# prepare data as a data.table
data <- as.data.table(source_data)

# set the data.table key for better indexing
setkey(data, subjid, param, agedays)

# generate new exclusion flag field using function
cleaned_data <- data[, clean_value:=cleangrowth(subjid, param, agedays, sex, measurement)]

# extract data limited only to values flagged for inclusion:
only_included_data <- cleaned_data[clean_vaue=='Include']
```

If our Example dataset above were named `source_data`, examining `cleaned_data`
would show:

| subjid | param | agedays | sex | measurement | clean_value |
| - | - | - | - | - | - |
| 1 | HEIGHTCM | 2790 | 0 | 118.5 | Include |
| 1 | HEIGHTCM | 3677 | 0 | 148.22 | Include |
| 1 | WEIGHTKG | 2790 | 0 | 23.8 | Include |
| 1 | WEIGHTKG | 3677 | 0 | 38.41 | Include |
| 2 | HEIGHTCM | 2112 | 1 | 118.2 | Include |
| 2 | HEIGHTCM | 2410 | 1 | 123.9 | Include |
| 2 | HEIGHTCM | 2708 | 1 | 128.03 | Include |
| 2 | HEIGHTCM | 3029 | 1 | 135.3 | Include |
| 2 | WEIGHTKG | 2112 | 1 | 25.84 | Include |
| 2 | WEIGHTKG | 2410 | 1 | 28.61 | Include |
| 2 | WEIGHTKG | 2708 | 1 | 28.61 | Exclude-Carried-Forward |
| 2 | WEIGHTKG | 3029 | 1 | 34.5 | Include |
| ... | ... | ... | ... | ... | ... |

In this sample output, we see that `growthcleanr` has indeed flagged the third
weight measurement for `subjid 2` as an apparent carried forward value. The
measurement remains in the dataset, and it is up to the researcher to determine
which exclusions to accept directly or investigate further.

## <a name="installation"></a>Installation

`growthcleanr has been developed and tested using R version 3 (most recently
3.6+). It should work using R on Windows, macOS, or Unix/Linux, although there
are some additional [platform-specific notes](#platform) you may wish to review.

In an up-to-date R environment such as RStudio, first install the dependencies:

```R
# Required packages for growthcleanr operation
install.packages(c("devtools", "dplyr", "data.table", "foreach", "doParallel", "Hmisc"))
# Optional additional packages referenced in this document
install.packages("argparse")
install.packages("bit64")
```

Note that this might take some time.

You can install the `growthcleanr` package directly from GitHub using `devtools`
in the R console with:

```R
devtools::install_github("carriedaymont/growthcleanr")
```

### Source-level install for developers

If you want to work with and potentially change the `growthcleanr` code itself,
you can download or clone the `growthcleanr` source code and then install it from
source. To clone the source using `git`:

```bash
git clone https://github.com/carriedaymont/growthcleanr.git
```

Either way, once you have the `growthcleanr` package source, open an R session
from the `growthcleanr` base directory. Then install growthcleanr using the R
`devtools` package:

```R
devtools::install(".")
```

You can also install the package from an installation file if one is obtained.

### <a name="platform"></a>Platform-specific notes

#### Windows

We have observed an issue with using `growthcleanr` with Windows in a large
agency where some R packages are installed by an administrator for a user who
may not have write access permissions on the package folder. Similar issues may
occur with some networked drives. This caused problems with `growthcleanr`'s
parallel processing option. If possible, install R and its packages in locations
hosted on the same local machine and folder(s) for which the primary user has
write permissions. These steps should help to avoid this problem.

Some users have reported that `growthcleanr` runs more slowly on Windows
compared with Linux or macOS.

#### macOS

`growthcleanr` uses the
[`data.table`](https://rdatatable.gitlab.io/data.table/) package for R
extensively. `data.table` provides a faster version of R's data frames, and is
used to improve `growthcleanr` performance. Typically `data.table` installs in a
manner that will be able to take advantage of multiple threads. You will know it
worked successfully if, when you load the `data.table` library in R, you see
something like the following:

```R
> library(data.table)
data.table 1.12.2 using 2 threads (see ?getDTthreads).  Latest news: r-datatable.com
```

That data.table reports "using 2 threads" indicates that installation has
succeeded. If the message reports using only one thread, see the advice under
the ["OpenMP enabled compiler for Mac"
instructions](https://github.com/Rdatatable/data.table/wiki/Installation#openmp-enabled-compiler-for-mac)
to re-install `data.table`.

Users have reported errors running multiple batches with `parallel = T` from
within RStudio. If this happens, the problem may be resolved by running from
RGui or from the command line using `Rscript`. An example standalone script that
may be used for this purpose is documented below under [Working with large data
sets](#largedata).

#### Docker

This package includes a `Dockerfile` that enables easy installation of R and
`growthcleanr` on a machine with [Docker](https://www.docker.com/) installed. It
requires an up-to-date Docker install, and a few command-line steps, but can save time
over installing R and `growthcleanr`'s dependencies manually.

To install and run `growthcleanr` using Docker:

```bash
% docker run -it mitre/growthcleanr:latest R
```

The first time this command is run, it might take a few minutes to download several
necessary components, but this should be fully automated. If successful, you should see
an R prompt, from which you can use `growthcleanr` as described below.

This R environment is virtualized inside Docker, however, and isolated from your local
machine. Because of this, you will need to map a local folder on your computer into the
Docker environment to work with your own data. For example, if your data is in
`/Users/exampleuser/analysis`, specify a folder mapping using the added `-v` step below:

```bash
% docker run -it -v /Users/exampleusers/analysis:/usr/src/app mitre/growthcleanr:latest R
```

Inside the Docker environment's R prompt, if you issue a command like `list.files()`,
you should see a list of the same files from that folder. You can now open and read data
files, run `cleangrowth()` and other analyses, and write result files to that same
directory.

Exit the Docker R environment with `quit()` as you normally would. Any new files you
saved will appear in the folder you mapped.

## <a name="usage"></a>Usage

With the `growthcleanr` package installed [as described above](#installation),
it can be loaded as a package:

```R
> library(growthcleanr)
Loading required package: data.table
data.table 1.12.8 using 2 threads (see ?getDTthreads).  Latest news: r-datatable.com
Loading required package: plyr
```

For convenience, `growthcleanr` ships with an example synthetic dataset created
using [Synthea](https://synthetichealth.github.io/synthea/), with support for
simulated pediatric growth measurement errors based on the protocol included as
supplementary material to the [Daymont et al. paper](#cite). This dataset is
automatically loaded with the `growthcleanr` library, and is called `syngrowth`.

```R
> dim(syngrowth)
[1] 56703     6
> head(syngrowth)
  id     subjid sex agedays    param measurement
1  1 1575544570   1    2177 HEIGHTCM      104.40
2  2 1575544570   1    2548 HEIGHTCM      117.70
3  3 1575544570   1    2548 HEIGHTCM      117.74
4  4 1575544570   1    2919 HEIGHTCM      117.70
5  5 1575544570   1    3290 HEIGHTCM      129.20
6  6 1575544570   1    3661 HEIGHTCM      133.90
```

It can be used as in the example above. Note that processing the example data
with `cleangrowth()` will likely take a few minutes to complete.

```R
> library(dplyr)
> data <- as.data.table(syngrowth)
> setkey(data, subjid, param, agedays)
> cleaned_data <- data[, clean_value:=cleangrowth(subjid, param, agedays, sex, measurement)]
> head(cleaned_data)
     id subjid sex agedays    param measurement             clean_value
1: 1510 775155   0     889 HEIGHTCM      84.900       Exclude-Duplicate
2: 1511 775155   0     889 HEIGHTCM      89.060                 Include
3: 1512 775155   0    1071 HEIGHTCM      92.500                 Include
4: 1513 775155   0    1253 HEIGHTCM      96.200                 Include
5: 1514 775155   0    1435 HEIGHTCM      96.200 Exclude-Carried-Forward
6: 1515 775155   0    1435 HEIGHTCM      99.692                 Include
> cleaned_data %>% group_by(clean_value) %>% tally(sort=TRUE)
# A tibble: 14 x 2
   clean_value                   n
   <ord>                     <int>
 1 Include                   38875
 2 Exclude-Duplicate         10546
 3 Exclude-Carried-Forward    6694
 4 Exclude-SD-Cutoff           168
 5 Exclude-EWMA-8              135
 6 Exclude-EWMA-Extreme         95
 7 Exclude-EWMA-9               93
 8 Exclude-Min-Height-Change    65
 9 Swapped-Measurements         16
10 Exclude-Too-Many-Errors       6
11 Exclude-EWMA-11               5
12 Exclude-EWMA-12               2
13 Exclude-Pair-Delta-18         2
14 Exclude-Max-Height-Change     1
```

If you are able to run these steps and see a similar result, you have the
`growthcleanr` package installed correctly. The resulting `cleaned_data` can be
reviewed, subsetted, and compared in more detail using all the tools R provides.

For complete information about the options that can be set on the
`cleangrowth()` function, see [Configuration](#configuration). Below are a
few additional examples:

This example shows three configuration options in use:

```R
cleaned_data <- data[, clean_value_both:=
                       cleangrowth(subjid, param, agedays, sex, measurement,
                                   lt3.exclude.mode = 'flag.both',
                                   ref.data.path = "inst/extdata/",
                                   quietly=F
                                   )]
```

  * `lt3.exclude.mode = "flag.both"` will exclude both measurements for a
    subject if they only have two unexcluded measurements of one parameter type
    with at least one implausible value and no same age-day measurements of the
    other parameter.

  * `ref.data.path = "inst/extdata"` shouldn't be necessary if you have the
    `growthcleanr` package installed, but if you are running it from its
    source directly you may need to specify its full path.

  * `quietly = F` enables verbose output, marking the progress of the algorithm
    through its many processing steps. This can be very helpful while testing.

This example shows built-in options for processing data in batches, which can
speed the process while working with large data sets:

```R
cleaned_data <- data[, clean_value_both:=
                       cleangrowth(subjid, param, agedays, sex, measurement,
                                   parallel=T,
                                   num.batches=4,
                                   log.path="logs"
                                   )]
```

The best `num.batches` value for your environment may vary depending on the
computing resources you have available. If you do not specify `num.batches`,
`growthcleanr` will estimate a batch count based on R functions for checking the
system hardware. If you are working with large datasets, you may need to
experiment with these options to determine the best settings for your needs.
Note that there are additional notes on
[working with large datasets](#largedata) below.

The default value of `log.path` is `"."`, the current working directory.
`growthcleanr` will write batch-specific log files to the `log.path` directory,
and will create the directory if necessary.

Note that if you run `growthcleanr` in parallel with multiple batches and see
warning errors such as the following, they can be ignored.

```R
Warning messages:
1: <anonymous>: ... may be used in an incorrect context: ‘.fun(piece, ...)’
2: <anonymous>: ... may be used in an incorrect context: ‘.fun(piece, ...)’
```

This set of warnings is coming from one of `growthcleanr`'s dependencies, and
does not indicate either failure or improper execution.

## <a name="configuration"></a>Configuration options

`growthcleanr` offers several options for configuration, with default values set
to address common cases. You may wish to experiment with the settings to
discover which work best for your research and your dataset.

All of the following options may be set as additional parameters in the call to
`cleangrowth()`.

### Algorithm-related options

The following options change the behavior of the growthcleanr algorithm.

* `recover.unit.error` - default `FALSE`; when `FALSE`, measurements identified
  as unit errors (e.g., apparent height values in inches instead of centimeters)
  will be flagged but not corrected, when `TRUE` these values will be corrected
  and included as valid measurements for cleaning.

* `sd.extreme` - default `25`; a very extreme value check on modified
  (recentered) Z-scores used as a first-pass elimination of clearly implausable
  values, often due to misplaced decimals.
* `z.extreme` - default `25`; similar usage as `sd.extreme`, for absolute
  Z-scores.

**NOTE**: many different steps in the `growthcleanr` algorithm use highly
refined techniques to identify implausible values that require more subtlety to
detect. These default values for `sd.extreme` and `z.extreme` are set very high
by design to eliminate completely implausible values early in the process; these
extreme values require none of `growthcleanr`'s additional, more subtle
approaches to be identified as exclusions. If `sd.extreme` and `z.extreme` are
configured to be lower values, measurements eliminated at the early step which
checks against these extreme values will not be further considered using later
techniques.

* `include.carryforward` - default `FALSE`; if set to `TRUE`, `growthcleanr`
  will skip algorithm step 9, which identifies carried forward measurements, and
  will not flag these values for exclusion.

* `ewma.exp` - default `-1.5`; the exponent used for weighting measurements when
  calculating exponentially weighted moving average (EWMA). This exponent should
  be negative to weight growth measurements closer to the measurement being
  evaluated more strongly. Exponents that are further from zero (e.g., `-3`)
  will increase the relative influence of measurements close in time to the
  measurement being evaluated compared to using the default exponent.

* `ref.data.path` - defaults to using CDC reference data from year 2000; supply
  a file path to use alternate reference data. Note that when running from an
  installed `growthcleanr` package (e.g. having called `library(growthcleanr)`),
  this path does not need to be specified. Developers testing the source code
  directly from the source directory will need to specify this as well.

* `error.load.mincount` - default `2`; minimum count of exclusions on parameter
  for one subject before considering excluding all measurements.

* `error.load.threshold` - default `0.5`; threshold to exceed for percentage
  of excluded measurement count relative to count of included other measurement
  (e.g., if 3 of 5 WTs are excluded, and 5 corresponding HTs are included, this
  exceeds the 0.5 threshold and the other two WTs will be excluded).

* `lt3.exclude.mode` - default `default`; determines type of exclusion procedure
  to use for 1 or 2 measurements of one type without matching same ageday
  measurements for the other parameter. Options include:
  * `default` - standard growthcleanr approach
  * `flag.both` - in case of two measurements with at least one beyond
    thresholds, flag both instead of one (as in default)

* `sd.recenter` - defaults to NA; data frame or table w/median SD-scores per day
  of life by gender and parameter. Columns must include param, sex, agedays, and
  sd.median (referred to elsewhere as "modified Z-score"). By default, median
  values will be calculated using growth data to be cleaned.

### Operational options

The following options change the execution of the program overall, with no
effect on the algorithm itself.

* `parallel` - default `FALSE`; set to `TRUE` to run `growthcleanr` in parallel.
  Running in parallel will split the input data into batches (while keeping all
  records for each subject together) and process each batch on a different
  processor/core to maximize throughput. Recommended for large datasets with more
  than 100K rows.

* `num.batches` - specifies the number of batches to run in parallel; only
  applies if `parallel` is set to `TRUE`. Defaults to the number of workers
  returned by the `getDoParWorkers()` function in the `foreach` package. Note
  that processing in parallel may affect overall system performance.

* `sdmedian.filename` - filename for optionally saving sd.median data calculated
  on the input dataset to as CSV. Defaults to `""`, for which this data will not
  be saved. Use for extracting medians for parallel processing scenarios other
  than the built-in parallel option. See notes on [large data sets](#largedata)
  for details.

* `sdrecentered.filename` - filename to save re-centered data to as CSV.
  Defaults to "", for which this data will not be saved. Useful for
  post-processing and debugging.

* `quietly` - default `TRUE`; when `TRUE`, displays function messages and will
  output log files when `parallel` is `TRUE`.

* `log.path` - default `"."`; sets directory for batch log file output when
  processing with `parallel = TRUE`. A new directory will be created if necessary.

## <a name="output"></a>Understanding output

As seen [in the example above](#example), `growthcleanr` creates a vector of
flags identifying which measurements should be included or excluded, and why.
This vector can be added to the original dataset for further review. The flags
`growthcleanr` sets are text labels, most of which are descriptive. The text
does not always map directly back to the precise step in the algorithm.

For convenience, the following table maps which `algorithm step` in the
algorithm, as defined by Supplemental File 3 within the
[Supplementary Material](https://academic.oup.com/jamia/article/24/6/1080/3767271#97610899)
of the Daymont et al. paper, define exclusions in the algorithm by the exclusion
identifier also prescribed in the algorithm ("`algorithm exc id`"). The
`exclusion type` column identifies the text of the exclusion flag set by
`growthcleanr`. The `notes` column at right identifies minor discrepancies between
the algorithm's step labels and labels used in comment text in `growthcleanr`.

| algorithm step | algorithm exc id | exclusion type | notes |
| - | - | - | - |
| 2d | 0 | Include | - |
| 2d | 1 | Missing | - |
| 5b | 2 | Exclude-Temporary-Duplicate | - |
| 7d | - | Swapped-Measurement | - |
| 8f | - | Unit-Error-High | - |
| 8f | - | Unit-Error-Low | - |
| 8f | - | Unit-Error-Possible | Not set in R code |
| 9c | 3 | Exclude-Carried-Forward | - |
| 10c | 4 | Exclude-SD-Cutoff | 10d, 10e |
| 11d | 5 | Exclude-EWMA-Extreme | 11e |
| 11f.ii | 6 | Exclude-EWMA-Extreme-Pair | 11i (R only) |
| 12d.i | 7 | Exclude-Duplicate | 12diii, 12ei, 12f |
| 14f.i | 8 | Exclude-EWMA-8 | Set in 14h (in R) |
| 14f.ii | 9 | Exclude-EWMA-9 | Set in 14h (in R) |
| 14f.iii | 10 | Exclude-EWMA-10 | Set in 14h (in R) |
| 14f.iv | 11 | Exclude-EWMA-11 | Set in 14h (in R) |
| 14f.v | 12 | Exclude-EWMA-12 | Set in 14h (in R) |
| 14f.vi | 13 | Exclude-EWMA-13 | Set in 14h (in R) |
| 14f.vii | 14 | Exclude-EWMA-14 | Set in 14h (in R) |
| 15s | 15 | Exclude-Min-Height-Change | Set in 15q (in R), 15s, 15t |
| 15s | 16 | Exclude-Max-Height-Change | Set in 15q (in R), 15s, 15t |
| 16b | 17 | Exclude-Pair-Delta-17 | - |
| 16b | 18 | Exclude-Pair-Delta-18 | - |
|  -  |  - | Exclude-Pair-Delta-19 | Added `flag.both` mode |
| 16d | 19 | Exclude-Single-Outlier | - |
| 17a | 20 | Exclude-Too-Many-Errors | - |
| 17a | 21 | Exclude-Too-Many-Errors-Other_Parameter | - |

A researcher wanting to compare two datasets -- for example, one including only
measurements `growthcleanr` flags for inclusion, and a second which also
includes carried forward values -- could extract only measurements where the
flag is set to "`Include`" for the first set and either "`Include`" or
"`Exclude-Carried-Forward`" for the second set. This flexibility gives each
researcher control over how to incorporate growthcleanr results into their own
research.

## <a name="largedata"></a> Working with large data sets

The nature of the `growthcleanr` algorithm is repetitive. It runs many checks
over each individual subject's measurements, in a specific sequence, and revises
its assessments as it goes. Because of this, it can take some time to process
even a modest dataset. For reference, the `syngrowth` synthetic data example
packaged with `growthcleanr` takes 2-3 minutes to process on a contemporary
laptop. `growthcleanr` uses optimized libraries (e.g. `data.table`) to go as
fast as possible, but there are limits given the repeated passes over data
required.

If you have under one million records, using `growthcleanr` may just
necessitate planning to wait a little while the job processes. There are
several strategies to improve performance you might consider:

 * The built-in `parallel` option should take advantage of multiple CPU cores
   when your hardware allows.
 * Running `growthcleanr` on a machine with more CPU cores and more RAM can
   improve this further.
 * Verify that supporting libraries like `data.table` are
   [installed properly](#platform) for your machine
 * Anecdotally, Windows users have sometimes seen slower performance when
   compared with similarly-appointed hardware running Linux or macOS.

For very large datasets, on the order of millions of records or more, cleaning
data can take several hours or more. When running this kind of job, the risk of
failures due to external factors such as power outages increases. The current
structure of `growthcleanr` does not checkpoint results as progress is made, so
if a job has to run overnight and fails in the middle, a researcher would likely
need to start over from scratch.

Because `growthcleanr` operates for the most part on individual subjects one at
a time, however, this issue might be mitigated by splitting the input data into
many small files, then running `growthcleanr` separately on each file, with
results re-combined at the end. The primary benefit of this strategy would be
the saving of the results of each file as it completes, allowing a researcher
whose job fails overnight to resume processing on only the as yet incomplete
input files. A secondary benefit might be a slight improvement of use of
available RAM.

Adopting this approach might require some custom code, and there are few
pitfalls to avoid. The following lays out a rough approach:

* Extract `sd.median` values for the entire dataset as a whole using the
  `sdmedian.filename` option. This is critical, as `growthcleanr` will generate
  these values itself for an input dataset if they are not provided. If a large
  dataset is split into 1,000 smaller files to be cleaned separately, each of
  those 1,000 `growthcleanr` jobs needs to use the same `sd.median` values to
  recenter or the results will be inconsistent. Pass this data in to
  `cleangrowth()` using the parameter `sd.recenter`.
* Split the data into many small files, but keep each individual subject's
  measurements together in one file. If, for example, a subject has 12
  measurements, all 12 should be in one and only one of the smaller files to
  maximize the longitudinal analysis of their values. The `splitinput()`
  function will perform this split safely, e.g.:

```R
> count = splitinput(syngrowth, fname="mydata")
> count
[1] 5
> list.files(".", pattern="mydata.*")
[1] "mydata.00000.csv" "mydata.00001.csv" "mydata.00002.csv" "mydata.00003.csv" "mydata.00004.csv" "mydata.00005.csv"
```

* Write a standalone driver script to execute `growthcleanr` on each separate
  file, then write out its results when complete. An example is below.
* Invoke the job using a tool like
  [GNU Parallel](https://www.gnu.org/software/parallel/), which can be run on
  all major platforms. An example invocation is below as well.
* If a failure occurs, set aside both the completed inputs and their results.
  Then re-run the job on only the remaining input data.
* When complete, re-combine the data into one file or as otherwise appropriate.

A standalone driver script might look like the following (note that it uses the
`argparse` library for processing options; this should make it easy to read, but
requires additional installation steps from what is recommended above under
[Installation](#installation)).

```R
library(argparse)
library(data.table)
library(growthcleanr)

parser <- ArgumentParser(description='Process some integers')
parser$add_argument('infile', metavar='INFILE', type="character", nargs=1,
                    help='input file')
parser$add_argument('outfile', metavar='OUTFILE', type='character', nargs=1,
                    help='output file')
parser$add_argument('--sdrecenter', type='character', nargs=1, default='',
                    help='sd.recenter data file')
parser$add_argument('--quietly', default=F, action='store_true',
                    help='Disable verbose output')
args <- parser$parse_args()

if (args$sdrecenter != '' ) {
  sdrecenter <- fread(args$sdrecenter)
} else {
  sdrecenter <- ''
}

df_in <- fread(args$infile)
df_out <- df_in[, clean_value:=
                cleangrowth(subjid, param, agedays, sex, measurement,
                            sd.recenter=sdrecenter,
                            quietly=args$quietly)
           ]
write.csv(df_out, args$outfile, row.names=FALSE)
```

A file like this could be saved as, for example, `gcdriver.R`, which can then
be invoked on a single input file using `Rscript`:

```bash
% Rscript gcdriver.R --quietly --sdrecenter sdmedians.csv mydata.csv mydata-cleaned.csv
```

If you have many small files saved, for example, as `mydata.00001.csv`,
`mydata.00002.csv`, etc., this driver can be invoked using `parallel`:

```bash
% ls mydata.?????.csv | parallel -j2 --eta \
  "Rscript gcdriver.R --quietly --sdrecenter sdmedians.csv {} {}-clean"
```

This lists your input files, and passes the filename list to parallel to use in
invoking the driver script, one file at a time, and to process as parallel jobs
as resources are available. Each run of the job then saves the cleaned output
with `-clean` appended to the input filename, and as each completes, the next
file on the list will be started until all are complete. A few things to note:

* The careful use of naming and wildcards when listing input files and naming
  output files can save accidental re-running of data from output files
* The `-j2` option specifies running two jobs at once; `-j0` would use as many
  CPU cores as a machine has available. This might need to vary to match
  specific hardware.
* The `--eta` option will report on progress.
* The `--quietly` option on the driver script will make it easier to monitor
  progress with less verbose output coming from `growthcleanr`.
* The `--sdrecenter` option on the driver script should be set to ensure each
  individual file is recentered using the same set for the entire input.
* If multiple cores are available, the job should proceed with speedup roughly
  similar to what can be gained using the parallel batching feature built in to
  `growthcleanr`, with the main difference being the saving of intermediate
  output as each smaller file completes.

If your data to be cleaned is very large, it might help to store it compressed,
for example with `gzip` and its corresponding `.gz` filename extension.
`growthcleanr` can read in `.gz` input, but you might need to install the
`R.utils` package first. `R` will provide a message if this is required.

## <a name="related"></a>Related tools

The CDC provides a
[SAS Program for the 2000 CDC Growth Charts](https://www.cdc.gov/nccdphp/dnpao/growthcharts/resources/sas.htm)
which can also be used to identify biologically implausible values using a
different approach.

[GrowthViz](https://github.com/mitre/GrowthViz) provides insights into how
`growthcleanr` assesses data, packaged in a Jupyter notebook. It ships with the
same `syngrowth` synthetic example dataset as `growthcleanr`, with results
included.
