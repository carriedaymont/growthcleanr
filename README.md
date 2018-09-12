# growthcleanr

R package for cleaning growth measurements

This package is used for cleaning height and weight measurements.

## Installation

You can install the package using `devtools` in the R console with:

```bash
devtools::install_github("carriedaymont/growthcleanr")
```

Or you can download the source code in case you want to make changes, then install it from source:

* Clone the github source:

```bash
git clone https://github.com/carriedaymont/growthcleanr.git
```

* Open an R session from the `growthcleanr` base directory.  Install using R `devtools` package:

```R
devtools::install(".")
```

You can also install the package from an installation file if one is obtained.

## Usage

Once the `growncleanr` package is installed, the minimal vectors needed are:

1. subjid - a unique identifier for each subject.
2. param - a string designating the measurement being 'HEIGHTCM' or 'WEIGHTKG'.
3. agedays - a numeric value for age in days at time of measurement.
3. sex - a numeric value of 0 (male) or 1 (female).
4. measurement - numeric value of measurement in cm for height or kg for weight.
5. include.carryforward - (optional), set to `TRUE` if values carried forward should not be excluded.
6. log.path - (optional), set path for log file output (e.g. `log/`).  **Will not produce log files if omitted.**
7. parallel - (optional), set to `TRUE` if program needs to run in parallel (recommended for large datasets, > 100K rows)
8. num.batches - If running in parallel, the number of batches and cores to use (caution, may effect overall system performance)

The `cleangrowth` function will return a for each measurement that will determine whether to `Include` the measurement, conclude that the measurement is `Missing`, or to `Exclude` the measurement (numerous variations of `Exclude` exist).

### Example

For a data.frame object `source_data` containing growth data:

```R
library(growthcleanr)

# save data as data.table
data<-as.data.table(source_data)

# set the data.table key for better indexing
setkey(data,subjid,param,agedays)

# generate new field using function
clean_data<-data[, clean_value:= cleangrowth(subjid, param, agedays, sex, measurement)]

# limit data only to values suggested to be included:
output_data<-clean_data[clean_vaue=='Include']
```