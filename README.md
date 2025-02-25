
<!-- README.md is generated from README.Rmd. Please edit that file -->

# growthcleanr

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/growthcleanr)](https://cran.r-project.org/package=growthcleanr)
[![R build
status](https://github.com/carriedaymont/growthcleanr/workflows/R-CMD-check/badge.svg)](https://github.com/carriedaymont/growthcleanr/actions)
[![Docker](https://github.com/mitre/growthcleanr/actions/workflows/build-and-publish-image-tag.yml/badge.svg)](https://github.com/mitre/growthcleanr/actions/workflows/build-and-publish-image-tag.yml)

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
observations and utilities for working with both functions. As of summer
2021, it also supports cleaning anthropometric measurements for adults
up to age 65. The adult algorithm has not yet been published in a
peer-reviewed publication, but is described in detail at [Adult
algorithm](https://carriedaymont.github.io/growthcleanr/articles/adult-algorithm.html).

## Installation

To install the stable version from CRAN:

``` r
install.packages("growthcleanr")
```
Or from Github directly:

``` r
library(devtools)
install_github("carriedaymont/growthcleanr")
```

## Summary

The `growthcleanr` package processes data prepared in a specific format
to identify biologically implausible height and weight measurements. It
bases these evaluations on techniques which use patient-specific
longitudinal analysis and variations from published growth trajectory
charts for pediatric subjects. These techniques are performed in a
specific order which refines and improves results throughout the
process.

Results from `growthcleanr` include a flag for each measurement
indicating whether it is to be included or excluded based on
plausibility, with a variety of specific types of exclusions identified
distinctly. These flags can be analyzed further by researchers studying
anthropometric EHR data to determine which measurements to include or
exclude in their own studies. No values are deleted or otherwise
removed; each is only flagged in a new column.

To start running `growthcleanr`, an R installation with a variety of
additional packages is required, as is a growth measurement dataset
prepared for use in `growthcleanr`.

The rest of this documentation includes:

### Getting started:

- [Quickstart](https://carriedaymont.github.io/growthcleanr/articles/quickstart.html),
  a brief tour of using growthcleanr, including data preparation
- [Installation](https://carriedaymont.github.io/growthcleanr/articles/installation.html),
  options for installing growthcleanr, with notes on specific platforms
  and source-level installation for developers
- [Usage](https://carriedaymont.github.io/growthcleanr/articles/usage.html),
  examples of cleaning data, multiple options, example data

### Advanced topics:

- [Configuration
  options](https://carriedaymont.github.io/growthcleanr/articles/configuration.html),
  changing growthcleanr operational settings
- [Understanding growthcleanr
  output](https://carriedaymont.github.io/growthcleanr/articles/output.html),
  the exclusion types growthcleanr identifies
- [Adult
  algorithm](https://carriedaymont.github.io/growthcleanr/articles/adult-algorithm.html),
  a detailed description of how growthcleanr assesses observations from
  adult subjects
- [Computing BMI percentiles and
  Z-scores](https://carriedaymont.github.io/growthcleanr/articles/utilities.html),
  additional functions for common data transforms and determining
  percentiles and Z-scores using the CDC method
- [Working with large
  datasets](https://carriedaymont.github.io/growthcleanr/articles/large-data-sets.html),
  notes and suggestions for running `growthcleanr` with large data
  sources
- [Next
  steps](https://carriedaymont.github.io/growthcleanr/articles/next-steps.html),
  notes on potential enhancements to the pediatric and adult algorithms
- [Developer
  guidelines](https://carriedaymont.github.io/growthcleanr/articles/developer-guidelines.html),
  advice for contributors to this package, including a CRAN release
  checklist

## Changes

For a detailed history of released versions, see the
[Changelog](https://carriedaymont.github.io/growthcleanr/news/index.html)
or`NEWS.md`. Tagged releases, starting with 1.2.3 in January 2021, are
listed [at
GitHub](https://github.com/carriedaymont/growthcleanr/releases).
