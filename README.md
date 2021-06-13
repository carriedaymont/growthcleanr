
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
observations and utilities for working with both functions. As of summer
2021, it also supports cleaning anthropometric measurements for adults
up to age 65. The adult algorithm has not yet been published in a
peer-reviewed publication, but is described in detail at [Adult
algorithm](articles/adult-algorithm.html).

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

-   [Quickstart](articles/quickstart.html) - a brief tour of using
    growthcleanr, including data preparation
-   [Installation](articles/installation.html) - including required
    dependencies and OS-specific notes
-   [Usage](articles/usage.html) - examples of cleaning data, multiple
    options, example data

### Advanced topics:

-   [Configuration options](articles/configuration.html) - changing
    growthcleanr operational settings
-   [Understanding growthcleanr output](articles/output.html) - the
    exclusion types growthcleanr identifies
-   [Adult algorithm](articles/adult-algorithm.html) - a detailed
    description of how growthcleanr assesses observations from adult
    subjects
-   [Computing BMI percentiles and Z-scores](articles/utilities.html) -
    additional functions for common data transforms and determining
    percentiles and Z-scores using the CDC method
-   [Working with large datasets](articles/large-data-sets.html) - notes
    and suggestions for running `growthcleanr` with large data sources
-   [Testing carried forward
    adjustments](articles/adjust-carry-forward.html) - notes on an
    experimental tool for refining how pediatric measurements are
    assessed

## Changes

For a detailed history of released versions, see the
[Changelog](news/index.html) or`NEWS.md`. Tagged releases, starting with
1.2.3 in January 2021, are listed [at
GitHub](https://github.com/carriedaymont/growthcleanr/releases).
