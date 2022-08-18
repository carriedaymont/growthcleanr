# CRAN submission growthcleanr 2.0.1 (3)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking installed package size ... NOTE
    installed size is  9.8Mb
    sub-directories of 1Mb or more:
      doc       5.7Mb
      extdata   3.3Mb

doc should not be included in build (listed in .Rbuildignore) and extdata has been compressed as much as possible -- data in that folder is necessary to run the package.

## CRAN pretest

We have addressed each of the notes in the pretest as below:

* Possibly misspelled words in DESCRIPTION:
    Anthropometric (3:25)
    anthropometric (22:56)

These are not misspelled.

*   Examples with CPU (user + system) or elapsed time > 10s
              user system elapsed
  cleangrowth 4.01    0.3    18.2

cleangrowth() example has been changed to \donttest{}

## Downstream dependencies
There are currently no downstream dependencies for this package.

# Previous Submissions

# CRAN submission growthcleanr 2.0.1 (2)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking installed package size ... NOTE
    installed size is  9.8Mb
    sub-directories of 1Mb or more:
      doc       5.7Mb
      extdata   3.3Mb

doc should not be included in build (listed in .Rbuildignore) and extdata has been compressed as much as possible -- data in that folder is necessary to run the package.

## CRAN pretest

We have addressed each of the notes in the pretest as below:

* Possibly misspelled words in DESCRIPTION:
    Anthropometric (3:25)
    anthropometric (22:56)

These are not misspelled.

*   Found the following (possibly) invalid URLs:
    URL: https://cran.r-project.org/web/packages/survey/index.html
      From: inst/doc/configuration.html
      Status: 200
      Message: OK
      CRAN URL not in canonical form
    The canonical URL of the CRAN page for a package is 
      https://CRAN.R-project.org/package=pkgname

This has been fixed.

*   The Description field should start with a capital letter.

This has been fixed ('growthcleanr' removed, includes capitalized).

*   Examples with CPU (user + system) or elapsed time > 10s
              user system elapsed
  cleangrowth 4.01    0.3    18.2

cleangrowth() example has been shortened (run on fewer subjects).

* Check: Overall checktime, Result: NOTE
  Overall checktime 13 min > 10 min

Tests have been shortened (elapsed time locally ~80s).

## Downstream dependencies
There are currently no downstream dependencies for this package.

# CRAN submission grwothcleanr 2.0.1 (1)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking installed package size ... NOTE
    installed size is  9.8Mb
    sub-directories of 1Mb or more:
      doc       5.7Mb
      extdata   3.3Mb

doc should not be included in build (listed in .Rbuildignore) and extdata has been compressed as much as possible -- data in that folder is necessary to run the package.

## Downstream dependencies
There are currently no downstream dependencies for this package.
