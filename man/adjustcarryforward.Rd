% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjustcarryforward.R
\name{adjustcarryforward}
\alias{adjustcarryforward}
\title{adjustcarryforward
\code{adjustcarryforward} Uses absolute height velocity to identify values
excluded as carried forward values for reinclusion.}
\usage{
adjustcarryforward(
  subjid,
  param,
  agedays,
  sex,
  measurement,
  orig.exclude,
  exclude_opt = 0,
  sd.recenter = NA,
  ewma.exp = -1.5,
  ref.data.path = "",
  quietly = TRUE,
  minfactor = 0.5,
  maxfactor = 2,
  banddiff = 3,
  banddiff_plus = 5.5,
  min_ht.exp_under = 2,
  min_ht.exp_over = 0,
  max_ht.exp_under = 0.33,
  max_ht.exp_over = 1.5
)
}
\arguments{
\item{subjid}{Vector of unique identifiers for each subject in the database.}

\item{param}{Vector identifying each measurement, may be 'WEIGHTKG', 'HEIGHTCM', or 'LENGTHCM'
'HEIGHTCM' vs. 'LENGTHCM' only affects z-score calculations between ages 24 to 35 months (730 to 1095 days).
All linear measurements below 731 days of life (age 0-23 months) are interpreted as supine length, and
all linear measurements above 1095 days of life (age 36+ months) are interpreted as standing height.
Note: at the moment, all LENGTHCM will be converted to HEIGHTCM. In the future, the algorithm will be updated to consider this difference.}

\item{agedays}{Numeric vector containing the age in days at each measurement.}

\item{sex}{Vector identifying the gender of the subject, may be 'M', 'm', or 0 for males, vs. 'F',
'f' or 1 for females.}

\item{measurement}{Numeric vector containing the actual measurement data.  Weight must be in
kilograms (kg), and linear measurements (height vs. length) in centimeters (cm).}

\item{orig.exclude}{Vector of exclusion assessment results from cleangrowth()}

\item{exclude_opt}{Number from 0 to 3 indicating which option to use to handle strings of carried-forwards:
0. no change.
\enumerate{
\item when deciding to exclude values, if we have a string of carried forwards,
drop the most deviant value, and all CFs in the same string, and move on as
normal.
\item when deciding to exclude values, if the most deviant in a
string of carried forwards is flagged, check all the CFs in that
string from 1:N. Exclude all after the first that is flagged for
exclusion when comparing to the Include before and after. Do not
remove things designated as include.
\item when deciding to exclude values, if the most deviant in a
string of carried forwards is flagged, check all the CFs in that
string from 1:N. Exclude all after the first that is flagged for
exclusion when comparing to the Include before and after. Make sure
remove things designated as include.
}}

\item{sd.recenter}{Data frame or table with median SD-scores per day of life}

\item{ewma.exp}{Exponent to use for weighting measurements in the exponentially weighted moving
average calculations. Defaults to -1.5. This exponent should be negative in order to weight growth
measurements closer to the measurement being evaluated more strongly. Exponents that are further from
zero (e.g. -3) will increase the relative influence of measurements close in time to the measurement
being evaluated compared to using the default exponent.}

\item{ref.data.path}{Path to reference data. If not supplied, the year 2000
Centers for Disease Control (CDC) reference data will be used.}

\item{quietly}{Determines if function messages are to be displayed and if log files (parallel only)
are to be generated. Defaults to TRUE.}

\item{minfactor}{Sweep variable for computing mindiff.next.ht in 15f, default 0.5}

\item{maxfactor}{Sweep variable for computing maxdiff.next.ht in 15f, default 2}

\item{banddiff}{Sweep variable for computing mindiff.next.ht in 15f, default 3}

\item{banddiff_plus}{Sweep variable for computing maxdiff.next.ht in 15, default 5.5}

\item{min_ht.exp_under}{Sweep variable for computing ht.exp in 15f, default 2}

\item{min_ht.exp_over}{Sweep variable for computing ht.exp in 15f, default 0}

\item{max_ht.exp_under}{Sweep variable for computing ht.exp in 15f, default 0.33}

\item{max_ht.exp_over}{Sweep variable for computing ht.exp in 15f, default 1.5}
}
\value{
Re-evaluated exclusion assessments based on height velocity.
}
\description{
adjustcarryforward
\code{adjustcarryforward} Uses absolute height velocity to identify values
excluded as carried forward values for reinclusion.
}
\examples{
\donttest{
# Run on a small subset of given data
df <- as.data.frame(syngrowth)
df <- df[df$subjid \%in\% unique(df[, "subjid"])[1:2], ]
clean_df <- cbind(df,
                  "gcr_result" = cleangrowth(df$subjid,
                                             df$param,
                                             df$agedays,
                                             df$sex,
                                             df$measurement))

# Adjust carry forward values in cleaned data
adj_clean <- adjustcarryforward(subjid = clean_df$subjid,
                                param = clean_df$param,
                                agedays = clean_df$agedays,
                                sex = clean_df$sex,
                                measurement = clean_df$measurement,
                                orig.exclude = clean_df$gcr_result)
}
}
