# Z-Score Calculation Reference — growthcleanr

This document describes every z-score computed by growthcleanr v3.0.0, how
each is calculated, which reference tables are used, and what is returned by
`cleangrowth()`.

growthcleanr uses **two fundamentally different z-score methods** for different
purposes. Understanding the distinction is essential for correctly interpreting
output or replicating calculations downstream.

---

## 1. Two Z-Score Methods

| Feature | CSD Method | LMS Method |
|---------|-----------|------------|
| **Used for** | Main cleaning algorithm (`sd.orig`, `tbc.sd`) | Diagnostic reference only (`z.orig`); extended BMI (`ext_bmiz`) |
| **Formula** | `(X - M) / csd` | `((X/M)^L - 1) / (L*S)` |
| **Effect at extremes** | Symmetric above/below median; less sensitive to skew at tails | Accounts for skewness via L; compressed z-scores at upper tail for right-skewed params (e.g., weight) |
| **Why used for cleaning** | Large weight in a very heavy child produces a larger z-change, making implausible values easier to detect |  |
| **Rounding** | Stata-style to 0.001 (half away from zero) | Stata-style to 0.001 |

---

## 2. CSD Method — `sd.orig`

The cleaning algorithm uses the **Conditional Standard Deviation (CSD)** method.
This divides the difference from the median by a half-SD derived from the
reference distribution, rather than the full LMS-adjusted SD.

### Formula

```
If measurement < M:   sd.orig = (measurement - M) / csd_neg
If measurement >= M:  sd.orig = (measurement - M) / csd_pos
```

Where:
- `M` = age/sex-specific median from reference table
- `csd_pos` = half the distance from M to the value at z = +2:
  `csd_pos = (M * (1 + 2*L*S)^(1/L) - M) / 2`
- `csd_neg` = half the distance from M to the value at z = -2:
  `csd_neg = (M - M * (1 - 2*L*S)^(1/L)) / 2`

The `csd_pos` and `csd_neg` values are pre-computed and stored in the
reference tables as columns `csdpos` and `csdneg`.

### Reference tables

| Table | Source | Used for |
|-------|--------|----------|
| `growthfile_who.csv.gz` | WHO 2006 | HEIGHTCM/WEIGHTKG age < 5y, all HEADCM |
| `growthfile_cdc_ext_infants.csv.gz` | CDC 2000 extended | HEIGHTCM/WEIGHTKG age >= 2y |

### Age blending — child algorithm (default, v3.0.0)

The child algorithm computes `sd.orig_who` and `sd.orig_cdc` separately, then
blends them based on age and parameter:

| Age range | Parameter | Formula |
|-----------|-----------|---------|
| < 2 years | any | `sd.orig = sd.orig_who` |
| 2–5 years | HEIGHTCM, WEIGHTKG | `sd.orig = (sd.orig_cdc * (age_yr - 2) + sd.orig_who * (5 - age_yr)) / 3` |
| > 5 years | HEIGHTCM, WEIGHTKG | `sd.orig = sd.orig_cdc` |
| any age | HEADCM | `sd.orig = sd.orig_who` |

At the boundaries (exactly 2y or 5y), the formula yields 100% WHO and 100%
CDC respectively (both arms equal at the boundary).

Both `sd.orig_who` and `sd.orig_cdc` are individually rounded to 0.001 before
blending. The blended result is then rounded to 0.001.

### Age blending — legacy pediatric algorithm (`use_legacy_algorithm = TRUE`)

The legacy algorithm uses the same WHO/CDC tables and formula but with a
narrower blending window (ages 2–4y) and divides by 2 instead of 3. For
corrected/uncorrected smoothing it also uses ages 2–4y divided by 2.

### Age blending — `calc_and_recenter_z_scores()` (CF rescue)

This internal helper — used during carry-forward rescue in Step 6 — uses
the same blending formula and WHO/CDC tables, but with a CDC-only cutoff
at age >= 4y (not > 5y):

| Age range | Parameter | Formula |
|-----------|-----------|---------|
| < 2 years | any | WHO only |
| 2–5 years | HEIGHTCM, WEIGHTKG | `(cdc*(age-2) + who*(5-age))/3` |
| >= 4 years | HEIGHTCM, WEIGHTKG | CDC only |
| any age | HEADCM | WHO only |

---

## 3. LMS Method — `z.orig`

`z.orig` is the standard Box-Cox z-score using the **CDC reference only**
(no WHO blending). It is included in the output as a reference but is
**not used by the cleaning algorithm**.

### Formula

```
If L != 0:  z.orig = ((X / M)^L - 1) / (L * S)
If L == 0:  z.orig = log(X / M) / S
```

Where `L`, `M`, `S` are the Box-Cox power, median, and coefficient of
variation from the CDC reference table.

**Note:** This is a CDC-only LMS z-score, not blended with WHO. It differs
from the smoothed LMS z-score used in downstream quality measure code
(which uses WHO/CDC blending with the same age windows as the CSD method).

---

## 4. Recentered Scores — `tbc.sd` and `ctbc.sd`

The cleaning algorithm works on **recentered** z-scores so that outlier
detection is relative to the dataset's own population rather than an
external reference standard.

### `tbc.sd` (recentered blended CSD)

```
tbc.sd = sd.orig - sd.median
```

Where `sd.median` is the **median of `sd.orig`** computed across all
non-excluded observations for a given `param` + `sex` combination within
the current dataset.

This means a subject near the population median of the input dataset will
have `tbc.sd` near 0, regardless of their position relative to the
national reference.

### `ctbc.sd` (corrected recentered CSD)

Same formula, but uses `sd.corr` (gestational-age-corrected z-score for
preterm infants) instead of `sd.orig`:

```
ctbc.sd = sd.corr - sd.median
```

Used only for subjects eligible for Fenton/gestational-age correction.

---

## 5. Gestational-Age-Corrected Scores — `sd.corr`, `sd.orig_uncorr`

For preterm infants, growthcleanr applies a gestational-age correction
using Fenton growth charts.

| Column | Description |
|--------|-------------|
| `sd.orig_uncorr` | Copy of `sd.orig` before gestational-age correction |
| `sd.corr` | CSD z-score after Fenton/gestational-age correction |

When no correction applies (term infants, older children), `sd.corr = sd.orig`
and `sd.orig_uncorr = sd.orig`.

---

## 6. EWMA (Exponentially Weighted Moving Average)

The EWMA is the central mechanism for detecting implausible values. It is
computed from **recentered CSD z-scores** (`tbc.sd`).

### Weight formula

```
weight_i = (5 + |agedays_current - agedays_i|) ^ exponent
```

- Default exponent: `-1.5`
- Measurements closer in time receive higher weight
- Default window: 15 observations on each side of the current measurement

### EWMA variants (internal — not in output)

| Variable | Excludes |
|----------|----------|
| `ewma.all` | Only the current observation |
| `ewma.before` | Current + immediately prior observation |
| `ewma.after` | Current + immediately subsequent observation |

The deviation from EWMA is `dewma = tbc.sd - ewma`.

EWMA values are **not included in `cleangrowth()` output**. They are
computed internally during each algorithm iteration and discarded.

---

## 7. Extended BMI Z-Scores — `ext_bmiz()`

The standalone `ext_bmiz()` function computes extended BMI z-scores for
children 2–19.9 years, following the CDC's extended BMI method.

For children with BMI at or below the 95th percentile, the extended z-score
equals the standard LMS z-score. For children above the 95th percentile, a
half-normal distribution with a fitted sigma parameter is used:

### Sigma formula (scale parameter)

```
Males:   sigma = 0.3728 + 0.5196 * age_yr - 0.0091 * age_yr^2
Females: sigma = 0.8334 + 0.3712 * age_yr - 0.0011 * age_yr^2
```

Where `age_yr` = age in years (= age in months / 12).

### Extended BMI percentile

```
For BMI <= 95th percentile:  ebp = 100 * pnorm(bz)        [standard LMS]
For BMI > 95th percentile:   ebp = 90 + 10 * pnorm((BMI - p95) / sigma)
```

The extended BMI z-score is then `qnorm(ebp / 100)`. If ebp > 99.99999999999999,
extended BMIz is capped at 8.21.

Reference: CDC SAS macro (updated December 2022); NCHS report doi:10.15620/cdc:121711.

---

## 8. Rounding

All CSD z-scores and LMS z-scores are rounded to **3 decimal places** using
Stata-style rounding (half away from zero). This differs from R's default
`round()`, which uses banker's rounding (round half to even).

`sd.orig_who` and `sd.orig_cdc` are rounded individually before blending.
The blended `sd.orig` is rounded again after combining.

---

## 9. Columns Returned by `cleangrowth()` (v3.0.0)

`cleangrowth()` returns a data.table with all original input columns plus:

| Column | Type | Description |
|--------|------|-------------|
| `exclude` | factor (then character) | Exclusion code. `"Include"` = retained by algorithm; `"Exclude-*"` = excluded with reason code. |
| `cf_rescued` | character | Carry-forward rescue reason (child algorithm only). Empty string `""` for non-rescued rows and adult rows. |
| `sd.orig` | numeric | Blended CSD z-score (see Section 2). The primary z-score used by the cleaning algorithm. |
| `tbc.sd` | numeric | Recentered blended CSD score (see Section 4). |
| `sd.orig_uncorr` | numeric | Pre-correction copy of `sd.orig` (see Section 5). For term infants/children: equals `sd.orig`. |
| `sd.corr` | numeric | Gestational-age-corrected CSD score (see Section 5). For term infants/children: equals `sd.orig`. |
| `ctbc.sd` | numeric | Recentered corrected CSD score (see Section 4). For term infants/children: equals `tbc.sd`. |

**Note:** `z.orig` (CDC-only LMS z-score, described in Section 3) is computed
during the algorithm but is **not included in the v3.0.0 output**. It is
retained in the code as a reference but removed from the returned data.table
via the checkpoint columns list.

**Note on `line`:** An internal row index column (`line`) is also present in
the output but is not meaningful to end users. It is used internally for
merging.

---

## 10. Key Distinctions for Downstream Use

### `sd.orig` vs a standard LMS z-score

`sd.orig` uses CSD, which is **not** the conventional z-score for growth
research. It is intentionally asymmetric to make the algorithm more sensitive
to implausible high values. Do not use `sd.orig` or `tbc.sd` as
population-referenced z-scores in analyses; compute LMS z-scores separately.

### `z.orig` vs a blended LMS z-score

`z.orig` uses CDC reference only. Standard research practice blends WHO
(younger children) and CDC (older children). For age-blended LMS z-scores,
use the same WHO/CDC blend as the cleaning algorithm (WHO-only < 2y,
weighted blend 2–5y ÷ 3, CDC-only > 5y) but applying the LMS formula.

### EWMA output

EWMA/DEWMA values are computed internally and **not returned**. If you
need DEWMA for quality measures, compute it independently using your
preferred z-score input (typically blended LMS z-scores).

---

## 11. Reference Table Summary

| File | Location | Content |
|------|----------|---------|
| `growthfile_who.csv.gz` | `inst/extdata/` | WHO 2006 L/M/S + `csdpos`/`csdneg` |
| `growthfile_cdc_ext_infants.csv.gz` | `inst/extdata/` | CDC 2000 extended L/M/S + csd columns |
| `CDCref_d.csv.gz` | `inst/extdata/` | CDC reference for `ext_bmiz()` |
| `who_hc_vel_3sd_infants.csv.gz` | `inst/extdata/` | WHO HC velocity thresholds (age < 3y) |
| `who_hc_maxvel_3sd_infants.csv.gz` | `inst/extdata/` | WHO HC max velocity thresholds (age < 3y) |
| Fenton tables | `inst/extdata/` | Gestational-age correction for preterm infants |
