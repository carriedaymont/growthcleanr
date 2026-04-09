# Adult growthcleanr Permissiveness Framework Specification

---

## Overview

The adult growthcleanr algorithm includes a permissiveness
framework that allows users to control how aggressively
measurements are excluded. This is controlled by a single
`permissiveness` parameter with four levels. The
`permissiveness` parameter sets defaults for all
sub-parameters, but individual sub-parameters can be
overridden. Explicit sub-parameter values always take
precedence over the permissiveness preset.

---

## Permissiveness Levels

**Current default looser** 

| Level | Intended Use | Weight | Height | RV handling |
|-------|-------------|--------|--------|-------------|
| `loosest` | Retain all plausible values | Allows for extreme weight changes, such as with bariatric surgery or rapid fluid status change | Allows for a large degree height measurement error and height gain in young adults as well as significant height loss | Considers repeated identical values as more likely to be copied from a prior value (not always copied) |
| `looser` | Retain most plausible values | Allows for slightly less extreme weight changes than loosest limits | Allows for a large degree of height measurement error and height gain in young adults, but not significant height loss | Considers repeated identical values as more likely to be copied from a prior value (not always copied) |
| `tighter` | Exclude some plausible values that likely reflect illness with large effect on body size | Excludes weights and weight changes that are likely due to conditions such as severe cachexia, very rapid weight change from bariatric surgery, rapid change in fluid status | Less permissive limits for height error, allows for height gain in young adults but not significant height loss | Considers repeated identical values to be independent, with subsequent measurements just as informative as other measurements |
| `tightest` | Exclude extreme measurements even if plausible | Exclude larger weight changes even if plausible | Less permissive for height error, no allowance for height change | Considers repeated identical values to be independent, with subsequent measurements just as informative as other measurements |

---

## BIV / Size Thresholds

These thresholds apply the same whether or not BMI is
available (i.e., whether or not same-day height and weight
are both present). BMI thresholds apply only when same-day
height and weight are available.

***Default (loosest) BIV limits***
| Parameter | Default | Unit | Notes |
|-----------|---------|------|-------|
| `overall_ht_min` | 50 | cm | |
| `overall_ht_max` | 244 | cm | |
| `overall_wt_min` | 20 | kg | Not pounds - specified because scale_max_lbs uses pounds|
| `overall_wt_max` | 500 | kg | Not pounds  - specified because scale_max_lbs uses pounds |
| `overall_bmi_min` | 5 | kg/m2 | |
| `overall_bmi_max` | 300 | kg/m2 | |


***BIV (Overall) Limits for other permissiveness***
| Level | BMI | Height (cm) | Weight (kg) |
|-------|-----|-------------|-------------|
| `loosest` | current defaults | current defaults | current defaults |
| `looser` | 12–65 | 120–230 | 30–270 |
| `tighter` | 16–45 | 142–213 | 36–159 |
| `tightest` | 18–40 | 147–208 | 39–136 |


***Default (loosest) 1D limits***
| Parameter | Default | Description |
|-----------|---------|-------------|
| `single_ht_min_bmi` | 60 cm | Height minimum when BMI available and within range |
| `single_ht_max_bmi` | 245 cm | Height maximum when BMI available and within range |
| `single_wt_min_bmi` | 12 kg | Weight minimum when BMI available and within range |
| `single_wt_max_bmi` | 350 kg | Weight maximum when BMI available and within range |
| `single_ht_min_nobmi` | 122 cm | Height minimum without BMI (tighter) |
| `single_ht_max_nobmi` | 245 cm | Height maximum without BMI |
| `single_wt_min_nobmi` | 30 kg | Weight minimum without BMI (tighter) |
| `single_wt_max_nobmi` | 350 kg | Weight maximum without BMI |
| `single_bmi_min` | 10 kg/m² | BMI below this is considered extreme |
| `single_bmi_max` | 250 kg/m² | BMI above this is considered extreme |


***1D (Single) Limits for other permissiveness***
| Level | BMI | Height (cm) regardless of BMI | Weight (kg) regardless of BMI |
|-------|-----|-------------|-------------|
| `loosest` | current defaults | current defaults | current defaults |
| `looser` | 12–65 | 120–230 | 30–270 |
| `tighter` | 16–45 | 142–213 | 36–159 |
| `tightest` | 18–40 | 147–208 | 39–136 |

---

## wtallow

### Overview

wtallow is a derived value representing the threshold above which an absolute weight difference between two adjacent measurements may result in exclusion. It depends on the time interval between measurements and the Upper Weight (UW) — the higher of the two weights being compared. A base UW of 120 kg receives no adjustment; UW > 120 scales up (PW-H only), UW < 120 scales down (PW-H and PW-L).

The complete specification including all formulas, caps, ET caps, and ceilings is in **`wtallow-formulas.md`**.

### Base formulas

| Level | wtallow formula | UW scaling |
|-------|-----------------|------------|
| `loosest` | PW-H (piecewise) | Scales up for UW>120 (capped at 180); scales down for UW<120 |
| `looser` | PW-H (piecewise) | Same as loosest |
| `tighter` | PW-L (piecewise-lower) | NOT scaled up for UW>120; scales down for UW<120 |
| `tightest` | allofus15 | Cap = min(40, effective PW-L at 12m for that UW) |

### UW-based wtallow scaling

**UW = 120 (base, no adjustment):** All formulas use their base values.

**UW > 120 (PW-H only, capped at UW=180):**
- Caps increase: `adj_cap = base_cap + 0.25 × (UW − 120)`
- 1m value = 25 (vs 20 at base)
- PW-L is NOT adjusted for UW > 120

**UW < 120 (PW-H and PW-L):**
- Caps scale down: `adj_cap = (base_cap − 20) × (UW / 120) + 20`
- Hard ceiling: wtallow ≤ UW × 2/3
- Ensures lighter patients get tighter limits

**allofus15:** Not directly UW-adjusted, but its 12m cap is limited to never exceed the effective PW-L value at 12m for that UW.

### ET caps (Evil Twins and Extreme EWMA thresholds)

ET caps are derived from the wtallow formula's caps:
- **PW-H and PW-L:** ET cap = wtallow cap at 6m/12m + 20 (UW-adjusted)
- **allofus15:** ET cap = allofus15-cap-12m (flat, regardless of interval)
- **Custom CSV:** ET caps fixed at 70/100

### wtallow tables

**PW-H (loosest/looser)** at selected UW values:
```
Months   @60     @120(base)  @160    @180(cap)
  1      20.0    20.0        25.0    25.0
  6      35.0    50.0        60.0    65.0
 12      50.0    80.0        90.0    95.0
 24      50.0    80.0        90.0    95.0
```

**PW-L (tighter)** at selected UW values:
```
Months   @60     @120(base)  @160
  1      20.0    20.0        20.0
  6      26.67   33.33       33.33
 12      36.67   53.33       53.33
 24      36.67   53.33       53.33
```
Note: PW-L is not scaled up for UW > 120.


---

## ET Caps (Evil Twins and Extreme EWMA Thresholds)

ET caps control the thresholds used by Evil Twins
(Step 9Wa) and Extreme EWMA (Step 9Wb). They are
derived from the wtallow formula's caps, not set
independently.

### Derivation from wtallow

For PW-H and PW-L formulas:
- ET cap ≤6m = wtallow-cap-6m + 20
- ET cap >6m = wtallow-cap-12m + 20

Both wtallow caps and ET caps are UW-adjusted:
- UW > 120 (PW-H only): caps increase by 0.25×(UW−120), capped at UW=180
- UW < 120 (PW-H and PW-L): caps scale down; ceiling of UW × 2/3

For allofus15: ET cap = allofus15-cap-12m (flat, regardless of interval).
For custom CSV: ET caps fixed at 70/100.

### ET cap tables

**PW-H (loosest/looser)** at selected UW values:
```
UW       ET ≤6m   ET >6m
 60 kg   55.0     70.0
120 kg   70.0     100.0
160 kg   80.0     110.0
180 kg   85.0     115.0
```

**PW-L (tighter)** at selected UW values:
```
UW       ET ≤6m   ET >6m
 60 kg   51.11    56.67
120 kg   53.33    73.33
160 kg   53.33    73.33
```
Note: PW-L is not scaled up for UW > 120.

### Gap between wtallow and ET cap

The ET cap always exceeds wtallow by exactly 20 at the
cap points (6m and 12m) for PW-H and PW-L formulas.
This 20 kg gap ensures that moderate EWMA outliers
(wtallow < dewma < ET cap) are caught by Step 11Wb,
while extreme outliers (dewma > ET cap) are caught
earlier by Step 9.

See `wtallow-formulas.md` for the complete specification.

---

## Perclimit (Percentage Criterion)

The perclimit criterion is an independent exclusion
pathway in the moderate EWMA step (Step 11Wb) and the
2D ordered step (Step 11Wa). It excludes values where
`wt/ewma < perclimit` for all three EWMA calculations
(main, bef, aft), regardless of whether the standard
EWMA or trajectory criteria are met.

This criterion targets very low weight values that may
not trigger the absolute wtallow threshold but are
implausibly small relative to the patient's other
measurements.

The three thresholds correspond to the parameters
`perclimit_low` (wt ≤45 kg), `perclimit_mid`
(45<wt≤80 kg), and `perclimit_high` (wt >80 kg).

### Perclimit by permissiveness level

| Level | wt ≤ 45 kg | 45 < wt ≤ 80 kg | wt > 80 kg |
|-------|-----------|-----------------|------------|
| `loosest` | 0.5 | 0.4 | 0 (disabled) |
| `looser` | 0.5 | 0.4 | 0 (disabled) |
| `tighter` | 0.7 | 0.4 | 0.4 |
| `tightest` | 0.7 | 0.4 | 0.4 |

### Application

In the moderate EWMA step (Step 11Wb), perclimit is
applied at the observation level: each weight
observation's perclimit is determined by that
observation's weight value.

In the 2D ordered step (Step 11Wa), perclimit is
applied at the subject level: the maximum perclimit
across all weight observations for the subject is used.

### Rationale

**Loosest/looser at ≤45 kg (0.5):** A perclimit of 0.5
means a weight must be at least 50% of its EWMA to
survive. This accommodates very underweight patients
(e.g., severe cachexia, anorexia) whose weights may
legitimately be quite low relative to their historical
measurements. The current value of 0.7 is too
restrictive at loosest/looser — a cachectic patient at
35 kg with EWMA of 52 (percewma = 0.67) would be
excluded under 0.7 but retained under 0.5.

**Tighter/tightest at ≤45 kg (0.7):** At these levels,
the goal is to exclude values likely reflecting severe
illness, so the more restrictive threshold is
appropriate.

**All levels at 45–80 kg (0.4):** At these weights, a
value less than 40% of its EWMA is almost certainly an
error regardless of permissiveness level.

**Loosest/looser: perclimit disabled (0) for wt > 80 kg:**
In the code, "disabled" is implemented as perclimit = 0,
meaning the check `percewma < perclimit` is never TRUE
(no ratio can be < 0). This is distinct from Inf, which
would incorrectly trigger on every observation.
For weights above 80 kg, UW-adjusted wtallow always
catches any value that would also be caught by
perclimit = 0.4. The mathematical reason: perclimit =
0.4 fires when dewma > 1.5 × wt. For wt > 80, this
requires dewma > 120. The UW-adjusted wtallow catches
these values first because wtallow scales with UW (the
upper weight of the pair) while the perclimit threshold
grows only with the observation's own weight. Above
80 kg the coverage gap vanishes entirely.

**Tighter/tightest: perclimit 0.4 for all wt > 45 kg:**
With PW-L or allofus15 wtallow (not scaled up for high
UW), the base caps are lower (53.33 for PW-L, 40 for
allofus15). wtallow still catches values that
perclimit = 0.4 would catch for wt > 45.
However, the 0.4 threshold is retained at tighter and
tightest as a belt-and-suspenders safeguard with
negligible cost.

---

## Error Load Threshold

The error load threshold (Step 14) controls when a
subject's remaining included values are excluded because
the ratio of errors to total non-SDE values exceeds the
threshold. A lower threshold means the algorithm is more
aggressive about escalating when multiple errors are
detected. The tolerance of the error load threshold addresses a different aspect of permissiveness than most other permissiveness thresholds, but follows the same pattern of more values being excluded in tighter/tightest.

| Level | error_load_threshold | if 10 obs, # errors needed to trigger |
|-------|---------------------| ------------|
| `loosest` | >0.41 | 5/10 |
| `looser` | >0.41 | 5/10 |
| `tighter` | >0.29 | 3/10 |
| `tightest` | >0.29 | 3/10 |

---

## Moderate EWMA Directional Factor

The `f` factor in the Moderate EWMA step (Step 11Wb)
controls how strict the directional EWMA check is. The
standard pathway requires `dewma_before > f * wtallow`
and `dewma_after > f * wtallow` (in addition to
`dewma_all > wtallow`). A lower `f` value makes the
directional check less strict, leading to more
exclusions.

| Level | mod_ewma_f |
|-------|-----------|
| `loosest` | 0.75 |
| `looser` | 0.75 |
| `tighter` | 0.60 |
| `tightest` | 0.60 |

---

## Height Band

| Level | ht_band |
|-------|---------|
| `loosest` | 3" |
| `looser` | 3" |
| `tighter` | 2" |
| `tightest` | 2" |

## Height Loss and Height Gain

| Level | htloss allowed | htgain allowed |
|-------|---------------|----------------|
| `loosest` | Yes | Yes |
| `looser` | No | Yes |
| `tighter` | No | Yes |
| `tightest` | No | No |

Note: htloss is not allowed in `looser` because allowing
height loss introduces a large number of additional
inclusions that are more likely to represent measurement
error than true height loss. Height gain is allowed in
`looser` and `tighter` because young adults may still be
growing, and measurement error in height is common.
Height gain rescue applies only when the earlier
measurement is at age ≤25 years.

---

## Repeated Value (RV) Handling

RVs are values identical to a prior value for the same
subject and parameter (on a different day). Two modes are
available:

### `independent` (default for `loosest` and `looser`)
RVs are treated as fully informative values:
- RVs contribute to EWMA calculations normally
- RVs are evaluated against that EWMA like any other
  value
- RVs are no more likely to be excluded than non-RV
  values

### `linked` (default for `tighter` and `tightest`)
RVs are treated as less informative values:
- RVs are excluded from EWMA calculations in all EWMA
  steps (extreme, moderate, firstRV, allRV, and any
  other EWMA-based steps)
- RVs do not have their own EWMA in firstRV; RV
  propagation serves as their exclusion criterion instead
- In allRV and other EWMA steps, RVs are evaluated
  against an EWMA calculated without them, making them
  more likely to be flagged
- RVs are also more likely to cause exclusion of
  surrounding non-RV values, because removing RVs from
  the EWMA can expose patterns that would otherwise be
  masked

### Downstream effects of `linked` mode

Because RVs are ghosted from the EWMA, their exclusion
(or survival) can affect how surrounding values are
evaluated. For example:

  Patient weights (close together): 60 100 95 60 60 60 70

- Under `independent`: the 60s contribute to the EWMA,
  pulling it toward 60; the 70 looks fine relative to
  that EWMA and survives.
- Under `linked`: the 60s are excluded from the EWMA;
  the 70 is evaluated against an EWMA that does not
  include them and may be flagged depending on the
  surrounding values.

This means `linked` mode can result in more exclusions
of non-RV values in addition to the RVs themselves.
Users choosing `tighter` or `tightest` should be aware
that this may affect subjects with runs of repeated
values followed by a return toward a prior level.

### RV handling in other steps

The `repval_handling` setting affects all steps where
RVs are evaluated, including:
- 2D non-ordered: treated like 3D under `linked`
- 3D non-distinct values: given same weight as distinct
  values under `linked`
- Evil Twins (ET): RVs no longer treated differently
  under `linked`

---

## User Overrides

The `permissiveness` parameter sets defaults for all
sub-parameters. Individual sub-parameters can be
overridden by passing them explicitly to `cleanadult()`.
Explicit values always take precedence over the
permissiveness preset.

For example, a user could select
`permissiveness = "tighter"` but override
`allow_ht_loss = TRUE` to get the tighter thresholds
while still allowing height loss.

Also: `wtallow_formula` accepts a custom `.csv` file
path if users need a fully custom weight allowance table.

---

## Summary Table

**Current default is looser**

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| BIV thresholds | current defaults | metric (wide) | metric (medium) | metric (narrow) |
| wtallow formula | PW-H (piecewise) | PW-H (piecewise) | PW-L (piecewise-lower) | allofus15 |
| UW scaling | UW-based (see wtallow-formulas.md) | UW-based | UW-based (down only) | cap limited by PW-L |
| ET caps | wtallow cap + 20 | wtallow cap + 20 | wtallow cap + 20 | allofus15-cap-12m |
| perclimit ≤45 kg | 0.5 | 0.5 | 0.7 | 0.7 |
| perclimit 45–80 kg | 0.4 | 0.4 | 0.4 | 0.4 |
| perclimit >80 kg | 0 (disabled) | 0 (disabled) | 0.4 | 0.4 |
| error_load_threshold | 0.41 | 0.41 | 0.29 | 0.29 |
| mod_ewma_f | 0.75 | 0.75 | 0.60 | 0.60 |
| ht_band | 3" | 3" | 2" | 2" |
| htloss | allowed | not allowed | not allowed | not allowed |
| htgain | allowed | allowed | allowed | not allowed |
| repval_handling | independent | independent | linked | linked |
