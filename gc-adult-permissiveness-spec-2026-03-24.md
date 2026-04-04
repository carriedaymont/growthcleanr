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

### Base formulas

| Level | wtallow formula | scaled for high weight |
|-------|----------------------|
| `loosest` | piecewise (current default) | scaled as below |
| `looser` | piecewise (current default) | scaled as below |
| `tighter` | piecewise-lower | not scaled |
| `tightest` | allofus15 (see below) | not scaled |

allofus 15 definition
1–2 days: 5 kg
3–7 days: 10 kg
8d to <6m: 15 kg
>=6m: linear increase from 15kg @ 6m to 40kg at 12m
>=12m: 40kg

### Weight-scaled wtallow (loosest and looser only)

For patients above 120 kg, the base wtallow is augmented
with a weight-scaled addition:

```
wtallow = base_formula(interval)
          + 0.50 * max(0, max(wt, ewma) - 120)
```

This accommodates the larger absolute weight changes that
are plausible for very heavy patients (e.g., bariatric
surgery, GLP-1 medications, fluid shifts). The 120 kg
threshold corresponds to BMI ~40–44 (class III obesity),
where bariatric surgery is standard of care and large
absolute weight changes become expected.

The 0.50 scaling factor is lower than the 0.70 factor
used for EWMA caps (see below), ensuring wtallow stays
below the ET cap.

Weight scaling is not applied at tighter or tightest
because those levels intentionally exclude extreme weight
changes associated with bariatric surgery, severe
cachexia, and rapid fluid status changes.

### Effective wtallow ceiling

At all permissiveness levels, the EWMA cap serves as a
ceiling on wtallow. That is:

```
effective_wtallow = min(wtallow, ewma_cap)
```

This is already implemented in the code. It ensures that
wtallow never exceeds the ET/extreme EWMA cap, which
would create a zone where moderate EWMA could flag a
value that evil twins and extreme EWMA would miss.

### wtallow tables

**Loosest** (piecewise + weight scaling):
```
Months   Base    @80    @120    @160    @200
  1      20.0   20.0    20.0    40.0    60.0
  3      32.0   32.0    32.0    52.0    72.0
  6      50.0   50.0    50.0    70.0    90.0
 12      80.0   80.0    80.0    100.0   120.0
 18      80.0   80.0    80.0    100.0   120.0
 24      80.0   80.0    80.0    100.0   120.0
```


---

## EWMA Caps

EWMA caps control the thresholds used by Evil Twins
(Step 9Wa) and Extreme EWMA (Step 9Wb). EWMA caps also
serve as the ceiling for wtallow (see above).

### Loosest and Looser: Weight-scaled caps

A 2-tier structure (<6m, ≥6m) with weight-dependent
scaling for heavy patients:

```
et_limit = baseline + 0.70 * max(0, maxwt - 120)

where baseline = 50  (if interval < 6 months)
                 80  (if interval >= 6 months)
      maxwt    = max(wt, ewma)
```

For patients ≤120 kg, this gives fixed caps of 50/80.
For heavier patients, caps scale:

```
MaxWt    <6m     ≥6m
 60 kg   50.0    80.0
100 kg   50.0    80.0
120 kg   50.0    80.0
140 kg   64.0    94.0
160 kg   78.0   108.0
180 kg   92.0   122.0
200 kg  106.0   136.0
```

Both loosest and looser use the same cap formula. The
differentiation between loosest and looser is handled
by wtallow and other parameters.

**Rationale for 2-tier structure:** The split at 6 months
reflects that peak weight loss from surgery and medication
typically occurs within the first year. A higher cap at
≥6m (80 vs 50) accommodates this period of large real
changes while keeping the short-interval cap tighter for
normal-weight patients.

### Tighter and Tightest: Fixed 2-tier caps

Fixed caps without weight scaling:

| Level | <6m cap | ≥6m cap |
|-------|---------|---------|
| `tighter` | 40 | 60 |
| `tightest` | 40 | 40 |

For tighter, the ≥6m cap (60) is higher than the <6m
cap (40). This is intentional: peak bariatric weight
loss occurs at 9–15 months, and a higher cap at longer
intervals accommodates patients with moderate obesity
who have real weight changes in this period. At the
tighter level, the goal is to exclude values reflecting
severe illness, not moderate-to-large real changes.

For tightest, both tiers are 40 kg.

No weight scaling is applied at tighter or tightest
because these levels intentionally do not accommodate
extreme bariatric-surgery-scale changes.

### Gap between wtallow and EWMA cap

For the system to work correctly, wtallow must remain
below the EWMA cap. At the extreme (loosest, 200 kg):

```
Months  wtallow  ET cap  Gap
  1      60.0    106.0   46.0
  3      76.0    106.0   30.0
  6     100.0    136.0   36.0
 12     110.0    136.0   26.0
 18     120.0    136.0   16.0
 24     120.0    136.0   16.0
```

The gap narrows to 16 kg at long intervals for very
heavy patients. This is adequate — if a 200 kg patient
has a dewma of 120–136 kg, they are flagged as moderate
(Step 11Wb) but not as extreme (Step 9). Values above
136 kg are flagged as extreme.

If a custom formula is used that grows past 80 at long
intervals, the EWMA cap ceiling prevents wtallow from
exceeding the ET cap.

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
For weights above 80 kg, wtallow always catches any
value that would also be caught by perclimit = 0.4.
The mathematical reason: perclimit = 0.4 fires when
dewma > 1.5 × wt. For wt > 80, this requires
dewma > 120. Even with weight-scaled wtallow, the
formula (base + 0.50 × (maxwt − 120)) catches these
values first because the wtallow scaling grows with
maxwt while the perclimit threshold grows only with wt.
The gap is at most a few kg wide in a narrow band for
wt 60–80 with very high EWMAs, and in those edge cases
the directional EWMA checks and trajectory criteria
provide additional redundant coverage. Above 80 kg the
gap vanishes entirely.

**Tighter/tightest: perclimit 0.4 for all wt > 45 kg:**
Without weight-scaled wtallow, the base wtallow caps
at 80 (piecewise) or 60 (piecewise-lower) or 40
(tightest). In this case, wtallow always catches
values that perclimit = 0.4 would catch for wt > 45.
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
| wtallow base | piecewise | piecewise | piecewise-lower | allofus15 |
| wtallow scaling | + 0.50×(maxwt−120)⁺ | + 0.50×(maxwt−120)⁺ | none | none |
| ewma_cap <6m | 50 + 0.70×(maxwt−120)⁺ | 50 + 0.70×(maxwt−120)⁺ | 40 | 40 |
| ewma_cap ≥6m | 80 + 0.70×(maxwt−120)⁺ | 80 + 0.70×(maxwt−120)⁺ | 60 | 40 |
| perclimit ≤45 kg | 0.5 | 0.5 | 0.7 | 0.7 |
| perclimit 45–80 kg | 0.4 | 0.4 | 0.4 | 0.4 |
| perclimit >80 kg | 0 (disabled) | 0 (disabled) | 0.4 | 0.4 |
| error_load_threshold | 0.41 | 0.41 | 0.29 | 0.29 |
| mod_ewma_f | 0.75 | 0.75 | 0.60 | 0.60 |
| ht_band | 3" | 3" | 2" | 2" |
| htloss | allowed | not allowed | not allowed | not allowed |
| htgain | allowed | allowed | allowed | not allowed |
| repval_handling | independent | independent | linked | linked |
