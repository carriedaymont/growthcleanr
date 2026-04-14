# CF Rescue Threshold Lookup Tables

**Date:** 2026-04-14
**Source:** `__Pipeline/CF-exploration/cf-threshold-schemes.md` (full methodology and rationale)

## Overview

CF rescue decides whether an identical-to-prior measurement is a true coincidence (I-T) or a carried-forward error (I-CF). The threshold T rescues (includes) values with |deltaZ| < T and excludes values with |deltaZ| >= T.

These thresholds replace the current fixed `cf_rescue_threshold` (0.05 standard, 0.10 imperial). They are looked up by **age bin**, **interval bin**, **param**, and **rounding type** (imperial vs other).

Three threshold levels: **0.05** (slow growth — short intervals or older ages), **0.20** (moderate), **0.40** (fast growth — longer intervals or younger ages). **NR** = no rescue (all identical values excluded). **--** = cell cannot occur.

---

## Lookup Dimensions

### Age bins (8)

| Bin | agedays range |
|-----|--------------|
| 0-3mo | 0 – 90 |
| 3-6mo | 91 – 182 |
| 6-12mo | 183 – 365 |
| 1-2y | 366 – 730 |
| 2-5y | 731 – 1826 |
| 5-10y | 1827 – 3652 |
| 10-15y | 3653 – 5478 |
| 15-20y | 5479+ |

Age bin is determined by the **current** measurement's agedays (the one being evaluated for CF).

### Interval bins (5)

| Bin | interval_days range |
|-----|-------------------|
| <1wk | 1 – 6 |
| 1wk-1mo | 7 – 29 |
| 1-6mo | 30 – 182 |
| 6mo-1y | 183 – 364 |
| >1y | 365+ |

Interval = agedays of current measurement minus agedays of the prior measurement (the one the current value is identical to).

### Rounding type (2)

- **imperial:** Rounding precision >= 1 inch (2.54 cm) for HT or >= 0.5 lb (0.227 kg) for WT. Applied only at age > 2y.
- **other:** All other rounding (metric, nearest-cm, etc.) and all ages <= 2y regardless of rounding.

Detection method: existing `rounding_rule` logic in growthcleanr (examines digit patterns).

---

## Threshold Tables

### HEIGHTCM — Other

| Age | <1wk | 1wk-1mo | 1-6mo | 6mo-1y | >1y |
|-----|:----:|:-------:|:-----:|:------:|:---:|
| 0-3mo | 0.40 | 0.40 | -- | -- | -- |
| 3-6mo | 0.20 | 0.40 | -- | -- | -- |
| 6-12mo | 0.05 | 0.40 | NR | -- | -- |
| 1-2y | 0.05 | 0.40 | 0.40 | -- | -- |
| 2-5y | 0.05 | 0.20 | 0.40 | -- | -- |
| 5-10y | 0.05 | 0.05 | 0.40 | 0.40 | -- |
| 10-15y | 0.05 | 0.05 | 0.20 | 0.40 | NR |
| 15-20y | 0.05 | 0.05 | 0.05 | 0.20 | 0.20 |

### HEIGHTCM — Imperial (>=2y only)

| Age | <1wk | 1wk-1mo | 1-6mo | 6mo-1y | >1y |
|-----|:----:|:-------:|:-----:|:------:|:---:|
| 2-5y | 0.05 | 0.20 | 0.40 | 0.40 | -- |
| 5-10y | 0.05 | 0.05 | 0.40 | 0.40 | -- |
| 10-15y | 0.05 | 0.05 | 0.20 | 0.40 | NR |
| 15-20y | 0.05 | 0.05 | 0.05 | 0.20 | 0.20 |

### WEIGHTKG — Other

| Age | <1wk | 1wk-1mo | 1-6mo | 6mo-1y | >1y |
|-----|:----:|:-------:|:-----:|:------:|:---:|
| 0-3mo | 0.40 | 0.40 | -- | -- | -- |
| 3-6mo | 0.20 | 0.40 | NR | -- | -- |
| 6-12mo | 0.05 | 0.20 | NR | -- | -- |
| 1-2y | 0.05 | 0.20 | 0.40 | -- | -- |
| 2-5y | 0.05 | 0.05 | 0.40 | 0.40 | -- |
| 5-10y | 0.05 | 0.05 | 0.20 | 0.40 | -- |
| 10-15y | 0.05 | 0.05 | 0.20 | 0.40 | NR |
| 15-20y | 0.05 | 0.05 | 0.20 | 0.20 | 0.20 |

### WEIGHTKG — Imperial (>=2y only)

| Age | <1wk | 1wk-1mo | 1-6mo | 6mo-1y | >1y |
|-----|:----:|:-------:|:-----:|:------:|:---:|
| 2-5y | 0.05 | 0.05 | 0.40 | 0.40 | -- |
| 5-10y | 0.05 | 0.05 | 0.20 | 0.40 | -- |
| 10-15y | 0.05 | 0.05 | 0.20 | 0.40 | NR |
| 15-20y | 0.05 | 0.05 | 0.20 | 0.20 | 0.20 |

---

## NR Rules

NR (no rescue) applies when:
1. Starting age < 12 months AND interval > 1 month
2. Starting age < 15 years AND interval > 1 year

Rationale: children growing fast enough that coincidental identical measurements essentially don't occur at these intervals.

## Implementation Notes

- **NR means T = 0**: all CFinit values in that cell are excluded (no deltaZ rescue attempted).
- **"--" cells** should not occur in practice (e.g., a 0-3mo child cannot have a 6mo interval). If encountered, treat as NR.
- **HEADCM:** Not yet determined. For now, use HEIGHTCM thresholds as a placeholder (HC growth rate is between HT and WT).
- **String length:** These thresholds are for CF-2 (second in a string of identical values). CF-3+ rules are not yet defined; current behavior (rescue with same threshold) can continue for now.
- **Rounding detection:** Use existing `rounding_rule` logic. If rounding type is unknown, default to "other" (more conservative — fewer false exclusions).
- **Parameter for backward compatibility:** The existing `cf_rescue_threshold` parameter could be retained as an override. When set to its default, use the lookup table. When explicitly set by user, use the uniform value.
