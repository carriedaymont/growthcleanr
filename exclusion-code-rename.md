# Exclusion Code Rename — Child to Adult Format

## Child Exclusion Codes (child algorithm only)

Current codes actively assigned in `child_clean.R` (the child algorithm path). Legacy-only codes (`pediatric_clean_legacy.R`) are not included — they will keep their current names.

| Current Child Code | New Long Code | Note |
|---|---|---|
| Include | Include | |
| Missing | Exclude-Missing | Assigned in `cleangrowth()` wrapper, not child/adult algorithm. Applies to both. |
| Not cleaned | Exclude-Not-Cleaned | |
| Exclude-Temporary-Extraneous-Same-Day | Exclude-C-Temp-Same-Day | |
| Exclude-Carried-Forward | Exclude-C-CF | |
| Exclude-EWMA1-Extreme | Exclude-C-Traj-Extreme | |
| Exclude-SDE-Identical | Exclude-C-Identical | |
| Exclude-SDE-All-Exclude | | No longer exists (dead code) |
| Exclude-SDE-All-Extreme | Exclude-C-Extraneous | |
| Exclude-SDE-EWMA | Exclude-C-Extraneous | Collapsed with SDE-All-Extreme and SDE-One-Day |
| Exclude-SDE-One-Day | Exclude-C-Extraneous | Collapsed with SDE-All-Extreme and SDE-EWMA |
| Exclude-EWMA2-middle | Exclude-C-Traj | |
| Exclude-EWMA2-birth-WT | Exclude-C-Traj | |
| Exclude-EWMA2-birth-WT-ext | Exclude-C-Traj | |
| Exclude-EWMA2-first | Exclude-C-Traj | |
| Exclude-EWMA2-first-ext | Exclude-C-Traj | |
| Exclude-EWMA2-last | Exclude-C-Traj | |
| Exclude-EWMA2-last-high | Exclude-C-Traj | |
| Exclude-EWMA2-last-ext | Exclude-C-Traj | |
| Exclude-EWMA2-last-ext-high | Exclude-C-Traj | |
| Exclude-EWMA2-birth-HT-HC | Exclude-C-Traj | WT not applicable for birth HT/HC step |
| Exclude-EWMA2-birth-HT-HC-ext | Exclude-C-Traj | WT not applicable for birth HT/HC step |
| Exclude-Min-diff | Exclude-C-Abs-Diff | Collapsed with Max-diff |
| Exclude-Max-diff | Exclude-C-Abs-Diff | Collapsed with Min-diff |
| Exclude-2-meas->1-year | Exclude-C-Pair | |
| Exclude-2-meas-<1-year | Exclude-C-Pair | |
| Exclude-1-meas | Exclude-C-Single | |
| Exclude-Error-load | Exclude-C-Too-Many-Errors | |
| Exclude-Absolute-BIV | Exclude-C-BIV | Collapsed with Standardized-BIV |
| Exclude-Standardized-BIV | Exclude-C-BIV | Collapsed with Absolute-BIV |
| Exclude-Evil-Twins | Exclude-C-Evil-Twins | |
| Exclude-SDE-EWMA-All-Extreme | | No longer exists (dead code) |
| Exclude-1-CF-deltaZ-<0.05 | Rescued | cf_rescued column only; main column = Include |
| Exclude-1-CF-deltaZ-<0.1-wholehalfimp | Rescued-Imperial | cf_rescued column only; main column = Include |
| Exclude-Teen-2-plus-CF-deltaZ-<0.05 | Rescued-Adol | cf_rescued column only; main column = Include |
| Exclude-Teen-2-plus-CF-deltaZ-<0.1-wholehalfimp | Rescued-Adol-Imperial | cf_rescued column only; main column = Include |

---

## Adult Exclusion Codes

Current codes from `adult_clean.R` and `adult_support.R`.

Changes: Add `Exclude-Missing` (assigned in `cleangrowth()` wrapper). Add parameter to Evil-Twins. Remove `-N` suffix from trajectory codes.

| Current Adult Code | New Adult Code | Note |
|---|---|---|
| Include | Include | |
| (Missing assigned by wrapper) | Exclude-Missing | Same code as child; assigned in `cleangrowth()` |
| Exclude-A-BIV | Exclude-A-BIV | Remove parameter |
| Exclude-A-BIV | Exclude-A-BIV | Remove parameter |
| Exclude-A-Scale-Max-Identical | Exclude-A-Scale-Max-Identical | Remove parameter |
| Exclude-A-Scale-Max | Exclude-A-Scale-Max | Remove parameter |
| Exclude-A-Scale-Max-RV-Propagated | Exclude-A-Scale-Max-RV-Propagated | Remove parameter, fix case |
| Exclude-A-Evil-Twins | Exclude-A-Evil-Twins | No change |
| Exclude-A-Identical | Exclude-A-Identical | Remove parameter |
| Exclude-A-Extraneous | Exclude-A-Extraneous | Remove parameter |
| Exclude-A-Identical | Exclude-A-Identical | Remove parameter |
| Exclude-A-Extraneous | Exclude-A-Extraneous | Remove parameter |
| Exclude-A-Traj-Ext-N | Exclude-A-Traj-Ext | Remove parameter and -N suffix |
| Exclude-A-Traj-Extreme-firstRV-N | Exclude-A-Traj-Extreme-firstRV | Remove parameter and -N suffix |
| Exclude-A-Traj-Extreme-allRV-N | Exclude-A-Traj-Extreme-allRV | Remove parameter and -N suffix |
| Exclude-A-Ord-Pair | Exclude-A-Ord-Pair | Remove parameter |
| Exclude-A-Ord-Pair-All | Exclude-A-Ord-Pair-All | Remove parameter |
| Exclude-A-Window | Exclude-A-Window | Remove parameter |
| Exclude-A-Window-All | Exclude-A-Window-All | Remove parameter |
| Exclude-A-2D-Ordered | Exclude-A-2D-Ordered | Remove parameter |
| Exclude-A-2D-Non-Ordered | Exclude-A-2D-Non-Ordered | Remove parameter |
| Exclude-A-Traj-Moderate-N | Exclude-A-Traj-Moderate | Remove parameter and -N suffix |
| Exclude-A-Traj-Moderate-allRV-N | Exclude-A-Traj-Moderate-allRV | Remove parameter and -N suffix |
| Exclude-A-Traj-Moderate-Error-Load-N | Exclude-A-Traj-Moderate-Error-Load | Remove parameter and -N suffix |
| Exclude-A-Traj-Moderate-Error-Load-RV-N | Exclude-A-Traj-Moderate-Error-Load-RV | Remove parameter and -N suffix |
| Exclude-A-Single | Exclude-A-Single | Remove parameter |
| Exclude-A-Single | Exclude-A-Single | Remove parameter |
| Exclude-A-Too-Many-Errors | Exclude-A-Too-Many-Errors | Remove parameter |
| Exclude-A-Too-Many-Errors | Exclude-A-Too-Many-Errors | Remove parameter |
