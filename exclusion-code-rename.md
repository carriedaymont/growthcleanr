# Exclusion Code Rename — Child to Adult Format

## Child Exclusion Codes (child algorithm only)

Current codes actively assigned in `child_clean.R` (the child algorithm path). Legacy-only codes (`pediatric_clean_legacy.R`) are not included — they will keep their current names.

Note: "-WT/HT/HC-" indicates there will be separate codes for WT/HT/HC.

| Current Child Code | New Long Code | Note |
|---|---|---|
| Include | Include | |
| Missing | Exclude-Missing | Assigned in `cleangrowth()` wrapper, not child/adult algorithm. Applies to both. |
| Not cleaned | Exclude-Not-Cleaned | |
| Exclude-Temporary-Extraneous-Same-Day | Exclude-C-Temp-Same-Day | |
| Exclude-Carried-Forward | Exclude-C-WT/HT/HC-CF | |
| Exclude-EWMA1-Extreme | Exclude-C-WT/HT/HC-Traj-Extreme | |
| Exclude-SDE-Identical | Exclude-C-WT/HT/HC-Identical | |
| Exclude-SDE-All-Exclude | | No longer exists (dead code) |
| Exclude-SDE-All-Extreme | Exclude-C-WT/HT/HC-Extraneous | |
| Exclude-SDE-EWMA | Exclude-C-WT/HT/HC-Extraneous | Collapsed with SDE-All-Extreme and SDE-One-Day |
| Exclude-SDE-One-Day | Exclude-C-WT/HT/HC-Extraneous | Collapsed with SDE-All-Extreme and SDE-EWMA |
| Exclude-EWMA2-middle | Exclude-C-WT/HT/HC-Traj | |
| Exclude-EWMA2-birth-WT | Exclude-C-WT-Traj | |
| Exclude-EWMA2-birth-WT-ext | Exclude-C-WT-Traj | |
| Exclude-EWMA2-first | Exclude-C-WT/HT/HC-Traj | |
| Exclude-EWMA2-first-ext | Exclude-C-WT/HT/HC-Traj | |
| Exclude-EWMA2-last | Exclude-C-WT/HT/HC-Traj | |
| Exclude-EWMA2-last-high | Exclude-C-WT/HT/HC-Traj | |
| Exclude-EWMA2-last-ext | Exclude-C-WT/HT/HC-Traj | |
| Exclude-EWMA2-last-ext-high | Exclude-C-WT/HT/HC-Traj | |
| Exclude-EWMA2-birth-HT-HC | Exclude-C-HT/HC-Traj | WT not applicable for birth HT/HC step |
| Exclude-EWMA2-birth-HT-HC-ext | Exclude-C-HT/HC-Traj | WT not applicable for birth HT/HC step |
| Exclude-Min-diff | Exclude-C-HT/HC-Abs-Diff | Collapsed with Max-diff |
| Exclude-Max-diff | Exclude-C-HT/HC-Abs-Diff | Collapsed with Min-diff |
| Exclude-2-meas->1-year | Exclude-C-WT/HT/HC-Pair | |
| Exclude-2-meas-<1-year | Exclude-C-WT/HT/HC-Pair | |
| Exclude-1-meas | Exclude-C-WT/HT/HC-Single | |
| Exclude-Error-load | Exclude-C-WT/HT/HC-Too-Many-Errors | |
| Exclude-Absolute-BIV | Exclude-C-WT/HT/HC-BIV | Collapsed with Standardized-BIV |
| Exclude-Standardized-BIV | Exclude-C-WT/HT/HC-BIV | Collapsed with Absolute-BIV |
| Exclude-Evil-Twins | Exclude-C-WT/HT/HC-Evil-Twins | |
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
| Exclude-A-HT-BIV | Exclude-A-HT-BIV | No change |
| Exclude-A-WT-BIV | Exclude-A-WT-BIV | No change |
| Exclude-A-WT-Scale-Max-Identical | Exclude-A-WT-Scale-Max-Identical | No change |
| Exclude-A-WT-Scale-Max | Exclude-A-WT-Scale-Max | No change |
| Exclude-A-WT-Scale-Max-RV-Propagated | Exclude-A-WT-Scale-Max-RV-Propagated | No change |
| Exclude-A-Evil-Twins | Exclude-A-WT-Evil-Twins | Add parameter |
| Exclude-A-HT-Identical | Exclude-A-HT-Identical | No change |
| Exclude-A-HT-Extraneous | Exclude-A-HT-Extraneous | No change |
| Exclude-A-WT-Identical | Exclude-A-WT-Identical | No change |
| Exclude-A-WT-Extraneous | Exclude-A-WT-Extraneous | No change |
| Exclude-A-WT-Traj-Ext-N | Exclude-A-WT-Traj-Ext | Remove -N suffix |
| Exclude-A-WT-Traj-Extreme-firstRV-N | Exclude-A-WT-Traj-Extreme-firstRV | Remove -N suffix |
| Exclude-A-WT-Traj-Extreme-allRV-N | Exclude-A-WT-Traj-Extreme-allRV | Remove -N suffix |
| Exclude-A-HT-Ord-Pair | Exclude-A-HT-Ord-Pair | No change |
| Exclude-A-HT-Ord-Pair-All | Exclude-A-HT-Ord-Pair-All | No change |
| Exclude-A-HT-Window | Exclude-A-HT-Window | No change |
| Exclude-A-HT-Window-All | Exclude-A-HT-Window-All | No change |
| Exclude-A-WT-2D-Ordered | Exclude-A-WT-2D-Ordered | No change |
| Exclude-A-WT-2D-Non-Ordered | Exclude-A-WT-2D-Non-Ordered | No change |
| Exclude-A-WT-Traj-Moderate-N | Exclude-A-WT-Traj-Moderate | Remove -N suffix |
| Exclude-A-WT-Traj-Moderate-allRV-N | Exclude-A-WT-Traj-Moderate-allRV | Remove -N suffix |
| Exclude-A-WT-Traj-Moderate-Error-Load-N | Exclude-A-WT-Traj-Moderate-Error-Load | Remove -N suffix |
| Exclude-A-WT-Traj-Moderate-Error-Load-RV-N | Exclude-A-WT-Traj-Moderate-Error-Load-RV | Remove -N suffix |
| Exclude-A-HT-Single | Exclude-A-HT-Single | No change |
| Exclude-A-WT-Single | Exclude-A-WT-Single | No change |
| Exclude-A-HT-Too-Many-Errors | Exclude-A-HT-Too-Many-Errors | No change |
| Exclude-A-WT-Too-Many-Errors | Exclude-A-WT-Too-Many-Errors | No change |
