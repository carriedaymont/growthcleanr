# R Child April-vs-January Comparison — Findings Index

Cross-session findings from `R-child-AprJanComparison-procedure.md`. One-liner per finding.

**Categories:** Bug fix · New error · Unclear · Intentional (other)
**Status:** open → approved → fixed → closed

---

## Session 1 — 2026-04-19 — Preprocessing (z-scores, GA correction, recentering)

Walkthrough note: [walkthrough-todo-2026-04-19.md](walkthrough-todo-2026-04-19.md) — see "R-vs-R comparison — Session 1" section.

- [AJ1] Preprocessing input: `"HEIGHIN"` typo in param validation list silently rendered HEIGHTIN data Missing — Bug fix — closed (current already correct, no change needed) — `Infants_Main.R:389-399` vs `child_clean.R:506-516`
- [AJ2] Step 2b: `!is.na(sd.orig)` guard added to potcorr_wt computation; reference would propagate NA into potcorr_wt for Missing first-WT rows — Bug fix — closed (current already correct, no change needed) — `Infants_Main.R:801-803` vs `child_clean.R:785-789`
- [AJ3] Step 2b: `intwt` low-weight floor bounds changed (250-560g→570g vs 100-499g→500g); lower bound widened, upper bound shifted from inclusive to exclusive — Intentional (other) — closed (confirmed intentional, no change needed) — `Infants_Main.R:821` vs `child_clean.R:822`
- [AJ4] Step 2b: `tmp[subjid %in% names(table(subjid) > 1),]` filter was effectively no-op in reference (`names()` of logical vector returns all names); current uses `which()` to filter correctly — Bug fix — closed (current already correct, no change needed) — `Infants_Main.R:1122` vs `child_clean.R:984`
- [AJ5] Step 2b: `c(orig_colnames, id)` ambiguous (variable vs string) in column-keep filter; current uses `c(orig_colnames, "id")` — Bug fix (cleanup) — closed (current already correct, no change needed) — `Infants_Main.R:1179-1180` vs `child_clean.R:1030-1031`
- [AJ6] Step 2b: At exactly agedays=1461 (4 years) for potcorr subjects, reference assigns sd.corr := sd.orig (override after smooth blend); current was assigning sd.corr := sd.c (smooth blend wins, no override). Resolved by tightening second `fcase()` branch from `ageyears_2b <= 4` to `ageyears_2b < 4` so the smooth-blend window is half-open `[2, 4)` and age=4 falls to default sd.orig — Bug fix — closed (fixed in current) — `child_clean.R:933-941`
- [AJ7] Recentering: HC ≥ 5y exclusion code changed from `'Missing'` to `'Exclude-Not-Cleaned'` (semantic recategorization, beyond rename) — Intentional (other) — closed (confirmed intentional, no change needed) — `Infants_Main.R:1391` vs `child_clean.R:1121`

**Session 1 status:** All 7 findings closed. AJ6 fixed in current code. AJ1/AJ2/AJ4/AJ5 confirmed correct in current as-is (bug fixes already present). AJ3/AJ7 confirmed intentional, no change needed.
