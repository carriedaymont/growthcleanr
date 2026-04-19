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

---

## Session 2 — 2026-04-19 — Step 5 (Temp SDE) + Early Step 13 (SDE-Identicals)

Walkthrough note: [walkthrough-todo-2026-04-19.md](walkthrough-todo-2026-04-19.md) — see "R-vs-R comparison — Session 2" section.

- [AJ8] Early Step 13: `keep_id` empty-group handling — current guards with explicit `length(incl_ids) == 0L` check; reference's `ifelse(..., min(...), max(...))` with `na.rm = TRUE` on an empty Include subset returns `Inf`/`-Inf` with an R warning. Behavior-equivalent downstream (`has_dup` is FALSE so filter can't fire), but reference leaks warnings on datasets with all-excluded same-day groups — Bug fix — closed (current already correct, no change needed) — `Infants_Main.R:2634-2636` vs `child_clean.R:2638-2643`. Pitfall: **NA / empty-set handling**.
- [AJ9] `.child_valid()` / `valid()`: reference assigns bare `"Missing"` and `"Not cleaned"` codes (no `Exclude-` prefix) at `Infants_Main.R:1207`, `:1212`, `:1387`, `:1391`. Reference's `!grepl("^Exclude", exclude)` returns TRUE for both, silently counting these rows as valid — most likely trigger is a subject with two HC measurements in the 3–5y window (both `"Not cleaned"`) where `extraneous.this.day` could re-label one as `Exclude-Temporary-Extraneous-Same-Day`. Current's 2026-04-14/16 code rename (`Missing` → `Exclude-Missing`, `Not cleaned` → `Exclude-Not-Cleaned`) fixes this as a side effect — Bug fix — closed (current already correct via rename; documented here so the bug-fix side effect is explicit) — `pediatric_support.R:21` vs `child_clean.R:5020`. Pitfalls: **Factor levels / exclusion codes** + **`.child_valid()` flag scope**.

**Session 2 status:** Both findings closed with no code change needed. AJ8 is a minor warning-avoidance improvement already in current. AJ9 documents a bug-fix side effect of the 2026-04-14/16 exclusion-code rename that was previously framed as cosmetic. No tests re-run — baseline unchanged at 63 / 48 / 28 / 41 / 13.
