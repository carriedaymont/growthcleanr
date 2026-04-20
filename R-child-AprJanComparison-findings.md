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

---

## Session 3 — 2026-04-19 — Step 6 (CF detection + rescue)

Walkthrough note: [walkthrough-todo-2026-04-19.md](walkthrough-todo-2026-04-19.md) — see "R-vs-R comparison — Session 3" section.

**Findings opened this session: 0.**

AJ## numbering does not advance. The non-intentional diffs surfaced during the diff are all housekeeping / cleanup with logic-equivalent behavior — logged under "Items NOT flagged (audit trail)" in the walkthrough note rather than as standalone findings (e.g., reference's dead `cf_string_length` variable; reference's no-op drop of never-created `is_eligible_include`; redundant `%% 1` OR branch in `wholehalfimp` for HEIGHTCM/HEADCM; Step 6 working columns persisting past Step 6 in reference; paste-based `not_single` collision risk isomorphic to Session 1's AJ4).

The big intentional change — CF rescue scheme redesign (lookup-table + `cf_rescue` + `cf_rescued` + `cf_detail`; rescued CFs → `Include` rather than different exclusion codes) — is documented in the procedure's "Known intentional changes" list and in `cf-rescue-thresholds.md`, so logged briefly and not analyzed further here.

**Session 3 status:** Session closed with 0 findings and no code change. Baseline unchanged at 63 / 48 / 28 / 41 / 13; no tests re-run. Next session candidate: **Session 4 — Step 7 (BIV) + Step 9 (Evil Twins)**.

---

## Session 4 — 2026-04-22 — Step 7 (BIV) + Step 9 (Evil Twins)

Walkthrough note: [walkthrough-todo-2026-04-22.md](walkthrough-todo-2026-04-22.md) — see "R-vs-R comparison — Session 4" section.

- [AJ10] Step 7: Absolute BIV weight age boundary — reference uses `agedays == 0` / `agedays != 0` for the 0.2 kg / 1 kg floor; current uses `agedays <= 365` / `agedays > 365` (extends 0.2 kg floor to entire first year, fixing exclusion of legitimate preterm weights 0.7–1.0 kg) — Intentional (other) — closed (confirmed via CLAUDE.md history 2026-03-16; no change needed) — `Infants_Main.R:3073–3075` vs `child_clean.R:3007–3010`. Pitfall: **Boundary changes**.
- [AJ11] Step 9: `any(start_df$otl, na.rm = TRUE)` vs reference `any(start_df$oob)` — current adds `na.rm = TRUE` preventing `if (NA)` error when `otl` contains NAs and no TRUE values — Bug fix — closed (current already correct, no change needed) — `Infants_Main.R:3187` vs `child_clean.R:3123`. Pitfall: **NA / empty-set handling**.

**Session 4 status:** 2 findings (AJ10 — Intentional (other); AJ11 — Bug fix). Both closed with no code change needed. Baseline unchanged at 63 / 48 / 28 / 41 / 13; no tests re-run. Next session candidate: **Session 5 — Child Step 11 (EWMA1)**.

---

## Session 5 — 2026-04-19 — Step 11 (EWMA1: Extreme EWMA)

Walkthrough note: [walkthrough-todo-2026-04-19.md](walkthrough-todo-2026-04-19.md) — see "R-vs-R comparison — Session 5" section.

- [AJ12] Step 11 / `ewma()` signature: EWMA window default changed from 25 (reference) to 15 (current `ewma_window` parameter); reference comment says 25 was "Changed … to 25 for better accuracy with minimal efficiency loss" — Intentional (other) — closed (confirmed intentional, no code change needed; all ewma() and ewma_cache_init() call sites already use `window = ewma_window` consistently) — `Infants_Main.R:2013–2015` vs `child_clean.R:1695`; call site `Infants_Main.R:3325` vs `child_clean.R:3285`. Pitfall: **Boundary changes**.

**Session 5 status:** 1 finding (AJ12 — Intentional (other)). Closed with no code change. Baseline unchanged at 63 / 48 / 28 / 41 / 13; no tests re-run.

---

## Session 6 — 2026-04-19 — Step 13 (Final SDE resolution, multi-phase B1/B2/B3)

Walkthrough note: [walkthrough-todo-2026-04-19.md](walkthrough-todo-2026-04-19.md) — see "R-vs-R comparison — Session 6" section.

- [AJ13] Step 13 Phase B3: `ewma_fill` and `spa_ewma` Inf guards — reference has no guard: when a `(subjid, param, agedays)` group has only temp SDE rows (stable Include rows excluded by BIV/EWMA1), `max(ewma.all[!was_temp_sde], na.rm = TRUE)` = `-Inf`, propagates to `spa_ewma = -Inf`, `absdewma = Inf`, `min_absdewma > 1` → all rows erroneously marked Extraneous; current adds `is.infinite → NA` conversions for both `ewma_fill` and `spa_ewma`, preventing erroneous marking and correctly keeping the best-tiebreaker row as Include — Bug fix — closed (current already correct, no change needed) — `Infants_Main.R:3722–3730` vs `child_clean.R:3593–3607`. Pitfall: **NA / empty-set handling**.

**Session 6 status:** 1 finding (AJ13 — Bug fix). Closed with no code change. Baseline unchanged at 63 / 48 / 28 / 41 / 13; no tests re-run. Next session candidate: **Session 7 — Child Step 15 (EWMA2) + Step 16 (Birth HT/HC)**.

---

## Session 7 — 2026-04-19 — Step 15 (EWMA2 Moderate) + Step 16 (Birth HT/HC)

Walkthrough note: [walkthrough-todo-2026-04-19.md](walkthrough-todo-2026-04-19.md) — see "R-vs-R comparison — Session 7" section.

**Findings opened this session: 0.**

AJ## numbering does not advance. All non-intentional diffs are known-intentional or equivalent cleanup — logged under "Items NOT flagged (audit trail)" in the walkthrough note. Key items: tbc_range pre-filter before the Step 15 while loop (safe efficiency optimization, range ≤ 1 groups cannot pass the `|dewma| > 1` addcrit threshold); DOP snapshot timing change (`dop_snap` keyed at iteration start vs. reference's live per-closure `data.df` scan — current is strictly more correct, eliminates implicit within-iteration ordering dependency); `include.temporary.extraneous = TRUE` removal from `include_counts` (Session 11 F74, confirmed intentional); sp_key-level Step 16 birth filter (Session 11 F85, confirmed intentional behavior-neutral tightening); EWMA caching (`ewma_cache_init`/`ewma_cache_update`/`ewma2_caches` — known intentional performance feature); `cols_to_drop_15_16` scope (Session 11 F81, drops p_plus/p_minus/tbc.p_plus/tbc.p_minus/first_meas after Step 16); all rounding-tolerance removals (per procedure, do NOT flag).

**Session 7 status:** 0 findings. Closed with no code change. Baseline unchanged at 63 / 48 / 28 / 41 / 13; no tests re-run. Next session candidate: **Session 8 — Child Step 17 (Height/HC Velocity)**.

---

## Session 8 — 2026-04-22 — Step 17 (Height/HC Velocity)

Walkthrough note: [walkthrough-todo-2026-04-22.md](walkthrough-todo-2026-04-22.md) — see "R-vs-R comparison — Session 8" section.

- [AJ14] Step 17 HC velocity: reference HC branch is effectively a no-op — for-loop condition `!is.na(df$whoagegrp.hc)` always FALSE because `whoagegrp.hc` is never set (reference sets `whoagegrp.ht` instead); loop body references `sub_m_who_ht_vel` (HT variable, undefined in HEADCM context) and `.ht`-suffix column names; result: only `mindiff = -1.5` fallback applied, no upper bound check. Current correctly implements HC velocity via pre-merged `.hc`-suffix columns and `fcase` lookup — Bug fix — closed (current already correct, no change needed) — `Infants_Main.R:4452–4523` vs `child_clean.R:4557–4609`. Pitfalls: **NA / empty-set handling** + **Parameter scope**.

**Session 8 status:** 1 finding (AJ14 — Bug fix). Closed with no code change. Baseline unchanged at 63 / 48 / 28 / 41 / 13; no tests re-run. Next session candidate: **Session 9 — Child Step 19 (Pairs/Singles) + Step 21 (Error Load) + Step 22 (Output)**.
