# GC Walk-Through — 2026-04-18b: Child Step 13 (Final SDE Resolution)

Dedicated detailed walkthrough of Main Child Step 13 (Final SDE Resolution) in `R/child_clean.R` (~lines 3361–3653 pre-edit). Adult algorithm has its own separate SDE resolution path in `adult_clean.R` / `adult_support.R` (Adult Steps 3, 9H, 10W) and is out of scope; this pass is child-scoped. Early Child Step 13 was already walked in `walkthrough-todo-2026-04-16d.md` Session 2 and is not re-walked here except where it affects Main Step 13 invariants.

Comment / narrative / small-cleanup only; no behavioral changes. Child test counts unchanged.

---

## Pre-session baseline tests

- test-cleangrowth.R: 63 PASS, 0 warnings
- test-child-regression.R: 48 PASS, 0 warnings
- test-child-edge-cases.R: 28 PASS, 0 warnings
- test-child-algorithms.R: 41 PASS (2 pre-existing codetools warnings from `parallel=TRUE` / `plyr::ddply` — baseline, unchanged)
- test-child-parameters.R: 13 PASS

Adult tests not run — no adult files touched.

---

## Scope

| Part | Lines (pre-edit) | What |
|---|---|---|
| Phase A: setup | `child_clean.R` 3361–3392 | Capture temp SDE ids, safety check, reset-to-Include + identify_temp_sde re-run, post safety check, `keep_cols_sde` capture |
| Pre-filter | 3393–3410 | Identify subjects with same-day groups; skip-branch if none |
| Construct `data.sde` | 3411–3425 | Subject-filter, `sde_this` flag, `.child_valid(..., include.temporary.extraneous = TRUE)` subset, subject-level `has_sde_subj` filter, sort |
| Phase B1: SDE-Identicals (defensive) | 3427–3448 | Whole-day `uniqueN(v) == 1L` check + partial `dup_count > 1L` check, both age-dependent |
| Restore temp SDEs | 3454–3457 | Capture `was_temp_sde`, reset TempSDE to Include, re-sort |
| Phase B2: One-Day SDE | 3459–3537 | `one_day_sde_flag`, median-based distance, DOP medians (`HT/WT/HC_dop_med`), SDE-All-Extreme `> 2`, one-day winner selection |
| Phase B3: SDE-EWMA | 3540–3632 | EWMA on non-temp-SDE rows, fill temp SDEs via `ewma_fill`, `spa_ewma`, `absdewma`, SDE-All-Extreme `> 1`, EWMA-based winner selection |
| Phase B4: Merge back | 3636–3653 | `sde_results` merge by `id`, Include-or-TempSDE-only overwrite, extra-col drop, `setkey` |

Inventory: `identify_temp_sde()` is called in Phase A (line 3382) with `exclude_from_dop_ids = temp_sde_ids_step13`. This is the only non-NULL caller across the 7 call sites (confirmed Session 8). `ewma()` is called in Phase B3 (line 3567). No other support functions are invoked from Step 13.

---

## Fix-now items

### F66. File header `identify_temp_sde()` bullet updated — FIXED
- **File/lines:** `R/child_clean.R:89` (pre-edit).
- **Issue.** `SUPPORTING FUNCTIONS` index bullet read `- identify_temp_sde(): Temporary SDE resolution (Step 5)`. Same stale-scope pattern fixed by F-bycatch2 in Session 9 for `calc_and_recenter_z_scores()`: Session 8 confirmed `identify_temp_sde()` is called from 7 sites, not just Child Step 5. The roxygen block at the function already lists all callers (F50).
- **Fix.** Bullet expanded to two lines: `Temporary SDE resolution — winner selection on each (subjid, param, agedays) group (Child Steps 5, 6, 7, 9, 11, 13)`. Step 11 is listed once even though there are two call sites (mid-loop and end-of-step) to match the other narrative conventions.

### F67. `# before dplyr chain` stale comment rewritten — FIXED
- **File/lines:** `R/child_clean.R:3412` (pre-edit).
- **Issue.** Comment read `# Filter to subjects with potential SDEs before dplyr chain`. The block immediately below uses only data.table (`data.df[subjid %in% ...]`, `data.df_sde_subset[, sde_this := ..., by = ...]`, `data.df_sde_subset[.child_valid(...)]`). No dplyr chain. Stale reference.
- **Fix.** Replaced with `# Filter to subjects with at least one same-day group before the per-row / # data.table work below.` Describes current behavior.

### F68. `# Include id for deterministic SDE order` comment clarified — FIXED
- **File/lines:** `R/child_clean.R:3450` (pre-edit).
- **Issue.** Comment said `Include id` but the `setkey(data.sde, subjid, param, agedays, internal_id)` on the next line uses `internal_id`, not the user's `id`. Loose wording — `id` is ambiguous now that user `id` is preserved separately from `internal_id`.
- **Fix.** Rewritten as: `Re-key so downstream per-row work uses internal_id as the within-SPA / tiebreaker.` (Line 3667's analogous comment at the start of Child Step 15/16 has the same loose wording; left alone here because it belongs to the Step 15/16 walkthrough.)

### F69. `# "Ba" value` stale Stata-era parenthetical removed — FIXED
- **File/lines:** `R/child_clean.R:3554` (pre-edit).
- **Issue.** Comment read `# Take the larger gap (the "Ba" value) and convert to years`. The `"Ba"` token is a Stata-era variable name that does not appear anywhere else in `child_clean.R` (grep confirms a single hit). A current R reader has no way to decode it.
- **Fix.** Rewritten as: `# Take the larger of the before/after gaps and convert to years; this drives / # the per-row EWMA exponent below (smaller gap -> less-negative exponent -> / # heavier weight on closer neighbors).` Drops the Stata token and states the downstream consequence.

### F70. Phase B1 narrative rewritten as defensive — FIXED
- **File/lines:** `child-algorithm-reference.md:535–542` (pre-edit).
- **Issue.** Pre-edit text read: `Runs again on data.sde to catch identicals that emerged after Child Steps 5–11 removed intervening values. Two checks: 1. Whole-day identical ... 2. Partial identical ... Both use the age-dependent id rule. Only Include rows are candidates (exclude == "Include" guard).`
  - The "emerged after Child Steps 5–11" framing is not accurate for the current code. After Early Child Step 13, no two same-day rows in `data.df` share the same `v` on their Include/TempSDE status; every later mechanism preserves that invariant:
    - Steps 5, 6, 7, 9, 11 never modify `v`.
    - Step 5 initial pass and Steps 6/7/9/11 post-step reruns of `identify_temp_sde()` preserve one-Include-per-SPA by construction — Step 5 selects a winner within each (subjid, param, agedays) group.
    - `cf_rescue = "all"` restored values sit on different days from their CF originators, so they do not reintroduce same-day same-value pairs.
    - The Step 11c and Step 13 pre-phase reruns at `child_clean.R:3358–3359` and `3381–3384` also preserve the same invariant.
  - The "Only Include rows are candidates (exclude == "Include" guard)" claim is only half-accurate: the whole-day check at line 3433–3435 has no `exclude == "Include"` guard and would demote a same-`v` TempSDE loser to `Exclude-C-Identical` if this block ever fired. Only the partial check at line 3446 has that guard.
- **Fix.** Heading changed to `**Phase B1: SDE-Identicals (defensive)**`. Prose rewritten to describe the block as a no-op safety net, walking through the invariant chain above (Early 13 + `v` not modified + every rerun preserves one-Include-per-SPA + `cf_rescue = "all"` doesn't break it) and explicitly noting that the whole-day check lacks the Include guard. Two-check list preserved with accurate per-check guard descriptions.
- **Code:** Phase B1 code itself left in place — removing it is behavior-adjacent (edge cases breaking the invariant would silently go through Phase B2/B3 as `Exclude-C-Extraneous` instead of `Exclude-C-Identical`, losing a diagnostic distinction). Safer to keep as defensive.

### F71. "Variables created and dropped" section rewritten — FIXED
- **File/lines:** `child-algorithm-reference.md:583–585` (pre-edit).
- **Issue.** Pre-edit list named `median.spz`, `median.dopz`, `absdmedian.spz`, `absdmedian.dopz` as variables "created on the `data.sde` subset." These columns are internal scratch in `identify_temp_sde()`'s *local copy* of its input df (lines 2119–2126 of `child_clean.R`) and never reach `data.sde`. The actual Main Step 13 working variables on `data.sde` / `ewma_df` were all missing.
- **Fix.** Section rewritten to enumerate the actual working variables grouped by phase (B1 defensive, B2 one-day, B3 EWMA), with code-location pointers. Added an explicit clarification that the `identify_temp_sde()` scratch columns live on a local copy and are discarded when the helper returns. Drop mechanism (Phase B4 `setdiff(names(data.df), keep_cols_sde)` at `child_clean.R:3651–3652`) is named.

### F72. Checklist finding 1 DOP mapping wording clarified — FIXED
- **File/lines:** `child-algorithm-reference.md:589` (pre-edit).
- **Issue.** Pre-edit text read `**DOP mapping consistent:** HC → HT in all three DOP median calculations.` Three DOP medians are computed (`HT_dop_med`, `WT_dop_med`, `HC_dop_med`), not just HC-related ones; "HC → HT in all three" is grammatically confusing and suggests the asymmetric HC rule applies everywhere.
- **Fix.** Rewritten to: `All three Phase B2 DOP median branches (HT_dop_med, WT_dop_med, HC_dop_med) follow the canonical mapping — WT↔HT is symmetric, HC→HT is asymmetric (HC uses HT as its anchor; HT does not use HC).` Matches the Phase B2 narrative at line 550 which already uses the `WT↔HT, HC→HT` convention.

---

## Checklist items applied (summary)

Applied the 20-item cross-cutting checklist from `algorithm-walkthrough-procedure.md`. Items 14–19 (wrapper-only) are n/a for Step 13. Items checked and their dispositions:

1. **Sort order determinism (both).** Every `order()` / `setkey()` in Step 13 includes `internal_id` as the final tiebreaker:
    - `data.sde[order(subjid, param, agedays, internal_id)]` at `child_clean.R:3425` and `:3457`
    - `setkey(data.sde, subjid, param, agedays, internal_id)` at `:3451`
    - `setkey(ewma_df, subjid, param, agedays, internal_id)` at `:3547`
    - `order(eligible_median, eligible_dop, eligible_tiebreaker)` at `:3529` (`tiebreaker_oneday` is `+internal_id` at `agedays == 0`, `-internal_id` otherwise)
    - `order(eligible_absdewma, eligible_tiebreaker)` at `:3621` (`tiebreaker_ewma` same convention)
    - `setkey(data.df, subjid, param, agedays, internal_id)` at `:3653`. All correct.
    - **Birth tiebreaking (child only).** At `agedays == 0` keep the lowest `internal_id`; otherwise keep the highest. Verified at four sites: Phase B1 whole-day `keep_id` (lines 3429–3432), Phase B1 partial `keep_id_dup` (3439–3442), Phase B2 `tiebreaker_oneday` (3517–3518), Phase B3 `tiebreaker_ewma` (3610–3611). All use the `+id / -id` ascending-sort convention or `min / max` conditional correctly.

2. **Z-score correctness (child).** All tbc / median / EWMA work uses `tbc.sd` (not `ctbc.sd`): Phase B2 medians at `:3469`, `:3490`, `:3496`, `:3502`; Phase B3 EWMA at `:3567`; `absdewma` at `:3590`. Consistent with `identify_temp_sde()`'s use of `tbc.sd`. `ctbc.sd` is not used in Step 13 (which predates Step 15/16's corrected-z pipeline).

3. **Dead code.** Phase B1 (lines 3427–3448) is a no-op safety net in current code — see F70 above. Kept as defensive. No other dead code found. `exclude %in% c("Exclude-C-Temp-Same-Day", "Include")` filters at lines 3594, 3605, 3615, 3627, 3630 include `"Exclude-C-Temp-Same-Day"` defensively, even though by those lines all TempSDEs have been converted to Include (`:3456`); not dead because the conversion is inside `data.sde` rather than the expression itself, and the defensive inclusion costs nothing.

4. **Exclusion code names (both).** All exclusion codes assigned in Step 13 match the canonical list: `Exclude-C-Identical` (via `.child_exc(param, "Identical")` at `:3434`, `:3447`), `Exclude-C-Extraneous` (via `.child_exc(param, "Extraneous")` at `:3484`, `:3537`, `:3607`, `:3632`), `Exclude-C-Temp-Same-Day` (at `:3384` literal), `Include` (at `:3381`, `:3456`, `:3629`). All in `exclude.levels` at `:537` and adjacent.

5. **Column names (child).** All current-state: `tbc.sd`, `internal_id`, `agedays`, `v`, `v.orig` only used elsewhere. No stale `subjid`/`v` or pre-rename column references in Step 13.

6. **`.child_valid()` call correctness (child).** Single call at `:3420`: `.child_valid(data.df_sde_subset, include.temporary.extraneous = TRUE)`. Correct — Step 13 needs TempSDE rows to participate in SDE resolution. No other `.child_valid()` calls in Step 13.

7. **`by`-group correctness.** All groupings verified correct. Notably, the DOP-median calculations group by `.(subjid, agedays)` (cross-parameter, same day) — correct for DOP logic. Within-param groupings use `.(subjid, param, agedays)`. Subject-level uses `.(subjid)`. Trajectory-level uses `.(subjid, param)`. See Part 1 analysis for the full table.

8. **Sort order assumptions.** Step 13 re-sorts `data.sde` whenever it might matter (lines 3425, 3451, 3457, 3547). The final `setkey(data.df, subjid, param, agedays, internal_id)` at `:3653` restores `data.df`'s invariant sort for downstream steps.

9. **data.table reference semantics.** All `:=` calls target `data.sde`, `ewma_df`, or `data.df` as intended. `ewma_df` is subsequently merged back into `data.sde` at `:3571` by `id`; new EWMA columns on `data.sde` do not clash with existing ones. `sde_results` is merged into `data.df` at `:3640` by `id` with `all.x = TRUE`.

10. **Factor level issues.** All exclusion-code assignments go through `.child_exc()` or use the canonical literal strings. No factor-level risks.

11. **Parameter scope (child).** All three params handled: HEIGHTCM, WEIGHTKG, HEADCM in Phase B2 DOP branches (`:3488`, `:3494`, `:3500`) and Phase B2 sort branches (`:3506–3511`). Phase B3 is param-agnostic and covers all three via grouping. HEADCM's DOP is HEIGHTCM (not reversed) — matches `get_dop()`.

12. **Interaction with later steps.** Phase B4 drops all Step-13 scratch columns (`:3651–3652`) via `setdiff(names(data.df), keep_cols_sde)`. Only the `exclude` column is updated on `data.df` at `:3648`. Step 15/16 (downstream) does not depend on any Step 13 scratch column.

13. **DOP logic correctness (child).** HEIGHTCM→WEIGHTKG, WEIGHTKG→HEIGHTCM, HEADCM→HEIGHTCM in Phase B2 DOP median branches. Matches canonical mapping.

20. **CRAN output compliance (wrapper/both).** All informational output in Step 13 is wrapped in `if (!quietly) message(...)` (lines 3363–3367, 3402–3404, 3409). `warning()` calls at `:3377`, `:3389` are not gated by `quietly` — correct for bug-indicator warnings that should always surface.

---

## Discussion points resolved during walkthrough

- **Phase B1 defensiveness.** Traced every possible path that could introduce a same-day same-`v` Include/TempSDE pair into `data.sde`:
    - Initial input: Early Child Step 13 resolves before Step 5 runs.
    - `v` modification: none of Child Steps 5, 6, 7, 9, 11 or their temp-SDE reruns modify `v` (only `exclude`).
    - Row creation: no step adds rows to `data.df`.
    - Step 5 / Step 11c / Step 13 pre-phase reruns of `identify_temp_sde()`: each selects exactly one winner per (subjid, param, agedays) group, preserving one-Include-per-SPA.
    - `cf_rescue = "all"`: restored CFs sit on different days than their originators; within a single day, CFs sharing the same `v` would already be Identical-resolved by Early Step 13, leaving at most one CF per (subjid, param, agedays, v) group for Step 6 to consider.
  - Net: Phase B1 fires on no current test input. Verified tests still pass after narrative rewrite confirming defensive framing — counts unchanged.
- **Phase B3 on one-day subject-params.** After Phase B2 selects a keeper for one-day SDE groups, Phase B3 re-processes the same groups. `n_available == 1` (only the keeper remains Include) so SDE-All-Extreme (`n_available >= 2` required) doesn't fire, and the single eligible row is trivially the `keep_id_ewma`, so no state changes. Phase B3 is idempotent on Phase B2 output.
- **`ewma_fill` semantics.** For temp SDE rows whose `ewma.all` is NA (they were excluded from `ewma_df` at `:3546`), `ewma_fill := max(ewma.all[!was_temp_sde], na.rm = TRUE)` per `.(subjid, param, agedays)` at `:3575–3579` assigns the day's non-temp-SDE EWMA. In practice a day has at most one non-temp-SDE Include row per SPA (Step 5/11c invariant), so `max` is effectively "copy from the keeper's EWMA." The `max` is a deterministic fallback if that invariant is ever broken.
- **Safety-check warnings (`:3375–3391`).** Pre- and post-`identify_temp_sde()` safety checks warn if `exclude == "Include"` duplicates exist per (subjid, param, agedays). These are not gated by `quietly` — intentional, since firing indicates a bug upstream (Step 11c or the Step 13 pre-phase rerun failed to resolve all SDE groups).

---

## Session 10 deferreds

None. All Fix-now items (F66–F72) were resolved this session. Items intentionally out of scope and not opened as new deferreds:

- **Phase B1 code removal.** As noted in F70, leaving the defensive block is the safer choice. Removing it would be a behavior-adjacent change requiring regression review.
- **Child Step 15/16 line 3667 `# Add id for consistent SDE order` comment.** Same `id` vs `internal_id` looseness as F68, but in the Step 15/16 block — belongs to the Step 15/16 walkthrough, not this one.
- **`suppressWarnings(min/max)` + `is.infinite()` pattern** at `child_clean.R:3583–3589` and `:3596–3601` and `:3473–3477` — already tracked in `CLAUDE.md → Known Issues → Robustness audit → Low priority cleanup`. Not a regression, not opening a new item here.
- **Phase A narrative completeness.** The current narrative overview at lines 528–533 does not explicitly call out the reset-to-Include step at `child_clean.R:3381` between the pre-dupe safety check and the `identify_temp_sde()` re-run. Minor clarity improvement, not misleading. Left as-is.
- **Overview wording (line 500) "using EWMA-based selection."** Technically Phase B2 (one-day) uses median-based selection, not EWMA. Narrative in the per-phase description is accurate; the Overview one-liner is a simplification. Not worth a separate fix — the per-phase text is what a reader will use.

---

## Session 10 — post-fix test results

After reinstalling the package and running the full child test suite:

- test-cleangrowth.R: 63 PASS, 0 warnings
- test-child-regression.R: 48 PASS, 0 warnings
- test-child-edge-cases.R: 28 PASS, 0 warnings
- test-child-algorithms.R: 41 PASS (2 baseline codetools warnings, unchanged)
- test-child-parameters.R: 13 PASS

All child counts identical to baseline. No regressions. Adult tests not re-run — no adult files touched in this session.
