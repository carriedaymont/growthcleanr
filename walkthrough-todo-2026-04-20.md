# GC Walk-Through — 2026-04-20: Child Step 17 (Height/HC Velocity)

Dedicated detailed walkthrough of Main Child Step 17 (raw-diff / velocity checks for HT and HC) in `R/child_clean.R` (~lines 4094–4706 pre-edit). Adult algorithm has no velocity / raw-diff step and uses raw-measurement-based limits instead, so this pass is child-scoped. Both the pre-filter (vectorized pass over all valid rows) and the per-(subjid, param) while-loop closure were reviewed.

Comment / narrative / small-cleanup with six small behavior-narrowing fixes (all passed tests; counts unchanged from baseline):

- F94: Tanner NA-out filter now matches the merge key — uses `tanner.months < 30` (midpoint-based, same quantity used in the Tanner merge) instead of this row's own agedays (`agedays/30.4375 < 30`). The previous check excluded pairs from Tanner based on starting age even when the midpoint was old enough for the reference. Same fix in the pre-filter and the inner loop.
- F95: The paired `< 365.25` / `> 365.25` formulas (mindiff, maxdiff — both Tanner and Tanner-to-WHO-fallback) left `d_agedays == 365.25` unassigned; the two formulas are continuous at gap = 1 year so switching `>` to `>=` is behavior-equivalent for integer agedays and covers the non-integer boundary. Same fix in pre-filter and inner loop.
- F97: The 2-row branch's strict `abs(tbc.sd) > shift(abs(tbc.sd), ...)` pair checks used `>`, which meant ties produced no exclusion at all; changed to `>=` so both rows flag, and the existing `order(-absval, internal_id)` candidate selection picks the lowest-`internal_id` tied row.
- F101: The pre-filter's `whoagegrp.ht > 24 | whoagegrp.ht.lead > 24` NA-out was dead (both branches impossible because the earlier cap at 24 months prevents either value from exceeding 24). The intended behavior — per Carrie — was to NA-out `whoagegrp.ht` when the next measurement crosses the 24-month WHO boundary, so a pair extending past 24 months does not extrapolate the WHO reference. Pre-filter now uses `(agedays + d_agedays) / 30.4375 > 24` directly; inner loop applies the same per-iteration check to `who_mindiff_ht` / `who_maxdiff_ht` after the fcase (whoagegrp.ht was pre-merged before the while loop on static per-row data, so boundary-crossing depends on the current iteration's d_agedays).
- F102 (collapse): `d_agedays < 20 → 1` and `20 ≤ d_agedays < 46 → 1` collapsed into a single `d_agedays < 46 → 1` clause in both pre-filter and inner loop; and the `153 ≤ d_agedays < 199 → 6` + `d_agedays >= 200 → 6` split collapsed into `d_agedays >= 153 → 6`, closing the former gap at `d_agedays == 199` exactly.
- **whoinc.age.ht NA init per iteration** (defense-in-depth): `df[, whoinc.age.ht := NA_integer_]` added at the top of the HT block inside the while loop so any future boundary change that reintroduces a gap cannot leave a stale value from the prior iteration. The pre-filter already had an explicit NA init.

Child test counts unchanged across all behavior changes (63 / 48 / 28 / 41 / 13 pass).

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
| Velocity table loading | `child_clean.R` 4107–4149 | Tanner pre-loaded in `cleangrowth()` and passed in; WHO HT / HC velocity tables loaded inside `cleanchild()` (see CLAUDE.md Known Issues for `gc_preload_refs()` extension candidate) |
| Pre-sort + valid_set | 4151–4164 | Re-sort `data.df` by (subjid, param, agedays, internal_id); compute HT/HC / non-single / `.child_valid()` mask |
| Step 17 pre-filter (HT) | 4166–4284 | Vectorized: `d_agedays` / `diff_prev` / `diff_next` per (subjid, param); Tanner merge + midpoint filter + max.ht.vel floors + mindiff/maxdiff; WHO merge + fcase + gap scaling + WHO-vs-Tanner choice + default fallback + birth adjustment |
| Step 17 pre-filter (HC) | 4286–4362 | Same skeleton as HT but WHO-only, 2-month smallest interval, ±1.5 tolerance, ±0.5 birth |
| Pre-filter wrap-up | 4364–4382 | `mindiff_prior` / `maxdiff_prior` via shift; `has_violation` flag; `sp17_to_process` set of groups that actually need the while loop |
| Per-group while loop | 4384–4694 | Inside-closure work: HT branch re-runs Tanner + WHO per iteration (with static pre-merge of whoagegrp.ht and who.ht.vel); HC branch re-runs WHO per iteration; 17O computes EWMA; 17P/R flag violating pair members (3+ rows) or compare abs(tbc.sd) directly (2 rows); exclude worst candidate; re-sort + repeat |
| Cleanup | 4695–4697 | Drop `sp_key` column |

**Inventory of support functions used:** `ewma()` called directly at line 4593 (not the cache API — EWMA is fully recomputed each iteration). `.child_valid()` at 4162 for the pre-step row mask. `.child_exc()` inside the closure for `Exclude-C-Abs-Diff`. No `identify_temp_sde()`, `calc_otl_evil_twins()`, `calc_and_recenter_z_scores()`, or `get_dop()`. (Narrative Code-location cell updated — F99.)

---

## Fix-now items

### F86. Stale "Bug fix:" comment on Tanner table alias — FIXED

- **File/lines:** `R/child_clean.R:4102–4105` (pre-edit).
- **Issue.** Comment `# Bug fix: was redundantly re-reading from disk each batch.` is stale changelog (F58/F75 pattern). The rest of the paragraph is current-state.
- **Fix.** Rewrote to current state without the "Bug fix:" framing.

### F87. Stale "# Chris updated this to >= from ==" comments — FIXED

- **File/lines:** `R/child_clean.R:4465`, `4470` (pre-edit); two sites sandwiching the inner-loop `whoinc.age.ht` edge-interval assignments.
- **Issue.** Stale changelog — described that a prior refactor changed `==` to `>=`, not what the code currently does.
- **Fix.** Removed both sites while doing the F102 collapse (the edge-interval assignments themselves were folded into the main clauses). A sanity grep in `R/*.R` confirmed no other "Chris updated / added / removed / fixed" remnants anywhere.

### F88. Stale "# Fix Min-diff tie-breaking" comment — FIXED

- **File/lines:** `R/child_clean.R:4663` (pre-edit).
- **Issue.** Stale changelog pattern (F58/F75); the current code is the Min-diff tie-breaking logic, not a fix to it.
- **Fix.** Rewritten as current-state rationale: which dewma side (before vs after) is paired with which exclusion code (1/3 vs 2/4) and why.

### F89. Stale line-number reference `lines 4030-4031` — FIXED

- **File/lines:** `R/child_clean.R:4499` (pre-edit).
- **Issue.** Comment read `For gap < 9 months: use WHO (already transformed at lines 4030-4031)`, but the transformation actually lives at the current lines 4495–4496 (the `who_mindiff_ht * .5 - 3` / `who_maxdiff_ht * 2 + 3` block just above). Line numbers drift with every edit, so naming them in a comment is fragile.
- **Fix.** Replaced the line-number pointer with a section-name pointer ("already transformed in 17G above"); 17G is an immediately-above labelled block that won't move around the way line numbers do.

### F90. `# Add id for consistent SDE order` stale comments at four sort sites — FIXED

- **File/lines:** `R/child_clean.R:4152`, `4421`, `4519`, `4570` (pre-edit).
- **Issue.** F68/F73/F80 pattern — comment says "id" but the `order()` key is `internal_id`. In Step 17 the tiebreaker is effectively unused (SDEs are resolved upstream so same-ageday rows don't enter Step 17) but the sort includes it for free determinism.
- **Fix.** All four sites rewritten to current state: name the downstream dependency (pre-loop shifts + merges, inner-loop pair checks) and note that `internal_id` is a deterministic tiebreaker against a rare same-ageday case that doesn't normally reach Step 17.

### F91. 3-line tiebreaker comment consolidation in Step 17 closure — FIXED

- **File/lines:** `R/child_clean.R:4683–4685` (pre-edit).
- **Issue.** Three comment lines expressing two thoughts (the "id tie-breaker" wording of F80 and the `order()` vs `which.max` rationale of F78):
  ```
  # choose the highest absval for exclusion
  # Use id tie-breaker for deterministic selection
  # which.max returns first tie, which depends on data order
  ```
- **Fix.** Collapsed to a single two-line claim + rationale matching the Session 11 F78 convention; uses `internal_id` wording and names "lowest internal_id wins" explicitly for ties.

### F92. Dead `id_all <- copy(df$id)` — FIXED

- **File/lines:** `R/child_clean.R:4390` (pre-edit).
- **Issue.** The line stashes a copy of the input `id` vector but nothing inside the closure ever reads `id_all`. (`ind_all` — the parallel copy for `index` — is used to write back exclusions by original row position, but `id_all` has no such consumer.)
- **Fix.** Removed; surrounding comment rewritten to name `ind_all`'s specific role.

### F93. HC narrative interval table did not match code — FIXED

- **File/lines:** `child-algorithm-reference.md:1389–1401` (pre-edit).
- **Issue.** Narrative listed `gap < 20 → 1`, `20 ≤ gap < 46 → 1`, and `gap ≥ 200 → 6` for HC, framed as "same gap-based mapping as HT but using HC-specific reference tables." None of those are in the HC code: the smallest HC interval is 2-month (46–75 days), and there is no `≥ 200` fallback — HC uses `153 ≤ gap < 200 → 6` only. HC pairs with `gap < 46` or `gap ≥ 200` (or HC itself above ~24 months) fall through to the HC default `mindiff = -1.5`, `maxdiff = NA`.
- **Fix.** Rewrote the HC subsection to match code: four intervals only (2/3/4/6 months; no 1-month), no `≥ 200` fallback, explicit note that pairs outside the reference band use only the decrease-side default (`mindiff = -1.5`) and have no upper bound — consistent with HC's narrower physiologic range. Also corrected the claim that HC has "the same gap-based mapping as HT."

### F94. Tanner NA-out filter uses this row's age instead of midpoint — FIXED

- **File/lines:** `R/child_clean.R:4197`, `4441` (pre-edit).
- **Issue.** `tanner.months` is computed from the midpoint age of the adjacent-measurement pair, and the merge uses `tanner.months` as the key. But the NA-out filter used `(agedays / 30.4375) < 30` — this row's own agedays — which is not the merge key. Pairs like (2y, 5y) with midpoint 3.5y merge fine on tanner.months but were NA'd out because the starting age < 30 months. Per Carrie this has been wrong for a long time.
- **Fix.** Changed both sites to `tanner.months < 30` (with `!is.na()` guard), matching the merge key. Side effect: pairs whose midpoint is ≥ 2.5y but whose starting age is younger now use Tanner — expected. Narrative updated to match.

### F95. `d_agedays == 365.25` gap at mindiff / maxdiff formula boundary — FIXED

- **File/lines:** `R/child_clean.R:4221`, `4225` (pre-filter, pre-edit); `4452`, `4456` (inner loop, pre-edit).
- **Issue.** The paired `< 365.25` / `> 365.25` clauses for mindiff and maxdiff (both Tanner and Tanner-WHO-fallback variants) left `d_agedays == 365.25` unassigned. Integer agedays can't hit 365.25 exactly, but the pattern is a silent footgun for any future non-integer inputs or for an analogous future pattern elsewhere.
- **Fix.** Changed the `> 365.25` branches to `>= 365.25` in both pre-filter and inner loop. The two formulas (exponents 2 / flat for mindiff, 0.33 / 1.5 for maxdiff) are continuous at gap = 1 year, so the boundary assignment matches either limit — behavior-equivalent for integer agedays and robust for non-integer. Narrative updated to match.

### F96. EWMA caching question (2-row case + overall Step 17) — DEFERRED

- **File/lines:** `R/child_clean.R:4593` (inner-loop EWMA call).
- **Issue.** Step 17 calls `ewma()` directly (the full-recompute path) rather than the incremental `ewma_cache_init()` / `ewma_cache_update()` used by Steps 11 / 15 / 16. And in the 2-row branch, `ewma()` is still invoked even though the result is not consumed (the 2-row tiebreaker uses `abs(tbc.sd)` directly).
- **Disposition.** Carrie asked to defer. Step 17's pre-filter already prunes non-violating groups, and the remaining groups are typically small, so the potential gain from caching or from gating the 2-row EWMA call is modest. Noted as a performance item for a later pass.

### F97. 2-row branch strict `>` → `>=` for tie-handling — FIXED

- **File/lines:** `R/child_clean.R:4637–4652` (pre-edit); eight comparison lines.
- **Issue.** All eight `abs(tbc.sd) > shift(abs(tbc.sd), ...)` checks used strict `>`. A 2-row violating pair with tied `abs(tbc.sd)` flagged NEITHER row, so the violation escaped Step 17 entirely — neither abs-diff exclusion was applied despite the pair being outside the mindiff/maxdiff envelope.
- **Fix.** Changed all eight comparisons to `>=`. In tied cases both rows are now flagged; the existing `candidates <- df[val_excl != "Include"]` + `order(-absval, internal_id)` selection then picks the lowest-internal_id tied row. Non-tied cases are unaffected because only one row satisfies `>=` when abs values differ. Narrative updated to match.

### F98. Narrative "Configurable parameters in scope: None" — FIXED

- **File/lines:** `child-algorithm-reference.md:1436–1438` (pre-edit).
- **Issue.** The narrative said Step 17 has no configurable parameters. Incorrect: `ewma_window` (default 15) is passed through to the `ewma()` call at line 4593 and controls how many neighbors on each side contribute to the 3+-measurement tiebreaker.
- **Fix.** Added a one-row "Configurable parameters" table listing `ewma_window` with its default and where it is consumed. The remaining sentence about hardcoded velocity thresholds and reference-table-derived values kept below the table.

### F99. Narrative Code-location cell missing support-function detail — FIXED

- **File/lines:** `child-algorithm-reference.md:1339` (pre-edit).
- **Issue.** Cell read `Inline in cleanchild() in child_clean.R; uses Tanner height velocity and WHO velocity reference tables loaded from inst/extdata/`. Same gap that F82 closed for Step 15/16 — specific support functions were not named.
- **Fix.** Rewrote to name `ewma()` (direct call, no cache), `.child_valid()`, `.child_exc()`, and to distinguish Tanner (pre-loaded once in `cleangrowth()` and passed in) from WHO HT / HC velocity tables (loaded inside `cleanchild()`; see Known Issues for the `gc_preload_refs()` extension candidate).

### F100. Tautology comment after the dead-branch collapse — FIXED

- **File/lines:** `R/child_clean.R:4694` (pre-edit).
- **Issue.** Comment `# count_exclude > 0 outer guard makes count_exclude >= 1 always true here.` was a vestige from the Session 11 pre-walkthrough collapse of a dead `if (count_exclude >= 1)` inner branch. The remaining comment described the tautology rather than the intent.
- **Fix.** Rewrote to describe why iteration continues (after excluding a candidate, new neighbors and d_agedays can produce fresh violations).

### F101. Dead WHO boundary-crossing check — FIXED (behavior change)

- **File/lines:** `R/child_clean.R:4231–4234` (pre-edit, pre-filter); inner loop had no equivalent check.
- **Issue.** Pre-filter had a two-branch OR intended to NA-out `whoagegrp.ht`:
  ```
  pf[ht_idx & ((!is.na(whoagegrp.ht) & whoagegrp.ht > 24) |
                (!is.na(whoagegrp.ht.lead) & whoagegrp.ht.lead > 24)),
     whoagegrp.ht := NA_integer_]
  ```
  Both branches are dead: (1) `whoagegrp.ht > 24` is impossible because the preceding assignment only writes when `agedays / 30.4375 ≤ 24`, capping the value at 24; (2) `whoagegrp.ht.lead` is a `shift()` of that same capped column, also never > 24. The intent — per Carrie — was "exclude if the second height is > 24 months", i.e., NA-out WHO for any pair whose next measurement crosses the 24-month WHO/CDC boundary.
- **Fix.**
    - Pre-filter: replaced the dead OR with a direct check on the next measurement's age, `(agedays + d_agedays) / 30.4375 > 24`. Dropped the now-unused `whoagegrp.ht.lead` column.
    - Inner loop: `whoagegrp.ht` is computed once pre-loop (agedays is static), but the boundary check depends on `d_agedays`, which changes each iteration. Solution: after the per-iteration `fcase` that populates `who_mindiff_ht` / `who_maxdiff_ht`, NA those two columns out for boundary-crossing pairs — achieves the same effect as NA-ing `whoagegrp.ht` without restructuring the pre-loop merge.
- **Behavior change.** For HT pairs with the current row ≤ 24 months and the next row > 24 months:
  - Previous: WHO velocity is used (either directly for `gap < 9 months`, or as Tanner fallback when `min.ht.vel` is NA for `gap ≥ 9 months`).
  - After fix: WHO is not used. If Tanner is available (now more often, after the F94 midpoint fix), Tanner applies. Otherwise the `mindiff = -3`, `maxdiff = NA` default fallback applies.
- **Test impact.** Child test counts unchanged (63 / 48 / 28 / 41 / 13).

### F102 (part 1 — collapse `< 20` + `20–45` into `< 46`) — FIXED

- **File/lines:** `R/child_clean.R:4229–4230` (pre-filter, pre-edit); `4458`, `4467` (inner loop).
- **Issue.** Two clauses both assign `whoinc.age.ht := 1L` but split the range as `< 20` and `20 ≤ gap < 46`. Per Carrie, `< 20` should not be a separate step; the correct behavior is a single `< 46` → 1 clause.
- **Fix.** Collapsed both the pre-filter and inner-loop assignments to a single `d_agedays < 46 → 1`. The inner-loop edge-interval block (which previously held the `< 20` line plus the stale "Chris updated" comments) disappears entirely. Narrative updated to match.

### F102 (part 2 — close `d_agedays == 199` gap) — FIXED

- **File/lines:** `R/child_clean.R:4234–4235` (pre-filter, pre-edit); `4462`, `4469` (inner loop).
- **Issue.** The `153 ≤ d_agedays < 199 → 6` and `d_agedays >= 200 → 6` split left `d_agedays == 199` unassigned, producing a single-integer gap where WHO did not apply. Combined with the inner loop's absence of an NA initializer, a row whose `d_agedays` transitioned to 199 between iterations kept its prior iteration's `whoinc.age.ht` value (and a knock-on compounding `-1.5` birth-adjustment bug via stale mindiff).
- **Fix.** Collapsed both sites to a single `d_agedays >= 153 → 6` clause. 199 is now covered. Narrative updated.

### F-bycatch 1. `iter_count` dead variable — FIXED

- **File/lines:** `R/child_clean.R:4415`, `4418` (pre-edit).
- **Issue.** `iter_count <- 0` and `iter_count <- iter_count + 1` — initialized and incremented but never read. Dead code left over from an earlier debug or tracing pattern.
- **Fix.** Removed both lines. Rewrote the surrounding sort comment at the same time (see F90 site 2).

### F-bycatch 2. Procedural addition: walkthrough checklist item 21 "Look for age gaps" — FIXED

- **File/lines:** `algorithm-walkthrough-procedure.md` (new item 21 in the 20-item cross-cutting checklist, plus a new row in the "Child-Specific Patterns to Watch For" table).
- **Rationale.** The F101 and F102-part-2 fixes both stemmed from the same latent category — threshold-binning clauses that leave an integer boundary value unassigned, potentially combined with a conditional `:=` pattern inside an iterative loop that lets stale prior-iteration values propagate. Added as a formal procedure item so future walkthroughs catch this category up front. Carrie is planning a focused retrospective pass over already-reviewed steps with this lens.

### whoinc.age.ht NA init per iteration (defense-in-depth) — FIXED

- **File/lines:** `R/child_clean.R:4477` (new line, inner loop HT block).
- **Context.** Even with F102-part-2 closing the `d_agedays == 199` gap, the underlying pattern (conditional `:=` assignments inside a loop with no NA reset) is a silent footgun — any future boundary change that reintroduces a gap would resurface the staleness. Per Carrie: "There should be no gaps, but I'd rather be safe."
- **Fix.** Added `df[, whoinc.age.ht := NA_integer_]` at the top of the inner-loop HT block, so any row whose `d_agedays` fails to match any subsequent clause gets NA rather than stale data.
- The inner-loop HC block already has an explicit NA init (`df[, whoinc.age.hc := NA_integer_]`) at the top of its assignment block, so no parallel change needed there.

---

## Checklist items applied (summary)

Applied the 20-item cross-cutting checklist from `algorithm-walkthrough-procedure.md` (plus the new item 21 introduced this session). Items 14–19 (wrapper-only) are n/a for Step 17.

1. **Sort order determinism (both).** `order()` / `setkey()` call sites:
    - `data.df[order(subjid, param, agedays, internal_id),]` at `child_clean.R:4153` (pre-step sort)
    - Pre-filter internal `pf <- pf[order(subjid, param, agedays, internal_id)]` at `:4173` (defensive; `pf` is a subset of the just-sorted `data.df`)
    - `setkey(df, sex, tanner.months)` at `:4442` (inner loop — merge key for Tanner)
    - `df <- df[order(agedays, internal_id),]` at three inner-loop sites (`:4427`, `:4527`, `:4578`)
    - `order(-candidates$absval, candidates$internal_id)` at `:4689` (worst-candidate selection, 2-row and 3+-row paths)
    - All sorts terminate with `internal_id`. Correct.
    - **Birth tiebreaking (child)** — not applicable in the SDE sense. Step 17 operates on distinct (subjid, param, agedays) rows (SDEs resolved upstream), so there are no ties on agedays to resolve as "keep lowest at agedays == 0 / highest elsewhere". The `internal_id` tiebreaker on `absval` is a deterministic-ordering device, not a clinical tiebreaking rule.

2. **Z-score correctness (child).** The only z-score used is `tbc.sd` (EWMA input at `:4593`, 2-row absval at `:4680`, 3+-row `dewma.*` at `:4597–4600`). No `ctbc.sd` / `sd.orig` / recentered perturbation variants in Step 17. Correct.

3. **Dead code.** Post-edit: no live dead code. Resolved: dead `id_all` (F92); dead `iter_count` (F-bycatch 1); dead `whoagegrp.ht > 24 | whoagegrp.ht.lead > 24` OR (F101, promoted from dead-code to wrong-code-fix). `dewma.all` is still computed at `:4605` but not consumed in the Step 17 closure (only `dewma.before` / `dewma.after` feed the tiebreaker); left in place for parallelism with the shared EWMA-output idiom — cost is one subtraction per row, clarity gain of the `dewma.*` trio outweighs the micro-efficiency.

4. **Exclusion code names (both).** Only one code assigned in Step 17: `"Exclude-C-Abs-Diff"`, via `.child_exc(param, "Abs-Diff")` at 12 sites in the inner loop (3+-row block and 2-row block). `Exclude-C-Abs-Diff` is in `exclude.levels.peds` at `:540`. Correct.

5. **Column names (child).** All current-state: `tbc.sd`, `ctbc.sd`, `internal_id`, `id`, `index`, `subjid`, `param`, `agedays`, `sex`, `v`, `exclude`, `d_agedays`, `diff_prev`, `diff_next`, `tanner.months`, `min.ht.vel`, `max.ht.vel`, `whoagegrp.ht`, `whoinc.age.ht`, `whoinc.1.ht`–`whoinc.6.ht`, `max.whoinc.1.ht`–`max.whoinc.6.ht`, `who_mindiff_ht`, `who_maxdiff_ht`, `mindiff`, `maxdiff`, `mindiff_prior`, `maxdiff_prior`, `pair`, `bef.g.aftm1`, `aft.g.aftm1`, `val_excl`, `val_excl_code`, `absval`, `ewma.all`, `ewma.before`, `ewma.after`, `dewma.all`, `dewma.before`, `dewma.after`, `exp_vals`, `whoagegrp.hc`, `whoinc.age.hc`, `who_mindiff_hc`, `who_maxdiff_hc`. HC variants use `.hc` suffix. No Stata-era or legacy names.

6. **`.child_valid()` call correctness (child).** Single call at `:4162`:
    - `.child_valid(data.df, include.temporary.extraneous = FALSE)` ANDed with `not_single` (sp-level row count > 1) and `param != "WEIGHTKG"` (weight excluded). Correct: Step 17 operates only on HT / HC rows that are Include and belong to subject-params with at least two measurements; singles are handled in Step 19, and weight can legitimately decrease so it is exempt entirely.
    - All other `include.*` flags default to FALSE — appropriate for velocity (extraneous-SDE losers and CF non-rescues are permanent exclusions and must not contribute raw-diff comparisons).

7. **`by`-group correctness.** Pre-filter uses `by = .(subjid, param)` for `d_agedays` / `diff_prev` / `diff_next` / `mindiff_prior` / `maxdiff_prior` shifts; `by = .(subjid)` for `tanner.months` midpoint calculation (the midpoint is across the HT trajectory, not per-param since Step 17 never runs on WT). Inner-loop closure is dispatched `by = .(subjid, param)`. All groupings verified. Correct.

8. **Sort order assumptions.** Pre-filter's per-group shifts and the inner loop's pairwise pair-check depend on (agedays, internal_id) order within each (subjid, param) group. Guaranteed by the pre-step sort at `:4153` and the inner-loop re-sorts at `:4427`, `:4527`, `:4578`. Correct.

9. **data.table reference semantics.** `:=` targets are `data.df` (outer step), `pf` (pre-filter local copy), or the closure-local `df` (from `copy(.SD)`). No unintended parent-df mutations. The inner-loop `df <- merge(df, who.ht.vel, ...)` / `df <- tanner.ht.vel.rev[df]` return new data.tables; the subsequent `:=` writes to those new objects (still inside the closure). Correct.

10. **Factor level issues.** `Exclude-C-Abs-Diff` is in `exclude.levels.peds` at `:540`. No other codes assigned. Safe.

11. **Parameter scope (child).** HT and HC only; weight exempt via the valid_set filter. HC has its own velocity tables, tighter tolerances, no Tanner. Narrative checklist item 4 matches.

12. **Interaction with later steps.** Columns created by Step 17 pre-filter all live on `pf` and are discarded by `rm(pf)` at `:4381`. Closure-local columns all discarded on closure return. `sp_key` (shared with Steps 15/16) is dropped at `:4697`. Only `exclude` is written back to `data.df`. Step 19 reads `exclude` and `tbc.sd` / `ctbc.sd`; neither is modified by Step 17 beyond exclusion flagging.

13. **DOP logic correctness (child).** Step 17 does not use DOP. Weight's-DOP-is-height etc. matters for Steps 6 / 11 / 13 / 15 / 19 but raw-diff velocity operates within a single (subjid, param) group.

21. **Age gaps (new item).** F94, F95, F101, F102 (parts 1 and 2) all identified and fixed under this lens. Post-edit, the whoinc.age.ht bin coverage is complete over all integer d_agedays values; the mindiff / maxdiff formula boundary at 365.25 is closed via `>=`; the Tanner NA filter matches the merge key; the WHO boundary-crossing NA-out applies the intended "next > 24 months" rule. The whoinc.age.ht NA init closes the staleness side-channel in case a future edit reintroduces a gap.

20. **CRAN output compliance (wrapper/both).** Informational output in Step 17 is `if (!quietly) message(...)` at two sites (`:4097–4100` entry banner and `:4378–4380` pre-filter diagnostic). No `cat()` / `print()` calls.

---

## Deferreds (noted, not opened as new todos unless flagged)

- **F96 (Step 17 EWMA caching).** Step 17 uses `ewma()` directly rather than `ewma_cache_*`. Pre-filter already prunes non-violating groups and per-group iteration counts are typically small, so the potential gain is modest. Additionally, the 2-row branch still invokes `ewma()` even though the result is not used (the 2-row tiebreaker keys on `abs(tbc.sd)` directly). Both are behavior-neutral efficiency items. Carrie explicitly deferred.
- **Known Issue (already in CLAUDE.md).** WHO velocity tables loaded inside `cleanchild()` rather than via `gc_preload_refs()`. Tiny files (< 1 kB each gzipped); correctness unaffected. Extension candidate called out in CLAUDE.md Known Issues (wrapper).
- **`dewma.all` computed but not consumed inside Step 17.** Kept for parallelism with the shared EWMA-output idiom (the `ewma.all` / `ewma.before` / `ewma.after` trio maps directly to the `dewma.all` / `dewma.before` / `dewma.after` trio). Removing only `dewma.all` would break symmetry; the cost of one subtraction per row is negligible.
- **Pre-filter redundant `pf <- pf[order(...)]` at `:4173`.** `pf` is a subset of the just-sorted `data.df`, so the re-sort is defensive but behavior-equivalent. Same pattern as Session 11's deferred redundant `order()` on `data.df`. Low-value cleanup.
- **`not_single` via `table(paste0(...))`.** Could be replaced by a data.table `by` aggregation for style consistency with the rest of Step 17. Low-value cleanup.

---

## Discussion points resolved during walkthrough

- **F94 directionality.** Carrie confirmed the narrative description (midpoint-based filter) was the intended behavior and the code (this-row-agedays-based filter) has been wrong "for a long time." The fix is `tanner.months < 30`, not any variant of the midpoint expression directly.
- **F101 interpretation.** The two-branch OR was originally presented as dead code; Carrie clarified the intended rule — "exclude if the second height is > 24 months." That made the code wrong (not dead), and the fix uses the next row's actual ageyears directly (`(agedays + d_agedays) / 30.4375 > 24`).
- **F102 <20 separation.** Carrie confirmed there is no separate clinical behavior for `d_agedays < 20`; the two clauses should both map to the 1-month WHO interval. Collapsed to a single `< 46 → 1` clause.
- **whoinc.age.ht staleness scope.** Clarified that the conditional `:=` in the inner loop only overwrites matching rows, so a row whose `d_agedays` fails to match any clause between iterations retains its prior iteration value. With the F102-part-2 gap-close plus the defense-in-depth NA init, staleness is fully neutralized under the current boundaries.
- **`dewma.all` consumption.** Verified that `dewma.all` is computed at line 4605 in Step 17 but never read downstream within the Step 17 closure; only `dewma.before` and `dewma.after` feed the 3+-row tiebreaker. Left in place for EWMA-trio parallelism.

---

## Post-fix test results

After reinstalling the package and running the full child test suite:

- test-cleangrowth.R: 63 PASS, 0 warnings
- test-child-regression.R: 48 PASS, 0 warnings
- test-child-edge-cases.R: 28 PASS, 0 warnings
- test-child-algorithms.R: 41 PASS (2 baseline codetools warnings, unchanged)
- test-child-parameters.R: 13 PASS

All child counts identical to baseline. No regressions, even though six of the fixes (F94, F95, F97, F101, F102 collapse + 199 gap close, plus the whoinc.age.ht NA init) were behavior changes in principle. The test data does not include pairs that exercise the affected boundary conditions — consistent with these being rare edge cases. Adult tests not re-run — no adult files touched in this session.
