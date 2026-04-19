# GC Walk-Through — 2026-04-21: Child Step 19 (Pairs and Singles)

Dedicated detailed walkthrough of Main Child Step 19 (pairs and singles evaluation against the designated other parameter) in `R/child_clean.R` (~lines 4722–4858 pre-edit). Adult algorithm has its own `Exclude-A-Single` / `Exclude-A-Ord-Pair` / `Exclude-A-2D-*` steps in Adult Steps 10H / 11Wa / 11Wa2 / 13 that operate on raw measurements with permissiveness-driven limits rather than DOP z-score comparison, so this pass is child-scoped.

Comment / narrative / small-cleanup pass with no behavior changes. Tests unchanged across the fix set (63 / 48 / 28 / 41 / 13 pass).

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
| Step header / banner | `child_clean.R` 4722–4727 | Console banner + quietly gate |
| Pre-step sort | 4729–4731 | Re-sort `data.df` by (subjid, param, agedays, internal_id); enables `df[1]` / `df[2]` indexing inside the closure |
| valid_set construction | 4733–4741 | Compute `sp_counts` (Include-only counts per (subjid, param)); valid_set = `.child_valid()` AND sp has <= 2 Includes |
| DOP snapshot | 4743–4749 | Freeze Include rows (with `tbc.sd` / `ctbc.sd`) keyed by (subjid, param, agedays) before the per-group closure runs |
| Per-(subjid, param) closure | 4751–4833 | For each (subjid, param) with 1–2 Includes: compute `comp_diff` vs DOP; pair rule (gap-split thresholds); then single rule |
| Return to data.df | 4832 | `exclude := closure_result` assigned only for valid_set rows |

**Inventory of support functions used:** `get_dop()` at line 4761 (scalar DOP name lookup). `.child_valid()` at 4740. `.child_exc()` at three sites (two pair, one single). No `ewma()` / `identify_temp_sde()` / `calc_otl_evil_twins()` / `calc_and_recenter_z_scores()` / cache API.

---

## Fix-now items

### F103. Consolidate stale sort comments — FIXED

- **File/lines:** `R/child_clean.R:4729–4730` (pre-edit).
- **Issue.** Two-line block:
  ```
  # order just for ease later
  # Add id for consistent SDE order
  ```
  First line is vague. Second line is the F68/F73/F80/F90 stale pattern: the sort is not about "SDE order" (Step 19 is not an SDE step and same-ageday rows are resolved upstream) and uses `internal_id`, not `id`.
- **Fix.** Rewrote as a 4-line current-state comment naming the specific downstream use (`df[1]` = earlier, `df[2]` = later in the 2-row pair case) and documenting the `internal_id` tiebreaker purpose.

### F104. Stale `# Compute valid_set AFTER sort, not before` — FIXED

- **File/lines:** `R/child_clean.R:4733` (pre-edit).
- **Issue.** F58/F75 stale-changelog pattern — describes a prior bug fix order rather than the current state.
- **Fix.** Replaced with a current-state comment describing what `valid_set` does: gates to Include rows whose (subjid, param) has 1 or 2 Include measurements, counting only Includes.

### F105. DOP snapshot Problem/Solution comment — FIXED

- **File/lines:** `R/child_clean.R:4743–4747` (pre-edit).
- **Issue.** 5-line "Problem:... Solution:..." format that reads like a debugging note rather than current-state rationale.
- **Fix.** Rewrote as 6-line current-state rationale: describes what the snapshot freezes (Include rows with z-scores), when it is taken (before the closure runs), and why it is needed (data.df is modified in place; without the snapshot, the DOP lookup for the second-processed param would miss rows excluded by the first, making the result depend on by-group order).

### F106. Clarify `# save initial exclusions to keep track` — FIXED

- **File/lines:** `R/child_clean.R:4752` (pre-edit).
- **Issue.** Imprecise — does not distinguish between the two vectors `ind_all` (indices) and `exclude_all` (the closure's return vector) being initialized.
- **Fix.** Rewrote to name both vectors and their specific roles (index-based back-write vs return-vector refresh).

### F107. Stata-era sub-labels `19A / 19B-C / 19D / 19E / 19F-G / 19H` — FIXED

- **File/lines:** `R/child_clean.R:4756`, `4763`, `4779`, `4785`, `4795`, `4811` (pre-edit); six sites.
- **Issue.** Six Stata-era letter-labelled sub-step markers that appear out of sequence in the code flow (19A → 19D → 19B/C → 19E → 19F/G → 19H) and do not map cleanly onto the narrative sections (which use numbered list items). The letter labels add navigation overhead without adding information.
- **Fix.** Replaced each with a descriptive current-state comment naming what the block does (e.g., "Pair case: compute z-score differences and the time gap between the two measurements"). The narrative's numbered list remains the authoritative walkthrough of the logic; the code comments now read as current-state rationale.

### F108. Stata-era `voi` terminology — FIXED

- **File/lines:** `R/child_clean.R:4763` (pre-edit).
- **Issue.** `# 19D: calculate the voi comparison` — `voi` (value-of-interest?) is undefined Stata-era terminology. Grep confirms no other `voi` / `VOI` occurrences in `R/*.R` outside this one comment.
- **Fix.** Rewrote as a current-state description of what `comp_diff` measures (disagreement with DOP; same-day if available, median fallback otherwise) and when it is NA (subject has no DOP Includes).

### F109. Dead `else` branch via unconditional partial keyed lookup — FIXED

- **File/lines:** `R/child_clean.R:4761`, `4775–4777` (pre-edit).
- **Issue.** The lookup `dop <- dop_snapshot[.(df$subjid[1], get_dop(df$param[1]))]` used the data.table default `nomatch = NA`, which returns **one NA row** when no match is found. As a result, `nrow(dop) > 0` was always TRUE and the else branch `df[, comp_diff := rep(NA, nrow(df))]` was dead code. (Behavior was still correct because the for-loop with a 1-NA-row dop produces `median(NA) = NA` and `abs(NA - x) = NA`, so `comp_diff` ends up NA either way.)
- **Fix.** Added `nomatch = NULL` to the keyed lookup. Now `dop` is 0 rows when the subject has no DOP Includes, the `nrow(dop) > 0` check guards correctly, and the else branch is reachable. Also avoids the fruitless per-row for-loop in the no-DOP case. Verified with a small data.table reproducer that `nomatch = NULL` returns a 0-row data.table, and that `median(numeric(0))` returns NA without warning.

### F110. `# Use id tie-breaker for deterministic selection` — FIXED

- **File/lines:** `R/child_clean.R:4786` (pre-edit).
- **Issue.** F80/F84 pattern — comment says "id" but the `order()` key is `internal_id`.
- **Fix.** Rolled into the F107 rewrite of the `19E` block; the new comment names `internal_id` explicitly and documents the direction (lowest internal_id wins — that row is the one excluded).

### F111. Vague `# save the results` — FIXED

- **File/lines:** `R/child_clean.R:4808` (pre-edit).
- **Issue.** Two-word comment above `exclude_all <- df$exclude` after the pair rule.
- **Fix.** Rewrote to describe the specific purpose: refresh `exclude_all` from df so the return vector reflects the pair decision even if the single rule does not fire.

### F112. Stale changelog `# Previous code at line 4529 unconditionally overwrote exclude_all, losing 2-meas results` — FIXED

- **File/lines:** `R/child_clean.R:4825–4826` (pre-edit).
- **Issue.** F58/F75 stale-changelog pattern — describes a prior bug that was fixed, rather than the current state, and references line 4529 which no longer corresponds to anything meaningful (the code has shifted).
- **Fix.** Rewrote as a 3-line current-state comment explaining why the write is index-keyed (to preserve any pair exclusion already recorded above).

### F113. Unused `id` in `.SDcols` — FIXED

- **File/lines:** `R/child_clean.R:4833` (pre-edit).
- **Issue.** Session 7 F77 pattern — `.SDcols` includes `'id'` but the closure body never reads `df$id`. Other uses (sorting, tiebreaking, writing back) are all via `internal_id` or `index`.
- **Fix.** Removed `'id'` from `.SDcols`. Remaining columns: `'index', 'internal_id', 'agedays', 'subjid', 'param', 'tbc.sd', 'ctbc.sd', 'exclude'`.

### F114. Spacing `>=365.25` → `>= 365.25` — FIXED

- **File/lines:** `R/child_clean.R:4800` (pre-edit).
- **Issue.** Minor style: no space between `>=` and `365.25` in the long-interval pair condition. Adjacent short-interval condition (3 lines below) has proper spacing. Behavior-neutral.
- **Fix.** Added space.

### F115. Narrative Code-location cell missing support-function detail — FIXED

- **File/lines:** `child-algorithm-reference.md:1464` (pre-edit).
- **Issue.** F99 pattern — cell reads `Inline in cleanchild() in child_clean.R; uses get_dop() for designated-other-parameter lookup`, but does not mention `.child_valid()` (row eligibility mask) or `.child_exc()` (exclusion-code constructor).
- **Fix.** Rewrote to name all three support functions and their specific role in Step 19.

### F116. Narrative "lowest id as tiebreaker" — FIXED

- **File/lines:** `child-algorithm-reference.md:1479` (pre-edit).
- **Issue.** F80/F84 pattern — narrative pair-selection step 3 says "lowest `id` as tiebreaker" but the code uses `internal_id`.
- **Fix.** Rewrote to name `internal_id` and clarify direction ("the row with the lowest internal_id is the one excluded"). Also expanded the step into a longer sentence that describes what `comp_diff` and the fallback are, matching the F107 code comments.

### F117. Narrative "Variables created and dropped" lists non-existent variables — FIXED

- **File/lines:** `child-algorithm-reference.md:1505` (pre-edit).
- **Issue.** Listed `diff_tbc.sd`, `diff_ctbc.sd`, `comp_diff`, `dop_tbc`, `dop_ctbc`. The last two are not variables in the code — they are notational shorthand for "the DOP's `tbc.sd` value" in the comp_diff formula. Only `comp_diff` is an actual column added to the closure-local copy of `.SD`; `diff_tbc.sd`, `diff_ctbc.sd`, and `diff_agedays` are local scalars. Also missed `diff_agedays`.
- **Fix.** Rewrote the section to distinguish between (a) local scalars in the pair branch (`diff_tbc.sd`, `diff_ctbc.sd`, `diff_agedays`) and (b) closure-local column (`comp_diff`). Explicitly noted that no columns are persisted to `data.df` for downstream steps.

### F118. Narrative Checklist finding 6 "All 3 codes exist" — FIXED

- **File/lines:** `child-algorithm-reference.md:1514` (pre-edit).
- **Issue.** Step 19 produces exactly two exclusion codes (`Exclude-C-Pair` and `Exclude-C-Single`); "All 3 codes" is a miscount.
- **Fix.** Rewrote to name both codes explicitly and point at `exclude.levels.peds`.

### F119. Narrative Overview "After all cleaning steps" — FIXED

- **File/lines:** `child-algorithm-reference.md:1468` (pre-edit).
- **Issue.** "After all cleaning steps" implies Step 19 runs last, but Step 21 (Error Load) runs after it. More precise framing is "after the prior individual-value cleaning steps."
- **Fix.** Rewrote the opening sentence accordingly, and added "for the same subject" to the DOP comparison framing.

### F120. Narrative Scope / Operates-on phrasing — FIXED

- **File/lines:** `child-algorithm-reference.md:1460` (pre-edit).
- **Issue.** "Include rows with 1–2 remaining measurements per subject-param" — "remaining measurements" is slightly ambiguous (could be read as "remaining rows" including Excludes).
- **Fix.** Changed to "1–2 remaining Include measurements per subject-param" to match the code's Include-only sp_count logic.

### F121. Narrative Single evaluation / Checklist finding 4 potcorr framing — FIXED

- **File/lines:** `child-algorithm-reference.md:1488`, `1512` (pre-edit).
- **Issue.** Pre-edit narrative:
  - Step 5 single rule was stated as "(|tbc.sd| > 5 & no DOP data)" without mentioning `comp_diff is NA` as the mechanism.
  - Checklist finding 4 said "Both tbc.sd and ctbc.sd checked" without clarifying that singles use only `tbc.sd`. The walkthrough also surfaced the implicit potcorr protection for singles (the `comp_diff` cross-check uses the DOP's `tbc.sd`, which for a potcorr subject is also on the uncorrected scale, so the disagreement check does not double-count the uncorrected extremeness).
- **Fix.** Clarified the single rule ("comp_diff is NA" as the explicit NA-check), and expanded checklist finding 4 to note the singles-only `tbc.sd` check and the implicit potcorr guard via the `comp_diff` cross-check.

### F122. Narrative DOP snapshot / Pair evaluation expanded — FIXED

- **File/lines:** `child-algorithm-reference.md:1472`, `1476–1483` (pre-edit).
- **Issue.** Pre-edit DOP snapshot subsection did not mention that `ctbc.sd` is also snapshotted (it is, via the `.(subjid, param, agedays, tbc.sd, ctbc.sd)` j expression). Pre-edit pair evaluation numbered list did not describe `comp_diff` explicitly (only via a formula) and did not name `diff_agedays`.
- **Fix.** Added `ctbc.sd` mention to the DOP snapshot subsection. Rewrote the pair-evaluation numbered list to: (1) explicitly name `diff_tbc.sd`, `diff_ctbc.sd`, `diff_agedays`; (2) describe `comp_diff` directly rather than via ambiguous notation; (3) call out the NA-fallback case for `comp_diff`; (4) update step 5 to describe what happens after a pair exclusion.

---

## Checklist items applied (summary)

Applied the 21-item cross-cutting checklist from `algorithm-walkthrough-procedure.md`. Items 14–19 (wrapper-only) are n/a for Step 19.

1. **Sort order determinism (both).** `order()` / `setkey()` call sites all terminate with `internal_id`:
    - `data.df[order(subjid, param, agedays, internal_id),]` at `child_clean.R:4733` (pre-step sort).
    - `order(-df$comp_diff, df$internal_id)` and `order(-abs(df$tbc.sd), df$internal_id)` at `:4804`, `:4807` (pair row selection).
    - `setkey(dop_snapshot, subjid, param, agedays)` at `:4752` — no `internal_id` because the DOP lookup is by keys, not relying on row order for determinism (partial keyed lookup returns all matching rows with `nomatch = NULL`, and the `median()` / `%in%` calls inside the closure are order-independent).
    - Correct.
    - **Birth tiebreaking (child)** — not applicable. Step 19 operates on distinct (subjid, param, agedays) rows (SDEs resolved upstream), so two rows in a pair cannot share `agedays`. The `internal_id` tiebreaker on tied `comp_diff` / `abs(tbc.sd)` is a deterministic-ordering device, not a clinical tiebreaking rule. "Keep lowest at agedays == 0 / highest elsewhere" applies to SDE steps, not to Step 19.

2. **Z-score correctness (child).** Uses `tbc.sd` (primary) for all decisions and `ctbc.sd` (in pair diff only) for the corrected-check gate that protects potcorr subjects. No `sd.orig` or other non-recentered variants. Correct.

3. **Dead code.** F109 resolves the `else { df[, comp_diff := ...] }` branch (partial keyed lookup without `nomatch = NULL` always returned 1 NA row, making the branch unreachable). F113 removes unused `'id'` from `.SDcols`. No other live dead code.

4. **Exclusion code names (both).** `Exclude-C-Pair` at two sites (pre-edit `:4801`, `:4805`; post-edit `:4821`, `:4825`) and `Exclude-C-Single` at two sites (`:4824`, `:4827` pre-edit; `:4848`, `:4852` post-edit). Both codes present in `exclude.levels.peds` at `:542–543`. Correct.

5. **Column names (child).** All current-state: `subjid`, `param`, `agedays`, `internal_id`, `id`, `index`, `tbc.sd`, `ctbc.sd`, `exclude`, `comp_diff`, `sp_n`. No Stata-era or legacy names (F108 removed the one `voi` occurrence).

6. **`.child_valid()` call correctness (child).** Single call at `:4742`:
    - `.child_valid(data.df, include.temporary.extraneous = FALSE)` ANDed with `only_single_pairs`.
    - `include.temporary.extraneous = FALSE` is the default; the explicit FALSE is defensive. Correct for Step 19 — temp SDEs are resolved in Step 13 and must not participate.
    - All other `include.*` flags default to FALSE. Correct — CF non-rescues and permanent extraneous losers must not participate in a plausibility check.

7. **`by`-group correctness.** `sp_counts` computed `by = .(subjid, param)`; closure dispatched `by = .(subjid, param)`. Both groupings are what we want — Step 19 is per subject-param. Correct.

8. **Sort order assumptions.** `data.df` sorted `(subjid, param, agedays, internal_id)` at the top of the step. Within each `by = .(subjid, param)` group, `.SD` preserves the parent order, so inside the closure `df[1]` is the earliest measurement and `df[2]` is the later (when nrow == 2). The pair check `diff_tbc.sd = df[2, tbc.sd] - df[1, tbc.sd]` and `diff_agedays = df[2, agedays] - df[1, agedays]` depends on this order. Correct.

9. **data.table reference semantics.** `data.df[valid_set, exclude := closure_result]` writes to `data.df` in place; the closure's `df <- copy(.SD)` takes a fresh copy so `df[, comp_diff := ...]` and `df[max_ind, exclude := ...]` inside the closure do not mutate `data.df`. `dop_snapshot` is read-only inside the closure. Correct.

10. **Factor level issues.** Both codes present in `exclude.levels.peds`. Safe.

11. **Parameter scope (child).** All 3 params handled uniformly. No HC exception (unlike Step 17 which excludes WT). Correct. Narrative Scope row matches.

12. **Interaction with later steps.** Only `exclude` is written back to `data.df`. `comp_diff` exists only on the closure-local `df`. All diff scalars (`diff_tbc.sd`, `diff_ctbc.sd`, `diff_agedays`) are local to the pair branch and vanish on closure return. Step 21 reads `exclude`; no other columns consumed downstream from Step 19.

13. **DOP logic correctness (child).** `get_dop()` assignments: `WEIGHTKG -> HEIGHTCM`, `HEIGHTCM -> WEIGHTKG`, `HEADCM -> HEIGHTCM`. All correct per the canonical DOP rule.

14–19. **Wrapper items — n/a for Step 19.**

20. **CRAN output compliance (wrapper/both).** Informational output at `:4724–4727` uses `if (!quietly) message(...)`. No `cat()` / `print()` calls.

21. **Age gaps (new item).** Checked:
    - Pair rule gap split: `diff_agedays >= 365.25` (long) and `diff_agedays < 365.25` (short). The two halves together cover the entire real line; `diff_agedays == 365.25` falls into the long branch. No gap.
    - Pair z-score thresholds: `> 4` (long) and `> 2.5` (short). These are strict inequalities by design; exact-boundary values (e.g., `diff_tbc.sd == 4.0` in the long branch) do not exclude. Deliberate behavioral choice, not an inadvertent gap.
    - Single thresholds: `abs(tbc.sd) > 3 & comp_diff > 5` (with DOP) and `abs(tbc.sd) > 5 & comp_diff is NA` (no DOP). Strict inequalities, same as pairs. Deliberate.
    - No `>= A & < B`, `>= B & < C`, ..., `>= Z` binning clauses in Step 19 (no age-binned or gap-binned reference tables).
    - No conditional `:=` inside a loop (no per-group iteration loop — Step 19 is one-pass per (subjid, param)).

---

## Deferreds (noted, not opened as new todos unless flagged)

- None. All identified items were fix-now. The pair and single z-score thresholds are deliberately strict (`>`), not a gap to close.

---

## Discussion points resolved during walkthrough

- **Partial keyed lookup returns 1-NA-row by default.** Verified with a minimal reproducer that `dop_snapshot[.("C", "HT")]` on an empty match returns a single-row data.table with NA columns (key fields NA too for non-key lookups). `nrow(dop) > 0` was therefore always TRUE, making the pre-edit `else { df[, comp_diff := rep(NA, nrow(df))] }` branch dead. The pre-edit behavior was still correct (the for-loop produces NA `comp_diff` when dop has no real data) but dead code is dead code. F109 switches the lookup to `nomatch = NULL` so the if/else guard is meaningful and the no-DOP case skips the per-row loop.
- **Non-potcorr `ctbc.sd` equals `tbc.sd`.** Confirmed by reading the wrapper recentering pipeline: `data.all[, ctbc.sd := sd.corr - sd.median]` at `:1095` and `data.all[, sd.corr := sd.orig]` at `:999` / `:1012` for non-potcorr subjects. So for non-potcorr, `ctbc.sd = sd.orig - sd.median = tbc.sd`, and the pair rule's `|diff_ctbc.sd|` check is equivalent to the `|diff_tbc.sd|` check. The `| is.na(diff_ctbc.sd)` disjunct is a defensive guard for any case where `ctbc.sd` happens to be NA, but for normal data it is not the mechanism that exempts non-potcorr subjects — they are already exempt via the equivalent check.
- **Single rule has no explicit `ctbc.sd` check but has an implicit potcorr guard.** Analyzed during the walkthrough: for a potcorr subject with one extreme-looking uncorrected `tbc.sd` (say, 6) but a normal corrected `ctbc.sd` (say, 0), the DOP is also likely potcorr and also has extreme-looking uncorrected values on the same trajectory. `comp_diff` is computed against the DOP's `tbc.sd` (same uncorrected scale), so for a well-behaved potcorr subject `comp_diff` is small and the `abs(tbc.sd) > 3 & comp_diff > 5` branch does not fire. The `abs(tbc.sd) > 5 & comp_diff is NA` branch only fires when the subject has no DOP data at all. Documented as a rationale point in the narrative checklist finding 4.
- **`19A–H` sub-labels are out of sequence and Stata-era.** Confirmed by reading the code flow: `19A` is "is it single or pair"; `19D` is "compute comp_diff"; `19B/C` is "pair differences"; `19E` is "pick worst"; `19F/G` is "apply pair rule"; `19H` is "re-evaluate remaining as single". Letter order does not match execution order, the labels do not correspond to narrative sections, and `voi` (in 19D) is undefined terminology. F107 removes all six labels and replaces them with descriptive current-state comments.
- **`dop_tbc`, `dop_ctbc` are not actual variables.** Verified by grep — they appear only in the narrative's formulae `|dop_tbc - tbc|` and `|median(dop_tbc) - tbc|` as notational shorthand, not as column names in the code. Only `comp_diff` is an actual column, and only on the closure-local copy of `.SD`. F117 removes them from the narrative's "Variables created and dropped" list.

---

## Post-fix test results

After reinstalling the package and running the full child test suite:

- test-cleangrowth.R: 63 PASS, 0 warnings
- test-child-regression.R: 48 PASS, 0 warnings
- test-child-edge-cases.R: 28 PASS, 0 warnings
- test-child-algorithms.R: 41 PASS (2 baseline codetools warnings, unchanged)
- test-child-parameters.R: 13 PASS

All child counts identical to baseline. No regressions, as expected since the fixes are comment-only, narrative-only, unused-column removal, or a dead-branch guard (behavior-equivalent). Adult tests not re-run — no adult files touched in this session.
