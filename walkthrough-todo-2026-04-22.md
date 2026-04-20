# GC Walk-Through — 2026-04-22: Child Step 21 (Error Load)

Dedicated detailed walkthrough of Main Child Step 21 (error-load escalation to all remaining Includes when the per-(subjid, param) error ratio exceeds a threshold) in `R/child_clean.R` (~lines 4860–4901 pre-edit). Adult algorithm has its own `eval_error_load()` support function in `adult_support.R` (Adult Step 14) with permissiveness-driven thresholds (0.41 / 0.29), param-specific error-code lists, and a `denom < 3` skip rule rather than a `mincount` floor, so this pass is child-scoped.

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
| Step header / banner | `child_clean.R` 4860–4866 | Console banner + quietly gate |
| valid_set + non_error_codes | 4868–4879 | `valid_set` (unused; all TRUE) + list of 5 codes removed from both numerator and denominator |
| By-group ratio computation | 4881–4892 | `data.df[, c("err_ratio", "n_errors") := {...}, by = .(subjid, param)]` computes the ratio per (subjid, param) |
| Exclusion assignment | 4894–4898 | `exclude := "Exclude-C-Too-Many-Errors"` for Include rows in groups where `err_ratio > threshold` AND `n_errors >= mincount` |
| Cleanup | 4900–4901 | Drop the two working columns |

**Inventory of support functions used:** `.child_exc()` at one site (exclusion code constructor). No `.child_valid()` / `identify_temp_sde()` / `get_dop()` / `ewma()` / any cache API. (Not using `.child_valid()` is correct — the step needs to see ALL rows, including BIV / Evil-Twins / etc. which must be counted as errors; the counting logic filters via the `non_error_codes` list instead.)

---

## Fix-now items

### F123. Dead `valid_set <- rep(TRUE, nrow(data.df))` — FIXED

- **File/lines:** `R/child_clean.R:4871`, plus `valid_set &` / `valid_set,` uses at `:4881` and `:4897` (pre-edit).
- **Issue.** Unlike Steps 7 / 9 / 11 / 17 / 19 where `valid_set` is computed via `.child_valid()` with meaningful filters, Step 21's `valid_set` is unconditionally TRUE for every row, so the three uses are no-ops. Step 21 is order-independent and correctly counts all rows (error / non-error classification happens inside the counting expression via the `non_error_codes` list).
- **Fix.** Removed the `valid_set <- rep(TRUE, nrow(data.df))` assignment and both uses. The two data.table calls now read `data.df[, ...]` (by-group computation) and `data.df[err_ratio > ... & ... & exclude == "Include", ...]` (exclusion write).

### F124. Stale changelog comment above `non_error_codes` — FIXED

- **File/lines:** `R/child_clean.R:4873–4874` (pre-edit).
- **Issue.** F58/F75 pattern — `# Non-error codes that should be excluded from both numerator AND denominator` followed by `# CF rescue codes removed — rescued CFs are now "Include" (stored in cf_rescued column)`. The second line describes what was removed during a prior fix rather than the current state.
- **Fix.** Rewrote as a current-state comment naming each code's role:
    - `Exclude-C-Identical` / `Exclude-C-Extraneous`: SDE housekeeping resolved in Child Step 13.
    - `Exclude-C-CF`: detected carried-forward that was not rescued (data-entry artifact). Rescued CFs are `Include`, tracked via `cf_rescued`, and count as Includes in the denominator.
    - `Exclude-Missing`: no measurement; cannot be a measurement error.
    - `Exclude-Not-Cleaned`: HC beyond cleaning age (> 3 × 365.25 days); out of scope.

### F125. Redundant short-form denominator comment — FIXED

- **File/lines:** `R/child_clean.R:4887` (pre-edit).
- **Issue.** `# Denominator is errors + includes (excludes SDE/CF/Missing)` — both redundant (the outer comment at `:4868–4870` already explains the ratio structure) and inaccurate (missing `Not-Cleaned`).
- **Fix.** Removed. The `# Count errors (not Include and not in non_error_codes)` label above the first counting line kept (shortened slightly), and the two assignment lines below it (`n_includes <- ...`, `denom <- ...`, `err_ratio <- ...`) now read as straight code without redundant commentary.

### F126. roxygen `error.load.mincount` description imprecise — FIXED

- **File/lines:** `R/child_clean.R:211–212` (pre-edit).
- **Issue.** "minimum count of exclusions on parameter before considering excluding all measurements. Defaults to 2." — "exclusions" is too broad (SDE/CF/Missing/Not-Cleaned exclusions do not count), and "on parameter" is ambiguous (the scope is per (subject, param) group, not per-param-across-all-subjects).
- **Fix.** Rewrote to name the exact exclusion codes that do NOT count ("real errors" = everything other than Identical / Extraneous / CF / Missing / Not-Cleaned) and to state the per-(subject, param) scope explicitly.

### F127. roxygen `error.load.threshold` description slightly misleading — FIXED

- **File/lines:** `R/child_clean.R:213–214` (pre-edit).
- **Issue.** "threshold of percentage of excluded measurement count to included measurement count that must be exceeded before excluding all measurements of either parameter" — "percentage" is imprecise (it's a ratio in [0, 1]); "excluded measurement count to included measurement count" reads like `errors : includes` but the actual formula is `errors / (errors + includes)`; "either parameter" suggests both params are evaluated together but the trigger is per-(subject, param).
- **Fix.** Rewrote to state the formula explicitly (`n_errors / (n_errors + n_includes)`), name the per-(subject, param) scope, note the SDE/CF/Missing/Not-Cleaned exclusion from both numerator and denominator, cross-reference `error.load.mincount`, and name the resulting exclusion code.

### F128. Narrative Code-location cell missing `.child_exc()` — FIXED

- **File/lines:** `child-algorithm-reference.md:1527` (pre-edit).
- **Issue.** F99/F115 pattern — cell reads `Inline in cleanchild() in child_clean.R` with no mention of the support function used.
- **Fix.** Added `.child_exc()` (builds `Exclude-C-Too-Many-Errors`).

### F129. Narrative Overview omits `Not-Cleaned` and narrows scope of exclusion — FIXED

- **File/lines:** `child-algorithm-reference.md:1531` (pre-edit).
- **Issue.** "excluding SDEs, CFs, and Missing from the denominator" — omits `Not-Cleaned`, and "from the denominator" understates the rule (codes are excluded from BOTH numerator and denominator). Also "all remaining Include values are also excluded" is missing the per-subject-param scope.
- **Fix.** Expanded to name all four categories (SDEs, CFs, Missing, Not-Cleaned), replaced "from the denominator" with "from both numerator and denominator", and added "for that subject-param" to the exclusion scope sentence.

### F130. Narrative Logic step 5 implicit per-group scope — FIXED

- **File/lines:** `child-algorithm-reference.md:1539` (pre-edit).
- **Issue.** "all Include rows → `Exclude-C-Too-Many-Errors`" — "all Include rows" is ambiguous; the trigger scope is per-(subjid, param) group (since `err_ratio` and `n_errors` are group-level).
- **Fix.** Changed to "all Include rows within that (subjid, param) group → `Exclude-C-Too-Many-Errors`".

### F131. Narrative "Variables created and dropped" incorrectly lists `n_includes` as a column — FIXED

- **File/lines:** `child-algorithm-reference.md:1559` (pre-edit).
- **Issue.** "Working columns (`n_errors`, `n_includes`, `err_ratio`) are computed via `data.table` by-group operations on `data.df` and not persisted." Only `err_ratio` and `n_errors` are assigned as columns on `data.df` (via `c("err_ratio", "n_errors") := {...}`); `n_includes` and `denom` are local scratch variables inside the `{}` block and never become columns.
- **Fix.** Rewrote to name the two actual columns (`err_ratio`, `n_errors`), note they are constant within a (subjid, param) group and dropped at the end of the step, and describe `n_includes` / `denom` as local scratch inside the expression.

### F132. Narrative Checklist finding 3 omits Missing / Not-Cleaned — FIXED

- **File/lines:** `child-algorithm-reference.md:1565` (pre-edit).
- **Issue.** "Denominator excludes SDEs/CFs: These are data-structure artifacts, not cleaning errors, so they are excluded from both numerator and denominator." Internally inconsistent (heading says "denominator", body says "both"); and the list of excluded categories omits Missing and Not-Cleaned.
- **Fix.** Rewrote as "Numerator and denominator exclude SDEs/CFs/Missing/Not-Cleaned: These are data-structure artifacts or out-of-scope rows, not cleaning errors, so they are excluded from both sides of the ratio (via `non_error_codes`)."

---

## Checklist items applied (summary)

Applied the 21-item cross-cutting checklist from `algorithm-walkthrough-procedure.md`. Items 14–19 (wrapper-only) are n/a for Step 21.

1. **Sort order determinism (both).** Step 21 does not sort; counts are order-independent within each (subjid, param) group. No `order()` / `setkey()` calls. No `internal_id` tiebreaking needed because the step produces one scalar per (subjid, param) group, not per-row decisions.
    - **Birth tiebreaking (child)** — not applicable. Step 21 has no row-level tiebreaking.

2. **Z-score correctness (child).** n/a. Step 21 uses only the `exclude` column and counting by (subjid, param).

3. **Dead code.** F123 (unused `valid_set`) fixed. No other dead code.

4. **Exclusion code names (both).** Only one code assigned: `Exclude-C-Too-Many-Errors`, via `.child_exc(param, "Too-Many-Errors")` at one site. Code present in `exclude.levels.peds` at `:544`. Correct.

5. **Column names (child).** All current-state: `err_ratio`, `n_errors`, `exclude`, `subjid`, `param`. No Stata-era or legacy names.

6. **`.child_valid()` call correctness (child).** Step 21 does NOT call `.child_valid()`, and this is correct by design. The step needs to see ALL rows (BIV / Evil-Twins / EWMA / velocity / pair / single exclusions all must contribute to `n_errors`); filtering first via `.child_valid()` would drop the real errors and miscount the ratio. Instead, the counting logic filters non-errors via the `non_error_codes` list explicitly.

7. **`by`-group correctness.** `by = c("subjid", "param")` on the ratio computation. Correct — Step 21 operates per (subject, param), and within a group `err_ratio` / `n_errors` are scalars shared by all rows.

8. **Sort order assumptions.** None. The by-group counts do not depend on row order.

9. **data.table reference semantics.** `:=` assigns to `data.df` in place (adds `err_ratio` and `n_errors` columns); the final `:= NULL` drops them. The exclusion assignment writes to `data.df$exclude` for Include rows in triggering groups. No accidental copies or joined-table mutations.

10. **Factor level issues.** `Exclude-C-Too-Many-Errors` is in `exclude.levels.peds` at `:544`. Safe.

11. **Parameter scope (child).** All 3 params (WT / HT / HC), grouped by `(subjid, param)`. No HC exception.

12. **Interaction with later steps.** Step 22 (output assembly) reads only `exclude` / z-score columns. `err_ratio` / `n_errors` are dropped in Step 21 cleanup, so no leakage.

13. **DOP logic correctness (child).** n/a. Step 21 does not use DOP.

14–19. **Wrapper items — n/a for Step 21.**

20. **CRAN output compliance (wrapper/both).** Informational output at `:4862–4866` uses `if (!quietly) message(...)`. No `cat()` / `print()` calls.

21. **Age gaps (new item).** n/a. No age-binning clauses, gap-to-interval mappings, or iterative `:=` assignments in Step 21.

---

## Edge cases verified during walkthrough

- **denom = 0 (all rows in `non_error_codes`).** `err_ratio` explicitly set to 0 → `err_ratio > 0.5` is FALSE → no trigger. Safe.
- **n_errors >= mincount but n_includes == 0.** Trigger fires at the by-group level, but the row-level filter `exclude == "Include"` matches 0 rows → no-op.
- **n_errors = 1 and n_includes = 1 (50% ratio with 2 measurements).** Blocked by the default `mincount = 2` guard. Behavior matches the narrative rationale.
- **HC rows beyond 3y.** Tagged `Exclude-Not-Cleaned` in preprocessing, excluded from both numerator and denominator via `non_error_codes`. Correct.
- **`Exclude-C-Temp-Same-Day` at Step 21.** Not in `non_error_codes`; would count as an error. Per Child Step 13 invariants, no temp SDE should survive to Step 21 — if one did it would indicate an upstream bug, and the current behavior (counting it) is correct: a diagnostic signal rather than a silent pass-through.

---

## Discussion points resolved during walkthrough

- **`valid_set` origin.** Git history shows `valid_set <- rep(TRUE, nrow(data.df))` has been in place since the initial v3.0.0 commit (`27bcc0f`). It does not parallel any earlier meaningful filter (it was always all-TRUE). Clean removal.
- **Why `.child_valid()` is not used.** Confirmed by reading the counting expression: `n_errors <- sum(!exclude %in% c("Include", non_error_codes))` needs rows with exclusions like `Exclude-C-BIV`, `Exclude-C-Evil-Twins`, `Exclude-C-Traj-Extreme`, `Exclude-C-Traj`, `Exclude-C-Abs-Diff`, `Exclude-C-Pair`, `Exclude-C-Single` to be counted as errors. `.child_valid()` filters those out (as "not currently Include"), so calling it would collapse `n_errors` to 0 (or very small) and defeat the step entirely.
- **Child vs adult design differences.** Child uses a single non-param-specific error list and a `mincount` error-count floor; adult uses param-specific error-code enumeration (height vs weight), different thresholds per permissiveness, and a `denom >= 3` sample-size floor. Both designs are valid for their respective contexts — documented in the fix-now writeups but not normalized.
- **Parameter table row audit (checklist item 6).** Confirmed: `error.load.mincount` (default 2) and `error.load.threshold` (default 0.5) defaults match the signature (`:343–344`), the `cleanchild()` internal signature (`:2526–2527`), the CLAUDE.md Configurable Parameters table, and the narrative Configurable parameters table at `:449–450`. All four sources consistent.

---

## Deferreds

- None. All identified items were fix-now.

---

## Post-fix test results

After reinstalling the package and running the full child test suite:

- test-cleangrowth.R: 63 PASS, 0 warnings
- test-child-regression.R: 48 PASS, 0 warnings
- test-child-edge-cases.R: 28 PASS, 0 warnings
- test-child-algorithms.R: 41 PASS (2 baseline codetools warnings, unchanged)
- test-child-parameters.R: 13 PASS

All child counts identical to baseline. No regressions, as expected — all fixes are dead-variable removal, comment rewrites, roxygen rewording, or narrative cleanup (no behavior changes). Adult tests not re-run — no adult files touched in this session.

---

# R-vs-R comparison — Session 4 — 2026-04-22

Scope: Child Step 7 (BIV) + Child Step 9 (Evil Twins). Per `R-child-AprJanComparison-procedure.md`. Pure R-vs-R diff of reference `Infants_Main.R` against current `child_clean.R`.

---

## Pre-session baseline tests

Baseline from today's Step 21 walkthrough session (same process instance):

- test-cleangrowth.R: 63 PASS
- test-child-regression.R: 48 PASS
- test-child-edge-cases.R: 28 PASS
- test-child-algorithms.R: 41 PASS (2 baseline codetools warnings)
- test-child-parameters.R: 13 PASS

---

## Scope

| Sub-area | Reference lines | Current lines | What |
|---|---|---|---|
| Step 7 absolute BIV | `Infants_Main.R:3049–3111` | `child_clean.R:2991–3039` | Absolute weight / HT / HC limits |
| Step 7 standardized BIV | `Infants_Main.R:3113–3158` | `child_clean.R:3041–3090` | Z-cutoffs + temp-SDE rerun + ageyears drop |
| `calc_oob_evil_twins` / `calc_otl_evil_twins` | `Infants_Main.R:2387–2429` | `child_clean.R:2254–2288` | OTL detection helper |
| Step 9 main | `Infants_Main.R:3160–3255` | `child_clean.R:3092–3204` | Evil Twins loop, exclusion, temp-SDE rerun |

---

## Known intentional changes encountered (logged briefly, not analyzed)

- **BIV exclusion codes unified:** `Exclude-Absolute-BIV` + `Exclude-Standardized-BIV` → `Exclude-C-BIV`; guard changed from `exclude != abs_biv` to `!grepl("^Exclude-C-BIV$", exclude)` as a side effect. Known code rename.
- **`biv.z.*` parametrization:** Eight per-cell parameters replace hardcoded `-25`/`-15`/`22`/`8`/`-15`/`15` z-thresholds. Listed in procedure's "Known intentional changes."
- **`oob` → `otl`, `calc_oob_evil_twins` → `calc_otl_evil_twins`:** Known OOB → OTL rename.
- **`id` → `internal_id` tiebreaker:** Known intentional throughout.
- **`valid()` → `.child_valid()`:** Known intentional rename.
- **`cat()` → `message()`:** Known intentional.
- **`Exclude-Evil-Twins` → `Exclude-C-Evil-Twins`:** Known code rename.
- **Rounding-tolerance removal in `calc_otl_evil_twins`:** Reference had `janitor::round_half_up(janitor::round_half_up(..., 3), 2)` before the `> 5` comparison; current omits this. Per procedure: do NOT flag.

---

## Findings

### AJ10 — Step 7: absolute BIV weight age boundary — Intentional (other) — closed

**Reference (`Infants_Main.R:3073–3075`):**
```r
data.df[valid_set & param == "WEIGHTKG" & v < 0.2 & agedays == 0,  exclude := exc_nam]
data.df[valid_set & param == "WEIGHTKG" & v < 1 & agedays != 0,    exclude := exc_nam]
```

**Current (`child_clean.R:3007–3010`):**
```r
data.df[valid_set & param == "WEIGHTKG" & v < 0.2 & agedays <= 365, exclude := "Exclude-C-BIV"]
data.df[valid_set & param == "WEIGHTKG" & v < 1   & agedays > 365,  exclude := "Exclude-C-BIV"]
```

Reference applies the strict 0.2 kg floor only at birth (agedays == 0), with the 1 kg floor for all other days. Current extends the 0.2 kg floor to the entire first year (agedays ≤ 365), with the 1 kg floor only after the first year.

Intentional bug fix documented in `__Pipeline/CLAUDE.md` (2026-03-16 GC archetype testing): "BIV threshold bug found and fixed: Absolute BIV weight limit was excluding legitimate preterm weights (0.7–1.0 kg). Changed to <0.2 kg for first year, <1 kg after that." Not listed in the procedure's "Known intentional changes" but confirmed intentional via CLAUDE.md history. **Pitfall: Boundary changes.**

No code change needed. Closed.

---

### AJ11 — Step 9: `any(start_df$otl, na.rm = TRUE)` vs `any(start_df$oob)` — Bug fix — closed

**Reference (`Infants_Main.R:3187`):**
```r
if (any(start_df$oob)) {
```

**Current (`child_clean.R:3123`):**
```r
if (any(start_df$otl, na.rm = TRUE)) {
```

Current adds `na.rm = TRUE`. Without it, if `otl` contains NAs (which occurs when `tbc.sd` or `ctbc.sd` is NA for any row in the valid set) and no TRUE values are present, `any(...)` returns NA; `if (NA)` errors with "missing value where TRUE/FALSE needed." In normal operation, Include rows at Step 9 should have valid z-scores (BIV excluded implausible z-scores; Missing/Not-Cleaned rows excluded from valid_set), so this edge case is rare, but the fix is correct. **Pitfall: NA / empty-set handling.**

Current already correct; no change needed. Closed.

---

## Items NOT flagged (audit trail)

**Step 7:**
- `v == 0` explicit tagging (`Infants_Main.R:3154`: `data.df[v == 0, exclude := exc_nam]`) absent from current: Dead code removal. The wrapper converts `measurement == 0` to NaN before dispatch; NaN rows get `Exclude-Missing` at init and are excluded from `valid_set`. No v==0 row can reach Step 7's valid_set. Removal is correct.
- `ageyears` column explicitly dropped after Step 7 in current (`child_clean.R:3090`: `data.df[, ageyears := NULL]`): Cleanup. Grep confirms `ageyears` in `data.df` is only created and used in Step 7 in both reference and current; later occurrences of `ageyears` are all local variables inside function scopes, not the `data.df` column. Reference leaves the column on `data.df` without dropping it (minor leak); current tidies it. No logic impact.
- Stale `# Removed nnte filter (nnte calculation removed)` comment in reference Step 7 (`Infants_Main.R:3066`): Already removed in current as part of cleanup.

**Step 9:**
- `paste0/table` → `.N by-group` for `not_single_pairs` (`Infants_Main.R:3177–3178` vs `child_clean.R:3114–3115`): Equivalent behavior. The paste-based table approach has the same theoretical underscore-concat collision risk as AJ4 in Session 1 — impossible in practice given fixed `WEIGHTKG/HEIGHTCM/HEADCM` param vocabulary (same audit-trail note as Sessions 1 and 3). Not re-flagged.
- Single data.table closure (reference) → per-group for-loop (current): Equivalent results. Evil Twins processing is fully independent per (subjid, param) group: OTL detection uses `same_sp_next` / `same_sp_prev` guards so cross-group comparisons never occur; the median tiebreaker is per-group; one exclusion per iteration within each group converges to the same final exclusion set regardless of whether groups are processed in parallel (reference's `by = .(subjid, param)` within the closure) or sequentially (current's for-loop). Result sets are identical.
- `sum_oob` / `any_oob` intermediate variables (reference) → direct `any(df$otl, na.rm = TRUE)` loop condition (current): Equivalent termination condition; `sum_oob >= 1` ↔ `any(otl == TRUE)` per group.
- `data.df[, sp_count_9 := NULL]` cleanup at end of Step 9 (`child_clean.R:3197`): New cleanup in current (reference leaves the column). Not a logic issue.
- Stale `# Removed nnte filter` and the old `# we only running carried forwards on valid values, non NNTE values` comment in reference Step 9 (`Infants_Main.R:3175–3176`): Already removed in current.

---

## Open questions

None.

---

## Session 4 status

2 findings (AJ10 — Intentional (other); AJ11 — Bug fix). Both closed with no code change needed. Baseline unchanged at 63 / 48 / 28 / 41 / 13; no tests re-run. Next session candidate: **Session 5 — Child Step 11 (EWMA1, complex closure and EWMA caching).**

---

# R-vs-R comparison — Session 8 — 2026-04-22

Scope: Child Step 17 (Height/HC Velocity). Per `R-child-AprJanComparison-procedure.md`. Pure R-vs-R diff of reference `Infants_Main.R` (~lines 4210–4715) against current `child_clean.R` (~lines 4115–4733).

---

## Pre-session baseline tests

Baseline from today's prior sessions (same process instance):

- test-cleangrowth.R: 63 PASS
- test-child-regression.R: 48 PASS
- test-child-edge-cases.R: 28 PASS
- test-child-algorithms.R: 41 PASS (2 baseline codetools warnings)
- test-child-parameters.R: 13 PASS

---

## Scope

| Sub-area | Reference lines | Current lines | What |
|---|---|---|---|
| Velocity table loading | `Infants_Main.R:4210–4283` | `child_clean.R:4115–4162` | tanner + WHO HT + WHO HC tables |
| Sort + valid_set + pre-filter | `Infants_Main.R:4285–4300` | `child_clean.R:4164–4405` | Pre-filter (perf opt in current) + valid_set |
| Step 17 closure | `Infants_Main.R:4301–4715` | `child_clean.R:4408–4733` | Per-(subjid, param) while loop |

---

## Known intentional changes encountered (logged briefly, not analyzed)

- **Pre-filter structure (performance):** Current builds a vectorized `pf` pre-filter over all valid rows, identifies groups with at least one raw-diff violation (`sp17_to_process`), and runs the closure only on those groups. Reference runs the closure on all valid groups. Equivalent results; groups without violations trivially exit the while loop in the reference.
- **WHO merge moved pre-loop:** Current pre-merges WHO HT velocity columns before the while loop (inside the closure, before the `testing <- TRUE / while` block), and uses `fcase` inside the loop to select the per-row value. Reference merges inside the while loop per iteration, causing `.x/.y` column suffix duplication on iteration 2+ when merge encounters columns already added in iteration 1. The pre-loop approach is both correct and more efficient. Reference merge at `Infants_Main.R:4387–4388`; current pre-merge at `child_clean.R:4420–4422`.
- **Tanner table loading pre-batch:** Reference loads `tanner.ht.vel.rev` fresh inside each per-batch `cleanbatch_infants()` call. Current loads it once in `cleangrowth()` and passes it as a parameter; aliased to `tanner.ht.vel.rev` at `child_clean.R:4118` for naming compatibility.
- **`dplyr::lead/lag` → `shift`:** Cosmetic refactor throughout. Reference used `dplyr::lead(agedays)` etc.; current uses `shift(agedays, n = 1L, type = "lead")` etc. Equivalent.
- **`.SDcols` narrowed:** Reference uses `.SDcols = colnames(data.df)` (all columns); current uses a specific 9-column list. Equivalent; current is more efficient.
- **`id` → `internal_id` as final tiebreaker:** Sort and candidate selection use `internal_id` in current, `id` in reference. Known rename.
- **`valid()` → `.child_valid()`:** Known rename.
- **`cat()` → `message()`:** Known rename.
- **Exclusion codes:** Reference uses `"Exclude-Min-diff"` / `"Exclude-Max-diff"`; current uses `.child_exc(param, "Abs-Diff")` → `"Exclude-C-Abs-Diff"`. Known code rename.
- **Rounding removal:** Reference applies `janitor::round_half_up(janitor::round_half_up(..., 3), 2)` to mindiff/maxdiff, diff_prev/diff_next, and dewma comparisons. Per procedure: do NOT flag.
- **17B Tanner months threshold:** F94 (walkthrough Session 12) — reference uses `(agedays/30.4375) < 30` (this row's age); current uses `tanner.months < 30` (midpoint-based). Known bug fix.
- **17C-b boundary `>` → `>=` for d_agedays == 365.25:** F95 — reference `d_agedays > 365.25` left exactly 365.25 unhandled; current `d_agedays >= 365.25`. Known bug fix.
- **17D WHO endpoint check (17D F101):** Reference used `whoagegrp.ht > 24 | dplyr::lead(whoagegrp.ht) > 24` (wrong: compares capped values, not actual endpoint); current uses `(agedays + d_agedays) / 30.4375 > 24`. Known bug fix (F101, Session 12).
- **17E whoinc.age.ht interval consolidation (F102):** Reference: `d_agedays >= 20 & < 46 → 1` and separate `d_agedays < 20 → 1`, with gap at d_agedays == 199; current: `d_agedays < 46 → 1`, `d_agedays >= 153 → 6` (gap closed). Known bug fix (F102, Session 12).
- **2-row branch `>` → `>=` for tied abs(tbc.sd) (F97):** Reference strict `>` meant tied-absval pairs flagged neither row, letting violations escape Step 17 entirely; current `>=` flags both, then `internal_id` picks the lowest. Known bug fix (F97, Session 12).
- **Dead `if (count_exclude >= 1)` else-branch collapsed:** Reference has `if (count_exclude > 0) { ... if (count_exclude >= 1) { testing <- TRUE; df <- df[index != idx, ] } else { testing <- FALSE } }` — inner if is always TRUE inside the outer if, making the else dead. Current collapses to `if (count_exclude > 0) { ... testing <- TRUE; df <- df[index != idx, ] } else { testing <- FALSE }`. Known pre-walkthrough cleanup.
- **Dead `id_all <- copy(df$id)` removed:** Reference line 4306 creates `id_all` but never reads it; already removed in current (confirmed in walkthrough Session 12). Not re-flagged.
- **Dead `iter_count` variable removed:** Reference lines 4311/4314 initialize and increment `iter_count` but never read it; already removed in current. Not re-flagged.
- **Reference debug block:** Lines 4644–4666 print to console when `unique(df$subjid) == "76234"`, plus commented-out code at lines 4683–4712. Debugging artifacts; already cleaned up in current.
- **`ewma()` now passes `window = ewma_window`:** AJ12 (Session 5) — known intentional parameter addition.

---

## Findings

### AJ14 — Step 17: HC velocity section in reference is effectively a no-op — Bug fix — closed

**Reference (`Infants_Main.R:4452–4523`, HEADCM branch inside closure):**

```r
# 17D
df[, whoagegrp.ht := round(agedays/30.4375)]           # sets .ht, not .hc
df[whoagegrp.ht > 24 | dplyr::lead(whoagegrp.ht) > 24,
   whoagegrp.ht := NA]

# 17J
df[d_agedays >= 46 & d_agedays < 76,  whoinc.age.hc := 2]
df[d_agedays >= 76 & d_agedays < 107, whoinc.age.hc := 3]
df[d_agedays >= 107 & d_agedays < 153,whoinc.age.hc := 4]
df[d_agedays >= 153 & d_agedays < 199,whoinc.age.hc := 6]  # gap at 199

# 17K
m_who_hc_vel <- merge(df, who.hc.vel, by = c("sex", "whoagegrp.ht"), ...)
for (i in unique(df$whoinc.age.hc[!is.na(df$whoinc.age.hc) &
                                  !is.na(df$whoagegrp.hc)])){  # whoagegrp.hc NEVER SET
  sub_m_who_hc_vel <- m_who_hc_vel[whoinc.age.ht == i,]  # wrong: uses .ht
  cn <- paste0("whoinc.", i, ".ht")                        # wrong: .ht suffix
  df[whoinc.age.hc == i,
     who_mindiff_hc := as.numeric(sub_m_who_ht_vel[, get(cn)])]  # wrong: ht df
  ...
}
```

The for-loop condition `!is.na(df$whoagegrp.hc)` is always FALSE because `whoagegrp.hc` is never set in the HC branch — the code sets `whoagegrp.ht` instead. The preallocate block at lines 4486–4492 runs (because `length(unique(...)) < 1`), setting `who_mindiff_hc = NA` and `who_maxdiff_hc = NA`. Then `df[, mindiff := who_mindiff_hc]` assigns all NA, and `df[is.na(mindiff), mindiff := -1.5]` applies the fixed fallback for every row. Result: HC velocity checking in the reference applies only `mindiff = -1.5` (or -2.0 at birth) with **no upper bound check**. Additionally, if the loop had somehow run, it would reference `sub_m_who_ht_vel` (HT data frame, undefined in HEADCM context) and `.ht`-suffix column names — causing an error or silently using HT data.

**Current (`child_clean.R:4419–4434`, pre-loop HC merge; `child_clean.R:4557–4609`, HEADCM in-loop):**

Pre-loop: copies `who.hc.vel`, renames all `whoinc.N.ht` / `max.whoinc.N.ht` / `whoagegrp.ht` columns to `.hc` suffix, creates `whoagegrp.hc` in `df` (NA for rows > 24 months), merges renamed table into `df`. Inside loop: computes `whoinc.age.hc` per-iteration (corrects the `d_agedays >= 153 & < 199` gap to `< 200`); uses `fcase` to extract `who_mindiff_hc` / `who_maxdiff_hc` from the pre-merged `.hc` columns; applies the 17L gap-scaling and 17M tolerance (±1.5 cm); falls back to `mindiff = -1.5` when no WHO data available.

Result: HC velocity checking in the current code correctly applies WHO HC velocity limits for pairs within the 0–24 month range, with a proper upper bound check and the standard gap-scaled tolerance. **Bug fix — closed (current already correct, no change needed).** `Infants_Main.R:4452–4523` vs `child_clean.R:4557–4609`. **Pitfalls: NA / empty-set handling** (loop condition always FALSE due to unset variable) + **Parameter scope** (HC-specific variable `whoagegrp.hc` never initialized).

---

## Items NOT flagged (audit trail)

- **HC 24-month endpoint boundary check not present for HC:** Unlike HT (where F101 added an explicit per-iteration NA-out of `who_mindiff_ht`/`who_maxdiff_ht` when `(agedays + d_agedays) / 30.4375 > 24`), the current HC section has no equivalent per-iteration endpoint check. In practice, all WHO HC velocity reference values correspond to `d_agedays < 200` (< 6.6 months), which is less than the 9-month threshold; a HC pair starting just before 24 months with a 6-month gap would end just before 30 months, still within the HC cleaning age (< 3 × 365.25 days). Edge case is rare and the practical impact is small. Not flagged as a regression since the reference HC section was entirely non-functional; noted here for completeness.
- **HC `whoinc.age.hc` boundary d_agedays == 199 gap:** Reference `d_agedays >= 153 & < 199` leaves d_agedays == 199 unhandled; current uses `< 200`. Analogous to F102 for HT — part of the broader HC fix (AJ14); no independent flag since reference HC was non-functional anyway.

---

## Open questions

None.

---

## Session 8 status

1 finding (AJ14 — Bug fix). Closed with no code change needed. Baseline unchanged at 63 / 48 / 28 / 41 / 13; no tests re-run since no code change. Next session candidate: **Session 9 — Child Step 19 (Pairs/Singles) + Step 21 (Error Load) + Step 22 (Output).**

---

# R-vs-R comparison — Session 9 — 2026-04-22

Pure R-vs-R diff of reference `Infants_Main.R` Steps 19, 21, and 22 (~lines 4725–4960) against current `child_clean.R` Steps 19, 21, and 22 (~lines 4735–4966). Per `R-child-AprJanComparison-procedure.md`.

---

## Pre-session baseline tests

63 / 48 / 28 / 41 / 13 — unchanged from Session 8.

---

## Scope

| Sub-area | Reference lines | Current lines | What |
|---|---|---|---|
| Step 19 (Pairs/Singles) | `Infants_Main.R:4725–4857` | `child_clean.R:4735–4871` | 1–2 remaining measurements evaluation |
| Step 21 (Error Load) | `Infants_Main.R:4859–4908` | `child_clean.R:4873–4918` | Error-load ratio and escalation |
| Step 22 (Output) | `Infants_Main.R:4918–4960` | `child_clean.R:4928–4966` | final_tbc assembly and return columns |

---

## Known intentional changes encountered (logged briefly)

- **`id` → `internal_id`** tiebreaker in Step 19 sort and `order()` calls — known rename.
- **`valid()` → `.child_valid()`** in Step 19 valid_set — known rename.
- **Rounding removed** (all `janitor::round_half_up()` calls in Steps 19 and 21) — per procedure, do NOT flag.
- **Exclusion codes renamed** — `"Exclude-2-meas->1-year"` / `"Exclude-2-meas-<1-year"` → `.child_exc(param, "Pair")`; `"Exclude-1-meas"` → `.child_exc(param, "Single")`; `"Exclude-Error-load"` → `.child_exc(param, "Too-Many-Errors")` — known rename.
- **`valid_set <- rep(TRUE, nrow(data.df))` removed** from Step 21 — F123 from Session 14 walkthrough; was unconditionally TRUE since v3.0.0 initial commit.
- **Error load threshold `.4` → `error.load.threshold`** — known intentional bug fix, documented in procedure's "Known intentional changes" list.
- **`cf_rescued` added to return_cols** in Step 22 — part of new CF rescue scheme (Step 6 intentional change).
- **`ewma1_it1.*` loop removed** from Step 22 return_cols — stale debug columns from `debug` parameter removed in Session 6 (F34).
- **`cf_detail` columns added** to Step 22 return — part of new CF rescue scheme.
- **Step 21 `non_error_codes` updated** to new exclusion code names (`Exclude-SDE-Identical` → `Exclude-C-Identical`, SDE sub-codes → `Exclude-C-Extraneous`, `Exclude-Carried-Forward` → `Exclude-C-CF`, `"Missing"` → `Exclude-Missing`) — known exclusion code rename.
- **Step 21 old CF rescue codes removed from `non_error_codes`** (`Exclude-1-CF-deltaZ-<0.05` etc.) — rescued CFs are now `Include` rows, so they no longer need to be listed as non-error exclusion codes; they count as Includes in the denominator.

---

## Findings

### AJ15 — Step 21: `"Not cleaned"` omitted from reference `non_error_codes`

**Category:** Bug fix

**Status:** closed (current already correct, no code change needed)

**Reference** (`Infants_Main.R:4873–4884`):
```r
non_error_codes <- c("Exclude-SDE-Identical",
                     "Exclude-SDE-All-Exclude",
                     "Exclude-SDE-All-Extreme",
                     "Exclude-SDE-EWMA",
                     "Exclude-SDE-One-Day",
                     "Exclude-Carried-Forward",
                     "Exclude-1-CF-deltaZ-<0.05",
                     "Exclude-1-CF-deltaZ-<0.1-wholehalfimp",
                     "Exclude-Teen-2-plus-CF-deltaZ-<0.05",
                     "Exclude-Teen-2-plus-CF-deltaZ-<0.1-wholehalfimp",
                     "Missing")
```

**Current** (`child_clean.R:4894–4898`):
```r
non_error_codes <- c("Exclude-C-Identical",
                     "Exclude-C-Extraneous",
                     "Exclude-C-CF",
                     "Exclude-Missing",
                     "Exclude-Not-Cleaned")
```

**Analysis:** The reference lists `"Missing"` as a non-error code (correctly excluding missing measurements from both numerator and denominator) but does NOT list `"Not cleaned"`. `"Not cleaned"` is assigned to HC rows with agedays > 3 × 365.25 days (and HC ≥ 5y) — rows that are out-of-scope for cleaning, not measurement errors. In the reference, these rows would be counted as `n_errors` per `(subjid, HEADCM)` group, inflating the error-load ratio and potentially triggering error-load escalation for subjects with many HC > 3y measurements. The current correctly adds `Exclude-Not-Cleaned` to `non_error_codes`, matching the intent already expressed by including `"Missing"`. The fix is a side effect of the exclusion code rename + deliberate addition to the non_error_codes list.

**Pitfall:** **Factor levels / exclusion codes** — adding `"Missing"` but forgetting `"Not cleaned"` when both are assigned in the same preprocessing step.

---

## Items NOT flagged (audit trail)

### Step 19

- **DOP keyed lookup + `nomatch = NULL`**: `dop_snapshot[.(df$subjid[1], get_dop(df$param[1])), nomatch = NULL]` vs reference filter syntax `dop_snapshot[subjid == ... & param == ...]`. Equivalent behavior (both return 0 rows on no match). The `nomatch = NULL` addition restores the dead else branch reachability — F109 fix from Session 13 walkthrough, confirmed intentional.
- **`isTRUE()` wrapping** of pair and single conditions: reference uses `if (round_half_up(...) > 4 & ...)`, current uses `if (isTRUE(abs(...) > 4 & ...))`. The `isTRUE()` makes the condition NA-safe (avoids `if (NA)` error if z-scores are unexpectedly NA). Rounding removal is per procedure; `isTRUE()` is defensive improvement.
- **`for (i in 1:nrow(df))` → `for (i in seq_len(nrow(df)))`**: cosmetic equivalence (both iterate over valid row indices when `nrow(df) > 0`).
- **Dead variables removed** from reference pair branch: `abs_tbd.sd`, `abs_ctbd.sd` (computed but never used — pair conditions use `abs(diff_tbc.sd)`, not these); `med_dop`, `med_cdop` (computed but never used — conditions use `abs(diff_tbc.sd)` directly). Correctly removed in current.
- **`.SDcols = colnames(data.df)` → trimmed `.SDcols`**: reference passes all columns to closure; current trims to `c('index', 'internal_id', 'agedays', 'subjid', 'param', 'tbc.sd', 'ctbc.sd', 'exclude')` — F113 from Session 13 walkthrough.
- **`paste0 + table` → data.table by-group `.N`** for `only_single_pairs`: equivalent computation, same results. Known pattern from prior sessions.
- **`exclude_all` refresh after pair rule** (`exclude_all <- df$exclude`): identical in both — keeps the pair exclusion recorded before trimming to single for the 19H re-evaluation.
- **Index-keyed single write** (`exclude_all[ind_all == df$index] <- ...`): identical logic in both — preserves any pair exclusion already recorded.

### Step 21

- **`valid_set &` filter removal**: reference had `data.df[valid_set, ...]` with `valid_set = rep(TRUE, ...)` (all-TRUE, dead); current drops `valid_set` entirely — F123 from Session 14 walkthrough.
- **CF rescue code removal from `non_error_codes`**: reference listed old CF rescue exclusion codes in non_error_codes; current does not (rescued CFs are now Include rows). The net effect is the same — neither rescued CFs (now Include) nor excluded CFs (`Exclude-C-CF`, listed in current) count toward `n_includes` or `n_errors` differently than intended.

### Step 22

- **Dead `checkpoint_data` merge block removed**: reference contains `if (!"potcorr" %in% colnames(data.df) && exists("checkpoint_data") && ...)` inside an outer `if (has_potcorr && ...)` — the inner condition is always FALSE (the outer block already confirmed `"potcorr" %in% colnames(data.df)`), making the merge a no-op. Correctly removed in current.

---

## Open questions

None.

---

## Session 9 status

1 finding (AJ15 — Bug fix). Closed with no code change needed. Baseline unchanged at 63 / 48 / 28 / 41 / 13; no tests re-run since no code change.
