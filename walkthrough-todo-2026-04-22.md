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
