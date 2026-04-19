# R Child April-vs-January Comparison Procedure

**Last updated:** 2026-04-22

Procedure for comparing the **current** child R algorithm (April 2026) against the **reference** R implementation (January 2026, where R and Stata were well-reconciled) to identify logic changes — intentional or accidental — introduced during the reconciliation/refactor work since January.

This is a **pure R-to-R comparison.** Do NOT consult Stata for this exercise. The reference R is the authoritative post-reconciliation baseline.

---

## Purpose

Catch logic regressions or unintended drift introduced while bug-fixing, refactoring, renaming, and reorganizing the child algorithm since January 2026. Some divergences are known and intentional (CF rescue scheme, BIV z-parameter framework, exclusion-code rename, etc.); the goal is to surface the rest and decide what to do with each.

---

## Reference and current paths

**Reference (post-reconciliation baseline, Jan 2026):**
- `__Pipeline/current-infant-R-growthcleanr/growthcleanr/R/Infants_Main.R` (5011 lines) — main algorithm + most support functions
- `__Pipeline/current-infant-R-growthcleanr/growthcleanr/R/pediatric_support.R` (166 lines) — `valid()`, `temporary_extraneous()`, `swap_parameters()`, `as_matrix_delta()`, `na_as_false()`

**Skip:** `__Pipeline/current-infant-R-growthcleanr/growthcleanr/R/pediatric_clean.R` — older non-infant pediatric path (`cleanbatch()`), removed in current.

**Current (April 2026):**
- `__Pipeline/gc-github-latest/R/child_clean.R` (5027 lines) — everything consolidated into one file
- `__Pipeline/gc-github-latest/R/utils.R` (368 lines) — shared helpers (does NOT contain `.child_valid()`)

---

## Function/name mapping (reference → current)

| Reference | Current | Notes |
|-----------|---------|-------|
| `cleanbatch_infants()` | `cleanchild()` | Per-batch main algorithm |
| `temporary_extraneous_infants()` | `identify_temp_sde()` | Shared SDE helper |
| `temporary_extraneous()` (pediatric_support.R) | `identify_temp_sde()` | Old generic path consolidated |
| `calc_oob_evil_twins()` | `calc_otl_evil_twins()` | "out of bounds" → "out of line" rename |
| `valid()` (pediatric_support.R) | `.child_valid()` | Renamed to prevent collision; gained `.` prefix |
| `swap_parameters()` (pediatric_support.R) | `get_dop()` | DOP lookup; reference also has a `get_dop()` in `Infants_Main.R` — confirm which the infant path called before walking each step |
| User `id` (sort key, tiebreaker) | `internal_id` | Sequential integer created at entry; user `id` preserved untouched |
| `cleanbatch()` (pediatric_clean.R) | (removed) | Legacy non-infant path |

---

## Known intentional changes — do NOT flag

These are documented design changes since reconciliation. If a diff matches one of these, log briefly as "intentional - {reason}" and move on without further analysis.

### Algorithm-level

- **CF rescue scheme (Step 6):** Lookup-table-based `cf_rescue` parameter (`standard`/`none`/`all`), `cf_rescued` output column, optional `cf_detail` columns (`cf_status`, `cf_deltaZ`). Replaces fixed 0.05 / 0.10 thresholds. Documented 2026-04-14 (`cf-rescue-thresholds.md`).
- **BIV parameters (Step 7):** Eight `biv.z.*` per-cell parameters replacing single `sd.extreme` / `z.extreme`. Documented 2026-04-17.
- **Error load threshold (Step 21):** Now uses `error.load.threshold` parameter (default 0.5) instead of hardcoded `.4`. Bug fix, documented 2026-04-12.
- **Adult algorithm integration:** Adult dispatch in `cleangrowth()` and adult-specific output columns (`mean_ht`, `bin_result`). Do not flag adult-related divergence.

### Naming and exclusion codes

- **Exclusion code rename:** All child codes now use `Exclude-C-{Reason}` format (not param-specific). `.child_exc()` helper generates codes. Documented 2026-04-14/16.
- **`Missing` → `Exclude-Missing`** and **`Not cleaned` → `Exclude-Not-Cleaned`**.
- **`valid()` → `.child_valid()`** rename (prevents collision with legacy `valid()`).
- **`internal_id` introduction:** Sequential integer for all internal sorting/tiebreaking. User `id` preserved untouched. Documented 2026-04-12.
- **`oob` → `otl`** rename in Evil Twins terminology (`calc_oob_evil_twins` → `calc_otl_evil_twins`).

### Removed (legacy / dead code)

- **Legacy pediatric path:** `pediatric_clean_legacy.R`, `pediatric_support_legacy.R`, `R/deprec/`, `R/modified_source_code/`. Documented 2026-04-16.
- **Legacy parameters:** `use_legacy_algorithm`, `prelim_infants`, `lt3.exclude.mode`, `ewma.exp`, `recover.unit.error`. Documented 2026-04-16.
- **Round-tolerance helpers:** `round_stata()` and intermediate z-score rounding removed (was a Stata-reconciliation artifact). **Per user instruction, do NOT flag rounding-tolerance removal.**
- **Fenton z-score data:** `fentlms_forz` deleted; switched to `fenton2025_ms_lookup_smoothed.csv` (CSD method). Documented 2026-04-16.

### Reorganization (cosmetic, not logic)

- **File consolidation:** Reference's `Infants_Main.R` + `pediatric_support.R` consolidated into current `child_clean.R`. Cosmetic; not a logic change.
- **`as_matrix_delta` made internal** (no longer exported).
- **`ewma_cache_init()` / `ewma_cache_update()` added** — incremental EWMA caching for batch performance. Walked in Session 7b (2026-04-18). Reference uses only `ewma()` directly.
- **`gc_preload_refs()` added** — pre-loads reference closures for repeated calls.
- **`cached_results` / `changed_subjids` added** — partial re-run support.

---

## Categorization scheme

For each diff that is NOT a known intentional change above, assign one of:

| Category | Definition | Action |
|---|---|---|
| **Bug fix** | Reference R didn't make sense; current R does. | Note briefly that it's a fix; no action needed. |
| **New error** | Reference R made sense; current R doesn't. | Flag for discussion → fix after approval. |
| **Unclear** | Could be either, or insufficient context to decide. | Flag for discussion → resolve with Carrie. |
| **Intentional (other)** | Documented elsewhere or discovered during walkthrough to be intentional. | Note briefly; no action. |

---

## Per-session workflow

1. **Pre-session git hygiene check.** Before doing anything else, run `git status` and `git log origin/<branch>..HEAD --oneline` to detect uncommitted modifications, untracked files, or local commits not yet pushed. Also check for accidental work in stale `/Users/Shared/` copies. If anything is found:
   - Surface to Carrie immediately with a one-screen summary (which files, what session/scope they likely belong to)
   - Confirm whether to commit/push first OR knowingly bundle into this session's commit
   - Do NOT silently inherit prior-session state — at minimum, document the state in the new session's walkthrough note so the eventual commit message can attribute correctly
   - This step exists because Sessions 12-14 of the algorithm walkthrough series silently went uncommitted; the catch-up was painful (commit `fa7f0df`)
2. **Pick step(s) for the session** (see Step Ordering below). Confirm with Carrie at start.
3. **Locate line ranges** in reference (`Infants_Main.R` / `pediatric_support.R`) and current (`child_clean.R`) for the step.
4. **Diff** the two ranges. Either:
   - `Read` both ranges and compare in head, or
   - Use `diff` via Bash (`diff -u <reference> <current> | less`) and walk hunks.
5. **For each diff,** explicitly walk the Logic Pitfalls checklist (below), then classify:
   - Function rename / file consolidation / known intentional change → log briefly, move on.
   - Otherwise → categorize (Bug fix / New error / Unclear / Intentional (other)) and record context, including which Logic Pitfalls category (if any) the finding falls into.
6. **Cross-check support functions** if the step calls them (e.g., `ewma()`, `identify_temp_sde()`, `get_dop()`). Either compare in line or defer to the dedicated support-function session (Session 10).
7. **Batch findings** at end of session in the dated walkthrough note (per `algorithm-walkthrough-procedure.md` pattern). Each session's walkthrough section should include:
   - Pre-session baseline test counts
   - Scope table (sub-area / reference lines / current lines / what)
   - Known intentional changes encountered (logged briefly, no further analysis)
   - Findings (numbered, categorized, with file:line refs and Logic Pitfalls category)
   - **Items NOT flagged** subsection — list items reviewed and explicitly chosen not to open as findings; provides an audit trail showing what coverage was attempted
   - Open questions for Carrie (if any)
8. **Update findings index** (`R-child-AprJanComparison-findings.md`) with one-liner per finding. Use `AJ##` prefix (Apr/Jan) to distinguish from the walkthrough series' `F##` numbering. Number monotonically across all sessions of this comparison series (Session 2 starts at AJ8).
9. **Carrie reviews and approves** findings.
10. **Fix approved findings in the same session.**
11. **Re-run tests** to confirm no regression. Reinstall package first; verify counts match the documented baseline (see "Tests to run after each fix" below).
12. **Update CLAUDE.md notes.** Prepend a one-paragraph session summary to the "Last updated:" entry in both `gc-github-latest/CLAUDE.md` and `__Pipeline/CLAUDE.md`. Keep concise; name the session number, scope, key findings (with status), and final test counts. The "Last updated:" date may be kept at the most recent walkthrough's date if older.
13. **Stage, commit, and push.** Single commit per session, with title `R-vs-R Apr/Jan Comparison Session N: <scope>`. Body should list any code changes (with line refs), the AJ## findings touched, and final test counts. Push to `origin/<branch>` immediately. Confirm via `git status` that the working tree is clean and `git log origin/<branch>..HEAD` is empty before ending the session — closes the loop on step 1.

---

## Step ordering (proposed)

| Session | Coverage | Complexity | Notes |
|---|---|---|---|
| 1 | Preprocessing: z-scores, GA correction (Step 2b), recentering | High | Affects all downstream steps; do first |
| 2 | Step 5 (Temp SDE) + Early Step 13 (SDE-Identicals) | Medium | Linked logic |
| 3 | Step 6 (CF) | High | Major intentional change (new CF rescue) — narrower comparison since most differences are intentional |
| 4 | Step 7 (BIV) + Step 9 (Evil Twins) | Low–Medium | Step 7 has known intentional change (`biv.z.*` params) |
| 5 | Step 11 (EWMA1) | High | Complex closure, EWMA caching |
| 6 | Step 13 (Final SDE) | High | Multi-phase (B1 / B2 / B3) |
| 7 | Step 15 (EWMA2) + Step 16 (Birth HT/HC) | High | Linked logic |
| 8 | Step 17 (Velocity) | High | Complex per-iteration loop |
| 9 | Step 19 (Pairs/Singles) + Step 21 (Error Load) + Step 22 (Output) | Low | Quick steps |
| 10 | Support functions: `ewma()`, `ewma_cache_*`, `as_matrix_delta()`, `identify_temp_sde()`, `calc_otl_evil_twins()`, `get_dop()`, `calc_and_recenter_z_scores()`, `.child_valid()`, `.child_exc()`, `sd_median()`, `read_anthro()` | Variable | Some already walked in earlier walkthrough sessions; reference back rather than re-walk in depth |

Approximate total: 8–10 sessions.

---

## Output files

**Per-session walkthrough notes:** `__Pipeline/gc-github-latest/walkthrough-todo-YYYY-MM-DD.md` — same pattern as existing notes (matches `algorithm-walkthrough-procedure.md`). Add a top-level `# R-vs-R comparison — Session N — YYYY-MM-DD` section with batched findings.

If the day's filename already exists (e.g., a regular walkthrough session was done earlier the same day), append the R-vs-R section to that file after a horizontal rule (`---`) separator. Do NOT create a parallel file with a suffix — keep one walkthrough-todo file per calendar day.

**Cross-session findings index:** `__Pipeline/gc-github-latest/R-child-AprJanComparison-findings.md` — single growing markdown file. One section per session, one bullet per finding:

```
- [AJ##] Step N: brief description — Category — Status — file:line refs (or commit ref if fixed)
```

Status values: `open` → `approved` → `fixed` → `closed`. Items can also go directly to `closed` when no change is needed (e.g., bug fixes already present in current code, intentional changes confirmed by Carrie).

**CLAUDE.md updates:** at end of each session, prepend a one-paragraph summary to the "Last updated:" entry in both `gc-github-latest/CLAUDE.md` and `__Pipeline/CLAUDE.md`. The walkthrough-todo file is the authoritative detail; the CLAUDE.md note is a pointer + key takeaways.

---

## Approval gate

**Carrie must approve every change before it is made**, even those that look like clear errors. Findings are batched at the end of each session for review. Do NOT make code changes during the comparison phase of the session.

---

## Tests to run after each fix

After fixing approved findings, **reinstall the package** then run the child test suite:

```bash
Rscript -e 'devtools::install_local(".", force=TRUE, upgrade="never")'

NOT_CRAN=true Rscript -e \
  'testthat::test_file("tests/testthat/test-child-regression.R")'
NOT_CRAN=true Rscript -e \
  'testthat::test_file("tests/testthat/test-child-algorithms.R")'
NOT_CRAN=true Rscript -e \
  'testthat::test_file("tests/testthat/test-child-parameters.R")'
NOT_CRAN=true Rscript -e \
  'testthat::test_file("tests/testthat/test-child-edge-cases.R")'
```

Adult tests do not need to be re-run unless an adult file was touched.

Baseline child counts (as of Session 14, 2026-04-22): **63 / 48 / 28 / 41 / 13** across the five child test files. If counts change unexpectedly, treat as a regression and investigate before continuing.

---

## Logic pitfalls to watch for at each diff

These are categories of subtle bugs that have surfaced in prior walkthroughs. When inspecting a diff, explicitly check each of these — most "New error" findings will fall into one of them. Per workflow step 5, name the pitfall category in each finding so the index can be cross-referenced over time.

### Boundary changes

- **`<` vs `<=`, `>` vs `>=`** — e.g., Step 17's `> 365.25` vs `>= 365.25` for the 1-year gap (F95); Step 17's tied-`abs(tbc.sd)` `>` vs `>=` (F97 — strict `>` flagged neither row in a tied pair, letting violations escape entirely).
- **Numeric thresholds** — e.g., `0.4` vs `0.41`, `0.5` vs `0.50`, hardcoded constants vs configurable parameters.
- **Age cutoffs** — `30 months` using `tanner.months` (midpoint-based) vs `agedays/30.4375` (this-row-only) (F94 — boundary-crossing pairs were keyed on the wrong age).
- **Closed vs open intervals** in lookup tables — interval `[20, 46)` vs `[20, 46]`; check both ends.

### Sort order

- **`setkey(data.df, subjid, param, agedays, ?)`** — `id` vs `internal_id` as final tiebreaker.
- **Missing `setkey` after row modification** — steps that assume rows are pre-sorted may silently break if a prior step added/removed rows.
- **Per-group order** within `by = ...` — relies on key being set before the operation.

### by-group scope

- **`by = subjid` vs `by = .(subjid, param)`** — silent wrong results.
- **`by = .(subjid)` over multi-param data** — collapses across params incorrectly.
- **`by = sp_key` vs `by = subjid`** — sp-keyed groups miss param boundaries; subj-keyed groups can leak across params.

### .child_valid() flag scope

- **Which `include.*` flags are passed** (`include.temporary.extraneous`, `include.extraneous`, `include.carryforward`).
- A step that should see CFs but doesn't pass `include.carryforward = TRUE` will silently miss them (or vice versa).
- A step that should NOT call `.child_valid()` at all (e.g., Step 21 needs to count error-coded rows) but does will silently undercount.

### Tiebreaker direction

- **"lowest internal_id wins"** vs **"highest internal_id wins"** — easy to flip.
- **At birth (`agedays == 0`)**: keep lowest id (pre-postnatal shift). All other ages: keep highest id. Reversal at the birth boundary is a known pattern; verify if present.

### Cache / snapshot freezing

- **DOP snapshots** taken before per-(subjid, param) closures must include both `tbc.sd` AND `ctbc.sd`.
- A snapshot taken AFTER the first param's exclusions invalidates the second param's lookup (F105 pattern in Step 19).
- **Static vs dynamic columns**: columns merged pre-loop on static data must not be re-used per-iteration if they depend on per-iteration variables (F101 pattern — `whoagegrp.ht` was static but `d_agedays` was per-iteration).

### Factor levels / exclusion codes

- Codes assigned to `exclude` must exist in `exclude.levels` — **unlisted strings become silent NAs.**
- Verify code names match between assignment and downstream filters: `Exclude-C-Pair` vs an older `Exclude-Carried-Forward-Pair` will silently fail to match.
- `.child_exc()` constructor: passes `param` argument but ignores it (codes are not param-specific in current).

### Per-iteration column lifecycle

- Columns set inside an iteration loop must be **reset to NA at top of each iteration** if they should not carry over (F102 — `whoinc.age.ht` reset added as defense-in-depth).
- Static columns (set pre-loop) vs dynamic columns (set per iteration) — easy to confuse, and easy to mistakenly use a static value in a dynamic context.

### NA / empty-set handling

- **`min()` / `max()` on empty subset** → `Inf` / `-Inf` with a warning. `suppressWarnings` + `is.infinite` is the established pattern.
- **`which.min` / `which.max` on empty vector** → `integer(0)`.
- **`median(NA)`** → `NA`; downstream comparisons silently produce `NA → FALSE`.
- **Partial keyed data.table lookup without `nomatch = NULL`** → returns 1 NA row when no match. This can mask dead `else` branches that look reachable (F109 pattern in Step 19).

### data.table reference semantics

- **`:=` modifies in place** — defensive `copy()` may have been omitted or added unnecessarily; either can be a bug.
- **Joined tables** — did the merge create new rows, modify the original, or return a copy?
- **Merge-by**: `merge(x, y, by = "subjid")` vs `x[y, on = "subjid"]` — different semantics around row preservation.

### Parameter scope

- HC has narrower applicable age range (≤ 3y); does the step correctly skip HC where intended?
- HT-only or WT-only logic correctly param-gated?
- Steps that loop over the three params: `setdiff` / `intersect` correctness when one param has no data for a subject.

### Interaction with later steps

- Variables stored on `data.df` for downstream use must persist past the step (no premature drop).
- Variables not needed downstream should be dropped before output (verify drop list).
- Order of step execution matters: a step that depends on a variable set by a later step is broken.

---

## Procedural pitfalls

- **Function renames cause spurious diffs.** Map names first (see table above) so you can recognize what's "the same function with a new name" vs. "different logic."
- **Known intentional changes WILL appear as diffs.** Resist the urge to investigate them in depth — log briefly and move on. The known list above should cover most.
- **Reference uses `valid()`, current uses `.child_valid()`.** Otherwise functionally equivalent; treat the rename as cosmetic.
- **Reference may have orphan / dead code paths** (e.g., second `valid()` in `Infants_Main.R` at line 4991 is likely orphaned — `pediatric_support.R`'s `valid()` is the live one). Don't compare orphan code to current.
- **Reference contains `pediatric_clean.R`** for the older non-infant path. **Skip it.** Comparison scope is `Infants_Main.R` + `pediatric_support.R` only.
- **Don't flag rounding-tolerance removal** (`round_stata()`, intermediate z-score rounding). Per user instruction.
- **Don't flag adult-related changes** (adult dispatch, output columns, etc.) — adult was not in reference.

---

## Related procedures

- `algorithm-walkthrough-procedure.md` — code/comments/narrative/tests walkthrough procedure (different scope; deeper but doesn't compare against reference).

The current procedure is **complementary** to the walkthrough procedure: walkthroughs find issues by inspecting current code in isolation; this procedure finds issues by diffing against a known-good baseline.
