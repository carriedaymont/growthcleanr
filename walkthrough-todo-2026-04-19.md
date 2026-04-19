# GC Walk-Through — 2026-04-19: Child Step 15/16 (EWMA2 Moderate + Birth HT/HC)

Dedicated detailed walkthrough of Main Child Steps 15 and 16 (EWMA2 / moderate trajectory outliers) in `R/child_clean.R` (~lines 3660–4088 pre-edit). Adult algorithm has no z-score pipeline and has its own separate moderate-EWMA logic in `adult_clean.R` / `adult_support.R` (Adult Step 11Wb); out of scope here. This pass is child-scoped.

Comment / narrative / small-cleanup with two targeted behavior-neutral tightenings:

- Removed `include.temporary.extraneous = TRUE` at `child_clean.R:3677` (all temp SDEs are resolved in Main Child Step 13; a straggler should be treated as ineligible, not counted — per Carrie).
- Tightened the Step 16 birth filter from subject-level (`subjid %in% subj_with_birth`) to per-sp_key (`sp_key %in% sp_with_birth`) at `child_clean.R:3948–3987`. A subject with birth HT but no birth HC (or vice versa) now has only the birth-carrying param processed in Step 16. The non-birth-carrying param was previously included in the filter but never produced an exclusion (all Step 16 rules require `agedays == 0`), so this is a pure efficiency improvement with no change in exclusion set.

Child test counts unchanged.

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
| Step 15/16 header | `child_clean.R` 3660–3668 | Narrative-style description of global iteration approach |
| Pre-loop setup | 3670–3679 | Re-sort; `sp_key` build; `sp_to_process_15` count filter |
| 15A | 3681–3686 | p_plus / p_minus = ±5% (WT) or ±1 cm (HT/HC) perturbations |
| 15B | 3688–3702 | tbc.p_plus / tbc.p_minus via `calc_and_recenter_z_scores()` ×2 |
| first_meas pre-calc | 3704–3718 | Per-param first-Include-row indicator used by Step 15's `first` rules |
| Step 15 iteration | 3720–3939 | Global while loop: DOP snapshot → recompute filter + first_meas → per-(subjid, param) closure (cache-hit or full ewma_cache_init) → addcrit → 9 rule families (middle / birth-WT ×2 / first ×2 / last ×6) → worst-candidate selection → `sp_with_new_excl` next-iter set |
| Step 16 pre-loop | 3941–3970 | Subject-param-level `sp_with_birth` filter; tbc.sd range > 1 pre-filter |
| Step 16 iteration | 3972–4078 | Same skeleton as Step 15; 4 rule families (birth-HT-HC ×2 gap tiers); no DOP lookup |
| Cleanup | 4080–4088 | Drop `sp_key`, `p_plus`, `p_minus`, `tbc.p_plus`, `tbc.p_minus`, `first_meas` |

**Inventory of support functions used:** `calc_and_recenter_z_scores()` (pre-loop only), `ewma_cache_init()` / `ewma_cache_update()` (inside both iteration loops), `get_dop()` (Step 15 `last-ext` rules only, line 3845), `.child_valid()` (filter recomputation), `.child_exc()` (exclusion-code generation). `ewma()` is NOT directly called — Steps 15 and 16 use the incremental cache API rather than the full-recompute path used by Steps 11/13/17. (Narrative Code-location cell updated — F82.)

---

## Fix-now items

### F73. `# Add id for consistent SDE order` stale comment at the pre-loop sort — FIXED
- **File/lines:** `R/child_clean.R:3670–3672` (pre-edit).
- **Issue.** Comment read:
  ```
  # Order data for processing
  # Add id for consistent SDE order
  data.df <- data.df[order(subjid, param, agedays, internal_id),]
  ```
  Two problems: (a) this is Step 15 (moderate EWMA), not SDE; (b) the comment says "id" but the `order()` key is `internal_id` — same loose wording fixed in Session 10 F68 for the Step 13 analogue. Session 10 left this specific site for the Step 15/16 walkthrough.
- **Fix.** Rewritten to name the downstream dependency:
  ```
  # Re-sort so the per-group closures below walk rows in
  # (subjid, param, agedays, internal_id) order; prior/next neighbor
  # differences and first_meas depend on this sort.
  ```
- Same `# Add id for consistent SDE order` pattern also appears at `child_clean.R:4148` (Step 17 sort) — out of scope here; belongs to the Step 17 walkthrough.

### F74. `# Include temp SDEs in count so they get p_plus/p_minus calculated` stale comment — FIXED
- **File/lines:** `R/child_clean.R:3674–3679` (pre-edit).
- **Issue.** Pre-edit block read:
  ```
  # Pre-identify subject-params that need processing: Include values, >2 values
  # Include temp SDEs in count so they get p_plus/p_minus calculated
  data.df[, sp_key := paste0(subjid, "_", param)]
  include_counts <- data.df[.child_valid(data.df, include.temporary.extraneous = TRUE),
                            .(n_include = .N), by = sp_key]
  sp_to_process_15 <- include_counts[n_include > 2, sp_key]
  ```
  The comment is misleading on two counts:
    1. Temp SDEs are all resolved in Main Child Step 13 (Phase B4 at `child_clean.R:3650–3653` overwrites every Include / Temp-Same-Day row with one of Include / `Exclude-C-Identical` / `Exclude-C-Extraneous`). By the time we reach Step 15, no row has `exclude == "Exclude-C-Temp-Same-Day"`, so the `include.temporary.extraneous = TRUE` flag is a no-op.
    2. Even if a temp SDE survived, the flag's effect would only be to raise the per-sp_key count by 1. The p_plus / p_minus assignment at lines 3684–3686 filters by `sp_key %in% sp_to_process_15 & param == <...>`, not by exclude status, so *every* row in a qualifying sp_key gets p_plus / p_minus regardless of the flag. The comment's causal claim ("so they get p_plus/p_minus calculated") is backwards.
- **Fix per Carrie's request.** Dropped the `include.temporary.extraneous = TRUE` flag entirely (use `.child_valid()` default). Rewrote the comment to describe current state:
  ```
  # Pre-identify subject-params that need processing: Include values only,
  # 3+ per subject-param. All temp SDEs were resolved in Main Child Step 13,
  # so the default `.child_valid()` (temp SDEs excluded) is all that is
  # needed here; a temp SDE leaking through should be treated as ineligible.
  ```
  Rationale: the invariant is enforced by Main Step 13, and Carrie preferred that any straggler *not* be counted (rather than kept as a defensive pad).

### F75. `# Fix first_meas logic` stale-changelog line removed — FIXED
- **File/lines:** `R/child_clean.R:3704–3718` (pre-edit).
- **Issue.** Pre-edit block led with:
  ```
  # Pre-calculate first_meas indicator
  # Fix first_meas logic
  # Handle HT/HC differently from WT
  # "Non-birth first" means: the first value, when that first value is not at birth
  # For WEIGHTKG: birth IS included in Step 15, so first_meas only if position 1 has agedays > 0
  # For HEIGHTCM/HEADCM: birth is EXCLUDED from Step 15, so calculate position among non-birth only
  ```
  `# Fix first_meas logic` is the F58–F60 stale-changelog pattern — it describes that the code was changed, not what the code does. The remaining lines are current-state descriptions but slightly redundant.
- **Fix.** Consolidated the comment block into a single current-state description of the two param-specific definitions and why they differ. Inline per-param one-liners above each `data.df[...]` assignment also removed as redundant with the consolidated block.

### F76. `valid()` stale reference in Step 15 inner-loop filter comment — FIXED
- **File/lines:** `R/child_clean.R:3758–3774` (pre-edit).
- **Issue.** Pre-edit comments read `# RECALCULATE filter each iteration - valid() depends on current exclusions / # Must recalculate filter inside loop`. `valid()` is the legacy name; the current helper is `.child_valid()`. Same inline per-param first_meas one-liners as F75.
- **Fix.** Stale `valid()` → `.child_valid()` in comment. Consolidated first_meas recalc comments to match the F75 structure and added a pointer that "the prior iteration may have excluded what was previously the first non-birth Include row" (which is why it's recomputed rather than cached).

### F77. `temporarily excluded before iteration starts` misleading phrasing — FIXED
- **File/lines:** `R/child_clean.R:3909–3911` (pre-edit).
- **Issue.** Pre-edit comment read `# All pot_excl candidates are eligible - birth HT/HC is already excluded / # from Step 15 processing (temporarily excluded before iteration starts)`. "Temporarily excluded" is misleading — birth HT/HC rows are not given a temporary exclusion code; they are simply filtered out of `step15_filter` via `!((param == "HEIGHTCM" | param == "HEADCM") & agedays == 0)` at lines 3725 and 3762, and so never enter `df` in the per-group closure.
- **Fix.** Rewritten to:
  ```
  # Exclude the worst candidate (highest abssum). Birth HT/HC rows are
  # already filtered out of step15_filter upstream, so they never enter
  # df and cannot be candidates here.
  ```

### F78. Tiebreaker comment consolidation in Step 15 closure — FIXED
- **File/lines:** `R/child_clean.R:3914–3918` (pre-edit).
- **Issue.** Three lines of comment explained the same thing:
  ```
  # Use internal_id tie-breaker for deterministic selection
  # which.max returns first tie, which depends on data order
  # Use order() with internal_id as tie-breaker for reproducibility
  ```
  Lines 1 and 3 are the same claim stated two ways; line 2 is the rationale for using `order()` instead of `which.max`.
- **Fix.** Collapsed to two lines stating claim + rationale once:
  ```
  # internal_id breaks ties on abssum (rare but possible — abssum
  # is continuous). order() is deterministic; which.max is not
  # when rows tie.
  ```

### F79. `valid()` stale reference in Step 16 inner-loop filter comment — FIXED
- **File/lines:** `R/child_clean.R:3981–3987` (pre-edit).
- **Issue.** Same `valid()` / "Must recalculate filter inside loop" stale-wording pair as F76, in the Step 16 iteration loop. Pre-edit filter still used `data.df$subjid %in% subj_with_birth` (the pre-fix subject-level set) — updated below as part of the Step 16 birth filter tightening.
- **Fix.** Updated comment to current-state and replaced the subject-level check with a note that `sp_to_process` already encodes the per-sp_key birth requirement:
  ```
  # Recompute step16_filter each iteration — .child_valid() depends on
  # current exclusions, which change as rows are flagged. sp_to_process
  # is already the surviving per-sp_key set, so we don't need to repeat
  # the sp_with_birth check here.
  ```

### F80. Step 16 `# Use id tie-breaker` stale comment — FIXED
- **File/lines:** `R/child_clean.R:4061` (pre-edit).
- **Issue.** Pre-edit comment: `# Use id tie-breaker for deterministic selection`. Same F68 pattern — `order(-candidates$abssum, candidates$internal_id)` on the line below uses `internal_id`, not user `id`.
- **Fix.** Rewritten to: `# internal_id breaks ties on abssum (same convention as Step 15).` Cross-refs Step 15 F78 so the reader doesn't need the rationale repeated.

### F81. Drop-columns comment at end of Steps 15/16 — FIXED
- **File/lines:** `R/child_clean.R:4083–4088` (pre-edit).
- **Issue.** Pre-edit comment read:
  ```
  # Drop columns no longer needed after Step 15
  # p_plus, p_minus, tbc.p_plus, tbc.p_minus: only used in EWMA2 +/-5% rule (Step 15)
  # first_meas: only used in EWMA2 first-measurement exclusion logic (Step 15)
  cols_to_drop_15 <- intersect(c("p_plus", "p_minus", ...
  ```
  Two problems: (a) the drop runs after Step 16 completes, not after Step 15 — Steps 15 and 16 share the addcrit perturbation columns; (b) "+/-5% rule (Step 15)" is a misnomer — the perturbation is ±5% for WT but ±1 cm for HT/HC, and it's the addcrit check used by both Steps 15 and 16, not a dedicated "rule". `tbc.p_plus` and `tbc.p_minus` are both explicitly in Step 16's `.SDcols` at line 4070, and the Step 16 closure reads them at lines 4034–4037.
- **Fix.** Rewritten comment to describe current-state role of each column:
  ```
  # Drop columns no longer needed after Steps 15-16.
  #   p_plus / p_minus / tbc.p_plus / tbc.p_minus: feed the addcrit
  #     perturbation check used by both Step 15 and Step 16.
  #   first_meas: only used by Step 15's first-row rules (Step 16 keys
  #     off agedays == 0 instead).
  ```
  Local variable renamed `cols_to_drop_15` → `cols_to_drop_15_16` for consistency with the block it belongs to.

### F82. Narrative Code-location cell — remove `ewma()`, add `get_dop()` — FIXED
- **File/lines:** `child-algorithm-reference.md:1190` (pre-edit).
- **Issue.** Pre-edit cell read `Inline in cleanchild() in child_clean.R; uses ewma(), ewma_cache_init(), and ewma_cache_update() (also defined in child_clean.R); z-score helpers calc_and_recenter_z_scores() for p_plus/p_minus conversion`. But `ewma()` is not directly called by Step 15 or Step 16 — both use the incremental cache API (`ewma_cache_init()` + `ewma_cache_update()`) exclusively. `get_dop()` is called at `child_clean.R:3845` for the Step 15 `last-ext` / `last-ext-high` DOP check and is not listed.
- **Fix.** Rewritten cell:
  ```
  Inline in cleanchild() in child_clean.R. Uses
  calc_and_recenter_z_scores() for the p_plus / p_minus pre-loop,
  ewma_cache_init() / ewma_cache_update() for the incremental EWMA
  inside the iteration loop, and get_dop() for the Step 15 last-ext
  DOP lookup.
  ```

### F83. "±5% rule" narrative misnomer — FIXED
- **File/lines:** `child-algorithm-reference.md:1205` (pre-edit).
- **Issue.** Pre-edit sentence read `The ±5% rule tests whether the measurement would still be an outlier if perturbed slightly — a robustness check that reduces false positives for values near exclusion boundaries.` But:
    1. ±5% is the WT-specific perturbation; HT/HC uses ±1 cm. Calling it "the ±5% rule" ignores the HT/HC case.
    2. It's not a specific rule — it's the addcrit perturbation check that wraps every Step 15/16 exclusion rule.
- **Fix.** Rewritten to: `The addcrit perturbation check (detailed below) tests whether the measurement would still be an outlier if nudged slightly toward its neighbors — a robustness check that reduces false positives for values near exclusion boundaries.`

### F84. "lowest `id` as tiebreaker" narrative — FIXED
- **File/lines:** `child-algorithm-reference.md:1285` (pre-edit).
- **Issue.** Pre-edit `### Worst-value selection` section read `Both steps use abs(tbc.sd + dewma.all) as the sort key (same as EWMA1), with lowest id as tiebreaker.` Two issues: (a) `id` vs `internal_id` (F68 pattern); (b) direction ambiguous — a reader can't tell whether "lowest id" is kept or excluded.
- **Fix.** Rewritten to:
  ```
  Both steps use abs(tbc.sd + dewma.all) as the sort key (same as
  EWMA1). Ties on the sort key are broken by lowest internal_id —
  the row with the lowest internal_id among the tied candidates
  becomes the exclusion. Ties are rare because abs(tbc.sd +
  dewma.all) is continuous, but internal_id makes the selection
  deterministic.
  ```

### F85. DOP snapshot refresh narrative — Step 15 only — FIXED
- **File/lines:** `child-algorithm-reference.md:1289` (pre-edit).
- **Issue.** Pre-edit `### Global iteration` closing line read `DOP snapshot refreshed each iteration.` But DOP is only used in Step 15 (`last-ext` / `last-ext-high` rules at `child_clean.R:3899–3906`). Step 16 has no DOP lookup (lines 4046–4056 are all agedays == 0 checks with no `tbc_dop` reference).
- **Fix.** Rewritten to: `Step 15's DOP snapshot (used by the last-ext and last-ext-high rules) is refreshed each iteration; Step 16 has no DOP lookup.`

---

## Behavior-neutral tightenings requested by user

### Removal of `include.temporary.extraneous = TRUE` at `child_clean.R:3677`
Covered by F74 above. Carrie's rationale: even if a temp SDE somehow survived Main Step 13, she would not want it counted toward the n ≥ 3 threshold for EWMA2 eligibility; an invariant break upstream should manifest as "too few Include rows", not as "counted as Include."

### Step 16 birth filter: subject-level → sp_key-level
- **File/lines:** `R/child_clean.R:3948–3987` (pre-edit).
- **Issue.** Pre-edit code identified `subj_with_birth` as the set of subjids with any birth HT or HC Include row, then filtered `step16_filter` to rows of those subjects that are HT or HC Include. Consequence: a subject with birth HT but no birth HC had both HT *and* HC included in step16_filter, so the HC group was handed to the per-group closure. All four Step 16 rules require `agedays == 0`, so the HC closure produced no exclusion (hit the "`candidates <- df[pot_excl != ""]`; `if (nrow(candidates) > 0)`" gate with zero candidates) and the HC sp_key would drop from `sp_to_process` on the very next iteration. The cost was one wasted `ewma_cache_init()` call per such HC group per first iteration — behavior-neutral but O(n²) in group size.
- **Fix.** Changed both initial and in-loop filters to use `sp_with_birth` (per-sp_key) instead of `subj_with_birth` (per-subject):
  - Initial filter (`sp_with_birth <- unique(data.df[... & agedays == 0, sp_key])`) at `child_clean.R:3945–3950` selects only HT/HC sp_keys that themselves have a birth measurement.
  - In-loop filter at `:3981–3987` no longer re-checks `sp_with_birth` directly — `sp_to_process` is already the surviving set after pre-filter and previous iterations, which inherits the birth requirement transitively. Filter now reads:
    ```
    step16_filter <- data.df$sp_key %in% sp_to_process &
      data.df$param %in% c("HEIGHTCM", "HEADCM") &
      .child_valid(data.df)
    ```
- **Result.** Identical exclusion set for all test cases (verified — child test counts unchanged), less work per first iteration.
- **Narrative alignment.** Narrative at `child-algorithm-reference.md:1281` ("Only processes HT/HC subject-params with a birth measurement...") already described the tightened intent correctly; the pre-fix code was looser than the narrative. Added a one-line clarification of the asymmetry case ("A subject with birth HT but no birth HC (or vice versa) has only the birth-carrying param processed here.").

---

## Checklist items applied (summary)

Applied the 20-item cross-cutting checklist from `algorithm-walkthrough-procedure.md`. Items 14–19 (wrapper-only) are n/a for Steps 15/16. Items checked and their dispositions:

1. **Sort order determinism (both).** `order()` / `setkey()` call sites in Steps 15/16:
    - `data.df[order(subjid, param, agedays, internal_id),]` at `child_clean.R:3672` (pre-loop sort)
    - `setkey(dop_snap, subjid, param)` at `:3756` (DOP snapshot — subjid+param key only, no agedays since the DOP lookup join key is agedays)
    - `order(-candidates$abssum, candidates$internal_id)` at `:3918` (Step 15 worst-candidate)
    - `order(-candidates$abssum, candidates$internal_id)` at `:4062` (Step 16 worst-candidate, post-F80)
    - All use `internal_id` as the final tiebreaker. Correct.
    - **Birth tiebreaking** — not applicable in the SDE sense. Step 15 explicitly filters birth HT/HC out of its scope (the birth-WT rules handle agedays == 0 for WT). Step 16 only processes agedays == 0 rows. No SDE-style "keep lowest at birth / keep highest otherwise" logic is needed because these are distinct-ageday EWMA candidates, not same-day duplicates.

2. **Z-score correctness (child).** All z-score work uses the correct column:
    - `tbc.sd` — primary z-score for dewma (`dewma.all = tbc.sd - ewma.all` at `:3823` and `:4029`), addcrit neighbor diffs, `tbc_diff_next` / `tbc_diff_prior`, last-row `tbc_prev`, DOP check in `last-ext`, worst-candidate `abssum`.
    - `tbc.p_plus` / `tbc.p_minus` — recentered perturbation z-scores used only in the addcrit cross-neighbor checks (`tbc_diff_plus_next` / `tbc_diff_minus_next` etc.). Computed once in the pre-loop (`:3696–3702`) via `calc_and_recenter_z_scores()` with shared `measurement.to.z_who` closure.
    - `ctbc.sd` — corrected z-score used for `c.dewma.all = ctbc.sd - c.ewma.all`. Fed into `ewma_cache_init()` alongside `tbc.sd`; the cache auto-detects `skip_ctbc` when `ctbc.sd == tbc.sd` and avoids duplicate work for non-potcorr subjects.
    - `sd.orig` / `sd.orig_cdc` / `sd.orig_who` — not used in Steps 15/16. (Used only in the wrapper's preprocessing and inside `calc_and_recenter_z_scores()`.)

3. **Dead code.** No live dead code in Steps 15/16 after edits. Minor observations kept as Deferred below (not dead, but low-value):
    - `df[, exp_vals := cache$exp_vals]` at `:3801` and `:4009` (cache-hit branch) writes `exp_vals` to the local `df` copy; the value is not read downstream in the closure because `ewma_cache_init()` is only called in the cache-miss branch. Kept for symmetry with the cache-miss branch, which does need `exp_vals` on `df` for its `ewma_cache_init()` call. Low-value cleanup; deferred.
    - `sex` and `v` in Step 15 closure's `.SDcols` (`:3926–3927`) and Step 16 closure's `.SDcols` (`:4070–4071`) — neither is read inside the closures. Analogous to Session 7 Step 11 `.SDcols` trim (9→6 cols) but low-value. Deferred.

4. **Exclusion code names (both).** Only one exclusion code is assigned in Steps 15/16: `"Exclude-C-Traj"`, via `.child_exc(param, "Traj")` at 10 sites in the Step 15 closure (`:3862, 3864, 3868, 3870, 3872, 3874, 3878, 3880, 3883, 3885, 3892, 3894, 3896, 3898, 3900, 3902, 3904, 3906`) and 4 sites in the Step 16 closure (`:4050, 4052, 4054, 4056`). `"Exclude-C-Traj"` is in `exclude.levels.peds` at `:540`. All consistent with the canonical list in `CLAUDE.md`.

5. **Column names (child).** All current-state: `tbc.sd`, `ctbc.sd`, `tbc.p_plus`, `tbc.p_minus`, `ewma.all`, `ewma.before`, `ewma.after`, `c.ewma.all`, `dewma.all`, `dewma.before`, `dewma.after`, `c.dewma.all`, `tbc_diff_next` / `tbc_diff_prior` / `tbc_diff_plus_next` / `tbc_diff_plus_prior` / `tbc_diff_minus_next` / `tbc_diff_minus_prior`, `addcrithigh`, `addcritlow`, `tbc_dop`, `pot_excl`, `abssum`, `first_meas`, `sp_key`, `maxdiff`, `ageyears`, `exp_vals`, `index`, `internal_id`, `subjid`, `param`, `agedays`, `sex`, `v`, `exclude`. No Stata-era or legacy names.

6. **`.child_valid()` call correctness (child).** Call sites in Steps 15/16 (post-edit):
    - `:3677` — default (temp SDEs excluded). Used for the `sp_to_process_15` count. Fixed by F74 — Carrie confirmed temp SDEs should never be counted here.
    - `:3724, :3761` — `include.temporary.extraneous = FALSE`. Used for Step 15 filter. Correct.
    - `:3948` — default (after post-fix). Used for `sp_with_birth` selection. Correct.
    - `:3953` — default (after post-fix). Used for initial `step16_filter`. Correct.
    - `:3986` — default (after post-fix). Used for in-loop `step16_filter`. Correct.
    - All `include.extraneous` / `include.carryforward` flags are at their default FALSE — appropriate for EWMA2 (extraneous-SDE losers and CF non-rescues are permanent exclusions and must not participate).

7. **`by`-group correctness.** All groupings verified:
    - `by = sp_key` at `:3678, :3726, :3737, :3956, :3966` (aggregates for count / tbc range pre-filters).
    - `by = .(subjid, param)` at `:3714, :3718, :3770, :3774, :3925, :4069` (first_meas recalc per trajectory; Step 15 and Step 16 closures dispatched per-trajectory).
    - `by = .(subjid)` — not used in Steps 15/16 (no subject-level aggregation needed; the sp_with_birth check is vectorized via `%in%`).
    - `setkey(dop_snap, subjid, param)` at `:3756` — keys for the `dop_snap[.(subjid, dop_param)]` lookup at `:3845`; the `on = "agedays"` join at `:3848` is the per-row attribute. Correct.

8. **Sort order assumptions.** Step 15 closure relies on `agedays`-order within each (subjid, param) group for `diff(agedays)`, `c(tbc.sd[2:.N], NA)`, `c(NA, tbc.sd[1:(.N-1)])`, `c(agedays[2:.N], NA)`, etc. The pre-loop sort at `:3672` plus the data.table-by-dispatch preserves within-group sort. The pre-loop sort at Step 17 (`:4149`) re-establishes the invariant for downstream steps — Steps 15/16 themselves do not need a post-loop re-sort because no row is re-ordered (only `exclude` is overwritten via `:=`).

9. **data.table reference semantics.** `:=` targets are all `data.df` or the closure-local `df` (from `copy(.SD)`). No unintended parent-df mutations. `sde_results`-style subtle join pitfalls aren't relevant here (no temp-SDE rerun after Step 15/16).

10. **Factor level issues.** `Exclude-C-Traj` is in `exclude.levels.peds` at `:540`. No other factor levels are assigned by these steps. Safe.

11. **Parameter scope (child).** Step 15: all 3 params (WT including birth; HT/HC excluding birth). Step 16: HT and HC only (WT not applicable — WT birth is in Step 15). Narrative checklist item 6 matches; no change needed.

12. **Interaction with later steps.** Columns created by Steps 15/16 and dropped before Step 17: `sp_key`, `p_plus`, `p_minus`, `tbc.p_plus`, `tbc.p_minus`, `first_meas`. Drop is confirmed at `:4081` (sp_key) and `:4086–4088` (others, post-F81). Step 17 and later steps do not reference any of these. Working variables inside closures (`diff_before`, `diff_after`, `maxdiff`, `ageyears`, `exp_vals`, `ewma.all`, etc.) live on the local `copy(.SD)` and are discarded when the closure returns — not persisted to `data.df`. Only `exclude` survives back to `data.df`.

13. **DOP logic correctness (child).** Step 15's DOP lookup at `:3845` uses `get_dop(df$param[1])` to map the current param to its DOP. Per the canonical mapping (WT↔HT, HC→HT), `get_dop()`'s scalar-input contract (Session 9 F65) is satisfied because all rows in the per-group closure share the same `param`. Step 16 has no DOP logic (no `tbc_dop` or `dop_snap` reference).

20. **CRAN output compliance (wrapper/both).** All informational output in Steps 15/16 is wrapped in `if (!quietly) message(...)` (lines 3667–3668, 3740–3742, 3750–3751, 3943–3944, 3968–3970, 3978–3980). No `cat()` / `print()` calls.

---

## Deferreds

Items intentionally out of scope for this session, not opened as new deferreds unless noted:

- **Dead `df[, exp_vals := cache$exp_vals]`** at `:3801` (Step 15) and `:4009` (Step 16). Cache-hit branch writes `exp_vals` to local `df` but never reads it downstream in that branch. Kept for symmetry with the cache-miss branch, which does use `df$exp_vals` as input to `ewma_cache_init()`. Removing would be behavior-neutral but low-value; documenting here rather than in `CLAUDE.md → Known Issues` because it is strictly local to the closure and does not affect any observable output or performance in a measurable way.
- **Step 15 and Step 16 `.SDcols` contain unused `sex` and `v`** (`:3926–3927`, `:4070–4071`). Analogous to Session 7 Step 11 trim. Low-value; same note as above.
- **Redundant `order()` at `:3672`.** `data.df` was already keyed at the end of Main Step 13 Phase B4 (`setkey(data.df, subjid, param, agedays, internal_id)` at `:3658`). No mutation between `:3658` and `:3672` changes sort order, so the re-sort is defensive. Kept.
- **Structural `15a`/`15b`/`15c` + `16a`/`16b` subsection headers** for readability — would match the Session 7 Step 11 `11a`/`11b`/`11c` convention. Cosmetic only; deferred.
- **`child_clean.R:4148` `# Add id for consistent SDE order` comment** at the Step 17 pre-loop sort — same F68/F73 pattern, but belongs to the Step 17 walkthrough, not this one. Not opened as a new deferred item.
- **File-header `ewma()` / `ewma_cache_init()` / `ewma_cache_update()` bullet at `child_clean.R:94`.** Currently reads `ewma() / ewma_cache_init() / ewma_cache_update(): EWMA machinery`. An F66-style expansion would name Steps 11, 13, 17 for `ewma()` and Steps 15, 16 for `ewma_cache_init()` / `ewma_cache_update()`. Deferred — the current wording is general enough not to mislead a reader, and the Support Functions section of `child-algorithm-reference.md` already carries the detailed per-caller information.

---

## Discussion points resolved during walkthrough

- **Temp SDEs at Step 15 entry.** Traced every code path that assigns / unsets `Exclude-C-Temp-Same-Day`:
    - Main Child Step 13 Phase B4 (`child_clean.R:3650–3653`) overwrites every row whose `exclude` is currently Include or Temp-Same-Day with one of `{Include, Exclude-C-Identical, Exclude-C-Extraneous}`. `sde_results` is derived from `data.sde$exclude` after Phase B3, which guarantees every data.sde row ends up in that set (no Temp-Same-Day in data.sde post-B3).
    - No code path between the end of Main Step 13 and the start of Step 15 modifies `exclude`.
    - Conclusion: by Step 15 start, `exclude == "Exclude-C-Temp-Same-Day"` is empty. The `include.temporary.extraneous = TRUE` flag at `:3677` was pure defensive padding; per Carrie's preference, a straggler should not be counted as Include, so the flag was removed (F74).
- **Phase 15B closure's `exp_vals` assignment in the cache-hit branch.** Verified that `df[, exp_vals := cache$exp_vals]` at `:3801` (and `:4009`) is written to the local `df` copy but not read after — the only downstream `df$exp_vals` consumer (`ewma_cache_init()` at `:3814`) lives inside the cache-miss branch. Symmetry-only; left alone.
- **Step 16 subject-level birth filter.** Traced that the pre-fix filter at `:3948–3955` included HT or HC groups of any subject with *any* birth HT/HC, even if the specific param had no birth. Closure ran `ewma_cache_init()` on those groups but no candidate was eligible (all four Step 16 rules require `agedays == 0`), so the group produced no exclusion and dropped from `sp_to_process` on the next iteration. Switching to `sp_with_birth` (per-sp_key) is behavior-neutral and saves one `ewma_cache_init()` per such group.
- **Step 15 pre-filter scope.** `tbc_range_15` at `:3733–3739` skips groups whose tbc.sd range is ≤ 1 because `|dewma|` is bounded by the tbc.sd range, so `addcrit` (which requires `dewma.before > 1` and `dewma.after > 1`) can never fire. Same logic at `:3961–3967` for Step 16. Both pre-filters run once before the iteration loop, not re-applied each iteration — acceptable because any group that loses its range > 1 property through an exclusion will simply produce no new candidate on its next iteration and drop out via `sp_with_new_excl`.
- **`addcrit` interpretation for addcrithigh.** For positive outliers, the `tbc.p_plus` (perturbed upward) check is less strict than the raw tbc.sd check (the gap to the neighbor widens), and the `tbc.p_minus` (perturbed downward) check is stricter (the gap narrows). `addcrithigh` requires *all three* (raw, plus, minus) to exceed 1 — the narrow-gap (minus) leg is the binding constraint for borderline positive outliers. Symmetric for `addcritlow`.

---

## Post-fix test results

After reinstalling the package and running the full child test suite:

- test-cleangrowth.R: 63 PASS, 0 warnings
- test-child-regression.R: 48 PASS, 0 warnings
- test-child-edge-cases.R: 28 PASS, 0 warnings
- test-child-algorithms.R: 41 PASS (2 baseline codetools warnings, unchanged)
- test-child-parameters.R: 13 PASS

All child counts identical to baseline. No regressions. Adult tests not re-run — no adult files touched in this session.
