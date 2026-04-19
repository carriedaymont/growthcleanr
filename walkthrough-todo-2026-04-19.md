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

---

# R-vs-R comparison — Session 1 — 2026-04-19

**Scope:** Preprocessing — z-scores, GA correction (Child Step 2b), recentering.

Per `R-child-AprJanComparison-procedure.md`. Pure R-vs-R comparison; reference (Jan 2026) is `__Pipeline/current-infant-R-growthcleanr/growthcleanr/R/Infants_Main.R` (5011 lines); current (Apr 2026) is `__Pipeline/gc-github-latest/R/child_clean.R` (5027 lines).

This is a **comparison-only** session per procedure — **no code edits made**. All findings batched for Carrie's review and approval before any fix.

## Pre-session baseline tests (re-confirmed for this comparison)

- test-cleangrowth.R: 63 PASS, 0 warnings
- test-child-regression.R: 48 PASS, 0 warnings
- test-child-edge-cases.R: 28 PASS, 0 warnings
- test-child-algorithms.R: 41 PASS, 2 warnings (pre-existing `parallel=TRUE` / `plyr::ddply` codetools — baseline, unchanged)
- test-child-parameters.R: 13 PASS, 0 warnings

Adult tests not run — no adult files in scope.

## Scope

| Sub-area | Reference lines | Current lines | What |
|---|---|---|---|
| Param validation + LENGTHCM relabel | `Infants_Main.R` 389–399, 703 | `child_clean.R` 506–516, 702 | Validates param names; relabels LENGTHCM → HEIGHTCM |
| Imperial conversion | `Infants_Main.R` 558–561 | `child_clean.R` 669–672 | HEIGHTIN→cm, WEIGHTLBS→kg |
| Z-score calculation | `Infants_Main.R` 705–761 | `child_clean.R` 705–753 | CDC + WHO closures, blending, sd.orig_uncorr |
| GA correction Child Step 2b | `Infants_Main.R` 765–1180 | `child_clean.R` 755–1031 | Fenton + corrected WHO/CDC, sd.corr assignment, uncorr revert |
| Post-Step-2b cleanup, setkey, index, Missing/Not-Cleaned | `Infants_Main.R` 1195–1212 | `child_clean.R` 1029–1046, 1051 | Drop temp cols, re-key, mark Missing & Not-Cleaned |
| Recentering | `Infants_Main.R` 1217–1330 | `child_clean.R` 1056–1104 | sd_median or supplied recenter file, tbc.sd, ctbc.sd |
| Post-recentering safety + HC ≥ 5y | `Infants_Main.R` 1387, 1391 | `child_clean.R` 1115, 1121 | Re-mark NA tbc.sd as Missing, HC ≥ 5y handling |

**Function/name renames already documented in procedure** (treated as cosmetic): `id` → `internal_id`, `swap_parameters()` → `get_dop()`, `valid()` → `.child_valid()`, `cleanbatch_infants()` → `cleanchild()`, `temporary_extraneous_infants()` → `identify_temp_sde()`, `calc_oob_evil_twins()` → `calc_otl_evil_twins()`.

## Known intentional changes encountered (logged briefly per procedure, no further analysis)

These match items in the procedure's "Known intentional changes" section. Logged to confirm visited; no action.

- **`prelim_infants` flag wrapper removed.** Reference wraps the entire z-score / Step 2b block in `if (prelim_infants) { ... }`; current runs unconditionally. Tied to legacy parameter removal (CLAUDE.md 2026-04-16).
- **`cat()` → `message()`** throughout preprocessing, recentering, and progress logging.
- **`read_anthro()` signature simplified** — removed `prelim_infants = TRUE` argument.
- **`ref_tables` preload optimization** — current uses `if (!is.null(ref_tables)) ref_tables$mtz_cdc_prelim else read_anthro(...)`; reference reads from disk every call.
- **`round_stata()` removed at z-score / sd.corr / tbc.sd / ctbc.sd steps.** Per procedure "Don't flag rounding-tolerance removal."
- **Fenton z-score formula: LMS → CSD method.** Reference uses `((v_fenton/fen_m)^fen_l - 1)/(fen_l * fen_s)` with L/M/S; current uses `(v_fenton − M) / (S_upper * M)` and `(v_fenton − M) / (S_lower * M)` with M/S_upper/S_lower. Tied to Fenton 2025 swap (CLAUDE.md 2026-04-16).
- **Fenton data file renames + content swap:** `fentlms_foraga.csv.gz` → `fent_foraga.csv.gz`; `fentlms_forz` deleted; switched to `fenton2025_ms_lookup_smoothed.csv` (CSD method).
- **Reference WHO/CDC merge eliminated** (lines ~907–922 reference) — current uses the `measurement.to.z` / `measurement.to.z_who` closures directly without the `growthfile_who` / `growthfile_cdc_ext_infants` merges. Performance optimization; the closures contain the same reference data. Verified equivalent.
- **`### CP MOD ###` / `### CP REPLACE BLOCK ###` / `### CP WORK ###` commented-out scratch blocks removed** in current.
- **`pc <- copy(data.all[subjid %in% pc_ids])` subset refactor** — current works on a `pc` subset (potcorr subjects only) inside the `if (has_potcorr)` fast path; reference works on the entire `data.all`. Performance optimization.
- **Single `fcase()` for `sd.c` and final `sd.corr` assignment** — current collapses reference's multiple sequential `data.all[...condition..., sd.c := ...]` overrides into one `fcase()` per variable. Verified the resulting branch logic is equivalent except at one boundary (see AJ6 below).
- **`agedays / 365.25` → cached `ageyears_2b`** column inside Step 2b. Cosmetic.
- **`id` → `internal_id`** for sort keys, tiebreaks, and `id_sort` formula throughout (multiple sites).
- **`internal_id` creation timing** — reference assigns `internal_id := seq_len(.N)` BEFORE any sort (line 380), so internal_id reflects input row order; current sets the canonical key first (`setkey(data.all.ages, subjid, param, agedays, id)` at line 496) then assigns `internal_id := seq_len(.N)` (line 497), so internal_id reflects sorted order. Reference never uses internal_id for tiebreaking, so the timing difference is moot for reference; for current, sorted-order internal_id is the documented contract.
- **`Missing` → `Exclude-Missing`** and **`Not cleaned` → `Exclude-Not-Cleaned`** code renames.
- **`cf_rescued` column initialization** at line 1051 (current) — new column tied to CF rescue scheme.
- **Recentering simplified to single path** — current keeps only the precomputed-rcfile path; reference had three paths (rcfile, NHANES, derive). NHANES + derive paths removed per CLAUDE.md 2026-04-16 legacy removal.
- **`sdmedian.filename` and `sdrecentered.filename` parameters removed** — tied to recentering simplification.
- **`include.carryforward` → `cf_rescue` / `cf_detail`** parameter swap.
- **8 `biv.z.*` parameters added** — Step 7 BIV parameters (used downstream in cleanchild).
- **`weight_cap` → `adult_scale_max_lbs`** rename, `adult_columns_filename` removed, `adult_permissiveness` added.
- **`ewma_window`, `ref_tables`, `cached_results`, `changed_subjids`, `tri_exclude`, `batch_size` parameters added.**
- **Outer batching loop restructured** — `### BATCHING ###` block at current line 525 onward is a structural change. Not a preprocessing logic change; not deep-dived in this session.
- **`nnte` cleanup** — reference's trailing `if (prelim_infants) { data.all[, nnte := FALSE] }` block (lines 1393–1401) removed.

## Findings batched for review

Per procedure, findings are logged with category. **No code changes made.**

### AJ1. "HEIGHIN" typo in param validation (Bug fix)

- **Reference:** `Infants_Main.R:389–399`
- **Current:** `child_clean.R:506–516`
- **Issue:** Reference's accepted-param validation list contains `"HEIGHIN"` (missing `T`); user-supplied `"HEIGHTIN"` rows therefore fail the `param %in% c(...)` check and have `v` and `v_adult` set to `NA` before the imperial conversion at line 558 ever runs. The conversion line `data.all[param == "HEIGHTIN", v := v*2.54]` then operates on rows whose `v` is already `NA`, silently producing `NA * 2.54 = NA`. Net effect: HEIGHTIN data in reference is silently rendered Missing.
- **Current:** uses correct spelling `"HEIGHTIN"`; validation passes; imperial conversion works.
- **Category:** Bug fix.
- **Likely test impact:** Test data likely uses HEIGHTCM/WEIGHTKG only, so the bug wouldn't surface in regression. Worth confirming with a targeted check (e.g., add a HEIGHTIN-only smoke test).

### AJ2. NA-safe `!is.na(sd.orig)` guard added to potcorr_wt (Bug fix)

- **Reference:** `Infants_Main.R:801–803`
  ```r
  data.all[param == "WEIGHTKG",
           potcorr_wt := (seq_len(.N) == 1L &
                          janitor::round_half_up(janitor::round_half_up(sd.orig, 3), 2) < -2 &
                          agemonths < 10),
           by = subjid]
  ```
- **Current:** `child_clean.R:785–789`
  ```r
  data.all[param == "WEIGHTKG",
           potcorr_wt := (seq_len(.N) == 1L &
                          !is.na(sd.orig) & sd.orig < -2 &
                          agemonths < 10),
           by = subjid]
  ```
- **Issue:** Without the `!is.na(sd.orig)` guard, `NA < -2` returns `NA`, which propagates into `potcorr_wt` as `NA`. The downstream `potcorr := any(potcorr_wt, na.rm = TRUE)` (current line 795) and the failed-merge reset `potcorr_wt := FALSE` patches paper over this in many cases, but `seq_len(.N) == 1L` evaluating against a row whose `v` is `NA` (Missing) leaves potcorr_wt as `NA` for that row — could falsely flag a subject as potcorr (if the NA propagates through a misordered subgroup) or fail to flag a potcorr subject whose first WT row was Missing.
- The double-rounding removal (`janitor::round_half_up(...)` chain) is tied to rounding-tolerance removal and per procedure is not flagged.
- **Current:** Adds explicit `!is.na(sd.orig)` and rephrases as a current-state comment "NA-safe: sd.orig is NA when measurement is NA (Missing rows)".
- **Category:** Bug fix.

### AJ3. `intwt` low-weight floor bounds changed (Intentional, but tied to undocumented bound shift)

- **Reference:** `Infants_Main.R:821` — `data.all[intwt >= 250 & intwt <= 560, intwt := 570]`
- **Current:** `child_clean.R:822` — `pc[intwt >= 100 & intwt < 500, intwt := 500]`
- **Issue:** The floor target changed from 570 g to 500 g (matching the new Fenton table's minimum entry — INTENTIONAL, tied to the Fenton table swap). But the floor's applicable input range also changed: lower bound `250 → 100`, and the upper bound went from inclusive `<= 560` (above the new floor target) to exclusive `< 500` (at the new floor target).
  - At weight = 500 g exactly: reference floors to 570 (since 500 ≤ 560); current keeps at 500 (since 500 < 500 is FALSE). Behavior change at the boundary.
  - At weight in [100, 250) g: reference does not floor (intwt < 250 fails the lower bound); current does floor to 500. Both are below the old reference's 250 lower bound. Per current's comment "users who configure lower BIV thresholds can still get Fenton-corrected z-scores for very low weights" — this is the intent.
- **Category:** Intentional (other) — tied to Fenton table swap, documented in current's inline comment.
- **Note:** Behavior at 500 g exactly is the only practical edge case. Verifying that the new Fenton table's 500 g entry is the right anchor (vs. 570 g in old table) is prudent.

### AJ4. `tmp[subjid %in% names(table(subjid) > 1),]` filter no-op (Bug fix — semantically wrong, behavior-equivalent in practice)

- **Reference:** `Infants_Main.R:1122` — `tmp <- tmp[subjid %in% names(table(subjid) > 1),]`
- **Current:** `child_clean.R:984` — `tmp <- tmp[subjid %in% names(which(table(subjid) > 1)),]`
- **Issue:** `table(subjid) > 1` returns a named logical vector. `names()` on a logical vector returns ALL names (regardless of TRUE/FALSE). So reference's filter is a no-op (keeps all subjects). Intent was to filter to subjects with > 1 measurement. Current correctly uses `which()` to extract indices first.
- **Practical impact:** None observable. For single-row subjects, `abs(sum(sd.corr[1] - sd.corr))` collapses to `abs(0)` = 0, so the downstream `sd.corr_abssumdiff > sd.orig_abssumdiff` check (`0 > 0`) is FALSE. Single-row subjects naturally fail the replacement criterion. The buggy filter was redundant defense.
- **Category:** Bug fix (semantic correctness).

### AJ5. `c(orig_colnames, id)` ambiguous in column-keep filter (Bug fix — cleanup)

- **Reference:** `Infants_Main.R:1179–1180`
  ```r
  data.all <- data.all[, colnames(data.all) %in% c(orig_colnames, id),
                       with = FALSE]
  ```
- **Current:** `child_clean.R:1030–1031`
  ```r
  data.all <- data.all[, colnames(data.all) %in% c(orig_colnames, "id"),
                       with = FALSE]
  ```
- **Issue:** In reference, `id` is the function parameter (a vector of user-provided id values), NOT the string `"id"`. `c(orig_colnames, id)` therefore concatenates orig_colnames with all id values (coerced to character). The `%in%` check then keeps any column whose name happens to match an id value (an unlikely coincidence). In practice, `orig_colnames` already contains `"id"` (assigned at reference line 769 = `copy(colnames(data.all))` after `data.all.ages[, id := id]`), so the `id` column is preserved either way. Reference's code is benign but confusing.
- **Category:** Bug fix (cleanup; no observable behavior change).

### AJ6. sd.corr at exactly 4 years (1461 days) for potcorr subjects — FIXED

- **Reference:** `Infants_Main.R:1047–1064` (multi-statement sequential override)
  ```r
  # Step 4: smooth blend for 2-4y
  data.all[(agedays/365.25 >= 2) & (agedays/365.25 <= 4) &
           !is.na(sd.c) & !is.na(sd.orig),
           sd.corr := (sd.orig * (4 - agedays/365.25) + sd.c * ((agedays/365.25) - 2)) / 2]
  # Step 5: <2y potcorr override
  data.all[agedays/365.25 <= 2 & potcorr == TRUE & !is.na(sd.c),
           sd.corr := sd.c]
  # Step 6: >=4y or non-potcorr fallback
  data.all[agedays/365.25 >= 4 | potcorr == FALSE | is.na(sd.corr),
           sd.corr := sd.orig]
  ```
- **Pre-fix current:** `child_clean.R:933–939` had `ageyears_2b > 2 & ageyears_2b <= 4` on the smooth-blend branch.
- **Issue:** At exactly `agedays = 1461` (4 years exactly) for a potcorr subject:
  - **Reference:** Step 4 fires (`<= 4` TRUE). Smooth blend at age=4 evaluates to `(sd.orig * 0 + sd.c * 2) / 2 = sd.c`. Then Step 6 fires (`>= 4` TRUE). Final: `sd.corr = sd.orig`.
  - **Pre-fix current:** Second fcase branch fires (`> 2 & <= 4` TRUE). Smooth blend evaluates to `sd.c`. Default branch doesn't fire because second branch matched. Final: `sd.corr = sd.c`.
- **Resolution (Carrie-approved):** At exactly agedays=1461, sd.corr should be sd.orig (matching reference). Tightened the second `fcase()` branch upper-bound from `ageyears_2b <= 4` to `ageyears_2b < 4`, so the smooth-blend window is half-open `[2, 4)` and age=4 falls to the default `sd.orig` branch. Added a clarifying inline comment naming the half-open convention and noting that reference achieved the same result via a separate `>= 4 -> sd.orig` override after the smooth blend.
- **Post-fix current:** `child_clean.R:933–941` (3-line comment block + 6-line fcase, total 9 lines)
  ```r
  pc[, sd.corr := fcase(
    ageyears_2b <= 2 & potcorr & !is.na(sd.c),
      sd.c,
    ageyears_2b > 2 & ageyears_2b < 4 & !is.na(sd.c) & !is.na(sd.orig),
      (sd.orig * (4 - ageyears_2b) + sd.c * (ageyears_2b - 2)) / 2,
    default = sd.orig
  )]
  ```
- **Practical scope:** Only affects integer agedays = 1461 for potcorr subjects (a single integer day per such subject). At ages 1460 and 1462 both versions agree (1460 → smooth blend; 1462 → sd.orig).
- **Category:** Bug fix — boundary now matches reference January-2026 baseline.
- **Tests:** Reinstalled package and re-ran the full child suite — 63 / 48 / 28 / 41 / 13 pass, identical to baseline (no regressions). The 2 baseline `parallel=TRUE` / `plyr::ddply` codetools warnings on `test-child-algorithms.R` are unchanged. Adult tests not re-run — no adult files touched.

### AJ7. HC ≥ 5 years exclusion code changed: 'Missing' → 'Exclude-Not-Cleaned' (Intentional, but worth visibility)

- **Reference:** `Infants_Main.R:1391` — `data.all[param == "HEADCM" & agedays >= 5*365.25, exclude := 'Missing']`
- **Current:** `child_clean.R:1121` — `data.all[param == "HEADCM" & agedays >= 5*365.25, exclude := 'Exclude-Not-Cleaned']`
- **Issue:** Beyond the rename `Missing → Exclude-Missing` and `Not cleaned → Exclude-Not-Cleaned`, this is a semantic recategorization. Reference treated HC ≥ 5y as 'Missing' (no z-score data); current treats as 'Exclude-Not-Cleaned' (out of cleaning scope).
- The HC > 3*365.25 rule earlier at `child_clean.R:1046` already marks all HC > 3y as `Exclude-Not-Cleaned`. The line 1115 safety check `data.all[is.na(tbc.sd), exclude := 'Exclude-Missing']` then transitions HC ≥ 5y rows to `Exclude-Missing` (because their tbc.sd is NA from the missing WHO HC reference). Line 1121 then re-sets them back to `Exclude-Not-Cleaned`. Per current's comment: "this keeps a single consistent code across both 'we don't clean HC >3y' and '>=5y has no reference'."
- **Category:** Intentional (other) — code consistency, documented in current's inline comment. Not in procedure's "Known intentional changes" list, so flagged for visibility.

## Items NOT flagged (verified equivalent or out of scope)

For the record, these were checked and found to be either equivalent in behavior or already covered by the procedure's known-intentional-changes list:

- WHO/CDC corrected z-score block: reference computes via merge of `growthfile_who` + `growthfile_cdc_ext_infants` and adds 40+ unused columns; current uses `measurement.to.z` / `measurement.to.z_who` closures directly. Same z-score values; current is faster.
- `sd.c` assignment via `fcase()` (current 913–923) vs sequential `data.all[...]` overrides (reference 952–994). Logic-equivalent at all boundaries (≤ 2, 2 < x < 5, ≥ 5).
- `sd.c_temp` save/restore for Fenton override fallback (reference 1003–1014) vs current's "only overwrite when Fenton non-NA" (current 928–929). Logic-equivalent.
- `keep_id` SDE-Identicals filter for the uncorr-evaluation `tmp` subset: reference uses `id`, current uses `internal_id`. Within a `(subjid, agedays)` group, internal_id ordering matches id ordering (because internal_id was assigned in canonically-sorted order in current). So `min(id)` ↔ `min(internal_id)` and `max(id)` ↔ `max(internal_id)` give the same row. Equivalent.
- `examine_only` working subset (reference 1076 on `data.all`, current 942 on `pc`): `pc` is restricted to potcorr subjects, so `pc$param == "WEIGHTKG" & pc$potcorr` selects the same rows as `data.all$param == "WEIGHTKG" & data.all$potcorr`.
- `sub_replace` selection logic (`abs(sum(diff))` calculation, `is_first` filter, threshold `>` comparison): identical between reference and current.
- LENGTHCM → HEIGHTCM relabeling: identical (single-line `:=` in both).
- Recentering merge mechanics (`sd.recenter[data.all]`, `sd.median` subtraction): identical.

## Open questions for Carrie — RESOLVED

1. **AJ3 (intwt floor bounds):** Carrie confirmed the change is correct as-is. No documentation update opened in this session.
2. **AJ7 (HC ≥ 5y code):** Carrie confirmed the recategorization to `Exclude-Not-Cleaned` is correct as-is.

**All findings closed:**
- **AJ6** — Approved and fixed in this session (smooth-blend window changed to half-open `[2, 4)`).
- **AJ1, AJ2, AJ4, AJ5** — Carrie confirmed current code is correct; bug fixes already present, no further change needed.
- **AJ3, AJ7** — Carrie confirmed intentional; no change needed.

## Approval and next steps

Per procedure, no code changes have been made in this session. After Carrie's review:

- Approved fixes will be applied in a separate edit pass with a re-test (full 5-file child suite, expected unchanged 63/48/28/41/13).
- Any "Unclear" items resolved as "fix" or "no-op" will be applied accordingly.
- Findings index will be updated with status transitions (`open → approved → fixed → closed`).

Next session candidate per `R-child-AprJanComparison-procedure.md`: **Session 2 — Step 5 (Temp SDE) + Early Step 13 (SDE-Identicals)**.

---

# R-vs-R comparison — Session 2 — 2026-04-19

**Scope:** Early Step 13 SDE-Identicals block, Step 5 initial temp-SDE caller, and the shared helper `temporary_extraneous_infants()` → `identify_temp_sde()`. Helper internals walked in-line since Step 5's logic lives entirely there; later recalc callers (Steps 6, 7, 9, 11, and Main Step 13 final SDE) are out of scope for this session.

Per `R-child-AprJanComparison-procedure.md`. Pure R-vs-R comparison; reference (Jan 2026) is `__Pipeline/current-infant-R-growthcleanr/growthcleanr/R/Infants_Main.R` (5011 lines) plus `__Pipeline/current-infant-R-growthcleanr/growthcleanr/R/pediatric_support.R` (166 lines); current (Apr 2026) is `__Pipeline/gc-github-latest/R/child_clean.R` (5027 lines).

This is a **comparison-only** session per procedure — **no code edits made**. Both findings closed with no change needed.

## Pre-session git hygiene

Clean working tree; up to date with `origin/efficiency-updates`; no untracked work in stale `/Users/Shared/` copies.

## Pre-session baseline tests

Session 1 baseline retained (nothing touched since): **63 / 48 / 28 / 41 / 13**. Adult tests not re-run — no adult files in scope.

## Scope

| Sub-area | Reference lines | Current lines | What |
|---|---|---|---|
| Early Step 13 (SDE-Identicals) | `Infants_Main.R` 2625–2643 | `child_clean.R` 2628–2649 | Flag same-day same-value Include duplicates before CF detection |
| Step 5 (Temp SDE initial pass) | `Infants_Main.R` 2656–2717 | `child_clean.R` 2651–2712 | First call to temp-SDE helper; sets `Exclude-C-Temp-Same-Day` on SDE losers |
| Temp SDE helper (`identify_temp_sde()`) | `Infants_Main.R` 4888–5002 (`temporary_extraneous_infants()`) + `pediatric_support.R` 59–145 (`temporary_extraneous()`) | `child_clean.R` 2050–2200 | DOP-based SDE resolution with optional `exclude_from_dop_ids` |
| `.child_valid()` / `valid()` | `pediatric_support.R` 10–32 | `child_clean.R` 4990–5033 | Eligibility predicate used throughout |

**Function/name renames already documented in procedure** (treated as cosmetic): `valid()` → `.child_valid()`, `temporary_extraneous_infants()` / `temporary_extraneous()` → `identify_temp_sde()`, `id` → `internal_id`, `Exclude-SDE-Identical` → `Exclude-C-Identical`, `Exclude-Temporary-Extraneous-Same-Day` → `Exclude-C-Temp-Same-Day`. Step 5's `cleanbatch_infants()` caller → `cleanchild()`.

## Known intentional changes encountered (logged briefly per procedure, no further analysis)

- **dplyr → data.table refactor** of Early Step 13 block (`arrange` / `group_by` / `mutate` / `ungroup` / `select` → `order()` + `:=` with `by =`). Cosmetic idiom change; logic equivalent.
- **`Exclude-SDE-Identical` → `Exclude-C-Identical`** (via `.child_exc(param, "Identical")`) — code rename (covered by 2026-04-14/16 rename).
- **`Exclude-Temporary-Extraneous-Same-Day` → `Exclude-C-Temp-Same-Day`** — code rename.
- **User `id` → `internal_id`** for sort key, tiebreaker, and `setkey()` throughout Early Step 13 and `identify_temp_sde()`.
- **`valid()` → `.child_valid()`** rename.
- **`cat("...\n")` → `message("...")`** for both Early Step 13 and Step 5 progress messages.
- **Step 5 debug block stripped to a warning** — reference `fwrite()`-ed three debug CSVs and `stop()`-ped; current just `warning()`s when the temp-SDE filter would incorrectly mark all rows in a group. Production-readiness cleanup.
- **Helper column-subset optimization** — current passes a narrow 7-column subset to `identify_temp_sde()` so the internal `copy(df)` is cheap; reference `copy()`-ed the full `data.df`.
- **Defensive `is.null(df$subjid)` / `is.null(df$param)` blocks removed** from helper (Session 8 F51 / D-a cleanup — branches never fired; all callers pass both columns).
- **`keyby = .(subjid, param, agedays, id)` → `by = ...` + explicit `order(...)`** in helper narrowing — cosmetic.
- **`extraneous_result & valid.rows` redundant mask dropped** at helper return (Session 8 cleanup — `extraneous` already gated by `valid.rows & extraneous.this.day` in its `:=` assignment).
- **Named-vector mapping → direct integer indexing** at helper return (`setNames(..., df$orig_row)` + `as.character(original_rows)` lookup → `result <- logical(nrow(df)); result[df$orig_row] <- df$extraneous`) — Session 8 cleanup, behavior equivalent.
- **`janitor::round_half_up(..., 3)` then `round_half_up(..., 2)` double-rounding of `absdmedian.spz` / `absdmedian.dopz`** removed in current — per procedure, "Don't flag rounding-tolerance removal."
- **Commented-out mid-Step-5 `saveRDS()` debug block** removed.
- **`data.table(data.df)` re-conversion** after reference's dplyr chain — no longer needed in current.

## Findings batched for review

Per procedure, findings are logged with category. **No code changes made.**

### AJ8. Early Step 13 `keep_id` empty-group handling (Bug fix — minor robustness / warning avoidance)

- **Reference:** `Infants_Main.R:2634–2636`
  ```r
  keep_id = ifelse(agedays == 0,
                   min(id[exclude == "Include"], na.rm = TRUE),
                   max(id[exclude == "Include"], na.rm = TRUE))
  ```
- **Current:** `child_clean.R:2638–2643`
  ```r
  data.df[, keep_id := {
    incl_ids <- internal_id[exclude == "Include"]
    if (length(incl_ids) == 0L) NA_integer_
    else if (agedays[1L] == 0L) min(incl_ids)
    else max(incl_ids)
  }, by = .(subjid, param, agedays, v)]
  ```
- **Issue:** Current guards against all-excluded `(subjid, param, agedays, v)` groups with an explicit `length(incl_ids) == 0L` check. Reference's `ifelse(..., min(...), max(...))` with `na.rm = TRUE` has no such guard — when a group has zero Include rows (e.g., all-Missing or all-already-excluded same-day duplicates), `min(integer(0), na.rm = TRUE)` / `max(integer(0), na.rm = TRUE)` returns `Inf` / `-Inf` along with an R warning: `no non-missing arguments to min/max; returning Inf`.
- **Practical impact:** Behavior downstream is equivalent. In an all-excluded group, `has_dup` is FALSE (because `n_same_value = sum(exclude == "Include") = 0`), so the filter `has_dup & exclude == "Include" & id != keep_id` cannot fire regardless of what `keep_id` contains. Reference does not produce wrong results — it just leaks spurious warnings to users whose datasets contain all-excluded same-day groups (most likely from all-Missing rows on a given day).
- **Category:** Bug fix (warning avoidance / robustness). Logic pitfall category: **NA / empty-set handling** (`min()` / `max()` on empty subset → `Inf` / `-Inf` with warning).
- **Status:** Closed — current already correct, no change needed.

### AJ9. Reference `valid()` regex does not catch `"Missing"` / `"Not cleaned"` (Bug fix — latent in reference; current corrects via documented code rename)

- **Reference:** `pediatric_support.R:21`
  ```r
  keep <- !grepl("^Exclude", exclude)
  ```
- **Current:** `child_clean.R:5020`
  ```r
  keep <- !grepl("^Exclude", exclude)
  ```
- **Issue:** The regex is the same between the two files, but the exclusion-code strings it matches against are different. Reference's preprocessing assigns bare `"Missing"` and `"Not cleaned"` codes (no `Exclude-` prefix) at:
  - `Infants_Main.R:1207` — `fifelse(is.na(v) | agedays < 0, 'Missing', 'Include')`
  - `Infants_Main.R:1212` — `data.all[param == "HEADCM" & agedays > (3*365.25), exclude := "Not cleaned"]`
  - `Infants_Main.R:1387` — `data.all[is.na(tbc.sd), exclude := 'Missing']`
  - `Infants_Main.R:1391` — `data.all[param == "HEADCM" & agedays >= 5*365.25, exclude := 'Missing']` (current's AJ7 recategorized this to `Exclude-Not-Cleaned`)

  Neither `"Missing"` nor `"Not cleaned"` starts with `Exclude`, so reference's `!grepl("^Exclude", exclude)` returns TRUE for both — silently including these ineligible rows in `valid.rows` at every Step that calls `valid()`.
- **Where it matters in Step 5:** At Step 5, the in-helper flag `extraneous.this.day := (.N > 1), by = .(subjid, param, agedays)` runs inside `valid.rows`. Reference counts Missing/Not-cleaned rows as valid, so a subject with (say) two HC measurements in the 3–5y window (both `"Not cleaned"`, both carrying real `tbc.sd` from the recentering pass before the 1387 safety mark overwrites NA tbc.sd) could have one row silently re-labeled `"Exclude-Temporary-Extraneous-Same-Day"`. The Missing/Not-cleaned row (with NA `tbc.sd` → NA `absdmedian.spz`) sorts last by `order(absdmedian.spz)`, so the SDE loser assignment may reassign codes in a way that depends on upstream row order.
- **Current:** Uses the same regex but pairs it with renamed codes (`Exclude-Missing`, `Exclude-Not-Cleaned`), both of which `^Exclude` correctly catches. `.child_valid()` returns FALSE for these rows as intended.
- **Category:** Bug fix in reference; current corrects via a side-effect of the 2026-04-14 / 2026-04-16 code rename. Logic pitfall categories: **Factor levels / exclusion codes** + **`.child_valid()` flag scope**.
- **Status:** Closed — current already correct; the rename's bug-fix side effect is now on record in this findings log. Worth noting even though no code change is needed, because the rename was framed as cosmetic in the procedure's "Known intentional changes" section and the bug-fix dimension was not previously explicit.

## Items NOT flagged (audit trail)

For the record, reviewed and explicitly chosen not to open as findings — all either listed above under "Known intentional changes" or verified logic-equivalent:

- Rename of `temporary_extraneous_infants()` → `identify_temp_sde()` and consolidation with `temporary_extraneous()` from `pediatric_support.R` — one helper in current vs. two near-duplicates in reference. Walked in detail in Session 8; no behavioral divergence.
- The `valid.rows & extraneous.this.day` gating inside the helper — same logic in reference and current.
- `keyby = .(subjid, param, agedays, id)` in helper subset → `by = ...` + explicit `order(...)` in current. Cosmetic.
- `setNames(...)` named-vector lookup vs. direct integer indexing at helper return. Equivalent.
- Step 5's progress-message style (`cat` → `message`, loss of explicit newline).
- Pre-session inventory confirms only Step 5's first call to the helper is "initial" (no `exclude_from_dop_ids`); Main Step 13 is the sole non-NULL caller. Later recalc callers (Steps 6, 7, 9, 11 mid-loop, 11 end-of-step) reuse the same signature but go out of scope for this session.
- `data.table(data.df)` re-conversion after reference's dplyr chain — current's pure-data.table refactor eliminates the need.

## Open questions for Carrie — RESOLVED

Carrie approved (in the turn that asked "OK to batch AJ8 + AJ9 as closed, write up the walkthrough note + update findings index + CLAUDE.mds, and commit/push?"). Both findings closed with no code change:

- **AJ8** — current code is correct as-is; reference emits spurious warnings but behavior downstream is equivalent.
- **AJ9** — current code is correct as-is; reference bug masked by the 2026-04-14/16 exclusion-code rename. Documenting here so the rename's bug-fix side effect is explicit.

## Approval and next steps

Per procedure, no code changes have been made in this session (both findings closed without change). No new test run needed — baseline unchanged at 63 / 48 / 28 / 41 / 13.

Next session candidate per `R-child-AprJanComparison-procedure.md`: **Session 3 — Child Step 6 (CF)**. Major intentional change (new CF rescue scheme with `cf_rescue` parameter + lookup thresholds); expect a narrower comparison since most differences are intentional.

---

# R-vs-R comparison — Session 3 — 2026-04-19

**Scope:** Main Child Step 6 (CF detection + rescue). Per procedure, this step has the largest intentional-change surface of the comparison series (CF rescue scheme redesign) — expected to yield fewer findings than Sessions 1/2, but with higher importance since CF is a well-trafficked step.

Per `R-child-AprJanComparison-procedure.md`. Pure R-vs-R comparison; reference (Jan 2026) is `__Pipeline/current-infant-R-growthcleanr/growthcleanr/R/Infants_Main.R` plus `pediatric_support.R`; current (Apr 2026) is `__Pipeline/gc-github-latest/R/child_clean.R`.

This is a **comparison-only** session per procedure — **no code edits made.** **No findings opened this session.** All non-intentional diffs are housekeeping / cleanup with logic-equivalent behavior in both reference and current; see "Items NOT flagged" below for the audit trail.

## Pre-session git hygiene

Clean working tree; up to date with `origin/efficiency-updates`; no untracked work in stale `/Users/Shared/` copies.

## Pre-session baseline tests

Baseline retained from prior sessions (nothing touched since): **63 / 48 / 28 / 41 / 13**. Adult tests not re-run — no adult files in scope.

## Scope

| Sub-area | Reference lines | Current lines | What |
|---|---|---|---|
| Step 6 outer frame + entry guard | `Infants_Main.R` 2677–2700 | `child_clean.R` 2673–2690 | `sde_identical_rows` init; `include.carryforward` wrapper (removed in current); `not_single` and `valid_set` construction |
| CF pre-filter by sp_key | `Infants_Main.R` 2702–2715 | `child_clean.R` 2692–2705 | `uniqueN(v.orig) < .N` SP-level dup filter |
| CF detection (prior-day single-value match) | `Infants_Main.R` 2717–2755 | `child_clean.R` 2707–2743 | Six-step data.table join detecting current-matches-prior-single pattern |
| CF flag writeback | `Infants_Main.R` 2760–2815 | `child_clean.R` 2747–2758 | Assign `Exclude-Carried-Forward` / `Exclude-C-CF`; `any_cf` short-circuit |
| Temp-SDE re-run post-CF | `Infants_Main.R` 2822–2836 | `child_clean.R` 2763–2773 | Reset temp SDEs, re-run helper now that CFs are excluded |
| `wholehalfimp` computation | `Infants_Main.R` 2838–2876 | `child_clean.R` 2775–2784 | Row-level whole / half imperial unit flag |
| SDE-Identical save & remove | `Infants_Main.R` 2879–2885 | `child_clean.R` 2786–2790 | Temporarily pull SDE-Identicals before string detection |
| `ageday_has_include` + positional string detection | `Infants_Main.R` 2899–2931 | `child_clean.R` 2795–2826 | Originator, `cf_string_num`, `originator_z` assignment |
| String-number / z propagation loop | `Infants_Main.R` 2932–2957 | `child_clean.R` 2827–2852 | `max_iterations` loop filling `cf_string_num` / `originator_z` across consecutive CFs |
| `seq_win` / `cs` / `absdiff` | `Infants_Main.R` 2959–2977 | `child_clean.R` 2854–2866 | Position-in-string + z-diff-from-originator |
| Column cleanup (first wave) | `Infants_Main.R` 2980–2981 | `child_clean.R` 2868–2869 | Drop intermediate string-detection columns |
| CF rescue (big intentional change) | `Infants_Main.R` 2990–3031 | `child_clean.R` 2871–2948 | Reference: fixed 0.05 / 0.1 thresholds, 4 rescue codes (still exclusions) — current: `cf_rescue = "none"/"standard"/"all"` + lookup table, rescued CFs → `Include` with `cf_rescued` tag |
| SDE-Identical restore + final sort | `Infants_Main.R` 3036–3040 | `child_clean.R` 2952–2956 | `rbind` + `order(subjid, param, agedays, id / internal_id)` |
| `cf_detail` columns (current only) | — | `child_clean.R` 2958–2979 | Optional `cf_status` / `cf_deltaZ` output columns |
| Post-Step-6 column cleanup (current only) | — | `child_clean.R` 2981–2989 | Drop `v.orig, wholehalfimp, seq_win, cs, absdiff, ageday_has_include, orig_ageday, cf_interval` |

**Out of scope per procedure:** `.cf_rescue_lookup()` implementation (`child_clean.R:2406` and `cf-rescue-thresholds.md`) — reference has nothing to compare against.

## Known intentional changes encountered (logged briefly per procedure, no further analysis)

- **CF rescue scheme redesign** (2026-04-14; documented in procedure's intentional-changes list and `cf-rescue-thresholds.md`):
  - New `cf_rescue = c("standard", "none", "all")` parameter replaces deprecated `include.carryforward`.
  - Rescued CFs become `"Include"` and are tagged via new `cf_rescued` column (`"Rescued"` for standard-mode lookup-based, `"Rescued-All"` for `cf_rescue = "all"` mode).
  - Reference's four rescue codes (`Exclude-1-CF-deltaZ-<0.05`, `Exclude-1-CF-deltaZ-<0.1-wholehalfimp`, `Exclude-Teen-2-plus-CF-deltaZ-<0.05`, `Exclude-Teen-2-plus-CF-deltaZ-<0.1-wholehalfimp`) are gone — rescued CFs are now Include; unrescued stay as `Exclude-C-CF`.
  - Age / interval / param / rounding-specific lookup thresholds (0.05 / 0.20 / 0.40 / NR) replace the fixed 0.05 / 0.1 thresholds. Reference's age/sex gating (teens only for multi-CF rescue) is gone; lookup table handles age sensitivity via its age-bin axis.
  - Optional `cf_detail = TRUE` adds `cf_status` (`"CF-NR"` / `"CF-Resc"` / `NA`) and `cf_deltaZ` columns to output.
  - `cf_rescue = "all"` rescues every detected CF including shared-SPA CFs (CF on a day that also has an Include); downstream Main Child Step 13 resolves any resulting multi-Include SPAs (Session 4a change).
- **`include.carryforward` deprecated** — no longer wraps Step 6. Reference's `if (!include.carryforward) { ... }` block is gone; Step 6 always runs. `include.carryforward` remains as an argument of `.child_valid()` only, not exposed through `cleangrowth()`.
- **Exclusion code rename**:
  - `Exclude-Carried-Forward` → `Exclude-C-CF` (via `.child_exc(param, "CF")`).
  - `Exclude-Temporary-Extraneous-Same-Day` → `Exclude-C-Temp-Same-Day`.
  - `Exclude-SDE-Identical` → `Exclude-C-Identical`.
- **User `id` → `internal_id`** in all sort keys (`order(subjid, param, agedays, internal_id)` at CF detection sort `child_clean.R:2711` and at SDE-Identical restore `child_clean.R:2955`).
- **`valid()` → `.child_valid()`** rename at the top of Step 6.
- **`temporary_extraneous_infants()` → `identify_temp_sde()`** at the temp-SDE re-run call. Helper walked in detail in Session 2 + Session 8; not re-walked here.
- **`cat(...)` → `message(...)`** for progress messages (Step 6 progress lines, CF-pre-filter count, CF-rescue-pre-filter log).
- **Reference's `not_single` via `table(paste0(subjid, "_", param))` + `%in% names(...)`** (`Infants_Main.R:2694–2695`) → **current's data.table by-group `.N > 1` + right-join** (`child_clean.R:2686–2687`). Behavior equivalent; current is cleaner and avoids the theoretical underscore-concat key collision. Session 1's AJ4 is isomorphic (`names(table(subjid) > 1)` no-op filter).
- **`janitor::round_half_up(abs(...), 3)` double-rounding of `absdiff`** removed in current (`absdiff` is now just `abs(sd.orig_uncorr - originator_z)`). Per procedure, "Don't flag rounding-tolerance removal."
- **Dead commented-out code** — reference `Infants_Main.R:2762–2808` has two large blocks of commented-out prior CF-detection implementations (dplyr / `map_lgl`, legacy data.table); removed in current.
- **`### CP DROP IN` / `### CP UP` Stata-era marker comments** and commented-out alternate wholehalfimp formulations (`Infants_Main.R:2838–2850`) — removed in current.
- **Stale-changelog comments** — reference's `CP TOGGLED TO FALSE but then back to true because nothing changed?`, `Removed nnte filter`, `Fix:`, `Fix CF rescue logic` blocks — removed in current (F58 / F75 pattern from walkthrough series).
- **`#### CP Untested, replicate of above code in data.table()`** orphan commented-out stub — removed in current.

## Findings batched for review

**No findings opened this session.** All non-intentional diffs are cleanup / housekeeping with equivalent behavior in both reference and current; see "Items NOT flagged" below for the audit trail.

## Items NOT flagged (audit trail)

For the record, reviewed and explicitly chosen not to open as findings — each is equivalent behavior or already covered by the intentional-changes list. No AJ## number assigned.

1. **Reference `cf_string_length` dead variable** (`Infants_Main.R:2960` creates the column `sum(cf_binary == TRUE)` by (subjid, param, cf_string_num); `:2981` drops it). Grep confirms no read sites in between — reference's rescue logic keys on `max(seq_win)` to distinguish single-CF vs multi-CF strings, not on `cf_string_length`. Current has removed the computation entirely. Equivalent behavior. Category: cosmetic cleanup.

2. **Reference `is_eligible_include` no-op drop** (`Infants_Main.R:2980–2981` drops this column in cleanup, but it is never created anywhere in Step 6 or earlier in the algorithm — data.table `:= NULL` on a non-existent column is silent). Current's cleanup list omits it. Equivalent behavior. Category: cosmetic cleanup.

3. **`wholehalfimp` for `HEIGHTCM` / `HEADCM` uses redundant `%% 1 | %% 0.5` OR in reference** (`Infants_Main.R:2869–2870` and `:2874–2875`); current uses only `%% 0.5 < 0.01` (`child_clean.R:2782, :2784`). The `%% 1 < 0.01` branch is redundant — for a whole integer N, both `N %% 1` and `N %% 0.5` equal 0, and for a half-integer N+0.5, `N+0.5 %% 0.5` equals 0 while `N+0.5 %% 1` equals 0.5. The `%% 0.5` check alone catches both whole and half cases. `WEIGHTKG` line uses only `%% 1 < 0.01` in both impls (reference `:2865`, current `:2780`) — whole-pound check only, correct in both. Equivalent behavior. Category: cosmetic cleanup.

4. **Step 6 working columns persist past Step 6 in reference** (`v.orig`, `wholehalfimp`, `seq_win`, `cs`, `absdiff`, `ageday_has_include`). Grep confirms no read sites past Step 6 in reference — all later uses of similar names are either Step 2b's own local `tmp[, seq_win := sequence(.N), by = subjid]` (on a different data.table, scoped out) or Step 13's own `absdiff_rel_to_median` / `absdiff_dop_med` (different variable names). Current drops all Step 6 working columns explicitly (`child_clean.R:2984–2989`). Orphan columns in reference are benign. Equivalent behavior. Category: output tidiness.

5. **Reference `not_single` paste-concat construction** — `paste0(subjid, "_", param) %in% names(tmp)[tmp > 1]` has theoretical underscore-concat key-collision risk, but the `WEIGHTKG` / `HEIGHTCM` / `HEADCM` param vocabulary is fixed and contains no underscores, so collision is impossible in practice. Current uses clean `by = .(subjid, param)` group-aggregation instead. Isomorphic to Session 1's AJ4. Not re-flagged here — kept as an audit-trail note that the same latent pattern was reviewed.

6. **`cf_subset <- cf_subset[order(subjid, param, agedays, internal_id)]` (current)** creates a re-ordered copy; reference's `setorder(cf_subset, subjid, param, agedays, id)` orders in place. Same final ordering; trivial efficiency nit only. (The `id` → `internal_id` rename is the only meaningful difference and is documented in the intentional-changes list.)

7. **Current's `data.df[!exclude == "Exclude-C-Identical"]`** parses as `!(exclude == "Exclude-C-Identical")`, i.e., equivalent to `exclude != "Exclude-C-Identical"` — matches reference's `data.df[exclude != "Exclude-SDE-Identical"]` (after rename). Same filter.

8. **Temp-SDE re-run helper call** — `temporary_extraneous_infants(data.df)` (ref) → `identify_temp_sde(data.df[, .(id, internal_id, subjid, param, agedays, tbc.sd, exclude)])` (current). Current passes a narrower column subset for efficient `copy()` inside the helper. Helper walked in detail in Sessions 2 + 8; no behavioral divergence in the re-run call path.

9. **Sort direction for SDE-Identical restore** — reference `order(subjid, param, agedays, id)` vs current `order(subjid, param, agedays, internal_id)`. `id` → `internal_id` rename only (documented intentional change). Same sort semantics.

10. **Originator selection**: both impls use `exclude == "Include" & nextcf == TRUE, originator := TRUE`. The comment "Any Include can be an originator - ageday_has_include only restricts CF rescue eligibility" is present in both. Identical logic.

11. **String propagation loop** — both use `for (i in 1:max_iterations)` with `shift(cf_string_num, type = "lag")` under the same `(ageday_has_include == FALSE | is.na(ageday_has_include))` guard. Identical structure and conditions.

12. **CF rescue outer filter** — both gate on `!is.na(seq_win) & (ageday_has_include == FALSE | is.na(ageday_has_include))`. Current adds an explicit `seq_win > 0` condition at the outer filter; reference folds this into the inner rescue-code branches via `seq_win != 0`. Equivalent effective filter — the outer gate plus inner filter in reference covers the same rows as current's combined outer gate.

## Open questions for Carrie — RESOLVED

No open questions. 0 findings to approve; the housekeeping items above are informational only.

Carrie approved (in the turn that asked "OK to proceed with writing up and committing?") the session conclusion: 0 AJ## findings, session closes with walkthrough note + CLAUDE.md updates + single commit.

## Approval and next steps

Per procedure, no code changes have been made in this session. No new test run needed — baseline unchanged at 63 / 48 / 28 / 41 / 13. Adult tests not re-run — no adult files in scope.

Next session candidate per `R-child-AprJanComparison-procedure.md`: **Session 4 — Child Step 7 (BIV) + Child Step 9 (Evil Twins)**. Step 7 has known intentional change (eight `biv.z.*` per-cell parameters replacing single `sd.extreme` / `z.extreme`); Step 9 expected low–medium complexity.
