# Adult GC Walk-Through — Deferred To-Dos

Items identified during the 2026-04-09 walk-through that are deferred for later resolution.
Focus: post-wtallow-redesign reconciliation of code, narrative, wtallow-formulas.md, and permissiveness spec.

---

## Pre-walk-through setup

(populated after tests pass)

---

## Steps walk-through (populated as walk-through proceeds)

### Step 9Wb (Extreme EWMA) — `remove_ewma_wt` lines 981–982

**Issue:** Comment says "NA → Inf for edge values (Stata missing-as-infinity)" but `adult_ewma_cache_init()` never produces NA for `ewma_before`/`ewma_after` — edge values fall back to `ewma_all`. The `ifelse(is.na(...), Inf, ...)` is a defensive no-op. If it ever did trigger, `Inf` would correctly confirm positive exclusion but would **not** confirm negative exclusion (would need `-Inf` for the negative case, matching Stata `missing < x` = TRUE).

**Why deferred:** No behavioral impact — NAs never occur, and the fallback-to-ewma_all behavior is correct in both directions. Misleading comment only.

**Fix needed:** (1) Update comment to explain that edge values fall back to `ewma_all` via the cache, making the Inf conversion unnecessary. (2) Optionally, if keeping the defensive code, split into `Inf` for positive and `-Inf` for negative to be correct if NAs ever did occur.

**Fix applied:** Removed `dewma_bef_safe`/`dewma_aft_safe` ifelse wrappers and updated comment to explain ewma_all fallback. Criteria now use `dewma_bef`/`dewma_aft` directly. **RESOLVED.**

### Step 11Wb (Moderate EWMA) — narrative trigger guard (line 979)

**Issue:** Narrative says "Skip if weight range ≤ min wtallow AND min percent ratio ≥ 0.7" but the code uses `max_perc` (the actual perclimit from the permissiveness level, which varies by weight tier and level), not a fixed 0.7. The "0.7" was approximately correct for tighter/tightest `perclimit_low` but is imprecise for other levels/tiers.

**Why deferred:** No behavioral impact — code is correct, narrative is imprecise. Not misleading enough to cause future work errors.

**Fix needed:** Update narrative line 979 to say "Skip if weight range ≤ min wtallow AND (perclimit disabled OR min percent ratio ≥ max perclimit)" to match the code's actual logic.

**Fix applied:** Narrative updated to match code logic. **RESOLVED.**
