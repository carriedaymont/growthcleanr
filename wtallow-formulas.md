# wtallow base formulas and adjustment

## General concepts and terms

***wtallow*** wtallow is a derived value used in multiple steps. Conceptually, it represents the threshold above which an absolute difference in weight between two values will result in exclusion of one of the two values. It is one of the most important values determining whether a value should be excluded. Importantly, most steps where wtallow is used include other constraints on exclusion or other paths to exclusion. Therefore, not every difference higher than wtallow will result in exclusion, and some differences below wtallow will result in exclusion. wtallow is dependent on two things - the interval between the two measurements being compared and the upper weight (UW) of the two measurements. 

***Upper weight (UW)*** UW is the higher of the two adjacent weights that are being evaluated when calculating wtallow. An UW of 120 is considered the base value. UW >120 are considered high and UW <120 are considered low.

***wtallow-et-cap*** The wtallow-et-cap is a value related to wtallow used only in the Evil Twins step. In the Evil Twins step, any time adjacent weights have an absolute difference larger than wtallow-et-cap, one of the two values will get excluded. The wtallow-et-cap-6m is the et-cap for intervals <=6m. The wtallow-et-cap-12m is the et-cap for intervals >6m.

***wtallow-cap*** These are values at 6m and 12m that represent a maximum value for wtallow even if the formula for wtallow would otherwise give a higher value. For intervals <=6m, no value can be above the wtallow-cap-6m. For intervals >6m, no value can be above the wtallow-cap-12m (this includes intervals >12m). The wtallow-cap values are also used in the wtallow formulas.

***Piecewise-Higher (PW-H)*** The formula used to determine wtallow in looser and loosest permissiveness levels. It is adjusted for high UW and low UW.
***Piecewise-Lower (PW-L)*** The formula used to determine wtallow in tighter permissiveness levels. It is adjusted for low UW (not high UW).
***allofus15*** The formula used to determine wtallow in tightest permissiveness levels. It is not adjusted based on UW.

***Suffixes*** For PW-H and PW-L, we will use suffixes -Base, -highUW, and -lowUW to designate formulas for the three UW conditions (=120, >120, and <120).

***wtallow-1-day*** The value of wtallow at one day. It is equal to `10 + 10 × log(1 + 5×1/30.4375) / log(6)`, approximately 10.85

***ceiling*** A maximum wtallow for PW-H-lowUW and PW-L-lowUW that limits wtallow beyond the limits provided by the wtallow-cap. A ceiling is only used for PW formulas adjusted for lowUW.

### Brief rationale for features of base formulas

Based on literature, input from experts, and clinical experience, the following concepts informed wtallow.

-Weight change can be large even in a single day or very few days from fluid shifts. After 1-2 weeks, the slope of potential weight change generally slows.
-Many of the largest weight changes occur after bariatric surgery or GLP1 agonist use. Literature shows that maximum weight loss after bariatric surgery or GLP1 agonist treatment occurs at approximately 12m, so wtallow does not rise after 12m. Of this weight change, a higher proportion occurs in the first 6m, so the slope is higher from 1m to 6m than 6m to 12m. Time courses for large weight loss with cancer can vary, but would generally not be larger or faster than those from weight loss surgery. 
-In general, PW-H was designed to include (i.e., not exclude) almost all physiologically plausible changes. PW-L was designed to include physiologically plausible changes for many, but not all, conditions. allofus15 was adapted from the Guide et al. algorithm to include changes that would be seen in people without severe or uncommon health problems.

### Rationale for adjustment

People with higher weight can have larger absolute changes in body weight. Similarly, people with lower weight will typically have smaller absolute changes in body weight. The relationship is complex, varies over time, and depends on whether you are considering weight gain or weight loss. Adjusting allowed wt change by the upper weight (UW) of a pair can account for some of these differences, although it does not account for all situations. Most important, adjusting the wtallow (max wt changed allowed) can increase the risk of inappropriately including erroneous high weights. For this reason, we cap the adjustment for upper weights at the values that would be used for an upper weight of 180.

## Formulas

### PW-H (Piecewise-Higher)

**PW-H-Base (UW = 120, no adjustment):**
wtallow-cap-H-6m-base = 50
wtallow-cap-H-12m-base = 80
wtallow-et-cap-H-6m-base = 70
wtallow-et-cap-H-12m-base = 100
ceiling = N/A

months ≤ 1: `10 + 10 × log(1 + 5×months) / log(6)`
1 < months ≤ 6: linear increase from 20 @1m to wtallow-cap-H-6m-base @6m; works out to `20 + 6 × (months − 1)`
6 < months ≤ 12: linear increase from wtallow-cap-H-6m-base @6m to wtallow-cap-H-12m-base @12m; works out to `50 + 5 × (months − 6)`
months > 12: flat at wtallow-cap-H-12m-base (80)

**PW-H-highUW (Scaling up, UW>120):**
***if UW is >180, use 180 as UW in formula***
wtallow-cap-H-6m-highUW = 50 + 0.25 x (UW-120)
wtallow-cap-H-12m-highUW = 80 + 0.25 x (UW-120)
wtallow-et-cap-H-6m-highUW = 70 + 0.25 x (UW-120)
wtallow-et-cap-H-12m-highUW = 100 + 0.25 x (UW-120)
ceiling = N/A

months ≤ 1: linear increase from wtallow-1day @ 1 day to 25 at 1m
1 < months ≤ 6: linear increase from 25 @1m to wtallow-cap-H-6m-highUW @6m
6 < months ≤ 12: linear increase from wtallow-cap-H-6m-highUW @6m to wtallow-cap-H-12m-highUW @12m
months > 12: flat at wtallow-cap-H-12m-highUW

**PW-H-lowUW (Scaling down, UW < 120):**
wtallow-cap-H-6m-lowUW = (50 − 20) × (UW / 120) + 20
wtallow-cap-H-12m-lowUW = (80 − 20) × (UW / 120) + 20
wtallow-et-cap-H-6m-lowUW = (50 − 20) × (UW / 120) + 40
wtallow-et-cap-H-12m-lowUW = (80 − 20) × (UW / 120) + 40
ceiling: UW × 2/3

months ≤ 1: `10 + 10 × log(1 + 5×months) / log(6)`
1 < months ≤ 6: linear increase from 20 @1m to wtallow-cap-H-6m-lowUW @6m
6 < months ≤ 12: linear increase from wtallow-cap-H-6m-lowUW @6m to wtallow-cap-H-12m-lowUW @12m
months > 12: flat at wtallow-cap-H-12m-lowUW

### PW-L (Piecewise-Lower)

**PW-L-Base (UW = 120, no adjustment):**
wtallow-cap-L-6m-base = 33.33
wtallow-cap-L-12m-base = 53.33
wtallow-et-cap-L-6m-base = 53.33
wtallow-et-cap-L-12m-base = 73.33

months ≤ 1: `10 + 10 × log(1 + 5×months) / log(6)`
1 < months ≤ 6: linear increase from 20 @1m to wtallow-cap-L-6m-base @6m
6 < months ≤ 12: linear increase from wtallow-cap-L-6m-base @6m to wtallow-cap-L-12m-base @12m
months > 12: flat at wtallow-cap-L-12m-base


**PW-L-highUW (UW>120):**           ***SAME AS BASE FORMULA***
wtallow-cap-L-6m-highUW = 33.33
wtallow-cap-L-12m-highUW = 53.33
wtallow-et-cap-L-6m-highUW = 53.33
wtallow-et-cap-L-12m-highUW = 73.33

months ≤ 1: `10 + 10 × log(1 + 5×months) / log(6)`
1 < months ≤ 6: linear increase from 20 @1m to wtallow-cap-L-6m-base @6m
6 < months ≤ 12: linear increase from wtallow-cap-L-6m-base @6m to wtallow-cap-L-12m-base @12m
months > 12: flat at wtallow-cap-L-12m-base

**PW-L-lowUW (Scaling down, UW < 120):**
wtallow-cap-L-6m-lowUW = (33.33 − 20) × (UW / 120) + 20
wtallow-cap-L-12m-lowUW = (53.33 − 20) × (UW / 120) + 20
wtallow-et-cap-L-6m-lowUW = (33.33 − 20) × (UW / 120) + 40
wtallow-et-cap-L-12m-lowUW = (53.33 − 20) × (UW / 120) + 40
ceiling: UW × 2/3

months ≤ 1: `10 + 10 × log(1 + 5×months) / log(6)`
1 < months ≤ 6: linear increase from 20 @1m to wtallow-cap-L-6m-lowUW @6m
6 < months ≤ 12: linear increase from wtallow-cap-L-6m-lowUW @6m to wtallow-cap-L-12m-lowUW @12m
months > 12: flat at wtallow-cap-L-12m-lowUW

---

### allofus15

allofus15 is not adjusted based on UW except that its 12m cap is limited to never exceed the PW-L 12m cap for that UW.

allofus15-cap-12m = min(40, effective PW-L at 12m for that UW)
- For UW ≥ 120: effective PW-L-12m = 53.33, so allofus15-cap-12m = 40
- For UW < 120: effective PW-L-12m = min(wtallow-cap-L-12m-lowUW, UW × 2/3) where wtallow-cap-L-12m-lowUW = (53.33 − 20) × (UW / 120) + 20; allofus15-cap-12m = min(40, that effective value)

allofus15 wtallow formula:
1–2 days: 5
3–7 days: 10
8d to <6m: 15
6m to 12m: linear increase from 15 @6m to allofus15-cap-12m @12m
>12m: flat at allofus15-cap-12m

allofus15-et-cap: allofus15-cap-12m
