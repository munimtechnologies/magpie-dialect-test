# READING-MAGPIE: do magpie fledglings use their own group's call sequences?

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-25. Made with AI agents working in MV2's name; no
human expert has checked it yet. Text CC BY 4.0. Every number below comes from
[results/magpie_results.json](results/magpie_results.json) and
[results/RESULTS-MAGPIE.md](results/RESULTS-MAGPIE.md), written unchanged by job magpie-dialect-v1
(run 2026-11-25 20:38-20:53 UTC, all unit tests passed, the plan and code were checked against
GitHub main before anything was computed). Plan: [PLAN.md](PLAN.md) (2026-11-23) and
[AMENDMENT-1.md](AMENDMENT-1.md) (2026-11-24), both public before the job ran.

Claim under test: Mason, King and Ridley (2026), Proc. R. Soc. B 293: 20251620,
doi:10.1098/rspb.2025.1620: Western Australian magpie fledglings produce the call sequences of
their own social group rather than those of the wider population, so the sequences are socially
learned. Data: the authors' figshare deposit doi:10.6084/m9.figshare.28251398 v2 (CC BY 4.0); all
five files matched figshare's MD5s.

## Short answer

1. **The group signal is real under a fair test.** A fledgling's repertoire overlaps its own
   group's record more than other groups' records: D = 0.1467 (mean own overlap 0.60 against 0.45
   for the other groups), p = 0.0124 against 9,999 random group assignments. This test fired in 1
   of 200 worlds in which each fledgling was paired with another group's records (0.5%), and it
   found a planted effect in 50 of 50 worlds. On the authors' robustness file (single occurrences
   removed) D = 0.1756, p = 0.0042.
2. **The published test cannot tell a fledgling's own group from any other group.** We reproduced
   the authors' model (binomial mixed model, fledgling presence of a type on group presence,
   AICc comparison): group coefficient 2.25, z 9.29, AICc 631.3 against 721.0 (population) and
   719.0 (intercept only). But the same decision, "group repertoire predicts", also fired in
   **200 of 200** donor worlds, where every fledgling was given the real records of a fledgling from
   a different group, a group it never lived in. So the published p value and AICc gap measure that
   common sequence types are common, not that a fledgling matches its own group. The model has no
   term for the sequence type.
3. **Order: no group-specific order was detected, and the test had little power.** Among types
   whose units also occur in another order, a fledgling used its group's order slightly more often
   (D_ord = 0.035, p = 0.161; 9 fledglings, 38 such types produced). Power in planted worlds was
   only 0.14, so this is absence of evidence. The authors' suggestion that "the way calls are
   combined" is what is learned is neither supported nor contradicted by these data.
4. **Social contacts do not explain repertoire size once recording effort is accounted for.**
   Repertoire size tracks the number of the fledgling's own recorded sequences (Spearman rho 0.94).
   The raw correlation of contacts with size is -0.004; with effort removed it is -0.88 (one-sided
   p for a positive effect 0.9988). With 9 fledglings this is descriptive.

## Registered verdicts

| test | statistic | p | reading |
|---|---|---|---|
| T0 reproduction | group b1 2.25, z 9.29; AICc 631.3 / 721.0 / 719.0 | 1.6e-20 | published-style decision fires |
| T1 own-group overlap excess (primary) | D 0.1467 (null mean 0.0421, 95th pct 0.1227) | 0.0124 | OWN_GROUP_EXCESS |
| T2 order beyond content | D_ord 0.0347 | 0.1612 | ORDER_ADDS_NOTHING_MEASURABLE |
| T3 donor worlds (200) | T0 decision fired 200/200; T1 fired 1/200 | | T0 MISCALIBRATED; T1 calibrated |
| Power (50 planted worlds each) | T1 50/50; T2 7/50 | | T1 1.00; T2 0.14 |
| T4 contacts vs size, effort removed | rho -0.879 (raw -0.004; effort-size 0.941) | 0.9988 | CONTACTS_EFFECT_NOT_SEPARABLE_FROM_EFFORT |

Secondary file (single occurrences removed): T0 fires (b1 2.43); T1 D 0.1756, p 0.0042; T3: T0
fired in 200 of 200 donor worlds, T1 in 0 of 200.

## Predictions (priors fixed in PLAN.md and AMENDMENT-1.md)

| id | prediction | prior | outcome | Brier |
|---|---|---|---|---|
| Q1 | T1 reads OWN_GROUP_EXCESS | 0.60 | true | 0.160 |
| Q2 | T2 reads ORDER_CARRIES_GROUP_SIGNAL | 0.25 | false | 0.063 |
| Q3 | the published-style decision is MISCALIBRATED (donor worlds) | 0.40 | true | 0.360 |
| Q4 | T4 reads CONTACTS_EFFECT_SURVIVES_EFFORT | 0.35 | false | 0.123 |
| Q5 | T2 is UNINFORMATIVE | 0.30 | false | 0.090 |
| Q6 | T0 decision fires in more than 50% of donor worlds | 0.55 | true | 0.203 |

Mean Brier score 0.166. T2 was not UNINFORMATIVE by the registered rule (38 multi-order types
produced, 9 fledglings entered), but its power of 0.14 means its negative reading carries little
weight.

## What it means

The headline, that fledglings match their own group, survives a calibrated test, with an effect
of about 15 percentage points of repertoire overlap. The evidence the paper published for it does
not measure that: the same model and decision rule declare "group repertoire predicts" for a group
the fledgling never belonged to, every time. For the ledger this is the pattern already seen in the
parrot duets: the effect survives, the published evidence does not.

## Limits

- 10 fledglings in 8 groups in the public file (the paper reports 11); one fledgling (NHFL) has an
  empty repertoire and drops out, leaving 9.
- Presence/absence only: there are no token counts or dates per type, so how often a type was used
  carries no weight.
- A fledgling's group record was captured during its own follows, so a shared event (the same
  intruder, the same alarm) puts a sequence in both. T1 cannot separate that from learning.
- Our model fit is our own (Gauss-Hermite maximum likelihood in numpy/scipy), not glmmTMB, so its
  coefficients may differ slightly from the paper's; the decision it feeds is the paper's.
- This was done with AI agents working in MV2's name and has not been checked by a human expert.
