# magpie-dialect-test: AMENDMENT-1 (2026-11-24, before any statistic was computed)

MV2 of Munim, Inc. Made with AI agents working in his name. Apache-2.0 (code), CC BY 4.0 (text).
This amendment changes the tests of [PLAN.md](PLAN.md) (v0, 2026-11-23) because the public files
do not contain what v0 assumed. It is written and published before the analysis job runs; the job
waits until this file and the code are public (gate.py).

## What was looked at before writing this

Exactly the following, recorded so readers can judge it (see LOADER-NOTES.md for the details):
the figshare API record (licence CC BY 4.0, version 2, file names, sizes, MD5s); the authors'
README.txt; the two R scripts (code only); the header line of each of the four CSV files; the first
two or three data rows of each CSV; the row counts of COMBINATION_repertoire.csv; the list of the
117 sequence-type labels; and, for each of the 10 fledglings, the metadata fields of the row for one
type label (ID, group, mother, closest contact, total number of contacts, proportion of time in
contact). No presence/absence column was tabulated, counted or compared for any fledgling, group or
type, and no model was fitted.

## What the files are

- COMBINATION_repertoire.csv: one row per fledgling and sequence type (10 fledglings x 117 types =
  1,170 rows). `call` is the ordered sequence of call units, separated by spaces (7 units: DS, LH,
  LHDS, NL, NLDS, SH, SHDS; the labels are not alphabetised, so order is recorded). `fl01` = the
  fledgling produced the type; `group01` = the type was recorded from its natal group (README:
  "group vocal data captured during fledgling follows"); `exgroup01` = the type was recorded in the
  broader population outside its group.
- There are **no token counts and no dates per type**: only presence/absence.
- 10 fledglings in 8 groups (JOG has 3 fledglings; RVD, GON, SCL, NH, PR, BWYb and MBG have 1 each).
  The paper says 11 fledglings; the public combination file has 10.
- The authors' test (Repertoire modelling.R): binomial GLMMs `fl01 ~ group01 + (1|ID)`,
  `fl01 ~ exgroup01 + (1|ID)` and `fl01 ~ 1 + (1|ID)` (glmmTMB), compared by AICc; the group model
  wins and its coefficient is reported. Rows are fledgling x type; **no term for the type**, so a type
  common everywhere is counted as evidence for group learning whenever it is in the group's record.

## Changes

1. **Units.** A fledgling's repertoire F_f is its set of types with fl01 = 1. Its follow record
   G_f is its set with group01 = 1; E_f its set with exgroup01 = 1. A group's record R_g is the union
   of G_f over its fledglings (for JOG, the union of 3 follow records).

2. **T0, reproduction (replaces v0 T0).** The three GLMMs above, fitted by maximum likelihood with a
   40-node Gauss-Hermite integral over the fledgling intercept (numpy/scipy, not glmmTMB). Reported:
   the group01 coefficient, its Wald z and p, and AICc of the three models. The **published-style
   decision** says "group repertoire predicts" when the group01 coefficient is > 0 with Wald p <= 0.05
   and the group model has the lowest AICc of the three.

3. **T3, calibration of the published-style decision (replaces v0 T3).** v0's pooled-token worlds
   cannot be built without tokens. Instead, 200 **donor worlds**: each fledgling f is given the real
   follow record and extra-group record (G_h, E_h) of a donor fledgling h from a **different** group
   (a random assignment in which no fledgling keeps a donor from its own group; the donor map is a
   permutation of the 10 fledglings). F_f is unchanged. In these worlds the "group" is a group the
   fledgling never lived in, so any "group repertoire predicts" decision is a false alarm. Report the
   share of worlds in which the T0 decision fires, and in which T1 (below, with 999 relabelings per
   world) reads OWN_GROUP_EXCESS. More than 10% = MISCALIBRATED (as registered).

4. **T1, own-group overlap excess (amended).** O(f, g) = |F_f n R_g| / |F_f|. D = mean over
   fledglings of O(f, own) - mean over the other 7 groups of O(f, g). No rarefaction is possible;
   every fledgling is compared with the same 8 records, and the null keeps record sizes. Null:
   reassign fledglings to groups at random keeping the number per group (604,800 distinct
   assignments, more than 100,000, so 9,999 random ones, as registered). Reading as v0.
   Fledglings with an empty F_f are dropped and counted.

5. **T2, order beyond content (amended).** v0's within-token order shuffle needs tokens and would in
   any case test whether order is shared population-wide, not whether it is group-specific. A
   **multi-order bag** is a multiset of units that appears in at least 2 distinct orders among the
   117 types. For fledgling f and record R, among f's types whose bag is multi-order and appears in R
   in some order, Q(f, R) = share that appear in R in exactly f's order. D_ord = mean over fledglings
   of Q(f, own) - mean over other groups g of Q(f, g) (terms with no bag match are skipped; a
   fledgling enters only if Q(f, own) and at least one Q(f, g) exist). Null: the same 9,999 group
   relabelings. Reading: D_ord > 0 with p <= 0.05 -> ORDER_CARRIES_GROUP_SIGNAL; otherwise
   ORDER_ADDS_NOTHING_MEASURABLE. UNINFORMATIVE if fewer than 10 distinct types produced by any
   fledgling have a multi-order bag (as registered), or if fewer than 5 fledglings enter.

6. **T4, contacts and repertoire size (amended).** No per-type token counts, so no rarefaction.
   Size = |F_f|; contacts = `total_soc_contacts` (Frequency_diversity_data.csv, constant per
   fledgling; the job checks this); effort = the sum over weeks of `combo_freq` (the fledgling's own
   recorded sequences). Statistic: Spearman rho between contacts and the residual of log(1 + size) on
   log(1 + effort) (least squares). Null 99,999 random permutations of contacts. Reading as v0
   (CONTACTS_EFFECT_SURVIVES_EFFORT if rho > 0 with p <= 0.05). The raw rho of contacts with size and
   the rho of effort with size are reported beside it. Ten fledglings: all four tests have low power.

7. **Power (amended).** For T1, 50 planted worlds: each fledgling type not in its own group's
   record is replaced, with probability max(0.10, observed D), by a random type from that record
   that the fledgling lacks. For T2, 50 planted worlds: each fledgling type whose multi-order bag is
   in its own record only in another order is switched to that order with probability 0.5. Each
   planted world is tested against 999 relabelings. Report the share detected at p <= 0.05.

8. **Secondary data.** T0, T1 and T3 are repeated on COMBINATION_repertoire_single_occurrences_removed.csv
   (the authors' own robustness file) and reported; they do not change any reading.

9. **Loader checks reported, not gating:** whether the follow records of the 3 JOG fledglings are
   identical; whether E_f equals the union of G_h over fledglings h of other groups; the share of
   F_f contained in G_f; whether every fledgling has the same 117 types.

Predictions stay as registered in PLAN.md (Q1 to Q5); Q3 now refers to the donor-world
calibration of the T0 decision. Added before the run: **Q6, the T0 decision fires in more than 50%
of donor worlds: prior 0.55.**

## Limits added

The group record of each fledgling was captured during that fledgling's own follows, so a shared
event (the same intruder, the same predator) puts a sequence in both F_f and G_f. T1 cannot separate
this from learning. Presence/absence data give no weight to how often a type was used.
