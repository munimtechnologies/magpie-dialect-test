# Magpie group-sequence re-test: results

Made with AI agents working in MV2's name (MV2, Munim, Inc.); not yet checked by a human expert.

## Primary file (COMBINATION_repertoire.csv)

- T0 group model: b1 2.251367486484622, z 9.286844079767961, p 1.5893080749713055e-20; AICc group 631.34, pop 721.00, base 719.00; published-style decision fires: True
- T1: D 0.1467, null mean 0.0421, null q95 0.1227, p 0.0124 -> OWN_GROUP_EXCESS
- T2: D_ord 0.034677139439044234 over 9 fledglings, p 0.1612, multi-order types produced 38 -> ORDER_ADDS_NOTHING_MEASURABLE
- T3 donor worlds: 200; T0 decision fired in 200 of 200 (MISCALIBRATED); T1 fired in 1
- Power: {"q_t1": 0.14667563477087284, "t1": {"worlds": 50, "hits": 50, "power": 1.0}, "t2": {"worlds": 50, "hits": 7, "power": 0.14}}
- T4: rho adjusted -0.8786687791935092, p 0.9988, raw rho -0.004219446843148565, effort-size rho 0.9412097036339401 -> CONTACTS_EFFECT_NOT_SEPARABLE_FROM_EFFORT

## Secondary file (COMBINATION_repertoire_single_occurrences_removed.csv)

- T0 decision fires: True (b1 2.4314373058152148, p 1.9747566280303524e-18)
- T1: D 0.1756, p 0.0042 -> OWN_GROUP_EXCESS
- T3: T0 fired in 200 of 200 donor worlds; T1 in 0

## Predictions

```
{
 "outcomes": {
  "Q1_T1_own_group_excess": true,
  "Q2_T2_order_carries_group_signal": false,
  "Q3_T0_decision_miscalibrated": true,
  "Q4_T4_contacts_survive_effort": false,
  "Q5_T2_uninformative": false,
  "Q6_T0_fires_in_over_half_of_donor_worlds": true
 },
 "priors": {
  "Q1_T1_own_group_excess": 0.6,
  "Q2_T2_order_carries_group_signal": 0.25,
  "Q3_T0_decision_miscalibrated": 0.4,
  "Q4_T4_contacts_survive_effort": 0.35,
  "Q5_T2_uninformative": 0.3,
  "Q6_T0_fires_in_over_half_of_donor_worlds": 0.55
 },
 "brier": {
  "Q1_T1_own_group_excess": 0.16000000000000003,
  "Q2_T2_order_carries_group_signal": 0.0625,
  "Q3_T0_decision_miscalibrated": 0.36,
  "Q4_T4_contacts_survive_effort": 0.12249999999999998,
  "Q5_T2_uninformative": 0.09,
  "Q6_T0_fires_in_over_half_of_donor_worlds": 0.20249999999999996
 },
 "mean_brier": 0.16624999999999998
}
```

## Loader checks

```
{
 "rows": 1170,
 "fledglings": 10,
 "types": 117,
 "groups": 8,
 "all_cells_present": true,
 "fledglings_per_group": {
  "BWYb": 1,
  "GON": 1,
  "JOG": 3,
  "MBG": 1,
  "NH": 1,
  "PR": 1,
  "RVD": 1,
  "SCL": 1
 },
 "multi_fledgling_group_records_identical": {
  "JOG": false
 },
 "exgroup_equals_union_of_other_groups_share": [
  1.0,
  0.9914529914529915,
  0.9658119658119658,
  1.0,
  1.0,
  1.0,
  1.0,
  0.9743589743589743,
  0.9829059829059829,
  0.9572649572649573
 ],
 "share_of_F_in_own_follow_record": [
  0.7222222222222222,
  0.3333333333333333,
  1.0,
  0.75,
  null,
  0.3409090909090909,
  0.5714285714285714,
  0.4,
  0.5,
  0.38461538461538464
 ],
 "repertoire_sizes": [
  18,
  6,
  6,
  4,
  0,
  44,
  21,
  10,
  4,
  13
 ],
 "follow_record_sizes": [
  41,
  3,
  11,
  22,
  15,
  22,
  28,
  12,
  13,
  10
 ]
}
```
