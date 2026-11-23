# magpie-dialect-test: pre-registered plan (v0, 2026-11-23)

MV2 of Munim, Inc. Made with AI agents working in his name. Apache-2.0 (code), CC BY 4.0 (text).

**Status: registered BEFORE any data was opened.** At the time this file was written, nobody
working on this repository had downloaded, opened or summarised any file of the data set below.
The only things known are the paper's public abstract, its press coverage and the file list on
figshare. Column names and file layout will be read next (headers and row counts only, no
statistic computed) and recorded in LOADER-NOTES.md before any test runs; anything decided
after that is an AMENDMENT with its own date.

## The claim

Mason, King and Ridley (2026), "Ontogenetic evidence of socially learned call sequences in
Western Australian magpies", Proc. R. Soc. B 293: 20251620, doi:10.1098/rspb.2025.1620
(CC BY 4.0). From the abstract: fledglings acquire their sequence repertoires gradually and
"produce the sequences of their social group - not the broader population", which the authors
read as social learning; fledglings with more social contacts produce more sequence types.
Press coverage called it the first evidence of an animal learning call-sequence structure, a
group-specific "grammar" or "dialect". A June 2026 preprint by the same group ("Innate building
blocks underly socially learned call sequences") reports that the individual call segments are
shared across groups from the first week, so the claim now rests on how calls are combined.

Data: Mason, Ridley and King, "Data and code: Ontogenetic evidence of socially learned call
sequences in Western Australian magpies", figshare doi:10.6084/m9.figshare.28251398 v2
(CC BY 4.0): COMBINATION_repertoire.csv, COMBINATION_repertoire_single_occurrences_removed.csv,
DISCRETE_repertoire.csv, Frequency_diversity_data.csv, total_repertoire.csv, two R scripts,
README.txt. 11 fledglings were followed for about 200 days.

We searched for an outside re-analysis, comment or reply (Proc B, bioRxiv, Sciety, news) on
2026-11-23 and found none.

## Why a re-test

Three features of the design can make sequences look group-specific without any learning of
combinations:

1. **Shared events.** Fledglings are recorded in the same sessions as their own group's adults,
   and the sequences are threat calls (intruders, foxes, aerial predators). Birds answering the
   same event give the same sequences.
2. **Sampling effort.** The number of sequence types an individual or group shows grows with the
   number of sequences recorded; unequal effort makes some repertoires look broader and some
   overlaps look larger.
3. **Content, not order.** If groups differ in which calls they use, sequences differ by group
   even when the order of calls is learned from nobody. The "grammar" reading needs order to
   carry group information beyond content.

## Units (fixed now, mapped to columns in LOADER-NOTES.md)

- A **token** is one recorded sequence (combination) by one bird.
- Its **type** is the authors' combination label (ordered string of call units).
- Its **bag** is the multiset of call units in the type, ignoring order.
- **Primary** data: COMBINATION_repertoire.csv as published. **Secondary**: the authors'
  single-occurrences-removed file. Adults and fledglings are the authors' own labels; group is
  the authors' group label. Birds of unknown group are dropped and counted.

## Tests

Random seeds fixed in code; all p values one-sided, (1 + hits) / (1 + worlds).

**T0, reproduction.** Run the authors' own group-specificity comparison as their R script does it
(the statistic and the decision rule are copied from the script, not guessed) and record its
value. If the script cannot be run in the container, T0 is reported as NOT_RUN and the
calibration in T3 uses the description in the paper.

**T1, own-group overlap excess (primary).** For each fledgling f, O_own(f) is the share of f's
types found among the adult tokens of f's own group; O_other(f) is the mean over the other groups
of the same share. Adult tokens of every group are rarefied to m = the smallest adult token count
among groups with a fledgling (mean over 200 random subsamples). D = mean over fledglings of
O_own - O_other. Null: reassign fledglings to groups at random, keeping the number of fledglings
per group (all assignments if at most 100,000, otherwise 9,999 random ones).
Reading: D > 0 with p <= 0.05 -> OWN_GROUP_EXCESS; otherwise NO_OWN_GROUP_EXCESS.

**T2, order beyond content.** D_type is T1's D; D_bag is the same statistic with every type
replaced by its bag. R = D_type - D_bag. Null: shuffle the order of call units inside every token
of every bird (999 worlds), recompute R. Reading: R above the null with p <= 0.05 ->
ORDER_CARRIES_GROUP_SIGNAL; otherwise ORDER_ADDS_NOTHING_MEASURABLE. If fewer than 10 fledgling
types have a bag that is seen in more than one order anywhere in the data, T2 is
UNINFORMATIVE (order is nearly fixed by content) and that count is reported.

**T3, calibration of the published-style test.** 200 null worlds in which every fledgling's
tokens are drawn from the pooled adult tokens of ALL groups (keeping each fledgling's token count
and, where dates exist, its sampling schedule). Report how often the T0 decision rule and T1 say
"own group" at p <= 0.05. If T0 is NOT_RUN, only T1 is calibrated. A test that fires in more
than 10% of these worlds is MISCALIBRATED.

**T4, social contacts and repertoire size.** For each fledgling, the number of distinct types in
a random subsample of k tokens (k = the smallest token count among fledglings with at least 20
tokens; mean over 200 subsamples; fledglings below 20 tokens are dropped and counted), against
the fledgling's number of social contacts as the authors measured it. Spearman rho; null: 99,999
random permutations of contacts over the fledglings kept (all orderings if there are 8 or
fewer). Reading: rho > 0 with p <= 0.05 -> CONTACTS_EFFECT_SURVIVES_EFFORT; otherwise
CONTACTS_EFFECT_NOT_SEPARABLE_FROM_EFFORT. The raw (unrarefied) correlation is reported beside it.

**Power.** For T1 and T2, 50 worlds with a planted own-group effect: each fledgling token is
replaced, with probability equal to the observed O_own - O_other (floored at 0.10), by a token
drawn from its own group's adults. Report the share detected.

## Predictions (priors, for a Brier score at grading)

| Q | Statement | Prior |
|---|---|---|
| Q1 | T1 reads OWN_GROUP_EXCESS | 0.60 |
| Q2 | T2 reads ORDER_CARRIES_GROUP_SIGNAL | 0.25 |
| Q3 | T3 finds the published-style test MISCALIBRATED | 0.40 |
| Q4 | T4 reads CONTACTS_EFFECT_SURVIVES_EFFORT | 0.35 |
| Q5 | T2 is UNINFORMATIVE | 0.30 |

## Limits stated in advance

Eleven fledglings is a small sample; every test here has low power, and an absence of evidence
will be reported as that. Shared events cannot be separated from learning with these files: a
positive T1 is consistent with both. Labels (calls, groups, contacts) are the authors'. No human
expert has checked this plan yet.
