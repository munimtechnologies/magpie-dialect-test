# magpie-dialect-test: loader notes (2026-11-24, before any statistic was computed)

MV2 of Munim, Inc. Made with AI agents working in his name. CC BY 4.0.

## Source and checks

figshare doi:10.6084/m9.figshare.28251398, version 2, "Data and code: Ontogenetic evidence of
socially learned call sequences in Western Australian magpies" (Mason, Ridley, King), licence
CC BY 4.0 (figshare API, read 2026-11-24). The job downloads from ndownloader.figshare.com and checks
each file against the size and MD5 figshare publishes:

| file | figshare file id | bytes | MD5 |
|---|---|---|---|
| COMBINATION_repertoire.csv | 51834599 | 78,432 | f60a6af46f0161225b5965eb2c917874 |
| COMBINATION_repertoire_single_occurrences_removed.csv | 55548710 | 47,142 | 85264331328a82127c4b2d807ce88a36 |
| Frequency_diversity_data.csv | 51834608 | 48,384 | 08feedc34b5446c0fbf778a262277361 |
| DISCRETE_repertoire.csv | 51834602 | 3,308 | 5be86be6a977aea606fbad023359db7f |
| total_repertoire.csv | 51834617 | 90,938 | 7585bf3be946d86b12871754a39a3ead |
| Repertoire modelling.R | 55548707 | 8,243 | 1d029c38bdbf48ba0ea01f9b60852503 |
| Frequency and diversity modelling.R | 55548704 | 41,050 | 9d7f7839cb4b2a5e4975a57d70d29996 |
| README.txt | 59052542 | 1,067 | e2bf57265a574bd433ea8feea058022d |

A copy saved by MV2's desktop agent on 2026-11-23 has the same byte sizes. The job records the
SHA-256 of what it downloads in fetch_report.json.

## Columns used

COMBINATION_repertoire.csv (and the single-occurrences-removed file, same layout):

| column | meaning (authors' README and R script) | used as |
|---|---|---|
| ID | fledgling | fledgling |
| group | natal group | group |
| mother, closest_contact, closest_contact_FL0_A1, closest_adult_contact | social partners | not used |
| total_number_of_contacts, proportion_spent_in_contact | social measures | reported only |
| call | sequence type: ordered call units separated by spaces | type; its sorted units = bag |
| fl01 | fledgling produced the type (0/1) | F |
| group01 | type recorded from the natal group during the fledgling's follows (0/1) | G |
| exgroup01 | type recorded in the population outside the group (0/1) | E |

Frequency_diversity_data.csv: one row per fledgling and week. Used: `ID`, `total_soc_contacts`
(expected constant per fledgling; checked), `combo_freq` (the fledgling's sequences that week; summed
as effort). All other columns unused.

## Shape (row counts and labels only)

- 1,170 data rows = 10 fledglings x 117 sequence types; every fledgling has the same 117 labels.
- Fledglings and groups: BITS (RVD), GONFL (GON), MBPX, TONY and YXBM (JOG), MONKEY (SCL), NHFL (NH),
  PRFL (PR), PXYM (BWYb), SOLID (MBG). 8 groups. The paper reports 11 fledglings.
- 7 call units: DS, LH, LHDS, NL, NLDS, SH, SHDS. Type lengths run from 2 to 15 units. Labels keep
  the order of production (for example both `DS NL` and `NL DS` exist), so order can be tested.
- No dates, no token counts, no per-type frequencies: presence/absence only. This is why
  AMENDMENT-1 replaces the token-based parts of PLAN.md.

## What was seen

Header lines; the first 2 or 3 data rows of each CSV; the 117 type labels; for each fledgling, the
metadata fields of its row for type `DS LH`. No 0/1 column was counted, summed or compared, and no
model was fitted, before AMENDMENT-1 and the code were written.
