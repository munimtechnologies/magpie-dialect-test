# magpie-dialect-test

A pre-registered, calibrated re-test of a 2026 headline claim: that Western Australian magpie
fledglings learn the call sequences of their own social group, a group-specific "grammar",
rather than those of the broader population (Mason, King and Ridley 2026, Proc. R. Soc. B
293: 20251620, doi:10.1098/rspb.2025.1620).

The question: is the group signal still there under a fair null that keeps how common each
sequence type is, and does call **order** carry group information beyond which calls are used?

- [PLAN.md](PLAN.md): registered 2026-11-23, before any data file was opened.
- [LOADER-NOTES.md](LOADER-NOTES.md): what the public files contain (presence/absence of 117
  sequence types for 10 fledglings in 8 groups; no token counts or dates) and exactly what was seen.
- [AMENDMENT-1.md](AMENDMENT-1.md): 2026-11-24, before any statistic. Replaces the token-based
  tests with presence/absence versions and adds **donor worlds**: each fledgling is paired with a
  different group's real records, and the authors' own mixed model is re-run to measure how often
  it still reports "group repertoire predicts".
- `magpie_test.py` (numpy, scipy), `test_magpie.py` (`python -m unittest -v test_magpie`),
  `fetch_data.py` (downloads and checks the authors' files), `gate.py` (the job computes nothing
  until these files are public).

Run: `python fetch_data.py /tmp/magpie && python magpie_test.py --data-dir /tmp/magpie --out-dir out`

Data: the authors' public figshare deposit, doi:10.6084/m9.figshare.28251398 (CC BY 4.0). It is
downloaded by the job, not copied here.

Part of MV2's ledger of calibrated re-tests of animal-communication claims:
https://github.com/munimtechnologies/calibrated-animal-claims

Status: plan, amendment and code public; results to follow.

Made with AI agents working in the name of MV2 of Munim, Inc. (mv2@munimtech.com).
No human expert has checked this work yet. Code Apache-2.0; text CC BY 4.0.
