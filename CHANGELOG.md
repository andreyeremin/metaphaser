# Changelog

## 0.2.0rc3 — 2026-11-25

- Add the optional partial-coverage sensitivity command for saved site calls.
- Report exact orientation-stratified permutation probabilities, coverage diagnostics and explicit untestable outcomes.
- Include historical study settings and frozen inputs reproducing all 383 primary tests and their BH correction.
- Preserve primary discovery, assignment and shared-position calculation behaviour.


## 0.2.0rc2 — 2026-11-17

- Separate the preliminary discovery threshold of two distinct reads from the
  final requirement of three validated P reads and three L or R reads.
- Refine candidate boundaries using the sequence between split alignments.
  Preserve the original hypotheses and add up to five improved boundary pairs.
- Keep final junction filters, modification extraction and the exact test unchanged.
- Add six regression cases. All 44 tests passed on Linux.
- Verify recovery of two previously missed SETUN06 intervals in a full-sample
  analysis and independently reconstruct all 242 tested modification comparisons.
- Provide installation instructions, a synthetic example and output documentation.

An initial full-sample attempt terminated with SIGSEGV. The cause was not
established. Three limited diagnostic repetitions and the resumed full-sample
analysis completed successfully. See [validation](docs/VALIDATION.md).

## 0.2.0rc1 — internal validation candidate

- Package the separate-channel modification analysis and exact
  orientation-stratified test.
- Add study-wide Benjamini–Hochberg correction with explicit family selection.
- Validate synthetic inputs, interruption/resume and saved statistical endpoints.
