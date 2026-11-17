# Validation of 0.2.0rc2

The source code in this release is the code tested on 17 November 2026.
Publication packaging adds documentation and installation metadata without
changing the analysis modules.

## Regression tests

All **44 tests passed**, with no failures or skips, in the Linux validation
environment. They include synthetic end-to-end runs, short intervals, all
four junction types, empty output, absent modification channels, MM/ML parsing,
source-copy ambiguity, parallel validation, exact probabilities, BH selection,
interruption/resume and the two discovery corrections.

## Full-sample check

A SETUN06 analysis processed **118,278 read identifiers**. Exported reads and
assembly mappings were verified and reused. Candidate discovery, synthetic
mapping, validation, modification extraction and statistics were rerun.

| Interval | Coordinates, zero-based half-open | L | R | P | B |
|---|---|---:|---:|---:|---:|
| contig_116 | [544868, 887996) | 11 | 10 | 3 | 0 |
| contig_44 | [72196, 636770) | 0 | 9 | 3 | 0 |

Both intervals had been missed by the preceding candidate. All nine previously
retained intervals remained present, giving eleven retained intervals in seven
overlapping regions. Three intervals met the read threshold at L, R and P.
None met it at all four junction types.

An independent integer-count reconstruction confirmed all **242 tested
p-values**, their exact numerators and denominators, the effects and the
numbers of reads. Thirty of these tests were P–L/P–R comparisons at retained
host-linked intervals. Changes to assignments at other intervals were checked:
one formerly simple P read also had an L junction, while three formerly
assigned L/R calls failed the unchanged edit-margin rule against added
alternative references.

These checks establish the reported behaviour on this sample and the test
fixtures. They are not a sensitivity/specificity benchmark across organisms
or sequencing models. The repository includes generated example data, not
the biological input files.

## Known runtime issue

The first full-sample attempt terminated with SIGSEGV during read recruitment.
The cause was not established. Three limited repetitions of the transition
from boundary refinement to BAM reading completed successfully. The resumed
full analysis and independent final checks also completed successfully.

If a similar crash occurs, retain the output directory and logs. Enable Python
fault traces with `PYTHONFAULTHANDLER=1`, and resume only with the same software,
inputs and analysis settings. Include the manifest, runtime versions and
relevant error trace when reporting the issue, after reviewing them for
private paths or sample metadata.
