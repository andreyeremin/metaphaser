# Partial-coverage sensitivity

This optional analysis uses saved `modifications/site_calls.tsv.gz` from a
completed run. It performs no mapping, basecalling or BAM access and needs no
GPU. It leaves the primary method and its multiple-testing family unchanged.

```bash
metaphaser partial-coverage \
  --site-calls RUN/modifications/site_calls.tsv.gz \
  --sample SAMPLE --locus LOCUS --modification 5mC --contrast P-R \
  --min-site-reads 6 --exact-state-limit 200000 \
  --output partial_coverage.json
```

Filters for sample, locus, modification and contrast are optional. Omitting
them processes every represented eligible comparison. Use the filters for a
bounded sensitivity check. Existing outputs are never overwritten.

## What is compared

Within each orientation, positions need usable calls in at least six pooled
reads, without using the state labels for selection. For each observed call,
subtract the pooled mean at that position. Each read receives the mean of its
observed centred values. Compare the state means, then give each orientation
equal weight. P–L uses the interval head and P–R its tail.

The position set, pooled means, whole-read values and group sizes stay fixed
during permutations within orientation. Exact arithmetic includes ties and
counts assignments of individual reads even when their values are equal.
The output effect is in percentage points. It differs from the primary
shared-position effect because both positions and per-read denominators can differ.

## Input and missing calls

Each input row identifies sample, locus, read, junction, role, modification,
orientation, contig, position, call status and high call. The current site-call
output and historical tables lacking `score_used` are supported. Channels are
processed separately. Unknown calls remain missing. They are never changed to zero.
Identical duplicate observations collapse; conflicting duplicates are errors.

All matching reads represented in the input are retained. If a represented
read has no usable selected observation, or a represented orientation lacks
a state, the comparison is reported as untestable. A read entirely absent
from the site-call table cannot be reconstructed by this command. A single
represented orientation can be tested and is explicitly flagged.

## Output and interpretation

The JSON report includes read counts, selected and informative sites, per-read
coverage and residuals, group-only sites, overlap measures, coverage-mask
differences, effect size and exact permutation numerator and denominator.
Results are nominal exploratory p-values without automatic BH correction.
Define the tested family explicitly before interpreting multiple comparisons.

The null requires exchangeability of whole observed profiles, including their
coverage patterns. Site centring does not correct state-dependent missingness.
The tests include an example with identical complete profiles in both states
but selective loss of low calls in one state. It produces an apparent
−50 percentage-point difference and p=1/462, demonstrating the limitation.
Coverage warnings are descriptive, and no warning is not proof of exchangeability.
Total contig depth alone does not justify preferring this mode.

The exact dynamic-programming state limit is enforced before allocating a
larger layer. Exceeding it produces an explicit failure and no substituted
p-value. Incomplete/untestable reports return exit status 2 and preserve
diagnostics. A completed report returns 0.

## Validation

Lightweight tests compare exact results with full enumeration on small
matrices and establish equivalence to the primary effect and p-value under
complete coverage. They also test fixed position selection, missing calls,
tied scores, orientation weights, duplicate conflicts and resource limits.
A bounded saved-data check for OI-2/J1213/5mC/P–R reproduced the historical
partial-coverage effect of −5.167800682955004 percentage points and p=1/72072.
This does not reclassify the primary 383 hypotheses or validate all missingness patterns.
