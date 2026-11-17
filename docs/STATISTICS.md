# Modification comparisons

## Observations and positions

The independent units used by the test are read identifiers. Each comparison
uses one modification type and one interval role: head for P–L, tail for P–R,
left flank for B–L, or right flank for B–R. Reads with ambiguous or complex
assignments do not enter the simple-state comparisons.

For each read orientation, usable positions are intersected across both states.
The same positions therefore contribute to every read fraction in that
orientation. Forward and reverse orientations can use different position sets.
Each read contributes the number of calls with ML ≥205 divided by the number
of shared positions. 5mC and 5hmC are never added together.

The parser retains explicit ML bytes. Omitted dot-mode positions use the
implemented low-call convention. Question-mode omissions remain unknown and
are excluded. Missing modification channels are not zero-valued measurements.
The parser currently handles forward-strand MM groups with letter codes
`C+m`, `C+h` and `A+a`. Other MM-strand encodings and numeric ChEBI codes do not
supply these named endpoints. Reference and read bases must agree at used sites.

## Effect and exact p-value

Within each eligible orientation, the effect is the mean alternative-state
fraction minus the mean baseline-state fraction. Eligible orientations need
at least one read in each state and at least one shared position. If both
orientations are eligible, their effects receive equal weights, regardless
of their read counts. `delta_pp` is this effect multiplied by 100.

Under the null hypothesis, complete read profiles are exchangeable between
states within each orientation. The test reallocates read labels, preserving
the state sizes within each orientation. It counts how many reallocations
give an absolute effect at least as large as the observed effect. This count
divided by the total number of reallocations is the exact two-sided p-value.
Ties are included. Repeated profiles retain their combinatorial multiplicities.

The implementation in `stratified_exact.py` uses integer counts and rational
arithmetic. It avoids Monte Carlo sampling. `p_numerator` and `p_denominator`
allow the reported probability to be checked independently. If the configured
dynamic-programming state limit is exceeded, the row has
`test_status=exact_state_limit` and no p-value. It is not replaced by p=1.

The main pipeline uses this exact implementation, not SciPy's default
two-sided tail convention. SciPy doubles the smaller one-sided p-value,
whereas metaphaser counts the absolute-effect tail.
[SciPy permutation-test documentation](https://docs.scipy.org/doc/scipy-1.16.1/reference/generated/scipy.stats.permutation_test.html).

## Multiple testing

`metaphaser adjust` selects rows with `test_status=tested`, `strand=combined`,
contrast P–L or P–R, and a host-linked structural class. It includes separate
modification types from every supplied sample table. It rejects duplicate
hypothesis keys, inconsistent exact probabilities and unresolved exact-limit
rows in that primary family.

BH is calculated with `scipy.stats.false_discovery_control(method="bh")`.
The output records the family size for each row. A correction using one
sample is a different family from a correction across an entire study.
BH's FDR guarantee depends on assumptions about the test dependence.
[SciPy BH documentation](https://docs.scipy.org/doc/scipy-1.16.1/reference/generated/scipy.stats.false_discovery_control.html).

## Scope

The fractions describe the positions shared by the included reads. Adding
reads can reduce that intersection. A large number of positions does not
increase the number of independent molecules. With small groups, exact
p-values can be coarse and uncertainty remains substantial. No bootstrap
confidence intervals are calculated by this version. `--bootstrap-replicates`
must be zero. `--seed` and `--permutations` are retained compatibility options
and do not change the exact calculation.
