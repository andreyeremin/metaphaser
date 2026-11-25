# Reproduce the 383 primary manuscript comparisons

This directory contains the frozen, study-derived inputs for all 383 primary
P–L/P–R comparisons: 152 for 5mC, 152 for 5hmC and 79 for 6mA. Nonsignificant
comparisons remain in the single Benjamini–Hochberg (BH) family.

## Run from a source checkout

Use a checkout containing this directory. The original `v0.2.0rc2` tag predates
the addition of the study tables. These commands run the statistical modules
from the current checkout explicitly; minimap2, modBAMs, assembly files and a
GPU are not needed for this saved-table calculation.

From the repository root, using Python 3.12:

```bash
python3.12 -m venv .venv-study
. .venv-study/bin/activate
python -m pip install numpy==1.26.4 scipy==1.16.1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PYTHONPATH=src python study/383/reproduce.py \
  --counts study/383/counts.tsv \
  --reference study/383/primary_comparisons.tsv \
  --output results/study-383
```

Use a fresh output directory. The script refuses existing output paths or
inputs whose SHA-256 hashes differ from the frozen study files. A successful
run writes `recomputed_primary.tsv` and `verification.json` into the specified
output directory and prints the verification summary. The JSON must report
`status: passed`, 383 matching numerators, denominators, effects, p-values and
BH values, 9,793 profiles and 661 orientation strata. It records the imported
metaphaser version and the hashes of the statistical source modules actually
used. The calculation preserves exact integer numerator/denominator pairs;
floating-point effects are checked within 1e-10 percentage points and p/q
values within relative tolerance 1e-14.

Expected summary: 10 nominal p-values below 0.05 and one q-value below 0.05.
OI-2, contig_3701, 5mC P–R has p = 1/72072, q =
0.005314130314130314 and a P-minus-R difference of
−8.170640834575261 percentage points.

## Files and units

| File | Content |
|---|---|
| `counts.tsv` | 9,793 per-read integer high-call counts and shared-site denominators at ML ≥205 |
| `primary_comparisons.tsv` | All 383 reference effects, exact probabilities and `q_BH_383` values |
| `reproduce.py` | Complete-family input checks, exact tests, BH correction and reference comparison |
| `historical_settings.json` | Source-checked historical settings and their differences from rc2 |
| `PROVENANCE.json` | Unchanged source-table derivation and source hashes from 13 November 2026 |
| `verification-20261113.json` | Archived successful rc1 reproduction; preserved as historical evidence |
| `verification.json` | Local verification of the current checkout against the frozen inputs |
| `LICENSE.md` | Scope of the code licence and the separately licensed study material |

`sample`, `locus`, `contrast` and `modification` identify a hypothesis.
`orientation` is `+` or `-`; `state` is P or its L/R baseline.
`profile_index` is one-based within hypothesis/orientation/state. Equal counts
are distinct molecule profiles and retain their permutation multiplicities.
Indices cannot link molecules between comparisons. Original read identifiers,
sequences, raw signals and per-position call matrices are not included.
The reference table retains authorised study sample/locus labels and interval
coordinates. `nominal_p_below_005` means nominal p <0.05; the study-wide
significance criterion is `q_BH_383 < 0.05`.

The OI inputs derive from 237 saved shared-position ML comparisons thresholded
at ML ≥205; SETUN contributes 146 saved integer-count comparisons. Within an
orientation, every included read has the same shared-site denominator. Eligible
orientations receive equal weights, independent of their read counts. The
whole-read exact test includes ties in the absolute-effect tail. See
[the statistical method](../../docs/STATISTICS.md).

## Historical settings and the public release

The [study reproduction guide](../../docs/STUDY_REPRODUCTION.md) separates
historical discovery settings from current defaults. `historical_settings.json`
is a descriptive record, **not** a file accepted by `metaphaser run --parameters`.
Running rc2 discovery again is not a reconstruction of the archived 383-test
family. These tables reproduce the downstream statistics after read assignment
and shared-site selection; they do not rerun assembly, discovery, mapping,
assignment or the selection of shared positions.

Count totals cannot reconstruct partial-coverage positions or missingness.
An experimental position-centred analysis needs a separate read-by-position
matrix and its own missingness checks and statistical family. It does not
replace this primary endpoint.

## Licence

The reproduction code and documentation retain the repository's MIT licence.
The original study package specified no licence for its derived tables.
[LICENSE.md](LICENSE.md) preserves that distinction without assigning a new
data licence or changing the existing software licence.
