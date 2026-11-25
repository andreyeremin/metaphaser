# metaphaser

Classify nanopore reads by genomic junction state and compare DNA modification
calls between those states using an assembly FASTA and a modBAM.

**Version:** 0.2.0rc3 · **License:** MIT · [Русская инструкция](README.ru.md)

metaphaser finds candidate intervals, constructs synthetic junction references,
assigns reads to L, R, P and B, and compares **5mC, 5hmC and 6mA separately** at
the same DNA positions. It writes read assignments, interval tables, per-read
modification profiles, exact p-values and an HTML report. A separate command
applies Benjamini–Hochberg correction across samples.

## Install

Requirements: Linux or macOS, Git, Conda or Mamba, Python 3.11 or later, and
minimap2. The commands below use Python 3.12 and minimap2 2.30. Native Windows
execution is not supported. Use a Linux environment such as WSL2 instead.

```bash
git clone --branch v0.2.0rc3 https://github.com/andreyeremin/metaphaser.git
cd metaphaser
conda env create -f environment.yml
conda activate metaphaser
python -m pip install -r requirements-tested.txt
python -m pip install --no-build-isolation --no-deps .
metaphaser --version
metaphaser doctor
```

`requirements-tested.txt` pins the Python dependencies used for validation.
SAMtools functionality is provided through pysam. Basecalling and assembly
are performed before running metaphaser.

If Python and minimap2 are already installed, you can use a virtual environment:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-tested.txt
python -m pip install --no-build-isolation --no-deps .
metaphaser doctor
```

minimap2 must be on `PATH`. Alternatively, pass its location with
`--minimap2 /path/to/minimap2` to `doctor`, `demo` and `run`.

## Try the included example

```bash
metaphaser demo --outdir demo_run --threads 2
```

This generates 56 synthetic reads and runs the complete analysis. The expected
result includes one interval with all four junction types. Open
`demo_run/results/report.html` and inspect `demo_run/results/elements.tsv`.
Use a new output directory for each demo invocation.

## Analyse your sample

Supply the assembly and **original reads with complete sequences and MM/ML
modification tags**. An unaligned modBAM is accepted. FASTQ alone cannot supply
modification measurements.

Check the inputs first:

```bash
metaphaser check --assembly assembly.fasta --modbam reads.mod.bam
```

Run directly in the terminal:

```bash
metaphaser run \
  --assembly assembly.fasta \
  --modbam reads.mod.bam \
  --sample SAMPLE \
  --cytosine-context all \
  --modifications 5mC,5hmC,6mA \
  --threads 4 --validation-workers 2 \
  --outdir results/SAMPLE
```

Choose the context and channels from the **basecalling model**, not the sample
name. For CpG-only 5mC/5hmC calls, use:

```bash
metaphaser run \
  --assembly assembly.fasta \
  --modbam reads.mod.bam \
  --sample SETUN06 \
  --cytosine-context cpg \
  --modifications 5mC,5hmC \
  --threads 4 --validation-workers 2 \
  --outdir results/SETUN06
```

`--modifications auto` uses available supported channels. An absent channel
is not reported as a zero modification level. `--cytosine-context` defaults
to `cpg`, so specify `all` for all-context cytosine comparisons.

If the BAM was aligned to a different assembly, add
`--ignore-original-alignments` to both `check` and `run`. The program then
uses fresh mappings for discovery. To inspect sequence-only BAMs, explicitly
use `--allow-no-modifications` with `run`.

## Run in the background and resume

The CLI normally runs in the foreground. To continue after disconnecting SSH:

```bash
nohup metaphaser run \
  --assembly assembly.fasta --modbam reads.mod.bam \
  --sample SAMPLE --cytosine-context all --modifications auto \
  --threads 4 --validation-workers 2 --outdir results/SAMPLE \
  > SAMPLE.log 2>&1 < /dev/null &
echo $! > SAMPLE.pid
tail -f SAMPLE.log
```

Pressing Ctrl+C during `tail` stops viewing the log, not the background job.
The program prints English progress messages, stage durations and the final
report path. Completion is also recorded in `manifest.json` and `summary.json`.

To resume an interrupted run, repeat its **same command** with `--resume`.
Input identities, software and scientific parameters must match. Do not
upgrade the software in the middle of a run. A completed output directory is
protected from accidental replacement. Mapping and validation parallelism are
controlled by `--threads` and `--validation-workers`, respectively.

## Read the results

| File | What it contains |
|---|---|
| `report.html` | Summary to open in a browser |
| `elements.tsv` | Intervals with at least 3 simple P reads and 3 simple L or R reads |
| `read_classifications.tsv.gz` | One row per input read identifier, including exclusions |
| `read_locus_classifications.tsv.gz` | Read assignments within each interval |
| `junction_calls.tsv.gz` | Junction calls and ambiguity or complexity flags |
| `loci.tsv`, `regions.tsv` | All candidates and regions formed by overlapping retained intervals |
| `element_intervals.fasta` | Assembly sequences of retained intervals |
| `modifications/comparisons.tsv` | Per-orientation summaries and combined-orientation exact tests |
| `modifications/read_profiles.tsv.gz` | Per-read modification fractions at shared positions |
| `modifications/site_calls.tsv.gz` | Individual calls and omission status |
| `modifications/panels.tsv.gz` | DNA positions used in each comparison |
| `modifications/epihaplotypes.tsv.gz` | Modification patterns counted separately by channel and state |
| `modifications/channels.json` | Available and measured modification channels |
| `manifest.json`, `input_metadata.json`, `logs/` | Settings, provenance, status and logs |

**L** crosses the interval start and **R** its end, along increasing assembly
coordinates. **P** connects the interval end to its start. **B** connects the
external flanks. These are junction labels, independent of read orientation.
A read crossing both original boundaries can contribute to both L and R.
The `*_simple` columns count reads that pass the final unambiguous, non-complex
classification. Coordinates are zero-based, with the right endpoint excluded.
Overlapping intervals are kept separately. Their count is not a count of
distinct biological elements.

## Modification comparisons and multiple testing

For each modification and orientation, both states are evaluated at positions
with usable calls in every included read. The endpoint is the per-read
fraction of calls with **ML ≥205**. The difference is alternative minus baseline:
for example, P minus L in `P-L`, expressed in **percentage points**.

The exact test permutes whole-read labels within each orientation, preserves
group sizes and gives available orientations equal weights. The two-sided
p-value counts absolute differences at least as large as observed, including
ties. A test needs at least one read per state and one shared position in each
included orientation. Small groups can have coarse attainable p-values.
See [statistical details](docs/STATISTICS.md).

The `run` command reports unadjusted p-values. Apply BH to a declared study
family by supplying **every sample** comparison file:

```bash
metaphaser adjust \
  --comparisons results/SAMPLE_A/modifications/comparisons.tsv \
                results/SAMPLE_B/modifications/comparisons.tsv \
  --output results/primary_BH.tsv
```

This selects tested combined-orientation P–L and P–R comparisons at retained
host-linked intervals, across the separate modification types. It excludes B
comparisons and descriptive orientation rows. Output includes `q_BH`,
`BH_family_size` and `BH_significant_005`. Do not select input tables by their
p-values. Use `--expected-tests N` to enforce a pre-established family size.

## Methods, validation and troubleshooting

- [Historical manuscript settings and reproduction of all 383 primary tests](docs/STUDY_REPRODUCTION.md)
- [Structural method and thresholds](docs/STRUCTURAL_METHOD.md)
- [Discovery changes in 0.2.0rc2](docs/DISCOVERY_FIX_20261117.md)
- [Validation and known limitations](docs/VALIDATION.md)
- [Troubleshooting](docs/TROUBLESHOOTING.md)
- [Release notes](CHANGELOG.md)

The package includes synthetic input data and regression tests. To run them
after installation, with minimap2 on `PATH`:

```bash
python -m pytest -q
```

An optional [GitHub Actions workflow](docs/github-actions-example.yml) builds
minimap2 and runs the tests in a clean Linux runner. To enable it in a fork,
copy it to `.github/workflows/tests.yml` with an account allowed to edit workflows.

Version 0.2.0rc3 remains a release candidate. The command-line workflow includes
structural classification, separate-channel modification comparisons and BH
correction. Read-quality sweeps, motif searches and gene annotation are not
implemented as CLI subcommands.

## Citation and license

Citation metadata are in [CITATION.cff](CITATION.cff). Cite the version you used.
The code is distributed under the [MIT license](LICENSE). The bundled example
contains generated sequences and synthetic modification calls.
The separate [study-derived tables](study/383/README.md) retain their
[original data-licensing status](study/383/LICENSE.md); the software licence
does not assign a new licence to those tables.

## Partial-coverage sensitivity

Use saved calls to examine positions covered by at least six pooled reads per
orientation. This optional mode reports an exploratory comparison separately
from the primary shared-position tests. It does not access BAM files or use a GPU.

```bash
metaphaser partial-coverage \
  --site-calls RUN/modifications/site_calls.tsv.gz \
  --sample SAMPLE --locus LOCUS --modification 5mC --contrast P-R \
  --output partial_coverage.json
```

Choose a new output file. See [the method, diagnostics and limits](docs/PARTIAL_COVERAGE.md).
