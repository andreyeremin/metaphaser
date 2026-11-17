# Troubleshooting

## `metaphaser` is not found

Activate the environment where the package was installed:

```bash
conda activate metaphaser
python -m metaphaser --version
```

Use `python -m pip` from that same environment when installing.

## Install on a Linux server without internet access

The target still needs Python 3.12 and minimap2 2.30. On an internet-connected
machine, download the dependencies for Linux x86-64 and CPython 3.12:

```bash
python3.12 -m pip download --only-binary=:all: \
  --platform manylinux_2_28_x86_64 --platform manylinux2014_x86_64 \
  --python-version 312 --implementation cp --abi cp312 \
  -r requirements-tested.txt --dest wheelhouse
```

Copy the source directory and `wheelhouse` to the server. From the copied
source directory, with minimap2 on `PATH`:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --no-index --find-links wheelhouse -r requirements-tested.txt
python -m pip install --no-build-isolation --no-deps .
metaphaser doctor
```

These wheel selectors target Linux x86-64 with glibc 2.28 or later. They are
not appropriate for ARM or an older glibc. Keep the Python version and wheel
platform consistent with the target server.

## minimap2 is not found or rejects an option

Use the supplied environment file, which selects minimap2 2.30, and run
`metaphaser doctor`. An older executable may appear earlier on `PATH`.
The `--minimap2` option accepts an explicit executable path.

## Input assembly does not match the original BAM reference

Confirm that the assembly and reads belong to the intended sample. If the BAM
was aligned to a different assembly, add `--ignore-original-alignments`.
The input must still retain complete original read sequences and MM/ML calls.
Hard-clipped primary sequences are excluded.

## No modification calls or no 6mA results

Inspect `input_metadata.json` and `modifications/channels.json`. A FASTQ or
sequence-only BAM cannot provide modification measurements. SETUN CpG-only
5mC/5hmC input has no 6mA channel. Set `--cytosine-context` according to the
basecalling model. `--allow-no-modifications` intentionally permits a
sequence-only structural analysis.

## No retained intervals

Check `loci.tsv`, `junction_calls.tsv.gz`, and the input exclusions in
`read_classifications.tsv.gz`. The preliminary screen requires two distinct
discovery reads, while `elements.tsv` requires three validated P reads and
three validated L or R reads. A candidate can fail final filters even when it
was discovered. Empty output tables are valid outputs.

## No p-value in a comparison row

Only `strand=combined` and `test_status=tested` identify completed tests.
Orientation rows are descriptive. `missing_state`, `no_shared_positions` and
`no_eligible_orientation` explain unavailable comparisons. An
`exact_state_limit` status requires a larger exact state budget and a new
output directory, because changing scientific parameters invalidates resume.
Do not replace missing p-values with one.

## Existing output directory or resume mismatch

Use a new output directory for a new analysis. For an interrupted run, repeat
the same command with `--resume`. Do not move or change inputs, update packages,
or alter scientific parameters if you want to reuse its checkpoints.

## Memory, CPU use or slow execution

Try `--threads 2 --validation-workers 1` for a smaller worker budget. Validation
workers hold reference and index data, so increasing their number can increase
memory consumption. Mapping and BAM conversion overlap. Read count, number of
candidates and alignment size affect runtime and storage needs. Keep enough
free space for intermediate FASTA and BAM files in the output directory.

## Crash without a Python exception

Preserve logs and checkpoints, then enable fault traces for the same command:

```bash
PYTHONFAULTHANDLER=1 metaphaser run ... --resume
```

Replace `...` with the original arguments. A SIGSEGV is a process crash, not
evidence of an SSH disconnect. The known validation incident is described in
[VALIDATION.md](VALIDATION.md). Report reproducible failures through
[GitHub Issues](https://github.com/andreyeremin/metaphaser/issues).
