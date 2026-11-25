# Historical manuscript settings and saved-table reproduction

The public software is [andreyeremin/metaphaser](https://github.com/andreyeremin/metaphaser),
under the [MIT licence](../LICENSE). The frozen primary study family is supplied
in [study/383](../study/383/README.md), with a separate data-licensing notice.
Use the exact commands there to recompute all 383 tests from the saved profiles.
The source snapshot under `v0.2.0rc2` predates this directory; use the `v0.2.0rc3` release, which includes this directory. No DOI is assigned.

## Historical versus release behaviour

| Setting | Historical manuscript analysis | rc2 default |
|---|---|---|
| Preliminary distinct discovery reads | 3 | 2 (`--min-discovery-reads`) |
| Observed boundary pairs per group | Up to 5 | Up to 5 |
| Additional positive-gap refined boundary pairs | Not used | Up to 5 |
| Final retained host-linked interval | ≥3 validated P reads and ≥3 validated L or R reads | Same |
| ML high-call threshold | ≥205 | ≥205 (`--ml-cutoff`) |
| Cytosine context | OI: all; SETUN: CpG | CpG (`--cytosine-context`) |
| Modification channels | OI: 5mC, 5hmC, 6mA; SETUN: 5mC, 5hmC | Available supported channels (`--modifications auto`) |
| Statistical endpoint | Fixed shared-position per-read fractions, equal orientation weights | Same exact implementation |
| Study-wide BH | One frozen family of 383 primary comparisons | Family determined by all supplied sample tables |

`--min-discovery-reads 3` restores the numerical preliminary threshold but
leaves rc2 boundary refinement active. There is no documented CLI switch to
disable that refinement. The parameter validator rejects
`refined_breakpoint_hypotheses: 0`. Consequently, neither that flag alone nor a
JSON file of historical numerical settings makes rc2 an exact historical
discovery replay. See [the rc2 discovery changes](DISCOVERY_FIX_20261117.md).

[historical_settings.json](../study/383/historical_settings.json) records the
historical structural parameters, mapping arguments, modification contexts,
statistical endpoint and family. The structural values agree between the
archived rc1 configuration and the historical SETUN configuration. Historical
discovery behaviour and sample-specific contexts are cross-checked against the
manuscript snapshot. SHA-256 hashes identify those source records without
publishing private filesystem paths or input manifests.

The shared structural settings include interval lengths 200–1,000,000 bp,
75-bp anchors, 85% discovery identity, at most 250-bp read gap/overlap,
25/200-bp boundary tolerances for intervals shorter than/at least 2000 bp,
2000-bp external flanks, up to 250 bp per side in the primary validation
window, 90% paired-base coverage, 95% acceptance identity, a maximum central
indel of 10 bp, and an edit advantage of max(3, ceil(0.02 × query length)).
Modification windows extend up to 1000 bp, limited by available internal arms.
The 500/1000-bp validation windows described in the manuscript are descriptive
extensions, not extra primary assignment requirements.

## What is reproduced

The full count and reference tables are frozen by SHA-256. Exact tests operate
on the saved integer counts without changing assignments, positions or
orientation membership. BH includes all 383 tests across the three separately
analysed modification channels. The existing exact-test and BH source modules
are byte-identical between the archived rc1 source and the public rc2 checkout.
The verification record identifies the actual imported source modules and
runtime, distinguishing this calculation from a fresh discovery run.

The historical receipt records Python 3.12.14, metaphaser 0.2.0rc1, NumPy
1.26.4 and SciPy 1.16.1. It remains separately named
`verification-20261113.json`. The current local receipt is `verification.json`.
Neither receipt establishes the runtime versions of upstream basecalling,
assembly or all original biological analyses.

## Partial coverage

Release v0.2.0rc3 adds an optional `partial-coverage` command operating on
saved site-call tables. See [its method and input contract](PARTIAL_COVERAGE.md).
This endpoint is exploratory and separate from the frozen 383 primary tests.
Its coverage diagnostics do not establish exchangeability or correct
state-dependent missingness. Primary discovery and shared-position comparison
code are unchanged from rc2.
