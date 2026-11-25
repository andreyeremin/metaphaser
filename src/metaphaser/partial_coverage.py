"""Exploratory, bounded partial-coverage sensitivity from saved site calls.

For each channel and orientation, select sites with >=6 pooled observed reads,
without consulting state labels. At selected site j, mu_j is the pooled binary
call mean. Read i has score r_i = mean_j(x_ij - mu_j) over its observed selected
sites. The effect is 100 * mean_orientation(mean_P(r) - mean_baseline(r)).
Whole reads, including their coverage masks, are permuted within orientation;
the panel, pooled means, scores and group sizes stay fixed. All arithmetic used
for selection, centring and exact tail membership is integer/Fraction arithmetic.

The same-profile null requires exchangeability of whole observed profiles,
including coverage, between structural states within orientation. Centring does
not repair state-dependent missingness. A zero-score or same-profile p=1 is an
algebraic result, not evidence for coverage exchangeability or no biological
effect. No automatic recommendation follows from total contig depth.

Only reads represented in the input can be accounted for. Represented reads
with no usable selected observations are retained and block their comparison;
an entirely absent read cannot be recovered from site_calls.tsv.gz alone.
"""
from __future__ import annotations

import argparse
import bisect
from collections import Counter
import csv
from fractions import Fraction
import gzip
import json
import math
import os
from pathlib import Path
import sys
import tempfile

from .stratified_exact import ExactLimit

CHANNELS = ("5mC", "5hmC", "6mA")
ROLES = {"head": "L", "tail": "R"}
DEFAULT_STATE_LIMIT = 200_000
REQUIRED = {"sample", "locus", "read_id", "junction", "role", "modification",
            "strand", "contig", "position", "call_status", "high_call"}
NULL_CAVEAT = (
    "The same-profile null requires whole observed profiles (calls and coverage "
    "masks) to be exchangeable between structural states within orientation. "
    "Site centring and pooled label-invariant selection do not correct "
    "group-dependent missingness. Same-profile or zero-effect p=1 does not "
    "validate exchangeability or establish biological equivalence."
)


def _positive_integer(value, name, minimum=1):
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def _mean(values):
    return sum(values, Fraction(0)) / len(values)


def _choose_distribution(values, k, limit):
    """Count distinct read assignments, grouping repeated exact scores.

    Check the bound before inserting a new state, rather than after building an
    unbounded layer. Feasibility pruning keeps only paths that can reach k.
    """
    dp = {(0, Fraction(0)): 1}
    remaining = len(values)
    for value, frequency in sorted(Counter(values).items()):
        remaining -= frequency
        nxt = {}
        for (used, score), weight in dp.items():
            lo, hi = max(0, k - used - remaining), min(frequency, k - used)
            for take in range(lo, hi + 1):
                key = (used + take, score + take * value)
                if key not in nxt:
                    if len(nxt) >= limit:
                        raise ExactLimit(
                            f"Partial-coverage exact subset-sum needs >{limit} states; "
                            "no p-value was substituted. Increase --exact-state-limit "
                            "only with an explicit resource budget."
                        )
                    nxt[key] = 0
                nxt[key] += weight * math.comb(frequency, take)
        dp = nxt
    result = {score: weight for (used, score), weight in dp.items() if used == k}
    if sum(result.values()) != math.comb(len(values), k):
        raise RuntimeError("Exact subset distribution does not normalize")
    return result


def exact_residual_test(strata, state_limit=DEFAULT_STATE_LIMIT):
    """Exact two-sided, equally weighted orientation test of rational scores.

    Each of one or two strata contains baseline_residuals and state_residuals.
    Floats are refused: callers must supply integers or Fractions. Equal scores
    can share a DP state but their individual-read binomial weights are retained.
    p_numerator/p_denominator are unreduced assignment counts, not rounded p.
    """
    _positive_integer(state_limit, "exact_state_limit")
    strata = list(strata)
    if not 1 <= len(strata) <= 2:
        raise ValueError("Expected one or two separate orientations")
    groups, observed = [], []
    for stratum in strata:
        a, b = [list(stratum[key]) for key in ("baseline_residuals", "state_residuals")]
        if not a or not b:
            raise ValueError("Both read groups must be nonempty in each orientation")
        if any(not isinstance(v, (int, Fraction)) or isinstance(v, bool) for v in a + b):
            raise ValueError("Exact residuals must be integers or Fraction values, not floats")
        a, b = list(map(Fraction, a)), list(map(Fraction, b))
        groups.append((a, b))
        observed.append(_mean(b) - _mean(a))
    effect = _mean(observed)
    total_space = math.prod(math.comb(len(a) + len(b), len(a)) for a, b in groups)
    result = dict(delta_pp=float(100 * effect), delta_pp_exact=100 * effect,
                  p_denominator=total_space)
    # With cutoff zero every assignment is extreme. This is an exact identity,
    # never a fallback from a failed or oversized DP calculation.
    if effect == 0:
        return dict(result, p_value=1.0, p_exact=Fraction(1), p_numerator=total_space)
    distributions = []
    for a, b in groups:
        na, nb, total = len(a), len(b), sum(a + b, Fraction(0))
        take_baseline = na <= nb
        raw = _choose_distribution(a + b, min(na, nb), state_limit)
        distributions.append({
            ((na * total - (na + nb) * v) if take_baseline else
             ((na + nb) * v - nb * total)) / (na * nb): weight
            for v, weight in raw.items()
        })
    cutoff = abs(sum(observed, Fraction(0)))
    if len(distributions) == 1:
        extreme = sum(w for value, w in distributions[0].items() if abs(value) >= cutoff)
    else:
        # Count the strictly interior interval with exact bisect boundaries;
        # do not build the Cartesian product of the orientation distributions.
        values = sorted(distributions[1])
        prefix = [0]
        for value in values:
            prefix.append(prefix[-1] + distributions[1][value])
        interior = 0
        for value, weight in distributions[0].items():
            left = bisect.bisect_right(values, -cutoff - value)
            right = bisect.bisect_left(values, cutoff - value)
            interior += weight * (prefix[right] - prefix[left])
        extreme = total_space - interior
    if not 0 < extreme <= total_space:
        raise RuntimeError("Exact tail count is inconsistent")
    probability = Fraction(extreme, total_space)
    return dict(result, p_value=float(probability), p_exact=probability, p_numerator=extreme)


def read_site_calls(path, *, sample=None, locus=None, modification=None, contrast=None):
    """Load relevant L/P head and R/P tail observations and preserve read rosters.

    Duplicate read-position rows must agree in call status/value and any saved
    score/query-position fields. Unknown calls register a read, never a zero.
    A read cannot belong to both states/orientations within the same role.
    """
    if modification is not None and modification not in CHANNELS:
        raise ValueError("Use a separate modification channel: 5mC, 5hmC or 6mA")
    if contrast not in (None, "P-L", "P-R"):
        raise ValueError("Use contrast P-L or P-R")
    comparisons, assignments, duplicates = {}, {}, {}
    counts = Counter(input_rows=0, selected_rows=0, duplicate_rows=0, unknown_observations=0)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not REQUIRED.issubset(reader.fieldnames or []):
            raise ValueError("Input must be a site_calls.tsv[.gz] table; missing columns: "
                             + ", ".join(sorted(REQUIRED - set(reader.fieldnames or []))))
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError("Duplicate input column names")
        for line, row in enumerate(reader, 2):
            counts["input_rows"] += 1
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Malformed site-call row at line {line}")
            if any(wanted is not None and row[key] != wanted for key, wanted in
                   (("sample", sample), ("locus", locus), ("modification", modification))):
                continue
            baseline = ROLES.get(row["role"])
            if baseline is None or row["junction"] not in (baseline, "P"):
                continue
            if contrast is not None and contrast != "P-" + baseline:
                continue
            if row["modification"] not in CHANNELS or row["strand"] not in ("+", "-"):
                raise ValueError(f"Invalid separate channel or orientation at line {line}")
            if any(not row[key].strip() for key in ("sample", "locus", "contig", "read_id")):
                raise ValueError(f"Empty sample/locus/contig/read identifier at line {line}")
            try:
                position = int(row["position"])
            except ValueError as exc:
                raise ValueError(f"Invalid position at line {line}") from exc
            if position < 0:
                raise ValueError(f"Negative position at line {line}")
            status, high = row["call_status"], row["high_call"]
            if status == "unknown" and high == "":
                value = None
            elif status in ("explicit", "implicit_low") and high in ("0", "1"):
                value = int(high)
                if status == "implicit_low" and value != 0:
                    raise ValueError(f"Conflicting implicit_low call at line {line}")
            else:
                raise ValueError(f"Invalid call_status/high_call pair at line {line}")
            key = (row["sample"], row["locus"], row["role"], row["modification"])
            read_key = (row["sample"], row["locus"], row["role"], row["read_id"])
            assignment = (row["junction"], row["strand"], row["contig"])
            if assignments.setdefault(read_key, assignment) != assignment:
                raise ValueError(f"Conflicting read group/orientation/contig at line {line}: {row['read_id']}")
            comparison = comparisons.setdefault(key, {"contig": row["contig"], "strata": {}})
            if comparison["contig"] != row["contig"]:
                raise ValueError(f"Conflicting contigs for one locus/role/channel at line {line}")
            groups = comparison["strata"].setdefault(row["strand"], {baseline: {}, "P": {}})
            observations = groups[row["junction"]].setdefault(row["read_id"], {})
            duplicate_key = (key, row["read_id"], position)
            signature = (status, value, *(row.get(k, "") for k in ("query_position", "stored_ml", "score_used")))
            counts["selected_rows"] += 1
            if duplicate_key in duplicates:
                if duplicates[duplicate_key] != signature:
                    raise ValueError(f"Conflicting duplicate read-position at line {line}: {row['read_id']}:{position}")
                counts["duplicate_rows"] += 1
                continue
            duplicates[duplicate_key] = signature
            if value is None:
                counts["unknown_observations"] += 1
            else:
                observations[position] = value
    counts["unique_read_position_rows"] = len(duplicates)
    counts["ignored_rows"] = counts["input_rows"] - counts["selected_rows"]
    return comparisons, dict(counts)


def prepare_orientation(groups, baseline, strand, min_site_reads=6):
    """Compute fixed panels/scores and descriptive coverage diagnostics, no test."""
    _positive_integer(min_site_reads, "min_site_reads", 6)
    if baseline not in ("L", "R") or strand not in ("+", "-"):
        raise ValueError("Expected baseline L/R and orientation +/-")
    if set(groups) != {baseline, "P"} or set(groups[baseline]) & set(groups["P"]):
        raise ValueError("Read groups must be separate and preserved")
    pooled, highs = Counter(), Counter()
    covered = {state: Counter() for state in groups}
    for state, reads in groups.items():
        for values in reads.values():
            for position, value in values.items():
                if value not in (0, 1) or not isinstance(value, int):
                    raise ValueError("Observed high calls must be binary integers")
                pooled[position] += 1
                highs[position] += value
                covered[state][position] += 1
    positions = sorted(p for p, n in pooled.items() if n >= min_site_reads)
    means = {p: Fraction(highs[p], pooled[p]) for p in positions}
    selected = set(positions)
    reads, scores, masks = [], {baseline: [], "P": []}, {baseline: Counter(), "P": Counter()}
    coverage = {}
    zero_reads = []
    for state in (baseline, "P"):
        sizes = []
        for rid, values in sorted(groups[state].items()):
            observed = sorted(selected.intersection(values))
            residual = _mean([Fraction(values[p]) - means[p] for p in observed]) if observed else None
            reads.append(dict(read_id=rid, state=state, observed_sites=len(values),
                              selected_observed_sites=len(observed), selected_positions=observed,
                              coverage_fraction=Fraction(len(observed), len(positions)) if positions else None,
                              residual=residual))
            sizes.append(len(observed))
            masks[state][tuple(observed)] += 1
            if residual is None:
                zero_reads.append(dict(read_id=rid, state=state))
            else:
                scores[state].append(residual)
        coverage[state] = dict(n_reads=len(sizes), zero_observation_reads=sizes.count(0),
                               min_selected_observations=min(sizes) if sizes else None,
                               max_selected_observations=max(sizes) if sizes else None,
                               mean_selected_observations=_mean(sizes) if sizes else None)
    na, nb = len(groups[baseline]), len(groups["P"])
    both = {p for p in positions if covered[baseline][p] and covered["P"][p]}
    complete = {p for p in positions if pooled[p] == na + nb}
    variable = {p for p in positions if 0 < highs[p] < pooled[p]}
    imbalanced = [p for p in positions if na and nb and
                  Fraction(covered[baseline][p], na) != Fraction(covered["P"][p], nb)]
    mask_tv = None
    warnings = []
    if na and nb:
        mask_tv = sum((abs(Fraction(masks[baseline][mask], na) - Fraction(masks["P"][mask], nb))
                       for mask in masks[baseline].keys() | masks["P"].keys()), Fraction(0)) / 2
        if mask_tv:
            warnings.append("group_dependent_missingness: observed coverage-mask distributions differ "
                            "between groups; descriptive warning, not a formal missingness test.")
        if positions and not both:
            warnings.append("no_between_group_site_overlap: centring cannot identify a state difference "
                            "at sites observed in only one group.")
    if zero_reads:
        warnings.append("zero_observation_reads: represented reads with no selected observations were "
                        "retained; the comparison is untestable, with no read or orientation dropped.")
    if not na or not nb:
        status = "missing_state"
    elif not positions:
        status = "no_selected_sites"
    elif zero_reads:
        status = "zero_observation_reads"
    else:
        status = "ready"
    site_rows = [dict(position=p, pooled_reads=pooled[p], pooled_mean=means[p],
                      baseline_reads=covered[baseline][p], state_reads=covered["P"][p],
                      variable=p in variable, observed_in_both_groups=p in both) for p in positions]
    effect = _mean(scores["P"]) - _mean(scores[baseline]) if status == "ready" else None
    return dict(strand=strand, test_status=status, n_baseline=na, n_state=nb,
                candidate_sites=len(pooled), selected_sites=len(positions), informative_sites=len(variable),
                informative_sites_with_group_overlap=len(variable & both), positions=positions, sites=site_rows,
                overlap=dict(sites_observed_in_both_groups=len(both), sites_observed_in_every_read=len(complete),
                             baseline_only_sites=sum(covered["P"][p] == 0 for p in positions),
                             state_only_sites=sum(covered[baseline][p] == 0 for p in positions),
                             union_jaccard=Fraction(len(both), len(positions)) if positions else None,
                             sites_with_unequal_group_coverage=len(imbalanced), coverage_mask_total_variation=mask_tv),
                read_coverage=coverage, reads=reads, zero_observation_reads=zero_reads, warnings=warnings,
                baseline_residuals=scores[baseline], state_residuals=scores["P"],
                delta_pp_exact=100 * effect if effect is not None else None)


def analyze_site_calls(path, *, min_site_reads=6, state_limit=DEFAULT_STATE_LIMIT, **filters):
    """Return an exploratory report; every present orientation must be usable."""
    _positive_integer(min_site_reads, "min_site_reads", 6)
    _positive_integer(state_limit, "exact_state_limit")
    comparisons, input_counts = read_site_calls(path, **filters)
    if not comparisons:
        raise ValueError("No matching L/P head or R/P tail site-call rows")
    results = []
    for (sample, locus, role, modification), data in sorted(comparisons.items()):
        baseline = ROLES[role]
        parts = [prepare_orientation(data["strata"][strand], baseline, strand, min_site_reads)
                 for strand in ("+", "-") if strand in data["strata"]]
        result = dict(sample=sample, locus=locus, contig=data["contig"], contrast="P-" + baseline,
                      role=role, modification=modification, exploratory=True,
                      test="partial_coverage_exact_within_orientation_equal_weights",
                      orientations=[part["strand"] for part in parts], strata=parts,
                      test_status="untestable", delta_pp=None, delta_pp_exact=None,
                      p_value=None, p_exact=None, p_numerator=None, p_denominator=None,
                      warnings=[f"{p['strand']}: {w}" for p in parts for w in p["warnings"]])
        if len(parts) == 1:
            result["warnings"].append("single_orientation: only one orientation is represented in the input.")
        if all(part["test_status"] == "ready" for part in parts):
            result["delta_pp_exact"] = _mean([p["delta_pp_exact"] for p in parts])
            result["delta_pp"] = float(result["delta_pp_exact"])
            try:
                result.update(exact_residual_test(parts, state_limit))
                result["test_status"] = "tested"
            except ExactLimit as exc:
                result["test_status"] = "exact_state_limit"
                result["error"] = str(exc)
        else:
            result["error"] = "No test: " + "; ".join(
                f"{p['strand']}: {p['test_status']}" for p in parts if p["test_status"] != "ready")
        results.append(result)
    return dict(schema_version=1, mode="exploratory_partial_coverage", exploratory=True,
                status="complete" if all(r["test_status"] == "tested" for r in results) else "incomplete",
                source=str(Path(path).resolve()), input_counts=input_counts,
                settings=dict(min_site_reads=min_site_reads, exact_state_limit=state_limit, filters=filters),
                selection="Fixed pooled observed-read count per site/channel/orientation; label-invariant.",
                effect="P minus baseline, in percentage points, with equal weights across represented orientations.",
                null_caveat=NULL_CAVEAT,
                recommendation="Exploratory sensitivity only. No automatic recommendation from total contig depth.",
                multiplicity="Nominal exploratory p-values; separate from the primary method and its BH family.",
                read_roster_scope="All matching reads represented in site_calls are retained, including unknown-only "
                                  "reads. Reads entirely absent from this table cannot be counted or recovered.",
                diagnostic_scope="Coverage imbalance/overlap warnings are descriptive, not a missingness test. "
                                 "Absence of warnings does not establish exchangeability.",
                site_metric_definitions="Candidate sites have at least one observed call; selected sites meet "
                                        "the pooled read minimum. Informative sites are selected sites with "
                                        "both 0 and 1 calls. Overlap and coverage-mask diagnostics use selected sites.",
                fraction_encoding="Rational quantities are stored as exact numerator/denominator strings in JSON.",
                comparisons=results)


def run_partial_coverage(site_calls, output, **options):
    """Write a complete diagnostic report atomically, without overwriting files."""
    output = Path(output)
    if output.exists():
        raise ValueError("Partial-coverage output already exists; choose a new filename")
    report = analyze_site_calls(site_calls, **options)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", dir=output.parent, delete=False, encoding="utf-8") as stream:
            temporary = Path(stream.name)
            json.dump(report, stream, indent=2, default=_json_fraction, allow_nan=False)
            stream.write("\n")
        os.link(temporary, output)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return dict(status=report["status"], exploratory=True, output=str(output.resolve()),
                comparisons=len(report["comparisons"]),
                tested=sum(r["test_status"] == "tested" for r in report["comparisons"]),
                errors=[dict(sample=r["sample"], locus=r["locus"], contrast=r["contrast"],
                             modification=r["modification"], error=r["error"])
                        for r in report["comparisons"] if "error" in r])


def _json_fraction(value):
    if isinstance(value, Fraction):
        return f"{value.numerator}/{value.denominator}"
    raise TypeError(f"Cannot serialize {type(value).__name__}")


def add_arguments(parser):
    """Shared by the standalone entry point and the separately supplied CLI patch."""
    def positive(value):
        try:
            number = int(value)
            _positive_integer(number, "value")
            return number
        except ValueError as exc:
            raise argparse.ArgumentTypeError(str(exc)) from exc

    def site_minimum(value):
        number = positive(value)
        if number < 6:
            raise argparse.ArgumentTypeError("pooled site-read minimum must be >= 6")
        return number

    parser.add_argument("--site-calls", required=True, help="Ready site_calls.tsv[.gz]; no BAM or GPU")
    parser.add_argument("--output", required=True, help="New exploratory JSON report, including coverage diagnostics")
    parser.add_argument("--min-site-reads", type=site_minimum, default=6, help="Pooled observed reads per orientation/site (>=6)")
    parser.add_argument("--exact-state-limit", type=positive, default=DEFAULT_STATE_LIMIT,
                        help="Maximum states per DP layer; fail explicitly above this bound (default: %(default)s)")
    parser.add_argument("--sample", help="Restrict to one saved sample")
    parser.add_argument("--locus", help="Restrict to one saved locus")
    parser.add_argument("--modification", choices=CHANNELS, help="Restrict to one separate channel")
    parser.add_argument("--contrast", choices=("P-L", "P-R"), help="L vs P at head, or R vs P at tail")


def run_from_args(args):
    result = run_partial_coverage(args.site_calls, args.output, min_site_reads=args.min_site_reads,
                                 state_limit=args.exact_state_limit, sample=args.sample, locus=args.locus,
                                 modification=args.modification, contrast=args.contrast)
    print(json.dumps(result))
    return 0 if result["status"] == "complete" else 2


def main(argv=None):
    parser = argparse.ArgumentParser(description="Exploratory partial-coverage sensitivity; no automatic depth recommendation.")
    add_arguments(parser)
    args = parser.parse_args(argv)
    try:
        return run_from_args(args)
    except (ValueError, OSError, EOFError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
