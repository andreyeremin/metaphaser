"""Recompute exact p and BH from study-derived molecule-level integer counts."""
import argparse
import collections
import csv
import hashlib
import importlib.metadata
import json
import math
import platform
from pathlib import Path

from metaphaser import __version__
from metaphaser import multiple_testing, stratified_exact
from metaphaser.stratified_exact import stratified
from metaphaser.multiple_testing import adjust_primary, KEY_FIELDS

EXPECTED_TESTS = 383
EXPECTED_PROFILES = 9793
EXPECTED_STRATA = 661
INPUT_SHA256 = {
    "counts": "f1b981433d5207a5dbd2828e4182153e687cb69ea1bfa6cea918adb6f70a6c3a",
    "reference": "6d4f211dbcae678faded33275f850f175982a9f452ab4f7cb82a657c04665f87",
}


def verify_inputs(counts_path, reference_path):
    """Require the complete frozen family, including nonsignificant tests."""
    for role, path in [("counts", counts_path), ("reference", reference_path)]:
        if hashlib.sha256(Path(path).read_bytes()).hexdigest() != INPUT_SHA256[role]:
            raise ValueError(f"The {role} file differs from the frozen 383-test study input.")


def read_table(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    output = Path(args.output)
    if output.exists():
        raise ValueError("Output already exists. Choose a fresh output directory.")
    verify_inputs(args.counts, args.reference)
    counts = read_table(args.counts)
    reference = read_table(args.reference)
    if len(reference) != EXPECTED_TESTS or len(counts) != EXPECTED_PROFILES:
        raise ValueError("The frozen study requires 383 comparisons and 9793 count profiles.")
    reference_by_key = {tuple(row[field] for field in KEY_FIELDS): row for row in reference}
    if len(reference_by_key) != len(reference):
        raise ValueError("Duplicate primary comparison key.")
    groups = collections.defaultdict(lambda: collections.defaultdict(list))
    seen = set()
    for row in counts:
        key = tuple(row[field] for field in KEY_FIELDS)
        strand, state = row["orientation"], row["state"]
        profile_key = key + (strand, state, int(row["profile_index"]))
        if profile_key in seen:
            raise ValueError("Duplicate count profile.")
        seen.add(profile_key)
        if key not in reference_by_key or strand not in {"+", "-"}:
            raise ValueError("Unrecognized comparison or orientation.")
        high, sites = int(row["high_calls"]), int(row["shared_positions"])
        if sites < 1 or not 0 <= high <= sites:
            raise ValueError("Invalid high-call count or denominator.")
        groups[key][strand].append((state, high, sites))
    if set(groups) != set(reference_by_key):
        raise ValueError("The count table does not cover the declared primary family.")
    if sum(len(strands) for strands in groups.values()) != EXPECTED_STRATA:
        raise ValueError("The frozen study requires 661 orientation strata.")
    results = []
    for key, row in reference_by_key.items():
        alternative, baseline = row["contrast"].split("-")
        strata = []
        observed_orientations = []
        for strand in ["+", "-"]:
            if strand not in groups[key]:
                continue
            observations = groups[key][strand]
            sites = {item[2] for item in observations}
            states = {item[0] for item in observations}
            if len(sites) != 1 or states != {baseline, alternative}:
                raise ValueError("Each orientation must use one shared denominator and both states.")
            strata.append(dict(baseline_counts=[high for state, high, _ in observations if state == baseline],
                               state_counts=[high for state, high, _ in observations if state == alternative],
                               shared_positions=sites.pop()))
            observed_orientations.append(strand)
        if set(observed_orientations) != set(row["orientations"].split(",")):
            raise ValueError("Orientation membership differs from the reference.")
        computed = stratified(strata)
        if computed["p_numerator"] != int(row["p_numerator"]) or computed["p_denominator"] != int(row["p_denominator"]):
            raise ValueError(f"Exact probability differs for {key}.")
        if abs(computed["delta_pp"] - float(row["delta_pp"])) > 1e-10:
            raise ValueError(f"Effect differs for {key}.")
        if not math.isclose(computed["p_value"], float(row["p_value"]), rel_tol=1e-14, abs_tol=0):
            raise ValueError(f"Stored p-value differs for {key}.")
        results.append(dict(row, **computed))
    adjusted = adjust_primary(results, expected_tests=EXPECTED_TESTS)
    for row in adjusted:
        if not math.isclose(row["q_BH"], float(row["q_BH_383"]), rel_tol=1e-14, abs_tol=0):
            raise ValueError("BH value differs from the reference.")
    output.mkdir(parents=True, exist_ok=False)
    with (output / "recomputed_primary.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(adjusted[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(adjusted)
    summary = dict(status="passed", primary_comparisons=len(adjusted), count_profiles=len(counts),
                   orientation_strata=EXPECTED_STRATA,
                   exact_numerators_matched=len(adjusted), exact_denominators_matched=len(adjusted),
                   effects_matched=len(adjusted), p_values_matched=len(adjusted), BH_values_matched=len(adjusted),
                   nominal_below_005=sum(float(row["p_value"]) < .05 for row in adjusted),
                   BH_below_005=sum(row["BH_significant_005"] for row in adjusted),
                   python=platform.python_version(),
                   packages={"metaphaser": __version__, **{name: importlib.metadata.version(name)
                             for name in ["numpy", "scipy"]}},
                   statistical_source_sha256={Path(module.__file__).name:
                       hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest()
                       for module in [stratified_exact, multiple_testing]},
                   input_sha256={Path(path).name: hashlib.sha256(Path(path).read_bytes()).hexdigest()
                                 for path in [args.counts, args.reference]})
    (output / "verification.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
