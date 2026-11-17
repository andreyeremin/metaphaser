"""Recheck the manuscript endpoints from author-supplied saved statistical inputs."""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import math
import platform
import tempfile
from pathlib import Path
from types import SimpleNamespace

from metaphaser.config import runtime_info
from metaphaser.modifications import compare
from metaphaser.multiple_testing import adjust_primary
from metaphaser.stratified_exact import stratified, controls


def read_table(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary-table", required=True)
    parser.add_argument("--oi-fixed-inputs", required=True)
    parser.add_argument("--setun-strata", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--expected-tests", type=int, default=383)
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=False)
    reference = read_table(args.primary_table)
    with gzip.open(args.oi_fixed_inputs, "rt") as stream:
        oi = json.load(stream)
    setun = collections.defaultdict(list)
    for row in read_table(args.setun_strata):
        if row["sample"].startswith("SETUN"):
            key = (row["sample"].replace("_barcode", ""), row["locus"], row["contrast"], row["modification"])
            setun[key].append(dict(strand=row["strand"], baseline_counts=json.loads(row["baseline_counts"]),
                                   state_counts=json.loads(row["state_counts"]), shared_positions=int(row["shared_positions"])))
    recalculated = []
    counts = collections.Counter()
    with tempfile.TemporaryDirectory(dir=output) as temporary:
        for row in reference:
            key = tuple(row[field] for field in ["sample", "locus", "contrast", "modification"])
            if key[0].startswith("SETUN"):
                result = stratified(setun[key])
                counts["SETUN_saved_count_exact_matches"] += 1
            else:
                name = "|".join(row[field] for field in ["sample", "locus", "contig", "contrast", "modification"])
                strata = oi[name]
                observations = {}
                state, baseline = row["contrast"].split("-")
                role = "head" if baseline == "L" else "tail"
                for stratum in strata:
                    positions = stratum["positions"]
                    for group, ids_key, matrix_key in [(baseline, "baseline_ids", "baseline_ml"), (state, "P_ids", "P_ml")]:
                        assert len(stratum[ids_key]) == len(stratum[matrix_key])
                        for read_id, scores in zip(stratum[ids_key], stratum[matrix_key]):
                            assert len(scores) == len(positions)
                            obs_key = (row["locus"], group, role, row["modification"], stratum["strand"], read_id)
                            assert obs_key not in observations
                            observations[obs_key] = dict(zip(positions, [int(score >= 205) for score in scores]))
                locus = dict(sample=row["sample"], locus=row["locus"], region="", structural_class="host_linked_P_one_boundary")
                results = compare(observations, [locus], Path(temporary), SimpleNamespace(exact_state_limit=2000000),
                                  SimpleNamespace(event=lambda *args, **kwargs: None))
                result = next(item for item in results if item["strand"] == "combined"
                              and item["contrast"] == row["contrast"] and item["modification"] == row["modification"])
                assert result["test_status"] == "tested"
                counts["OI_fixed_ML_compare_exact_matches"] += 1
            assert result["p_numerator"] == int(row["p_numerator"]), key
            assert result["p_denominator"] == int(row["p_denominator"]), key
            assert abs(result["delta_pp"] - float(row["delta_pp"])) < 1e-10, key
            recalculated.append(dict(row, **{field: result[field] for field in ["p_numerator", "p_denominator", "p_value", "delta_pp"]}))
    adjusted = adjust_primary(recalculated, args.expected_tests)
    for row in adjusted:
        assert math.isclose(row["q_BH"], float(row["q_BH_383"]), rel_tol=1e-14, abs_tol=0)
    with (output / "reconstructed_primary.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(adjusted[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(adjusted)
    significant = [{field: row[field] for field in ["sample", "contig", "contrast", "modification", "delta_pp", "p_value", "q_BH"]}
                   for row in adjusted if row["BH_significant_005"]]
    verification = dict(status="passed", primary_comparisons=len(adjusted),
                        nominal_below_005=sum(float(row["p_value"]) < .05 for row in adjusted),
                        BH_below_005=len(significant), significant=significant,
                        checks=dict(counts), controls=controls(),
                        runtime=runtime_info("unused in saved-endpoint reconstruction"),
                        input_sha256={Path(path).name: hashlib.sha256(Path(path).read_bytes()).hexdigest()
                                      for path in [args.primary_table, args.oi_fixed_inputs, args.setun_strata]})
    (output / "verification.json").write_text(json.dumps(verification, indent=2) + "\n")
    print(json.dumps(verification, indent=2))


if __name__ == "__main__":
    main()
