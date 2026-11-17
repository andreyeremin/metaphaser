import csv
import json
import subprocess
import sys
from types import SimpleNamespace

import pytest

from metaphaser.modifications import compare
from metaphaser.multiple_testing import adjust_primary, adjust_files
from metaphaser.stratified_exact import stratified, controls, ExactLimit


def test_existing_stratified_controls():
    assert controls()["status"] == "passed"


def test_exact_limit_never_substitutes_a_sampled_p_value():
    with pytest.raises(ExactLimit):
        stratified([dict(baseline_counts=[0, 1, 2], state_counts=[3, 4, 5], shared_positions=5)], limit=1)


def test_comparison_keeps_one_read_and_one_site_and_equal_orientation_weights(tmp_path):
    observations = {}
    for strand, baseline, alternative in [("+", [1], [0]), ("-", [0, 1], [1])]:
        for state, values in [("L", baseline), ("P", alternative)]:
            for index, high in enumerate(values):
                observations[("J1", state, "head", "5mC", strand, f"{state}{strand}{index}")] = {10: high}
    loci = [dict(sample="synthetic", locus="J1", region="region1", structural_class="host_linked_P_one_boundary")]
    rows = compare(observations, loci, tmp_path, SimpleNamespace(exact_state_limit=2000000),
                   SimpleNamespace(event=lambda *args, **kwargs: None))
    combined = next(row for row in rows if row["contrast"] == "P-L" and row["strand"] == "combined")
    expected = stratified([dict(baseline_counts=[1], state_counts=[0], shared_positions=1),
                           dict(baseline_counts=[0, 1], state_counts=[1], shared_positions=1)])
    assert combined["test_status"] == "tested"
    assert combined["delta_pp"] == -25
    assert combined["p_numerator"] == expected["p_numerator"]
    assert combined["p_denominator"] == 6
    assert combined["orientations"] == "+,-"


def primary_row(locus="J1", numerator=1, denominator=100):
    return dict(sample="synthetic", locus=locus, contrast="P-R", modification="5mC",
                strand="combined", structural_class="host_linked_P_one_boundary", test_status="tested",
                test="exact_within_orientation_equal_weights", p_numerator=numerator,
                p_denominator=denominator, p_value=numerator / denominator)


def test_bh_expected_family_duplicates_and_exact_probabilities():
    rows = [primary_row(), primary_row("J2", 1, 10)]
    adjusted = adjust_primary(rows, expected_tests=2)
    assert [row["q_BH"] for row in adjusted] == pytest.approx([.02, .1])
    with pytest.raises(ValueError, match="Expected"):
        adjust_primary(rows, expected_tests=383)
    with pytest.raises(ValueError, match="Duplicate"):
        adjust_primary([rows[0], rows[0]])
    with pytest.raises(ValueError, match="disagrees"):
        adjust_primary([dict(rows[0], p_value=.5)])


def test_bh_cli_excludes_orientation_and_B_rows_and_rejects_incomplete_exact_test(tmp_path):
    source = tmp_path / "comparisons.tsv"
    rows = [primary_row(), dict(primary_row("J2"), strand="+"),
            dict(primary_row("J3"), contrast="B-R")]
    with source.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    output = tmp_path / "primary.tsv"
    result = subprocess.run([sys.executable, "-m", "metaphaser", "adjust", "--comparisons", str(source),
                             "--output", str(output), "--expected-tests", "1"], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["tests"] == 1
    with pytest.raises(ValueError, match="already exists"):
        adjust_files([source], output)
    text = source.read_text().replace("\ttested\t", "\texact_state_limit\t")
    source.write_text(text)
    with pytest.raises(ValueError, match="incomplete"):
        adjust_files([source], tmp_path / "incomplete.tsv")
