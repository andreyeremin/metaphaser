"""Small CPU-only controls; never open BAMs or run study-size permutations."""
import csv
from fractions import Fraction as F
import gzip
import itertools
import json
import math
import os
from pathlib import Path
import random
import subprocess
import sys

import pytest

from metaphaser import partial_coverage as pc
from metaphaser.stratified_exact import ExactLimit, stratified


def rows_for(matrix, *, n_baseline=3, baseline="L", strand="+", modification="5mC", role=None):
    rows = []
    for index, values in enumerate(matrix):
        for pos, value in enumerate(values):
            rows.append(dict(sample="fixture", locus="J1", read_id=f"{strand}r{index}",
                             junction=baseline if index < n_baseline else "P",
                             role=role or ("head" if baseline == "L" else "tail"),
                             modification=modification, strand=strand, contig="ctg", position=10 + pos,
                             call_status="unknown" if value is None else "explicit",
                             high_call="" if value is None else value))
    return rows


def save_rows(tmp_path, rows, name="site_calls.tsv.gz"):
    path = tmp_path / name
    opener = gzip.open if name.endswith(".gz") else open
    with opener(path, "wt", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=sorted(pc.REQUIRED), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def groups_for(matrix, n_baseline=3, baseline="L"):
    return {state: {f"r{i}": {p: v for p, v in enumerate(values) if v is not None}
                    for i, values in enumerate(matrix) if (i < n_baseline) == (state == baseline)}
            for state in (baseline, "P")}


def brute(strata):
    observed = sum((sum(s["state_residuals"], F(0)) / len(s["state_residuals"])
                    - sum(s["baseline_residuals"], F(0)) / len(s["baseline_residuals"])
                    for s in strata), F(0))
    distributions = []
    for s in strata:
        a, b = s["baseline_residuals"], s["state_residuals"]
        values, na, nb = a + b, len(a), len(b)
        total = sum(values, F(0))
        diffs = []
        for ids in itertools.combinations(range(na + nb), na):
            selected = sum((values[i] for i in ids), F(0))
            diffs.append((total - selected) / nb - selected / na)
        distributions.append(diffs)
    extreme = total = 0
    for assignment in itertools.product(*distributions):
        extreme += abs(sum(assignment, F(0))) >= abs(observed)
        total += 1
    return F(extreme, total), 100 * observed / len(strata), extreme, total


def test_full_coverage_matches_primary_for_every_label_assignment():
    matrix = [[0, 1, 1, 0], [1, 1, 0, 0], [0, 1, 1, 1],
              [0, 0, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0]]
    fixed = None
    for ids in itertools.combinations(range(6), 3):
        groups = {state: {f"r{i}": dict(enumerate(values)) for i, values in enumerate(matrix)
                          if (i in ids) == (state == "L")} for state in ("L", "P")}
        part = pc.prepare_orientation(groups, "L", "+")
        result = pc.exact_residual_test([part])
        expected = stratified([dict(baseline_counts=[sum(matrix[i]) for i in ids],
                                   state_counts=[sum(matrix[i]) for i in range(6) if i not in ids],
                                   shared_positions=4)])
        assert result["p_numerator"] == expected["p_numerator"]
        assert result["p_denominator"] == expected["p_denominator"]
        assert result["delta_pp"] == pytest.approx(expected["delta_pp"], abs=1e-13)
        panel = (part["positions"], [s["pooled_mean"] for s in part["sites"]])
        assert fixed is None or panel == fixed
        fixed = panel


def test_partial_selection_and_scores_are_label_invariant_and_match_brute():
    matrix = [[1, None, 0], [1, 1, None], [0, 1, 0], [1, 0, 1],
              [None, 0, 0], [0, None, 1], [1, 1, 1], [0, 0, 0]]
    fixed = None
    for ids in itertools.combinations(range(8), 4):
        groups = {state: {f"r{i}": {p: v for p, v in enumerate(values) if v is not None}
                          for i, values in enumerate(matrix) if (i in ids) == (state == "L")}
                  for state in ("L", "P")}
        part = pc.prepare_orientation(groups, "L", "+")
        assert part["selected_sites"] == 3
        invariant = (part["positions"], [s["pooled_mean"] for s in part["sites"]],
                     sorted((r["read_id"], r["residual"]) for r in part["reads"]))
        assert fixed is None or invariant == fixed
        fixed = invariant
        expected_p, effect, num, den = brute([part])
        actual = pc.exact_residual_test([part])
        assert (actual["p_exact"], actual["delta_pp_exact"], actual["p_numerator"], actual["p_denominator"]) == (
            expected_p, effect, num, den)


def test_fractional_two_orientation_counts_match_exhaustive_enumeration():
    rng = random.Random(163)
    for _ in range(32):
        strata = []
        for _ in range(2):
            strata.append({key: [F(rng.randint(-3, 3), rng.choice((3, 7, 11)))
                                  for _ in range(rng.randint(1, 3))]
                           for key in ("baseline_residuals", "state_residuals")})
        expected_p, effect, num, den = brute(strata)
        actual = pc.exact_residual_test(strata)
        assert (actual["p_exact"], actual["delta_pp_exact"], actual["p_numerator"], actual["p_denominator"]) == (
            expected_p, effect, num, den)
        swapped = [dict(baseline_residuals=s["state_residuals"], state_residuals=s["baseline_residuals"])
                   for s in reversed(strata)]
        other = pc.exact_residual_test(swapped)
        assert other["p_exact"] == actual["p_exact"]
        assert other["delta_pp_exact"] == -actual["delta_pp_exact"]


def test_reverse_orientation_has_equal_weight_despite_read_imbalance(tmp_path):
    forward = rows_for([[1]] * 3 + [[0]] * 3)
    reverse = rows_for([[0]] * 4 + [[1], [1], [0], [0]], n_baseline=4, strand="-")
    result = pc.analyze_site_calls(save_rows(tmp_path, forward + reverse))["comparisons"][0]
    assert result["orientations"] == ["+", "-"]
    assert result["delta_pp_exact"] == -25
    primary = stratified([dict(baseline_counts=[1] * 3, state_counts=[0] * 3, shared_positions=1),
                          dict(baseline_counts=[0] * 4, state_counts=[1, 1, 0, 0], shared_positions=1)])
    assert result["p_numerator"] == primary["p_numerator"]
    assert result["p_denominator"] == primary["p_denominator"] == 1400
    assert result["p_exact"] == brute(result["strata"])[0]


def test_unknowns_excluded_from_pooled_selection_and_not_imputed(tmp_path):
    matrix = [[1, 1], [1, 0], [0, 1], [0, 0], [1, None], [0, None]]
    report = pc.analyze_site_calls(save_rows(tmp_path, rows_for(matrix)))
    part = report["comparisons"][0]["strata"][0]
    assert part["positions"] == [10]
    assert part["sites"][0]["pooled_mean"] == F(1, 2)
    assert part["n_baseline"] == part["n_state"] == 3
    assert report["input_counts"]["unknown_observations"] == 2
    assert all(r["selected_observed_sites"] == 1 for r in part["reads"])


@pytest.mark.parametrize("extra", [None, 1])
def test_zero_observation_read_is_retained_and_blocks_entire_comparison(tmp_path, extra):
    forward = rows_for([[1]] * 3 + [[0]] * 3)
    row = dict(forward[0], read_id="zero", position=90, high_call="" if extra is None else 1,
               call_status="unknown" if extra is None else "explicit")
    reverse = rows_for([[0]] * 3 + [[1]] * 3, strand="-")
    result = pc.analyze_site_calls(save_rows(tmp_path, forward + [row] + reverse))["comparisons"][0]
    assert result["test_status"] == "untestable"
    assert result["p_value"] is result["delta_pp"] is None
    assert result["orientations"] == ["+", "-"]
    part = result["strata"][0]
    assert part["n_baseline"] == 4
    assert len(part["reads"]) == 7
    assert part["zero_observation_reads"] == [dict(read_id="zero", state="L")]
    assert next(r for r in part["reads"] if r["read_id"] == "zero")["residual"] is None
    assert "zero_observation_reads" in result["error"]


def test_same_profile_null_with_disjoint_group_coverage_warns(tmp_path):
    matrix = [[1, None]] * 6 + [[None, 1]] * 6
    report = pc.analyze_site_calls(save_rows(tmp_path, rows_for(matrix, n_baseline=6)))
    result = report["comparisons"][0]
    assert result["delta_pp_exact"] == 0 and result["p_exact"] == 1
    assert result["p_numerator"] == result["p_denominator"] == math.comb(12, 6)
    part = result["strata"][0]
    assert part["overlap"]["sites_observed_in_both_groups"] == 0
    assert part["overlap"]["coverage_mask_total_variation"] == 1
    assert part["informative_sites"] == 0
    assert any("group_dependent_missingness" in w for w in result["warnings"])
    assert any("no_between_group_site_overlap" in w for w in result["warnings"])
    assert "exchangeable" in report["null_caveat"] and "does not" in report["null_caveat"]
    assert "No automatic recommendation" in report["recommendation"]
    assert report["exploratory"] is result["exploratory"] is True


def test_equal_marginal_coverage_but_different_masks_warns():
    # Every site has exactly 50% coverage in each group, but joint masks differ.
    matrix = [[1, 1, None, None]] * 3 + [[None, None, 1, 1]] * 3
    matrix += [[1, None, 1, None]] * 3 + [[None, 1, None, 1]] * 3
    part = pc.prepare_orientation(groups_for(matrix, 6), "L", "+")
    assert part["overlap"]["sites_with_unequal_group_coverage"] == 0
    assert part["overlap"]["coverage_mask_total_variation"] == 1
    assert any("group_dependent_missingness" in w for w in part["warnings"])


def test_state_dependent_missingness_can_create_an_effect_under_same_complete_profiles(tmp_path):
    # Both groups have exactly three (1,0) and three (0,1) complete profiles.
    # Selectively missing low calls in L creates an apparent P-L difference.
    # This is a counterexample to missingness correction, not a calibration test.
    full_group = [[1, 0]] * 3 + [[0, 1]] * 3
    complete = pc.prepare_orientation(groups_for(full_group * 2, 6), "L", "+")
    assert pc.exact_residual_test([complete])["delta_pp_exact"] == 0
    assert pc.exact_residual_test([complete])["p_exact"] == 1
    observed = [[1, None]] * 3 + [[None, 1]] * 3 + full_group
    report = pc.analyze_site_calls(save_rows(tmp_path, rows_for(observed, n_baseline=6)))
    result = report["comparisons"][0]
    part = result["strata"][0]
    assert part["n_baseline"] == part["n_state"] == 6
    assert part["informative_sites_with_group_overlap"] == 2
    assert [s["pooled_reads"] for s in part["sites"]] == [9, 9]
    assert [s["pooled_mean"] for s in part["sites"]] == [F(2, 3), F(2, 3)]
    assert part["baseline_residuals"] == [F(1, 3)] * 6
    assert part["state_residuals"] == [F(-1, 6)] * 6
    assert result["delta_pp_exact"] == -50
    assert result["p_exact"] == F(1, 462)
    assert part["overlap"]["coverage_mask_total_variation"] == 1
    assert any("group_dependent_missingness" in w for w in result["warnings"])
    assert "do not correct group-dependent missingness" in report["null_caveat"]


def test_original_fourteen_column_site_calls_shape(tmp_path):
    rows = rows_for([[1, 0], [1, 0], [0, 1], [0, None], [1, None], [0, 1]])
    for row in rows:
        row["query_position"] = 100 + row["position"]
        score = "" if row["call_status"] == "unknown" else 255 * row["high_call"]
        row["stored_ml"] = row["score_used"] = score
    rows[2].update(call_status="implicit_low", high_call=0, stored_ml="", score_used=0)
    fields = ["sample", "locus", "read_id", "junction", "role", "modification", "strand",
              "contig", "position", "query_position", "stored_ml", "score_used", "call_status", "high_call"]
    source = tmp_path / "real_shape.tsv.gz"
    with gzip.open(source, "wt", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows + [rows[0]])
    report = pc.analyze_site_calls(source)
    assert report["status"] == "complete"
    assert report["input_counts"]["duplicate_rows"] == 1
    assert report["input_counts"]["unknown_observations"] == 2
    assert report["comparisons"][0]["strata"][0]["positions"] == [10]


def test_duplicates_collapsed_before_threshold_and_conflicts_rejected(tmp_path):
    rows = rows_for([[1]] * 3 + [[0]] * 3)
    expected = pc.analyze_site_calls(save_rows(tmp_path, rows))
    actual = pc.analyze_site_calls(save_rows(tmp_path, rows + rows))
    assert actual["comparisons"] == expected["comparisons"]
    assert actual["input_counts"]["duplicate_rows"] == 6
    with pytest.raises(ValueError, match="Conflicting duplicate"):
        pc.analyze_site_calls(save_rows(tmp_path, rows + [dict(rows[0], high_call=0)]))
    # Repeating five unique observations must not manufacture pooled coverage 6.
    report = pc.analyze_site_calls(save_rows(tmp_path, rows[:5] * 3))
    assert report["comparisons"][0]["strata"][0]["selected_sites"] == 0


@pytest.mark.parametrize("change", [dict(junction="P"), dict(strand="-"), dict(contig="other")])
def test_conflicting_molecule_assignment_rejected(tmp_path, change):
    rows = rows_for([[1]] * 3 + [[0]] * 3)
    with pytest.raises(ValueError, match="Conflicting read group"):
        pc.analyze_site_calls(save_rows(tmp_path, rows + [dict(rows[0], **change)]))


def test_channels_roles_samples_and_loci_are_separate(tmp_path):
    rows = rows_for([[1]] * 3 + [[0]] * 3)
    rows += rows_for([[0]] * 3 + [[1]] * 3, modification="5hmC")
    rows += rows_for([[0]] * 3 + [[1]] * 3, baseline="R", modification="6mA")
    rows += [dict(r, sample="second", locus="J2") for r in rows[:6]]
    rows += [dict(rows[0], role="host_left", junction="B")]
    path = save_rows(tmp_path, rows)
    report = pc.analyze_site_calls(path)
    assert len(report["comparisons"]) == 4
    results = {(r["sample"], r["contrast"], r["modification"]): r["delta_pp_exact"] for r in report["comparisons"]}
    assert results["fixture", "P-L", "5mC"] == -100
    assert results["fixture", "P-L", "5hmC"] == 100
    assert results["fixture", "P-R", "6mA"] == 100
    selected = pc.analyze_site_calls(path, sample="fixture", locus="J1", modification="6mA", contrast="P-R")
    assert len(selected["comparisons"]) == 1
    assert selected["comparisons"][0]["role"] == "tail"


def test_sparse_channel_does_not_borrow_reads_from_another(tmp_path):
    rows = rows_for([[1]] * 3 + [[0]] * 3)
    rows += [dict(rows[i], modification="5hmC") for i in (0, 1, 3, 4)]
    results = pc.analyze_site_calls(save_rows(tmp_path, rows))["comparisons"]
    assert {r["modification"]: r["test_status"] for r in results} == {"5mC": "tested", "5hmC": "untestable"}


def test_repeated_profiles_keep_combinatorial_weights_without_enumeration():
    result = pc.exact_residual_test([dict(baseline_residuals=[F(0)] * 20,
                                         state_residuals=[F(1, 7)] * 20)], state_limit=21)
    assert result["p_numerator"] == 2
    assert result["p_denominator"] == math.comb(40, 20)
    assert result["p_exact"] == F(2, math.comb(40, 20))


def test_state_limit_is_explicit_and_never_p_one(tmp_path):
    with pytest.raises(ExactLimit, match="no p-value was substituted"):
        pc.exact_residual_test([dict(baseline_residuals=[F(0)] * 3,
                                    state_residuals=[F(1, 7)] * 3)], state_limit=1)
    source = save_rows(tmp_path, rows_for([[0]] * 3 + [[1]] * 3))
    output = tmp_path / "incomplete.json"
    summary = pc.run_partial_coverage(source, output, state_limit=1)
    report = json.loads(output.read_text())
    result = report["comparisons"][0]
    assert summary["status"] == "incomplete"
    assert result["test_status"] == "exact_state_limit"
    assert result["p_value"] is result["p_numerator"] is result["p_denominator"] is None
    assert result["delta_pp_exact"] == "100/1"
    assert ">1 states" in result["error"]


def test_true_zero_effect_is_exact_without_dp():
    result = pc.exact_residual_test([dict(baseline_residuals=[F(1, 7)] * 50,
                                         state_residuals=[F(1, 7)] * 50)], state_limit=1)
    assert result["p_numerator"] == result["p_denominator"] == math.comb(100, 50)


def test_tiny_nonzero_rational_effect_is_not_rounded_to_zero():
    result = pc.exact_residual_test([dict(baseline_residuals=[F(0)] * 3,
                                         state_residuals=[F(1, 10**30)] * 3)])
    assert result["p_exact"] == F(1, 10)
    assert result["delta_pp_exact"] == F(1, 10**28)


@pytest.mark.parametrize("value", [0.5, float("nan"), "1/2", None])
def test_exact_test_refuses_inexact_input(value):
    with pytest.raises(ValueError, match="Fraction"):
        pc.exact_residual_test([dict(baseline_residuals=[value], state_residuals=[F(0)])])


@pytest.mark.parametrize("changes", [dict(high_call=""), dict(high_call="2"),
                                      dict(call_status="unknown", high_call="0"),
                                      dict(call_status="implicit_low", high_call="1"),
                                      dict(modification="5modC"), dict(position=-1), dict(strand="?")])
def test_invalid_site_calls_fail_clearly(tmp_path, changes):
    rows = rows_for([[1]] * 3 + [[0]] * 3)
    rows[0].update(changes)
    with pytest.raises(ValueError):
        pc.analyze_site_calls(save_rows(tmp_path, rows))


def test_missing_state_no_sites_empty_filter_and_invalid_bounds(tmp_path):
    rows = rows_for([[1]] * 6, n_baseline=6)
    path = save_rows(tmp_path, rows)
    result = pc.analyze_site_calls(path)["comparisons"][0]
    assert result["strata"][0]["test_status"] == "missing_state"
    assert result["p_value"] is None
    for options in (dict(min_site_reads=5), dict(state_limit=0), dict(state_limit=1.5), dict(sample="absent")):
        with pytest.raises(ValueError):
            pc.analyze_site_calls(path, **options)


def test_gzip_plain_order_invariance_and_exclusive_output(tmp_path):
    rows = rows_for([[1, 0], [1, 1], [0, 1], [0, 0], [1, 1], [0, 1]])
    compressed = save_rows(tmp_path, rows)
    plain = save_rows(tmp_path, list(reversed(rows)), "site_calls.tsv")
    assert pc.analyze_site_calls(compressed)["comparisons"] == pc.analyze_site_calls(plain)["comparisons"]
    output = tmp_path / "report.json"
    summary = pc.run_partial_coverage(compressed, output)
    assert summary["status"] == "complete"
    before = output.read_bytes()
    report = json.loads(before)
    result = report["comparisons"][0]
    assert F(result["p_exact"]) == F(result["p_numerator"], result["p_denominator"])
    assert all(isinstance(r["residual"], str) for r in result["strata"][0]["reads"])
    with pytest.raises(ValueError, match="already exists"):
        pc.run_partial_coverage(compressed, output)
    assert output.read_bytes() == before


def test_module_is_usable_without_scientific_runtime(tmp_path):
    source = save_rows(tmp_path, rows_for([[0]] * 3 + [[1]] * 3))
    # -S excludes site-packages, proving no numpy/pysam/scipy/rich/GPU dependency.
    result = subprocess.run([sys.executable, "-S", "-m", "metaphaser.partial_coverage",
                             "--site-calls", str(source), "--output", str(tmp_path / "standalone.json")],
                            text=True, capture_output=True,
                            env=dict(os.environ, PYTHONPATH=str(Path(pc.__file__).resolve().parents[1])))
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["tested"] == 1


@pytest.fixture
def integrated_cli():
    from metaphaser import cli
    return cli


def test_integrated_cli_parsing_and_dispatch(integrated_cli, tmp_path, capsys, monkeypatch):
    source = save_rows(tmp_path, rows_for([[0]] * 3 + [[1]] * 3))
    output = tmp_path / "cli.json"
    argv = ["partial-coverage", "--site-calls", str(source), "--output", str(output)]
    args = integrated_cli.parser().parse_args(argv)
    assert args.min_site_reads == 6 and args.exact_state_limit == pc.DEFAULT_STATE_LIMIT
    assert args.command == "partial-coverage"
    monkeypatch.setitem(sys.modules, "metaphaser.pipeline", None)
    assert integrated_cli.main(argv) == 0
    assert json.loads(capsys.readouterr().out)["tested"] == 1
    assert output.is_file()
    assert integrated_cli.main(argv[:-1] + [str(tmp_path / "limited.json"), "--exact-state-limit", "1"]) == 2
    assert json.loads(capsys.readouterr().out)["errors"]
    # Existing primary run parsing retains its original defaults.
    primary = integrated_cli.parser().parse_args(["run", "--assembly", "a", "--modbam", "b", "--outdir", "o"])
    assert primary.exact_state_limit == 2_000_000
    assert primary.bootstrap_replicates == 0


@pytest.mark.parametrize("extra", [["--min-site-reads", "5"], ["--exact-state-limit", "0"],
                                    ["--modification", "5modC"], ["--contrast", "B-L"],
                                    ["--modbam", "forbidden.bam"]])
def test_integrated_cli_rejects_invalid_options(integrated_cli, extra):
    with pytest.raises(SystemExit) as error:
        integrated_cli.parser().parse_args(["partial-coverage", "--site-calls", "saved.tsv.gz",
                                        "--output", "report.json"] + extra)
    assert error.value.code == 2
