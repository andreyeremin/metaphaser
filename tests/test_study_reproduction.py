"""Lightweight checks of the frozen publication inputs and overwrite guard."""
import importlib.util
from pathlib import Path
import subprocess
import sys

import pytest


STUDY = Path(__file__).resolve().parents[1] / "study" / "383"
spec = importlib.util.spec_from_file_location("study_reproduce", STUDY / "reproduce.py")
study = importlib.util.module_from_spec(spec)
spec.loader.exec_module(study)


def test_frozen_inputs_keep_the_complete_primary_family():
    study.verify_inputs(STUDY / "counts.tsv", STUDY / "primary_comparisons.tsv")
    assert len(study.read_table(STUDY / "counts.tsv")) == 9793
    assert len(study.read_table(STUDY / "primary_comparisons.tsv")) == 383


def test_truncated_family_is_rejected_before_calculation(tmp_path):
    truncated = tmp_path / "primary.tsv"
    lines = (STUDY / "primary_comparisons.tsv").read_text().splitlines(True)
    truncated.write_text("".join(lines[:-1]))
    with pytest.raises(ValueError, match="frozen 383-test"):
        study.verify_inputs(STUDY / "counts.tsv", truncated)


def test_existing_output_is_preserved_before_inputs_are_read(tmp_path):
    sentinel = tmp_path / "existing.txt"
    sentinel.write_text("preserve this output")
    result = subprocess.run(
        [sys.executable, str(STUDY / "reproduce.py"), "--counts", "absent.tsv",
         "--reference", "absent.tsv", "--output", str(tmp_path)],
        capture_output=True, text=True, timeout=15,
    )
    assert result.returncode != 0
    assert "Output already exists" in result.stderr
    assert sentinel.read_text() == "preserve this output"
