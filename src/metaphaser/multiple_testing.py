"""Study-wide BH correction of completed host-linked P comparisons."""
from __future__ import annotations

import csv
import math
import os
import tempfile
from fractions import Fraction
from pathlib import Path

from scipy.stats import false_discovery_control

HOST_CLASSES = {"host_linked_P_both_boundaries", "host_linked_P_one_boundary"}
KEY_FIELDS = ("sample", "locus", "contrast", "modification")


def adjust_primary(rows, expected_tests=None):
    """Correct a caller-selected primary family, checking exact p provenance."""
    rows = [dict(row) for row in rows]
    if not rows:
        raise ValueError("The primary family is empty.")
    keys = [tuple(row[field] for field in KEY_FIELDS) for row in rows]
    if len(set(keys)) != len(keys):
        raise ValueError("Duplicate sample/locus/contrast/modification in the primary family.")
    if expected_tests is not None and len(rows) != expected_tests:
        raise ValueError(f"Expected {expected_tests} tests, found {len(rows)}.")
    probabilities = []
    for row in rows:
        if row["contrast"] not in {"P-L", "P-R"} or row["modification"] not in {"5mC", "5hmC", "6mA"}:
            raise ValueError("Primary correction requires separate-channel P-L/P-R comparisons.")
        numerator, denominator = int(row["p_numerator"]), int(row["p_denominator"])
        if denominator <= 0 or not 0 <= numerator <= denominator:
            raise ValueError("Invalid exact probability numerator or denominator.")
        probability = float(Fraction(numerator, denominator))
        if not math.isclose(probability, float(row["p_value"]), rel_tol=1e-12, abs_tol=0):
            raise ValueError("Stored p-value disagrees with the exact probability.")
        probabilities.append(probability)
    adjusted = false_discovery_control(probabilities, method="bh")
    for row, probability in zip(rows, adjusted):
        row["q_BH"] = float(probability)
        row["BH_family_size"] = len(rows)
        row["BH_significant_005"] = bool(probability < 0.05)
    return rows


def adjust_files(paths, output, expected_tests=None):
    """Select completed primary comparisons across every supplied CLI run."""
    output = Path(output)
    if output.exists():
        raise ValueError("BH output already exists. Choose a new output filename.")
    selected = []
    for path in paths:
        with Path(path).open(newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            required = set(KEY_FIELDS) | {"strand", "structural_class", "test_status", "test",
                                          "p_value", "p_numerator", "p_denominator"}
            if not required.issubset(reader.fieldnames or []):
                raise ValueError("BH input must be a metaphaser modifications/comparisons.tsv table.")
            for row in reader:
                primary = (row["strand"] == "combined"
                           and row["structural_class"] in HOST_CLASSES
                           and row["contrast"] in {"P-L", "P-R"}
                           and row["modification"] in {"5mC", "5hmC", "6mA"})
                if not primary:
                    continue
                if row["test_status"] == "no_eligible_orientation":
                    continue
                if row["test_status"] != "tested" or row["test"] != "exact_within_orientation_equal_weights":
                    raise ValueError("A primary exact test is incomplete. Complete it before BH correction.")
                selected.append(row)
    adjusted = adjust_primary(selected, expected_tests)
    output.parent.mkdir(parents=True, exist_ok=True)
    # Publish only a complete table, with exclusive creation at the final name.
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", newline="", dir=output.parent, delete=False) as stream:
            temporary = Path(stream.name)
            fields = list(dict.fromkeys(key for row in adjusted for key in row))
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(adjusted)
        os.link(temporary, output)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return {"tests": len(adjusted), "nominal_below_005": sum(float(row["p_value"]) < .05 for row in adjusted),
            "BH_below_005": sum(row["BH_significant_005"] for row in adjusted)}
