"""Defaults and identities for the tested short-junction protocol."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import shutil
import subprocess
from pathlib import Path

DEFAULTS = {
    "min_length": 200, "max_length": 1_000_000, "min_anchor": 75,
    "discovery_identity": .85, "max_query_gap": 250, "min_support": 3,
    "short_max": 2000, "short_tolerance": 25, "long_tolerance": 200,
    "breakpoint_hypotheses": 5, "external_flank": 2000, "base_window": 250,
    "mod_window": 1000, "paired_fraction": .90, "accept_identity": .95,
    "max_core_indel": 10, "min_edit_margin": 3, "edit_margin_fraction": .02,
    "ml_cutoff": 205, "min_discovery_reads": 2, "refined_breakpoint_hypotheses": 5,
}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for part in iter(lambda: stream.read(2**20), b""):
            h.update(part)
    return h.hexdigest()


def code_identity():
    root = Path(__file__).parent
    h = hashlib.sha256()
    for file in sorted(root.rglob("*.py")):
        h.update(str(file.relative_to(root)).encode())
        h.update(file.read_bytes())
    return h.hexdigest()


def file_identity(path, full_hash=False):
    path = Path(path).resolve(strict=True)
    st = path.stat()
    value = {"path": str(path), "bytes": st.st_size, "mtime_ns": st.st_mtime_ns}
    if full_hash:
        value["sha256"] = sha256(path)
    return value


def find_executable(value, name):
    candidate = value or shutil.which(name)
    if candidate is None:
        raise ValueError(f"{name} was not found. Put it on PATH or use --{name} /path/to/{name}.")
    path = Path(candidate).expanduser().resolve()
    if not path.is_file() or not os.access(path, os.X_OK):
        raise ValueError(f"Executable is not available: {path}")
    proc = subprocess.run([str(path), "--version"], capture_output=True, text=True, timeout=20)
    if proc.returncode:
        raise ValueError(f"Cannot run {name}: {proc.stderr.strip()}")
    return str(path), (proc.stdout or proc.stderr).splitlines()[0]


def runtime_info(minimap2):
    import platform
    versions = {k: importlib.metadata.version(k) for k in ("pysam", "numpy", "scipy", "edlib", "rich")}
    return {"python": platform.python_version(), "platform": platform.platform(),
            "packages": versions, "minimap2": minimap2}


def parameters(args):
    p = dict(DEFAULTS)
    if args.parameters:
        changes = json.loads(Path(args.parameters).read_text())
        unknown = set(changes)-set(p)
        if unknown:
            raise ValueError("Unknown analysis parameters: " + ", ".join(sorted(unknown)))
        p.update(changes)
    for field in ("min_length", "max_length", "min_anchor", "max_query_gap", "min_support", "ml_cutoff", "min_discovery_reads"):
        value = getattr(args, field, None)
        if value is not None:
            p[field] = value
    integers = set(DEFAULTS)-{"discovery_identity", "paired_fraction", "accept_identity", "edit_margin_fraction"}
    for key in integers:
        if isinstance(p[key], bool) or not isinstance(p[key], int) or p[key] < (0 if key in {"max_query_gap", "ml_cutoff"} else 1):
            raise ValueError(f"{key} must be a valid integer.")
    for key in ("discovery_identity", "paired_fraction", "accept_identity", "edit_margin_fraction"):
        if not 0 < p[key] <= 1:
            raise ValueError(f"{key} must be greater than zero and at most one.")
    if not 0 <= p["ml_cutoff"] <= 255 or p["min_length"] > p["max_length"]:
        raise ValueError("Invalid ML cutoff or interval length range.")
    if p["base_window"] < p["min_anchor"] or p["external_flank"] < p["min_anchor"]:
        raise ValueError("The test window and external flanks must be at least min_anchor bases long.")
    return p
