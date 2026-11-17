from __future__ import annotations

import contextlib
import csv
import fcntl
import gzip
import json
import os
from pathlib import Path


def json_write(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name+".tmp")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(tmp, "wt") as f:
        json.dump(data, f, allow_nan=False, default=str, separators=(",", ":"))
        f.write("\n")
    tmp.replace(path)


def json_read(path):
    path = Path(path)
    with (gzip.open(path, "rt") if path.suffix == ".gz" else path.open()) as f:
        return json.load(f)


def write_tsv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name+".tmp")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(tmp, "wt", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k:json.dumps(v,separators=(",", ":")) if isinstance(v,(list,dict,tuple)) else
                             "" if v is None else v for k,v in row.items()})
    tmp.replace(path)


def read_fasta(path):
    result, name, parts = {}, None, []
    with (gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)) as f:
        for number, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    result[name] = "".join(parts).upper()
                words = line[1:].split()
                if not words:
                    raise ValueError(f"Empty FASTA identifier at line {number}.")
                name = words[0]
                if name in result:
                    raise ValueError(f"Duplicate FASTA identifier: {name}")
                parts = []
            elif name is None:
                raise ValueError("FASTA must start with a >identifier header.")
            else:
                if set(line.upper())-set("ACGTRYSWKMBDHVN"):
                    raise ValueError(f"Invalid DNA sequence on FASTA line {number}.")
                parts.append(line)
    if name is not None:
        result[name] = "".join(parts).upper()
    if not result or any(not seq for seq in result.values()):
        raise ValueError("The assembly contains no sequence or an empty contig.")
    return result


@contextlib.contextmanager
def run_lock(folder):
    path = Path(folder)/".metaphaser.lock"
    with path.open("a+") as f:
        try:
            fcntl.flock(f.fileno(), fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError(f"Another metaphaser process owns this output directory: {folder}") from None
        f.seek(0)
        f.truncate()
        f.write(str(os.getpid())+"\n")
        f.flush()
        try:
            yield
        finally:
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)
