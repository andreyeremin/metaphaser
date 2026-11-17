import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from metaphaser.demo import make_demo


def command(*args):
    return subprocess.run([sys.executable,"-m","metaphaser",*map(str,args)],capture_output=True,text=True)


@pytest.mark.integration
@pytest.mark.parametrize("length",[200,430,800,1500,3000])
def test_cli_classifies_all_states_and_phases_modifications(tmp_path,length):
    mm=shutil.which("minimap2")
    if not mm:pytest.skip("minimap2 is not installed")
    assembly,bam=make_demo(tmp_path/"inputs",length=length)
    out=tmp_path/"result"
    args=["run","--assembly",assembly,"--modbam",bam,"--sample","test","--outdir",out,"--threads","1","--minimap2",mm,"--no-color"]
    result=command(*args)
    assert result.returncode==0,result.stderr
    summary=json.loads((out/"summary.json").read_text())
    assert summary["host_linked_P_groups"]>=1
    assert summary["all_four_groups"]>=1
    with gzip.open(out/"read_classifications.tsv.gz","rt") as f:rows=list(csv.DictReader(f,delimiter="\t"))
    assert len(rows)==56 and len({r["read_id"] for r in rows})==56
    byid={r["read_id"]:r for r in rows}
    assert byid["demo_P_forward_00"]["classification"]=="P"
    assert byid["demo_B_forward_00"]["classification"]=="B"
    assert byid["demo_L_forward_00"]["classification"]=="integrated"
    assert byid["ordinary_00"]["classification"]=="unassigned"
    with (out/"modifications/comparisons.tsv").open() as f:comparisons=list(csv.DictReader(f,delimiter="\t"))
    hits=[r for r in comparisons if r["contrast"]=="P-L" and r["modification"]=="5mC" and r["strand"]=="combined"]
    assert hits and any(r["test_status"]=="tested" and float(r["delta_pp"])<0 for r in hits)
    assert (out/"report.html").is_file()
    resumed=command(*args,"--resume")
    assert resumed.returncode==0,resumed.stderr
    assert "Reused:" in resumed.stderr
    assert "Running minimap2" not in resumed.stderr
    # A missing mapping must invalidate that stage and downstream validation.
    (out/"work/synthetic.bam").unlink()
    repaired=command(*args,"--resume")
    assert repaired.returncode==0,repaired.stderr
    assert "Running minimap2" in repaired.stderr
    assert json.loads((out/"summary.json").read_text())["all_four_groups"]==summary["all_four_groups"]
    rejected=command(*args,"--resume","--ml-cutoff","200")
    assert rejected.returncode==2 and "Resume identity changed" in rejected.stderr
    # An output collision must not overwrite the completed run.
    rejected=command(*args)
    assert rejected.returncode==2 and "not empty" in rejected.stderr


@pytest.mark.integration
def test_no_junction_input_finishes_with_empty_well_formed_tables(tmp_path):
    mm=shutil.which("minimap2")
    if not mm:pytest.skip("minimap2 is not installed")
    fa,bam=make_demo(tmp_path/"inputs",ordinary_only=True)
    out=tmp_path/"result"
    r=command("run","--assembly",fa,"--modbam",bam,"--outdir",out,"--minimap2",mm,"--threads","1","--no-color")
    assert r.returncode==0,r.stderr
    assert json.loads((out/"summary.json").read_text())["host_linked_P_groups"]==0
    assert (out/"elements.tsv").read_text().startswith("sample\tlocus")
    assert "<html" in (out/"report.html").read_text()


@pytest.mark.integration
def test_missing_6ma_channel_does_not_become_a_zero_measurement(tmp_path):
    mm=shutil.which("minimap2")
    if not mm:pytest.skip("minimap2 is not installed")
    fa,bam=make_demo(tmp_path/"inputs",include_6ma=False)
    out=tmp_path/"result"
    r=command("run","-a",fa,"-b",bam,"-o",out,"--minimap2",mm,"-t","1","--no-color")
    assert r.returncode==0,r.stderr
    channels=json.loads((out/"modifications/channels.json").read_text())
    assert "6mA" not in channels["measured_channels"]
    assert "6mA" not in (out/"modifications/comparisons.tsv").read_text()
    assert json.loads((out/"summary.json").read_text())["all_four_groups"]==1


def test_cli_invalid_options_have_nonzero_exit_status():
    r=command("run","--assembly","missing.fa","--modbam","missing.bam","--outdir","unused","--bootstrap-replicates","2")
    assert r.returncode==2
    assert "must be 0" in r.stderr
