import json
import random
import sys

import pysam
import pytest

from metaphaser._vendor import junctions_v17 as core
from metaphaser.config import DEFAULTS


@pytest.mark.skipif(sys.platform=="win32",reason="The validation pool uses POSIX fork")
def test_parallel_validation_preserves_calls_and_reuses_chunks(tmp_path):
    rng=random.Random(42)
    seq="".join(rng.choices("ACGT",k=4000))
    group=dict(id="J1",contig="ctg",a=1000,b=1800,hypotheses=[(1000,1800)])
    refs=core.build_references([group],{"ctg":seq},DEFAULTS)
    target=next(r for r in refs if r["aliases"][0]["kind"]=="P" and any(a["model"]=="no_repeat" for a in r["aliases"]))
    header=pysam.AlignmentHeader.from_dict({"HD":{"VN":"1.6"},"SQ":[{"SN":r["id"],"LN":len(r["seq"])} for r in refs]})
    bam=tmp_path/"validation.bam"
    with pysam.AlignmentFile(bam,"wb",header=header) as f:
        for i in range(8):
            a=pysam.AlignedSegment(header)
            a.query_name=f"P{i:04}"
            a.query_sequence=target["seq"]
            a.reference_id=header.get_tid(target["id"])
            a.reference_start=0
            a.cigartuples=[(7,len(target["seq"]))]
            a.set_tag("NM",0)
            f.write(a)
        # These complete read groups exercise actual multi-chunk partitioning.
        for i in range(1800):
            a=pysam.AlignedSegment(header)
            a.query_name=f"noise{i:05}"
            a.query_sequence="".join(rng.choices("ACGT",k=500))
            a.flag=4
            f.write(a)
    cfg=dict(indexed_competition=True,ranked_candidate_search=True,sparse_modification_pairs=True,
             validation_workers=2,validation_chunks=4,validation_checkpoint_dir=str(tmp_path/"checkpoints"))
    serial,rejected=core.validate_synthetic_ranked(cfg,refs,bam,DEFAULTS)
    assert len(serial)==8 and {c["kind"] for c in serial}=={"P"}
    parallel,parallel_rejected=core.validate_synthetic_parallel(cfg,refs,bam,DEFAULTS)
    assert serial==parallel and rejected==parallel_rejected
    info=json.loads((tmp_path/"checkpoints/summary.json").read_text())
    assert info["chunks"]>=2 and info["workers"]==2
    again,rejected_again=core.validate_synthetic_parallel(cfg,refs,bam,DEFAULTS)
    assert again==serial and rejected_again==rejected
    info=json.loads((tmp_path/"checkpoints/summary.json").read_text())
    assert info["reused_chunks"]==info["chunks"]
