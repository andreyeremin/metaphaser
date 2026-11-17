import math
import random
import shutil

import pytest

from metaphaser._vendor import junctions_v17 as core
from metaphaser.config import DEFAULTS


@pytest.mark.integration
def test_copy_without_discovery_seeds_prevents_unique_P_assignment(tmp_path):
    mm=shutil.which("minimap2")
    if not mm:pytest.skip("minimap2 is not installed")
    rng=random.Random(158)
    duplicated="".join(rng.choices("ACGT",k=5000))
    unique="".join(rng.choices("ACGT",k=5000))
    seqs={"seeded_copy":duplicated,"unseeded_copy":duplicated,"unique":unique}
    groups=[dict(id=lid,contig=ctg,a=2000,b=2800,hypotheses=[(2000,2800)]) for lid,ctg in [("duplicate","seeded_copy"),("unique","unique")]]
    refs=core.build_references(groups,seqs,DEFAULTS)
    fa=tmp_path/"assembly.fasta"
    core.fasta_write(fa,seqs.items())
    calls=[]
    for ref in refs:
        for al in ref["aliases"]:
            if al["kind"]!="P" or al["model"]!="no_repeat":continue
            w=min(250,al["h"])
            calls.append(dict(read_id=al["locus"],ref_id=ref["id"],locus=al["locus"],contig=al["contig"],kind="P",
                model=al["model"],a=al["a"],b=al["b"],h=al["h"],window=w,query_slice=ref["seq"][al["join"]-w:al["join"]+w],
                edit_distance=0,required_margin=max(3,math.ceil(.02*2*w)),locus_ambiguous=False,simple=True))
    cfg=dict(minimap2=mm,junction_seed_length=11)
    result=core.audit_source_copies(cfg,{"fasta":str(fa)},calls,refs,seqs,tmp_path,1,DEFAULTS)
    assert result["tested_P_calls"]==2
    assert next(c for c in calls if c["locus"]=="unique")["simple"]
    assert not next(c for c in calls if c["locus"]=="duplicate")["simple"]
    assert next(c for c in calls if c["locus"]=="duplicate")["source_copy_ambiguous"]


def test_overlapping_occupied_labels_are_distinct_from_another_source_copy():
    aliases=[dict(locus="overlap",kind="L",contig="ctg",arms=[dict(start=500,end=1000),dict(start=1000,end=2000)]),
             dict(locus="other",kind="L",contig="other_ctg",arms=[dict(start=500,end=1000),dict(start=1000,end=2000)])]
    call=dict(kind="L",locus="main",contig="ctg",genomic_query_min=900,genomic_query_max=1100,
              competing_loci=["overlap"],simple=False,complex_read=False,genomic_contiguous_alternative=False)
    other=dict(call,competing_loci=["other"])
    assert core.normalize_occupied_labels([call,other],[dict(aliases=aliases)])==1
    assert call["simple"] and not other["simple"]
