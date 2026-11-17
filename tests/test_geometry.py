import random

import pytest

from metaphaser._vendor import junctions_v17 as core
from metaphaser.config import DEFAULTS
from metaphaser.io import read_fasta
from metaphaser.reporting import classify_loci


@pytest.mark.parametrize("reverse",[False,True])
def test_200_base_backward_interval_is_discovered(reverse):
    if reverse:
        coords=[(1000,1100),(1100,1200)]
    else:
        coords=[(1100,1200),(1000,1100)]
    pieces=[dict(contig="ctg",r0=a,r1=b,q0=i*100,q1=(i+1)*100,reverse=reverse,
                 identity=1.,paired=100,mapq=60,route="test") for i,(a,b) in enumerate(coords)]
    rows=core.split_events("read",pieces,DEFAULTS,"test")
    assert len(rows)==1 and rows[0]["length"]==200


def test_short_reference_arms_do_not_repeat_the_interval():
    rng=random.Random(5)
    seq="".join(rng.choices("ACGT",k=4000))
    refs=core.build_references([dict(id="J1",contig="ctg",a=1000,b=1200,hypotheses=[(1000,1200)])],{"ctg":seq},DEFAULTS)
    r=next(r for r in refs if any(a["kind"]=="P" and a["model"]=="no_repeat" for a in r["aliases"]))
    assert len(r["seq"])==200
    assert r["seq"]==seq[1100:1200]+seq[1000:1100]


def test_grouping_cannot_chain_distant_boundaries():
    events=[]
    for i,a in enumerate([1000,1020,1040,1060]):
        for j in range(3):events.append(dict(read_id=f"r{i}_{j}",contig="ctg",a=a,b=a+800,length=800,gap=0,route="test"))
    groups=core.group_events(events,DEFAULTS)
    assert len(groups)==2
    for g in groups:assert g["amax"]-g["amin"]<=50


def test_identical_sequence_keeps_source_aliases():
    rng=random.Random(4)
    seq="".join(rng.choices("ACGT",k=4000))
    gs=[dict(id=k,contig=ctg,a=1000,b=1800,hypotheses=[(1000,1800)]) for k,ctg in [("J1","c1"),("J2","c2")]]
    refs=core.build_references(gs,{"c1":seq,"c2":seq},DEFAULTS)
    assert any({a["locus"] for a in r["aliases"]}=={"J1","J2"} for r in refs if r["aliases"][0]["kind"]=="P")


def test_b_only_and_whole_contig_join_are_not_host_linked_elements():
    def locus(name,a,b,L=0,R=0,P=0,B=0):
        return dict(sample="x",locus=name,contig="ctg",a=a,b=b,L_simple=L,R_simple=R,P_simple=P,B_simple=B)
    rows=[locus("P",100,900,L=3,P=3),locus("B",1000,1800,L=10,R=10,B=8),locus("whole",0,3000,P=5)]
    rows,regions=classify_loci(rows,{"ctg":"A"*3000},3)
    assert [r["structural_class"] for r in rows]==["host_linked_P_one_boundary","B_without_P","contig_end_join"]
    assert len(regions)==1


def test_fasta_duplicate_names_are_rejected(tmp_path):
    p=tmp_path/"bad.fa"
    p.write_text(">a\nACGT\n>a\nACGT\n")
    with pytest.raises(ValueError,match="Duplicate"):
        read_fasta(p)
