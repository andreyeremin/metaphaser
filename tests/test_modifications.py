import array

import pysam
import pytest
import csv
import gzip
from types import SimpleNamespace

from metaphaser.modifications import parse_modifications,call_at,extract
from metaphaser.config import DEFAULTS


def read(mm,ml,seq="CACACG"):
    a=pysam.AlignedSegment()
    a.query_name="test"
    a.query_sequence=seq
    a.flag=4
    a.set_tag("MM",mm)
    a.set_tag("ML",array.array("B",ml))
    a.set_tag("MN",len(seq))
    return a


def test_implicit_low_and_unknown_are_distinct():
    a=read("C+m.,0;C+h?,0;",[230,20])
    _,modes,explicit=parse_modifications(a)
    assert call_at("C","m",0,modes,explicit)==(230,"explicit")
    assert call_at("C","m",2,modes,explicit)==(0,"implicit_low")
    assert call_at("C","h",2,modes,explicit)==(None,"unknown")
    assert ("A","a") not in modes


def test_reverse_stored_read_uses_forward_query_coordinates():
    a=read("C+m?,0;",[230],seq="CGTTGC")
    a.flag=16
    seq,modes,explicit=parse_modifications(a)
    assert seq=="GCAACG"
    # MM positions are specified in the original sequenced orientation.
    assert explicit["C","m"]=={1:230}


def test_mn_mismatch_is_not_silently_interpreted_as_low_calls():
    a=read("C+m.,0;",[230])
    a.set_tag("MN",999)
    with pytest.raises(ValueError,match="MN"):
        parse_modifications(a)


def test_missing_ml_has_no_usable_modification_channel():
    a=read("C+m.,0;",[230])
    a.set_tag("ML",None)
    _,modes,explicit=parse_modifications(a)
    assert not modes and not explicit


def test_only_separate_modification_channels_are_emitted(tmp_path):
    bam=tmp_path/"mods.bam"
    with pysam.AlignmentFile(bam,"wb",header={"HD":{"VN":"1.6"},"SQ":[{"SN":"c","LN":3}]}) as f:
        a=read("C+m?,0;C+h?,0;",[120,100],seq="CGA")
        f.write(a)
    calls=[dict(read_id="test",simple=True,locus="J1",contig="c",kind="L",reverse=False,pairs=[(0,0,"head")])]
    ob,info=extract({"sample":"x","bam":str(bam)},calls,{"c":"CGA"},tmp_path,DEFAULTS,SimpleNamespace(event=lambda *a,**k:None))
    assert ob["J1","L","head","5mC","+","test"][0]==0
    assert ob["J1","L","head","5hmC","+","test"][0]==0
    assert all(key[3] in {"5mC", "5hmC"} for key in ob)
    with gzip.open(tmp_path/"site_calls.tsv.gz","rt") as f:rows=list(csv.DictReader(f,delimiter="\t"))
    assert {row["modification"] for row in rows} == {"5mC", "5hmC"}
    assert "6mA" not in info["available_channels"]
