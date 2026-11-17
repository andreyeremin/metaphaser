"""Small deterministic input files for installation and end-to-end checks."""
from __future__ import annotations

import array
import json
import random
from pathlib import Path

import pysam


def make_demo(folder,length=800,per_state=6,ordinary_only=False,include_6ma=True):
    folder=Path(folder)
    folder.mkdir(parents=True,exist_ok=True)
    if any(folder.iterdir()):raise ValueError(f"Demo input directory is not empty: {folder}")
    rng=random.Random(20261016+length)
    seq="".join(rng.choices("ACGT",k=12000+length))
    a,b=6000,6000+length
    h=min(2000,length//2)
    states={"L":seq[a-1000:a+h],"R":seq[b-h:b+1000],
            "P":seq[b-h:b]+seq[a:a+h],"B":seq[a-1000:a]+seq[b:b+1000]}
    reads=[]
    if not ordinary_only:
        for state,s in states.items():
            for orientation in ["+","-"]:
                for n in range(per_state):
                    read=s if orientation=="+" else s.translate(str.maketrans("ACGT","TGCA"))[::-1]
                    reads.append((f"demo_{state}_{'forward' if orientation=='+' else 'reverse'}_{n:02}",read,state,orientation))
    for n in range(8):
        reads.append((f"ordinary_{n:02}",seq[1000+n*40:2600+n*40],"ordinary","+"))
    assembly=folder/"assembly.fasta"
    assembly.write_text(">demo_contig\n"+seq+"\n")
    modbam=folder/"reads.mod.bam"
    truth=[]
    with pysam.AlignmentFile(modbam,"wb",header={"HD":{"VN":"1.6"},"SQ":[{"SN":"demo_contig","LN":len(seq)}]}) as bam:
        for rid,read,state,orientation in reads:
            r=pysam.AlignedSegment(bam.header)
            r.query_name=rid
            r.query_sequence=read
            r.flag=4
            r.query_qualities=pysam.qualitystring_to_array("I"*len(read))
            cpos=[i for i,v in enumerate(read) if v=="C"]
            apos=[i for i,v in enumerate(read) if v=="A"]
            headers=[]
            ml=[]
            channels=[("C","m",cpos),("C","h",cpos)]
            if include_6ma:channels.append(("A","a",apos))
            for base,code,positions in channels:
                headers.append(base+"+"+code+"?,"+",".join("0" for _ in positions)+";")
                for i,pos in enumerate(positions):
                    high=state in {"L","R","B"} and code in {"m","a"} and i%3==0
                    ml.append(230 if high else 0)
            r.set_tag("MM","".join(headers))
            r.set_tag("ML",array.array("B",ml))
            r.set_tag("MN",len(read))
            r.set_tag("qs",40.0)
            bam.write(r)
            truth.append(dict(read_id=rid,state=state,orientation=orientation))
    (folder/"truth.json").write_text(json.dumps(dict(interval=[a,b],length=length,reads=truth),indent=2))
    return assembly,modbam
