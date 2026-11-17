"""Managed I/O workers. The parent can stop a blocked native reader safely."""
from __future__ import annotations

import argparse
from pathlib import Path

import pysam

from .io import json_read,json_write,read_fasta


def convert(args):
    out=Path(args.output)
    unsorted=out.with_name(out.name+".unsorted.bam")
    status=Path(args.status)
    count=0
    json_write(status,dict(operation="convert",records=0))
    with pysam.AlignmentFile("-","r",check_sq=False) as src,pysam.AlignmentFile(unsorted,"wb",template=src,format_options=[b'level=1']) as dst:
        for a in src:
            dst.write(a)
            count+=1
            if count%100_000==0:json_write(status,dict(operation="convert",records=count))
    json_write(status,dict(operation="sort",records=count))
    pysam.sort("-n","-@",str(args.threads),"-m","256M","-o",str(out),str(unsorted))
    unsorted.unlink()
    json_write(status,dict(operation="complete",records=count))


def audit(args):
    from ._vendor import junctions_v17 as core
    request=json_read(args.request)
    work=Path(request["work"])
    calls=json_read(work/"background_calls.json.gz")
    refs=json_read(work/"references.json.gz")
    seqs=read_fasta(request["sample"]["fasta"])
    result=core.audit_source_copies(request["config"],request["sample"],calls,refs,seqs,work,
                                   request["threads"],request["parameters"])
    json_write(work/"copy_audit_summary.json",result)
    json_write(work/"final_calls.json.gz",calls)


def main():
    parser=argparse.ArgumentParser()
    sub=parser.add_subparsers(dest="mode",required=True)
    c=sub.add_parser("convert")
    c.add_argument("--output",required=True)
    c.add_argument("--status",required=True)
    c.add_argument("--threads",type=int,default=1)
    a=sub.add_parser("audit")
    a.add_argument("--request",required=True)
    args=parser.parse_args()
    convert(args) if args.mode=="convert" else audit(args)


if __name__=="__main__":main()
