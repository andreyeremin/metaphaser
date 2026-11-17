"""Modification profiles from original MM/ML tags and accepted junctions."""
from __future__ import annotations

import collections
import csv
import gzip
import hashlib
import json
import re

import pysam

from .io import write_tsv
from .stratified_exact import stratified, ExactLimit

CHANNELS = {"6mA":("A","a"), "5mC":("C","m"), "5hmC":("C","h")}
KINDS = ["6mA","5mC","5hmC"]
CONTRASTS = [("L","P","head"),("R","P","tail"),("L","B","host_left"),("R","B","host_right")]


def parse_modifications(a):
    seq = a.get_forward_sequence().upper()
    modes,explicit = {},collections.defaultdict(dict)
    if not a.has_tag("MM") or not a.has_tag("ML"):
        return seq,modes,explicit
    if a.has_tag("MN") and a.get_tag("MN") != a.query_length:
        raise ValueError("MN length differs from the stored sequence length")
    for group in a.get_tag("MM").split(";"):
        if not group:
            continue
        m = re.fullmatch(r"([ACGTUN])([+-])([a-z]+|\d+)([.?]?)",group.split(",")[0])
        if not m:
            raise ValueError(f"Malformed MM group: {group.split(',')[0]}")
        base,strand,codes,mode=m.groups()
        if strand == "+":
            for code in ([codes] if codes.isdigit() else codes):
                modes[base,code]=mode or "."
    decoded=a.modified_bases_forward
    if decoded is None and any("," in group and group.split(",",1)[1] for group in a.get_tag("MM").split(";") if group):
        raise ValueError("MM/ML tags could not be decoded")
    for (base,strand,code),values in (decoded or {}).items():
        if strand != 0:
            continue
        for pos,score in values:
            if score >= 0:
                explicit[base,str(code)][pos]=score
    return seq,modes,explicit


def call_at(base,code,pos,modes,explicit):
    value=explicit.get((base,code),{}).get(pos)
    if value is not None:
        return value,"explicit"
    if modes.get((base,code)) == ".":
        return 0,"implicit_low"
    return None,"unknown"


def extract(sample,calls,seqs,folder,p,logger,context="cpg",selection="auto"):
    selected=collections.defaultdict(list)
    for c in calls:
        if c["simple"]:
            selected[c["read_id"]].append(c)
    observations=collections.defaultdict(dict)
    seen=set()
    measured=set()
    available=set()
    diagnostics=collections.Counter()
    wanted=set(CHANNELS) if selection=="auto" else set(selection.split(","))
    output=folder/"site_calls.tsv.gz"
    fields=["sample","locus","read_id","junction","role","modification","strand",
            "contig","position","query_position","stored_ml","score_used","call_status","high_call"]
    with gzip.open(output.with_suffix(".gz.tmp"),"wt",newline="") as target:
        writer=csv.DictWriter(target,fieldnames=fields,delimiter="\t",lineterminator="\n")
        writer.writeheader()
        if selected:
            with pysam.AlignmentFile(sample["bam"],"rb",check_sq=False) as bam:
                for number,a in enumerate(bam.fetch(until_eof=True),1):
                    if number%100_000==0:
                        logger.event(f"Modification scan: {number:,} records, {len(seen):,}/{len(selected):,} selected reads")
                    rid=a.query_name
                    if a.is_secondary or a.is_supplementary or a.is_qcfail or a.is_duplicate or rid not in selected or rid in seen:
                        continue
                    if not a.query_sequence or any(op==5 for op,n in a.cigartuples or []):
                        continue
                    seen.add(rid)
                    try:
                        query,modes,explicit=parse_modifications(a)
                    except (ValueError,TypeError) as exc:
                        diagnostics["malformed_modification_reads"]+=1
                        raise ValueError(f"Cannot decode MM/ML tags for read {rid}: {exc}") from exc
                    available.update(k for k,(base,code) in CHANNELS.items() if (base,code) in modes)
                    dedup=set()
                    for c in selected[rid]:
                        strand="-" if c["reverse"] else "+"
                        sequence=seqs[c["contig"]]
                        for q,r,role in c["pairs"]:
                            if not 0<=q<len(query) or not 0<=r<len(sequence):
                                diagnostics["out_of_range_pairs"]+=1
                                continue
                            base=query[q]
                            ref=sequence[r] if strand=="+" else sequence[r].translate(str.maketrans("ACGT","TGCA"))
                            if base!=ref:
                                diagnostics["base_mismatches"]+=1
                                continue
                            values=[]
                            if base=="A" and "6mA" in wanted and ("A","a") in modes:
                                v,status=call_at("A","a",q,modes,explicit)
                                values.append(("6mA",v,status))
                            is_cpg=(sequence[r:r+2]=="CG") if strand=="+" else (r>0 and sequence[r-1:r+1]=="CG")
                            if base=="C" and (context=="all" or is_cpg):
                                m,ms=call_at("C","m",q,modes,explicit)
                                h,hs=call_at("C","h",q,modes,explicit)
                                if "5mC" in wanted and ("C","m") in modes:values.append(("5mC",m,ms))
                                if "5hmC" in wanted and ("C","h") in modes:values.append(("5hmC",h,hs))
                            for kind,value,status in values:
                                key=(c["locus"],c["kind"],role,kind,strand,rid)
                                item=(key,r)
                                if item in dedup:
                                    continue
                                dedup.add(item)
                                high=None if value is None else int(value>=p["ml_cutoff"])
                                if high is not None:
                                    observations[key][r]=high
                                    if kind in CHANNELS:measured.add(kind)
                                diagnostics[status]+=1
                                # Omitted bytes are not printed as measured probability zero.
                                writer.writerow(dict(sample=sample["sample"],locus=c["locus"],read_id=rid,
                                    junction=c["kind"],role=role,modification=kind,strand=strand,contig=c["contig"],
                                    position=r,query_position=q,stored_ml=value if status=="explicit" else "",score_used="" if value is None else value,
                                    call_status=status,high_call="" if high is None else high))
                    if len(seen)==len(selected):
                        break
    output.with_suffix(".gz.tmp").replace(output)
    return observations,{"selected_reads":len(selected),"found_reads":len(seen),
                         "available_channels":sorted(available),"selected_channels":sorted(available&wanted),
                         "measured_channels":sorted(measured),"diagnostics":dict(diagnostics)}


def compare(observations,loci,folder,options,logger):
    """Separate channels, fixed positions per orientation, exact stratified test."""
    index=collections.defaultdict(list)
    for key,values in observations.items():
        if key[3] not in CHANNELS:
            raise ValueError('Only separate modification channels are permitted')
        if values:index[key[:5]].append((key[5],values))
    rows,profiles,panels,haplotypes=[],[],[],[]
    for number,locus in enumerate(loci,1):
        lid=locus['locus']
        for baseline,state,role in CONTRASTS:
            for kind in KINDS:
                strata=[];parts=[];observed=False
                common=dict(sample=locus['sample'],locus=lid,region=locus['region'],
                    structural_class=locus['structural_class'],contrast=state+'-'+baseline,
                    role=role,modification=kind)
                for strand in ['+','-']:
                    groups=[index.get((lid,c,role,kind,strand),[]) for c in [baseline,state]]
                    if not any(groups):continue
                    observed=True
                    ids=[{rid for rid,_ in g} for g in groups]
                    if any(len(g)!=len(ii) for g,ii in zip(groups,ids)) or ids[0]&ids[1]:
                        raise ValueError('Duplicate or overlapping molecules within a comparison')
                    shared=set.intersection(*(set(v) for g in groups for _,v in g)) if all(groups) else set()
                    panel=sorted(shared)
                    tag=hashlib.sha256(json.dumps([lid,baseline,state,role,kind,strand,panel]).encode()).hexdigest()[:16]
                    counts=[[sum(v[pos] for pos in panel) for _,v in group] for group in groups]
                    means=[100*sum(c)/(len(c)*len(panel)) if c and panel else None for c in counts]
                    eligible=all(groups) and bool(panel)
                    status='descriptive' if eligible else 'missing_state' if not all(groups) else 'no_shared_positions'
                    row=dict(common,strand=strand,panel_id=tag,n_baseline=len(groups[0]),n_state=len(groups[1]),
                        shared_positions=len(panel),baseline_percent=means[0],state_percent=means[1],
                        delta_pp=means[1]-means[0] if eligible else None,test_status=status,p_value=None,
                        test='orientation_summary',orientations=strand,p_numerator=None,p_denominator=None)
                    rows.append(row)
                    if eligible:
                        strata.append(dict(baseline_counts=counts[0],state_counts=counts[1],shared_positions=len(panel)))
                        parts.append(row)
                    panels.append(dict(panel_id=tag,locus=lid,contrast=row['contrast'],role=role,modification=kind,strand=strand,positions=panel))
                    for config,group,highs in zip([baseline,state],groups,counts):
                        patterns=collections.defaultdict(list)
                        for (rid,values),high in zip(group,highs):
                            pattern=''.join(str(values[pos]) for pos in panel)
                            profiles.append(dict(sample=locus['sample'],locus=lid,panel_id=tag,contrast=row['contrast'],
                                modification=kind,strand=strand,state=config,read_id=rid,shared_positions=len(panel),
                                high_calls=high if panel else None,high_call_fraction=high/len(panel) if panel else None,
                                pattern=pattern,test_status=status))
                            if panel:patterns[pattern].append(rid)
                        for pattern,rids in sorted(patterns.items()):
                            haplotypes.append(dict(sample=locus['sample'],locus=lid,panel_id=tag,contrast=row['contrast'],
                                modification=kind,strand=strand,state=config,pattern=pattern,n_reads=len(rids),read_ids=sorted(rids)))
                if observed:
                    baseline_mean=sum(r['baseline_percent'] for r in parts)/len(parts) if parts else None
                    state_mean=sum(r['state_percent'] for r in parts)/len(parts) if parts else None
                    row=dict(common,strand='combined',panel_id='',n_baseline=sum(r['n_baseline'] for r in parts),
                        n_state=sum(r['n_state'] for r in parts),shared_positions=sum(r['shared_positions'] for r in parts),
                        baseline_percent=baseline_mean,state_percent=state_mean,
                        delta_pp=state_mean-baseline_mean if parts else None,test_status='no_eligible_orientation',
                        p_value=None,test='exact_within_orientation_equal_weights',
                        orientations=','.join(r['strand'] for r in parts),p_numerator=None,p_denominator=None)
                    if strata:
                        try:
                            row.update(stratified(strata,limit=options.exact_state_limit));row['test_status']='tested'
                        except ExactLimit as exc:
                            row['test_status']='exact_state_limit';row['test_note']=str(exc)
                    rows.append(row)
        if number%25==0:logger.event(f'Modification comparisons: {number:,}/{len(loci):,} loci')
    write_tsv(folder/'comparisons.tsv',rows,['sample','locus','region','structural_class','contrast','role','modification','strand',
        'panel_id','n_baseline','n_state','shared_positions','baseline_percent','state_percent','delta_pp',
        'test_status','p_value','test','orientations','p_numerator','p_denominator','test_note'])
    write_tsv(folder/'read_profiles.tsv.gz',profiles,['sample','locus','panel_id','contrast','modification','strand','state',
        'read_id','shared_positions','high_calls','high_call_fraction','pattern','test_status'])
    write_tsv(folder/'panels.tsv.gz',panels,['panel_id','locus','contrast','role','modification','strand','positions'])
    write_tsv(folder/'epihaplotypes.tsv.gz',haplotypes,['sample','locus','panel_id','contrast','modification','strand','state','pattern','n_reads','read_ids'])
    return rows
