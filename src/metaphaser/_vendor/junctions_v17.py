#!/usr/bin/env python3
"""Adaptive junction analysis of existing ONT modBAMs; see METHODS.md."""
from __future__ import annotations

import argparse
import bisect
import collections
import gzip
import hashlib
import importlib.metadata
import itertools
import json
import math
import multiprocessing
import os
import random
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import edlib
import numpy as np
import pysam

VERSION = '1.7.0'
RC_TABLE = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
CONSUME_REF = {0, 2, 3, 7, 8}
CONSUME_QUERY = {0, 1, 4, 7, 8}


def log(*args):
    print(time.strftime('%Y-%m-%dT%H:%M:%S'), *args, flush=True)


def rc(s):
    return s.translate(RC_TABLE)[::-1]


def dump(path, value):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_name(p.name + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, default=str) + '\n')
    tmp.replace(p)


def tsv(path, rows):
    rows = list(rows)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w') as f:
        f.write('\t'.join(keys) + '\n')
        for row in rows:
            f.write('\t'.join(json.dumps(row.get(k, ''), separators=(',', ':')) if isinstance(row.get(k), (dict, list, tuple)) else str(row.get(k, '')) for k in keys) + '\n')


def fasta_read(path):
    seqs, name, parts = {}, None, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(parts).upper()
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name is not None:
        seqs[name] = ''.join(parts).upper()
    return seqs


def fasta_write(path, rows):
    with open(path, 'w') as f:
        for name, seq in rows:
            f.write('>' + name + '\n' + seq + '\n')


def checksum(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(2 ** 20), b''):
            h.update(block)
    return h.hexdigest()


def mm_argv(cfg, threads, small=False):
    return [cfg['minimap2'], '-a', '-x', 'map-ont', '-k', str(cfg.get('junction_seed_length',9)) if small else '11', '-w', '5', '-m', '20', '-s', '40', '-n', '2', '-N', '50', '-p', '0.5', '--secondary=yes', '--secondary-seq', '-Y', '--eqx', '--mask-len', '75', '--no-long-join', '-r', '100,200', '-t', str(threads)]


def map_reads(cfg, reference, reads, out, logfile, threads, small=False):
    cmd = mm_argv(cfg, threads, small) + [str(reference), str(reads)]
    log('map', json.dumps(cmd))
    with open(logfile, 'w') as err:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err)
        try:
            with pysam.AlignmentFile(proc.stdout, 'r') as src, pysam.AlignmentFile(out, 'wb', template=src, format_options=[b'level=1']) as dst:
                for a in src:
                    dst.write(a)
        finally:
            if proc.stdout:
                proc.stdout.close()
        code = proc.wait()
        if code:
            raise RuntimeError(f'minimap2 failed ({code}): {logfile}')
    return cmd


def parse_cigar(text):
    codes = dict(zip('MIDNSHP=X', range(9)))
    return [(codes[op], int(n)) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', text)]


def segment(contig, start, cigar, reverse, nm, mapq, route):
    total = sum(n for op, n in cigar if op in CONSUME_QUERY or op == 5)
    left = sum(n for op, n in itertools.takewhile(lambda x: x[0] in (4, 5), cigar))
    right = sum(n for op, n in itertools.takewhile(lambda x: x[0] in (4, 5), reversed(cigar)))
    q0, q1 = left, total - right
    if reverse:
        q0, q1 = total - q1, total - q0
    aln = sum(n for op, n in cigar if op in (0, 1, 2, 7, 8))
    paired = sum(n for op, n in cigar if op in (0, 7, 8))
    return dict(contig=contig, r0=start, r1=start + sum(n for op, n in cigar if op in CONSUME_REF), q0=q0, q1=q1, reverse=reverse, identity=1 - nm / max(1, aln), paired=paired, mapq=mapq, route=route)


def segments_of(a, include_sa=False):
    if a.is_unmapped or not a.cigartuples:
        return []
    result = [segment(a.reference_name, a.reference_start, a.cigartuples, a.is_reverse, a.get_tag('NM') if a.has_tag('NM') else 0, a.mapping_quality, 'alignment')]
    if include_sa and a.has_tag('SA'):
        for item in a.get_tag('SA').strip(';').split(';'):
            if not item:
                continue
            ctg, pos, strand, cig, mq, nm = item.split(',')
            result.append(segment(ctg, int(pos) - 1, parse_cigar(cig), strand == '-', int(nm), int(mq), 'original_SA'))
    return result


def split_events(read_id, segments, p, route, query_windows=False):
    good = {}
    for x in segments:
        if x['paired'] >= p['min_anchor'] and x['identity'] >= p['discovery_identity']:
            good[(x['contig'], x['r0'], x['r1'], x['q0'], x['q1'], x['reverse'])] = x
    segs = sorted(good.values(), key=lambda x: (x['q0'], x['q1']))
    events = []
    for i, x in enumerate(segs):
        for y in segs[i + 1:]:
            gap = y['q0'] - x['q1']
            if gap > p['max_query_gap']:
                break
            if abs(gap) > p['max_query_gap'] or x['contig'] != y['contig'] or x['reverse'] != y['reverse']:
                continue
            # Two alternative descriptions of the same query bases are not a junction.
            if min(x['q1'], y['q0']) - x['q0'] < p['min_anchor'] or y['q1'] - max(y['q0'], x['q1']) < p['min_anchor']:
                continue
            a, b = (x['r0'], y['r1']) if x['reverse'] else (y['r0'], x['r1'])
            length = b - a
            if p['min_length'] <= length <= p['max_length']:
                event=dict(read_id=read_id, contig=x['contig'], a=a, b=b, length=length, gap=gap, identity=min(x['identity'], y['identity']), route=route, event='P')
                if query_windows and gap > 0:
                    event.update(query_window_start=max(x['q0'],x['q1']-p['base_window']),
                                 query_window_end=min(y['q1'],y['q0']+p['base_window']),
                                 query_reverse=x['reverse'])
                events.append(event)
    return events


def cigar_events(a, seqs, p):
    if a.is_unmapped or not a.cigartuples or a.query_sequence is None:
        return []
    seq = seqs.get(a.reference_name)
    if seq is None:
        return []
    q, r, result = 0, a.reference_start, []
    ops = a.cigartuples
    for i, (op, n) in enumerate(ops):
        if op in (1, 2) and p['min_length'] <= n <= p['max_length']:
            left = sum(k for o, k in ops[:i] if o in (0, 7, 8))
            right = sum(k for o, k in ops[i + 1:] if o in (0, 7, 8))
            if min(left, right) >= p['min_anchor']:
                if op == 2:
                    result.append(dict(read_id=a.query_name, contig=a.reference_name, a=r, b=r+n, length=n, gap=0, identity=1.0, route='CIGAR_deletion', event='B'))
                elif n <= p['short_max']:
                    inserted = a.query_sequence[q:q+n]
                    for aa, bb in [(r-n, r), (r, r+n)]:
                        if aa < 0 or bb > len(seq):
                            continue
                        d = edlib.align(inserted, seq[aa:bb], mode='NW', task='distance', k=math.ceil(.15*n))['editDistance']
                        if d >= 0:
                            result.append(dict(read_id=a.query_name, contig=a.reference_name, a=aa, b=bb, length=n, gap=0, identity=1-d/n, route='CIGAR_insertion_copy', event='P'))
        if op in CONSUME_QUERY:
            q += n
        if op in CONSUME_REF:
            r += n
    return result


def export_source(sample, seqs, reads_path, p):
    counts, seen, events, mm_codes = collections.Counter(), set(), [], collections.Counter()
    with pysam.AlignmentFile(sample['bam'], 'rb') as f, open(reads_path, 'w') as out:
        header = f.header.to_dict()
        indexed = {'mapped_records': f.mapped, 'unmapped_records': f.unmapped} if f.has_index() else {}
        for a in f.fetch(until_eof=True):
            counts['records'] += 1
            if a.is_secondary or a.is_supplementary:
                continue
            counts['primary_records'] += 1
            if a.is_qcfail or a.is_duplicate:
                counts['flag_filtered'] += 1
                continue
            if a.query_name in seen:
                counts['duplicate_primary_identifier'] += 1
                continue
            seen.add(a.query_name)
            if not a.query_sequence or any(op == 5 for op, n in a.cigartuples or []):
                counts['missing_or_hardclipped_sequence'] += 1
                continue
            seq = a.get_forward_sequence()
            out.write('>' + a.query_name + '\n' + seq + '\n')
            counts['exported_reads'] += 1
            counts['exported_bases'] += len(seq)
            if a.has_tag('MM'):
                for group in a.get_tag('MM').split(';'):
                    if group:
                        mm_codes[group.split(',')[0]] += 1
            events.extend(split_events(a.query_name, segments_of(a, True), p, 'original_SA'))
            events.extend(cigar_events(a, seqs, p))
    return dict(counts=counts, index=indexed, header=header, modification_groups=mm_codes), events


def discover_fresh(path, seqs, p):
    events = []
    with pysam.AlignmentFile(path, 'rb') as f:
        for name, group in itertools.groupby(f, key=lambda a: a.query_name):
            rows = list(group)
            events.extend(split_events(name, list(itertools.chain.from_iterable(segments_of(a) for a in rows)), p, 'sensitive_remap'))
            for a in rows:
                if not a.is_secondary:
                    events.extend(cigar_events(a, seqs, p))
    return events


def group_events(events, p):
    # Bins accelerate lookup; centre and total range constraints prevent chaining.
    groups, bins = [], collections.defaultdict(list)
    for e in sorted(events, key=lambda x: (x['contig'], x['a'], x['b'], x['read_id'])):
        short = e['length'] < p['short_max']
        tol = p['short_tolerance'] if short else p['long_tolerance']
        ix, iy = e['a'] // tol, e['b'] // tol
        choices = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for k in bins[(e['contig'], short, ix+dx, iy+dy)]:
                    g = groups[k]
                    if abs(e['a']-g['a']) <= tol and abs(e['b']-g['b']) <= tol and max(g['amax'],e['a'])-min(g['amin'],e['a']) <= 2*tol and max(g['bmax'],e['b'])-min(g['bmin'],e['b']) <= 2*tol:
                        choices.append(k)
        if choices:
            k = min(choices, key=lambda k: abs(e['a']-groups[k]['a'])+abs(e['b']-groups[k]['b']))
            g = groups[k]
        else:
            k = len(groups)
            g = dict(contig=e['contig'], a=e['a'], b=e['b'], amin=e['a'], amax=e['a'], bmin=e['b'], bmax=e['b'], events=[])
            groups.append(g)
            bins[(e['contig'], short, ix, iy)].append(k)
        g['events'].append(e)
        for key, value in [('amin',e['a']),('bmin',e['b'])]:
            g[key] = min(g[key],value)
        for key, value in [('amax',e['a']),('bmax',e['b'])]:
            g[key] = max(g[key],value)
    retained = []
    for g in groups:
        ids = {e['read_id'] for e in g['events']}
        if len(ids) < p['min_support']:
            continue
        exact = collections.defaultdict(set)
        rank = {}
        for e in g['events']:
            key = (e['a'],e['b'])
            exact[key].add(e['read_id'])
            rank[key] = min(rank.get(key, 10**6), abs(e['gap']))
        keys = sorted(exact, key=lambda k: (-len(exact[k]),rank[k],k))[:p['breakpoint_hypotheses']]
        g.update(id=f'J{len(retained)+1:04d}', a=keys[0][0], b=keys[0][1], hypotheses=keys, seed_reads=sorted(ids), support=len(ids), routes=dict(collections.Counter(e['route'] for e in g['events'])))
        retained.append(g)
    return retained


def build_references(groups, seqs, p):
    refs, by_seq = [], {}
    def add(sequence, alias):
        if len(sequence) < 2*p['min_anchor']:
            return
        key = (alias['kind'], sequence)
        if key in by_seq:
            refs[by_seq[key]]['aliases'].append(alias)
        else:
            by_seq[key] = len(refs)
            refs.append(dict(id=f'R{len(refs)+1:06d}', seq=sequence, aliases=[alias]))
    for g in groups:
        s = seqs.get(g['contig'])
        if s is None:
            continue
        for hi,(a,b) in enumerate(g['hypotheses']):
            if a < 0 or b > len(s) or b <= a:
                continue
            models = [('no_repeat',a,b,0)]
            outside = [d for d in range(4,min(100,len(s)-b,(b-a)//2)+1) if s[a:a+d]==s[b:b+d] and set(s[a:a+d]) <= set('ACGT')]
            inside = [d for d in range(4,min(100,(b-a)//3)+1) if s[a:a+d]==s[b-d:b] and set(s[a:a+d]) <= set('ACGT')]
            if outside:
                d=max(outside)
                models.append(('repeat_copy_shift',a+d,b+d,d))
            if inside:
                d=max(inside)
                models.extend([('remove_right_repeat',a,b-d,d),('remove_left_repeat',a+d,b,d)])
            for model,u,v,d in models:
                h=min(p['external_flank'],(v-u)//2)
                if h < p['min_anchor']:
                    continue
                f0=max(0,a-p['external_flank']);f1=min(len(s),b+p['external_flank'])
                common=dict(locus=g['id'],contig=g['contig'],a=a,b=b,model=model,hypothesis=hi,repeat_length=d,h=h)
                def alias(kind, arms):
                    return dict(common,kind=kind,join=arms[0][2]-arms[0][1],arms=[dict(role=role,start=start,end=end) for role,start,end in arms])
                # Occupied arms always preserve the genomic junctions at a and b.
                for kind,arms in [('L',[('host_left',f0,a),('head',a,a+h)]),('R',[('tail',b-h,b),('host_right',b,f1)]),('P',[('tail',v-h,v),('head',u,u+h)]),('B',[('host_left',f0,a),('host_right',b,f1)])]:
                    if min(end-start for role,start,end in arms) < p['min_anchor']:
                        continue
                    add(''.join(s[start:end] for role,start,end in arms),alias(kind,arms))
                if inside and model=='no_repeat':
                    d=max(inside)
                    arms=[('host_left',f0,a+d),('host_right',b,f1)]
                    al=alias('B',arms);al['model']='B_retained_left_repeat';al['repeat_length']=d
                    add(''.join(s[start:end] for role,start,end in arms),al)
    return refs


def alignment_metric_reference(a, reference, alias, window, p):
    c=alias['join']; n=len(reference)
    w=min(window,c,n-c)
    if w < p['min_anchor']:
        return None
    lo,hi=c-w,c+w
    pairs=[(q,r) for q,r in a.get_aligned_pairs(matches_only=True) if lo<=r<hi]
    sides=[sum(r<c for q,r in pairs),sum(r>=c for q,r in pairs)]
    if min(sides) < math.ceil(p['paired_fraction']*w):
        return None
    core=min(25,w)
    if sum(c-core<=r<c+core for q,r in pairs) < math.ceil(p['paired_fraction']*2*core):
        return None
    ins=dele=0;qpos=0;rpos=a.reference_start;core_gap=0
    for op,k in a.cigartuples:
        if op == 1 and lo<=rpos<hi:
            ins+=k
            if c-core<=rpos<c+core:core_gap=max(core_gap,k)
        if op in (2,3):
            dele+=max(0,min(hi,rpos+k)-max(lo,rpos))
            if max(rpos,c-core)<min(rpos+k,c+core):core_gap=max(core_gap,k)
        if op in CONSUME_QUERY:qpos+=k
        if op in CONSUME_REF:rpos+=k
    matches=sum(a.query_sequence[q].upper()==reference[r] for q,r in pairs)
    identity=matches/max(1,len(pairs)+ins+dele)
    if identity < p['accept_identity'] or core_gap > p['max_core_indel']:
        return None
    q0=min(q for q,r in pairs);q1=max(q for q,r in pairs)+1
    length=a.query_length
    f0,f1=(length-q1,length-q0) if a.is_reverse else (q0,q1)
    return dict(window=w,paired_left=sides[0],paired_right=sides[1],identity=identity,core_gap=core_gap,q0=q0,q1=q1,forward_q0=f0,forward_q1=f1,query_slice=a.query_sequence[q0:q1])


def alignment_metric(a, reference, alias, window, p):
    c=alias['join'];w=min(window,c,len(reference)-c)
    if w<p['min_anchor']:return None
    need=math.ceil(p['paired_fraction']*w)
    start,end=a.reference_start,a.reference_end
    upper_left=max(0,min(c,end)-max(c-w,start))
    upper_right=max(0,min(c+w,end)-max(c,start))
    if min(upper_left,upper_right)<need or a.query_alignment_length<2*need:return None
    return alignment_metric_reference(a,reference,alias,window,p)


def recruit_reads(fresh_bam, groups, out_path, flank):
    by_contig=collections.defaultdict(list)
    seeds=set()
    for g in groups:
        by_contig[g['contig']].extend([g['a'],g['b']]);seeds.update(g['seed_reads'])
    for values in by_contig.values():values.sort()
    count=0
    with pysam.AlignmentFile(fresh_bam,'rb') as f,open(out_path,'w') as out:
        for rid,group in itertools.groupby(f,key=lambda a:a.query_name):
            rows=list(group);keep=rid in seeds
            for a in rows:
                if a.is_unmapped:continue
                points=by_contig.get(a.reference_name,[])
                k=bisect.bisect_left(points,a.reference_start-flank)
                if k<len(points) and points[k]<=a.reference_end+flank:keep=True
            if keep:
                a=next((a for a in rows if a.query_sequence and not a.is_secondary and not a.is_supplementary),None)
                if a is None:a=next((a for a in rows if a.query_sequence),None)
                if a is not None:
                    out.write('>'+rid+'\n'+a.get_forward_sequence()+'\n');count+=1
    return count


def validate_synthetic_reference(cfg, refs, ref_bam, p):
    index={r['id']:r for r in refs}
    reverse_refs=[rc(r['seq']) for r in refs]
    postings=None;kmer_size=11
    if cfg.get('indexed_competition'):
        postings=collections.defaultdict(list)
        for j,ref in enumerate(refs):
            for word in {ref['seq'][k:k+kmer_size] for k in range(max(0,len(ref['seq'])-kmer_size+1))}:
                postings[word].append(j)
        log('lossless competitor index ready',len(refs),'references',len(postings),'distinct 11-mers')
    calls=[];rejections=collections.Counter();query_cache={}
    with pysam.AlignmentFile(ref_bam,'rb') as f:
        for rid,group in itertools.groupby(f,key=lambda a:a.query_name):
            if cfg.get('indexed_competition'):query_cache={}
            rows=list(group);possible=[]
            for a in rows:
                if a.is_unmapped or not a.query_sequence:continue
                ref=index[a.reference_name]
                for al in ref['aliases']:
                    window=min(p['base_window'],al['h'])
                    m=alignment_metric(a,ref['seq'],al,window,p)
                    if m is None:continue
                    qs=m['query_slice'];qskey=min(qs,rc(qs))
                    key=(qskey,al['kind'],ref['id'])
                    if key not in query_cache:
                        current=edlib.align(qs,ref['seq'],mode='HW',task='distance')['editDistance']
                        margin=max(p['min_edit_margin'],math.ceil(p['edit_margin_fraction']*len(qs)))
                        best_other=10**9;other_loci=set()
                        eligible=range(len(refs))
                        if postings is not None:
                            bound=current+margin
                            block=len(qs)//(bound+1)
                            if block>=kmer_size:
                                found=set()
                                # <=bound edits can damage at most bound disjoint
                                # query blocks. One of bound+1 exact 11-mers must
                                # survive. Include both reference orientations.
                                for j in range(bound+1):
                                    word=qs[j*block:j*block+kmer_size]
                                    found.update(postings.get(word,()))
                                    found.update(postings.get(rc(word),()))
                                eligible=sorted(found)
                        for competitor_index in eligible:
                            competitor=refs[competitor_index]
                            # Compare the same complete query slice, including its RC.
                            limit=current+margin
                            ds=[edlib.align(qs,z,mode='HW',task='distance',k=limit)['editDistance'] for z in (competitor['seq'],reverse_refs[competitor_index])]
                            ds=[v for v in ds if v>=0]
                            if not ds:continue
                            d=min(ds)
                            kinds={x['kind'] for x in competitor['aliases']}
                            equivalent=kinds <= ({'L','R'} if al['kind'] in ('L','R') else {al['kind']})
                            if equivalent:
                                if d<=current+margin:other_loci.update(x['locus'] for x in competitor['aliases'])
                            else:best_other=min(best_other,d)
                        query_cache[key]=(current,margin,best_other,sorted(other_loci))
                    current,margin,best_other,other_loci=query_cache[key]
                    if best_other-current<margin:
                        rejections['indistinguishable_configuration']+=1;continue
                    paired=[];genomic_query_positions=[]
                    for q,r in a.get_aligned_pairs(matches_only=True):
                        offset=0
                        for arm in al['arms']:
                            length=arm['end']-arm['start']
                            if offset<=r<offset+length:
                                near=(al['join']-r if r<al['join'] else r-al['join'])
                                qf=a.query_length-1-q if a.is_reverse else q
                                genomic_r=arm['start']+r-offset
                                if m['forward_q0']<=qf<m['forward_q1']:genomic_query_positions.append(genomic_r)
                                if near<=min(p['mod_window'],al['h']):
                                    keep=True
                                    if cfg.get('sparse_modification_pairs'):
                                        canonical=rc(a.query_sequence[q]) if a.is_reverse else a.query_sequence[q]
                                        keep=canonical=='A' and cfg.get('keep_A_modification_pairs',True)
                                        if canonical=='C':
                                            rp=r-offset
                                            keep=(rp==0 or ref['seq'][r-1:r+1]=='CG') if a.is_reverse else (rp==length-1 or ref['seq'][r:r+2]=='CG')
                                    if keep:paired.append((qf,genomic_r,arm['role']))
                                break
                            offset+=length
                    tiers={}
                    for w in (500,1000):
                        tiers[str(w)]=None if min(al['h'],al['join'],len(ref['seq'])-al['join'])<w else alignment_metric(a,ref['seq'],al,w,p) is not None
                    possible.append(dict(read_id=rid,ref_id=ref['id'],locus=al['locus'],contig=al['contig'],kind=al['kind'],model=al['model'],a=al['a'],b=al['b'],h=al['h'],repeat_length=al['repeat_length'],identity=m['identity'],window=m['window'],forward_q0=m['forward_q0'],forward_q1=m['forward_q1'],genomic_query_min=min(genomic_query_positions),genomic_query_max=max(genomic_query_positions)+1,query_slice=qs,edit_distance=current,required_margin=margin,competing_distance=None if best_other==10**9 else best_other,competing_loci=other_loci,reverse=a.is_reverse,pairs=paired,tiers=tiers))
            # Collapse repeat-copy/breakpoint alternatives of the same event and slice.
            unique={}
            for call in sorted(possible,key=lambda x:(x['edit_distance'],-x['identity'],x['repeat_length'],x['model'])):
                key=(call['locus'],call['kind'])
                if key not in unique:unique[key]=call
            calls.extend(unique.values())
    return calls,dict(rejections)


def prepare_ranked_index(cfg, refs):
    index={r['id']:r for r in refs}
    reverse_refs=[rc(r['seq']) for r in refs]
    postings=None;kmer_size=11
    if cfg.get('indexed_competition'):
        postings=collections.defaultdict(list)
        for j,ref in enumerate(refs):
            for word in {ref['seq'][k:k+kmer_size] for k in range(max(0,len(ref['seq'])-kmer_size+1))}:
                postings[word].append(j)
        log('lossless competitor index ready',len(refs),'references',len(postings),'distinct 11-mers')
    return index,reverse_refs,postings,kmer_size


def bounded_bam_records(handle, stop):
    while stop is None or handle.tell()<stop:
        try:yield next(handle)
        except StopIteration:return


def validate_synthetic_ranked(cfg, refs, ref_bam, p, prepared=None, start=None, stop=None):
    index,reverse_refs,postings,kmer_size=prepared if prepared is not None else prepare_ranked_index(cfg,refs)
    calls=[];rejections=collections.Counter();nreads=0
    with pysam.AlignmentFile(ref_bam,'rb') as f:
        if start is not None:f.seek(start)
        records=bounded_bam_records(f,stop) if stop is not None else f
        for rid,group in itertools.groupby(records,key=lambda a:a.query_name):
            rows=list(group);candidates=[];distance_cache={};query_cache={};unique={}
            for a in rows:
                if a.is_unmapped or not a.query_sequence:continue
                ref=index[a.reference_name]
                for al in ref['aliases']:
                    window=min(p['base_window'],al['h'])
                    m=alignment_metric(a,ref['seq'],al,window,p)
                    if m is None:continue
                    qs=m['query_slice'];qskey=min(qs,rc(qs));key=(qskey,al['kind'],ref['id'])
                    if key not in distance_cache:distance_cache[key]=edlib.align(qs,ref['seq'],mode='HW',task='distance')['editDistance']
                    current=distance_cache[key]
                    rank=(current,-m['identity'],al['repeat_length'],al['model'],len(candidates))
                    candidates.append((rank,a,ref,al,m,current))
            for rank,a,ref,al,m,current in sorted(candidates,key=lambda x:x[0]):
                group_key=(al['locus'],al['kind'])
                if group_key in unique:
                    rejections['duplicate_hypotheses_skipped']+=1;continue
                qs=m['query_slice'];qskey=min(qs,rc(qs));key=(qskey,al['kind'],ref['id'])
                if key not in query_cache:
                    margin=max(p['min_edit_margin'],math.ceil(p['edit_margin_fraction']*len(qs)))
                    best_other=10**9;other_loci=set()
                    eligible=range(len(refs))
                    if postings is not None:
                        bound=current+margin
                        block=len(qs)//(bound+1)
                        if block>=kmer_size:
                            found=set()
                            # <=bound edits can damage at most bound disjoint
                            # query blocks. One of bound+1 exact 11-mers must
                            # survive. Include both reference orientations.
                            for j in range(bound+1):
                                word=qs[j*block:j*block+kmer_size]
                                found.update(postings.get(word,()))
                                found.update(postings.get(rc(word),()))
                            eligible=sorted(found)
                    for competitor_index in eligible:
                        competitor=refs[competitor_index]
                        # Compare the same complete query slice, including its RC.
                        limit=current+margin
                        ds=[edlib.align(qs,z,mode='HW',task='distance',k=limit)['editDistance'] for z in (competitor['seq'],reverse_refs[competitor_index])]
                        ds=[v for v in ds if v>=0]
                        if not ds:continue
                        d=min(ds)
                        kinds={x['kind'] for x in competitor['aliases']}
                        equivalent=kinds <= ({'L','R'} if al['kind'] in ('L','R') else {al['kind']})
                        if equivalent:
                            if d<=current+margin:other_loci.update(x['locus'] for x in competitor['aliases'])
                        else:best_other=min(best_other,d)
                    query_cache[key]=(current,margin,best_other,sorted(other_loci))
                current,margin,best_other,other_loci=query_cache[key]
                if best_other-current<margin:
                    rejections['indistinguishable_configuration']+=1;continue
                paired=[];genomic_query_positions=[]
                for q,r in a.get_aligned_pairs(matches_only=True):
                    offset=0
                    for arm in al['arms']:
                        length=arm['end']-arm['start']
                        if offset<=r<offset+length:
                            near=(al['join']-r if r<al['join'] else r-al['join'])
                            qf=a.query_length-1-q if a.is_reverse else q
                            genomic_r=arm['start']+r-offset
                            if m['forward_q0']<=qf<m['forward_q1']:genomic_query_positions.append(genomic_r)
                            if near<=min(p['mod_window'],al['h']):
                                keep=True
                                if cfg.get('sparse_modification_pairs'):
                                    canonical=rc(a.query_sequence[q]) if a.is_reverse else a.query_sequence[q]
                                    keep=canonical=='A' and cfg.get('keep_A_modification_pairs',True)
                                    if canonical=='C':
                                        rp=r-offset
                                        keep=(rp==0 or ref['seq'][r-1:r+1]=='CG') if a.is_reverse else (rp==length-1 or ref['seq'][r:r+2]=='CG')
                                if keep:paired.append((qf,genomic_r,arm['role']))
                            break
                        offset+=length
                tiers={}
                for w in (500,1000):
                    tiers[str(w)]=None if min(al['h'],al['join'],len(ref['seq'])-al['join'])<w else alignment_metric(a,ref['seq'],al,w,p) is not None
                unique[group_key]=dict(read_id=rid,ref_id=ref['id'],locus=al['locus'],contig=al['contig'],kind=al['kind'],model=al['model'],a=al['a'],b=al['b'],h=al['h'],repeat_length=al['repeat_length'],identity=m['identity'],window=m['window'],forward_q0=m['forward_q0'],forward_q1=m['forward_q1'],genomic_query_min=min(genomic_query_positions),genomic_query_max=max(genomic_query_positions)+1,query_slice=qs,edit_distance=current,required_margin=margin,competing_distance=None if best_other==10**9 else best_other,competing_loci=other_loci,reverse=a.is_reverse,pairs=paired,tiers=tiers)
            calls.extend(unique.values());nreads+=1
            if nreads%5000==0 and not cfg.get('_chunk_worker'):log('ranked validation progress',nreads,'read groups;',len(calls),'retained calls')
    if '_expected_read_groups' in cfg:assert nreads==cfg['_expected_read_groups'],'Chunk read-group count mismatch'
    return calls,dict(rejections)


_FORK_VALIDATION_CONTEXT=None


def checkpoint_identity(cfg,refs,ref_bam,p):
    path=Path(ref_bam).resolve();st=path.stat();h=hashlib.sha256()
    for ref in refs:
        h.update(ref['id'].encode());h.update(b'\0');h.update(ref['seq'].encode());h.update(b'\0')
        h.update(json.dumps(ref['aliases'],sort_keys=True).encode());h.update(b'\0')
    return dict(code_sha256=checksum(__file__),bam=str(path),bam_bytes=st.st_size,bam_mtime_ns=st.st_mtime_ns,references_sha256=h.hexdigest(),parameters=p,indexed_competition=bool(cfg.get('indexed_competition')),sparse_modification_pairs=bool(cfg.get('sparse_modification_pairs')),keep_A_modification_pairs=bool(cfg.get('keep_A_modification_pairs',True)))


def partition_bam(ref_bam,nchunks,identity,path):
    path=Path(path)
    if path.exists():
        saved=json.loads(path.read_text());assert saved['identity']==identity,'Checkpoint input identity changed'
        return saved['chunks']
    target=max(65536,Path(ref_bam).stat().st_size//max(1,nchunks))
    chunks=[];groups=records=0;previous=None;first_read=last_read=None
    with pysam.AlignmentFile(ref_bam,'rb') as f:
        start=f.tell()
        while True:
            offset=f.tell()
            try:a=next(f)
            except StopIteration:break
            new_group=a.query_name!=previous
            if new_group and groups and (offset>>16)-(start>>16)>=target:
                chunks.append(dict(index=len(chunks),start=start,stop=offset,read_groups=groups,records=records,first_read=first_read,last_read=last_read))
                start=offset;groups=records=0;first_read=last_read=None
            if new_group:
                groups+=1
                if first_read is None:first_read=a.query_name
                last_read=a.query_name
            previous=a.query_name;records+=1
        if groups:chunks.append(dict(index=len(chunks),start=start,stop=None,read_groups=groups,records=records,first_read=first_read,last_read=last_read))
    assert chunks and chunks[0]['start']>=0
    for left,right in zip(chunks,chunks[1:]):
        assert left['stop']==right['start'] and left['last_read']!=right['first_read']
    dump(path,dict(identity=identity,chunks=chunks,read_groups=sum(c['read_groups'] for c in chunks),alignment_records=sum(c['records'] for c in chunks)))
    log('BAM partition ready',len(chunks),'chunks;',sum(c['read_groups'] for c in chunks),'whole read groups')
    return chunks


def validation_chunk_worker(chunk):
    cfg,refs,bam,p,prepared,identity,folder=_FORK_VALIDATION_CONTEXT
    path=Path(folder)/('chunk_%04d.json.gz'%chunk['index'])
    if path.exists():
        with gzip.open(path,'rt') as f:saved=json.load(f)
        assert saved['identity']==identity and saved['chunk']==chunk,'Chunk checkpoint identity mismatch'
        return dict(index=chunk['index'],path=str(path),read_groups=chunk['read_groups'],calls=len(saved['calls']),reused=True,elapsed_seconds=0)
    started=time.monotonic()
    local=dict(cfg,_chunk_worker=True,_expected_read_groups=chunk['read_groups'])
    calls,rejections=validate_synthetic_ranked(local,refs,bam,p,prepared,chunk['start'],chunk['stop'])
    tmp=path.with_name(path.name+'.tmp')
    with gzip.open(tmp,'wt',compresslevel=1) as f:json.dump(dict(identity=identity,chunk=chunk,calls=calls,rejections=rejections),f,separators=(',',':'))
    tmp.replace(path)
    return dict(index=chunk['index'],path=str(path),read_groups=chunk['read_groups'],calls=len(calls),reused=False,elapsed_seconds=round(time.monotonic()-started,3))


def validate_synthetic_parallel(cfg,refs,ref_bam,p):
    global _FORK_VALIDATION_CONTEXT
    workers=int(cfg['validation_workers']);folder=Path(cfg['validation_checkpoint_dir']);folder.mkdir(parents=True,exist_ok=True)
    identity=checkpoint_identity(cfg,refs,ref_bam,p)
    chunks=partition_bam(ref_bam,cfg.get('validation_chunks',workers*4),identity,folder/'partitions.json')
    workers=min(workers,len(chunks));prepared=prepare_ranked_index(cfg,refs)
    _FORK_VALIDATION_CONTEXT=(cfg,refs,str(ref_bam),p,prepared,identity,str(folder))
    finished=[];started=time.monotonic()
    try:
        with multiprocessing.get_context('fork').Pool(workers) as pool:
            for result in pool.imap_unordered(validation_chunk_worker,chunks,chunksize=1):
                finished.append(result)
                log('parallel validation checkpoint',len(finished),'of',len(chunks),'chunks;',sum(v['read_groups'] for v in finished),'read groups;',sum(v['calls'] for v in finished),'calls;',workers,'workers')
    finally:_FORK_VALIDATION_CONTEXT=None
    calls=[];rejections=collections.Counter()
    for result in sorted(finished,key=lambda r:r['index']):
        with gzip.open(result['path'],'rt') as f:record=json.load(f)
        assert record['identity']==identity and record['chunk']==chunks[result['index']]
        for c in record['calls']:c['pairs']=[tuple(v) for v in c['pairs']]
        calls.extend(record['calls']);rejections.update(record['rejections'])
    dump(folder/'summary.json',dict(status='complete',workers=workers,chunks=len(chunks),read_groups=sum(c['read_groups'] for c in chunks),alignment_records=sum(c['records'] for c in chunks),calls=len(calls),reused_chunks=sum(v['reused'] for v in finished),elapsed_seconds=round(time.monotonic()-started,3),chunks_detail=sorted(finished,key=lambda r:r['index']),identity=identity))
    return calls,dict(rejections)


def validate_synthetic(cfg, refs, ref_bam, p):
    if cfg.get('ranked_candidate_search'):
        if cfg.get('validation_workers',1)>1:return validate_synthetic_parallel(cfg,refs,ref_bam,p)
        return validate_synthetic_ranked(cfg,refs,ref_bam,p)
    return validate_synthetic_reference(cfg,refs,ref_bam,p)


def background_check(cfg, sample, calls, seqs, out, threads):
    fasta=out/'background_queries.fasta';bam=out/'background_queries.bam'
    # Each validated slice has an independent identifier: no full-read score bias.
    for i,c in enumerate(calls):c['call_id']=f'C{i+1:07d}';c['background_competitors']=[]
    if not calls:return []
    fasta_write(fasta,((c['call_id'],c['query_slice']) for c in calls))
    cmd=map_reads(cfg,sample['fasta'],fasta,bam,out/'background_mapping.log',threads,True)
    by_id={c['call_id']:c for c in calls}
    with pysam.AlignmentFile(bam,'rb') as f:
        for a in f:
            if a.is_unmapped:continue
            c=by_id[a.query_name];qlen=len(c['query_slice'])
            if a.reference_end-a.reference_start>2*qlen+50:continue
            if a.query_alignment_length<.9*qlen:continue
            ref=seqs[a.reference_name][max(0,a.reference_start-30):a.reference_end+30]
            qs=a.query_sequence
            limit=c['edit_distance']+c['required_margin']
            d=edlib.align(qs,ref,mode='HW',task='distance',k=limit)['editDistance']
            if d>=0 and d-c['edit_distance']<c['required_margin']:
                c['background_competitors'].append(dict(contig=a.reference_name,start=a.reference_start,end=a.reference_end,edits=d,mapq=a.mapping_quality))
    for c in calls:
        other=set(c['competing_loci'])-{c['locus']}
        c['locus_ambiguous']=bool(other)
        c['genomic_contiguous_alternative']=c['kind'] in ('P','B') and bool(c['background_competitors'])
        c['simple']=not c['locus_ambiguous'] and not c['genomic_contiguous_alternative']
    by_read=collections.defaultdict(list)
    for c in calls:by_read[(c['locus'],c['read_id'])].append(c)
    for key,values in by_read.items():
        kinds={c['kind'] for c in values}
        complex_read=('P' in kinds and bool(kinds-{'P'})) or ('B' in kinds and bool(kinds-{'B'}))
        for c in values:c['complex_read']=complex_read;c['simple']=c['simple'] and not complex_read
    return cmd


def input_provenance(sample):
    result={}
    for key in ('bam','fasta','assembly_info'):
        path=Path(sample[key]);st=path.stat()
        result[key]=dict(path=str(path),bytes=st.st_size,mtime_ns=st.st_mtime_ns)
        if key!='bam':result[key]['sha256']=checksum(path)
    return result


def audit_source_copies(cfg,sample,calls,origins,seqs,out,threads,p):
    """Find other possible source loci even if they produced no discovery seeds.

    A P-spanning query cannot map contiguously to an occupied genomic copy.
    Consequently the whole-junction background check alone is insufficient:
    map its two reference arms separately, pair their possible genomic origins,
    and then compare the SAME observed query slice to each reconstructed P.
    """
    aliases={}
    for ref in origins:
        for al in ref['aliases']:
            aliases[(ref['id'],al['locus'],al['kind'],al['model'],al['a'],al['b'])]=al
    arm_queries={};call_keys={};arm_info={}
    for i,c in enumerate(calls):
        if c['kind']!='P':continue
        key=(c['ref_id'],c['locus'],c['kind'],c['model'],c['a'],c['b'])
        al=aliases[key];arms=al['arms'];w=min(c['window'],al['h'])
        signature=(al['contig'],tuple((x['start'],x['end'],x['role']) for x in arms),w)
        if signature not in arm_queries:
            tag=f'U{len(arm_queries)+1:06d}';arm_queries[signature]=tag
            for role,arm in zip(('tail','head'),arms):
                start,end=(arm['end']-w,arm['end']) if role=='tail' else (arm['start'],arm['start']+w)
                arm_info[tag+'_'+role]=dict(tag=tag,role=role,seq=seqs[al['contig']][start:end],w=w)
        call_keys[i]=(arm_queries[signature],al)
    if not arm_info:return dict(tested_P_calls=0,ambiguous_P_calls=0,commands=[])
    out=Path(out);fa=out/'source_arms.fasta';bam=out/'source_arms.bam'
    fasta_write(fa,((name,x['seq']) for name,x in arm_info.items()))
    cmd=mm_argv(cfg,threads,True);cmd[cmd.index('-N')+1]='200';cmd += [sample['fasta'],str(fa)]
    with open(out/'source_arms_mapping.log','w') as err:
        proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=err)
        with pysam.AlignmentFile(proc.stdout,'r') as src,pysam.AlignmentFile(bam,'wb',template=src) as dst:
            for a in src:dst.write(a)
        proc.stdout.close()
        if proc.wait():raise RuntimeError('Source-arm mapping failed')
    hits=collections.defaultdict(list)
    with pysam.AlignmentFile(bam,'rb') as f:
        for a in f:
            if a.is_unmapped:continue
            info=arm_info[a.query_name];w=info['w']
            if a.query_alignment_length<.9*w:continue
            pairs=a.get_aligned_pairs(matches_only=True)
            if not pairs:continue
            seq=seqs[a.reference_name];lo=max(0,a.reference_start-15);hi=min(len(seq),a.reference_end+15)
            q=a.query_sequence
            d=edlib.align(q,seq[lo:hi],mode='HW',task='distance',k=math.ceil(.05*w))['editDistance']
            if d<0:continue
            # Reference arms have no sequencing error. Permit only small endpoint
            # extrapolation, and retain secondary/supplementary placements.
            start=pairs[0][1]-pairs[0][0]
            end=pairs[-1][1]+1+(w-1-pairs[-1][0])
            near_start=pairs[0][0]<=10;near_end=w-1-pairs[-1][0]<=10
            if info['role']=='tail':
                edge=start if a.is_reverse else end
                covered=near_start if a.is_reverse else near_end
            else:
                edge=end if a.is_reverse else start
                covered=near_end if a.is_reverse else near_start
            if not covered:continue
            value=dict(contig=a.reference_name,edge=edge,reverse=a.is_reverse,edits=d)
            if value not in hits[a.query_name]:hits[a.query_name].append(value)
    audit=[]
    for i,(tag,al) in call_keys.items():
        c=calls[i];w=min(c['window'],al['h']);query=c['query_slice']
        candidates=[];seen=set()
        for tail in hits[tag+'_tail']:
            for head in hits[tag+'_head']:
                if tail['contig']!=head['contig'] or tail['reverse']!=head['reverse']:continue
                a,b=(tail['edge'],head['edge']) if tail['reverse'] else (head['edge'],tail['edge'])
                if b-a<p['min_length'] or b-a>p['max_length']:continue
                s=seqs[tail['contig']]
                if a<0 or b>len(s) or a+w>b-w:continue
                # The seed's repeat-equivalent boundary shifts are one source.
                if tail['contig']==al['contig'] and abs(a-al['arms'][1]['start'])<=max(50,al['repeat_length']+25) and abs(b-al['arms'][0]['end'])<=max(50,al['repeat_length']+25):continue
                signature=(tail['contig'],a,b,tail['reverse'])
                if signature in seen:continue
                seen.add(signature)
                alternative=s[b-w:b]+s[a:a+w]
                if tail['reverse']:alternative=rc(alternative)
                distance=edlib.align(query,alternative,mode='HW',task='distance',k=c['edit_distance']+c['required_margin'])['editDistance']
                if distance>=0 and distance-c['edit_distance']<c['required_margin']:
                    candidates.append(dict(contig=tail['contig'],a=a,b=b,length=b-a,reverse=tail['reverse'],edit_distance=distance))
        saturated=len(hits[tag+'_tail'])>=200 or len(hits[tag+'_head'])>=200
        c['source_copy_alternatives']=candidates[:50];c['source_search_saturated']=saturated
        c['source_copy_ambiguous']=bool(candidates) or saturated
        if c['source_copy_ambiguous']:
            c['locus_ambiguous']=True;c['simple']=False
        audit.append(dict(read_id=c['read_id'],locus=c['locus'],source_copy_ambiguous=c['source_copy_ambiguous'],saturated=saturated,alternative_count=len(candidates),alternatives=candidates[:50]))
    dump(out/'source_copy_audit.json',audit)
    tsv(out/'source_copy_audit.tsv',audit)
    return dict(tested_P_calls=len(audit),ambiguous_P_calls=sum(r['source_copy_ambiguous'] for r in audit),saturated_calls=sum(r['saturated'] for r in audit),commands=[cmd])


def audit_existing(cfg,sample,out,threads):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    old=Path(cfg['reuse_results'])/sample['sample']
    deadline=time.monotonic()+3600
    if not (old/'summary.json').exists():log('waiting for original approved sample stage',sample['sample'])
    while not (old/'summary.json').exists():
        if time.monotonic()>deadline:raise TimeoutError('Original approved sample stage did not finish within one hour')
        time.sleep(5)
    controls_summary=json.loads((Path(cfg['reuse_results'])/'controls/summary.json').read_text())
    assert controls_summary['status']=='passed'
    previous=json.loads((old/'summary.json').read_text())
    assert previous['status']=='complete'
    seqs=fasta_read(sample['fasta']);calls=json.loads((old/'calls.json').read_text());groups=json.loads((old/'discovery_groups.json').read_text());origins=json.loads((old/'reference_origins.json').read_text())
    audit=audit_source_copies(cfg,sample,calls,origins,seqs,out,threads,cfg['parameters'])
    if cfg.get('normalize_occupied_labels'):
        audit['occupied_labels_normalized']=normalize_occupied_labels(calls,origins)
    loci=loci_summary(groups,calls,sample['sample'])
    # Add assembly context without converting Flye's circular flag to a topology assay.
    flye={}
    for line in Path(sample['assembly_info']).read_text().splitlines():
        if line.startswith('#'):continue
        parts=line.split()
        if len(parts)>=4:flye[parts[0]]=parts
    for row in loci:
        n=len(seqs[row['contig']]);row['contig_length']=n
        row['near_contig_ends']=row['a']<=50 and n-row['b']<=50
        row['flye_circular']=flye.get(row['contig'],['','','','unknown'])[3]
    mod_summary,contrasts=modifications(sample,calls,seqs,out,cfg['parameters'])
    dump(out/'calls.json',calls);dump(out/'loci.json',loci);tsv(out/'loci.tsv',loci)
    dump(out/'reference_origins.json',origins)
    summary=dict(previous,version=VERSION,copy_audit=audit,modifications=mod_summary,loci=loci,source_results=str(old),software=software(cfg))
    for kind in ('P','B'):
        summary['confirmed_simple_'+kind]=sum(r['confirmed_simple_'+kind] for r in loci)
        summary['confirmed_short_'+kind]=sum(r['confirmed_simple_'+kind] and r['short'] for r in loci)
    dump(out/'summary.json',summary)
    log('copy audit complete',sample['sample'],audit)


def normalize_occupied_labels(calls,origins):
    """A genomic occupied read may support several overlapping locus tests.

    Distinguish a second label for the SAME genomic query placement from a
    different possible genomic copy. Only the latter is source ambiguity.
    P/B assignments and actual duplicate-copy ambiguity are unchanged.
    """
    placements=collections.defaultdict(list)
    for ref in origins:
        for al in ref['aliases']:
            if al['kind'] in ('L','R') and al['arms'][0]['end']==al['arms'][1]['start']:
                placements[al['locus']].append((al['contig'],al['arms'][0]['start'],al['arms'][1]['end']))
    restored=0
    for c in calls:
        if c['kind'] not in ('L','R'):continue
        if 'genomic_query_min' in c:lo,hi=c['genomic_query_min'],c['genomic_query_max']
        else:
            positions=[r for q,r,role in c['pairs'] if c['forward_q0']<=q<c['forward_q1']]
            if not positions:continue
            lo,hi=min(positions),max(positions)+1
        distinct=[];same=[]
        for locus in set(c['competing_loci'])-{c['locus']}:
            if any(ctg==c['contig'] and start<=lo and end>=hi for ctg,start,end in placements[locus]):same.append(locus)
            else:distinct.append(locus)
        was_simple=c['simple']
        c['overlapping_occupied_locus_labels']=sorted(same)
        c['distinct_source_locus_labels']=sorted(distinct)
        c['locus_ambiguous']=bool(distinct)
        c['simple']=not c['locus_ambiguous'] and not c['complex_read'] and not c['genomic_contiguous_alternative']
        restored+=c['simple'] and not was_simple
    return restored


def software(cfg):
    result=dict(pipeline=VERSION,python=sys.version,code_sha256=checksum(__file__))
    result['packages']={x:importlib.metadata.version(x) for x in ('pysam','numpy','edlib','scipy')}
    for name in ('minimap2','samtools'):
        r=subprocess.run([cfg[name],'--version'],capture_output=True,text=True,check=True)
        result[name]=r.stdout.splitlines()[0]
    return result


def loci_summary(groups,calls,sample):
    result=[];by_locus=collections.defaultdict(list)
    for c in calls:by_locus[c['locus']].append(c)
    for g in groups:
        cc=by_locus[g['id']]
        pbest=collections.Counter((c['a'],c['b'],c['model'],c['repeat_length']) for c in cc if c['kind']=='P' and c['simple'])
        if pbest:a,b,model,d=pbest.most_common(1)[0][0]
        else:a,b,model,d=g['a'],g['b'],'no_simple_P',0
        row=dict(sample=sample,locus=g['id'],contig=g['contig'],a=a,b=b,length=b-a,short=b-a<2000,model=model,repeat_length=d,discovery_reads=g['support'],discovery_routes=g['routes'])
        for kind in ('L','R','P','B'):
            row[kind+'_all']=len({c['read_id'] for c in cc if c['kind']==kind})
            row[kind+'_simple']=len({c['read_id'] for c in cc if c['kind']==kind and c['simple']})
            row[kind+'_500']=len({c['read_id'] for c in cc if c['kind']==kind and c['simple'] and c['tiers']['500'] is True})
            row[kind+'_1000']=len({c['read_id'] for c in cc if c['kind']==kind and c['simple'] and c['tiers']['1000'] is True})
        row['complex_reads']=len({c['read_id'] for c in cc if c['complex_read']})
        row['ambiguous_reads']=len({c['read_id'] for c in cc if c['locus_ambiguous']})
        row['genomic_alternative_reads']=len({c['read_id'] for c in cc if c['genomic_contiguous_alternative']})
        row['confirmed_simple_P']=row['P_simple']>=3
        row['confirmed_simple_B']=row['B_simple']>=3
        row['P_read_ids']=sorted({c['read_id'] for c in cc if c['kind']=='P' and c['simple']})
        row['B_read_ids']=sorted({c['read_id'] for c in cc if c['kind']=='B' and c['simple']})
        result.append(row)
    return result


def mod_values(a):
    seq=a.get_forward_sequence().upper()
    explicit=collections.defaultdict(dict);modes={}
    if not a.has_tag('MM') or not a.has_tag('ML'):return seq,explicit,modes
    for group in a.get_tag('MM').split(';'):
        if not group:continue
        header=group.split(',')[0]
        match=re.fullmatch(r'([ACGTUN])([+-])([a-z]+|\d+)([.?]?)',header)
        if match:
            canonical,strand,codes,mode=match.groups()
            if strand=='+':
                for code in codes:modes[(canonical,code)]=mode or '.'
    for (canonical,strand,code),values in (a.modified_bases_forward or {}).items():
        if strand!=0:continue
        for pos,ml in values:
            if ml>=0:explicit[(canonical,str(code))][pos]=ml
    return seq,explicit,modes


def permutation(a,b,seed):
    values=np.array(a+b,dtype=float);na=len(a)
    observed=float(np.mean(b)-np.mean(a));space=math.comb(len(values),na)
    extreme=0
    if space<=100000:
        total=space
        for ids in itertools.combinations(range(len(values)),na):
            selected=values[list(ids)]
            stat=(values.sum()-selected.sum())/(len(values)-na)-selected.mean()
            extreme+=abs(stat)>=abs(observed)-1e-12
        return observed,extreme/total,total,'exact'
    rng=np.random.default_rng(seed);total=19999
    for _ in range(total):
        ids=rng.choice(len(values),na,replace=False);selected=values[ids]
        stat=(values.sum()-selected.sum())/(len(values)-na)-selected.mean()
        extreme+=abs(stat)>=abs(observed)-1e-12
    return observed,(extreme+1)/(total+1),total,'Monte_Carlo'


def modifications(sample,calls,seqs,out,p):
    selected=collections.defaultdict(list)
    for c in calls:
        if c['simple']:selected[c['read_id']].append(c)
    observations=collections.defaultdict(dict);modes_seen=collections.Counter();bad=0
    with pysam.AlignmentFile(sample['bam'],'rb') as f:
        for a in f.fetch(until_eof=True):
            if a.is_secondary or a.is_supplementary or a.query_name not in selected:continue
            try:query,explicit,modes=mod_values(a)
            except Exception:bad+=1;continue
            modes_seen.update(str(k)+v for k,v in modes.items())
            def value(base,code,q):
                if q in explicit.get((base,code),{}):return explicit[(base,code)][q]
                return 0 if modes.get((base,code))=='.' else None
            for c in selected[a.query_name]:
                s=seqs[c['contig']]
                for q,r,role in c['pairs']:
                    if q<0 or q>=len(query) or r<0 or r>=len(s):continue
                    base=query[q];strand='-' if c['reverse'] else '+'
                    refbase=s[r] if strand=='+' else rc(s[r])
                    if base!=refbase:continue
                    candidates=[]
                    if base=='A' and ('A','a') in modes:candidates.append(('6mA',value('A','a',q)))
                    cpg=(strand=='+' and s[r:r+2]=='CG') or (strand=='-' and r>0 and s[r-1:r+1]=='CG')
                    if base=='C' and cpg:
                        m=value('C','m',q);h=value('C','h',q)
                        candidates.extend([('5mC_CpG',m),('5hmC_CpG',h),('modified_C_CpG',None if m is None or h is None else min(255,m+h))])
                    for kind,v in candidates:
                        if v is not None:
                            key=(c['locus'],c['kind'],role,kind,strand,a.query_name)
                            observations[key][r]=int(v>=p['ml_cutoff'])
    rows=[];details=[];observation_groups=collections.defaultdict(list)
    for key,values in observations.items():
        if values:observation_groups[key[:5]].append((key[-1],values))
    for locus in sorted({c['locus'] for c in calls}):
        for baseline,state,role in [('L','P','head'),('R','P','tail'),('L','B','host_left'),('R','B','host_right')]:
            for kind in ('6mA','5mC_CpG','5hmC_CpG','modified_C_CpG'):
                for strand in ('+','-'):
                    groups=[]
                    for config in (baseline,state):
                        groups.append(observation_groups.get((locus,config,role,kind,strand),[]))
                    if not all(groups):continue
                    panel=set.intersection(*(set(vals) for group in groups for rid,vals in group))
                    if not panel:continue
                    fractions=[[sum(vals[pos] for pos in panel)/len(panel) for rid,vals in group] for group in groups]
                    delta=float(np.mean(fractions[1])-np.mean(fractions[0]))
                    row=dict(sample=sample['sample'],locus=locus,contrast=state+'-'+baseline,role=role,modification=kind,strand=strand,n_baseline=len(groups[0]),n_state=len(groups[1]),shared_positions=len(panel),baseline_percent=100*float(np.mean(fractions[0])),state_percent=100*float(np.mean(fractions[1])),delta_pp=100*delta,p_value=None,test='descriptive')
                    if min(len(groups[0]),len(groups[1]))>=3 and len(panel)>=3:
                        seed=int(hashlib.sha256(json.dumps(row,sort_keys=True).encode()).hexdigest()[:8],16)
                        delta,pv,nperm,method=permutation(fractions[0],fractions[1],seed)
                        row.update(p_value=pv,test=method,permutations=nperm)
                    rows.append(row)
                    for config,group,fracs in zip((baseline,state),groups,fractions):
                        for (rid,vals),frac in zip(group,fracs):details.append(dict(sample=sample['sample'],locus=locus,contrast=state+'-'+baseline,role=role,kind=kind,strand=strand,configuration=config,read_id=rid,shared_positions=len(panel),fraction=frac,panel=sorted(panel)))
    tsv(out/'modification_read_fractions.tsv',details)
    dump(out/'modification_contrasts.json',rows)
    tsv(out/'modification_contrasts.tsv',rows)
    return dict(contrasts=len(rows),tested=sum(r['p_value'] is not None for r in rows),malformed_mod_reads=bad,modes=dict(modes_seen)),rows


def sample_run(cfg,sample,out,threads):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    work=Path('work')/sample['sample'];work.mkdir(parents=True,exist_ok=True)
    p=cfg['parameters'];started=time.time()
    log('start',sample['sample'])
    seqs=fasta_read(sample['fasta']);provenance=input_provenance(sample)
    reads=work/'reads.fasta'
    source,events=export_source(sample,seqs,reads,p)
    dump(out/'input_provenance.json',dict(files=provenance,source=source,software=software(cfg)))
    log('exported',source['counts'])
    fresh=work/'sensitive_assembly.bam'
    command=map_reads(cfg,sample['fasta'],reads,fresh,out/'sensitive_mapping.log',threads)
    events.extend(discover_fresh(fresh,seqs,p))
    tsv(out/'discovery_events.tsv',events)
    groups=group_events(events,p)
    dump(out/'discovery_groups.json',groups)
    log('discovery',len(events),'events',len(groups),'groups')
    refs=build_references(groups,seqs,p)
    fasta_write(out/'junction_references.fasta',((r['id'],r['seq']) for r in refs))
    dump(out/'reference_origins.json',[dict(id=r['id'],length=len(r['seq']),sha256=hashlib.sha256(r['seq'].encode()).hexdigest(),aliases=r['aliases']) for r in refs])
    recruited=work/'recruited.fasta'
    nrecruit=recruit_reads(fresh,groups,recruited,p['external_flank'])
    calls=[];reject={};ref_command=[];bg_command=[]
    if refs and nrecruit:
        mapped=work/'synthetic.bam'
        ref_command=map_reads(cfg,out/'junction_references.fasta',recruited,mapped,out/'synthetic_mapping.log',threads,True)
        calls,reject=validate_synthetic(cfg,refs,mapped,p)
        log('synthetic passing',len(calls),'checking complete assembly')
        bg_command=background_check(cfg,sample,calls,seqs,out,threads)
    dump(out/'calls.json',calls)
    tsv(out/'junction_calls.tsv',({k:v for k,v in c.items() if k not in ('pairs','query_slice')} for c in calls))
    loci=loci_summary(groups,calls,sample['sample'])
    dump(out/'loci.json',loci);tsv(out/'loci.tsv',loci)
    mod_summary,contrasts=modifications(sample,calls,seqs,out,p)
    summary=dict(sample=sample['sample'],status='complete',elapsed_seconds=round(time.time()-started,1),input=source['counts'],assembly_contigs=len(seqs),assembly_bases=sum(map(len,seqs.values())),discovery_events=len(events),discovery_groups=len(groups),short_discovery_groups=sum(g['b']-g['a']<2000 for g in groups),synthetic_references=len(refs),recruited_reads=nrecruit,passing_calls=len(calls),rejection_counts=reject,confirmed_simple_P=sum(r['confirmed_simple_P'] for r in loci),confirmed_simple_B=sum(r['confirmed_simple_B'] for r in loci),confirmed_short_P=sum(r['confirmed_simple_P'] and r['short'] for r in loci),confirmed_short_B=sum(r['confirmed_simple_B'] and r['short'] for r in loci),modifications=mod_summary,commands=[command,ref_command,bg_command],loci=loci)
    dump(out/'summary.json',summary)
    log('complete',sample['sample'],summary['confirmed_short_P'],'short P',summary['confirmed_short_B'],'short B')
    return summary


def mutate(seq,rate,rng):
    # Substitutions and one-base indels; deterministic independent technical errors.
    out=[]
    for base in seq:
        u=rng.random()
        if u<rate/3:continue
        if u<2*rate/3:out.append(rng.choice('ACGT'))
        out.append(rng.choice([x for x in 'ACGT' if x!=base]) if 2*rate/3<=u<rate else base)
    return ''.join(out)


def controls(cfg,out,threads):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);p=cfg['parameters'];rng=random.Random(181605)
    seqs={};truth=[];reads=[];groups=[]
    for length in (200,430,800,1500,3000):
        s=''.join(rng.choices('ACGT',k=10000+length));ctg='control_'+str(length)
        a=5000;b=a+length;seqs[ctg]=s
        g=dict(id='T'+str(length),contig=ctg,a=a,b=b,hypotheses=[(a,b)],seed_reads=[],support=5,routes={})
        groups.append(g)
        h=min(2000,length//2)
        examples={'P':s[b-h:b]+s[a:a+h],'L':s[a-1000:a]+s[a:a+h],'R':s[b-h:b]+s[b:b+1000],'B':s[a-1000:a]+s[b:b+1000],'ordinary':s[a-1000:b+1000],'tandem':s[a-1000:b]+s[a:b]+s[b:b+1000]}
        for rate in (0,.01,.03):
            for kind,seq in examples.items():
                for rep in range(5):
                    rid=f'{ctg}_{kind}_{rate}_{rep}';reads.append((rid,mutate(seq,rate,rng)));truth.append(dict(read_id=rid,length=length,kind=kind,error_rate=rate))
        for rep in range(5):
            reads.append((ctg+'_chimera_'+str(rep),s[100:500]+''.join(rng.choices('ACGT',k=400))))
            truth.append(dict(read_id=ctg+'_chimera_'+str(rep),length=length,kind='unrelated_chimera',error_rate=0))
    # Exact duplicate source locus: end-to-start family evidence must retain both aliases.
    duplicate=dict(groups[2],id='T800_copy',contig='control_800_copy')
    seqs[duplicate['contig']]=seqs['control_800'];groups.append(duplicate)
    refs=build_references(groups,seqs,p)
    duplicate_refs=[r for r in refs if any(x['locus']=='T800' and x['kind']=='P' for x in r['aliases'])]
    assert any({x['locus'] for x in r['aliases']} >= {'T800','T800_copy'} for r in duplicate_refs)
    # Both repeat-copy rules and the no-repeat model must coexist.
    repeat='ACGTACGTTG';s=''.join(rng.choices('ACGT',k=1400));s=s[:500]+repeat+s[510:890]+repeat+s[900:]
    repeat_refs=build_references([dict(id='repeat',contig='repeat',a=500,b=900,hypotheses=[(500,900)])],{'repeat':s},p)
    models={a['model'] for r in repeat_refs for a in r['aliases']}
    assert {'no_repeat','remove_right_repeat','remove_left_repeat'} <= models
    assert all(r['aliases'][0]['join']<=len(r['seq'])-75 for r in refs)
    fasta_write(out/'assembly.fasta',seqs.items());fasta_write(out/'reads.fasta',reads);fasta_write(out/'references.fasta',((r['id'],r['seq']) for r in refs))
    map_reads(cfg,out/'references.fasta',out/'reads.fasta',out/'synthetic.bam',out/'mapping.log',threads,True)
    calls,reject=validate_synthetic(cfg,refs,out/'synthetic.bam',p)
    background_check(cfg,{'fasta':str(out/'assembly.fasta')},calls,seqs,out,threads)
    by_read=collections.defaultdict(set)
    for c in calls:
        if not c['genomic_contiguous_alternative']:by_read[c['read_id']].add(c['kind'])
    rows=[]
    for length in (200,430,800,1500,3000):
        for rate in (0,.01,.03):
            for kind in ('P','L','R','B','ordinary','tandem'):
                selected=[t for t in truth if t['length']==length and t['error_rate']==rate and t['kind']==kind]
                detected=sum((kind in by_read[t['read_id']]) if kind in ('P','L','R','B') else bool(by_read[t['read_id']] & {'P','B'}) for t in selected)
                complex_count=len({c['read_id'] for c in calls if c['read_id'] in {t['read_id'] for t in selected} and c['complex_read']})
                rows.append(dict(length=length,error_rate=rate,kind=kind,total=len(selected),detected=detected,complex_reads=complex_count))
    # Test the discovery parser against real mapper output, separately from reference validation.
    map_reads(cfg,out/'assembly.fasta',out/'reads.fasta',out/'assembly.bam',out/'discovery_mapping.log',threads)
    events=discover_fresh(out/'assembly.bam',seqs,p)
    discovered={e['read_id'] for e in events if e['event']=='P'}
    for row in rows:
        selected=[t for t in truth if all(t[k]==row[k] for k in ('length','error_rate','kind'))]
        row['discovery_P_reads']=sum(t['read_id'] in discovered for t in selected)
    # Discovery on reverse-strand segments is explicitly checked.
    toy=[dict(contig='x',r0=1000,r1=1200,q0=0,q1=200,reverse=True,identity=1,paired=200,mapq=60,route='test'),dict(contig='x',r0=1200,r1=1400,q0=200,q1=400,reverse=True,identity=1,paired=200,mapq=60,route='test')]
    assert split_events('reverse',toy,p,'unit')[0]['length']==400
    perfect=[r for r in rows if r['error_rate']==0 and r['kind'] in ('P','B')]
    failed=[r for r in perfect if r['detected']<r['total']]
    false_positive=[c for c in calls if ('_ordinary_' in c['read_id'] or '_chimera_' in c['read_id']) and c['kind'] in ('P','B') and not c['genomic_contiguous_alternative']]
    summary=dict(status='passed' if not failed and not false_positive else 'failed',software=software(cfg),parameters=p,rows=rows,failed_perfect_controls=failed,false_positive_junction_calls=len(false_positive),duplicate_alias_preserved=True,repeat_models=sorted(models),rejection_counts=reject,discovery_events=len(events))
    tsv(out/'control_results.tsv',rows);dump(out/'calls.json',calls);dump(out/'summary.json',summary)
    if failed or false_positive:raise RuntimeError('Synthetic control gate failed; inspect results/controls/summary.json')
    log('controls passed')


def aggregate(cfg,out):
    out=Path(out);summaries=[];loci=[];contrasts=[]
    for sample in cfg['samples']:
        folder=out/sample['sample'];summary=json.loads((folder/'summary.json').read_text());summaries.append(summary)
        loci.extend(json.loads((folder/'loci.json').read_text()));contrasts.extend(json.loads((folder/'modification_contrasts.json').read_text()))
    tested=[(i,r['p_value']) for i,r in enumerate(contrasts) if r['p_value'] is not None]
    tested.sort(key=lambda x:x[1]);minimum=1.0
    for rank in range(len(tested),0,-1):
        i,pv=tested[rank-1];minimum=min(minimum,pv*len(tested)/rank);contrasts[i]['BH_q']=minimum
    tsv(out/'loci.tsv',loci);tsv(out/'modification_contrasts.tsv',contrasts)
    dump(out/'modification_contrasts.json',contrasts)
    dump(out/'summary.json',dict(version=VERSION,status='complete',samples=summaries,controls=json.loads((out/'controls/summary.json').read_text()),tested_contrasts=len(tested),BH_q_lt_005=sum(r.get('BH_q',1)<.05 for r in contrasts)))


def copy_controls(cfg,out,threads):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    old=Path(cfg['reuse_results'])/'controls'
    summary=json.loads((old/'summary.json').read_text());assert summary['status']=='passed'
    seqs=fasta_read(old/'assembly.fasta');groups=[]
    for length in (430,800):groups.append(dict(id='copy_test_'+str(length),contig='control_'+str(length),a=5000,b=5000+length,hypotheses=[(5000,5000+length)]))
    refs=build_references(groups,seqs,cfg['parameters']);calls=[]
    for ref in refs:
        al=ref['aliases'][0]
        if al['kind']!='P':continue
        w=min(250,al['h']);c=al['join'];qs=ref['seq'][c-w:c+w]
        calls.append(dict(read_id=al['locus'],ref_id=ref['id'],locus=al['locus'],contig=al['contig'],kind='P',model=al['model'],a=al['a'],b=al['b'],h=al['h'],window=w,query_slice=qs,edit_distance=0,required_margin=max(3,math.ceil(.02*len(qs))),locus_ambiguous=False,simple=True))
    result=audit_source_copies(cfg,dict(fasta=str(old/'assembly.fasta')),calls,refs,seqs,out,threads,cfg['parameters'])
    unique=[c for c in calls if c['locus']=='copy_test_430'];duplicate=[c for c in calls if c['locus']=='copy_test_800']
    assert unique and duplicate
    assert all(not c['source_copy_ambiguous'] for c in unique),'Unique source falsely declared ambiguous'
    assert all(c['source_copy_ambiguous'] for c in duplicate),'Unseeded identical copy was missed'
    summary['copy_control']=dict(status='passed',unique_source_retained=True,unseeded_duplicate_detected=True,audit=result)
    if cfg.get('normalize_occupied_labels'):
        aliases=[dict(locus='same_place',kind='L',contig='c',arms=[dict(start=500,end=1000),dict(start=1000,end=2000)]),dict(locus='different_copy',kind='L',contig='other',arms=[dict(start=500,end=1000),dict(start=1000,end=2000)])]
        test=dict(kind='L',locus='original',contig='c',pairs=[(q,900+q,'host_left') for q in range(200)],forward_q0=0,forward_q1=200,competing_loci=['same_place'],simple=False,complex_read=False,genomic_contiguous_alternative=False)
        test_other=dict(test,competing_loci=['different_copy'])
        restored=normalize_occupied_labels([test,test_other],[dict(aliases=aliases)])
        assert restored==1 and test['simple'] and not test_other['simple']
        summary['occupied_placement_control']=dict(status='passed',same_genomic_placement_retained=True,different_copy_remains_ambiguous=True)
    dump(out/'summary.json',summary)
    log('unseeded source-copy controls passed')


def fast_controls(cfg,out,threads):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if cfg.get('remap_junctions') or cfg.get('verify_seed_controls'):
        controls(cfg,out/'seed_controls',threads)
    old=Path(cfg['reuse_results'])/'controls';seqs=fasta_read(old/'assembly.fasta')
    groups=[dict(id='T'+str(n),contig='control_'+str(n),a=5000,b=5000+n,hypotheses=[(5000,5000+n)]) for n in (200,430,800,1500,3000)]
    groups.append(dict(groups[2],id='T800_copy',contig='control_800_copy'))
    refs=build_references(groups,seqs,cfg['parameters'])
    if cfg.get('verify_interval_bounds'):
        lookup={r['id']:r for r in refs};checked=0
        with pysam.AlignmentFile(old/'synthetic.bam','rb') as f:
            for a in f:
                if a.is_unmapped or not a.query_sequence:continue
                ref=lookup[a.reference_name]
                for al in ref['aliases']:
                    w=min(cfg['parameters']['base_window'],al['h'])
                    assert alignment_metric(a,ref['seq'],al,w,cfg['parameters'])==alignment_metric_reference(a,ref['seq'],al,w,cfg['parameters'])
                    checked+=1
        log('interval-bound equivalence checks passed',checked)
    exhaustive,reject1=validate_synthetic(dict(cfg,indexed_competition=False,ranked_candidate_search=False),refs,old/'synthetic.bam',cfg['parameters'])
    indexed,reject2=validate_synthetic(dict(cfg,indexed_competition=True),refs,old/'synthetic.bam',cfg['parameters'])
    assert exhaustive==indexed,'Indexed/ranked competitor search differs from exhaustive search'
    if not cfg.get('ranked_candidate_search'):assert reject1==reject2
    parallel_control=None
    if cfg.get('parallel_validation_workers',1)>1:
        parallel_cfg=dict(cfg,validation_workers=2,validation_chunks=7,validation_checkpoint_dir=str(out/'parallel_control'))
        parallel,parallel_rejections=validate_synthetic_parallel(parallel_cfg,refs,old/'synthetic.bam',cfg['parameters'])
        assert parallel==indexed and parallel_rejections==reject2,'Parallel validation differs from serial validation'
        resumed,resumed_rejections=validate_synthetic_parallel(parallel_cfg,refs,old/'synthetic.bam',cfg['parameters'])
        assert resumed==indexed and resumed_rejections==reject2,'Resumed checkpoints changed validation'
        state=json.loads((out/'parallel_control/summary.json').read_text())
        assert state['reused_chunks']==state['chunks'] and state['chunks']>1
        parallel_control=dict(status='passed',selected_calls=len(parallel),serial_parallel_calls_equal=True,rejection_counts_equal=True,resume_equal=True,all_checkpoints_reused=True,chunks=state['chunks'])
        # Bounded actual-data regression: use the unchanged complete NOSH BAM
        # and compare both implementations before any large sample starts.
        nosh=next(s for s in cfg['samples'] if s['sample']=='NOSH_AS_18')
        origin=Path(cfg['reuse_results'])/nosh['sample']
        sequences=fasta_read(origin/'junction_references.fasta')
        origins=json.loads((origin/'reference_origins.json').read_text())
        actual_refs=[dict(id=r['id'],seq=sequences[r['id']],aliases=r['aliases']) for r in origins]
        actual_bam=Path(cfg['source_mapped_results'])/nosh['sample']/cfg['source_mapped_filename']
        actual_cfg=dict(cfg,keep_A_modification_pairs=True)
        t0=time.monotonic();serial,serial_rejections=validate_synthetic_ranked(actual_cfg,actual_refs,actual_bam,cfg['parameters']);serial_seconds=time.monotonic()-t0
        parallel_actual=dict(actual_cfg,validation_workers=2,validation_chunks=5,validation_checkpoint_dir=str(out/'parallel_nosh_regression'))
        t0=time.monotonic();parallel,parallel_rejections=validate_synthetic_parallel(parallel_actual,actual_refs,actual_bam,cfg['parameters']);parallel_seconds=time.monotonic()-t0
        assert serial==parallel and serial_rejections==parallel_rejections,'Actual NOSH serial/parallel mismatch'
        previous=json.loads((Path(cfg['reuse_completed_results'])/nosh['sample']/'calls.json').read_text())
        keys=list(serial[0]) if serial else []
        normalize=lambda x:json.loads(json.dumps(x))
        assert normalize(serial)==[{k:c[k] for k in keys} for c in previous],'Previously completed NOSH calls changed'
        parallel_control['actual_NOSH']=dict(status='passed',calls=len(serial),serial_parallel_calls_equal=True,previous_complete_calls_equal=True,serial_seconds=round(serial_seconds,3),parallel_seconds=round(parallel_seconds,3))
    # Sparse storage retains every A/CpG observation used by the modification
    # code; compare it to dense coordinate storage on actual control alignments.
    dense,_=validate_synthetic(dict(cfg,indexed_competition=True,sparse_modification_pairs=False),refs,old/'synthetic.bam',cfg['parameters'])
    def projection(c):return {k:v for k,v in c.items() if k!='pairs'}
    assert [projection(c) for c in dense]==[projection(c) for c in indexed]
    for full,sparse in zip(dense,indexed):assert set(sparse['pairs'])<=set(full['pairs'])
    copy_controls(cfg,out,threads)
    summary=json.loads((out/'summary.json').read_text())
    summary['indexed_competition_control']=dict(status='passed',exhaustive_call_count=len(exhaustive),indexed_call_count=len(indexed),identical_selected_calls=True,rejection_counters_equal=reject1==reject2,sparse_coordinates_preserve_classification=True,seed_length=11,ranked=bool(cfg.get('ranked_candidate_search')))
    if parallel_control is not None:
        summary['parallel_validation_control']=parallel_control
        summary['validation_software']=software(cfg)
    if cfg.get('remap_junctions') or cfg.get('verify_seed_controls'):
        summary['junction_seed_control']=json.loads((out/'seed_controls/summary.json').read_text())
        summary['junction_seed_length']=cfg['junction_seed_length']
    dump(out/'summary.json',summary)
    log('indexed/exhaustive equivalence controls passed',len(indexed),'calls')


def complete_bam(path):
    try:
        with open(path,'rb') as f:
            f.seek(-28,os.SEEK_END)
            return f.read()==bytes.fromhex('1f8b08040000000000ff0600424302001b0003000000000000000000')
    except (OSError,ValueError):return False


def fast_sample(cfg,sample,out,threads):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);started=time.time()
    old=Path(cfg['reuse_results'])/sample['sample']
    native_results=Path(cfg['reuse_results']).parent
    mapped=native_results/'work'/sample['sample']/'synthetic.bam'
    if cfg.get('source_mapped_results'):mapped=Path(cfg['source_mapped_results'])/sample['sample']/cfg.get('source_mapped_filename','synthetic_seed11.bam')
    recruited=native_results/'work'/sample['sample']/'recruited.fasta'
    deadline=time.monotonic()+cfg.get('alignment_wait_seconds',21600)
    log('waiting for completed existing synthetic BAM',str(mapped))
    required=[old/'discovery_groups.json',old/'reference_origins.json',old/'junction_references.fasta',old/'input_provenance.json']
    while not (all(p.is_file() for p in required) and (recruited.is_file() if cfg.get('remap_junctions') else complete_bam(mapped))):
        if time.monotonic()>deadline:raise TimeoutError('Existing approved mapping stage did not complete')
        time.sleep(5)
    provenance=json.loads((old/'input_provenance.json').read_text());source=provenance['source']
    new_map_command=[]
    if cfg.get('remap_junctions'):
        mapped=out/'synthetic_seed11.bam'
        new_map_command=map_reads(cfg,old/'junction_references.fasta',recruited,mapped,out/'synthetic_mapping.log',threads,True)
    log('reusing completed mapped BAM',sample['sample'])
    assert complete_bam(mapped)
    dump(out/'alignment_reuse_status.json',dict(source_bam=str(mapped),BGZF_EOF_present=True,source_results=str(old),source_provenance=provenance['files'],recomputed_junction_mapping=bool(cfg.get('remap_junctions')),source_recruited_fasta=str(recruited)))
    seqs=fasta_read(sample['fasta']);groups=json.loads((old/'discovery_groups.json').read_text());origins=json.loads((old/'reference_origins.json').read_text());refseq=fasta_read(old/'junction_references.fasta')
    refs=[dict(id=r['id'],seq=refseq[r['id']],aliases=r['aliases']) for r in origins]
    local=dict(cfg,keep_A_modification_pairs=any(x.startswith('A+a') for x in source['modification_groups']),validation_workers=min(threads,int(cfg.get('parallel_validation_workers',1))),validation_checkpoint_dir=str(out/'validation_chunks'))
    calls,reject=validate_synthetic(local,refs,mapped,cfg['parameters'])
    log('indexed synthetic validation complete',len(calls),'calls')
    regression={}
    comparison=Path(cfg.get('comparison_results',str(Path(cfg['reuse_results']))))/sample['sample']/'calls.json'
    if comparison.exists():
        previous=json.loads(comparison.read_text())
        keys=['read_id','ref_id','locus','contig','kind','model','a','b','h','repeat_length','identity','window','forward_q0','forward_q1','query_slice','edit_distance','required_margin','competing_distance','competing_loci','reverse','tiers']
        project=lambda rows:[{k:c[k] for k in keys} for c in rows]
        equal=project(calls)==project(previous)
        if not cfg.get('remap_junctions'):assert equal,'Indexed calls differ from the original exhaustive calls'
        regression['original_exhaustive_calls_equal']=equal
    bg=background_check(cfg,sample,calls,seqs,out,threads)
    audit=audit_source_copies(cfg,sample,calls,origins,seqs,out,threads,cfg['parameters'])
    audit['occupied_labels_normalized']=normalize_occupied_labels(calls,origins)
    loci=loci_summary(groups,calls,sample['sample']);flye={}
    for line in Path(sample['assembly_info']).read_text().splitlines():
        if line.startswith('#'):continue
        parts=line.split()
        if len(parts)>=4:flye[parts[0]]=parts
    for row in loci:
        n=len(seqs[row['contig']]);row['contig_length']=n
        row['near_contig_ends']=row['a']<=50 and n-row['b']<=50
        row['flye_circular']=flye.get(row['contig'],['','','','unknown'])[3]
    mod_summary,contrasts=modifications(sample,calls,seqs,out,cfg['parameters'])
    other=Path(cfg.get('regression_results','/nonexistent'))/sample['sample']/'modification_contrasts.json'
    if other.exists():
        equal=contrasts==json.loads(other.read_text())
        if not cfg.get('remap_junctions'):assert equal,'Sparse modification coordinates changed a completed contrast'
        regression['dense_modification_contrasts_equal']=equal
    with open(old/'discovery_events.tsv') as f:nevents=max(0,sum(1 for _ in f)-1)
    with open(recruited,'rb') as f:nrecruit=sum(line.startswith(b'>') for line in f)
    dump(out/'calls.json',calls);dump(out/'loci.json',loci);tsv(out/'loci.tsv',loci)
    tsv(out/'junction_calls.tsv',({k:v for k,v in c.items() if k not in ('pairs','query_slice')} for c in calls))
    dump(out/'reference_origins.json',origins);dump(out/'input_provenance.json',dict(provenance,validation_software=software(cfg)))
    summary=dict(sample=sample['sample'],version=VERSION,status='complete',elapsed_seconds=round(time.time()-started,1),input=source['counts'],assembly_contigs=len(seqs),assembly_bases=sum(map(len,seqs.values())),discovery_events=nevents,discovery_groups=len(groups),short_discovery_groups=sum(g['b']-g['a']<2000 for g in groups),synthetic_references=len(refs),recruited_reads=nrecruit,passing_calls=len(calls),rejection_counts=reject,copy_audit=audit,modifications=mod_summary,source_results=str(old),reused_alignment=str(mapped),regression=regression,commands=[new_map_command,bg],loci=loci)
    for kind in ('P','B'):
        summary['confirmed_simple_'+kind]=sum(r['confirmed_simple_'+kind] for r in loci)
        summary['confirmed_short_'+kind]=sum(r['confirmed_simple_'+kind] and r['short'] for r in loci)
    dump(out/'summary.json',summary)
    log('fast validation complete',sample['sample'],summary['confirmed_short_P'],'short P',summary['confirmed_short_B'],'short B',regression)


def reuse_complete_sample(cfg,sample,out):
    source=Path(cfg['reuse_completed_results'])/sample['sample'];out=Path(out)
    summary=json.loads((source/'summary.json').read_text());assert summary['status']=='complete'
    out.mkdir(parents=True,exist_ok=True)
    copy_names={'loci.json','modification_contrasts.json'}
    for entry in source.iterdir():
        if not entry.is_file() or entry.name=='summary.json':continue
        target=out/entry.name
        if entry.name in copy_names:shutil.copy2(entry,target)
        else:target.symlink_to(entry.resolve())
    dump(out/'reuse_provenance.json',dict(status='complete',source=str(source),source_pipeline_version=summary.get('version'),reused_without_recalculation=True))
    shutil.copy2(source/'summary.json',out/'summary.json')
    log('reused complete sample unchanged',sample['sample'],str(source))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['sample','controls','aggregate','audit','copy_controls','fast_controls','fast_sample','reuse_complete']);ap.add_argument('--config',required=True);ap.add_argument('--sample');ap.add_argument('--output',required=True);ap.add_argument('--threads',type=int,default=8);args=ap.parse_args()
    cfg=json.loads(Path(args.config).read_text())
    if args.mode=='controls':controls(cfg,args.output,args.threads)
    elif args.mode=='fast_controls':fast_controls(cfg,args.output,args.threads)
    elif args.mode=='fast_sample':fast_sample(cfg,next(x for x in cfg['samples'] if x['sample']==args.sample),args.output,args.threads)
    elif args.mode=='reuse_complete':reuse_complete_sample(cfg,next(x for x in cfg['samples'] if x['sample']==args.sample),args.output)
    elif args.mode=='copy_controls':copy_controls(cfg,args.output,args.threads)
    elif args.mode=='aggregate':aggregate(cfg,args.output)
    elif args.mode=='audit':audit_existing(cfg,next(x for x in cfg['samples'] if x['sample']==args.sample),args.output,args.threads)
    else:sample_run(cfg,next(x for x in cfg['samples'] if x['sample']==args.sample),args.output,args.threads)


if __name__=='__main__':
    main()
