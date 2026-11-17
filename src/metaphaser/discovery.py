"""Candidate discovery and read-sequence refinement before strict validation."""
from __future__ import annotations

import collections
import itertools

import edlib
import pysam

from ._vendor import junctions_v17 as core


def discover_fresh(path, seqs, parameters):
    events = []
    with pysam.AlignmentFile(path, 'rb') as handle:
        for name, records in itertools.groupby(handle, key=lambda a: a.query_name):
            rows = list(records)
            segments = list(itertools.chain.from_iterable(core.segments_of(a) for a in rows))
            proposed = core.split_events(name, segments, parameters, 'sensitive_remap', query_windows=True)
            complete = next((a for a in rows if a.query_sequence and
                             not any(op == 5 for op, _ in a.cigartuples or [])), None)
            sequence = complete.get_forward_sequence() if complete is not None else None
            for event in proposed:
                if sequence is not None and 'query_window_start' in event:
                    start, end = event['query_window_start'], event['query_window_end']
                    if 0 <= start < end <= len(sequence):
                        window = sequence[start:end]
                        event['junction_query'] = core.rc(window) if event['query_reverse'] else window
            events.extend(proposed)
            for record in rows:
                if not record.is_secondary:
                    events.extend(core.cigar_events(record, seqs, parameters))
    return events


def junction_score(query, contig, a, b, seqs, parameters):
    group = dict(id='refinement', contig=contig, a=a, b=b, hypotheses=[(a, b)])
    best = None
    for reference in core.build_references([group], seqs, parameters):
        aliases = [al for al in reference['aliases'] if al['kind'] == 'P']
        if not aliases:
            continue
        match = edlib.align(query, reference['seq'], mode='HW', task='locations')
        crosses = any(start <= al['join']-parameters['min_anchor'] and
                      end+1 >= al['join']+parameters['min_anchor']
                      for al in aliases for start, end in match['locations'])
        if crosses and (best is None or match['editDistance'] < best):
            best = match['editDistance']
    return best


def group_and_refine(events, seqs, parameters):
    # Discovery is a screen. Only final validated read counts retain an element.
    seed_min = min(parameters['min_support'], parameters['min_discovery_reads'])
    groups = core.group_events(events, dict(parameters, min_support=seed_min))
    summary = dict(discovery_min_reads=seed_min, final_min_reads=parameters['min_support'],
                   groups=len(groups), groups_below_final_seed_count=0,
                   refined_groups=0, added_hypotheses=0, tested_gap_hypotheses=0)
    for group in groups:
        if group['support'] < parameters['min_support']:
            summary['groups_below_final_seed_count'] += 1
        original = [list(pair) for pair in group['hypotheses']]
        original_set = {tuple(pair) for pair in original}
        candidates = collections.defaultdict(dict)
        details = []
        cache = {}
        seen = set()
        for event in group['events']:
            gap = event['gap']
            query = event.get('junction_query')
            if event.get('event') != 'P' or not query or not 0 < gap <= parameters['max_query_gap']:
                continue
            signature = (event['read_id'], event['a'], event['b'], query)
            if signature in seen:
                continue
            seen.add(signature)
            a, b = event['a'], event['b']
            def score(pair):
                key = (query, *pair)
                if key not in cache:
                    cache[key] = junction_score(query, group['contig'], *pair, seqs, parameters)
                return cache[key]
            baseline = [score(pair) for pair in original_set]
            baseline = min((s for s in baseline if s is not None), default=None)
            if baseline is None:
                continue
            proposed = {(a-step, b) for step in range(1, gap+1)}
            proposed.update((a, b+step) for step in range(1, gap+1))
            for pair in sorted(proposed-original_set):
                aa, bb = pair
                if aa < 0 or bb > len(seqs[group['contig']]) or not parameters['min_length'] <= bb-aa <= parameters['max_length']:
                    continue
                summary['tested_gap_hypotheses'] += 1
                distance = score(pair)
                if distance is not None and distance < baseline:
                    previous = candidates[pair].get(event['read_id'])
                    evidence = dict(read_id=event['read_id'], edits=distance,
                                    original_edits=baseline, improvement=baseline-distance,
                                    query_length=len(query), observed_a=a, observed_b=b, read_gap=gap)
                    if previous is None or evidence['improvement'] > previous['improvement']:
                        candidates[pair][event['read_id']] = evidence
        ranked = sorted(candidates, key=lambda pair: (-len(candidates[pair]),
                        -sum(v['improvement'] for v in candidates[pair].values()),
                        min(v['edits'] for v in candidates[pair].values()), pair))
        added = ranked[:parameters['refined_breakpoint_hypotheses']]
        for pair in added:
            details.append(dict(a=pair[0], b=pair[1], evidence=list(candidates[pair].values())))
        group['observed_hypotheses'] = original
        group['refined_hypotheses'] = details
        group['hypotheses'] = original + [list(pair) for pair in added]
        summary['refined_groups'] += bool(added)
        summary['added_hypotheses'] += len(added)
        for event in group['events']:
            event.pop('junction_query', None)
    return groups, summary
