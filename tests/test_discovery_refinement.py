import random

import pytest

from metaphaser._vendor import junctions_v17 as core
from metaphaser.config import DEFAULTS
from metaphaser.discovery import group_and_refine
from metaphaser.reporting import classify_loci


def event(read, a=1000, b=5000, gap=0, query=None):
    row = dict(read_id=read, contig='ctg', a=a, b=b, length=b-a, gap=gap,
               identity=1., route='test', event='P')
    if query is not None:
        row['junction_query'] = query
    return row


def sequence():
    rng = random.Random(173)
    return ''.join(rng.choices('ACGT', k=7000))


def test_two_discovery_reads_do_not_become_three_validated_reads():
    seq = sequence()
    groups, stats = group_and_refine([event('a'), event('b'), event('a')], {'ctg': seq}, DEFAULTS)
    assert len(groups) == 1 and groups[0]['support'] == 2
    locus = dict(sample='x', locus=groups[0]['id'], contig='ctg', a=1000, b=5000,
                 L_simple=8, R_simple=6, P_simple=2, B_simple=0)
    rows, regions = classify_loci([locus], {'ctg': seq}, DEFAULTS['min_support'])
    assert not rows[0]['host_linked_P'] and not regions
    locus['P_simple'] = 3
    rows, regions = classify_loci([locus], {'ctg': seq}, DEFAULTS['min_support'])
    assert rows[0]['host_linked_P'] and len(regions) == 1


def test_refinement_recovers_a_boundary_inside_an_unaligned_query_gap():
    seq = sequence()
    true_a, b, gap = 1000, 5000, 20
    query = seq[b-250:b] + seq[true_a:true_a+250+gap]
    events = [event('a', true_a+gap, b, gap, query), event('b', true_a+gap, b)]
    groups, stats = group_and_refine(events, {'ctg': seq}, DEFAULTS)
    assert [true_a+gap, b] in groups[0]['observed_hypotheses']
    assert [true_a, b] in groups[0]['hypotheses']
    assert groups[0]['support'] == 2
    assert all(v['edits'] < v['original_edits'] for h in groups[0]['refined_hypotheses'] for v in h['evidence'])
    assert len(groups[0]['refined_hypotheses']) <= DEFAULTS['refined_breakpoint_hypotheses']


def test_a_positive_gap_does_not_require_a_reference_boundary_shift():
    seq = sequence()
    query = seq[4750:5000] + 'N'*12 + seq[1000:1250]
    groups, _ = group_and_refine([event('a', gap=12, query=query), event('b')], {'ctg': seq}, DEFAULTS)
    assert groups[0]['hypotheses'] == [[1000, 5000]]


def test_duplicate_events_cannot_satisfy_discovery_threshold():
    groups, _ = group_and_refine([event('a'), event('a')], {'ctg': sequence()}, DEFAULTS)
    assert not groups


@pytest.mark.parametrize('reverse', [False, True])
def test_query_window_keeps_original_read_coordinates(reverse):
    coords = [(1000, 1200), (4800, 5000)] if reverse else [(4800, 5000), (1000, 1200)]
    pieces = [dict(contig='ctg', r0=a, r1=b, q0=q0, q1=q1, reverse=reverse,
                   identity=1., paired=200, mapq=60, route='test')
              for (a, b), (q0, q1) in zip(coords, [(50, 250), (270, 470)])]
    row, = core.split_events('x', pieces, DEFAULTS, 'test', query_windows=True)
    assert (row['a'], row['b'], row['gap']) == (1000, 5000, 20)
    assert (row['query_window_start'], row['query_window_end']) == (50, 470)
    assert row['query_reverse'] == reverse
