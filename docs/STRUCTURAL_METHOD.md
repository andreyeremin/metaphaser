# Algorithm and implementation boundaries

The structural implementation is in `src/metaphaser/_vendor/junctions_v17.py`,
with candidate refinement in `src/metaphaser/discovery.py`. The pipeline calls
their discovery, reference-building, validation, background and source-copy
routines. Version 0.2.0rc2 adds the changes described in
[DISCOVERY_FIX_20261117.md](DISCOVERY_FIX_20261117.md).

## Discovery

Complete primary sequences are exported after QC/duplicate-flag filtering and
identifier deduplication. Missing sequences and hard-clipped primary records
are recorded as input-filtered. Original primary/SA evidence is combined with
fresh full-read mappings unless original alignments were explicitly disabled.

Both discovery and synthetic mapping use minimap2:

```
-a -x map-ont -k 11 -w 5 -m 20 -s 40 -n 2 -N 50 -p 0.5
--secondary=yes --secondary-seq -Y --eqx --mask-len 75
--no-long-join -r 100,200
```

Original read coordinates are used for ordering, including reverse mappings.
All eligible same-contig, same-orientation piece pairs are examined. Each must
contribute 75 read bases outside the mutual overlap. The gap or overlap is at
most 250 bases. Backward jumps of 200 to 1,000,000 bases propose P intervals.
CIGAR deletions in that range propose B. Insertions of 200 to 2000 bases propose
P when they match adjacent contig sequence within a 15% edit bound. Each side
must contain at least 75 paired bases in total.

Discovery segments require at least 75 paired bases and 85% identity.

Grouping uses both distance from the initial anchor and a full boundary-range
bound of twice the tolerance. It cannot chain neighbouring events indefinitely.
Observed boundary pairs are ranked by distinct read count, then absolute read
gap, then coordinates. The first five are retained.

The boundary tolerance is 25 bp for intervals shorter than 2000 bp and 200 bp
for longer intervals. The preliminary screen requires two distinct discovery
read identifiers. Positive gaps are checked for boundary extensions, retaining
up to five improved pairs in addition to the original hypotheses. Final
host-linked intervals still require three simple P reads and three L or R reads.

## References and acceptance

L and R preserve the original boundaries. P joins tail to head and B joins host
flanks. The original four references are always tested when both arms have at
least 75 bases. Exact A/C/G/T matches of 4 to 100 bases also supply local
boundary alternatives. A match at the interval start and immediately after its
end shifts both P cuts right. Matching copies at the two internal ends supply
P models omitting either copy and a B model retaining the left copy. Shift
length is limited to half the interval and available right sequence. Removal
length is limited to a third of the interval. Identical references keep all
locus aliases.

External flanks extend up to 2000 bp and stop at contig ends. Internal arms
are capped at 2000 bp and half the interval length after boundary adjustment.
References with fewer than 75 bases on either side are omitted. The default
validation window is up to 250 bases on each side of the join, limited by the
available arms. At least 90% of each side and 90% of the central 50 bases must
be paired. Identity must be at least 95%, with no central indel longer than
10 bases. The required edit-distance advantage over alternative junction
configurations is the larger of three edits and 2% of the query-window length,
rounded up. Modification extraction extends up to 1000 bp from the junction,
limited by the internal arm length.

Comparisons use the same query slice in both orientations. Compatible occupied
labels describing the same genomic placement do not create source ambiguity.
P and B must beat uninterrupted genomic explanations. Source-copy checks for P
also map its head and tail separately, including loci with no discovery seeds.
Arm placements need 90% coverage, at most 5% edits and an endpoint within 10
bases of the junction-facing end. Mapping intervals are extended by 15 bases
for arm edit checks. The uninterrupted-background check uses 90% query coverage,
at most twice the query length plus 50 bases, and 30-base interval extensions.
At least 200 arm hits marks the source search as saturated and prevents unique
P assignment. Other source labels within max(50 bases, matching-copy length +
25 bases) are equivalent boundary descriptions.

Validation processes work on complete read groups. Their checkpoint identities
include input alignment metadata, reference sequences/aliases, code and
parameters. The indexed competitor shortlist is covered by the inherited regression tests. Every shortlisted reference still receives the edit-distance test.
