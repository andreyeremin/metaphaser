# Discovery changes in 0.2.0rc2

The preliminary count of split-read events is now separate from the final
count of validated junction reads. The candidate screen uses two distinct
read identifiers by default (`--min-discovery-reads 2`). Retained host-linked
intervals still require at least three validated P reads and at least three
validated L or R reads. The discovery threshold never increases the number
of validated reads and does not change the identity, anchor, competition,
assembly-background, or source-copy checks.

For a positive gap between two split alignments, the sequence in that gap
may extend the interval head or tail beyond the alignment endpoint. The
program tests extensions of either cut by up to the observed read-gap
length, within the existing interval-length and contig limits. It compares
the same read window to each candidate P reference with Edlib. A candidate
must improve on the best original boundary hypothesis and the alignment
must cross the join with at least 75 reference bases on either side.

The original observed hypotheses remain present. Up to five additional
boundary pairs are retained per group, ranked by the number of distinct
reads giving an improvement, total edit-distance improvement, minimum edit
distance and coordinates. These are proposals only. Reads are subsequently
remapped and subjected to the unchanged full junction validation. No sample
names, read identifiers or known locus coordinates are present in the
implementation.

`work/groups.json` records observed and refined hypotheses and the read-level
evidence for each refinement. `work/discovery_refinement.json` records the
candidate and final thresholds and the number of added hypotheses. The
software version and source digest change, so older runs cannot be silently
resumed with the new algorithm.

This changes candidate generation. It is not a claim of identical results to
0.2.0rc1 or of exhaustive sensitivity. The exact modification test and BH
implementation are unchanged.
