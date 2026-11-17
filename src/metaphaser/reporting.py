"""Stable public tables, interval sequences and a self-contained HTML report."""
from __future__ import annotations

import collections
import html
import json
import sqlite3
from pathlib import Path

from .io import write_tsv


def classify_loci(rows,seqs,threshold):
    for r in rows:
        L,R,P,B=(r[k+"_simple"]>=threshold for k in "LRPB")
        if P and L and R:
            name="host_linked_P_both_boundaries"
        elif P and (L or R):
            name="host_linked_P_one_boundary"
        elif P and r["a"]<=25 and len(seqs[r["contig"]])-r["b"]<=25:
            name="contig_end_join"
        elif P:
            name="P_without_host_boundary"
        elif B:
            name="B_without_P"
        else:
            name="below_read_threshold"
        r.update(structural_class=name,empty_site_detected=B,
                 both_host_boundaries=L and R,host_linked_P=P and (L or R),
                 all_four_junctions=L and R and P and B,region="",
                 contig_length=len(seqs[r["contig"]]))
    regions=[]
    by_contig=collections.defaultdict(list)
    for r in rows:
        if r["host_linked_P"]:by_contig[r["contig"]].append(r)
    for contig,members in sorted(by_contig.items()):
        current=None
        for r in sorted(members,key=lambda x:(x["a"],x["b"],x["locus"])):
            if current is None or r["a"]>=current["end"]:
                current=dict(region=f"REG{len(regions)+1:04}",contig=contig,start=r["a"],end=r["b"],loci=[])
                regions.append(current)
            current["end"]=max(current["end"],r["b"])
            current["loci"].append(r["locus"])
            r["region"]=current["region"]
    return rows,regions


def public_tables(folder,rows,regions,calls,seqs,sample,inventory):
    fields=["sample","locus","region","contig","a","b","length","short","structural_class",
        "L_simple","R_simple","P_simple","B_simple","L_all","R_all","P_all","B_all",
        "complex_reads","ambiguous_reads","genomic_alternative_reads","host_linked_P","both_host_boundaries",
        "empty_site_detected","all_four_junctions","discovery_reads","contig_length","model","repeat_length",
        "L_500","R_500","P_500","B_500","L_1000","R_1000","P_1000","B_1000"]
    write_tsv(folder/"loci.tsv",rows,fields)
    elements=[r for r in rows if r["host_linked_P"]]
    write_tsv(folder/"elements.tsv",elements,fields)
    write_tsv(folder/"regions.tsv",regions,["region","contig","start","end","loci"])
    with (folder/"element_intervals.fasta").open("w") as f:
        for r in elements:
            seq=seqs[r["contig"]][r["a"]:r["b"]]
            f.write(f'>{sample}:{r["locus"]} contig={r["contig"]} start={r["a"]} end={r["b"]} coordinates=zero_based_half_open\n')
            for start in range(0,len(seq),80):f.write(seq[start:start+80]+"\n")
    columns=["read_id","locus","contig","kind","simple","complex_read","locus_ambiguous",
        "genomic_contiguous_alternative","source_copy_ambiguous","source_search_saturated","a","b",
        "model","repeat_length","reverse","identity","window","edit_distance","required_margin",
        "competing_distance","competing_loci","forward_q0","forward_q1","tiers"]
    columns += ["source_copy_alternatives","background_competitors","overlapping_occupied_locus_labels","distinct_source_locus_labels"]
    write_tsv(folder/"junction_calls.tsv.gz",calls,columns)
    by_read=collections.defaultdict(list)
    for c in calls:by_read[c["read_id"]].append(c)
    per_locus=[]
    for rid,values in sorted(by_read.items()):
        grouped=collections.defaultdict(list)
        for c in values:grouped[c["locus"]].append(c)
        for lid,cc in grouped.items():
            good={c["kind"] for c in cc if c["simple"]}
            reasons=set()
            if any(c["locus_ambiguous"] for c in cc):reasons.add("source_ambiguous")
            if any(c["complex_read"] for c in cc):reasons.add("complex_junctions")
            if any(c["genomic_contiguous_alternative"] for c in cc):reasons.add("assembly_explained")
            state="integrated" if good and good<={"L","R"} else "+".join(sorted(good)) if good else "unassigned"
            per_locus.append(dict(sample=sample,read_id=rid,locus=lid,state=state,
                                  junctions=sorted(good),observed_junctions=sorted({c["kind"] for c in cc}),reasons=sorted(reasons)))
    write_tsv(folder/"read_locus_classifications.tsv.gz",per_locus,
        ["sample","read_id","locus","state","junctions","observed_junctions","reasons"])
    locus_by_read=collections.defaultdict(list)
    for row in per_locus:locus_by_read[row["read_id"]].append(row)
    classes=collections.Counter()
    def read_rows():
        with sqlite3.connect(inventory) as db:
            for rid,length,quality,eligible,reason,mm,ml,records in db.execute(
                    "SELECT read_id,length,quality,eligible,reason,has_mm,has_ml,records FROM reads ORDER BY rowid"):
                entries=locus_by_read.get(rid,[])
                accepted=[x for x in entries if x["junctions"]]
                if not eligible:label="input_filtered"
                elif accepted:
                    states={x["state"] for x in accepted}
                    label=next(iter(states)) if len(states)==1 else "multiple_states"
                elif any("complex_junctions" in x["reasons"] for x in entries):label="complex"
                elif any("source_ambiguous" in x["reasons"] for x in entries):label="ambiguous"
                elif any("assembly_explained" in x["reasons"] for x in entries):label="assembly_explained"
                else:label="unassigned"
                classes[label]+=1
                yield dict(sample=sample,read_id=rid,classification=label,length=length,mean_q=quality,
                    exported=bool(eligible),input_filter=reason,has_MM=bool(mm),has_ML=bool(ml),
                    primary_records=records,locus_assignments=entries)
    write_tsv(folder/"read_classifications.tsv.gz",read_rows(),
        ["sample","read_id","classification","length","mean_q","exported","input_filter","has_MM","has_ML","primary_records","locus_assignments"])
    return dict(classes)


def html_report(path,summary,loci,contrasts):
    esc=lambda x:html.escape(str(x),quote=True)
    elements=[r for r in loci if r["host_linked_P"]]
    tests=[r for r in contrasts if r["test_status"]=="tested" and r["locus"] in {e["locus"] for e in elements}]
    tests.sort(key=lambda r:r["p_value"])
    def table(headers,rows):
        return '<div class="scroll"><table><thead><tr>'+''.join('<th>'+esc(h)+'</th>' for h in headers)+ \
            '</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+esc(v)+'</td>' for v in row)+'</tr>' for row in rows)+'</tbody></table></div>'
    overview=table(["Locus","Region","Contig","Start","End","Length","L","R","P","B","Class"],
        [[r[k] for k in ["locus","region","contig","a","b","length","L_simple","R_simple","P_simple","B_simple","structural_class"]] for r in elements[:500]])
    comparisons=table(["Locus","Contrast","Modification","Strand","Baseline reads","State reads","Sites","Difference (pp)","p","Test"],
        [[r["locus"],r["contrast"],r["modification"],r["strand"],r["n_baseline"],r["n_state"],r["shared_positions"],
          f'{r["delta_pp"]:.3g}',f'{r["p_value"]:.3E}',r["test"]] for r in tests[:100]])
    links=[("elements.tsv","Host-linked P loci"),("loci.tsv","All discovered loci"),("regions.tsv","Overlapping regions"),
        ("read_classifications.tsv.gz","All read classifications"),("read_locus_classifications.tsv.gz","Assignments by read and locus"),
        ("junction_calls.tsv.gz","Junction-level calls"),("element_intervals.fasta","Assembly interval sequences"),
        ("modifications/comparisons.tsv","Modification comparisons"),("modifications/read_profiles.tsv.gz","Per-read modification profiles"),
        ("modifications/site_calls.tsv.gz","Per-position calls"),("modifications/epihaplotypes.tsv.gz","Modification patterns"),
        ("modifications/panels.tsv.gz","Shared genomic positions"),("manifest.json","Parameters and provenance"),("logs/run.log","Run log")]
    cards=[("Input reads",summary["unique_input_reads"]),("Host-linked P groups",len(elements)),
           ("Both host boundaries + P",sum(r["both_host_boundaries"] for r in elements)),
           ("All four junctions",sum(r["all_four_junctions"] for r in elements))]
    body=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>metaphaser | {esc(summary['sample'])}</title><style>
body{{font:16px system-ui,sans-serif;color:#193c4a;max-width:1400px;margin:36px auto;padding:0 24px;background:#fafcfd}}
h1{{font-size:36px;margin-bottom:6px}}h2{{margin-top:38px}}.sub{{color:#5e7681}}.cards{{display:flex;flex-wrap:wrap;gap:16px;margin:28px 0}}
.card{{background:white;border:1px solid #d7e2e7;padding:20px;min-width:180px;flex:1}}.num{{font-size:32px;color:#007f83;font-weight:700}}
.scroll{{overflow:auto}}table{{border-collapse:collapse;width:100%;background:white}}td,th{{padding:10px 12px;border-bottom:1px solid #d7e2e7;text-align:left;white-space:nowrap}}
th{{background:#eaf5f3}}a{{color:#006f79}}li{{margin:8px 0}}code{{background:#edf2f5;padding:2px 5px}}footer{{margin:40px 0;color:#5e7681}}
</style><h1>metaphaser</h1><div class="sub">{esc(summary['sample'])} | version {esc(summary['version'])} | complete</div>
<div class="cards">{''.join('<div class="card"><div class="num">'+esc(v)+'</div>'+esc(k)+'</div>' for k,v in cards)}</div>
<p><b>L</b>: left host boundary. <b>R</b>: right host boundary. <b>P</b>: end-to-start connection. <b>B</b>: flank-to-flank connection.</p>
<p>Counts use simple, uniquely assigned read identifiers. A read spanning L and R appears in both counts. Overlapping locus groups share a region identifier.</p>
<h2>Host-linked P loci</h2>{overview}<p class="sub">Showing up to 500 groups. The tables contain all results.</p>
<h2>Modification comparisons</h2><p>One observation per read, identical shared positions, state labels permuted within orientation, with equal orientation weights. p-values are unadjusted. A missing test means the read or position requirement was not met.</p>
{comparisons}<p class="sub">Showing up to 100 tested comparisons, ordered by p-value. The full table also includes untested comparisons and their reasons.</p>
<h2>How to read the modification patterns</h2><p>Each pattern uses the ordered positions in <code>panels.tsv.gz</code>. A 1 denotes a high-confidence call and a 0 a usable call below the cutoff. Unknown calls are excluded from the shared panel. The separate 5mC and 5hmC patterns are retained.</p>
<h2>Files</h2><ul>{''.join('<li><a href="'+esc(name)+'">'+esc(label)+'</a></li>' for name,label in links)}</ul>
<footer>Coordinates are zero-based and half-open. Interval FASTA sequences come from the assembly. Structural classes describe observed junctions and do not assign a mobile-element gene family or physical topology.</footer></html>'''
    Path(path).write_text(body)
