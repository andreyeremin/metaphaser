"""Command-line orchestration around the preserved junction validator."""
from __future__ import annotations

import collections
import datetime as dt
import hashlib
import json
import multiprocessing
import os
import signal
import sqlite3
import subprocess
import sys
import time
from pathlib import Path

import pysam

from . import __version__
from ._vendor import junctions_v17 as core
from .config import code_identity, file_identity, find_executable, runtime_info
from .io import json_read, json_write, read_fasta, run_lock
from .logging import RunLogger
from .modifications import extract, compare
from .reporting import classify_loci, public_tables, html_report
from .discovery import discover_fresh, group_and_refine


class Processes:
    """Own every mapping child, including the preserved source-copy audit."""
    PIPE=subprocess.PIPE
    def __init__(self):self.children=[]
    def Popen(self,*args,**kwargs):
        kwargs.setdefault("start_new_session",True)
        proc=subprocess.Popen(*args,**kwargs)
        self.children.append(proc)
        return proc
    @staticmethod
    def run(*args,**kwargs):return subprocess.run(*args,**kwargs)
    def stop(self):
        for proc in self.children:
            if proc.poll() is None:
                try:os.killpg(proc.pid,signal.SIGTERM)
                except ProcessLookupError:pass
                try:proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid,signal.SIGKILL)
                    proc.wait()


def preflight(assembly,modbam,ignore_original=False):
    seqs=read_fasta(assembly)
    if not core.complete_bam(modbam):
        raise ValueError("The modBAM is empty or missing its BGZF end marker. Finish or repair the BAM before running metaphaser.")
    with pysam.AlignmentFile(modbam,"rb",check_sq=False) as bam:
        mismatch=[(name,length) for name,length in zip(bam.references,bam.lengths)
                  if name not in seqs or len(seqs[name])!=length]
        if mismatch and not ignore_original:
            raise ValueError("BAM reference names or lengths do not match the assembly. Use the matching assembly, or --ignore-original-alignments to discover junctions from a fresh mapping. First mismatch: "+str(mismatch[0]))
        tagged=0
        primary=0
        headers=set()
        for n,a in enumerate(bam.fetch(until_eof=True),1):
            if not a.is_secondary and not a.is_supplementary:
                primary+=1
                if a.has_tag("MM") and a.has_tag("ML"):
                    tagged+=1
                    headers.update(x.split(",")[0] for x in a.get_tag("MM").split(";") if x)
            if n>=10_000:break
    return seqs,dict(assembly_contigs=len(seqs),assembly_bases=sum(map(len,seqs.values())),
                    sampled_primary_reads=primary,sampled_tagged_reads=tagged,
                    modification_groups=sorted(headers),ignored_original_alignments=ignore_original)


def export_reads(sample,seqs,work,p,logger,ignore_original=False,allow_no_modifications=False):
    inventory=work/"read_inventory.sqlite"
    if inventory.exists():inventory.unlink()
    db=sqlite3.connect(inventory)
    db.execute("CREATE TABLE reads(read_id TEXT PRIMARY KEY,length INTEGER,quality REAL,eligible INTEGER,reason TEXT,has_mm INTEGER,has_ml INTEGER,records INTEGER)")
    counts=collections.Counter()
    events=[]
    seen=set()
    modes=collections.Counter()
    with pysam.AlignmentFile(sample["bam"],"rb",check_sq=False) as bam,(work/"reads.fasta").open("w") as fasta:
        bam_header=bam.header.to_dict()
        for a in bam.fetch(until_eof=True):
            counts["records"]+=1
            if a.is_secondary or a.is_supplementary:continue
            counts["primary_records"]+=1
            rid=a.query_name
            if not rid:raise ValueError("A primary BAM record has no read identifier.")
            if any(c.isspace() for c in rid):raise ValueError(f"Read identifier contains whitespace: {rid!r}")
            mm,ml=a.has_tag("MM"),a.has_tag("ML")
            quality=float(a.get_tag("qs")) if a.has_tag("qs") else \
                sum(a.query_qualities)/len(a.query_qualities) if a.query_qualities is not None and len(a.query_qualities) else None
            reason=""
            if a.is_qcfail or a.is_duplicate:reason="qc_or_duplicate_flag"
            elif rid in seen:reason="duplicate_primary_identifier"
            else:
                seen.add(rid)
                if not a.query_sequence:reason="missing_sequence"
                elif any(op==5 for op,n in a.cigartuples or []):reason="hard_clipped_sequence"
            eligible=not reason
            # Prefer the record actually exported if an earlier record was flag-filtered.
            db.execute("""INSERT INTO reads VALUES(?,?,?,?,?,?,?,1) ON CONFLICT(read_id) DO UPDATE SET
                length=CASE WHEN excluded.eligible>reads.eligible THEN excluded.length ELSE reads.length END,
                quality=CASE WHEN excluded.eligible>reads.eligible THEN excluded.quality ELSE reads.quality END,
                reason=CASE WHEN excluded.eligible>reads.eligible THEN excluded.reason ELSE reads.reason END,
                has_mm=CASE WHEN excluded.eligible>reads.eligible THEN excluded.has_mm ELSE reads.has_mm END,
                has_ml=CASE WHEN excluded.eligible>reads.eligible THEN excluded.has_ml ELSE reads.has_ml END,
                eligible=MAX(reads.eligible,excluded.eligible),records=reads.records+1""",
                (rid,a.query_length or 0,quality,int(eligible),reason,int(mm),int(ml)))
            if reason:
                counts[reason]+=1
                continue
            seq=a.get_forward_sequence()
            fasta.write(">"+rid+"\n"+seq+"\n")
            counts["exported_reads"]+=1
            counts["exported_bases"]+=len(seq)
            if mm and ml:
                counts["tagged_exported_reads"]+=1
                for group in a.get_tag("MM").split(";"):
                    if group:modes[group.split(",")[0]]+=1
            if not ignore_original:
                events.extend(core.split_events(rid,core.segments_of(a,True),p,"original_SA"))
                events.extend(core.cigar_events(a,seqs,p))
            if counts["exported_reads"]%25_000==0:
                db.commit()
                logger.event(f"Exported {counts['exported_reads']:,} reads ({counts['exported_bases']:,} bases)")
    db.commit()
    counts["unique_input_reads"]=db.execute("SELECT COUNT(*) FROM reads").fetchone()[0]
    db.close()
    if not counts["exported_reads"]:raise ValueError("No complete eligible primary read sequences were found.")
    if not counts["tagged_exported_reads"] and not allow_no_modifications:
        raise ValueError("No exported read contains both MM and ML tags. Supply a modBAM or explicitly use --allow-no-modifications for sequence-only analysis.")
    json_write(work/"original_events.json.gz",events)
    value=dict(counts=dict(counts),modification_groups=dict(modes),bam_header=bam_header)
    json_write(work/"source.json",value)
    json_write(work.parent/"input_metadata.json",dict(sample=sample["sample"],bam_header=bam_header,
        modification_groups=dict(modes),assembly_contigs=len(seqs),assembly_bases=sum(map(len,seqs.values()))))
    return value


class Stages:
    def __init__(self,folder,manifest,logger,resume):
        self.folder=folder
        self.manifest=manifest
        self.logger=logger
        self.resume=resume
        self.invalidated=False

    def run(self,index,title,products,action):
        self.logger.stage(index,12,title)
        key=str(index)
        saved=self.manifest["stages"].get(key)
        reuse=self.resume and not self.invalidated and saved and saved.get("status")=="complete"
        if reuse:
            for name,identity in saved["products"].items():
                path=self.folder/name
                if not path.is_file() or file_identity(path,full_hash=bool(identity.get("sha256")))!=identity:
                    reuse=False
                    break
                if path.suffix==".bam" and not core.complete_bam(path):reuse=False;break
        if reuse:
            self.logger.end_stage(reused=True)
            return
        self.invalidated=True
        if saved:self.logger.event("Stage output changed or is incomplete. Rebuilding this and following stages.","WARNING")
        self.manifest["current_stage"]=title
        self.manifest["stages"][key]={"status":"running","started":dt.datetime.now(dt.timezone.utc).isoformat()}
        json_write(self.folder/"manifest.json",self.manifest)
        action()
        identities={}
        for name in products:
            path=self.folder/name
            if not path.is_file():raise RuntimeError(f"Stage did not produce {path}")
            identities[name]=file_identity(path,full_hash=path.suffix not in {".bam",".sqlite"})
        self.manifest["stages"][key]={"status":"complete","products":identities}
        json_write(self.folder/"manifest.json",self.manifest)
        self.logger.end_stage()


def run(options,p):
    folder=Path(options.outdir).expanduser().resolve()
    assembly=Path(options.assembly).expanduser().resolve(strict=True)
    bam=Path(options.modbam).expanduser().resolve(strict=True)
    if assembly.is_relative_to(folder) or bam.is_relative_to(folder):
        raise ValueError("Input files must be outside the output directory.")
    if folder.exists() and any(folder.iterdir()) and not options.resume:
        raise ValueError(f"Output directory is not empty: {folder}. Choose a new directory or add --resume.")
    if options.resume and not (folder/"manifest.json").is_file():
        raise ValueError("--resume needs an existing metaphaser manifest.json.")
    minimap2,mm_version=find_executable(options.minimap2,"minimap2")
    runtime=runtime_info(mm_version)
    sample=dict(sample=options.sample,bam=str(bam),fasta=str(assembly))
    identity=dict(inputs={"assembly":file_identity(assembly,True),"modbam":file_identity(bam,options.hash_inputs)},
        sample=options.sample,parameters=p,modifications=options.modifications,cytosine_context=options.cytosine_context,
        ignore_original_alignments=options.ignore_original_alignments,allow_no_modifications=options.allow_no_modifications,
        random_seed=options.seed,permutations=options.permutations,exact_state_limit=options.exact_state_limit,
        bootstrap_replicates=options.bootstrap_replicates,code_sha256=code_identity(),
        minimap2={"path":minimap2,"version":mm_version,"sha256":hashlib.sha256(Path(minimap2).read_bytes()).hexdigest()},
        runtime=runtime)
    folder.mkdir(parents=True,exist_ok=True)
    with run_lock(folder):
        if options.resume:
            manifest=json_read(folder/"manifest.json")
            if manifest["identity"]!=identity:
                raise ValueError("Resume identity changed: inputs, code, dependencies or analysis parameters differ. Use a new output directory. Threads may be changed without changing the analysis identity.")
        else:
            manifest={"tool":"metaphaser","version":__version__,"identity":identity,"stages":{},"status":"created",
                      "created":dt.datetime.now(dt.timezone.utc).isoformat(),"coordinates":"zero_based_half_open"}
        logger=RunLogger(folder/"logs",options.no_color,options.debug)
        owner=Processes()
        previous_term=signal.getsignal(signal.SIGTERM)
        def terminate(signum,frame):raise KeyboardInterrupt("SIGTERM received")
        signal.signal(signal.SIGTERM,terminate)
        if hasattr(os,"register_at_fork"):
            os.register_at_fork(after_in_child=lambda:signal.signal(signal.SIGTERM,signal.SIG_DFL))
        core.log=logger.core
        core.subprocess=owner
        cfg={"minimap2":minimap2,"junction_seed_length":11,"parameters":p,"indexed_competition":True,
             "ranked_candidate_search":True,"verify_interval_bounds":True,
             "sparse_modification_pairs":options.cytosine_context=="cpg",
             "validation_workers":options.validation_workers or min(options.threads,12),
             "validation_checkpoint_dir":str(folder/"work/validation_checkpoints")}
        if "fork" not in multiprocessing.get_all_start_methods():cfg["validation_workers"]=1
        work=folder/"work"
        work.mkdir(exist_ok=True)
        mods=folder/"modifications"
        mods.mkdir(exist_ok=True)
        success=False
        manifest["status"]="running"
        json_write(folder/"manifest.json",manifest)
        started=time.monotonic()
        try:
            logger.event(f"metaphaser {__version__} | sample {options.sample}")
            logger.event(f"Mapping threads: {options.threads}. Validation workers: {cfg['validation_workers']}.")
            logger.event("Checking FASTA, BAM reference coordinates and modification tags")
            seqs,input_check=preflight(assembly,bam,options.ignore_original_alignments)
            logger.event(f"Assembly: {len(seqs):,} contigs, {sum(map(len,seqs.values())):,} bases")
            if not input_check["sampled_tagged_reads"]:
                logger.event("No MM/ML pair was found in the first 10,000 records. The export stage will check all eligible primary reads.","WARNING")
            stages=Stages(folder,manifest,logger,options.resume)

            def mapper(config,reference,reads,out,logfile,threads,small=False):
                cmd=core.mm_argv(config,threads,small)+[str(reference),str(reads)]
                logger.event("Running minimap2",command=cmd)
                progress=Path(str(out)+".status.json")
                if progress.exists():progress.unlink()
                with open(logfile,"w") as err:
                    proc=owner.Popen(cmd,stdout=subprocess.PIPE,stderr=err)
                    converter=owner.Popen([sys.executable,"-m","metaphaser.worker","convert","--output",str(out),
                        "--status",str(progress),"--threads",str(max(1,min(threads,4)))],stdin=proc.stdout,stderr=err)
                    proc.stdout.close()
                    last=None
                    while proc.poll() is None or converter.poll() is None:
                        if (proc.poll() not in [None,0]) or (converter.poll() not in [None,0]):
                            raise RuntimeError(f"Alignment/conversion failed. See {logfile}")
                        if progress.exists():
                            current=json_read(progress)
                            if current!=last:
                                if current["operation"]=="sort":logger.event(f"Grouping {current['records']:,} alignments by read identifier")
                                elif current["records"]:logger.event(f"Converted {current['records']:,} alignments")
                                last=current
                        time.sleep(.1)
                    if proc.returncode or converter.returncode:raise RuntimeError(f"Alignment/conversion failed. See {logfile}")
                return cmd
            core.map_reads=mapper

            stages.run(1,"Export complete reads and inspect original alignments",
                ["work/source.json","work/reads.fasta","work/read_inventory.sqlite","work/original_events.json.gz","input_metadata.json"],
                lambda:export_reads(sample,seqs,work,p,logger,options.ignore_original_alignments,options.allow_no_modifications))
            source=json_read(work/"source.json")
            cfg["keep_A_modification_pairs"]=(options.modifications=="auto" or "6mA" in options.modifications.split(",")) and \
                any(name.startswith("A+a") for name in source["modification_groups"])
            stages.run(2,"Map reads to the metagenome assembly",["work/assembly.bam"],
                lambda:mapper(cfg,assembly,work/"reads.fasta",work/"assembly.bam",folder/"logs/assembly_mapping.log",options.threads))

            def discovery():
                events=json_read(work/"original_events.json.gz")
                events.extend(discover_fresh(work/"assembly.bam",seqs,p))
                groups,discovery_summary=group_and_refine(events,seqs,p)
                json_write(work/"groups.json",groups)
                json_write(work/"discovery_counts.json",dict(events=len(events),groups=len(groups)))
                json_write(work/"discovery_refinement.json",discovery_summary)
                logger.event(f"Discovered {len(events):,} connections in {len(groups):,} retained groups")
                logger.event(f"Candidate screen: {discovery_summary['discovery_min_reads']} discovery reads; final structural threshold: {p['min_support']} validated reads")
                logger.event(f"Refined positive read gaps in {discovery_summary['refined_groups']} groups; added {discovery_summary['added_hypotheses']} boundary hypotheses")
            stages.run(3,"Find and group sequence connections",["work/groups.json","work/discovery_counts.json","work/discovery_refinement.json"],discovery)
            groups=json_read(work/"groups.json")

            def references():
                refs=core.build_references(groups,seqs,p)
                json_write(work/"references.json.gz",refs)
                core.fasta_write(folder/"junction_references.fasta",((r["id"],r["seq"]) for r in refs))
                json_write(folder/"reference_origins.json.gz",[dict(id=r["id"],length=len(r["seq"]),aliases=r["aliases"]) for r in refs])
                logger.event(f"Built {len(refs):,} synthetic junction references")
            stages.run(4,"Build L, R, P and B reference sequences",
                ["work/references.json.gz","junction_references.fasta","reference_origins.json.gz"],references)
            refs=json_read(work/"references.json.gz")

            def recruit():
                n=core.recruit_reads(work/"assembly.bam",groups,work/"recruited.fasta",p["external_flank"])
                if refs and n:mapper(cfg,folder/"junction_references.fasta",work/"recruited.fasta",work/"synthetic.bam",folder/"logs/junction_mapping.log",options.threads,True)
                else:
                    with pysam.AlignmentFile(work/"synthetic.bam","wb",header={"HD":{"VN":"1.6"},"SQ":[]}):pass
                json_write(work/"recruitment.json",dict(reads=n))
                logger.event(f"Recruited {n:,} reads around the discovered boundaries")
            stages.run(5,"Recruit reads and map to synthetic junctions",["work/recruited.fasta","work/recruitment.json","work/synthetic.bam"],recruit)
            recruited=json_read(work/"recruitment.json")["reads"]

            def validation():
                checkpoint=Path(cfg["validation_checkpoint_dir"])
                if checkpoint.exists():
                    # Existing chunk checkpoints are reused only if their complete identity matches.
                    current=core.checkpoint_identity(cfg,refs,work/"synthetic.bam",p) if refs and recruited else None
                    saved=checkpoint/"partitions.json"
                    if saved.exists() and json_read(saved).get("identity")!=current:
                        checkpoint.rename(checkpoint.with_name("validation_checkpoints.stale."+str(time.time_ns())))
                calls,rejections=core.validate_synthetic(cfg,refs,work/"synthetic.bam",p) if refs and recruited else ([],{})
                json_write(work/"validated_calls.json.gz",calls)
                json_write(work/"rejections.json",rejections)
                logger.event(f"Retained {len(calls):,} junction calls before genomic-copy checks")
            stages.run(6,"Validate junction windows and compare alternative configurations",
                ["work/validated_calls.json.gz","work/rejections.json"],validation)

            def background():
                cc=json_read(work/"validated_calls.json.gz")
                core.background_check(cfg,sample,cc,seqs,work,options.threads)
                restored=core.normalize_occupied_labels(cc,refs)
                json_write(work/"background_calls.json.gz",cc)
                logger.event(f"Checked uninterrupted assembly alternatives. Retained {restored:,} same-placement L/R labels.")
            stages.run(7,"Check assembly alternatives and compatible host-boundary labels",["work/background_calls.json.gz"],background)

            def copy_audit():
                request=work/"copy_audit_request.json"
                json_write(request,dict(config=cfg,sample=sample,work=str(work),threads=options.threads,parameters=p))
                logfile=folder/"logs/source_copy_audit.log"
                with logfile.open("w") as log:
                    proc=owner.Popen([sys.executable,"-m","metaphaser.worker","audit","--request",str(request)],stdout=log,stderr=log)
                    while proc.poll() is None:time.sleep(.1)
                    if proc.returncode:raise RuntimeError(f"Source-copy audit failed. See {logfile}")
                audit=json_read(work/"copy_audit_summary.json")
                logger.event(f"Source-copy check: {audit['tested_P_calls']:,} P calls examined, {audit['ambiguous_P_calls']:,} ambiguous")
            stages.run(8,"Check other possible source copies",["work/final_calls.json.gz","work/copy_audit_summary.json"],copy_audit)
            calls=json_read(work/"final_calls.json.gz")
            loci=core.loci_summary(groups,calls,options.sample)
            loci,regions=classify_loci(loci,seqs,p["min_support"])

            def tables():
                classes=public_tables(folder,loci,regions,calls,seqs,options.sample,work/"read_inventory.sqlite")
                json_write(work/"read_classes.json",classes)
                json_write(folder/"loci.json",loci)
                logger.event(f"Classified {sum(classes.values()):,} input read identifiers")
            stages.run(9,"Write read classifications and structural-state tables",
                ["loci.tsv","elements.tsv","regions.tsv","loci.json","element_intervals.fasta","read_classifications.tsv.gz",
                 "read_locus_classifications.tsv.gz","junction_calls.tsv.gz","work/read_classes.json"],tables)

            def modification_calls():
                ob,info=extract(sample,calls,seqs,mods,p,logger,options.cytosine_context,options.modifications)
                records=[dict(key=list(k),values=v) for k,v in ob.items()]
                json_write(work/"modification_observations.json.gz",records)
                json_write(mods/"channels.json",info)
                logger.event("Measured modification channels: "+(", ".join(info["measured_channels"]) or "none"))
            stages.run(10,"Link original modification calls to structural states",
                ["work/modification_observations.json.gz","modifications/site_calls.tsv.gz","modifications/channels.json"],modification_calls)

            def statistics():
                records=json_read(work/"modification_observations.json.gz")
                ob={tuple(r["key"]):{int(k):v for k,v in r["values"].items()} for r in records}
                contrasts=compare(ob,loci,mods,options,logger)
                json_write(work/"contrasts.json.gz",contrasts)
                logger.event(f"Wrote {len(contrasts):,} modification comparisons, {sum(r['test_status']=='tested' for r in contrasts):,} tested")
            stages.run(11,"Phase modification patterns and compare matched positions",
                ["modifications/comparisons.tsv","modifications/read_profiles.tsv.gz","modifications/epihaplotypes.tsv.gz",
                 "modifications/panels.tsv.gz","work/contrasts.json.gz"],statistics)
            contrasts=json_read(work/"contrasts.json.gz")

            def report():
                elements=[r for r in loci if r["host_linked_P"]]
                summary=dict(tool="metaphaser",version=__version__,sample=options.sample,status="complete",
                    unique_input_reads=source["counts"]["unique_input_reads"],input=source["counts"],
                    host_linked_P_groups=len(elements),both_boundaries_P_groups=sum(r["both_host_boundaries"] for r in elements),
                    all_four_groups=sum(r["all_four_junctions"] for r in elements),regions=len(regions),
                    all_loci=len(loci),structural_classes=dict(collections.Counter(r["structural_class"] for r in loci)),
                    read_classes=json_read(work/"read_classes.json"),modifications=json_read(mods/"channels.json"),
                    tested_comparisons=sum(r["test_status"]=="tested" for r in contrasts),
                    p_values="unadjusted",parameters=p,elapsed_seconds=round(time.monotonic()-started,2))
                json_write(folder/"summary.json",summary)
                html_report(folder/"report.html",summary,loci,contrasts)
            stages.run(12,"Write the summary and HTML report",["summary.json","report.html"],report)
            summary=json_read(folder/"summary.json")
            manifest.update(status="complete",finished=dt.datetime.now(dt.timezone.utc).isoformat(),current_stage=None,
                            last_requested_threads=options.threads)
            json_write(folder/"manifest.json",manifest)
            logger.event(f"Complete: {summary['host_linked_P_groups']:,} host-linked P groups, {summary['both_boundaries_P_groups']:,} with both host boundaries, {summary['all_four_groups']:,} with all four junctions", "SUCCESS")
            logger.event(f"Report: {folder/'report.html'}", "SUCCESS")
            logger.display_summary(summary)
            success=True
            return summary
        except BaseException as exc:
            owner.stop()
            manifest.update(status="interrupted" if isinstance(exc,KeyboardInterrupt) else "failed",error=str(exc) or type(exc).__name__)
            json_write(folder/"manifest.json",manifest)
            logger.event(f"{manifest['status'].capitalize()}: {manifest['error']}","ERROR")
            logger.event("Completed stage outputs were preserved. Fix the problem and repeat the same command with --resume.")
            raise
        finally:
            owner.stop()
            logger.close(success)
            signal.signal(signal.SIGTERM,previous_term)
