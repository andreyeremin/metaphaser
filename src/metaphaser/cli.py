"""metaphaser command-line interface."""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

from rich.console import Console
from rich.table import Table

from . import __version__
from .config import DEFAULTS, find_executable, parameters, runtime_info


def positive(value):
    n=int(value)
    if n<1:raise argparse.ArgumentTypeError("must be greater than zero")
    return n


def parser():
    root=argparse.ArgumentParser(prog="metaphaser",description="Link metagenomic junction states and nanopore DNA modification calls.")
    root.add_argument("--version",action="version",version="metaphaser "+__version__)
    commands=root.add_subparsers(dest="command")
    run=commands.add_parser("run",help="Classify reads and phase modification profiles")
    run.add_argument("--assembly","-a",required=True,help="Metagenome assembly, FASTA or FASTA.gz")
    run.add_argument("--modbam","-b",required=True,help="Original BAM with complete read sequences and MM/ML tags")
    run.add_argument("--outdir","-o",required=True,help="New output directory, or an existing run with --resume")
    run.add_argument("--sample","-s",help="Sample identifier, derived from the BAM name if omitted")
    run.add_argument("--threads","-t",type=positive,default=2,help="Mapping CPU threads (default: %(default)s)")
    run.add_argument("--validation-workers",type=positive,help="Validation processes (default: min(threads, 12))")
    run.add_argument("--minimap2",help="Path to minimap2, otherwise use PATH")
    run.add_argument("--resume",action="store_true",help="Reuse completed stages after checking inputs, code and settings")
    run.add_argument("--hash-inputs",action="store_true",help="Also hash the entire BAM for stronger resume identity")
    run.add_argument("--ignore-original-alignments",action="store_true",help="Use fresh assembly mappings only for discovery")
    run.add_argument("--allow-no-modifications",action="store_true",help="Allow sequence-only BAMs")
    run.add_argument("--modifications",default="auto",help="auto, or comma-separated 5mC,5hmC,6mA")
    run.add_argument("--cytosine-context",choices=["cpg","all"],default="cpg",help="Cytosine positions to compare (default: cpg)")
    run.add_argument("--bootstrap-replicates",type=int,default=0,help="Compatibility option, only 0 is supported by the exact stratified method")
    run.add_argument("--permutations",type=positive,default=19_999,help="Unused compatibility option, exact testing never falls back to sampling")
    run.add_argument("--exact-state-limit",type=positive,default=2_000_000,help="Maximum exact dynamic-programming states")
    run.add_argument("--seed",type=int,default=20261016,help="Recorded compatibility value, the exact test is deterministic")
    run.add_argument("--parameters",help="JSON object overriding documented analysis parameters")
    for name in ["min-length","max-length","min-anchor","max-query-gap","min-support","ml-cutoff"]:
        run.add_argument("--"+name,type=int,help=f"Override {name.replace('-','_')} (default: {DEFAULTS[name.replace('-','_')]})")
    run.add_argument("--no-color",action="store_true",help="Disable terminal colors")
    run.add_argument("--min-discovery-reads",type=positive,help="Preliminary candidate read count, distinct from final --min-support (default: 2)")
    run.add_argument("--debug",action="store_true",help="Show debugging details and Python tracebacks")
    doctor=commands.add_parser("doctor",help="Check Python dependencies and minimap2 without analysing reads")
    doctor.add_argument("--minimap2")
    doctor.add_argument("--json",action="store_true")
    check=commands.add_parser("check",help="Check input format, reference identity and initial modification tags")
    check.add_argument("--assembly",required=True)
    check.add_argument("--modbam",required=True)
    check.add_argument("--ignore-original-alignments",action="store_true")
    demo=commands.add_parser("demo",help="Generate small labelled inputs and run the complete pipeline")
    demo.add_argument("--outdir",required=True)
    demo.add_argument("--threads",type=positive,default=2)
    demo.add_argument("--minimap2")
    demo.add_argument("--no-color",action="store_true")
    demo.add_argument("--length",type=positive,default=800)
    adjust = commands.add_parser("adjust", help="Apply study-wide BH to exact host-linked P comparisons")
    adjust.add_argument("--comparisons", nargs="+", required=True, help="All sample modifications/comparisons.tsv files in the study family")
    adjust.add_argument("--output", required=True, help="New primary comparison TSV with BH q-values")
    adjust.add_argument("--expected-tests", type=positive, help="Fail when the family has a different number of tests")
    return root


def main(argv=None):
    root=parser()
    args=root.parse_args(argv)
    console=Console(stderr=True,highlight=False)
    try:
        if args.command is None:
            root.print_help()
            return 0
        if args.command == "adjust":
            from .multiple_testing import adjust_files
            result = adjust_files(args.comparisons, args.output, args.expected_tests)
            print(json.dumps(result))
            return 0
        if args.command=="doctor":
            mm,version=find_executable(args.minimap2,"minimap2")
            info=runtime_info(version)
            info.update(metaphaser=__version__,minimap2_path=mm)
            if args.json:print(json.dumps(info,indent=2))
            else:
                table=Table(title="metaphaser | Runtime check")
                table.add_column("Component");table.add_column("Version / path")
                for key,value in [("metaphaser",__version__),("Python",info["python"]),*info["packages"].items(),("minimap2",version),("executable",mm)]:table.add_row(key,str(value))
                console.print(table)
                console.print("Ready to run.",style="green")
            return 0
        if args.command=="check":
            from .pipeline import preflight
            _,result=preflight(args.assembly,args.modbam,args.ignore_original_alignments)
            print(json.dumps(result,indent=2))
            return 0
        if args.command=="demo":
            from .demo import make_demo
            folder=Path(args.outdir).expanduser().resolve()
            assembly,bam=make_demo(folder/"inputs",args.length)
            command=["run","--assembly",str(assembly),"--modbam",str(bam),"--outdir",str(folder/"results"),"--sample","demo","--threads",str(args.threads)]
            if args.minimap2:command.extend(["--minimap2",args.minimap2])
            if args.no_color:command.append("--no-color")
            return main(command)
        if args.validation_workers and args.validation_workers>args.threads:
            raise ValueError("--validation-workers cannot exceed --threads.")
        if args.bootstrap_replicates != 0:
            raise ValueError("--bootstrap-replicates must be 0 for the exact stratified method.")
        if not 0<=args.seed<2**32:raise ValueError("--seed must be between 0 and 4294967295.")
        if args.modifications!="auto":
            values=[v.strip() for v in args.modifications.split(",")]
            if not values or set(values)-{"5mC","5hmC","6mA"}:raise ValueError("Use auto, or a comma-separated selection of 5mC,5hmC,6mA.")
            args.modifications=",".join(sorted(set(values)))
        if args.sample is None:args.sample=Path(args.modbam).name.removesuffix(".bam")
        if not re.fullmatch(r"[A-Za-z0-9_.-]+",args.sample):raise ValueError("The sample identifier may contain only letters, digits, dots, underscores and hyphens.")
        p=parameters(args)
        # Avoid BLAS pools multiplying the explicitly requested worker budget.
        for key in ["OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"]:os.environ[key]="1"
        from .pipeline import run
        result=run(args,p)
        print(json.dumps({"status":result["status"],"sample":result["sample"],"outdir":str(Path(args.outdir).resolve()),
                          "host_linked_P_groups":result["host_linked_P_groups"],"all_four_groups":result["all_four_groups"]}))
        return 0
    except KeyboardInterrupt:
        console.print("Interrupted. Completed stages were preserved for --resume.",style="yellow")
        return 130
    except Exception as exc:
        console.print(f"ERROR: {exc}",style="red",markup=False)
        if getattr(args,"debug",False):console.print_exception(show_locals=False)
        return 2


if __name__=="__main__":sys.exit(main())
