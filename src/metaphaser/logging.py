"""English terminal progress plus plain and structured persistent logs."""
from __future__ import annotations

import datetime as dt
import json
import threading
import time
from pathlib import Path

from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, TimeElapsedColumn
from rich.table import Table


class RunLogger:
    def __init__(self, folder, no_color=False, debug=False):
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        self.console = Console(stderr=True, no_color=no_color, highlight=False)
        self.plain = (folder/"run.log").open("a", buffering=1)
        self.jsonl = (folder/"events.jsonl").open("a", buffering=1)
        self.debug = debug
        self.lock = threading.RLock()
        self.current = None
        self.progress = None
        self.last_event = time.monotonic()
        self.started = time.monotonic()
        self.stopping = threading.Event()
        self.thread = threading.Thread(target=self._heartbeat, daemon=True)
        self.thread.start()

    def event(self, message, level="INFO", **details):
        stamp = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
        with self.lock:
            self.last_event = time.monotonic()
            self.plain.write(f"{stamp} {level:<7} {message}\n")
            self.jsonl.write(json.dumps({"time":stamp, "level":level, "message":message,
                                        "elapsed_seconds":round(time.monotonic()-self.started,2),
                                        **details}, default=str)+"\n")
            if level == "DEBUG" and not self.debug:
                return
            style = {"INFO":"cyan", "SUCCESS":"green", "WARNING":"yellow", "ERROR":"red", "DEBUG":"dim"}.get(level, "")
            self.console.print(f"{stamp[11:19]} {level:<7} {message}", style=style, markup=False)

    def core(self, *parts):
        if parts and parts[0] == "map":
            self.event("Mapping reads with minimap2", command=parts[1] if len(parts)>1 else "")
        else:
            self.event(" ".join(map(str, parts)))

    def stage(self, index, total, title):
        self.end_stage()
        self.current = (index, total, title, time.monotonic())
        self.event(f"[{index:02}/{total:02}] {title}", stage=index)
        if self.console.is_terminal:
            self.progress = Progress(SpinnerColumn(style="cyan"), TextColumn("{task.description}"),
                                     TimeElapsedColumn(), console=self.console, transient=True)
            self.progress.start()
            self.progress.add_task(title, total=None)

    def end_stage(self, reused=False):
        if self.progress:
            self.progress.stop()
            self.progress = None
        if self.current:
            _,_,title,started = self.current
            self.event(f"{'Reused' if reused else 'Finished'}: {title} ({time.monotonic()-started:.1f} s)", "SUCCESS")
            self.current = None

    def _heartbeat(self):
        while not self.stopping.wait(5):
            if self.current and time.monotonic()-self.last_event >= 30:
                self.event(f"Still working: {self.current[2]} ({time.monotonic()-self.current[3]:.0f} s)")

    def display_summary(self,result):
        if not self.console.is_terminal:
            return
        table=Table(title="metaphaser | Run complete",border_style="green")
        table.add_column("Result",style="bold")
        table.add_column("Count",justify="right")
        for label,key in [("Input read identifiers","unique_input_reads"),("Host-linked P groups","host_linked_P_groups"),
                          ("Both host boundaries + P","both_boundaries_P_groups"),("All four junctions","all_four_groups"),
                          ("Overlapping regions","regions"),("Tested modification contrasts","tested_comparisons")]:
            table.add_row(label,f"{result[key]:,}")
        self.console.print(table)

    def close(self, success=True):
        self.stopping.set()
        self.thread.join(timeout=1)
        if self.progress:
            self.progress.stop()
        if success:
            self.end_stage()
        self.plain.close()
        self.jsonl.close()
