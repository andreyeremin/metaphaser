import json
import os
import shutil
import signal
import subprocess
import sys
import time

import pytest

from metaphaser.demo import make_demo


@pytest.mark.integration
def test_sigterm_preserves_completed_stages_and_allows_resume(tmp_path):
    mm=shutil.which("minimap2")
    if not mm:pytest.skip("minimap2 is not installed")
    assembly,bam=make_demo(tmp_path/"inputs")
    marker=tmp_path/"mapping_started"
    delay=tmp_path/"slow_minimap2"
    delay.write_text(f'''#!/usr/bin/env python3
import os,sys,time
from pathlib import Path
if '--version' not in sys.argv:
    Path({str(marker)!r}).write_text(str(os.getpid()))
    if not os.environ.get('METAPHASER_TEST_FAST'):
        time.sleep(30)
os.execv({mm!r}, [{mm!r}]+sys.argv[1:])
''')
    delay.chmod(0o755)
    out=tmp_path/"results"
    args=[sys.executable,"-m","metaphaser","run","--assembly",str(assembly),"--modbam",str(bam),
          "--outdir",str(out),"--threads","1","--minimap2",str(delay),"--no-color"]
    proc=subprocess.Popen(args,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    for _ in range(100):
        if marker.exists() or proc.poll() is not None:break
        time.sleep(.05)
    assert marker.exists(),proc.communicate(timeout=5)
    child_pid=int(marker.read_text())
    proc.send_signal(signal.SIGTERM)
    _,err=proc.communicate(timeout=10)
    assert proc.returncode==130,err
    manifest=json.loads((out/"manifest.json").read_text())
    assert manifest["status"]=="interrupted"
    assert manifest["stages"]["1"]["status"]=="complete"
    with pytest.raises(ProcessLookupError):os.kill(child_pid,0)
    env=dict(os.environ,METAPHASER_TEST_FAST="1")
    again=subprocess.run(args+["--resume"],env=env,capture_output=True,text=True,timeout=30)
    assert again.returncode==0,again.stderr
    assert "Reused: Export complete reads" in again.stderr
