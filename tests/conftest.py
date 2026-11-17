import os
import sys
from pathlib import Path

import pytest


@pytest.fixture(autouse=True)
def use_test_runtime_executables(monkeypatch):
    monkeypatch.setenv("PATH",str(Path(sys.executable).parent)+os.pathsep+os.environ.get("PATH",""))
