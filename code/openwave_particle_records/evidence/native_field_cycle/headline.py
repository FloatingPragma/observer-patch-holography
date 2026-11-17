#!/usr/bin/env python3
"""Read-only full replay; includes original frozen theorem documentation."""
from pathlib import Path
import runpy
runpy.run_path(str(Path(__file__).with_name('replay.py')),run_name='__main__')
