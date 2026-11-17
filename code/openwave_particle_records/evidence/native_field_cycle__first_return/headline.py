#!/usr/bin/env python3
"""Read-only complete amended first-return replay."""
from pathlib import Path
import runpy
import sys
if len(sys.argv)!=1:raise SystemExit('Read-only replay takes no options')
runpy.run_path(str(Path(__file__).with_name('probe.py')),run_name='__main__')
