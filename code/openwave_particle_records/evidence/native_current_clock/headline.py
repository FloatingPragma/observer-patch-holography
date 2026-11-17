#!/usr/bin/env python3
"""Read-only full replay with fixed linear-algebra threading and custody."""
import hashlib
import os
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
manifest=HERE/'MANIFEST.sha256'
if manifest.exists():
    for line in manifest.read_text().splitlines():
        expected,name=line.split('  ',1)
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==expected,name
env=dict(os.environ)
env['OPENBLAS_NUM_THREADS']='1'
env['OMP_NUM_THREADS']='1'
for script in ('probe.py','interval_audit.py'):
    result=subprocess.run([sys.executable,str(HERE/script)],env=env)
    if result.returncode:raise SystemExit(result.returncode)
