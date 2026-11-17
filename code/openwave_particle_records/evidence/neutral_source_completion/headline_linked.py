#!/usr/bin/env python3
"""Read-only complete replay with custody checks."""
import hashlib, os, subprocess, sys
from pathlib import Path
p=Path(__file__).resolve().parent
m=p/'MANIFEST.sha256'
if m.exists():
    for row in m.read_text().splitlines():
        digest,name=row.split('  ',1)
        assert hashlib.sha256((p/name).read_bytes()).hexdigest()==digest,name
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
raise SystemExit(subprocess.run([sys.executable,str(p/'probe_linked.py')],env=env).returncode)
