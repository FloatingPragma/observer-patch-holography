#!/usr/bin/env python3
"""Run the two self-contained offline components in separate interpreters."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
HERE=Path(__file__).resolve().parent
initial=json.loads((HERE/'PREREG.json').read_text())
for name,digest in initial['files'].items():
    assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,name
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
for folder in(HERE,HERE/'berry'):
    subprocess.run([sys.executable,str(folder/'headline.py')],cwd=folder,env=env,check=True)
