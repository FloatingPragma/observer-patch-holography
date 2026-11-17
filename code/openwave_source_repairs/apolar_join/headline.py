#!/usr/bin/env python3
"""Compact replay of complete frozen source, observer and transfer stages."""
from pathlib import Path
import subprocess,sys,time
HERE=Path(__file__).resolve().parent
started=time.perf_counter()
for name in ('route.py','field.py','transfer.py'):
    subprocess.run([sys.executable,'-B',str(HERE/name)],cwd=HERE,check=True)
print('PASS_COMPACT_SOURCE_APOLAR_JOIN seconds=%.6f'%(time.perf_counter()-started))
