#!/usr/bin/env python3
"""Read-only offline headline; Python3 and NumPy, no network or producer imports."""
from pathlib import Path
import hashlib,os,subprocess,sys,time
HERE=Path(__file__).resolve().parent
def main():
    started=time.monotonic();manifest=HERE/'MANIFEST.sha256'
    if manifest.exists():
        for line in manifest.read_text().splitlines():
            expected,name=line.split('  ',1)
            assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==expected,name
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    result=subprocess.run([sys.executable,str(HERE/'integrate.py')],cwd=HERE,env=env)
    print('Read-only common-source headline seconds:',time.monotonic()-started)
    raise SystemExit(result.returncode)
if __name__=='__main__':main()
