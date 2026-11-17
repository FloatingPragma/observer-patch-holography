#!/usr/bin/env python3
"""Read-only reconstruction of all active native and common-input results."""
import json,os,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
start=time.perf_counter()
env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',PYTHONINTMAXSTRDIGITS='0',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
for script in ['run.py','check.py','join_common.py']:
    process=subprocess.run([sys.executable,'-B',str(HERE/script)],cwd=HERE,env=env,text=True,capture_output=True)
    if process.returncode:
        print(process.stdout);print(process.stderr,file=sys.stderr);raise SystemExit(process.returncode)
result=json.loads((HERE/'RECEIPT.json').read_text())
print(json.dumps(dict(all_gates_pass=result['all_gates_pass'],literal_native_calls=result['literal_native_calls'],
    complete_outcomes=result['complete_outcomes'],preparation_inconclusive_upper=result['preparation_inconclusive_upper'],
    controller_seed_independent_of_K=result['controller_seed_independent_of_K'],
    replay_changes_evidence=False,seconds=time.perf_counter()-start),sort_keys=True))
