"""Reproduce deterministic records and receipt offline; no remote imports."""
import hashlib
import time
import json
from experiment import HERE, run, canonical
from trajectories import run as trajectory_run

start=time.monotonic()
receipt,records=run()
for name,obj in [('RECEIPT.json',receipt),('RECORDS.json',records)]:
    encoded=canonical(obj).encode()
    if (HERE/name).read_bytes()!=encoded:
        raise SystemExit('Replay differs: '+name)
trajectory_receipt,trajectory_records=trajectory_run()
for name,obj in [('TRAJECTORY_RECEIPT.json',trajectory_receipt),('TRAJECTORIES.json',trajectory_records)]:
    if (HERE/name).read_text()!=json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n':
        raise SystemExit('Replay differs: '+name)
print(canonical({'component_checks':receipt['all_component_checks'],
                 'trajectory_checks':all(trajectory_receipt['gates'].values()),
                 'physical_rest_energy_established':False,
                 'census_native_calls':receipt['direct_native_calls'],
                 'trajectory_native_calls':trajectory_receipt['native_event_records'],
                 'receipt_sha256':hashlib.sha256(canonical(receipt).encode()).hexdigest(),
                 'seconds':time.monotonic()-start}))
raise SystemExit(0 if receipt['all_component_checks'] and all(trajectory_receipt['gates'].values()) else 1)
