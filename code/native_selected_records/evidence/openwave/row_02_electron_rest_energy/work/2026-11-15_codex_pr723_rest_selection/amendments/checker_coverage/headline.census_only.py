"""Reproduce deterministic records and receipt offline; no remote imports."""
import hashlib
import time
from experiment import HERE, run, canonical

start=time.monotonic()
receipt,records=run()
for name,obj in [('RECEIPT.json',receipt),('RECORDS.json',records)]:
    encoded=canonical(obj).encode()
    if (HERE/name).read_bytes()!=encoded:
        raise SystemExit('Replay differs: '+name)
print(canonical({'component_checks':receipt['all_component_checks'],
                 'physical_rest_energy_established':False,
                 'native_calls':receipt['direct_native_calls'],
                 'receipt_sha256':hashlib.sha256(canonical(receipt).encode()).hexdigest(),
                 'seconds':time.monotonic()-start}))
raise SystemExit(0 if receipt['all_component_checks'] else 1)
