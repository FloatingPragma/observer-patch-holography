import argparse,json
from .certificate import certificate,HERE
from model.openwave_larmor_dynamics.certificate import equivalent
p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--check',action='store_true');a=p.parse_args()
value=certificate();path=HERE/'results.json'
if a.write:path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
if a.check:assert equivalent(json.loads(path.read_text()),value),'receipt mismatch'
print(value['verdict'])
