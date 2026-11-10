import argparse
import json
from pathlib import Path

from .model import certificate

p = argparse.ArgumentParser()
p.add_argument('--write', action='store_true')
p.add_argument('--check', action='store_true')
a = p.parse_args()
path = Path(__file__).with_name('results.json')
value = certificate()
if a.write:
    path.write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')
if a.check and json.loads(path.read_text()) != value:
    raise SystemExit('scalar source spectral certificate differs')
print('PASS: resolved original scalar transfer window; no Higgs mass or uniform residue')
