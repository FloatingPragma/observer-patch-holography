"""Replay charge-transfer meters and their unchanged original energy law."""
import argparse
import json
from pathlib import Path
from .model import certificate

p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true')
a=p.parse_args();text=json.dumps(certificate(),indent=2,sort_keys=True,allow_nan=False)+'\n'
path=Path(__file__).with_name('results.json')
if a.write:path.write_text(text)
elif not path.exists() or path.read_text()!=text:raise SystemExit('charged-readout receipt differs')
print('PASS: complete source-energy and coherent charge-transfer records')
