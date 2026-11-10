import argparse
import json
from pathlib import Path
from . import model

p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args()
path=Path(__file__).with_name('results.json');expected=json.dumps(model.certificate(),indent=2,sort_keys=True,allow_nan=False)+'\n'
if a.write:path.write_text(expected)
elif path.read_text()!=expected:raise ValueError('energy-response receipt differs from actual source')
print('PASS: original complete source energy-current response; geometric attachment explicit')
