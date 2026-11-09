import argparse
import json
from pathlib import Path
from model.particle_common_contract import contract
from . import model

parser=argparse.ArgumentParser()
mode=parser.add_mutually_exclusive_group(required=True)
mode.add_argument('--write',action='store_true')
mode.add_argument('--check',action='store_true')
args=parser.parse_args()
path=Path(__file__).with_name('results.json')
expected=json.dumps(model.certificate(),indent=2,sort_keys=True,allow_nan=False)+'\n'
if args.write:path.write_text(expected)
elif not path.is_file() or path.read_text()!=expected:
    raise ValueError('boundary correlation receipt differs from its source')
print('PASS: Gauss-compatible boundary correlation and polynomial error')
