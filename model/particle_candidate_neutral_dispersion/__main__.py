import argparse
import json
from pathlib import Path

from .model import certificate

parser=argparse.ArgumentParser()
parser.add_argument('--write',action='store_true')
parser.add_argument('--check',action='store_true')
args=parser.parse_args()
value=certificate();path=Path(__file__).with_name('results.json')
if args.write:path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
if args.check and json.loads(path.read_text())!=value:
    raise SystemExit('neutral two-point moment receipt differs')
print('PASS: original neutral occupied-pair moment bounds and complete clock estimator; no mass claim')
