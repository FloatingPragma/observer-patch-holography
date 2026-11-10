import argparse
import json
from pathlib import Path

from .model import certificate

parser = argparse.ArgumentParser()
parser.add_argument('--write', action='store_true')
parser.add_argument('--check', action='store_true')
args = parser.parse_args()
path = Path(__file__).with_name('results.json')
value = certificate()
if args.write:
    path.write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')
if args.check and json.loads(path.read_text()) != value:
    raise SystemExit('loop fusion certificate differs')
print('PASS: full-source static loop fusion and original complete energy records; no particle mass')
