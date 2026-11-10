import argparse
import json
from pathlib import Path
from .certificate import certificate

parser = argparse.ArgumentParser()
mode = parser.add_mutually_exclusive_group(required=True)
mode.add_argument('--write',action='store_true')
mode.add_argument('--check',action='store_true')
args = parser.parse_args()
value = certificate()
path = Path(__file__).with_name('results.json')
if args.write:
    path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
elif json.loads(path.read_text()) != value:
    raise SystemExit('centre/current correspondence differs from its complete recomputation')
print('PASS: full exterior weights, native fifteen currents, centre kernel and parity witnesses')
