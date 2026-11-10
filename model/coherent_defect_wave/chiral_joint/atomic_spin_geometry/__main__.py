import argparse
import json
from pathlib import Path
from .certificate import certificate

parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--write', action='store_true')
group.add_argument('--check', action='store_true')
args = parser.parse_args()
result = certificate()
path = Path(__file__).with_name('results.json')
if args.write:
    path.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
elif json.loads(path.read_text()) != result:
    raise SystemExit('complete interacting spin geometry receipt differs')
print('PASS: full cube spectrum, actual source momenta and controlled interacting ground')
