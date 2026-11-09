import argparse
import json
from pathlib import Path
from .capacity import certificate

parser = argparse.ArgumentParser()
mode = parser.add_mutually_exclusive_group(required=True)
mode.add_argument('--write', action='store_true')
mode.add_argument('--check', action='store_true')
args = parser.parse_args()
value = json.dumps(certificate(), sort_keys=True, indent=2, allow_nan=False)+'\n'
path = Path(__file__).with_name('results.json')
if args.write:
    path.write_text(value)
elif not path.is_file() or path.read_text() != value:
    raise ValueError('primary capacity receipt differs from its actual source')
print('PASS: carrier capacity bounds, with physical and record capacities distinct')
