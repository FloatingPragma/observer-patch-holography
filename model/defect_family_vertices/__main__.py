"""Build or check the complete native family scalar-response receipt."""
import argparse
import json
from pathlib import Path
from .certificate import certificate

parser = argparse.ArgumentParser(description=__doc__)
mode = parser.add_mutually_exclusive_group(required=True)
mode.add_argument('--write', action='store_true')
mode.add_argument('--check', action='store_true')
args = parser.parse_args()
value = certificate()
path = Path(__file__).with_name('results.json')
if args.write:
    path.write_text(json.dumps(value, sort_keys=True, indent=2)+'\n')
elif value != json.loads(path.read_text()):
    raise ValueError('native family scalar response differs')
print('PASS: full native contacts, complete scalar records, charged pulses and neutral channels')
