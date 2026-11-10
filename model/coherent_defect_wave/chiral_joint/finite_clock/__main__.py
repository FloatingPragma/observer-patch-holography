"""Write or independently replay the exact finite-clock receipt."""
import argparse
import json
from pathlib import Path
from ..receipt_replay import check_receipt
from .certificate import certificate

parser = argparse.ArgumentParser(description=__doc__)
mode = parser.add_mutually_exclusive_group(required=True)
mode.add_argument('--write', action='store_true')
mode.add_argument('--check', action='store_true')
args = parser.parse_args()
result = certificate()
rendered = json.dumps(result, indent=2, sort_keys=True)+'\n'
path = Path(__file__).with_name('results.json')
if args.write:
    path.write_text(rendered)
else:
    check_receipt(path, rendered, 'finite rational clock')
print('PASS: exact finite unitary, complete ring records and controlled neutral arrival')
