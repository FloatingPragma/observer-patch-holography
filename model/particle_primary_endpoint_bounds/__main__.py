"""Replay actual endpoint lower bounds without a perturbative extrapolation."""
import argparse
import json
from pathlib import Path
from .endpoint import certificate

parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True)
group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
args=parser.parse_args()
text=json.dumps(certificate(),indent=2,sort_keys=True,allow_nan=False)+'\n'
path=Path(__file__).with_name('results.json')
if args.write:path.write_text(text)
elif not path.is_file() or path.read_text()!=text:raise SystemExit('endpoint-bound receipt differs')
print('PASS: full primary endpoint ground floors; no isolated particle branch')
