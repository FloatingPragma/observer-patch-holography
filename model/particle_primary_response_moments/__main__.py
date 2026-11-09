"""Replay the original selected-state response bound."""
import argparse
import json
from pathlib import Path
from .model import certificate

parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True)
group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
args=parser.parse_args()
text=json.dumps(certificate(),indent=2,sort_keys=True,allow_nan=False)+'\n'
path=Path(__file__).with_name('results.json')
if args.write:path.write_text(text)
elif not path.is_file() or path.read_text()!=text:raise SystemExit('response-moment receipt differs')
print('PASS: original populated finite excitation band, complete records; no mass inference')
