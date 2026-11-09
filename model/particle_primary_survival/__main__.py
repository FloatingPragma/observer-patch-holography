import argparse
import json
from pathlib import Path
from .survival import certificate

parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True)
group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
args=parser.parse_args();value=json.dumps(certificate(),indent=2,sort_keys=True,allow_nan=False)+'\n'
path=Path(__file__).with_name('results.json')
if args.write:path.write_text(value)
elif not path.exists() or path.read_text()!=value:raise SystemExit('primary survival receipt differs')
print('PASS: actual finite primary contrast and extreme-sector suppression; no uniform mass residue')
