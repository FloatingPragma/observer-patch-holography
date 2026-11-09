import argparse
import json
from pathlib import Path
from .filter import certificate

parser=argparse.ArgumentParser();mode=parser.add_mutually_exclusive_group(required=True)
mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true')
args=parser.parse_args();value=json.dumps(certificate(),indent=2,sort_keys=True,allow_nan=False)+'\n'
path=Path(__file__).with_name('results.json')
if args.write:path.write_text(value)
elif not path.exists() or path.read_text()!=value:raise SystemExit('filtered probe receipt differs')
print('PASS: complete finite filter apparatus and original Gram scores; no primary survival claim')
