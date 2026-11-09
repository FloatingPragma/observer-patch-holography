import argparse
import json
from pathlib import Path
from .model import certificate

p=argparse.ArgumentParser();mode=p.add_mutually_exclusive_group(required=True)
mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true')
args=p.parse_args();path=Path(__file__).with_name('results.json')
value=json.dumps(certificate(),sort_keys=True,indent=2,allow_nan=False)+'\n'
if args.write:path.write_text(value)
elif not path.is_file() or path.read_text()!=value:raise ValueError('charged-string receipt differs')
print('PASS: original charged endpoint-pair energy response; species and rest mass unproved')
