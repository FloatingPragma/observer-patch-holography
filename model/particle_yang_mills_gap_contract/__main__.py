import argparse
import json
from pathlib import Path
from .transfer import certificate

parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True)
group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
args=parser.parse_args();value=json.dumps(certificate(),indent=2,sort_keys=True)+'\n'
path=Path(__file__).with_name('results.json')
if args.write:path.write_text(value)
elif not path.exists() or path.read_text()!=value:raise SystemExit('Euclidean gap transfer receipt differs')
print('PASS: exact finite Euclidean correlations and same-clock gap criterion')
