"""Reproduce complete prepared charged-line observer readouts."""
import argparse
import json
from .readouts import HERE, certificate, validate_pins

parser=argparse.ArgumentParser(description=__doc__)
mode=parser.add_mutually_exclusive_group(required=True)
mode.add_argument('--write',action='store_true')
mode.add_argument('--check',action='store_true')
args=parser.parse_args()
path=HERE/'results.json'
if args.check:
    saved=json.loads(path.read_text())
    validate_pins(saved)
result=certificate()
if args.write:
    path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
elif result!=saved:
    raise ValueError('fresh prepared-line readout certificate differs')
print(result['verdict'])
