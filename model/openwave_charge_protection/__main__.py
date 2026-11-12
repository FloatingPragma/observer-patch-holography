"""Recompute compact degree protection and its current-attachment test."""
import argparse
import json
from .protection import HERE,certificate,validate_pins

parser=argparse.ArgumentParser(description=__doc__)
group=parser.add_mutually_exclusive_group(required=True)
group.add_argument('--write',action='store_true')
group.add_argument('--check',action='store_true')
args=parser.parse_args();path=HERE/'results.json'
if args.check:
    saved=json.loads(path.read_text());validate_pins(saved)
result=certificate()
if args.write:path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
elif result!=saved:raise ValueError('fresh winding/current certificate differs')
print(result['verdict'])
