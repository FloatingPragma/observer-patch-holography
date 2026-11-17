"""Recompute the free native beat and its complete source records."""
import argparse
import json
import hashlib
from .clock import HERE,ROOT,certificate
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

parser=argparse.ArgumentParser(description=__doc__)
mode=parser.add_mutually_exclusive_group(required=True)
mode.add_argument('--write',action='store_true')
mode.add_argument('--check',action='store_true')
args=parser.parse_args();path=HERE/'results.json'
if args.check:
    saved=json.loads(path.read_text())
    for name,digest in bound_pins(saved['source_sha256']).items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:raise ValueError('changed source: '+name)
value=certificate()
if args.write:path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
elif value!=saved:raise ValueError('fresh native beat certificate differs')
print(value['verdict'])
