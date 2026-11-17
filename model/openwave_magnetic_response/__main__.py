"""Reproduce the finite external-field magnetic response."""
import argparse
import hashlib
import json
from pathlib import Path
from .response import run
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def produce():
    value=run()
    paths=[*sorted(HERE.glob('*.py')),HERE/'PROOF.md',ROOT/'model/native_mass_channels/algebra.py']
    value['source_sha256']=bound_pins({str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    return value


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    mode=ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true')
    args=ap.parse_args();value=produce();path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    elif json.loads(path.read_text())!=value:raise SystemExit('source response differs from the stored result')
    print(value['verdict'])


if __name__=='__main__':main()
