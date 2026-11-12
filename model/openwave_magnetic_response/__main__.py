"""Reproduce the finite external-field magnetic response."""
import argparse
import hashlib
import json
from pathlib import Path
from .response import run

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def produce():
    value=run()
    paths=[*sorted(HERE.glob('*.py')),HERE/'PROOF.md',ROOT/'model/native_mass_channels/algebra.py']
    value['source_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
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
