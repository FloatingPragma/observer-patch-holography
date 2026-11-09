"""Replay the actual selected source's critical-profile work bounds."""
import argparse
import json
from pathlib import Path
from . import model

OUTPUT=Path(__file__).with_name('results.json')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    flags=parser.add_mutually_exclusive_group(required=True)
    flags.add_argument('--write',action='store_true')
    flags.add_argument('--check',action='store_true')
    args=parser.parse_args()
    text=json.dumps(model.certificate(),sort_keys=True,indent=2,allow_nan=False)+'\n'
    if args.write:OUTPUT.write_text(text)
    elif not OUTPUT.exists() or OUTPUT.read_text()!=text:
        raise ValueError('critical-profile source receipt differs')
    print('PASS: selected critical-phase work has controlled spatial refinement; no mass inferred')


if __name__=='__main__':main()
