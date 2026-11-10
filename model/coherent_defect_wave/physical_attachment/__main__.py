"""Replay spin, currents and exact chiral mass conditions."""
import argparse
import json
from pathlib import Path
from .certificate import certificate


def main():
    output=Path(__file__).with_name('results.json')
    parser=argparse.ArgumentParser(description=__doc__)
    choice=parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--write',action='store_true')
    choice.add_argument('--check',nargs='?',const=output,type=Path)
    args=parser.parse_args()
    result=json.dumps(certificate(),sort_keys=True,indent=2)+'\n'
    if args.write:
        output.write_text(result)
    elif not args.check.exists() or args.check.read_text()!=result:
        raise ValueError('stale spin/current/chiral mass receipt')
    print('PASS exact routing spin/current attachment and chiral mass constraints')


if __name__=='__main__':
    main()
