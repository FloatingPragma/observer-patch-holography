"""Check complete native transport, shared-link response and scaling boundaries."""
import argparse
import json
from pathlib import Path
from .certificate import certificate


def main():
    output=Path(__file__).with_name('results.json')
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true')
    group.add_argument('--check',nargs='?',const=output,type=Path)
    args=parser.parse_args();text=json.dumps(certificate(),sort_keys=True,indent=2)+'\n'
    if args.write:output.write_text(text)
    elif not args.check.exists() or args.check.read_text()!=text:
        raise ValueError('stale propagating-native join receipt')
    print('PASS native CAR mobility, common-clock Dirac controls, complete interaction bounds and scaling tests')


if __name__=='__main__':main()
