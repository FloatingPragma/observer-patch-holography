"""Reproduce the finite native-routing join and its controlled detector bound."""
import argparse
import json
from pathlib import Path
from .certificate import certificate


def main():
    path=Path(__file__).with_name('results.json')
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true')
    group.add_argument('--check',nargs='?',const=path,type=Path)
    args=parser.parse_args();text=json.dumps(certificate(),sort_keys=True,indent=2)+'\n'
    if args.write:path.write_text(text)
    elif not args.check.exists() or args.check.read_text()!=text:
        raise ValueError('stale native-routing join receipt')
    print('PASS finite shared-clock native interaction, Gauss code, complete records and exact free restriction')


if __name__=='__main__':main()
