"""Replay the actual-source energy forms and mixed-word history bounds."""
import argparse
import json
from pathlib import Path
from . import model

OUTPUT=Path(__file__).with_name('results.json')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    flags=parser.add_mutually_exclusive_group(required=True)
    flags.add_argument('--write',action='store_true');flags.add_argument('--check',action='store_true')
    args=parser.parse_args();value=json.dumps(model.certificate(),sort_keys=True,indent=2,allow_nan=False)+'\n'
    if args.write:OUTPUT.write_text(value)
    elif not OUTPUT.exists() or OUTPUT.read_text()!=value:raise ValueError('stale primary energy-form receipt')
    print('PASS: common primary Weyl/filtered histories; nonelastic survival unproved')


if __name__=='__main__':main()
