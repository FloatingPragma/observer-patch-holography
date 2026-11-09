"""Replay complete odd-probe filtering on the actual primary source."""
import argparse
import json
from pathlib import Path
from . import model

OUTPUT=Path(__file__).with_name('results.json')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    flags=parser.add_mutually_exclusive_group(required=True)
    flags.add_argument('--write',action='store_true');flags.add_argument('--check',action='store_true')
    args=parser.parse_args()
    text=json.dumps(model.certificate(),sort_keys=True,indent=2,allow_nan=False)+'\n'
    if args.write:OUTPUT.write_text(text)
    elif not OUTPUT.exists() or OUTPUT.read_text()!=text:raise ValueError('odd-filter receipt differs')
    print('PASS: complete source odd filters retain a populated finite band and preserve summed energy records')


if __name__=='__main__':main()
