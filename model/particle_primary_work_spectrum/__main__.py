"""Replay source modular fluctuations and critical-phase spectral bounds."""
import argparse
import json
from pathlib import Path
from . import model

OUTPUT=Path(__file__).with_name('results.json')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    action=parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--write',action='store_true');action.add_argument('--check',action='store_true')
    args=parser.parse_args();text=json.dumps(model.certificate(),sort_keys=True,indent=2,allow_nan=False)+'\n'
    if args.write:OUTPUT.write_text(text)
    elif not OUTPUT.exists() or OUTPUT.read_text()!=text:raise ValueError('work-spectrum receipt differs')
    print('PASS: actual modular/affinity and work records bounded; fixed-energy survival unproved')


if __name__=='__main__':main()
