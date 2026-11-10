"""Replay an equal-work complete apparatus on the unchanged selected source."""
import argparse
import json
from pathlib import Path
from . import model


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    action=parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--write',action='store_true')
    action.add_argument('--check',action='store_true')
    args=parser.parse_args();path=Path(__file__).with_name('results.json')
    text=json.dumps(model.certificate(),sort_keys=True,indent=2,allow_nan=False)+'\n'
    if args.write:path.write_text(text)
    elif not path.exists() or path.read_text()!=text:raise ValueError('equal-work receipt differs')
    print('PASS: same source and work; distinct complete apparatus has elastic limit')


if __name__=='__main__':main()
