"""Replay exact composite identity tests and the complete charge counter."""
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
    data=json.dumps(model.certificate(),sort_keys=True,indent=2,allow_nan=False)+'\n'
    if args.write:path.write_text(data)
    elif not path.exists() or path.read_text()!=data:raise ValueError('composite identity receipt differs')
    print('PASS: source-current candidates and complete finite transfer detector')


if __name__=='__main__':main()
