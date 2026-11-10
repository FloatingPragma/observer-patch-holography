"""Write or freshly verify the finite affine-Gauss receipt."""
import argparse
import json
from pathlib import Path
from .certificate import certificate,validate_pins,write


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true')
    mode.add_argument('--check',action='store_true')
    args=parser.parse_args()
    if args.write:value=write()
    else:
        expected=json.loads(Path(__file__).with_name('results.json').read_text())
        validate_pins(expected['input_sha256'])
        value=certificate()
        if expected!=value:raise SystemExit('fresh exact receipt differs')
    print(value['schema'])
    print('Finite source, Gauss constraints and complete observer records verified.')


if __name__=='__main__':main()
