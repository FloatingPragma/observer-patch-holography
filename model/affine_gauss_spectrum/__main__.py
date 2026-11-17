"""Write or independently replay the full charged-sector energy enclosures."""
import argparse
import json
from pathlib import Path
from .model import certificate, validate_pins
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true')
    mode.add_argument('--check',action='store_true')
    args=parser.parse_args();path=Path(__file__).with_name('results.json')
    if args.check:
        expected=json.loads(path.read_text());validate_pins(bound_pins(expected['input_sha256']))
    fresh=certificate()
    if args.write:path.write_text(json.dumps(fresh,indent=2)+'\n')
    elif expected!=fresh:raise ValueError('actual charged spectrum replay differs')
    print('PASS: actual charged spectral intervals, strict field dressing and complete preparation')


if __name__=='__main__':main()
