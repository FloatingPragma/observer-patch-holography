"""Write or replay full-native refined charged-energy windows."""
import argparse,json
from pathlib import Path
from .certificate import certificate,validate
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true');args=parser.parse_args()
    path=Path(__file__).with_name('results.json')
    if args.check:
        expected=json.loads(path.read_text());validate(bound_pins(expected['input_sha256']))
    fresh=certificate()
    if args.write:path.write_text(json.dumps(fresh,indent=2,sort_keys=True)+'\n')
    elif fresh!=expected:raise ValueError('fresh native resolution replay differs')
    print('PASS: exact native roots, nested full-H weights, neutral shift and explicit refined selection')

if __name__=='__main__':main()
