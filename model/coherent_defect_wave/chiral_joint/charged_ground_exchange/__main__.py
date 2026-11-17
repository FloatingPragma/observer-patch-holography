"""Write exact evidence or replay it independently without mutations."""
import argparse
import json
from pathlib import Path
from .certificate import certificate,write,validate_pins
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[4]))
    from tools.pin_scope import bound_pins

parser=argparse.ArgumentParser(description=__doc__)
mode=parser.add_mutually_exclusive_group(required=True)
mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true')
args=parser.parse_args()
if args.write:result=write()
else:
    expected=json.loads(Path(__file__).with_name('results.json').read_text())
    validate_pins(bound_pins(expected['input_sha256']));result=certificate()
    if expected!=json.loads(json.dumps(result)):raise SystemExit('exact charged response replay differs')
print(result['schema'])
print('Complete charged responses and the exact pair-denominator exchange verified.')
