"""Explicit receipt generation or independent read-only replay."""
import argparse
import json
from pathlib import Path
from .certificate import certificate,write,validate_pins,replay_differences

parser=argparse.ArgumentParser(description=__doc__)
mode=parser.add_mutually_exclusive_group(required=True)
mode.add_argument('--write',action='store_true',help='regenerate the complete receipt')
mode.add_argument('--check',action='store_true',help='independently replay without writing')
args=parser.parse_args()
if args.write:
    result=write()
else:
    expected=json.loads(Path(__file__).with_name('results.json').read_text())
    validate_pins(expected['input_sha256'])
    result=certificate()
    differences=replay_differences(expected,result)
    if differences:raise SystemExit('receipt replay differs: '+', '.join(differences[:8]))
print(result['schema'])
print('Exact full-domain EOM/current/spin identities and complete same-H residues verified.')
