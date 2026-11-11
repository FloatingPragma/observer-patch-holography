"""Compute the exact complete native scalar-cost response receipt."""
import argparse
import json
from pathlib import Path
from .response import certificate

parser = argparse.ArgumentParser(description=__doc__)
mode = parser.add_mutually_exclusive_group(required=True)
mode.add_argument('--write', action='store_true')
mode.add_argument('--check', action='store_true')
args = parser.parse_args()
value = certificate()
path = Path(__file__).with_name('results.json')
if args.write:
    path.write_text(json.dumps(value, sort_keys=True, indent=2)+'\n')
elif value != json.loads(path.read_text()):
    raise ValueError('exact scalar-cost response differs from receipt')
for family, row in value['channels'].items():
    result = row['complete_unconditioned_second_energy_transfer_moment']
    print(family, 'M2', result['value_interval'], 'dM2/dmu', result['derivative_interval'])
