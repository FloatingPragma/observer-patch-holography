import argparse
import json
from .packet import build, validate, OUTPUT

parser = argparse.ArgumentParser(description='Replay quantitative spin packet readout')
action = parser.add_mutually_exclusive_group(required=True)
action.add_argument('--write', action='store_true')
action.add_argument('--check', action='store_true')
args = parser.parse_args()
value = build()
if args.write:
    OUTPUT.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n')
else:
    from tools import strict_json
    validate(strict_json.load(OUTPUT), value)
print('PASS quantitative spin packet readout; no physical mass or new native execution claimed')
