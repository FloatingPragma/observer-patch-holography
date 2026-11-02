import argparse
import json
from pathlib import Path
from . import packet

parser = argparse.ArgumentParser(description=__doc__)
mode = parser.add_mutually_exclusive_group(required=True)
mode.add_argument('--write', action='store_true')
mode.add_argument('--check', action='store_true')
args = parser.parse_args()
path = Path(__file__).with_name('results.json')
if args.write:
    path.write_text(json.dumps(packet.build(), sort_keys=True, indent=2)+'\n')
else:
    packet.validate(json.loads(path.read_text()))
print('ballistic free-boundary arithmetic and custody PASS')
