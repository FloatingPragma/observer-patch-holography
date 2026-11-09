"""Write or replay exact native charged-residue certificates."""
import argparse
import json
from pathlib import Path

from .residue import certificate


parser = argparse.ArgumentParser(description=__doc__)
group = parser.add_mutually_exclusive_group()
group.add_argument('--write', action='store_true')
group.add_argument('--check', action='store_true')
args = parser.parse_args()
value = certificate()
encoded = json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'
path = Path(__file__).with_name('certificate.json')
if args.write:
    path.write_text(encoded)
    print(path)
elif args.check:
    if not path.exists() or json.loads(path.read_text()) != value:
        raise SystemExit('native charged-residue certificate differs')
    print('native charged-residue certificate verified')
else:
    print(encoded, end='')
