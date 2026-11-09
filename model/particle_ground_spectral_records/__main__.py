"""Exact receipt and separately labelled floating spectral diagnostics."""
import argparse
import json
from pathlib import Path

from .records import certificate, spectral_diagnostics


parser = argparse.ArgumentParser(description=__doc__)
group = parser.add_mutually_exclusive_group()
group.add_argument('--write', action='store_true')
group.add_argument('--check', action='store_true')
group.add_argument('--diagnostics', action='store_true')
args = parser.parse_args()
value = spectral_diagnostics() if args.diagnostics else certificate(include_diagnostics=False)
encoded = json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'
path = Path(__file__).with_name('certificate.json')
if args.write:
    path.write_text(encoded)
    print(path)
elif args.check:
    if not path.exists() or json.loads(path.read_text()) != value:
        raise SystemExit('native spectral record certificate differs')
    print('native spectral record certificate verified')
else:
    print(encoded, end='')
