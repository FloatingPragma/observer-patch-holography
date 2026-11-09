import argparse
import json
from pathlib import Path
from .model import certificate

p = argparse.ArgumentParser()
g = p.add_mutually_exclusive_group(required=True)
g.add_argument('--write', action='store_true')
g.add_argument('--check', action='store_true')
a = p.parse_args()
path = Path(__file__).with_name('results.json')
value = json.dumps(certificate(), indent=2, sort_keys=True, allow_nan=False)+'\n'
if a.write:
    path.write_text(value)
elif not path.is_file() or path.read_text() != value:
    raise ValueError('geometric cut receipt differs')
print('PASS: covariant regional cut, boundary entropy and complete records')
