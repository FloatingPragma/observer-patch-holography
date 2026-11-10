import argparse
import json
from pathlib import Path

from .model import build, render

p = argparse.ArgumentParser()
p.add_argument('--write', action='store_true')
p.add_argument('--check', action='store_true')
a = p.parse_args()
root = Path(__file__).parent
value = build()
outputs = {'results.json': json.dumps(value, indent=2, sort_keys=True) + '\n',
           'RESULTS.md': render(value)}
for name, expected in outputs.items():
    path = root / name
    if a.write:
        path.write_text(expected)
    if a.check and path.read_text() != expected:
        raise SystemExit('candidate test surface differs: ' + name)
print('PASS: candidate identity tests retain their source and physical domains')
