import argparse
import json
from pathlib import Path

from .model import build, render
from .steps import render as render_steps
from .catalogue import render as render_catalogue
from .paths import render as render_paths

p = argparse.ArgumentParser()
p.add_argument('--write', action='store_true')
p.add_argument('--check', action='store_true')
a = p.parse_args()
root = Path(__file__).parent
value = build()
outputs = {'results.json': json.dumps(value, indent=2, sort_keys=True) + '\n',
           'RESULTS.md': render(value),
           'CONFIRMATION_STEPS.md': render_steps(value),
           'CATALOGUE.md': render_catalogue(value),
           'MASS_PATHS.md': render_paths(value)}
for name, expected in outputs.items():
    path = root / name
    if a.write:
        path.write_text(expected)
    if a.check and path.read_text() != expected:
        raise SystemExit('candidate test surface differs: ' + name)
print('PASS: candidate identity tests retain their source and physical domains')
