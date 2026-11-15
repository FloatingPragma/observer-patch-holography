"""Capture read-only producer inputs; never invoke producer run/freeze functions."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import shutil
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
seen = set()


def audit(event, args):
    if event == 'open' and isinstance(args[0], (str, bytes)):
        try:
            p = Path(args[0])
            if str(p).startswith(str(ROOT)) and not str(p).startswith(str(HERE)) and p.is_file():
                seen.add(p)
        except (OSError, ValueError):
            pass


sys.addaudithook(audit)
sys.path[:0] = [str(ROOT/'openwave-tuning/common_dynamics/interacting_gaussian'),
               str(ROOT/'openwave-tuning/next/CHAIN2/local12_limit_audit')]
import branches
import certified_bounds
import fock_transfer
import refinement_adapter as native

fock_transfer.parameters(1, 1, 8)
native.source_constants(48)
native.reset_repetition_plan(1, F(1, 100))
branches.field.geometry()
for name in ('SPECTRAL_FRAME.json', 'FOCK.json', 'MANIFEST.json'):
    (ROOT/'openwave-tuning/common_dynamics/quadratic_field'/name).read_bytes()
for p in ('openwave-tuning/common_dynamics/quantum_field/gaussian_bounds.py',
          'openwave-tuning/common_dynamics/quantum_field/PROOF.md',
          'reverse-engineering-reality/code/source_selection_model/response.json',
          'reverse-engineering-reality/code/m1_source_realization/receipt.json'):
    (ROOT/p).read_bytes()
for module in list(sys.modules.values()):
    path = getattr(module, '__file__', None)
    if path:
        p = Path(path).resolve()
        if p.is_relative_to(ROOT) and not p.is_relative_to(HERE):
            seen.add(p)
seen = {Path(importlib.util.source_from_cache(str(p))) if p.suffix == '.pyc' else p for p in seen}
pins = {}
for p in sorted(seen):
    data = p.read_bytes()
    rel = str(p.relative_to(ROOT))
    target = HERE/'inputs'/rel
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(data)
    pins[rel] = hashlib.sha256(data).hexdigest()
assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest() == h for p,h in pins.items())
(HERE/'INPUTS.json').write_text(json.dumps({'schema':'oph.native_field.inputs.v1',
    'files_sha256':pins, 'capture':'read-only imports and parameter evaluation; no producer execution or mutation'}, indent=2)+'\n')
print(json.dumps({'frozen_files':len(pins), 'quadratic_manifest':pins['openwave-tuning/common_dynamics/quadratic_field/MANIFEST.json']}))
