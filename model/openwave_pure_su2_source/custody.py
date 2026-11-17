"""Strict source custody and finite numerical replay for source packages."""
import hashlib
import json
import math
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]


def pin_package(here, extra=()):
    paths = {str(p.relative_to(ROOT)) for p in here.iterdir()
             if p.suffix in ('.py', '.md', '.tex')}
    paths.update(extra)
    return bound_pins({p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
            for p in sorted(import_closure(paths))})


def equivalent(a, b):
    if isinstance(b, dict) and isinstance(a, dict):
        # Recorded pin maps are compared without prose and test entries.
        a = {k: bound_pins(v) if k == 'source_sha256' else v for k, v in a.items()}
        b = {k: bound_pins(v) if k == 'source_sha256' else v for k, v in b.items()}
    if isinstance(b, dict):
        return isinstance(a, dict) and set(a) == set(b) and all(equivalent(a[k], v) for k, v in b.items())
    if isinstance(b, list):
        return isinstance(a, list) and len(a) == len(b) and all(equivalent(x, y) for x, y in zip(a, b, strict=True))
    if isinstance(b, float):
        return type(a) is float and math.isfinite(a) and math.isfinite(b) and abs(a-b) <= 2e-11*max(1., abs(b))
    return type(a) is type(b) and a == b


def execute(here, pins, run):
    import argparse
    p = argparse.ArgumentParser()
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write', action='store_true')
    g.add_argument('--check', action='store_true')
    args = p.parse_args()
    before = pins()
    value = run()
    assert pins() == before, 'source changed during execution'
    value['source_sha256'] = bound_pins(before)
    target = here/'results.json'
    if args.write:
        target.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n')
    else:
        assert equivalent(json.loads(target.read_text()), value), 'complete receipt replay differs'
    print(value['verdict'])
