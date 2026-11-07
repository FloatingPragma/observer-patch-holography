"""Reproduce the native canonical restriction and large-population receipt."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json
from ..native_program import receipt as program
from . import restriction, product


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'
SOURCE_FILES = tuple(restriction.owning.SOURCE_PATHS)+(
    'model/completion/program_class_v2.json',
    'model/local_repair_family/record_clock_campaign.json',
    'Lean/ObserverPatchHolography/Locality/NoSignalling.lean')+tuple(
    'model/native_capacity_restriction/'+name for name in
    ('__init__.py', 'restriction.py', 'packet.py', 'test_restriction.py',
     'product.py', 'test_product.py', 'PRODUCT.md', 'README.md', 'PROOF.md'))


def source_bindings():
    result = {name: digest for name, digest in program.source_bindings().items()
              if not name.startswith('docs/')}
    result.update({name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                   for name in SOURCE_FILES if not name.startswith('docs/')})
    return dict(sorted(result.items()))


def produce():
    bindings = source_bindings()
    result = restriction.certificate()
    result['independent_product_dynamic_restriction'] = product.certificate()
    result['scope']['product_theorem_is_separate_from_fair_canonical_family'] = True
    result['source_bindings'] = bindings
    if source_bindings() != bindings:
        raise ArithmeticError('canonical restriction sources changed during replay')
    return result


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def verify(path=OUTPUT):
    if strict_json.load(Path(path)) != produce():
        raise ValueError('native capacity restriction receipt differs from exact replay')
    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--write', action='store_true')
    actions.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(produce()))
        print('Native canonical restriction receipt written.')
    else:
        verify()
        print('Native canonical restriction receipt verified.')
