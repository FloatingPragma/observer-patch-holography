"""Reproduce the complete native independent-command comparison."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json
from ..native_program import receipt as owning
from . import causality


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'
SOURCE_FILES = tuple(causality.dataflow.SOURCE_FILES)+tuple(
    'model/native_record_causality/'+name for name in
    ('__init__.py', 'causality.py', 'packet.py', 'test_causality.py', 'README.md', 'PROOF.md'))


def source_bindings():
    # Bind executable sources and their retained receipts. Documentation notes
    # inherited by a legacy receipt are not evidence for this theorem.
    result = {name: digest for name, digest in owning.source_bindings().items()
              if not name.startswith('docs/')}
    result.update({name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                   for name in SOURCE_FILES})
    return dict(sorted(result.items()))


def produce():
    bindings = source_bindings()
    result = causality.certificate()
    result['source_bindings'] = bindings
    if source_bindings() != bindings:
        raise ArithmeticError('native comparison sources changed during replay')
    return result


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def verify(path=OUTPUT):
    if strict_json.load(Path(path)) != produce():
        raise ValueError('native independent-command receipt differs from exact replay')
    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', nargs='?', const=OUTPUT, type=Path)
    group.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    args = parser.parse_args()
    if args.write is not None:
        args.write.write_text(canonical(produce()))
        print('Native independent-command receipt written.')
    else:
        verify(args.check)
        print('Native complete-record comparison verified.')
