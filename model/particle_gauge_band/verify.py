"""Generate or independently replay complete band and quantum-link receipts."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import tempfile

from . import gauge_band, quantum_link_completion, native_link_projector_witness

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RECEIPTS = ('gauge_band.json', 'quantum_link_completion.json', 'native_link_projector_witness.json')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def produce(directory):
    gauge_band.execute(directory)
    quantum_link_completion.execute(directory)
    native_link_projector_witness.execute(directory)


def compare(expected, actual, path='root'):
    """Exact discrete data, numerical tolerance only for floating diagnostics.

    Eigenvector phases and degenerate bases are not compared. Every replay
    validates its complete matrices, bounds, covariance and CAR internally.
    """
    if isinstance(expected, dict):
        if set(expected) != set(actual):
            raise AssertionError(f'{path}: field set differs')
        for key in expected:
            if key not in {'source_bindings', 'witness_sha256', 'operator_artifact_sha256'}:
                compare(expected[key], actual[key], f'{path}/{key}')
    elif isinstance(expected, list):
        if len(expected) != len(actual):
            raise AssertionError(f'{path}: list length differs')
        for i, (left, right) in enumerate(zip(expected, actual, strict=True)):
            compare(left, right, f'{path}/{i}')
    elif isinstance(expected, float):
        if not math.isfinite(expected) or not math.isfinite(actual) or abs(expected-actual) > 5e-8+2e-10*abs(expected):
            raise AssertionError(f'{path}: numerical invariant differs')
    elif type(expected) is not type(actual) or expected != actual:
        raise AssertionError(f'{path}: exact value differs')


def check():
    stored = {name: json.loads((HERE/name).read_text()) for name in RECEIPTS}
    for name, packet in stored.items():
        for path, sha in packet['source_bindings'].items():
            if digest(ROOT/path) != sha:
                raise AssertionError(f'{name}: source binding changed: {path}')
    if digest(HERE/'gauge_band_witnesses.npz') != stored['gauge_band.json']['witness_sha256']:
        raise AssertionError('band binary witness hash differs')
    if digest(HERE/'completed_link_operators.json') != stored['quantum_link_completion.json']['operator_artifact_sha256']:
        raise AssertionError('completed link operators hash differs')
    with tempfile.TemporaryDirectory(prefix='oph-gauge-band-') as scratch:
        replay = Path(scratch)
        produce(replay)
        for name, packet in stored.items():
            compare(packet, json.loads((replay/name).read_text()), name)
        if json.loads((HERE/'completed_link_operators.json').read_text()) != json.loads((replay/'completed_link_operators.json').read_text()):
            raise AssertionError('exact link operators differ from independent replay')
    print('Complete gauge-band, rotor-coordinate and native CAR replays passed.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.check:
        check()
    else:
        produce(HERE)


if __name__ == '__main__':
    main()
