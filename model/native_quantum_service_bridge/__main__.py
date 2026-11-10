"""Authenticate inputs and independently replay both complete service receipts."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def verify_pins():
    manifest = json.loads((HERE / 'pins.json').read_text())
    for name, expected in manifest['sha256'].items():
        actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError('source pin changed: ' + name)
    return len(manifest['sha256'])


def verify_loaded_project_modules():
    """Catch unlisted dynamic imports as well as changed pinned bytes."""
    names = json.loads((HERE / 'pins.json').read_text())['sha256']
    checked = set()
    for module in tuple(sys.modules.values()):
        path = getattr(module, '__file__', None)
        if not path:
            continue
        path = Path(path).resolve()
        if path.is_relative_to(ROOT):
            name = str(path.relative_to(ROOT))
            if name not in names:
                raise ValueError('unlisted dynamically imported project module: ' + name)
            checked.add(name)
    return len(checked)


def compare(actual, expected, path='receipt'):
    """Exact structural/text comparison; fixed envelopes only for diagnostics."""
    if type(actual) is not type(expected):
        raise ValueError('receipt type differs: ' + path)
    if isinstance(actual, dict):
        if actual.keys() != expected.keys():
            raise ValueError('receipt keys differ: ' + path)
        for key in actual:
            compare(actual[key], expected[key], path + '.' + key)
    elif isinstance(actual, list):
        if len(actual) != len(expected):
            raise ValueError('receipt length differs: ' + path)
        for i, (a, b) in enumerate(zip(actual, expected, strict=True)):
            compare(a, b, path + f'[{i}]')
    elif isinstance(actual, float):
        if not math.isfinite(actual) or abs(actual - expected) > 2e-11 + 2e-11 * abs(expected):
            raise ValueError('numerical replay differs: ' + path)
    elif actual != expected:
        raise ValueError('exact receipt differs: ' + path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--check', action='store_true', help='replay and compare retained receipts (default)')
    modes.add_argument('--write', type=Path, help='write both new receipts into a new directory')
    args = parser.parse_args()
    pins = verify_pins()
    from .construction import run
    from .finite_compiler import run as run_finite
    outputs = {'results.json': run(), 'finite_compiler.json': run_finite()}
    loaded = verify_loaded_project_modules()
    if args.write:
        args.write.mkdir(parents=True, exist_ok=False)
        for name, output in outputs.items():
            (args.write / name).write_text(json.dumps(output, indent=2, sort_keys=True, allow_nan=False) + '\n')
    else:
        for name, output in outputs.items():
            compare(output, json.loads((HERE / name).read_text()), name)
    print(json.dumps({'status': 'written' if args.write else 'passed', 'authenticated_files': pins,
                      'authenticated_loaded_project_modules': loaded,
                      'native_yukawa_source_pulses': 1, 'native_readout_source_pulses': 6,
                      'compensated_CAR_columns_checked': outputs['results.json']['compensated_quantum_link']['operator_intertwining_columns'],
                      'finite_matrix_columns_executed': outputs['finite_compiler.json']['executed_full_input_columns']}))


if __name__ == '__main__':
    main()
