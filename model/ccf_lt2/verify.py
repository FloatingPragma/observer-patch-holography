#!/usr/bin/env python3
"""Execute and strictly verify the official finite common construction."""
from __future__ import annotations

if not __debug__:
    raise RuntimeError("CCF-LT-2 verification requires enabled assertions; remove -O/PYTHONOPTIMIZE")
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
if __package__ in (None, ''):
    sys.path.insert(0, str(HERE.parents[1]))
    __package__ = 'model.ccf_lt2'


def scientific_inputs():
    paths = [p for p in HERE.glob('*.py') if not p.name.startswith('test_')]
    paths += [HERE/name for name in ('source_snapshot.json', 'calibration_inputs.json',
                                    'representation_snapshot.json')]
    inputs = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
    for name in ('SPEC.md','AUTHORITY.json'):
        inputs['../'+name] = hashlib.sha256((HERE.parent/name).read_bytes()).hexdigest()
    return inputs


def read_result(path):
    """Reject ambiguous keys and nonfinite data before result comparison."""
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError('Duplicate result key: '+key)
            result[key] = value
        return result
    def reject(value):
        raise ValueError('Nonfinite JSON value: '+value)
    data = json.loads(Path(path).read_text(), object_pairs_hook=unique, parse_constant=reject)
    def finite(value):
        if isinstance(value, float) and not math.isfinite(value):
            raise ValueError('Nonfinite JSON number')
        if isinstance(value, dict):
            for item in value.values():
                finite(item)
        elif isinstance(value, list):
            for item in value:
                finite(item)
    finite(data)
    return data


def compare(expected, actual, path='root'):
    """Require all exact data and all numerical values, without dropped gates."""
    if isinstance(expected, dict):
        assert isinstance(actual, dict) and expected.keys() == actual.keys(), path
        for key in expected:
            if key not in ('elapsed_seconds', 'environment'):
                compare(expected[key], actual[key], path+'.'+key)
    elif isinstance(expected, list):
        assert isinstance(actual, list) and len(expected) == len(actual), path
        for i, (a, b) in enumerate(zip(expected, actual)):
            compare(a, b, path+f'[{i}]')
    elif isinstance(expected, float):
        assert isinstance(actual, (float, int)) and math.isfinite(actual), path
        assert math.isfinite(expected) and abs(expected-actual) <= 1e-10, (path, expected, actual)
    else:
        assert type(expected) is type(actual) and expected == actual, (path, expected, actual)


def run(full=False):
    from . import checks
    import numpy, scipy, sympy
    before = scientific_inputs()
    result = checks.run(full)
    assert scientific_inputs() == before, 'Scientific inputs changed during execution'
    result['model'] = 'CCF-LT-2'
    result['specification'] = 'A-ccf2-local-capacity-v2 finite matter, field and response construction'
    result['scientific_inputs_sha256'] = before
    result['environment'] = {'python': sys.version.split()[0], 'numpy': numpy.__version__,
                             'scipy': scipy.__version__, 'sympy': sympy.__version__}
    return json.loads(json.dumps(result))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full', action='store_true', help='Execute all complete-source and per-vertex capacity checks')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--output', type=Path, help='Regenerate a current result file')
    group.add_argument('--verify', type=Path, help='Strictly compare an explicit result file')
    args = parser.parse_args(argv)
    result = run(args.full)
    reference = args.verify or (None if args.output else HERE/('RESULTS.json' if args.full else 'HEADLINE.json'))
    if reference is not None:
        compare(read_result(reference), result)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'all_passed': True, 'model': 'CCF-LT-2',
                      'full_source_included': args.full,
                      'strict_result_comparison': reference is not None,
                      'elapsed_seconds': result['elapsed_seconds']}, indent=2))
    return result


if __name__ == '__main__':
    main()
