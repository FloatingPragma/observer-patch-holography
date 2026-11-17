"""Replay the same-source current comparison and its explicit repair controller.

One receipt separates exact algebra from compiled apparatus observations.
Checking owns a fresh source replay; hashes and self-consistent claims alone
are insufficient. Numerical instrument diagnostics are checked against
retained factors; fixed admission budgets are checked separately from the
platform-dependent observations.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from . import native_response_comparison as comparison
from . import native_response_reduction as reduction
from . import native_response_reference as reference
from . import native_response_biased as biased
from . import native_response_dynamics as dynamics
from . import native_response_cooling as cooling

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = 'oph.native_response_compatibility.v2'
ARRAY_OWNERS = ('reference', 'cooling', 'controller')


def _controller():
    from . import native_response_controller
    return native_response_controller


def source_bindings():
    paths = (set(comparison.source_bindings()) | set(reference.source_bindings())
             | set(reduction.source_bindings()) | set(biased.source_bindings())
             | set(dynamics.source_bindings()) | set(cooling.source_bindings())
             | set(_controller().source_bindings()))
    paths.update('model/common_source/' + name for name in (
        'native_response_compatibility.py', 'test_native_response_compatibility.py',
        'native_response_comparison.py', 'test_native_response_comparison.py',
        'native_response_reduction.py', 'test_native_response_reduction.py',
        'native_response_reference.py', 'test_native_response_reference.py',
        'native_response_biased.py', 'test_native_response_biased.py',
        'native_response_dynamics.py', 'test_native_response_dynamics.py',
        'native_response_cooling.py', 'test_native_response_cooling.py',
        'native_response_controller.py', 'test_native_response_controller.py'))
    paths.add('tools/strict_json.py')
    return {path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
            for path in sorted(paths)}


def execute():
    exact = comparison.exact_certificate()
    limits = reduction.report()
    ground = biased.report()
    feedback = dynamics.source_certificate(bias=cooling.BIAS)
    reference_report, reference_arrays = reference.execute()
    cooling_report, cooling_arrays = cooling.execute()
    controller_report, controller_arrays = _controller().execute(reference_arrays, cooling_arrays,
        reference_report=reference_report, cooling_report=cooling_report)
    arrays = _pack_arrays(reference=reference_arrays, cooling=cooling_arrays, controller=controller_arrays)
    return {'schema': SCHEMA, 'comparison': exact, 'reduction': limits,
            'biased': ground, 'dynamics': feedback, 'reference': reference_report,
            'cooling': cooling_report, 'controller': controller_report,
            'source_composition': {
                'controller_reference_arrays': 'reference__',
                'controller_clock_and_bath_arrays': 'cooling__',
                'independent_surrogate_clock_samples_used': False,
                'reference_and_cooling_source_executions_reused_by_controller': True}}, arrays


def _pack_arrays(**owners):
    if owners.keys() != set(ARRAY_OWNERS):
        raise ValueError('all three source array owners are required')
    result = {}
    for owner in ARRAY_OWNERS:
        values = owners[owner]
        if type(values) is not dict or not values:
            raise ValueError('a nonempty source array mapping is required: '+owner)
        for key, value in values.items():
            if type(key) is not str or not key or '__' in key:
                raise ValueError('an unambiguous local source array name is required')
            result[owner+'__'+key] = value
    return result


def split_arrays(arrays):
    """Validate ownership prefixes without coercing or copying original arrays."""
    if type(arrays) is not dict:
        raise ValueError('an original complete array dictionary is required')
    result = {owner: {} for owner in ARRAY_OWNERS}
    for key, value in arrays.items():
        if type(key) is not str:
            raise ValueError('original source array names must be strings')
        owner, separator, name = key.partition('__')
        if not separator or owner not in result or not name or '__' in name:
            raise ValueError('unknown or ambiguous source array owner: '+key)
        result[owner][name] = value
    if any(not values for values in result.values()):
        raise ValueError('a complete source array owner is missing')
    return result


def _exact(expected, actual, path='receipt'):
    """Exact original types, no Boolean/integer or rounded-number aliases."""
    if type(actual) is not type(expected):
        raise ValueError('response compatibility type differs: ' + path)
    if type(expected) is dict:
        if expected.keys() != actual.keys():
            raise ValueError('response compatibility keys differ: ' + path)
        for key in expected:
            _exact(expected[key], actual[key], path + '.' + key)
    elif type(expected) is list:
        if len(expected) != len(actual):
            raise ValueError('response compatibility length differs: ' + path)
        for i, (left, right) in enumerate(zip(expected, actual, strict=True)):
            _exact(left, right, path + '[' + str(i) + ']')
    elif expected != actual:
        raise ValueError('response compatibility value differs: ' + path)


def compare_execution(report, arrays, fresh, fresh_arrays):
    if type(report) is not dict or report.keys() != fresh.keys():
        raise ValueError('response compatibility complete report differs')
    for key in fresh:
        if key not in ARRAY_OWNERS:
            _exact(fresh[key], report[key], key)
    retained = split_arrays(arrays)
    replayed = split_arrays(fresh_arrays)
    reference.compare_execution(report['reference'], retained['reference'],
                                fresh['reference'], replayed['reference'])
    cooling.compare_execution(report['cooling'], retained['cooling'],
                              fresh['cooling'], replayed['cooling'])
    _controller().compare_execution(report['controller'], retained['controller'],
        fresh['controller'], replayed['controller'],
        reference_arrays=retained['reference'], cooling_arrays=retained['cooling'],
        fresh_reference_arrays=replayed['reference'], fresh_cooling_arrays=replayed['cooling'])
    # Historical primitive leakage diagnostics were removed by the owners.
    # The remaining resources and admission budgets are word-derived values,
    # so a broad array-replay tolerance must not forge their small entries.
    for owner in ARRAY_OWNERS:
        for section in ('resources', 'error_budget'):
            if section in fresh[owner]:
                _fixed_numbers(fresh[owner][section], report[owner][section], owner+'.'+section)


def _fixed_numbers(expected, actual, path):
    if type(actual) is not type(expected):
        raise ValueError('fixed source value type differs: '+path)
    if type(expected) is dict:
        if expected.keys() != actual.keys():
            raise ValueError('fixed source value keys differ: '+path)
        for key in expected:
            _fixed_numbers(expected[key], actual[key], path+'.'+key)
    elif type(expected) is list:
        if len(expected) != len(actual):
            raise ValueError('fixed source value length differs: '+path)
        for index, (a, b) in enumerate(zip(expected, actual, strict=True)):
            _fixed_numbers(a, b, path+'['+str(index)+']')
    elif type(expected) is float:
        if not math.isfinite(actual) or not math.isclose(expected, actual, rel_tol=8e-15, abs_tol=0.):
            raise ValueError('fixed source numerical value differs: '+path)
    elif expected != actual:
        raise ValueError('fixed source exact value differs: '+path)


def check(path):
    path = Path(path)
    pins = source_bindings()
    report = strict_json.load(path)
    if type(report) is not dict or report.get('source_bindings') != pins:
        raise ValueError('response compatibility source bindings differ')
    artifact = path.with_suffix('.npz')
    expected = {'filename': artifact.name,
                'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}
    _exact(expected, report.get('artifact'), 'artifact')
    with np.load(artifact, allow_pickle=False) as stored:
        if len(stored.files) != len(set(stored.files)):
            raise ValueError('duplicate response artifact keys')
        arrays = {key: stored[key] for key in stored.files}
    fresh, fresh_arrays = execute()
    compare_execution({key: value for key, value in report.items()
                       if key not in ('source_bindings', 'artifact')},
                      arrays, fresh, fresh_arrays)
    if source_bindings() != pins:
        raise ValueError('response sources changed during replay')


def write(path):
    path = Path(path)
    artifact = path.with_suffix('.npz')
    if path.exists() or artifact.exists():
        raise ValueError('fresh response receipt and artifact paths required')
    pins = source_bindings()
    report, arrays = execute()
    if source_bindings() != pins:
        raise ValueError('response sources changed during execution')
    # Serialize before publishing either artifact; reject NaN or invalid JSON.
    report = json.loads(json.dumps(report, allow_nan=False))
    with artifact.open('xb') as output:
        np.savez_compressed(output, **arrays)
    report['source_bindings'] = pins
    report['artifact'] = {'filename': artifact.name,
                          'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}
    with path.open('x', encoding='utf-8', newline='\n') as output:
        json.dump(report, output, sort_keys=True, indent=2, allow_nan=False)
        output.write('\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--output', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    if args.check:
        check(args.check)
        print('Native response comparison, biased reference, cooling and controller: PASS')
    else:
        write(args.output)


if __name__ == '__main__':
    main()
