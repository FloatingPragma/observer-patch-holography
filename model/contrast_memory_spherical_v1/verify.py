"""Recompute the A-contrast-memory-spherical-v1 model and compare it with the stored results.

    python -m model.contrast_memory_spherical_v1.verify                  # headline subset, a few minutes at most
    python -m model.contrast_memory_spherical_v1.verify --full           # everything, including 80 and 320 cells and the events
    python -m model.contrast_memory_spherical_v1.verify --output PATH    # write a fresh result file (no comparison)
    python -m model.contrast_memory_spherical_v1.verify --verify PATH    # compare with an explicit result file

The comparison uses tolerances: `|a - b| <= 1e-8 max(1, |a|)` for floating-point values, exact equality for
integers, strings and booleans, and wider bands for the few values sampled from floating-point probabilities.
Results must replay on different platforms and BLAS libraries, so bytes are never compared. The default mode
compares the keys it recomputes with the stored full result; `--full` also requires the same set of keys. The
exit status is non-zero on any mismatch, and the failing keys are named.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import platform
import sys
import time
from pathlib import Path

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')

HERE = Path(__file__).resolve().parent
RESULTS = HERE / 'RESULTS.json'
SCHEMA = 'A-contrast-memory-spherical-v1 results 1'
SKIP = {'environment', 'seconds', 'mode'}
# keys whose values come from samples drawn with floating-point probabilities: absolute tolerances
SAMPLED = {
    'particle_stability.single_readings.no_pair_runs': 5,
    'particle_stability.single_readings.runs_on_first_seam': 5,
    'particle_stability.single_readings.distinct_seams_read': 1,
}
SAMPLED_PREFIX = {'collar.reading.sampled.': 1e-4}


def compute(full=False):
    from . import carrier, catalogue, collar, events, extras, law, port_response, reader, stability, theorem
    from .common import Clock, clean
    import numpy, scipy, sympy
    clock = Clock()
    out = {'schema': SCHEMA, 'model': 'A-contrast-memory-spherical-v1', 'law': 'A-contrast-memory-spherical-v1'}
    out['local_law'] = law.exact_checks()
    clock.mark('local_law')
    out['carrier'] = carrier.run(levels=(0, 1, 2) if full else (0,))
    clock.mark('carrier')
    out['theorem'] = theorem.run(levels=(0, 1, 2) if full else (0,))
    clock.mark('theorem')
    st = stability.run(full)
    st['verdict'] = stability.verdict(st)
    out['particle_stability'] = st
    clock.mark('particle_stability')
    out['reader'] = reader.run(full)
    clock.mark('reader')
    out['events'] = events.run(full)
    clock.mark('events')
    out['catalogue'] = catalogue.run(full)
    clock.mark('catalogue')
    out['port_response'] = port_response.run(full)
    clock.mark('port_response')
    out['collar'] = collar.run(full)
    clock.mark('collar')
    out['extras'] = extras.run(full)
    clock.mark('extras')
    out['mode'] = 'full' if full else 'default'
    out['environment'] = {'python': platform.python_version(), 'numpy': numpy.__version__,
                          'scipy': scipy.__version__, 'sympy': sympy.__version__, 'platform': platform.platform()}
    out['seconds'] = clock.sections
    return clean(out)


def checks(result) -> list[str]:
    """Statements every run must satisfy, independent of the stored file."""
    from . import carrier
    bad = ['carrier.' + x for x in carrier.check(result['carrier'])]
    if not all(result['local_law']['checks'].values()):
        bad.append('local_law.checks')
    for key, ok in result['particle_stability']['verdict'].items():
        if not ok:
            bad.append('particle_stability.verdict.' + key)
    for key, entry in result['theorem'].items():
        if key.startswith('level') and not (entry['amax_below_6'] and entry['beta_below_sqrt10']):
            bad.append(f'theorem.{key}')
    if not result['reader']['finite_readout']['within_bounds']:
        bad.append('reader.finite_readout.within_bounds')
    for name, entry in result['events']['same_order_for_every_observer'].items():
        if not (entry['every_observer_order_equal'] and entry['sequences_differ'] and entry['largest_seam_record_difference'] == 0):
            bad.append('events.same_order_for_every_observer.' + name)
    return bad


def compare(expected, actual, path, failures, full):
    if path in SKIP or path.split('.')[-1] in SKIP:
        return
    if isinstance(actual, dict):
        if not isinstance(expected, dict):
            failures.append(f'{path}: type')
            return
        missing = [k for k in actual if k not in expected]
        if missing:
            failures.extend(f'{path}.{k}: absent in the stored result' for k in missing)
        if full:
            failures.extend(f'{path}.{k}: absent in this run' for k in expected if k not in actual and k not in SKIP)
        for k in actual:
            if k in expected:
                compare(expected[k], actual[k], f'{path}.{k}' if path else k, failures, full)
        return
    if isinstance(actual, list):
        if not isinstance(expected, list) or len(expected) != len(actual):
            failures.append(f'{path}: list length')
            return
        for i, (e, a) in enumerate(zip(expected, actual)):
            compare(e, a, f'{path}[{i}]', failures, full)
        return
    if any(isinstance(value, float) and not math.isfinite(value) for value in (expected, actual)):
        failures.append(f'{path}: non-finite value')
        return
    if path in SAMPLED:
        if type(expected) is not int or type(actual) is not int or abs(expected - actual) > SAMPLED[path]:
            failures.append(f'{path}: stored {expected}, computed {actual}')
        return
    for prefix, tol in SAMPLED_PREFIX.items():
        if path.startswith(prefix) and isinstance(actual, float):
            if type(expected) not in (int, float) or abs(expected - actual) > tol:
                failures.append(f'{path}: stored {expected}, computed {actual}')
            return
    if isinstance(actual, bool) or isinstance(expected, bool) or actual is None or isinstance(actual, str):
        if expected != actual or type(expected) is not type(actual):
            failures.append(f'{path}: stored {expected!r}, computed {actual!r}')
        return
    if isinstance(actual, int):
        if type(expected) is not int or expected != actual:
            failures.append(f'{path}: stored {expected}, computed {actual}')
        return
    if not isinstance(expected, (int, float)) or not math.isfinite(expected):
        failures.append(f'{path}: stored {expected!r}, computed {actual!r}')
        return
    if abs(expected - actual) > 1e-8 * max(1.0, abs(expected)):
        failures.append(f'{path}: stored {expected!r}, computed {actual!r}')


def main(argv=None):
    ap = argparse.ArgumentParser(description='Recompute the A-contrast-memory-spherical-v1 model and compare with stored results.')
    ap.add_argument('--full', action='store_true', help='recompute everything (long)')
    group = ap.add_mutually_exclusive_group()
    group.add_argument('--output', type=Path, help='write a fresh result file')
    group.add_argument('--verify', type=Path, help='compare with this result file')
    args = ap.parse_args(argv)
    t0 = time.monotonic()
    result = compute(args.full)
    failures = checks(result)
    reference = None if args.output else (args.verify or RESULTS)
    if reference is not None:
        stored = json.loads(Path(reference).read_text())
        if stored.get('schema') != SCHEMA:
            failures.append('schema')
        compare(stored, result, '', failures, args.full)
    if args.output:
        args.output.write_text(json.dumps(result, indent=1) + '\n')
    for f in failures:
        print('MISMATCH ' + f, file=sys.stderr)
    print(json.dumps({'verified': not failures and reference is not None, 'mode': result['mode'],
                      'compared_with': str(reference.name) if reference is not None else None,
                      'written': str(args.output) if args.output else None, 'failures': len(failures),
                      'seconds': round(time.monotonic() - t0, 1)}))
    if failures:
        raise SystemExit(1)
    return result


if __name__ == '__main__':
    main()
