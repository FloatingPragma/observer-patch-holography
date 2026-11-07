"""Replay exact controls for the native twelve-direction transition theorem."""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from . import transition as t


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_FILES = (
    'model/native_lorentz_transition/transition.py',
    'model/native_lorentz_transition/certificate.py',
    'model/native_lorentz_transition/PROOF.md',
    'model/common_source/null_tomography.py',
    'model/common_source/slow_geometry.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py',
    'evidence/native_geometric_source_20260925/spec.json',
    'Lean/Screen/PrimitivePortFrameQuotient.lean',
    'Lean/Screen/PortFrameGram.lean',
    'Lean/ObserverPatchHolography/EinsteinBranch/Tensor.lean',
    'Lean/Geometry/LorentzOverlapCocycle.lean',
    'Lean/Geometry/CelestialSoldering.lean',
)


def _encode(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {key: _encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_encode(item) for item in value]
    return value


def build():
    source = t.source_contract()
    examples = t.controls()
    classified = {name: t.classify(rows) for name, rows in examples.items()
                  if name != 'null_collapse'}
    for name in ('identity', 'boost', 'rotation', 'native_cyclic_rotation', 'dilation'):
        t.proper_future(examples[name])
    if classified['parity']['proper'] or not classified['parity']['future_preserving']:
        raise AssertionError('parity was conflated with a proper transition')
    if classified['time_reversal']['future_preserving']:
        raise AssertionError('time reversal was accepted as future preserving')
    collapse = examples['null_collapse']
    if any(value != t.q(0) for value in t.null_charges(collapse)):
        raise AssertionError('the singular countercontrol is not null')
    try:
        t.classify(collapse)
    except ValueError as exc:
        if 'invertible' not in str(exc):
            raise
    else:
        raise AssertionError('a singular null collapse passed the transition gate')
    first = t.scaled(Fraction(3, 2), examples['boost'])
    second = t.scaled(Fraction(2, 3), examples['rotation'])
    composed = t.compose(first, second)
    if composed['metric_factor'] != t.q(1):
        raise AssertionError('inverse conformal scales did not cancel')
    frame = t.raw_frame()
    permutation = []
    cyclic = examples['native_cyclic_rotation']
    for axis in frame:
        image = tuple(t._dot(row[1:], axis) for row in cyclic[1:])
        permutation.append(frame.index(image))
    if sorted(permutation) != list(range(12)):
        raise AssertionError('the native rotation did not permute all ports')
    perturbed = [list(row) for row in examples['boost']]
    perturbed[1][2] = t.q(Fraction(1, 1000))
    residual = t.metric_residual(perturbed)
    if residual['maximum_metric_residual'] == t.q(0):
        raise AssertionError('the perturbed map should have a nonzero metric defect')
    return _encode({
        'schema': 'oph.native_lorentz_transition.v1',
        'source_bindings': {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                            for name in SOURCE_FILES},
        'exact_scalar_encoding': 'rational pair [a,b] denotes a+b*sqrt(5)',
        'native_design_rank': source['exact_design_rank'],
        'native_kernel_dimension': source['exact_kernel_dimension'],
        'native_metric': t.metric(),
        'controls': classified,
        'noncommuting_composition_with_cancelled_scales': composed,
        'native_cyclic_port_permutation': permutation,
        'perturbed_boost_metric_residual': residual,
        'checks': {'native_kernel_replayed': True,
                   'twelve_null_tests_imply_conformal_identity': True,
                   'properness_and_future_orientation_separate': True,
                   'singular_null_collapse_rejected': True,
                   'conformal_factors_compose': True,
                   'metric_error_bound_exact': True},
        'scope': {'linear_transition_supplied': True,
                  'physical_null_lift_supplied': True,
                  'source_selected_clock_or_conformal_scale': False,
                  'source_selected_boost': False,
                  'continuous_spinor_group_isomorphism_formalized': False,
                  'physical_spacetime_dynamics_derived': False},
        'all_pass': True,
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', nargs='?', const=OUTPUT, type=Path)
    mode.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    args = parser.parse_args()
    packet = build()
    if args.write is not None:
        args.write.write_text(json.dumps(packet, indent=2, sort_keys=True) + '\n')
        print('native Lorentz transition certificate written')
    else:
        if json.loads(args.check.read_text()) != packet:
            raise SystemExit('native Lorentz transition certificate differs')
        print('native Lorentz transition certificate matches')


if __name__ == '__main__':
    main()
