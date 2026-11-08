"""Replayed source energy and regional entropy response certificate."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from . import native, response

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).parent


def source_bindings():
    paths = set(HERE.glob('*.py'))
    queue = list(paths)
    legacy = ROOT/'evidence/observer_dynamics_20260925/oph-physics-sim'

    def include(path):
        if path.is_file() and path not in paths:
            paths.add(path)
            queue.append(path)

    include(legacy/'oph_exact/carrier.py')
    include(legacy/'oph_exact/__init__.py')
    while queue:
        path = queue.pop()
        root = legacy if path.is_relative_to(legacy) else ROOT
        package = '.'.join(path.relative_to(root).with_suffix('').parts[:-1])
        for node in ast.walk(ast.parse(path.read_text())):
            names = []
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                name = '.'*node.level+(node.module or '')
                base = importlib.util.resolve_name(name, package) if node.level else name
                names = [base]+[base+'.'+a.name for a in node.names]
            for name in names:
                root = legacy if name.split('.')[0] in ('oph_fpe', 'oph_exact') else ROOT
                segments = name.split('.')
                include(root.joinpath(*segments).with_suffix('.py'))
                for j in range(1, len(segments)+1):
                    include(root.joinpath(*segments[:j])/'__init__.py')
    paths.update(ROOT/p for p in (
        'code/source_selection_model/response.json',
        'evidence/native_geometric_source_20260925/spec.json',
        'code/m1_source_realization/receipt.json',
        'model/native_spectral_detector/contract.cpp',
        'model/native_spin_clock_join/README.md',
        'model/native_spin_observer/README.md',
        'model/native_gravity_identification/PROOF.md',
        'model/native_regional_modular_clock/README.md',
        'Lean/ObserverPatchHolography/EinsteinBranch/Entropy.lean',
        'Lean/ObserverPatchHolography/EinsteinBranch/SmallBall.lean',
        'Lean/ObserverPatchHolography/EinsteinBranch/Composition.lean'))
    paths.update((HERE/'README.md', HERE/'PROOF.md'))
    return {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(paths)}


def claims():
    return {
        'actual_source_preparation_pulse_and_complete_readout_executed': True,
        'same_selected_input_and_observation_clock_for_all_pulses': True,
        'general_finite_driven_energy_continuity': True,
        'current_source_derivative_and_work_term': True,
        'regional_modular_first_law_on_actual_source_state': True,
        'nonzero_entropy_response_realized': True,
        'uniform_faithfulness_floor_for_two_spin_control': True,
        'general_finite_retarded_response_and_error_bound': True,
        'native_success_postselected_or_renormalized': False,
        'arbitrary_infinite_ground_preparation_modulus': False,
        'uniform_floor_for_arbitrary_growing_regions': False,
        'source_weights_identified_as_geometric_lapse': False,
        'physical_modular_boost_or_symmetric_Lorentz_stress': False,
        'spatial_momentum_Ward_identity': False,
        'geometric_generalized_entropy_stationarity': False,
        'responsive_area_and_universal_Newton_scale': False,
        'Einstein_dynamics_or_relativistic_particle_masses': False}


def produce():
    pins = source_bindings()
    examples = response.response_bounds(F(3, 4), F(1), F(1), F(1, 16))
    native_runs = native.check()
    for run in native_runs:
        error = F(run['complete_isometry_error_bound'])
        run['decoded_actual_modular_bounds'] = {
            k: str(v) for k, v in response.modular_stability(F(1, 5), error).items()}
    value = {
        'schema': 'oph.spin_stress_response.v1', 'all_pass': True,
        'exact_driven_energy_checks': response.exact_currents(),
        'same_source_modular_response': response.exact_response(),
        'finite_bound_example': {k: str(v) for k, v in examples.items()},
        'relative_entropy_remainder_example': str(response.first_law_bound(F(1, 5), F(3, 64))),
        'complete_native_runs': native_runs, 'claims': claims(),
        'finite_theorem_attachment': {
            'module': 'OPH.EinsteinBranch.Entropy',
            'predicate': 'FiniteFirstLawPremises',
            'supplied_conjunct': 'entropyDifferential = modularPairing',
            'physical_2pi_bulk_plus_edge_split_supplied': False,
            'SmallBallPremises_supplied': False,
            'ContinuumEinsteinPremises_supplied': False},
        'source_bindings': pins}
    if pins != source_bindings():
        raise RuntimeError('source changed during response verification')
    return value


def compare(actual, expected, path='receipt'):
    if type(actual) is not type(expected):
        raise ValueError(path+': type changed')
    if type(expected) is dict:
        if actual.keys() != expected.keys():
            raise ValueError(path+': keys changed')
        for key in expected:
            compare(actual[key], expected[key], path+'.'+key)
    elif type(expected) is list:
        if len(actual) != len(expected):
            raise ValueError(path+': complete outcome carrier changed')
        for j, (a, e) in enumerate(zip(actual, expected, strict=True)):
            compare(a, e, path+'.'+str(j))
    elif type(expected) is float:
        if not math.isfinite(actual) or actual < 0 or abs(actual-expected) > 2e-11:
            raise ValueError(path+': numeric diagnostic changed')
    elif actual != expected:
        raise ValueError(path+': exact value changed')


def validate(value):
    expected = produce()
    compare(value, expected)
    for a, e in zip(value['complete_native_runs'], expected['complete_native_runs'], strict=True):
        # One global complement-amplitude tolerance, never one per small record.
        aw = np.asarray(a['native_first_failure_weights'])
        ew = np.asarray(e['native_first_failure_weights'])
        if np.linalg.norm(np.sqrt(aw)-np.sqrt(ew)) > 2e-11:
            raise ValueError('complete first-failure amplitude changed')
    return expected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--check', type=Path, nargs='?', const=HERE/'results.json')
    group.add_argument('--write', type=Path, nargs='?', const=HERE/'results.json')
    args = parser.parse_args()
    if args.write:
        args.write.write_text(json.dumps(produce(), indent=2, sort_keys=True, allow_nan=False)+'\n')
    else:
        validate(strict_json.load(args.check))
    print('PASS: same-source driven energy and regional modular response')


if __name__ == '__main__':
    main()
