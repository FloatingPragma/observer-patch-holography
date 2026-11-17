"""Complete selected-spin modular and response geometry certificate."""
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
from . import native, pulse, regional

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
        'model/native_ordered_ground/witness.npz',
        'model/native_ordered_ground/README.md',
        'model/native_spin_clock_join/README.md',
        'model/native_spin_observer/README.md',
        'model/native_gravity_identification/PROOF.md',
        'model/native_regional_modular_clock/README.md',
        'Lean/ObserverPatchHolography/EinsteinBranch/Entropy.lean',
        'Lean/ObserverPatchHolography/EinsteinBranch/SmallBall.lean',
        'Lean/ObserverPatchHolography/EinsteinBranch/Composition.lean',
        'Lean/ObserverPatchHolography/EinsteinBranch/RegisterSurface.lean'))
    paths.update((HERE/'README.md', HERE/'PROOF.md'))
    return {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(paths)}


def claims():
    return {
        'same_complete_eight_spin_source_reference_and_H_B_energy_convention': True,
        'full_256_eigenvalue_ground_certificate_rechecked': True,
        'actual_regional_faithfulness_and_certified_logarithm': True,
        'complete_native_regional_modular_flow_executed': True,
        'complete_native_blank_preparation_and_source_current_pulse_executed': True,
        'complete_native_source_current_pointer_contrast_certified': True,
        'complete_native_16_setting_readout_executed': True,
        'positive_regional_detector_Fisher_metric': True,
        'isotropic_internal_bond_modular_boost_ansatz_excluded': True,
        'fixed_graph_area_universal_stationarity_refuted': True,
        'current_filter_exact_Ward_bound': True,
        'native_output_postselected_or_renormalized': False,
        'Gibbs_reference_substituted_for_source_ground': False,
        'infinite_ordered_state_marginal_identified': False,
        'physical_modular_boost_normalization_derived': False,
        'controlled_modular_or_current_pulse_is_autonomous_H_B_evolution': False,
        'nonzero_spatial_momentum_or_symmetric_stress_derived': False,
        'regional_Fisher_metric_is_Lorentzian_spacetime_metric': False,
        'responsive_area_Newton_constant_or_stationarity_derived': False,
        'continuum_Einstein_premises_discharged': False,
        'particle_rest_masses_or_Einstein_dynamics_derived': False}


def rational(value):
    if type(value) is F:
        if max(value.numerator.bit_length(), value.denominator.bit_length()) > 6000:
            return {'outward_lower': str(-regional.upper(-value)),
                    'outward_upper': str(regional.upper(value))}
        return str(value)
    if type(value) is dict:
        return {k: rational(v) for k, v in value.items()}
    if type(value) in (list, tuple):
        return [rational(v) for v in value]
    return value


def produce():
    pins = source_bindings()
    state = regional.state()
    value = {
        'schema': 'oph.spin_modular_geometry.v1', 'all_pass': True,
        'source_ground_certificate': state['full_source_certificate'],
        'source_reference': rational({k: state[k] for k in
            ('rho', 'tangent', 'ground_vector_error', 'reference_floor', 'actual_floor',
             'current_norm_upper')}),
        'rational_sylvester_principal_minors': [str(v) for v in state['principal_minors']],
        'certified_modular_logarithm': rational(regional.logarithm()),
        'same_source_response': rational(regional.source_response()),
        'finite_source_pulse_bound': rational(regional.finite_pulse()),
        'current_filter': rational(regional.current_filter()),
        'native_modular_flow': native.modular_flow(),
        'native_complete_detector': native.detector(),
        'native_blank_preparation_and_current_pulse': pulse.run(),
        'native_complete_current_readout': pulse.readout(),
        'claims': claims(),
        'finite_theorem_attachment': {
            'Entropy.FiniteFirstLawPremises_first_conjunct': True,
            'physical_2pi_bulk_plus_center_split': False,
            'center_is_geometric_edge_variation': False,
            'SmallBallPremises': False, 'ContinuumEinsteinPremises': False},
        'source_bindings': pins}
    if pins != source_bindings():
        raise ValueError('source changed during complete regional replay')
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
            raise ValueError(path+': complete record carrier changed')
        for j, (a, e) in enumerate(zip(actual, expected, strict=True)):
            compare(a, e, path+'.'+str(j))
    elif type(expected) is float:
        if not math.isfinite(actual) or actual < 0 or abs(actual-expected) > 2e-11:
            raise ValueError(path+': numeric diagnostic changed')
    elif actual != expected:
        raise ValueError(path+': exact value changed')


def instrument_rows(value):
    return ([value['native_modular_flow'], value['native_blank_preparation_and_current_pulse'],
             value['native_complete_detector']['selector']]
            + value['native_complete_detector']['settings']
            + value['native_complete_current_readout']['experiments'])


def validate(value):
    expected = produce()
    compare(value, expected)
    for actual, fresh in zip(instrument_rows(value), instrument_rows(expected), strict=True):
        a = np.asarray(actual['native_first_failure_weights'])
        b = np.asarray(fresh['native_first_failure_weights'])
        if np.linalg.norm(np.sqrt(a)-np.sqrt(b)) > 2e-11:
            raise ValueError('complete joint failure amplitude changed')
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
    print('PASS: selected-source regional modular flow and response geometry')


if __name__ == '__main__':
    main()

