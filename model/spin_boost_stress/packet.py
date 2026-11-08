"""Same-source weighted-energy stress and finite boost response certificate."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np

from tools import strict_json
from model.spin_modular_geometry.packet import compare, rational
from . import algebra, response, native

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
        'model/protected_defect_catalog/PROOF.md',
        'model/spin_modular_geometry/PROOF.md',
        'model/source_energy_effort/PROOF.md',
        'model/native_spin_clock_join/energy.py',
        'Lean/QFT/CoarseEnergyWard.lean',
        'Lean/QFT/LatticeMomentumBalance.lean',
        'Lean/QFT/PeriodicMomentum.lean',
        'Lean/Dynamics/PrivateInner.lean',
        'Lean/Dynamics/StoneConverse.lean',
        'Lean/Dynamics/WardLimitManifest.lean',
        'paper/tex_fragments/OPERATIONAL_CAUSAL_SELECTION.tex',
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
    return {'actual_selected_source_energy_polarization_and_flux': True,
            'same_energy_flux_operator_under_H_B_and_H_J': True,
            'same_charge_preserved_by_source_candidate': True,
            'three_local_forced_stress_columns_constructed': True,
            'finite_coarse_Ward_premises_instantiated': True,
            'actual_finite_charge_preserving_translation_logs_constructed': True,
            'exact_subcell_energy_dipole_refinement_rule': True,
            'created_channel_energy_flux_acceleration_certified': True,
            'complete_native_weighted_energy_flow_executed_from_blank': True,
            'selected_ground_noninvariance_certified': True,
            'finite_boost_commutator_and_history_error_bounds': True,
            'candidate_speed_choice_removes_ground_defect': False,
            'source_energy_flux_is_conserved_spatial_momentum': False,
            'local_forced_stress_is_symmetric_Lorentz_stress': False,
            'subcell_energy_dipoles_may_be_discarded': False,
            'small_rapidity_steps_remove_second_commutator_defect': False,
            'operational_spacetime_boost_covariance_derived': False,
            'continuum_Ward_manifest_inhabited': False,
            'interacting_relativistic_dispersion_or_particle_mass_derived': False,
            'native_success_postselected_or_renormalized': False}


def produce():
    pins = source_bindings()
    controls = algebra.source_controls()
    value = {'schema': 'oph.spin_boost_stress.v1', 'all_pass': True,
             'source_controls': rational(controls),
             'coarse_energy_and_dipole_transport': algebra.coarse_transport(),
             'source_forced_stress': rational(algebra.forced_stress()),
             'actual_finite_translation_readouts': algebra.finite_translations(),
             'candidate_speed_controls': rational([algebra.speed_normalization(c)
                                                    for c in (F(1, 2), F(1), F(2))]),
             'complete_history_example': rational(response.history_bound(
                 F(1, 32), controls[0]['boost_defect_norm_upper'], [F(1, 4), F(1, 3)], [F(1, 8)])),
             'native_source_candidate': native.run(), 'claims': claims(),
             'reused_theorems': {
                 'native_spin_clock_join.energy.weighted_identity': 'actual eight-source matrices checked',
                 'QFT.CoarseEnergyWard.coarse_continuity': 'matrix additive-group premises checked',
                 'QFT.LatticeMomentumBalance.coarse_forced_balance': 'all real matrix pairings',
                 'QFT.PeriodicMomentum.polynomial_conserved': 'actual even source shifts checked',
                 'Dynamics.PrivateInner.hasDerivAt_realInnerExponentialFlow': 'actual Hermitian B source',
                 'Dynamics.StoneConverse': 'fixed source clock and ground scalar retained',
                 'Dynamics.WardLimitManifest': 'source refinement/chart/limit premises not supplied',
                 'operational_causal_selection': 'one-boost complete-history covariance not supplied'},
             'source_bindings': pins}
    if pins != source_bindings():
        raise ValueError('source changed during boost/stress certificate')
    return value


def validate(value):
    expected = produce()
    compare(value, expected)
    a, b = (np.asarray(row['native_source_candidate']['native_first_failure_weights'])
            for row in (value, expected))
    if np.linalg.norm(np.sqrt(a)-np.sqrt(b)) > 2e-11:
        raise ValueError('global complete failure amplitude changed')
    return expected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', nargs='?', type=Path, const=HERE/'results.json')
    group.add_argument('--check', nargs='?', type=Path, const=HERE/'results.json')
    args = parser.parse_args()
    if args.write:
        args.write.write_text(json.dumps(produce(), indent=2, sort_keys=True, allow_nan=False)+'\n')
    else:
        validate(strict_json.load(args.check))
    print('PASS: source weighted energy, forced stress and complete boost response')


if __name__ == '__main__':
    main()

