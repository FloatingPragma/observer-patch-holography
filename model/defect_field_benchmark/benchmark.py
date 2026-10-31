"""Reproducible invariant receipt for the same-action charged-particle bridge."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from ..common_source import calibrated_matter as matter
from .quasiparticles import (ABSOLUTE_OPERATOR_GUARD, MODE_COUNT, NAMBU_COUNT,
                            ORTHOGONALITY_GUARD, current_probe_data, moving_symbol,
                            quasiparticles, validate_chart)

ROOT = Path(__file__).resolve().parents[2]
RECEIPT = Path(__file__).with_name('results.json')
SOURCES = (
    'model/defect_field_benchmark/quasiparticles.py',
    'model/defect_field_benchmark/benchmark.py',
    'model/common_source/calibrated_matter.py',
    'model/common_source/calibrated_matter_parameters.json',
    'model/common_source/calibrated_matching.py',
    'model/common_source/calibrated_matching_results.json',
    'model/common_source/chiral_candidate/action.py',
    'model/common_source/chiral_candidate/kinetic.py',
    'model/observer_particles/results/representation_identification_operators.json',
    'model/observer_particles/results/unified_yukawa_operators.json',
    'model/defect_field_bridge/colored.py',
    'Lean/Geometry/InternalEnergyInertia.lean',
)


def charged_mass_comparison():
    """Compare diagonalized outputs to independently evaluated input masses."""
    chart = quasiparticles()
    channels = {'up': (4, 6), 'down': (2, 6), 'charged_lepton': (6, 2)}
    rows = []
    for channel, (charge6, multiplicity) in channels.items():
        expected = [float(value) for value in matter.mass_coordinates(channel)]
        for sign in (-1, 1):
            blocks = [b for b in chart['projectors'] if b['charge6'] == sign*charge6]
            if len(blocks) != 3 or any(b['multiplicity'] != multiplicity for b in blocks):
                raise ArithmeticError('charged eigenspace multiplicity mismatch')
            for family, (block, target) in enumerate(zip(blocks, expected, strict=True)):
                if abs(block['energy']-target) > ABSOLUTE_OPERATOR_GUARD:
                    raise ArithmeticError('independent running input/action energy comparison failed')
                rows.append({'channel': channel, 'family_order': family,
                             'electric_charge_sixths': sign*charge6,
                             'multiplicity': multiplicity,
                             'diagonalized_energy_GeV': block['energy'],
                             'independent_y_times_v_over_sqrt2_GeV': target})
    return rows


def mass_shell_certificate():
    """Check the same mass kernel and Weyl symbol without identifying native p."""
    rest = quasiparticles()['bdg']
    energies = quasiparticles()['energies']
    axes = np.eye(3)
    derivative = [moving_symbol(p)-rest for p in axes]
    identity = np.eye(NAMBU_COUNT)
    for i, left in enumerate(derivative):
        if np.linalg.norm(left@rest+rest@left) > ABSOLUTE_OPERATOR_GUARD:
            raise ArithmeticError('mass/kinetic anticommutator failed')
        for j, right in enumerate(derivative):
            target = 2*identity if i == j else 0*identity
            if np.linalg.norm(left@right+right@left-target) > ORTHOGONALITY_GUARD:
                raise ArithmeticError('same-action Weyl Clifford identity failed')
    direction = np.array([2., 3., 6.])/7
    for magnitude in (.002, .13, 1.7):
        actual = np.linalg.eigvalsh(moving_symbol(direction*magnitude))
        wanted = np.sort(np.r_[np.hypot(energies, magnitude), -np.hypot(energies, magnitude)])
        if np.max(abs(actual-wanted)) > ABSOLUTE_OPERATOR_GUARD:
            raise ArithmeticError('complete same-mass finite-momentum spectrum failed')
    return {'same_mass_kernel_and_all_180_directions': True,
            'kinetic_anticommutator_and_Clifford_identities_checked': True,
            'supplied_speed_c': 1,
            'energy_squared': 'E(p)^2=E(0)^2+c^2|p|^2',
            'charged_inertial_curvature': 'd^2 E/dp_i^2 at 0 = c^2/E(0); m=E(0)/c^2',
            'momentum_pair_convention': 'Nambu particle p, conjugate hole -p; lower block -h(-p).T',
            'native_port_transport_is_this_momentum': False,
            'finite_momentum_source_execution': False,
            'rest_quasiparticle_basis_or_zero_half_continues_to_finite_p': False}


def certificate():
    chart = quasiparticles()
    validate_chart(chart['positive_basis'], chart['energies'], chart['charge6'])
    currents = current_probe_data()
    masses = charged_mass_comparison()
    parameters = matter.parameters()
    first, second = chart['ramsey_indices']
    positives = np.flatnonzero(chart['charge6'] == 6)
    partner = int(positives[0])
    # Only invariant quantities are serialized. Neither eigenvector entries
    # nor arbitrary vectors inside a degenerate subspace are receipt data.
    return {
        'schema': 'oph.defect_field_benchmark.static_charged_quasiparticles.v1',
        'role': 'specified executable calibrated matter law; same-input recovery, not held-out particle prediction',
        'source_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCES},
        'domain': {'original_annihilators': MODE_COUNT, 'all_Nambu_directions': NAMBU_COUNT,
                   'positive_charged_quasiparticles': 84, 'neutral_zero_quasiparticles': 6,
                   'neutral_zero_particle_hole_directions': 12,
                   'zero_half_convention': list(chart['zero_annihilator_indices']),
                   'spatial_extension_ports': 12,
                   'spatial_extension_CAR_modes': 1080,
                   'colored_sector_is_original_Fock90_sector': False,
                   'static_occupation_number': 2,
                   'per_port_hard_core_colored_dimension': 66*90**2,
                   'same_port_CAR_collision_complement_dimension': 12*(90*89//2)},
        'canonical_action': {'particle_hole': 'C(u,v)=(conjugate(v),conjugate(u))',
                             'full_basis': '[W,CW]',
                             'normal_form': 'H=E_vac+sum_a E_a gamma_a^dagger gamma_a',
                             'vacuum_offset_GeV': chart['vacuum_energy'],
                             'electric_charge': 'Q=sum_a q_a gamma_a^dagger gamma_a; neutral vacuum',
                             'all_CAR_action_and_charge_columns_checked': True,
                             'degenerate_chart_is_physical_selection': False,
                             'operator_absolute_residual_guard': ABSOLUTE_OPERATOR_GUARD,
                             'CAR_charge_absolute_residual_guard': ORTHOGONALITY_GUARD},
        'calibration': {'scheme_scale': parameters['common_scale'],
                        'input_role': parameters['role'],
                        'all_nine_Yukawas_v_and_CKM_are_inputs': True,
                        'source': parameters['sources']['SMDR_forward_inverse_calibration'],
                        'new_independent_mass_predictions': 0,
                        'running_action_energies_are_pole_masses': False,
                        'external_matching_is_native_interacting_evolution': False,
                        'empirical_fit_uncertainties_propagated_here': False},
        'charged_energy_blocks': masses,
        'currents': {'all_12_full_180_by_180_Nambu_probes_retained_by_API': True,
                     'primitive_basis_anomalous_block_norms': [p['pair_block_norm'] for p in currents['probes']],
                     'exact_neutral_Higgs_orbit_rank': currents['exact_higgs_orbit_rank'],
                     'exact_unbroken_span_dimension': currents['exact_stabilizer_dimension'],
                     'unbroken_span': 'color su(3) plus electric u(1)',
                     'unbroken_source_basis_coefficients': currents['exact_stabilizer_coefficients'],
                     'broken_anomalous_map_nonzero_singular_values': currents['anomalous_singular_values'][:3].tolist(),
                     'broken_probe_pair_blocks_dropped': False,
                     'all_12_currents_act_inside_static_Nqp2': False},
        'same_object_repair_action': {
            'native_branch': 'K_h=(1/sqrt(2)) sum_x |f_h(x)><x| tensor (-SWAP)^p(x,h) (U_h tensor U_h)',
            'internal_phase_action': 'U_h=exp[-i tau(1+j_h/4) (diag(E)+E_vac/2)+i kappa j_h diag(q)]',
            'common_clock': 'one supplied tau in inverse GeV for every color, hbar=1',
            'scheduler_rule': 'the chosen positive factor 1+j_h/4 rescales the entire static action',
            'spatial_action': '12 copies of the internal cell give H_spatial=12 E_vac+sum_(x,a) E_a n_(x,a)',
            'chosen_common_reference_subtraction': '-11 E_vac times identity, giving one global E_vac in the compiled sector',
            'exchange_action': 'pi times the +1 eigenspace projector of SWAP gives exp(-i A)=-SWAP',
            'action_is_derived_from_scalar_attempt_or_pulse_counts': False,
            'common_global_phase_is_absolute_energy_measurement': False},
        'same_charge_Ramsey_chart': {'first': first, 'second': second, 'opposite_charge_partner': partner,
                                   'electric_charge_sixths': int(chart['charge6'][first]),
                                   'gap_GeV': float(chart['energies'][second]-chart['energies'][first]),
                                   'initial_total_charge': 0,
                                   'common_electric_phase_cancels': True,
                                   'electric_charge_inferred_by_Ramsey_marginal': False,
                                   'separate_charge_probe': 'companion neutral QND pointer gives p0=(1-sin(q))/2',
                                   'absolute_common_offset_observed': False},
        'same_mass_kinematic_constraint': mass_shell_certificate(),
        'remaining_physical_joins': [
            'Measured charged masses enter through a replayed common-scale empirical calibration, rather than arising from repair counts.',
            'The compiled rest dynamics and hard-core true exchange must be read through the companion complete native/field instrument receipt.',
            'Bare pairing does not preserve physical N=2; only the static Bogoliubov quasiparticle N=2 is used here.',
            'All six neutral modes are massless here; measured neutrino oscillations require the existing optional Weinberg law and its own canonical join.',
            'Broken electroweak and time-dependent Higgs probes can change quasiparticle number and require a larger process domain.',
            'Same-port collisions, other occupation sectors, local gauge link dynamics and interacting poles are not supplied by this hard-core static law.',
            'Native spatial transport is not yet identified with the momentum and speed in the checked Weyl symbol.'],
    }


def compare_receipt(actual, expected, path='receipt'):
    """Portable invariant comparison; hashes, counts and text stay exact."""
    if isinstance(actual, dict) and isinstance(expected, dict):
        if actual.keys() != expected.keys():
            raise ValueError(f'{path}: receipt keys differ')
        for key in actual:
            compare_receipt(actual[key], expected[key], path+'.'+key)
    elif isinstance(actual, (tuple, list)) and isinstance(expected, (tuple, list)):
        if len(actual) != len(expected):
            raise ValueError(f'{path}: receipt lengths differ')
        for i, (left, right) in enumerate(zip(actual, expected, strict=True)):
            compare_receipt(left, right, path+f'[{i}]')
    elif type(actual) is float and type(expected) in (float, int):
        if not np.isfinite(actual) or not np.isfinite(expected) or not np.isclose(actual, expected, rtol=2e-11, atol=2e-12):
            raise ValueError(f'{path}: numerical invariant differs')
    elif actual != expected or type(actual) is not type(expected):
        raise ValueError(f'{path}: exact invariant differs')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = certificate()
    if args.write:
        RECEIPT.write_text(json.dumps(result, indent=2)+'\n')
        print(f'wrote {RECEIPT}')
    elif args.check:
        compare_receipt(result, json.loads(RECEIPT.read_text()))
        print('static charged quasiparticle receipt OK')
    else:
        print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
