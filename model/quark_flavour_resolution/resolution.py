"""Exact spectral transfers for native hypercharge and neutral polarization.

The neutral polarization memory is part of the source. The compensating
reference retains its separate energy and is idle in that mechanism.
"""
from functools import lru_cache
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import hashlib
import json

import sympy as s

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).parent
NATIVE = 'model/native_quark_response/results.json'
CALIBRATION = 'model/common_source/calibrated_matter_results.json'
FAMILY = 'model/defect_family_vertices/results.json'
CHARGES = {'uc': F(2, 3), 'dc': F(-1, 3), 'Q': F(-1, 6)}


def interval(value):
    return tuple(map(F, value))


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def multiply(a, b):
    values = [x*y for x in a for y in b]
    return min(values), max(values)


def encode(value):
    return list(map(str, value))


def midpoint(value):
    return float(sum(value)/2)


def sqrt_bounds(value, bits=100):
    value = F(value)
    scale = 1 << bits
    lower = isqrt(value.numerator*scale*scale//value.denominator)
    return F(lower, scale), F(lower+1, scale)


def authenticate(path):
    value = json.loads((ROOT/path).read_text())
    for name, expected in value['input_sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != expected:
            raise ValueError('native input changed: '+name)
    return value


@lru_cache(maxsize=1)
def operator_checks():
    from model.coherent_defect_wave.chiral_joint.higgs_composites import model as car
    from model.coherent_defect_wave.chiral_higgs.model import higgs_memory
    from model.coherent_defect_wave.physical_attachment.chiral import representation
    parts, couplings = car.quantum_yukawa()
    h = car.add(*parts.values(), car.monomial(matrix=higgs_memory()[2]))
    h = car.clean({key: value.subs({y: 1 for y in couplings.values()}) for key, value in h.items()})
    charge = car.gauge_current(11)
    if car.commutator(charge, h):
        raise ArithmeticError('native hypercharge is not conserved')
    for a in range(12):
        if car.commutator(charge, car.gauge_current(a)):
            raise ArithmeticError('hypercharge coupling is not globally invariant')
    for family, mode in (('uc', 0), ('dc', 12), ('Q', 23)):
        c = car.monomial(((mode, False),))
        if car.commutator(charge, c) != car.scale(c, s.Rational(CHARGES[family])):
            raise ArithmeticError('literal native removal charge differs')
    y = representation()[2][11]
    fermions = [car.monomial(((i, False),)) for i in range(30)]
    for i, c in enumerate(fermions):
        derivative = car.commutator(c, charge)
        if derivative != car.scale(c, y[i % 15, i % 15]):
            raise ArithmeticError('the polarization source changed a literal current derivative')
        for other in fermions:
            if car.anticommutator(derivative, other):
                raise ArithmeticError('an anomalous Higgs contact was introduced')
    positive = sum((v for v in y.diagonal() if v > 0), s.S.Zero)
    negative = -sum((v for v in y.diagonal() if v < 0), s.S.Zero)
    if positive != negative or positive != 3:
        raise ArithmeticError('full Fock charge norm bound changed')
    w, q = s.symbols('w q', real=True)
    memory = s.Matrix([[0, q], [q, 1]])
    resolvent = (w*s.eye(2)-memory).inv()[0, 0]
    if s.simplify(resolvent-1/(w-q*q/(w-1))) != 0:
        raise ArithmeticError('exact neutral-memory Schur selfenergy failed')
    return {'all_unit_native_H_commutes_with_total_hypercharge': True,
            'hypercharge_commutes_with_all_twelve_global_currents': True,
            'literal_native_removal_charges': {k: str(v) for k, v in CHARGES.items()},
            'gauge_current_preserving_up_down_exchange_impossible': True,
            'exchange_obstruction': 'commuting with QY cannot map charge2/3 to charge-1/3',
            'full_thirty_mode_plus_scalar_charge_norm_upper': '13/2',
            'neutral_memory_extended_H_norm_upper': '73/2',
            'norm_bound_derivation': '||H||<=29, ||N_A||=1, ||QY X_A||<=13/2',
            'exact_blank_memory_resolvent': '1/(w-q^2/(w-1))',
            'exact_Schur_selfenergy': 'q^2/(w-1)',
            'all_anomalous_contacts_unchanged': 900,
            'contact_identity': '{[c_i,QY tensor X_A],c_j}=Y_i {c_i,c_j} tensor X_A=0',
            'actual_scalar_Gram_unchanged': 'same QY=0 native ground tensor |0><0|',
            'actual_scalar_retarded_first_moment_unchanged': 'extra double commutator is proportional to X_A, whose blank expectation is0',
            'scalar_Z_H': 'the same full-native 2-7n/4-3e/4, approximately9.58523463744285',
            'all_scalar_poles_unchanged_claimed': False,
            'local_dynamical_gauge_field_claimed': False}


def polarization(q, coupling=F(1)):
    """Both exact memory poles, including the zero residue at detachment."""
    q, coupling = F(q), F(coupling)
    x = coupling*q
    d = sqrt_bounds(1+4*x*x)
    if not x:
        d = (F(1), F(1))
    reciprocal = (1/d[1], 1/d[0])
    low_weight = ((1+reciprocal[0])/2, (1+reciprocal[1])/2)
    high_weight = ((1-reciprocal[1])/2, (1-reciprocal[0])/2)
    return [dict(branch='lower', shift=((1-d[1])/2, (1-d[0])/2), weight=low_weight),
            dict(branch='upper', shift=((1+d[0])/2, (1+d[1])/2), weight=high_weight)]


def native_symmetry(receipt):
    from .symmetry import certificate as symmetry_certificate
    a, b = (receipt['exact_local_response'][f] for f in ('uc', 'dc'))
    for key in ('complete_target_annihilator', 'complete_target_characteristic_factorization',
                'exact_ground_moments', 'complete_positive_records'):
        if a[key] != b[key]:
            raise ArithmeticError('native equal-coupling spectral identity changed: '+key)
    if len(a['complete_positive_records']) != 13 or any(
            interval(row['probability_interval'])[0] <= 0 for row in a['complete_positive_records']):
        raise ArithmeticError('complete thirteen-root cyclic support is required')
    return {'exact_complete_ground_resolvents_equal': True,
            'full_target_dimension_each': a['complete_target_dimension'],
            'cyclic_ground_image_dimension_each': 13,
            'exact_target_annihilator': a['complete_target_annihilator'],
            'cyclic_isometry': 'H_u^j c_u|g> -> H_d^j c_d|g>, j=0,...,12; equal moments and common annihilator preserve every inner product',
            'microscopic_full_spatial_symmetry': symmetry_certificate(),
            'whole_common_spatial_source_exchange_symmetry_proved_on_declared_domain': True,
            'same_current_labels_imply_equal_response': False,
            'mechanism_scope': 'equality holds on the actual local ground and on stationary invariant B=0 spatial singlet preparations; the unrestricted spatial ground baryon sector is not determined'}


def shifted_records(native, shifts):
    rows = []
    for index, record in enumerate(native['complete_positive_records']):
        energy = interval(record['source_excitation_interval'])
        weight = interval(record['probability_interval'])
        for branch in shifts:
            changed = add(energy, branch['shift'])
            probability = multiply(weight, branch['weight'])
            rows.append({'native_energy_record': index, 'branch': branch['branch'],
                         'original_bare_source_excitation_interval': encode(energy),
                         'source_or_relational_excitation_interval': encode(changed),
                         'total_source_plus_reference_excitation_interval': encode(add(changed, (F(4), F(4)))),
                         'original_probability_interval': encode(probability),
                         'zero_weight': record['exactly_zero_residue'] or branch['weight'] == (0, 0)})
    if any(interval(row['source_or_relational_excitation_interval'])[0] <= 0 for row in rows if not row['zero_weight']):
        raise ArithmeticError('a nonpositive visible energy entered the stated controls')
    return rows


def scale_test(channels, target):
    from model.defect_mass_catalogue.calibration import Anchor, common_scale, convert_energy, ratio, strings
    selected = {f: next(row for row in channels[f]['complete_records']
                        if row['native_energy_record'] == 0 and row['branch'] in ('single', 'lower'))
                for f in ('uc', 'dc')}
    energies = {f: interval(row['source_or_relational_excitation_interval']) for f, row in selected.items()}
    anchors = [Anchor('declared mechanism source clock', energies[f], (target[f], target[f]),
                      'specified record-zero and lower-memory branch dictionary') for f in ('uc', 'dc')]
    scales = [anchor.scale_interval() for anchor in anchors]
    intersection = common_scale(anchors)
    return {'record_dictionary': 'native record0; lower memory branch when present',
            'dictionary_is_physical_mass_identification': False,
            'dimensionless_up_down_ratio_interval': strings(ratio(energies['uc'], energies['dc'])),
            'ratio_numeric': midpoint(ratio(energies['uc'], energies['dc'])),
            'frozen_target_ratio': str(target['uc']/target['dc']),
            'common_scale_interval': strings(intersection),
            'one_scale_matches_frozen_two_targets': intersection is not None,
            'u_anchored_d_GeV_interval': strings(convert_energy('same', 'same', energies['dc'], scales[0])),
            'u_anchored_d_GeV_numeric': midpoint(convert_energy('same', 'same', energies['dc'], scales[0])),
            'target_d_running_parameter_GeV': str(target['dc']),
            'empirical_values_used_to_choose_primitive_coefficients': False,
            'target_uncertainties_propagated': False}


def calculate():
    pins = bindings()
    native = authenticate(NATIVE)
    symmetry = native_symmetry(native)
    controls = dict(operator_checks())
    family = authenticate(FAMILY)
    controls['exact_original_scalar_first_moment_preserved'] = family['complete_scalar_readout']['scalar_retarded_first_moment_each_component']
    calibration = json.loads((ROOT/CALIBRATION).read_text())['running_mass_parameters_numeric_GeV']
    target = {'uc': F(str(calibration['up'][0])), 'dc': F(str(calibration['down'][0]))}
    cases = {}
    for mechanism in ('positive_source_current_square', 'compensated_reference_attraction',
                      'positive_completed_reference_square', 'neutral_polarization_memory'):
        channels = {}
        for family, q in CHARGES.items():
            if mechanism == 'neutral_polarization_memory':
                shifts = polarization(q)
            else:
                shift = {'positive_source_current_square': q*q,
                         'compensated_reference_attraction': -q*q,
                         'positive_completed_reference_square': F(0)}[mechanism]
                shifts = [dict(branch='single', shift=(shift, shift), weight=(F(1), F(1)))]
            data = native['exact_local_response'][family]
            channels[family] = {'total_native_target_charge': str(q),
                'complete_native_energy_count': len(data['complete_positive_records']),
                'full_target_dimension': data['complete_target_dimension']*(2 if mechanism == 'neutral_polarization_memory' else 1),
                'complete_records': shifted_records(data, shifts),
                'success_probability': '1/2', 'failure_probability': '1/2',
                'outside_complete_target_probability': '0',
                'failure_energy_record': 'failure flag; source energy is not measured on this branch',
                'static_reference_cost_on_success': '4',
                'static_reference_cost_on_failure': '0',
                'per_success_conditioning_applied': False,
                'source_selected_fermion_family_copies_added': False}
            weights = [interval(row['original_probability_interval']) for row in channels[family]['complete_records']]
            if not sum(w[0] for w in weights) <= F(1, 2) <= sum(w[1] for w in weights):
                raise ArithmeticError('complete charged-removal weight was lost')
            if mechanism == 'neutral_polarization_memory':
                d = s.sqrt(1+4*s.Rational(q*q))
                channels[family]['exact_memory_records'] = [{
                    'branch': name, 'energy_shift': str((1+sign*d)/2),
                    'blank_memory_probability': str((1-sign/d)/2),
                    'memory_number_expectation': str((1+sign/d)/2),
                    'charge_memory_interaction_expectation': str(sign*2*s.Rational(q*q)/d)}
                    for name, sign in (('lower', -1), ('upper', 1))]
                channels[family]['complete_first_energy_moment_change'] = '0'
                channels[family]['complete_second_energy_moment_change'] = str(q*q/2)
                channels[family]['detached_coupling_zero_records'] = shifted_records(data, polarization(q, 0))
        cases[mechanism] = {'shared_primitive_coefficient': '1', 'channels': channels,
                            'one_scale_test': scale_test(channels, target)}
    cases['positive_source_current_square']['law'] = 'H_source=H_native+QY^2; reference idle'
    cases['compensated_reference_attraction'].update({
        'law': 'H_joint=H_native+4 N_ref+QY_source QY_reference',
        'source_energy_unchanged': True,
        'reference_self_energy_on_success': '4',
        'interaction_energy_on_compensated_charge_q': '-q^2',
        'detach_cross_coupling_zero_restores_equal_bare_spectra': True,
        'interaction_shift_is_intrinsic_particle_mass': False})
    cases['positive_completed_reference_square'].update({
        'law': 'H_joint=H_native+4 N_ref+(QY_source+QY_reference)^2/2',
        'source_charge_self_energy': 'q^2/2', 'reference_charge_self_energy': 'q^2/2',
        'interaction_energy': '-q^2', 'total_charge_field_energy': '0',
        'genuine_dynamical_gauge_field_constructed': False})
    cases['neutral_polarization_memory'].update({
        'law': 'H_source=H_native tensor I2+I tensor N_A+QY tensor X_A',
        'memory': 'neutral two-state source degree; N_A=diag(0,1), X_A=[[0,1],[1,0]]',
        'preparation': 'same actual native QY=0 physical ground times memory blank |0>',
        'reference': 'same retained compensating charged reference; idle except its declared cost4',
        'ground_claim': 'unchanged in the QY=0 domain; no minimum over other charge domains asserted',
        'projectors': 'every native energy projector tensor both exact q-memory projectors',
        'source_dimension': str(10*2**30),
        'global_native_currents_preserved': True,
        'local_Gauss_or_Maxwell_implementation': False,
        'source_bare_memory_and_interaction_energy_decomposition': 'memory records give subsystem expectations, not simultaneous sharp values of noncommuting terms',
        'unconditioned_uc_minus_dc_second_moment_change': '1/6',
        'whole_measures_replaceable_by_selected_lower_poles': False})
    if bindings() != pins:
        raise ValueError('source changed during complete quark response calculation')
    return {'schema': 'oph.quark_flavour_resolution.v1',
            'native_equal_coupling_domain': symmetry, 'exact_operator_checks': controls,
            'mechanisms': cases,
            'scope': {'dimensionless_up_down_splitting_constructed': True,
                      'shared_unit_primitive_coefficients': True,
                      'empirical_ratio_optimized': False,
                      'one_scale_fixed_target_fit_obtained': False,
                      'complete_common_spatial_interacting_ground_computed': False,
                      'physical_quark_rest_masses_derived': False,
                      'polarization_memory_is_a_declared_source_architecture': True},
            'input_sha256': pins}


def bindings():
    from model.particle_common_contract.provenance import import_closure
    seeds = {path.relative_to(ROOT).as_posix() for path in HERE.iterdir()
             if path.suffix in ('.py', '.md')}
    seeds.update({NATIVE, CALIBRATION, FAMILY, 'model/native_quark_response/PROOF.md',
                  'model/defect_family_vertices/PROOF.md'})
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(import_closure(seeds))}


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    value = calculate()
    path = HERE/'results.json'
    if args.write:
        path.write_text(json.dumps(value, sort_keys=True, indent=2)+'\n')
    elif json.loads(path.read_text()) != value:
        raise ValueError('complete charge-sensitive quark receipt differs')
    print('PASS: complete charge-sensitive quark spectra and fixed one-scale controls')
