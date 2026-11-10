"""Same-source elementary/composite filtered states and complete records."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from math import isqrt
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from model.native_weak_current_interface import model as weak
from model.native_higgs_spectral_join import model as higgs
from model.native_charge_resolution import certificate as resolution
from model.native_block_dirac import algebra as dirac
from model.particle_common_contract.provenance import import_closure

PARENTS = {
    'model/native_weak_current_interface/results.json': 'db5ae5c7b9f43762a4cdb7128761d8f2014cf3c60d0c9f324569d27d1d7b3f66',
    'model/native_higgs_spectral_join/results.json': '5749af81e1612f0658733b51acd52c2d08a1a69d174a3707989094d66d149973',
    'model/native_charge_resolution/results.json': 'f22c9a266a43de8aaae776b1293894081c3dc20b7d8e5678e30ef6af993dd19e',
    'model/native_block_dirac/results.json': '355dcc45b548203733d1ec6de5e9ffbd82512696976959464d2e1a0c8a70df14'}
SCALE = 2**128
RESEARCH = HERE/'research_reference.json'
RESEARCH_SHA = '13200ce38d74e97c6fee01f985032f6d5a3a3c18a0beb8d6aa8aa399d34561a8'


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def rounded(value, upper=False):
    value = F(value); n = value.numerator*SCALE; d = value.denominator
    return F(-(-n//d) if upper else n//d, SCALE)


def root(value, upper=False):
    value = F(value)
    if value < 0: raise ValueError('nonnegative radicand required')
    n = isqrt(value.numerator*SCALE**2//value.denominator)
    if upper and n*n*value.denominator < value.numerator*SCALE**2: n += 1
    return F(n, SCALE)


@lru_cache(None)
def sources():
    data = {}
    for name, expected in PARENTS.items():
        if digest(ROOT/name) != expected: raise ValueError('frozen parent differs: '+name)
        data[name] = json.loads((ROOT/name).read_text())
        higgs.validate(data[name]['input_sha256'])
    for name, replay in (
            ('model/native_weak_current_interface/results.json', weak.certificate),
            ('model/native_higgs_spectral_join/results.json', higgs.certificate),
            ('model/native_charge_resolution/results.json', resolution.certificate)):
        if replay() != data[name]: raise ValueError('fresh source replay differs: '+name)
    return data


def preparation(depth=3):
    if not isinstance(depth, int) or depth < 0: raise ValueError('nonnegative finite filter depth required')
    _, _, _, p = higgs.source()
    assert p['actual_local_preparation']['preparation_success_probability'] == '1'
    overlap = F(p['atomic_cube_trial']['exact_interacting_ground_trial_overlap_lower'])
    mixed = F(p['finite_energy_agreement_filter']['mixed_ground_excited_amplitude_upper'])
    M = int(p['finite_energy_agreement_filter']['explicit_integer_M'])
    r0 = (1-overlap)/overlap; c = r0+2*mixed*mixed
    assert 0 < r0 < c < 1
    error2 = 2*r0*c**depth
    Q = F(1)
    for j in range(depth): Q /= (1+r0*c**j)**(2**(depth-j))
    return {'binary_filter_depth': depth, 'native_trial_copies': 2**depth,
        'agreement_filter_count': 2**depth-1, 'each_filter_pointer_dimension': M,
        'initial_excited_odds_upper': str(r0), 'geometric_odds_contraction': str(c),
        'retained_purification_error_squared_upper': str(rounded(error2, True)),
        'retained_purification_error_upper': str(root(error2, True)),
        'complete_preparation_success_lower': str(rounded(Q)),
        'success_lower_formula': 'product_(j=0..d-1)(1+r0*c^j)^(-2^(d-j))',
        'deterministic_initial_trial': True, 'base_trial_success': '1',
        'filter_depth_counts_all_rounds': True, 'all_partners_and_filter_flags_retained': True,
        'exact_interacting_ground_projection_assumed': False}


def filtered_bounds(energy, weight, epsilon, source_remainder=F(720,2**34)):
    """Rational enclosures for a genuine spectral bin, never its energy alias."""
    L, U = map(F, energy); wlo, whi = map(F, weight); epsilon = F(epsilon)
    if not (0 < L <= U and 0 < wlo <= whi <= 1 and epsilon >= 0):
        raise ValueError('positive source energy and elementary weight required')
    e = (L+U)/2; delta = (U-L)/2; s = root(wlo)
    q = rounded(delta+source_remainder/s, True)
    if not (epsilon < s and q < e): raise ValueError('nonzero filtered-state margin required')
    actual_q = rounded(q+(2+e+q)*epsilon/(s-epsilon), True)
    if actual_q >= e: raise ValueError('actual preparation loses the alignment margin')
    a_low = s-epsilon; b_low = (e-q)*s-2*epsilon
    assert b_low > 0
    ideal_trace = rounded(q/e, True); actual_trace = rounded(actual_q/e, True)
    a_to_ideal = rounded(epsilon/s, True)
    b_to_ideal = rounded(2*epsilon/((e-q)*s), True)
    gamma = L*wlo-source_remainder
    assert gamma > 0
    out = {'energy_center': str(e), 'energy_halfwidth': str(delta),
        'source_EOM_remainder_norm_upper': str(source_remainder),
        'ideal_relative_vector_residual_upper': str(q),
        'actual_prepared_relative_vector_residual_upper': str(actual_q),
        'ideal_normalized_state_trace_distance_upper': str(ideal_trace),
        'actual_normalized_state_trace_distance_upper': str(actual_trace),
        'ideal_normalized_state_infidelity_upper': str(rounded((q/e)**2, True)),
        'actual_normalized_state_infidelity_upper': str(rounded((actual_q/e)**2, True)),
        'actual_elementary_vs_ideal_normalized_trace_upper': str(a_to_ideal),
        'actual_composite_vs_ideal_normalized_trace_upper': str(b_to_ideal),
        'ideal_elementary_in_bin_probability_lower': str(wlo),
        'ideal_composite_B_over_2_in_bin_probability_lower': str(rounded(((e-q)*s/2)**2)),
        'actual_elementary_in_bin_probability_lower': str(rounded(a_low*a_low)),
        'actual_composite_B_over_2_in_bin_probability_lower': str(rounded(b_low*b_low/4)),
        'actual_same_history_total_variation_upper': str(actual_trace),
        'actual_dressed_W_expectation_difference_upper': str(rounded(4*actual_trace, True)),
        'positive_cross_residue_lower': str(rounded(gamma)),
        'bound_applies_to_source_and_common_retained_preparation_environment': True}
    out['diagnostic_binary64'] = {key: float(F(value)) for key, value in out.items() if isinstance(value,str)}
    return out


def row(multiplier, prep):
    parent = sources()['model/native_charge_resolution/results.json']
    window = next(x for x in parent['nested_window_frontier']
                  if x['radius_in_units_of_full_hopping_norm'] == multiplier)
    energy = window['source_excitation_interval_relative_actual_native_ground']
    absolute = [F(x)-4 for x in window['absolute_source_plus_reference_bin']]
    weight = window['unnormalized_bin_weight_interval']
    result = filtered_bounds(energy, weight, prep['retained_purification_error_upper'])
    Q = F(prep['complete_preparation_success_lower'])
    for pump in ('elementary', 'composite_B_over_2'):
        p = F(result[f'actual_{pump}_in_bin_probability_lower'])
        result[f'unconditional_{pump}_in_bin_probability_lower'] = str(rounded(Q*p))
        result[f'joint_fair_setting_unconditional_{pump}_in_bin_probability_lower'] = str(rounded(Q*p/2))
    result.update({'radius_in_units_of_360_lambda': multiplier,
        'absolute_native_source_bin': list(map(str,absolute)),
        'excitation_enclosure_relative_actual_ground': energy,
        'elementary_unnormalized_bin_weight_interval': weight,
        'retains_original_certified_elementary_probability_lower': window['retains_original_certified_weight_lower'],
        'spectral_projector_uses_absolute_bin_not_excitation_enclosure': True,
        'exact_finite_spectral_service': True, 'implemented_bin_QPE_claimed': False,
        'unconditional_A_and_B_over_2_complete_laws_are_close_claimed': False,
        'subsequent_histories_may_be_adaptive_but_are_the_same_policy': True})
    return result


def selector_transfer(rows):
    algebra = dict(dirac.four_mode_code())
    spinors = dirac.spinor_checks()
    assert algebra['complete_four_mode_Fock_code_dimension'] == 16
    assert spinors['both_covariant_spin_sums_checked']
    result = {'selection': 'E_H_sharp maximizes Re <A_up g,P_E B_up g> in the radius2v absolute source bin; lower energy breaks ties',
        'actual_positive_cross_G1_eigen_doublet_exists': True,
        'source_excitation_interval': rows[0]['excitation_enclosure_relative_actual_ground'],
        'is_frozen_maximum_elementary_E_star': False,
        'old_Dirac_receipt_or_old_spectral_couplers_modified': False,
        'declared_free_architecture': 'four native copies per cell with selected G1 singlet occupation code and positive sqrt(m^2+L_a)-m kinetic matrix',
        'mass_parameter': 'm=E_H_sharp-Eg in source clock units',
        'same_constituent_h_ground_Gauss_and_full_reference': True,
        'spectral_matrix_units_and_full_range_kinetic_operator_must_be_redeclared': True,
        'parametric_four_mode_CAR_and_rotation_algebra_replayed': algebra,
        'both_frequency_spinor_sum_controls_replayed': True,
        'source_composite_equals_code_field_on_all_occupations': False,
        'interacting_native_Lorentz_or_SM_chirality_transfer': False}
    # The helper's literal E_star label describes its frozen instantiation.
    # The replay establishes only its selector-independent polynomial algebra.
    result['parametric_four_mode_CAR_and_rotation_algebra_replayed']['mass_term_on_code'] = 'm times total canonical particle-plus-antiparticle number'
    return result


def calculate():
    sources(); prep = preparation(3); rows = [row(m,prep) for m in (2,128)]
    native_norm = F(480)+F(360,2**34)
    Hnorm = 2**prep['binary_filter_depth']*(native_norm+120)
    return {'source': 'unchanged full native eight-cell h=H0+2^-34 T, actual ground g and complete old/dual/port reference',
        'current_interface': {'W_norm': '2', 'A_norm': '1', 'B_norm_upper': '2',
            'same_full_domain_EOM_replayed': True, 'scalar_eigenbackground_used': False},
        'physical_preparation': prep, 'filtered_state_settings': rows,
        'complete_instruments': {'settings': ['A', 'B/2'], 'each_setting_probability': '1/2',
            'each_setting_Kraus_maps': ['P C', '(I-P) C', 'sqrt(I-C^dagger C)'],
            'each_setting_probabilities': ['||P C Psi||^2', '||(I-P) C Psi||^2', '1-||C Psi||^2'],
            'retained_pump_spin_doublet_dimension': 4, 'energy_flag_dimension': 2,
            'setting_dimension': 2, 'regular_binary_cubic_reference_dimension': 48,
            'full_setting_pump_energy_alphabet_size': 16,
            'potentially_nonzero_labels_on_declared_empty_pump_input': 6,
            'preparation_failure_probability': '1-Q', 'full_joint_success_labels': 'Q times setting/pump/energy probabilities',
            'every_fine_preparation_record_and_partner_retained': True,
            'history_comparison_conditions_only_on_success_in_selected_bin': True,
            'source_only_marginals_are_regular_reference_twirls': True},
        'same_clock_apparatus': {'raw_full_kinematic_H_ext_norm_upper': str(Hnorm),
            'three_pulse_durations': ['1','1','1'], 'bare_wait': 't>=0', 'total_elapsed_source_time': '3+t',
            'driven_Hermitian_log_service': 'V=i Log(exp(-i H_ext tau) U)/tau-H_ext; phases in (-pi,pi]',
            'actual_pulse': 'exp(-i H_ext tau) U', 'each_control_norm_upper': str(Hnorm+F(22,7)),
            'net_source_work_absolute_upper': str(2*Hnorm),
            'gross_three_pulse_switching_work_upper': str(6*(Hnorm+F(22,7))),
            'finite_coefficient_error_complete_record_TV': 'sum_j tau_j ||delta L_j||',
            'autonomous_local_controller_or_efficient_native_circuit_claimed': False,
            'full_register_capacity_nats_formula': '2^d*(240 ln2+8 ln5+24 ln119+ln(2^31+1)+16 ln2)+(2^d-1)lnM+ln48+4ln2',
            'controller_coefficient_workspace_in_capacity_formula': False},
        'typed_free_Dirac_corollary': selector_transfer(rows)}


def certificate():
    value = calculate(); seeds = set(PARENTS)
    for parent in sources().values(): seeds.update(parent['input_sha256'])
    seeds.add('model/particle_common_contract/provenance.py')
    seeds.update(str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md'))
    seeds.add(str(RESEARCH.relative_to(ROOT)))
    value.update({'schema':'oph.native_higgs_filtered_state.v1',
        'input_sha256':{name:digest(ROOT/name) for name in sorted(import_closure(seeds))},
        'owned_source_sha256':{p.name:digest(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.md')}})
    if digest(RESEARCH)!=RESEARCH_SHA:raise ValueError('immutable filtered-state research reference differs')
    ignored={'schema','input_sha256','owned_source_sha256'}
    science=lambda row:{key:item for key,item in row.items() if key not in ignored}
    if science(value)!=science(json.loads(RESEARCH.read_text())):
        raise ValueError('fresh public filtered-state science differs from the frozen source calculation')
    value['independent_research_comparison']={
        'reference_sha256':RESEARCH_SHA,'all_scientific_fields_equal':True,
        'ignored_only_provenance_fields':sorted(ignored),'fresh_source_and_record_bounds_replayed':True}
    return value


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true'); group.add_argument('--check',action='store_true')
    args=parser.parse_args(); path=HERE/'results.json'
    if args.check:
        saved=json.loads(path.read_text()); higgs.validate(saved['input_sha256'])
        for name,sha in saved['owned_source_sha256'].items():
            if digest(HERE/name)!=sha: raise ValueError('owned source changed: '+name)
    actual=certificate()
    if args.write:path.write_text(json.dumps(actual,indent=2,sort_keys=True)+'\n')
    elif actual!=saved:raise ValueError('filtered-state source receipt differs')
    print('PASS: same-H filtered elementary/composite states, complete imperfect-preparation records and typed free code')


if __name__=='__main__':main()
