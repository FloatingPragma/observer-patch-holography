"""Source-word polynomial preparation, full records and Wilson-gap bounds.

The source is a periodic cubic gauge-coordinate field. Bounds are finite
and volume dependent. A prepared gap is not a vacuum or continuum theorem.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import itertools
import json
from pathlib import Path

import sympy as s

from model.particle_observer_outputs import field_route


ROOT = Path(__file__).resolve().parents[2]
SOURCE_INPUTS = (
    'model/particle_observer_outputs/field_route.py',
    'model/particle_observer_outputs/FIELD_ROUTE.md',
    'model/interacting_extension/chiral_lattice/source_data.py',
    'model/interacting_extension/chiral_lattice/light_branch.py',
    'model/interacting_extension/chiral_lattice/PROOF.md',
    'model/interacting_extension/chiral_electric/PROOF.md',
    'model/common_source/overlap_regulator.py',
    'model/observer_particles/results/representation_identification_operators.json',
    'model/observer_particles/results/representation_identification.json',
    'model/particle_observer_outputs/READOUT_CONTRACT.md',
)


def integer(value, name, minimum=1):
    if type(value) is not int or value < minimum:
        raise ValueError(name+' must be an integer at least '+str(minimum))
    return value


def rational(value, name, positive=False):
    if type(value) not in (int, str, F):
        raise ValueError(name+' requires exact rational data')
    result = F(value)
    if result < 0 or positive and result == 0:
        raise ValueError(name+' has the wrong sign')
    return result


def clean(matrix):
    return s.ImmutableMatrix(matrix).applyfunc(s.simplify)


@lru_cache(maxsize=1)
def native_casimirs():
    data = field_route.source_data.load()
    packet = json.loads((ROOT/SOURCE_INPUTS[7]).read_text())['even_identification_operators']
    def diagonal(name):
        answer = [s.S.Zero]*15
        for i, j, value in packet[name]['entries']:
            if i != j:
                raise ArithmeticError('native quadratic Casimir is not diagonal')
            answer[i] = s.Rational(value)
        return answer
    color, weak = diagonal('color_C2'), diagonal('weak_C2')
    q = data['charges']
    # Conventional fundamental normalization gives the actual M15 trace
    # metric 2 on both nonabelian factors and 120 on q=6Y.
    if sum(color)/8 != 2 or sum(weak)/3 != 2 or sum(x*x for x in q) != 120:
        raise ArithmeticError('the retained native trace metric differs')
    rows = []
    for block in data['blocks']:
        indices = block['indices']
        c3 = {color[i] for i in indices}
        c2 = {weak[i] for i in indices}
        if len(c3) != 1 or len(c2) != 1:
            raise ArithmeticError('a native charge block has split Casimirs')
        c3, c2 = c3.pop(), c2.pop()
        casimir = c3/2+c2/2+s.Rational(block['q6Y']**2, 120)
        rows.append({'q6Y': block['q6Y'], 'dimension': block['dimension'],
                     'color_C2': c3, 'weak_C2': c2, 'trace_metric_C2': casimir})
    if sum(row['trace_metric_C2'] for row in rows) != s.Rational(33, 10):
        raise ArithmeticError('source trace-metric Casimir sum differs')
    return rows


def polynomial_gap_plan(side, square_root_degree, target_bits=None):
    """Exact finite certificate for gauge-averaged source words, L=n².

    Haar normalization and suppression are kept as factored integer data;
    no large polynomial coefficient table is materialized.
    """
    side = integer(side, 'periodic side')
    n = integer(square_root_degree, 'square root of per-link degree')
    sites, edges, degree = side**3, 3*side**3, n*n
    # ceil(log2(1+8n)) without floating logarithms.
    cover_bits = (8*n).bit_length()
    exponent = 2*edges+900*edges*cover_bits-(n*n)//1080
    energy = F(29, 10)*edges*degree*degree
    if target_bits is not None:
        target_bits = integer(target_bits, 'target gap bits')
    return {'periodic_side': side, 'sites': sites, 'links': edges,
            'per_link_word_degree': degree, 'total_word_degree': edges*degree,
            'root_degree': n,
            'normalized_gap_bad_event': 'sigma_min(H_W)<1/2',
            'Haar_ball_lower': {'base': 1+8*n, 'inverse_power': 450},
            'unnormalized_state_norm_squared_lower': {
                'rational_base': '59/60', 'rational_power': 2*edges,
                'integer_base': 1+8*n, 'inverse_power': 900*edges},
            'bad_probability_upper': {'base': 2, 'exponent': min(0, exponent)},
            'uncapped_binary_exponent': exponent,
            'electric_energy_upper': energy,
            'code_dimension_upper': {'base': 451, 'exponent': edges*degree},
            'logical_preparation_qubits_upper': 9*edges*degree,
            'target_bits': target_bits,
            'target_met': target_bits is None or exponent <= -target_bits,
            'finite_source_word_code': True,
            'gauge_invariant_amplitude': True,
            'prepared_state_is_claimed_vacuum': False,
            'large_preparation_coefficients_materialized': False}


def choose_preparation(side, target_bits):
    integer(target_bits, 'target gap bits')
    n = 1
    attempts = 0
    while True:
        value = polynomial_gap_plan(side, n, target_bits)
        attempts += 1
        if value['target_met']:
            return {**value, 'doubling_steps': attempts}
        n *= 2


def same_clock_interval(plan, degree=128, spacing='1', plaquette='1/100',
                        half_gap_bits=10, occupied_modes=60):
    """Short-time gap leakage and polynomial/exact-sign record comparison.

    The initial state is the normalized finite word polynomial, with all 60
    modes at one cell filled. No intermediate measurements are inserted.
    The energy bound uses the same positive electric Casimir and H_m.
    """
    b = integer(half_gap_bits, 'half gap bits')
    count = integer(occupied_modes, 'fixed conserved occupation count', 0)
    if count > 60*plan['sites']:
        raise ValueError('occupation count exceeds the full source Fock carrier')
    if plan['uncapped_binary_exponent'] > -2*b:
        raise ValueError('preparation does not certify the requested initial gap weight')
    a = rational(spacing, 'spacing', True)
    kappa = rational(plaquette, 'plaquette coefficient')
    inherited = field_route.polynomial_bounds(degree, a)
    epsilon = inherited['sign_error_on_admitted_gap']
    bounded = count*inherited['global_one_body_H_norm_upper']+kappa*plan['links']
    positive_energy = plan['electric_energy_upper']+2*bounded
    eta = F(1, 2**(2*b))
    duration = eta/(2*positive_energy)
    good_bound = 4*eta
    comparison = F(count, 1)*duration/a*(epsilon+(2+epsilon)*F(2, 2**b))
    return {'Legendre_degree': degree, 'spacing': a,
            'conserved_fermion_number': count,
            'bounded_action_norm_on_fixed_number_sector': bounded,
            'positive_shifted_H_energy_upper': positive_energy,
            'initial_small_gap_probability_upper': eta,
            'common_duration': duration,
            'small_gap_probability_upper_through_interval': good_bound,
            'complete_final_record_TV_poly_vs_exact_sign_upper': min(F(1), comparison),
            'comparison_changes_only_sign_polynomial': True,
            'same_electric_plaquette_state_and_clock': True,
            'intermediate_measurements_covered': False,
            'preparation_service_error_included': False,
            'uniform_positive_time_as_volume_grows_claimed': False}


def plaquette_moment_plan(side, target_bad_bits):
    """Polynomial amplitude concentrating the positive plaquette cost.

    The chosen moment allowance is tailored to the inherited HJL barrier
    Bad => S>=1/48000. All large resource powers remain factored.
    """
    side = integer(side, 'periodic side')
    bits = integer(target_bad_bits, 'target probability bits')
    edges = 3*side**3
    target = F(1, 2**bits)
    threshold_bits = 1
    while F(edges, 2*4**threshold_bits) > target/(4*48000):
        threshold_bits += 1
    tail_bits = 1
    while F(2*edges, 2**tail_bits) > target/(4*48000):
        tail_bits += 1
    n = 1
    while True:
        plan = polynomial_gap_plan(side, n)
        exponent = (2*edges+900*edges*(8*n).bit_length()
                    -(n*n)//(480*4**threshold_bits))
        if exponent <= -tail_bits:
            break
        n *= 2
    threshold = F(1, 2**threshold_bits)
    tail_upper = F(1, 2**tail_bits)
    moment = edges*threshold**2/2+2*edges*tail_upper
    if moment > target/(2*48000):
        raise ArithmeticError('the prepared plaquette moment exceeds its allowance')
    return {**plan, 'requested_all_time_bad_weight': target,
            'plaquette_norm_threshold': threshold,
            'large_plaquette_probability_upper': tail_upper,
            'large_plaquette_uncapped_binary_exponent': exponent,
            'positive_plaquette_sum_mean_upper': moment,
            'HJL_bad_gap_plaquette_sum_lower': F(1, 48000)}


def energy_confinement_plan(side=1, target_bad_bits=40, degree=128,
                            spacing='1', occupied_modes=60):
    """A distinct declared high-plaquette-cost member with all-time gap control.

    H_m=H_E+K_m-kappa sum ReTr(U_p)/15; adding kappa*plaquette_count
    yields H'_m=H_E+K_m+kappa*S. Energy conservation, not repeated
    filtering or time rescaling, supplies the bound.
    """
    plan = plaquette_moment_plan(side, target_bad_bits)
    a = rational(spacing, 'spacing', True)
    count = integer(occupied_modes, 'fixed fermion number', 0)
    if count > 60*plan['sites']:
        raise ValueError('occupation count exceeds the source carrier')
    poly = field_route.polynomial_bounds(degree, a)
    kinetic = count*poly['global_one_body_H_norm_upper']
    target = plan['requested_all_time_bad_weight']
    energy = plan['electric_energy_upper']
    kappa = 2*48000*(energy+2*kinetic)/target
    bound = 48000*(plan['positive_plaquette_sum_mean_upper']
                   +(energy+2*kinetic)/kappa)
    if bound > target:
        raise ArithmeticError('conserved energy does not certify the requested gap weight')
    small_kappa = min(F(1), 48000*(plan['positive_plaquette_sum_mean_upper']
                                   +100*(energy+2*kinetic)))
    return {**plan, 'Legendre_degree': degree, 'spacing': a,
            'conserved_fermion_number': count,
            'kinetic_operator_norm_upper': kinetic,
            'declared_plaquette_coefficient': kappa,
            'source_sign_approximation_error': poly['sign_error_on_admitted_gap'],
            'all_time_small_gap_probability_upper': bound,
            'all_time_target_upper': target,
            'same_preparation_action_and_clock_during_observation': True,
            'plaquette_coefficient_differs_from_1_100_control': True,
            'no_repeated_gap_filtering': True,
            'same_bound_for_compressed_bad_effect_on_every_containing_word_code': True,
            'compressed_bad_effect_is_not_truncated_Wilson_gap': True,
            'kappa_1_100_energy_bound_on_this_preparation': small_kappa,
            'small_kappa_failure_is_inconclusive_not_a_dynamical_no_go': True,
            'uniform_finite_coupling_refinement_claimed': False}


def all_time_final_record_bound(plan, duration):
    """One final complete detector after any finite declared model time."""
    t = rational(duration, 'duration')
    eta = plan['all_time_target_upper']
    # The planner target is exactly dyadic. This is an upper square root
    # even when its exponent is odd.
    bits = eta.denominator.bit_length()-1
    sqrt_upper = F(1, 2**(bits//2))
    epsilon = plan['source_sign_approximation_error']
    error = (plan['conserved_fermion_number']*t/plan['spacing']
             *(epsilon+(2+epsilon)*sqrt_upper))
    return {'duration': t,
            'complete_final_record_TV_poly_vs_exact_sign_upper': min(F(1), error),
            'all_time_gap_bound_used': True,
            'same_electric_plaquette_state_and_clock': True,
            'spatial_locality_or_cutoff_error_included': False,
            'intermediate_measurement_backaction_covered': False}


@lru_cache(maxsize=1)
def one_cell_seed():
    """Exact L=1 preparation on actual normalized source characters.

    The 11-character code is a subspace of L²(G); its orthogonal complement
    is retained by identity in preparation and by a detector output.
    """
    irreps = [{'label': 'trivial', 'dimension': 1, 'q6Y': 0,
               'trace_metric_C2': s.S.Zero, 'probability': s.Rational(90, 91)}]
    for row in native_casimirs():
        for conjugate in (False, True):
            irreps.append({**row,
                'label': str(row['q6Y'])+('_conjugate' if conjugate else '_native'),
                'q6Y': -row['q6Y'] if conjugate else row['q6Y'],
                'conjugate': conjugate, 'probability': s.Rational(1, 910)})
    state = s.Matrix([30]+[1]*10)/s.sqrt(910)
    blank = s.eye(11)[:, 0]
    a = (blank.T*state)[0]
    difference = state-blank
    unitary = clean(s.eye(11)-difference*difference.T/(1-a))
    if clean(unitary.T*unitary) != s.eye(11) or clean(unitary*blank) != state:
        raise ArithmeticError('complete finite character preparation failed')
    probabilities = [s.simplify(abs(x)**2) for x in state]
    if probabilities != [r['probability'] for r in irreps]:
        raise ArithmeticError('original character probabilities changed')
    one_energy = sum(r['probability']*r['trace_metric_C2']/2 for r in irreps)
    if one_energy != s.Rational(33, 9100):
        raise ArithmeticError('source trace-metric electric energy differs')
    return {'irreps': irreps, 'state': state, 'unitary': unitary,
            'one_link_energy': one_energy}


def seed_records():
    """All 12³ electric-irrep records, including the outside-menu outcomes."""
    data = one_cell_seed()
    menu = data['irreps']+[{'label': 'outside_menu', 'probability': s.S.Zero,
                           'trace_metric_C2': None, 'dimension': None}]
    records = []
    mean = s.S.Zero
    for indices in itertools.product(range(12), repeat=3):
        rows = [menu[i] for i in indices]
        probability = s.prod(row['probability'] for row in rows)
        energy = None if any(i == 11 for i in indices) else sum(
            row['trace_metric_C2']/2 for row in rows)
        if energy is not None:
            mean += probability*energy
        records.append({'outcomes': [row['label'] for row in rows],
                        'probability': probability, 'electric_energy': energy})
    if sum(row['probability'] for row in records) != 1 or mean != s.Rational(99, 9100):
        raise ArithmeticError('complete seed detector lost its source energy law')
    return {'per_link_menu': menu, 'character_code_dimension': 11**3,
            'complete_record_count': len(records),
            'positive_record_count': sum(r['probability'] != 0 for r in records),
            'zero_record_count': sum(r['probability'] == 0 for r in records),
            'original_records': records, 'mean_electric_energy': mean,
            'postselection_used': False,
            'preparation_extended_by_identity_on_orthogonal_complement': True,
            'entire_irrep_multiplicity_in_each_detector_projector': True,
            'small_gap_weight_certified_by_L1_seed': False,
            'prepared_subspace_claimed_invariant_under_full_field_H': False}


@lru_cache(maxsize=1)
def native_gap_closure():
    """Exact native weak-doublet Wilson symbol, with every M15 multiplicity.

    This is a coordinate symbol, not a normalizable point-coordinate state.
    Its neighborhood proves that Haar preparation does not ensure a gap.
    """
    data = field_route.source_data.load()
    raw_gamma, _ = data['gamma']
    gamma = tuple(s.Matrix(g).applyfunc(lambda x: s.nsimplify(x, rational=True))
                  for g in raw_gamma)
    beta = gamma[3]
    pauli = (s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]))
    links = tuple((s.eye(2)+s.I*s.sqrt(3)*p)/2 for p in pauli)+(s.eye(2),)
    shifts = []
    for mu, link in enumerate(links):
        pminus, pplus = (s.eye(4)-gamma[mu])/2, (s.eye(4)+gamma[mu])/2
        shifts.append(s.kronecker_product(pminus, link)
                      +s.kronecker_product(pplus, link.H))
    kernel = clean(s.kronecker_product(beta, s.eye(2))*(2*s.eye(8)-sum(shifts, s.zeros(8))))
    if kernel != kernel.H:
        raise ArithmeticError('actual source Wilson symbol is not Hermitian')
    squared = clean(kernel*kernel)
    if squared.eigenvals() != {s.S.Zero: 4, s.Integer(3): 4}:
        raise ArithmeticError('native weak-link gap-closure control changed')
    weak_dimensions = sum(r['dimension'] for r in native_casimirs() if r['weak_C2'] != 0)
    if weak_dimensions != 8:
        raise ArithmeticError('the native weak multiplicities changed')
    spectral = {0: 2*weak_dimensions, 1: 4*(15-weak_dimensions), 3: 2*weak_dimensions}
    return {'spatial_beta': 'gamma_tuple[3], not gamma5',
            'weak_links': ['(I+i sqrt3 sigma1)/2', '(I+i sqrt3 sigma2)/2', 'I'],
            'full_Dirac_M15_dimension': 60,
            'H_W_squared_spectral_ranks': spectral,
            'symbol_auxiliary_trace_probe_probabilities': {k: s.Rational(v, 60) for k, v in spectral.items()},
            'actual_quantum_point_coordinate_preparation_claimed': False,
            'bare_Haar_reference_guarantees_admitted_gap': False,
            'gap_is_not_a_gauge_choice': True}


def serial(value):
    if isinstance(value, dict):
        return {str(k): serial(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serial(v) for v in value]
    if isinstance(value, (F, s.Basic)):
        return str(value)
    return value


def certificate():
    plan = choose_preparation(1, 20)
    clock = same_clock_interval(plan)
    step = max(2*clock['Legendre_degree']+1, 4)
    integrated = clock['bounded_action_norm_on_fixed_number_sector']*clock['common_duration']
    # Preserve a rigorous simple dyadic ceiling rather than raising the
    # long exact time denominator to every term in the factorial tail.
    dyadic = F(1)
    for _ in range(128):
        if dyadic/2 < integrated:
            break
        dyadic /= 2
    if integrated > dyadic:
        raise ArithmeticError('declared dyadic interaction bound is too small')
    cutoff = field_route.cutoff_bound(plan['total_word_degree']+4*step,
        plan['total_word_degree'], step, dyadic)
    confined = energy_confinement_plan()
    return serial({'schema': 'oph.gauge_coordinate_preparation.v1',
        'source_program': 'same periodic cubic Dirac/M15 coordinate field',
        'electric_metric': '-Tr_M15(T_a T_b), electric coefficient 1/2',
        'native_Casimirs': native_casimirs(),
        'one_cell_complete_seed_records': seed_records(),
        'exact_native_gap_closure': native_gap_closure(),
        'certified_small_gap_preparation': plan,
        'common_clock_accuracy_interval': clock,
        'same_source_total_degree_cutoff': cutoff,
        'distinct_high_plaquette_cost_preparation': confined,
        'all_time_configuration_final_record_at_time_one': all_time_final_record_bound(confined, 1),
        'finite_unitary_preparation_not_retry_postselection': True,
        'higher_degree_polynomial_coefficients_or_tape_executed': False,
        'coordinate_parent_is_analytic_comparison': True,
        'finite_evolution_uses_compressed_parent_polynomial': True,
        'small_gap_event_is_not_gap_of_truncated_Wilson_matrix': True,
        'periodic_boundary_required_for_free_gap_bound': True,
        'physical_vacuum_selected': False,
        'fixed_coupling_uniform_refinement_state_claimed': False,
        'interacting_relativistic_continuum_claimed': False,
        'physical_particle_assignments': [], 'source_predicted_GeV_masses': []})
