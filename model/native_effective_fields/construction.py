"""Source exchange fields, controlled transfer logarithm and actual currents.

The field carrier is the existing fixed twelve-port binary total-two sector.
Euclidean native-attempt parameter, analytic unitary parameter and observer
clock readings are separate. No continuum, physical dimension, or empirical
agreement follows from this finite construction.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import combinations
import math
from pathlib import Path

import numpy as np
import sympy as sp

from ..native_clock import clock
from ..native_fluctuations import experiment as fluctuations
from ..native_fluctuations import source as currents
from ..observer_reconstruction import native_geometry
from ..observer.scalar_response.experiment import canonical
from ..primitive_observer import runtime, transfer
from ..federated_quantum_sector import construction as quantum
from . import charged_reference, charged_reference_source


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*clock.SOURCE_PATHS, *fluctuations.SOURCE_PATHS,
    'model/observer_reconstruction/native_geometry.py',
    'model/native_effective_fields/__init__.py', 'model/native_effective_fields/__main__.py',
    'model/native_effective_fields/construction.py',
    'model/native_effective_fields/charged_reference.py',
    'model/native_effective_fields/charged_reference_source.py',
    'model/observer/source_scalar/compiler.py',
    'model/observer/source_discrimination/circuit.py',
    'code/source_selection_model/response.json',
    'code/m1_source_realization/receipt.json',
    'model/federated_quantum_sector/construction.py',
    'Lean/SourceLaw/TransferMatrixLift.lean', 'Lean/Dynamics/PrivateInner.lean'))))
CURRENT_VALUES = (-1, 0, 1)


@lru_cache(maxsize=1)
def source_kernel():
    result = transfer.native_counts()
    transfer.exact_audit(result)
    return result


def _order(value):
    if type(value) is not int or not 0 <= value <= 24:
        raise ValueError('an exact logarithm order from zero through twenty-four required')
    return value


def seam_term(seam):
    """Restriction of n_a+n_b-2 n_a n_b-(raise_a lower_b+h.c.)."""
    kernel = source_kernel()
    if (type(seam) is not tuple or len(seam) != 2
            or any(type(port) is not int for port in seam) or seam not in kernel.seams):
        raise ValueError('one declared oriented source seam required')
    a, b = seam
    index = {state: i for i, state in enumerate(kernel.states)}
    matrix = [[0]*66 for _ in range(66)]
    for i, state in enumerate(kernel.states):
        matrix[i][i] = state[a]+state[b]-2*state[a]*state[b]
        if state[a] != state[b]:
            changed = list(state)
            changed[a], changed[b] = state[b], state[a]
            matrix[i][index[tuple(changed)]] = -1
    return tuple(map(tuple, matrix))


@lru_cache(maxsize=1)
def generator():
    """Immutable exact integer exchange generator G=60(I-K), not a fit."""
    kernel = source_kernel()
    terms = [seam_term(seam) for seam in kernel.seams]
    result = tuple(tuple(sum(term[i][j] for term in terms) for j in range(66)) for i in range(66))
    if any(result[i][j] != 60*int(i == j)-kernel.counts[i][j]
           for i in range(66) for j in range(66)):
        raise ArithmeticError('local occupation/exchange terms do not reproduce the native kernel')
    return result


def transfer_matrix():
    """Immutable exact positive native configuration transfer K on 66 states."""
    return tuple(tuple(F(value, 60) for value in row) for row in source_kernel().counts)


@lru_cache(maxsize=25, typed=True)
def log_polynomial(order):
    """Exact sum_{k=1..order}(G/60)^k/k, with no floating matrix products."""
    order = _order(order)
    g = np.asarray(generator(), dtype=object)
    power = np.eye(66, dtype=object)
    total = [[F(0)]*66 for _ in range(66)]
    for k in range(1, order+1):
        power = power@g
        denominator = k*60**k
        for i in range(66):
            for j in range(66):
                total[i][j] += F(int(power[i, j]), denominator)
    return tuple(map(tuple, total))


def log_remainder(order):
    """Operator bound on H-log_polynomial from 0<=G/60<=I/3.

    Positive spectral calculus gives 0<=H-H_m<=r_m I with
    r_m=(1/3)^(m+1)/((m+1)(1-1/3)). This does not certify a later
    floating matrix product or scalar conversion of H_m.
    """
    order = _order(order)
    return F(1, 2*(order+1)*3**order)


def local_field_audit():
    kernel, g = source_kernel(), np.asarray(generator(), dtype=np.int64)
    first = kernel.seams[0]
    sharing = next(e for e in kernel.seams if len(set(e)&set(first)) == 1)
    disjoint = next(e for e in kernel.seams if not set(e)&set(first))
    a, b, c = (np.asarray(seam_term(e), dtype=np.int64) for e in (first, sharing, disjoint))
    if np.any(a@c-c@a) or not np.any(a@b-b@a):
        raise ArithmeticError('the exact seam locality/noncommutativity controls changed')
    adjacency = [sum(state[a]*state[b] for a, b in kernel.seams) for state in kernel.states]
    deleted = g+2*np.diag(adjacency)
    deleted_row_sums = tuple(map(int, deleted.sum(axis=1)))
    if set(deleted_row_sums) != {0, 2} or sum(deleted_row_sums) != 60:
        raise ArithmeticError('deleting the density interaction lost its exact native counterexample')
    return {'occupation_operators': 'n_i|x>=x_i|x>, with x_i in{0,1} and sum_i x_i=2',
        'seam_generator': 'n_a+n_b-2*n_a*n_b-(raise_a*lower_b+raise_b*lower_a)',
        'sum_of_30_local_terms_equals_60_I_minus_native_counts': True,
        'full_matrix_entries_checked': 66**2, 'declared_ports': 12,
        'disjoint_seam_commutator': '0', 'overlapping_seam_commutator_nonzero': True,
        'commutator_witness_seams': [list(first), list(sharing), list(disjoint)],
        'source_coefficient_of_density_product_in_G': '-2',
        'deleted_density_interaction_constant_vector_is_ground': False,
        'deleted_density_one_attempt_total_mass_from_uniform': str(1-F(sum(deleted_row_sums), 66*60)),
        'deletion_control_is_not_a_stochastic_replacement_law': True,
        'local_generator_operator_bounds': '0 <= G <= 20 I',
        'exact_transfer_operator_bounds': '(2/3) I <= K <= I'}


def logarithm_audit():
    matrix = np.asarray(generator(), dtype=np.int64)
    h, vacuum, spectral = transfer.spectral_reconstruction(source_kernel())
    rows = []
    witness = None
    for order in (1, 2, 4, 8, 16):
        polynomial, remainder = log_polynomial(order), log_remainder(order)
        residual = float(np.linalg.norm(h-np.asarray(polynomial, float), 2))
        if residual > float(remainder)+2e-14:
            raise ArithmeticError('diagnostic transfer logarithm disagrees with the exact analytic enclosure')
        rows.append({'order': order, 'exact_operator_remainder_upper': str(remainder),
            'polynomial_sha256': hashlib.sha256(canonical([[str(x) for x in row] for row in polynomial]).encode()).hexdigest(),
            'floating_spectral_residual_diagnostic': residual,
            'unitary_operator_error_for_parameter_t': f'abs(t)*({remainder})',
            'positive_Euclidean_semigroup_error_for_parameter_s': f's*({remainder}), s>=0'})
        if order == 16:
            for i in range(66):
                for j in range(i+1, 66):
                    if matrix[i, j] == 0 and polynomial[i][j] > remainder:
                        witness = {'state_indices': [i, j], 'one_attempt_exchange_entry': 0,
                            'H_entry_strictly_positive_lower': str(polynomial[i][j]-remainder),
                            'H_entry_upper': str(polynomial[i][j]+remainder)}
                        break
                if witness:
                    break
    if witness is None:
        raise ArithmeticError('the certified nonlocal logarithm witness disappeared')
    trace = sum(transfer_matrix()[i][i] for i in range(66))
    if trace != 56:
        raise ArithmeticError('the transfer-reference normalization changed')
    return {'polynomial_enclosures': rows, 'nonlocal_logarithm_witness': witness,
        'spectral_diagnostics': spectral,
        'faithful_transfer_reference': {'definition': 'sigma_K=K/56', 'trace_K': str(trace),
            'modular_generator': '-log(sigma_K)=H+log(56) I',
            'state_action_convention': 'sigma_K^(-is) rho sigma_K^(is)=exp(is H) rho exp(-is H)',
            'reference_is_native_stationary_density_I66_over66': False,
            'reference_physically_prepared_by_native_records': False,
            'reference_uniquely_selected_by_OPH': False},
        'unitary_continuation_is_native_CPTP_repair': False,
        'finite_positive_generator_is_a_relativistic_mass_gap': False}


def spatial_locality_audit():
    """Join the record-recovered finite metric to the native field terms."""
    geometry = native_geometry.recover(native_geometry.execute_calibration())
    if not geometry['rank_three_unit_gram_recovered']:
        raise ArithmeticError('record-recovered rank-three metric is unavailable')
    encoded = geometry['gram_qsqrt5']
    kernel = source_kernel()
    edge_set = {tuple(sorted(edge)) for edge in kernel.seams}
    shell_counts = {'edge': 0, 'other': 0}
    for i in range(12):
        for j in range(i+1, 12):
            dot = encoded[i][j]
            if (i, j) in edge_set:
                if dot != ['0', '1/5'] or encoded[j][i] != ['0', '1/5']:
                    raise ArithmeticError('a native exchange seam is not a shortest-shell metric pair')
                shell_counts['edge'] += 1
            else:
                if dot not in (['0', '-1/5'], ['-1', '0']) or encoded[j][i] != dot:
                    raise ArithmeticError('a non-edge is not on a longer recovered spatial shell')
                shell_counts['other'] += 1
    if shell_counts != {'edge': 30, 'other': 36}:
        raise ArithmeticError('the rank-three metric does not separate source edges from other shells')
    # Independently verify that each local contribution is exactly the
    # identity-minus-swap term on that geometrically shortest edge.
    index = {state: i for i, state in enumerate(kernel.states)}
    identity = np.eye(66, dtype=np.int64)
    summed = np.zeros((66, 66), dtype=np.int64)
    for seam in kernel.seams:
        swap = np.zeros((66, 66), dtype=np.int64)
        a, b = seam
        for column, state in enumerate(kernel.states):
            changed = list(state)
            changed[a], changed[b] = changed[b], changed[a]
            swap[index[tuple(changed)], column] = 1
        term = np.asarray(seam_term(seam), dtype=np.int64)
        if not np.array_equal(term, identity-swap):
            raise ArithmeticError('a native local field term escaped its shortest spatial edge')
        summed += term
    if not np.array_equal(summed, np.asarray(generator(), dtype=np.int64)):
        raise ArithmeticError('the full native field generator is not the sum of its metric-edge terms')
    return {'metric_record_count': geometry['record_count'],
        'rank_three_unit_gram_certified': True,
        'edge_squared_chord_length': '2-2*sqrt(5)/5 in recovered unit Gram',
        'pair_shell_counts': shell_counts,
        'native_exchange_terms_checked': len(kernel.seams),
        'every_native_exchange_term_is_geometrically_edge_local': True,
        'native_generator_is_exact_sum_of_local_edge_terms': True,
        'logarithm_can_have_longer_range_terms': True,
        'physical_length_scale_or_continuum_claimed': False}


def _ports(ports):
    if (type(ports) is not tuple or not 1 <= len(ports) <= 4
            or any(type(p) is not int or not 0 <= p < 12 for p in ports)):
        raise ValueError('one through four declared integer ports required')
    return ports


def exact_current_law(ports):
    """All signed words and terminal integer counts; no initial load access."""
    ports = _ports(ports)
    kernel = source_kernel()
    rows = {(): np.ones(66, dtype=np.int64)}
    for step, port in enumerate(ports, start=1):
        marks = tuple(np.asarray(m, dtype=np.int64) for m in currents.marked_counts(kernel, port))
        rows = {word+(value,): row@mark for word, row in rows.items()
                for value, mark in zip(CURRENT_VALUES, marks, strict=True)}
        # The integer total is <=66*60^4<2^30, so int64 has ample exact range.
        if sum(int(row.sum()) for row in rows.values()) != 66*60**step:
            raise ArithmeticError('complete current words lost native path mass')
    return {word: tuple(map(int, row)) for word, row in rows.items()}, 66*60**len(ports)


def current_moments(ports):
    ports = _ports(ports)
    if len(ports) != 4:
        raise ValueError('the connected four-point observable requires four ports')
    rows, denominator = exact_current_law(ports)
    def expectation(indices):
        return sum((F(sum(row)*math.prod(word[i] for i in indices), denominator)
                    for word, row in rows.items()), F(0))
    means = tuple(expectation((i,)) for i in range(4))
    if any(means):
        raise ArithmeticError('stationary recorded increments must have zero mean')
    covariances = {ij: expectation(ij) for ij in combinations(range(4), 2)}
    fourth = expectation(range(4))
    gaussian = (covariances[0, 1]*covariances[2, 3]+covariances[0, 2]*covariances[1, 3]
                +covariances[0, 3]*covariances[1, 2])
    return {'ports': list(ports), 'means': list(map(str, means)),
        'pair_covariances': {f'{i},{j}': str(v) for (i, j), v in covariances.items()},
        'four_point': str(fourth), 'Gaussian_matching_two_point_prediction': str(gaussian),
        'connected_four_point': str(fourth-gaussian),
        'complete_word_count': len(rows), 'nonzero_word_count': sum(bool(sum(row)) for row in rows.values()),
        'initial_load_disclosure_needed': False, 'hidden_initial_configuration_measured': False}


def current_program(ports):
    ports = _ports(ports)
    instruments = {p: currents.current_instrument(source_kernel(), p) for p in set(ports)}
    return runtime.Program('native_effective_field_four_current_probe', 66, np.eye(66)/66,
        tuple(runtime.Event(f'signed_current_port_{port}_at_attempt_{step+1}', instruments=(instruments[port],))
              for step, port in enumerate(ports)),
        binding={'ports': list(ports), 'parameter': 'successive native attempts',
                 'initial_law': 'stationary I66/66; actual configuration not disclosed',
                 'record': 'delta n_port from selected seam and committed native current',
                 'coherent_reference_or_pulse_execution_claimed': False})


def execute_currents(ports, *, checkpoint=False):
    ports = _ports(ports)
    if type(checkpoint) is not bool:
        raise ValueError('Boolean checkpoint option required')
    p = current_program(ports)
    history = runtime.execute(p)
    expected, denominator = exact_current_law(ports)
    max_state = max_probability = 0.
    laws = []
    for branch in history.branches:
        word = tuple(CURRENT_VALUES[value] for value in branch.record)
        row = expected[word]
        target = np.diag([float(F(x, denominator)) for x in row])
        max_state = max(max_state, float(np.linalg.norm(branch.density-target, 'fro')))
        max_probability = max(max_probability, abs(branch.probability-float(F(sum(row), denominator))))
        if dict(branch.resources) != {'native_attempts': len(ports), 'retained_signed_current_records': len(ports)}:
            raise ArithmeticError('complete current resource costs changed')
        laws.append({'word': list(word), 'exact_probability': str(F(sum(row), denominator)),
                     'executed_probability_diagnostic': branch.probability})
    if len(laws) != 3**len(ports):
        raise ArithmeticError('a zero or nonzero current word was discarded')
    resume = None
    if checkpoint:
        prefix = runtime.execute(p, stop_after=min(2, len(ports)))
        restored = runtime.checkpoint_from_json(p, runtime.checkpoint_to_json(prefix))
        resumed = runtime.execute(p, checkpoint=restored)
        resume = currents.bounded(runtime.complete_distance(history, resumed), 'current checkpoint')
    return {'ports': list(ports), 'complete_record_law': laws,
        'full_terminal_dimension': 66, 'all_zero_outcome_words_retained': True,
        'exact_terminal_counts_sha256': hashlib.sha256(canonical([[list(word), list(row)]
            for word, row in expected.items()]).encode()).hexdigest(),
        'full_terminal_density_residual_diagnostic': currents.bounded(max_state, 'four-current state'),
        'complete_probability_residual_diagnostic': currents.bounded(max_probability, 'four-current probability'),
        'serialized_checkpoint_complete_distance_diagnostic': resume,
        'native_C_CP_disturbance_retained': True, 'source_pulse_words_emitted': False}


def record_spectral_reduction():
    """Minimal linear transfer representation of accessible absolute-current words."""
    reduced = clock.lumping(source_kernel())
    p, d = sp.Matrix(reduced['matrix']), sp.diag(*reduced['stationary'])
    marks = [sp.Matrix(m) for m in reduced['marked']]
    one, prior = sp.ones(3, 1), sp.Matrix([reduced['stationary']])
    h = marks[1]*one
    mean = sp.Rational(5, 33)
    v = h-mean*one
    lower, upper = (25-sp.sqrt(3))/30, (25+sp.sqrt(3))/30
    minus, plus = sp.Rational(1, 7260)-sp.sqrt(3)/16335, sp.Rational(1, 7260)+sp.sqrt(3)/16335
    if p.eigenvals() != {sp.Integer(1): 1, lower: 1, upper: 1}:
        raise ArithmeticError('accessible transfer poles changed')
    if (prior*v)[0] != 0 or not all(sp.simplify(x) == 0 for x in
            (p-lower*sp.eye(3))*(p-upper*sp.eye(3))*v):
        raise ArithmeticError('centered current does not obey the exact two-pole recurrence')
    if not (minus > 0 and plus > 0):
        raise ArithmeticError('current spectral weights must be strictly positive')
    for n in range(9):
        if sp.simplify((v.T*d*p**n*v)[0]-minus*lower**n-plus*upper**n) != 0:
            raise ArithmeticError('native current covariance is not the two-pole spectral law')
    hankel = sp.Matrix(3, 3, lambda i, j: (prior*marks[1]**(i+j)*one)[0])
    if hankel.det() != sp.Rational(1, 9702990000):
        raise ArithmeticError('the three-dimensional record realization lost minimal rank')
    # Exact reversible detailed balance is the positive inner product behind
    # the symmetric transfer S=sqrt(D) P sqrt(D)^(-1), not a new stochastic law.
    if d*p != p.T*d or any(d*m != m.T*d for m in marks):
        raise ArithmeticError('marked current reflection symmetry failed')
    return {'classes': ['adjacent occupied ports', 'distance two occupied ports', 'opposite occupied ports'],
        'class_stationary_law': list(map(str, reduced['stationary'])),
        'exact_transfer_matrix': [[str(x) for x in row] for row in reduced['matrix']],
        'exact_marked_matrices': [[[str(x) for x in row] for row in m] for m in reduced['marked']],
        'symmetric_transfer': 'sqrt(D) P sqrt(D)^(-1), D=diag(5/11,5/11,1/11)',
        'nonconstant_eigenvalues': ['(25-sqrt(3))/30', '(25+sqrt(3))/30'],
        'positive_spectral_weights': ['1/7260-sqrt(3)/16335', '1/7260+sqrt(3)/16335'],
        'bit_covariance_with_n_unread_attempts': 'w_minus*r_minus^n + w_plus*r_plus^n',
        'minimal_linear_record_dimension': 3, 'all_one_word_Hankel_determinant': str(hankel.det()),
        'minimality_scope': 'stationary absolute-current words; not full native coherent outputs',
        'marked_insertions_are_POVM_effects': False}


def record_decay_sanity():
    """Recover dimensionless relaxation poles from finite current records alone.

    The first and last moving bits are retained; the intervening native
    attempts are unread, including zero-current outcomes. This is a declared
    uniform-scheduler protocol, not the fixed two-event E/F controller or a
    particle rest-energy measurement.
    """
    reduced = clock.lumping(source_kernel())
    transfer = sp.Matrix(reduced['matrix'])
    moving = sp.Matrix(reduced['marked'][1])
    prior = sp.Matrix([reduced['stationary']])
    one = sp.ones(3, 1)
    mean = (prior*moving*one)[0]
    covariance = tuple(sp.factor((prior*moving*transfer**n*moving*one)[0]-mean**2)
                       for n in range(4))
    c0, c1, c2, c3 = covariance
    hankel = sp.factor(c0*c2-c1*c1)
    if hankel <= 0:
        raise ArithmeticError('two accessible current decay modes are not identified')
    pole_sum = sp.factor((c0*c3-c1*c2)/hankel)
    pole_product = sp.factor((c1*c3-c2*c2)/hankel)
    discriminant = sp.factor(pole_sum**2-4*pole_product)
    if discriminant <= 0:
        raise ArithmeticError('the record decay poles are not distinct and real')
    slow, fast = ((pole_sum+sp.sqrt(discriminant))/2,
                  (pole_sum-sp.sqrt(discriminant))/2)
    if not 0 < fast < slow < 1:
        raise ArithmeticError('the two recovered record poles are not decays')
    # This protocol-conditional readout is dimensionless: no supplied seconds, energy
    # constant, physical mass operator, or particle label enters this ratio.
    ratio = sp.log(1/fast)/sp.log(1/slow)

    kernel = source_kernel()
    index = {state: i for i, state in enumerate(kernel.states)}
    control = []
    for occupied in ((4, 7), (4, 8)):
        state = tuple(int(i in occupied) for i in range(12))
        source_class = reduced['classes'][index[state]]
        prepared = tuple(F(int(i == source_class)) for i in range(3))
        one_bit = sum(clock.word_law(reduced, (1,), prepared))
        two_bits = sum(clock.word_law(reduced, (1, 1), prepared))
        control.append({'occupied_ports': list(occupied),
                        'one_moving_bit_probability': str(one_bit),
                        'two_moving_bits_probability': str(two_bits)})
    if (control[0]['one_moving_bit_probability'] != control[1]['one_moving_bit_probability']
            or control[0]['two_moving_bits_probability'] == control[1]['two_moving_bits_probability']):
        raise ArithmeticError('the wider native record protocol lost its separation control')
    return {'parameter': 'unread native attempt count between two retained moving bits',
        'preparation_and_schedule': 'stationary binary total-two law; independent uniform thirty-seam choices',
        'single_moving_bit_probability': str(mean),
        'four_centered_lag_covariances': list(map(str, covariance)),
        'record_Hankel_determinant': str(hankel),
        'recovered_pole_sum': str(pole_sum),
        'recovered_pole_product': str(pole_product),
        'recovered_slow_and_fast_poles': list(map(str, (sp.simplify(slow), sp.simplify(fast)))),
        'dimensionless_log_decay_ratio_fast_over_slow': str(ratio),
        'floating_ratio_diagnostic': float(ratio.evalf(17)),
        'prepared_state_record_control': control,
        'fixed_two_event_E_F_controller_supplies_this_long_lag_law': False,
        'source_selects_uniform_scheduler_or_stationary_preparation': False,
        'particle_mass_ratio_or_physical_energy_identified': False}


def charged_coherence_scope():
    """Separate the fixed-sector record poles from a full-algebra charge block.

    The latter is an exact native CP coherence decay. Its two-point expression
    requires a vacuum/one-token coherence and an analyzer; neither is among
    the retained current records, and a CP decay is not a rest energy.
    """
    kernel = source_kernel()
    states = set(kernel.states)
    if len(states) != 66 or {sum(state) for state in states} != {2}:
        raise ArithmeticError('the native current transfer left its total-two sector')
    projected_ladder_entries = 0
    for state in states:
        for port in range(12):
            changed = list(state)
            changed[port] = 1-changed[port]
            projected_ladder_entries += tuple(changed) in states
    if projected_ladder_entries:
        raise ArithmeticError('a number-changing ladder remained inside the fixed sector')

    for endpoint in range(4):
        for current in CURRENT_VALUES:
            if quantum.native_unit('C', endpoint, endpoint, (0, 1), current) != \
                    quantum.native_unit('E', endpoint, endpoint, (0, 1), current):
                raise ArithmeticError('classical endpoint/current records distinguish C from E')

    # The full native C/E channels also act on cross-number matrix units.
    # Check every port directly from the endpoint rule, without constructing
    # a fitted Hamiltonian or assigning a particle species to a mode.
    rates = {}
    for law in ('C', 'E'):
        expected = F(11, 12) if law == 'C' else F(5, 6)
        for port in range(12):
            row = {}
            for seam in kernel.seams:
                for key, weight in quantum.native_unit(law, 1 << port, 0, seam).items():
                    row[key] = row.get(key, F(0))+weight/F(len(kernel.seams))
            if row != {(1 << port, 0): expected}:
                raise ArithmeticError('the exact charge-one coherence block changed')
        rates[law] = str(expected)
    return {'record_transfer_sector': 'binary total-two; N=2 I_66',
        'projected_single_port_ladder_nonzero_entries': projected_ladder_entries,
        'all_operators_within_record_transfer_neutral_under_occupation_U1': True,
        'full_native_CP_vacuum_one_token_coherence_eigenvalues': rates,
        'degenerate_port_modes_per_coherent_law': 12,
        'classical_current_records_select_C_or_E_coherent_law': False,
        'formal_charge_one_two_point': 'Tr(b_i Phi^n(b_i^dagger |vac><vac|))=r^n',
        'two_point_parameter': 'unread native attempt count n',
        'vacuum_one_token_coherence_and_analyzer_executed_by_current_observer': False,
        'hard_core_ladder_is_selected_charged_lepton_CAR_operator': False,
        'CP_coherence_decay_identified_with_particle_rest_energy': False}


def run():
    kernel = source_kernel()
    adjacency = {i: {b if a == i else a for a, b in kernel.seams if i in (a, b)} for i in range(12)}
    neighbor = min(adjacency[0])
    opposite, = (i for i in range(1, 12) if i not in adjacency[0] and not adjacency[0]&adjacency[i])
    probes = ((0, 0, 0, 0), (0, neighbor, 0, neighbor), (0, opposite, 0, opposite))
    predictions = [current_moments(ports) for ports in probes]
    if any(F(row['connected_four_point']) == 0 for row in predictions):
        raise ArithmeticError('a native non-Gaussian four-point witness vanished')
    return {'schema': 'oph.native_effective_fields.v1', 'all_pass': True,
        'source_sha256': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS},
        'local_source_field': local_field_audit(), 'controlled_logarithm': logarithm_audit(),
        'finite_spatial_support': spatial_locality_audit(),
        'observer_record_spectral_reduction': record_spectral_reduction(),
        'observer_record_decay_sanity': record_decay_sanity(),
        'charged_coherence_scope': charged_coherence_scope(),
        'charged_reference_readout': charged_reference.reference_readout_scope(kernel.seams),
        'charged_spatial_pole_gate': charged_reference.charged_spatial_pole_gate(kernel.seams),
        'charged_reference_source_word': charged_reference_source.source_word_readout(kernel.seams),
        'four_current_predictions': predictions,
        'complete_executions': [execute_currents(ports, checkpoint=(i == 0)) for i, ports in enumerate(probes)],
        'scope': {'one_native_carrier_no_independent_replica_limit': True,
            'uniform_seam_scheduler_fair_coin_and_stationary_preparation_supplied': True,
            'field_generator_coefficients_follow_native_counts': True,
            'exact_logarithm_is_strictly_graph_local': False,
            'hidden_occupation_field_readouts_granted': False,
            'nonGaussianity_alone_establishes_quantum_interaction': False,
            'native_CP_evolution_identified_with_unitary_continuation': False,
            'spacetime_continuum_physical_dimension_or_relativistic_QFT_derived': False,
            'new_empirical_observations': False,
            'floating_propagation_or_physical_service_errors_certified': False}}
