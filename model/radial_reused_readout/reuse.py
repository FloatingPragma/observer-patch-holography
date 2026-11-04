"""Exact complete histories, native mean evolution and finite classification."""
from fractions import Fraction as F
from hashlib import sha256
from math import comb, factorial

import numpy as np

from ..matter_geometry_response import source as scheduler
from ..observer import native
from ..radial_geometry_join import bounds as parent_bounds, source


def count(value, maximum=1024):
    if type(value) is not int or not 0 <= value <= maximum:
        raise ValueError('a bounded nonnegative integer round count is required')
    return value


def history_coefficients(rounds, ones):
    """One ordered history; its quantum cross term is retained explicitly."""
    rounds, ones = count(rounds), count(ones)
    if ones > rounds:
        raise ValueError('pointer count exceeds the complete history length')
    denominator = 50**rounds
    return {'parity_zero': F(49**(rounds-ones), denominator),
            'parity_one': F(49**ones, denominator),
            'cross_parity': F(7**rounds, denominator),
            'ordered_histories_with_this_count': comb(rounds, ones)}


def history_probability(q, rounds, ones):
    if type(q) is not F or not 0 <= q <= 1:
        raise ValueError('an exact parity expectation in [0,1] is required')
    row = history_coefficients(rounds, ones)
    return (1-q)*row['parity_zero']+q*row['parity_one']


def unread_coefficients(rounds):
    rounds = count(rounds)
    diagonal = [sum((comb(rounds, k)*history_coefficients(rounds, k)[key]
                     for k in range(rounds+1)), F())
                for key in ('parity_zero', 'parity_one')]
    cross = sum((comb(rounds, k)*history_coefficients(rounds, k)['cross_parity']
                 for k in range(rounds+1)), F())
    if diagonal != [1, 1] or cross != F(7, 25)**rounds:
        raise ArithmeticError('complete repeated instrument lost a parity block')
    return {'diagonal_blocks': diagonal, 'parity_cross_multiplier': cross,
            'origin_hopping_energy_injection_coefficient': -(1-cross)}


def matmul(a, b):
    return [[sum((x*y for x, y in zip(row, col, strict=True)), F())
             for col in zip(*b, strict=True)] for row in a]


def matrix_power(a, rounds):
    rounds = count(rounds)
    result = [[F(i == j) for j in range(len(a))] for i in range(len(a))]
    while rounds:
        if rounds & 1:
            result = matmul(result, a)
        a = matmul(a, a)
        rounds //= 2
    return result


def reused_mean(q, rounds):
    if type(q) is not F or not 0 <= q <= 1:
        raise ValueError('an exact parity expectation in [0,1] is required')
    count(rounds)
    zero = matrix_power(source.mean_matrix(F(0)), rounds)
    one = matrix_power(source.mean_matrix(F(1)), rounds)
    return [[(1-q)*a+q*b for a, b in zip(left, right, strict=True)]
            for left, right in zip(zero, one, strict=True)]


def native_transition_control():
    """Independent complete native transition enumeration on all 78 inputs."""
    states = native.states()
    index = {state: i for i, state in enumerate(states)}
    transfers = [np.zeros((len(states), len(states)), dtype=np.int64) for _ in (0, 1)]
    digest = sha256()
    histories = 0
    for initial in states:
        for bit in (0, 1):
            for word in range(256):
                edge = scheduler.selected_seam(bit, word)
                for coin in (0, 1):
                    if edge is None:
                        output, current = initial, 0
                    else:
                        output, (_, current) = native.native(initial, edge, coin)
                    recovered = list(output)
                    if edge is not None:
                        recovered[edge[0]] += current
                        recovered[edge[1]] -= current
                    if tuple(recovered) != initial:
                        raise ArithmeticError('complete native current record is not reversible')
                    for parity in (0, 1):
                        weight = 49 if bit == parity else 1
                        transfers[parity][index[output], index[initial]] += weight
                    digest.update(repr((initial, bit, word, coin, edge, output, current)).encode())
                    histories += 1
    denominator = 50*512
    loads = np.array(states, dtype=np.int64).T
    for parity, transfer in enumerate(transfers):
        if not np.all(transfer.sum(axis=0) == denominator):
            raise ArithmeticError('idle or zero-current histories were deleted')
        actual = loads@transfer
        matrix = source.mean_matrix(F(parity))
        for j, state in enumerate(states):
            for i in range(12):
                expected = sum((matrix[i][k]*state[k] for k in range(12)), F())
                if F(int(actual[i, j]), denominator) != expected:
                    raise ArithmeticError('native conditional mean differs from its endpoint source')
    # A repeated-field source has one latent parity. Exact two-step propagation
    # on the complete 78-state carrier differs from fresh independent parity.
    mixture2 = loads@(transfers[0]@transfers[0]+transfers[1]@transfers[1])
    exact = reused_mean(F(1, 2), 2)
    independent = matrix_power(source.mean_matrix(F(1, 2)), 2)
    for j, state in enumerate(states):
        for i in range(12):
            expected = sum((exact[i][k]*state[k] for k in range(12)), F())
            if F(int(mixture2[i, j]), 2*denominator**2) != expected:
                raise ArithmeticError('complete two-step native mixture is incorrect')
    difference = exact[0][0]-independent[0][0]
    if difference <= 0:
        raise ArithmeticError('correlated reuse was silently replaced by fresh preparations')
    return {'complete_native_input_states': len(states),
            'complete_one_round_native_histories_checked': histories,
            'complete_histories_sha256': digest.hexdigest(),
            'exact_endpoint_transition_denominator': denominator,
            'all_initial_load_mean_identities_checked': True,
            'full_78_state_two_round_mean_checked': True,
            'two_round_reused_minus_fresh_matrix_00': difference,
            'reused_mean': '(1-q)*M_0^N+q*M_1^N',
            'fresh_independent_mean_substituted': False}


def classify(scores):
    if not isinstance(scores, (tuple, list)) or len(scores) != 32:
        raise ValueError('all 32 declared pointer scores are required')
    if any(type(bit) is not int or bit not in (0, 1) for bit in scores):
        raise ValueError('binary scores including retained failure scores are required')
    return 'charged' if sum(scores) >= 16 else 'ground'


def score_records(records):
    """Total score retaining failures; records are never removed or rewritten."""
    if not isinstance(records, (list, tuple)) or len(records) != 32:
        raise ValueError('all 32 complete record frames including a halted suffix are required')
    scores, halted = [], False
    for record in records:
        if type(record) is not dict:
            raise ValueError('a complete labeled record frame is required')
        kind = record.get('kind')
        if halted:
            if kind != 'halted':
                raise ValueError('only a declared halted suffix may follow a failure')
            scores.append(0)
        elif kind == 'pointer':
            value = record.get('value')
            if type(value) is not int or not 0 <= value < 6:
                raise ValueError('the complete six-outcome physical pointer is required')
            scores.append(int(value == 1))
            halted = value not in (0, 1)
        elif kind == 'failure' and type(record.get('native_label')) is str and record['native_label']:
            scores.append(0)
            halted = True
        else:
            raise ValueError('invalid pointer, native failure, or halted frame')
    return tuple(scores)


def classification_bound():
    rounds = 32
    exponent = 2*rounds*F(12, 25)**2
    exp_lower = sum((exponent**j/factorial(j) for j in range(21)), F())
    if exponent != F(9216, 625) or exp_lower <= 10**6:
        raise ArithmeticError('rational exponential tail certificate failed')
    tails = [sum((F(comb(rounds, k)*49**(rounds-k), 50**rounds)
                  for k in range(start, rounds+1)), F()) for start in (16, 17)]
    if max(tails) >= F(1, 10**6):
        raise ArithmeticError('exact binomial countercheck exceeds the stated tail')
    parent = parent_bounds.joined_bounds()
    _, word_bound = source.word_and_bound()
    attempt_error = F(word_bound['nine_coin_attempt_trace_distance_upper'])
    if attempt_error >= F(1, 10**9):
        raise ArithmeticError('complete parent control guarantee is absent')
    e = parent['ground_D0_expectation_interval'][1]
    total = e+F(1, 10**6)+rounds*F(1, 10**9)
    if total >= F(1, 2500):
        raise ArithmeticError('single-preparation classifier is not certified')
    return {'field_preparations': 1, 'rounds': rounds, 'tie_rule': 'charged at count >=16',
            'promised_input': 'actual G/2 or either unique lowest neighboring charge state',
            'conditional_record1_probabilities': [F(1, 50), F(49, 50)],
            'Hoeffding_exponent': exponent, 'exp_lower_terms_included': 21,
            'exp_lower_exceeds_one_million': True,
            'conditional_wrong_parity_probability_strict_upper': F(1, 10**6),
            'exact_ground_parity_conditional_error': tails[0],
            'exact_charged_parity_conditional_error': tails[1],
            'actual_spectral_input_wrong_parity_upper': e,
            'one_round_complete_half_trace_error_upper': attempt_error,
            'full_32_round_half_trace_error_upper': rounds*attempt_error,
            'each_spectral_input_misclassification_strict_upper': total,
            'simple_misclassification_strict_upper': F(1, 2500),
            'failure_protocol': 'retain first failure and prefix; score it and all unused slots zero',
            'spectral_sector_preservation_by_measurement_assumed': False}


def resources():
    word, bounds = source.word_and_bound()
    # The native fair-coin word has one preparation pulse and a three-pulse
    # copied pointer. Every round uses eight scheduler bits and one placement bit.
    return {'field_preparations': 1, 'rounds': 32,
            'native_parity_pulses_per_round': len(word),
            'native_fair_coin_pulses_per_coin': 4,
            'fresh_independent_fair_coins_per_round': 9,
            'native_pulses_total_upper': 32*(len(word)+9*4),
            'native_integer_repair_attempts_upper': 32,
            'fresh_probe_work_bits_per_round': 5,
            'fresh_two_bit_coin_registers_per_round': 9,
            'fresh_apparatus_logical_bits_without_recycling': 32*(5+9*2),
            'initial_complete_field_Fock_modes': 54,
            'native_work_states': 78,
            'retained_classical_record_frames': 32,
            'frame_fields': ['physical pointer or first-failure label', 'scheduler word',
                            'placement coin', 'seam or idle', 'current', 'terminal native loads'],
            'first_failure_label': 'round, pulse offset, native complement slot; padded halted suffix',
            'all_nonlogical_pointer_outputs_scored_zero_and_retained': True,
            'source_word_sha256': bounds['word_sha256'],
            'complete_32_round_native_word_executed': False,
            'physical_apparatus_preparation_and_drift_compensation_are_inputs': True}
