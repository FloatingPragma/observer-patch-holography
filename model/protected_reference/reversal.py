"""Weighted Petz reversal of the actual complete 132-state repair instrument.

Every nonzero reverse is an existing raw branch, with its original biased
coin probability. The full record involution also retains all zero/abort
slots. History balance uses the stationary state and the fixed selector.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product

from . import reference


SOURCE_FILES = (
    'model/protected_reference/reversal.py',
    'model/protected_reference/REVERSAL.md',
    'model/protected_reference/test_reversal.py',
)


def label(row):
    return row.law, row.seam, row.coin, row.raw


def reverse_label(value):
    """An involution of all 720 labels, including impossible raw records.

    Moving slots exchange coin and current. Every holding slot is its own
    reverse. Empty and abort slots are fixed, an explicit convention on
    zero operators rather than a claim of a unique operational reverse.
    """
    if type(value) is not tuple or len(value) != 4:
        raise ValueError('one complete native raw label is required')
    law, seam, coin, raw = reference.exchange._label(*value)
    moving = {
        ('C', 0, 1): (1, 2), ('C', 1, 2): (0, 1),
        ('E', 0, 5): (1, 6), ('E', 1, 6): (0, 5),
    }
    other_coin, other_raw = moving.get((law, coin, raw), (coin, raw))
    return law, seam, other_coin, other_raw


def reverse_history(history):
    if type(history) is not tuple:
        raise ValueError('a tuple of complete raw labels is required')
    return tuple(reverse_label(value) for value in reversed(history))


def _index(rows):
    result = {label(row): row for row in rows}
    if len(result) != len(rows):
        raise ArithmeticError('duplicate complete raw labels')
    return result


@dataclass(frozen=True)
class History:
    probability: F
    pairs: tuple[tuple[int, int], ...]


def _compose(history, indexed):
    mapping, weight = {i: i for i in range(132)}, F(1)
    for value in history:
        reverse_label(value)  # Validate the full label even after a zero slot.
        if value not in indexed:
            raise ValueError('the complete raw label is absent from the instrument')
        row = indexed[value]
        step = dict(row.cover_pairs())
        if len(step) != len(row.cover_pairs()) or len(set(step.values())) != len(step):
            raise ArithmeticError('a raw branch is not a partial injection')
        mapping = {i: step[j] for i, j in mapping.items() if j in step}
        weight *= row.probability
    return History(weight, tuple(sorted(mapping.items())))


def history_operator(history, rows=None):
    """Exact partial permutation and squared amplitude of a complete history.

    Its action is defined on the full 132-state carrier; no outcome is
    normalized or removed. The empty history is the identity operator.
    """
    if type(history) is not tuple:
        raise ValueError('a tuple of complete raw labels is required')
    indexed = _index(reference.branches() if rows is None else rows)
    return _compose(history, indexed)


def _history_balance(history, indexed, pi):
    forward = _compose(history, indexed)
    backward = _compose(reverse_history(history), indexed)
    if backward.pairs != tuple(sorted((j, i) for i, j in forward.pairs)):
        raise ArithmeticError('the reversed complete history is not the adjoint partial map')
    for i, j in forward.pairs:
        if forward.probability*pi[i % 66] != backward.probability*pi[j % 66]:
            raise ArithmeticError('the weighted complete-history identity failed')
    first = forward.probability*sum((pi[i % 66]/2 for i, _ in forward.pairs), F())
    second = backward.probability*sum((pi[i % 66]/2 for i, _ in backward.pairs), F())
    if first != second:
        raise ArithmeticError('stationary complete-history probabilities differ')
    return {'length': len(history), 'cover_rank': len(forward.pairs),
            'stationary_probability': str(first),
            'forward_squared_amplitude': str(forward.probability),
            'reverse_squared_amplitude': str(backward.probability)}


def branch_certificate(rows=None):
    """Exact weighted adjoint identities on all source and orientation columns.

    The rational identity w_h*pi_in = w_reverse*pi_out, together with the
    existing ordered-cover adjoint, proves equality of the Petz Kraus
    operators. The memory density commutes with every orientation factor,
    including for input coherences and an arbitrary untouched reference.
    """
    rows = tuple(reference.branches() if rows is None else rows)
    indexed, pi = _index(rows), reference.local.reference(reference.parameters())
    expected = {label(row) for row in reference.branches()}
    if set(indexed) != expected or len(rows) != 720:
        raise ArithmeticError('the complete 720-slot native instrument is required')
    counts, ratios = Counter(), Counter()
    effects, stationary = [F()]*132, [F()]*132
    for row in rows:
        value = label(row)
        reverse = reverse_label(value)
        if reverse_label(reverse) != value or reverse not in indexed:
            raise ArithmeticError('the raw-record reversal is not a complete involution')
        other = indexed[reverse]
        pairs, inverse = tuple(sorted(row.cover_pairs())), tuple(sorted(other.cover_pairs()))
        if len({i for i, _ in pairs}) != len(pairs) or len({j for _, j in pairs}) != len(pairs):
            raise ArithmeticError('a complete cover branch is not a partial injection')
        if inverse != tuple(sorted((j, i) for i, j in pairs)):
            raise ArithmeticError('the reverse raw branch lost the ordered-cover adjoint')
        if row.abort != other.abort or row.abort and pairs:
            raise ArithmeticError('abort-slot reversal or zero support changed')
        if row.probability <= 0 or other.probability <= 0:
            raise ArithmeticError('positive native selector weights are required')
        if not pairs:
            counts['zero_slots'] += 1
            counts['abort_slots'] += int(row.abort)
            if reverse != value:
                raise ArithmeticError('zero slots must retain the declared fixed-label convention')
        elif all(i == j for i, j in pairs):
            counts['holding_slots'] += 1
            if reverse != value:
                raise ArithmeticError('a holding raw branch is not fixed by reversal')
        else:
            counts['moving_slots'] += 1
        for i, j in pairs:
            if row.probability*pi[i % 66] != other.probability*pi[j % 66]:
                raise ArithmeticError('the weighted Petz branch identity failed')
            # Each branch intertwines the memory flip, so sigma_alpha's
            # square roots cancel for arbitrary faithful alpha.
            flip_i, flip_j = (i+66) % 132, (j+66) % 132
            if (flip_i, flip_j) not in pairs:
                raise ArithmeticError('the branch changed the protected memory factor')
            effects[i] += row.probability
            stationary[j] += row.probability*pi[i % 66]/2
            ratios[str(pi[i % 66]/pi[j % 66])] += 1
            counts['cover_columns'] += 1
    if effects != [F(1)]*132:
        raise ArithmeticError('the original complete instrument is not trace preserving')
    if stationary != [pi[i % 66]/2 for i in range(132)]:
        raise ArithmeticError('the original complete instrument lost its faithful reference')
    return {'raw_slots': len(rows), **dict(sorted(counts.items())),
            'exp_modular_transfer_columns': dict(sorted(ratios.items())),
            'moving_label_pairs': {'C': [[0, 1], [1, 2]], 'E': [[0, 5], [1, 6]]},
            'other_labels_fixed_including_zero_and_abort': True,
            'zero_operator_reverse_unique': False,
            'same_seam_and_law': True,
            'weighted_identity': 'w_h*pi_in = w_reverse*pi_out',
            'petz_kraus_identity': 'rho_alpha^(1/2) K_h^dagger rho_alpha^(-1/2) = K_reverse(h)',
            'holds_self_dual': True,
            'reverse_modular_transfer': '-Delta_h',
            'complete_reverse_channel_equals_original': True,
            'heisenberg_channel_kms_self_adjoint': True,
            'identity_extends_to_arbitrary_untouched_reference': True}


def history_certificate():
    rows = reference.branches()
    indexed, pi = _index(rows), reference.local.reference(reference.parameters())
    first = tuple(label(row) for row in rows if row.seam == (0, 5))
    second = tuple(label(row) for row in rows if row.seam == (0, 11))
    checked, zeros = 0, 0
    for values in product(first, second):
        result = _history_balance(values, indexed, pi)
        zeros += result['cover_rank'] == 0
        checked += 1
    loops = {}
    for law in ('C', 'E'):
        loop = reference.protection.exchange_loop(law)
        values = tuple((name, tuple(seam), coin, raw) for name, seam, coin, raw in loop['history'])
        loops[law] = _history_balance(values, indexed, pi)
        loops[law]['orientation_action'] = 'X'
        loops[law]['modular_transfer'] = '0'
    return {'empty_history': _history_balance((), indexed, pi),
            'all_two_event_raw_histories_on_two_declared_seams': checked,
            'zero_two_event_histories_retained': zeros,
            'odd_exchange_loops': loops,
            'all_finite_history_proof': 'iterate the exact Petz branch identity with reversed label order',
            'stationary_complete_history_probabilities_equal': True,
            'initial_state': 'rho_alpha; 0<alpha<1, or its invariant supported corner at alpha=0,1',
            'selector': 'fixed native state-independent seam/law/biased-coin selector',
            'arbitrary_input_forward_reverse_record_probabilities_equal': False,
            'arbitrary_record_adaptive_policy_included': False}


def contraction_control():
    """Two orthogonal inputs lose distinguishability in the unread channel.

    Their complete raw records are retained by the instrument. This control
    traces those records and shows that Petz self-duality is no inverse
    dynamics on arbitrary source states.
    """
    x = reference.exchange.INDEX[(0, 5)]+66
    y = reference.exchange.INDEX[(5, 11)]
    outputs = [reference.apply_exact({(index, index): F(1)}) for index in (x, y)]
    difference = {key: outputs[0].get(key, F())-outputs[1].get(key, F())
                  for key in outputs[0].keys() | outputs[1].keys()}
    if any(i != j for i, j in difference):
        raise ArithmeticError('a basis-state unread output gained a coherence')
    distance = sum((abs(value) for value in difference.values()), F())/2
    if not 0 < distance < 1:
        raise ArithmeticError('the unread distinguishability contraction witness changed')
    return {'input_half_trace_distance': '1', 'output_half_trace_distance': str(distance),
            'initial_cover_indices': [x, y],
            'unread_channel_has_a_cptp_left_inverse_on_all_states': False,
            'petz_self_duality_is_samplewise_or_unitary_rewind': False}


def certificate():
    return {'schema': 'oph.protected_reference.reversal.v1',
            'branches': branch_certificate(), 'histories': history_certificate(),
            'irreversibility_control': contraction_control(),
            'scope': {
                'same_weighted_132_state_source_and_raw_instrument': True,
                'full_reference_requires_0_lt_alpha_lt_1': True,
                'alpha_endpoints_use_invariant_support_corner': True,
                'quantum_operator_identity_not_only_population_balance': True,
                'protected_charge_preserved': True,
                'detailed_balance_selects_coherent_hamiltonian': False,
                'petz_reversal_is_physical_time_reversal': False,
                'modular_parameter_identified_with_native_attempt_or_laboratory_time': False,
                'physical_stress_split_or_interacting_qft_derived': False}}
