from dataclasses import replace
from fractions import Fraction as F

import numpy as np
import pytest

from . import reference as r
from . import reversal as v


def _dense(row):
    matrix = np.zeros((132, 132))
    for i, j in row.cover_pairs():
        matrix[j, i] = np.sqrt(float(row.probability))
    return matrix


def _reference_roots(alpha):
    # Construct the roots independently in the orientation basis. Distinct
    # charge weights test coherences that a diagonal population check omits.
    pi = np.array([float(weight) for weight in r.local.reference(r.parameters())])
    plus = np.ones((2, 2))/2
    minus = np.array([[1, -1], [-1, 1]])/2
    root = np.sqrt(alpha)*plus+np.sqrt(1-alpha)*minus
    inverse = plus/np.sqrt(alpha)+minus/np.sqrt(1-alpha)
    return np.kron(root, np.diag(np.sqrt(pi))), np.kron(inverse, np.diag(1/np.sqrt(pi)))


def test_complete_weighted_reversal_retains_raw_zeros_holds_and_aborts():
    result = v.branch_certificate()
    assert result['raw_slots'] == 720
    assert result['cover_columns'] == 15840
    assert result['moving_slots'] == 120
    assert result['holding_slots'] == 180
    assert result['zero_slots'] == 420
    assert result['abort_slots'] == 240
    assert result['exp_modular_transfer_columns'] == {'1': 15200, '25/8': 320, '8/25': 320}
    assert result['complete_reverse_channel_equals_original']
    assert not result['zero_operator_reverse_unique']
    for row in r.branches():
        label = v.label(row)
        assert v.reverse_label(v.reverse_label(label)) == label
        if not row.inputs or row.inputs == row.outputs:
            assert v.reverse_label(label) == label


@pytest.mark.parametrize('alpha', [1/3, 1/2, 5/7])
def test_petz_operator_identity_on_noncommuting_input_and_entangled_reference(alpha):
    root, inverse = _reference_roots(alpha)
    rows = r.branches()
    index = {v.label(row): row for row in rows}
    rng = np.random.default_rng(1741)
    # A 132-by-3 amplitude matrix represents a pure state entangled with an
    # untouched three-dimensional reference, without a product-state input.
    entangled = rng.normal(size=(132, 3))+1j*rng.normal(size=(132, 3))
    entangled /= np.linalg.norm(entangled)
    tested = 0
    for row in rows:
        if row.seam not in ((0, 5), (0, 11)):
            continue
        actual = root@_dense(row).conj().T@inverse
        expected = _dense(index[v.reverse_label(v.label(row))])
        assert np.linalg.norm(actual-expected) < 5e-15
        assert np.linalg.norm(actual@entangled-expected@entangled) < 5e-15
        tested += 1
    assert tested == 48


@pytest.mark.parametrize('charge', [-1, 1])
def test_endpoint_reference_uses_its_invariant_charge_corner(charge):
    rows = r.branches()
    indexed = {v.label(row): row for row in rows}
    pi = np.array([float(weight) for weight in r.local.reference(r.parameters())])
    root, inverse = np.sqrt(pi), 1/np.sqrt(pi)
    for row in rows:
        matrix, reverse = np.zeros((66, 66)), np.zeros((66, 66))
        for i, j, sign in zip(row.inputs, row.outputs, row.signs, strict=True):
            matrix[j, i] = np.sqrt(float(row.probability))*(sign if charge == -1 else 1)
        other = indexed[v.reverse_label(v.label(row))]
        for i, j, sign in zip(other.inputs, other.outputs, other.signs, strict=True):
            reverse[j, i] = np.sqrt(float(other.probability))*(sign if charge == -1 else 1)
        assert np.max(np.abs(root[:, None]*matrix.T*inverse[None, :]-reverse)) < 1e-15


def test_full_heisenberg_kms_symmetry_including_cross_charge_coherences():
    root, _ = _reference_roots(1/3)
    rows = r.branches()
    rng = np.random.default_rng(817)
    a, b = [rng.normal(size=(132, 132))+1j*rng.normal(size=(132, 132)) for _ in range(2)]
    a, b = a/np.linalg.norm(a), b/np.linalg.norm(b)

    def heisenberg(operator):
        result = np.zeros((132, 132), dtype=complex)
        for row in rows:
            pairs = row.cover_pairs()
            if pairs:
                inputs, outputs = zip(*pairs, strict=True)
                result[np.ix_(inputs, inputs)] += float(row.probability)*operator[np.ix_(outputs, outputs)]
        return result

    inner = lambda first, second: np.trace(root@first.conj().T@root@second)
    assert abs(inner(a, heisenberg(b))-inner(heisenberg(a), b)) < 1e-15


def test_using_fair_coin_in_the_biased_source_breaks_weighted_reversal():
    rows = list(r.branches())
    index = next(i for i, row in enumerate(rows)
                 if v.label(row) == ('C', (0, 11), 0, 1))
    rows[index] = replace(rows[index], probability=F(1, 120))
    with pytest.raises(ArithmeticError, match='weighted Petz'):
        v.branch_certificate(rows)


def test_erasing_one_ordering_flip_breaks_the_full_quantum_adjoint():
    rows = list(r.branches())
    index = next(i for i, row in enumerate(rows)
                 if v.label(row) == ('C', (0, 11), 0, 1))
    row = rows[index]
    signs = list(row.signs)
    signs[signs.index(-1)] = 1
    rows[index] = replace(row, signs=tuple(signs))
    with pytest.raises(ArithmeticError, match='ordered-cover adjoint'):
        v.branch_certificate(rows)


def test_dropping_a_zero_abort_address_is_rejected():
    rows = list(r.branches())
    rows.pop(next(i for i, row in enumerate(rows) if row.abort))
    with pytest.raises(ArithmeticError, match='complete 720-slot'):
        v.branch_certificate(rows)


def test_odd_exchange_histories_and_all_two_event_outcomes_reverse():
    result = v.history_certificate()
    assert result['empty_history']['cover_rank'] == 132
    assert result['empty_history']['stationary_probability'] == '1'
    assert result['all_two_event_raw_histories_on_two_declared_seams'] == 576
    assert result['zero_two_event_histories_retained'] == 494
    for value in result['odd_exchange_loops'].values():
        assert value['cover_rank'] == 2
        assert value['stationary_probability'] == '1/88267806'
        assert value['orientation_action'] == 'X'
        assert value['modular_transfer'] == '0'
    assert not result['arbitrary_input_forward_reverse_record_probabilities_equal']
    assert not result['arbitrary_record_adaptive_policy_included']


def test_stationary_history_balance_does_not_equal_arbitrary_input_probabilities():
    word = (('C', (0, 11), 0, 1),)
    forward = v.history_operator(word)
    reverse = v.history_operator(v.reverse_history(word))
    initial = r.exchange.INDEX[(0, 5)]
    assert initial in dict(forward.pairs)
    assert initial not in dict(reverse.pairs)
    assert forward.probability > 0
    # Keeping the same nonstationary input gives positive forward and zero
    # reverse probability, despite the exact weighted operator identity.
    assert forward.probability*F(initial in dict(forward.pairs)) != reverse.probability*F(initial in dict(reverse.pairs))


def test_zero_history_keeps_later_label_validation_and_no_postselection():
    word = (('C', (0, 11), 0, 3), ('E', (0, 5), 1, 6))
    assert v.history_operator(word).pairs == ()
    assert v.history_operator(v.reverse_history(word)).pairs == ()
    assert v.reverse_history(v.reverse_history(word)) == word
    with pytest.raises(ValueError, match='record address'):
        v.history_operator((word[0], ('E', (0, 5), 1, 8)))


def test_petz_self_duality_cannot_be_used_as_an_inverse_on_all_states():
    result = v.contraction_control()
    assert result['input_half_trace_distance'] == '1'
    assert result['output_half_trace_distance'] == '949/990'
    assert not result['unread_channel_has_a_cptp_left_inverse_on_all_states']


@pytest.mark.parametrize('value', [[], ('C', (0, 11), 0), ('E', (0, 11), True, 6),
                                  ('C', (0, 11), 0, 4), ('E', (0, 11), 1, 8)])
def test_only_complete_typed_native_labels_are_accepted(value):
    with pytest.raises(ValueError):
        v.reverse_label(value)


def test_public_certificate_keeps_the_physical_boundary():
    result = v.certificate()
    assert result['scope']['same_weighted_132_state_source_and_raw_instrument']
    assert result['scope']['quantum_operator_identity_not_only_population_balance']
    assert result['scope']['alpha_endpoints_use_invariant_support_corner']
    assert not result['scope']['detailed_balance_selects_coherent_hamiltonian']
    assert not result['scope']['petz_reversal_is_physical_time_reversal']
    assert not result['scope']['physical_stress_split_or_interacting_qft_derived']
