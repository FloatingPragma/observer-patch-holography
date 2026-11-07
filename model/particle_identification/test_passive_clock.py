"""Source block coverage and independent small channel/Choi witnesses."""
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm

from .passive_clock import (native_mixing_upper, passive_error_lower,
                            source_contract_audit, unitary_error_lower)


def _bell():
    vector = np.array([1., 0., 0., 1.])/np.sqrt(2)
    return np.outer(vector, vector)


def _trace_distance(a, b):
    return np.abs(np.linalg.eigvalsh(a-b)).sum()/2


def _apply_first(kraus, state):
    matrices = [np.kron(k, np.eye(2)) for k in kraus]
    return sum(k@state@k.conj().T for k in matrices)


def test_source_audit_covers_every_quantum_direction_and_exact_extremal_root():
    audit = source_contract_audit()
    assert audit['carrier_dimension'] == 66
    assert audit['nonconstant_quantum_spectral_radius'] == '(55+sqrt(5))/60'
    assert audit['candidates']['C']['off_diagonal_spectral_radius_numerator_upper_over_60'] == 50
    assert audit['candidates']['E']['off_diagonal_spectral_radius_numerator_upper_over_60'] == 42
    for candidate, lower in (('C', 40), ('E', 22)):
        row = audit['candidates'][candidate]
        assert row['full_quantum_gershgorin_lower_numerator_over_60'] == lower
        assert row['quantum_block_size_census'] == {'1': 2970, '10': 132, '66': 1}
        assert row['complete_operator_directions'] == 66**2


def test_outward_rational_bound_and_tensor_cost():
    assert F(229, 240) > (55+np.sqrt(5))/60
    assert native_mixing_upper(0) == 1
    assert passive_error_lower(0, 2) == 0
    assert F(74, 100) < passive_error_lower(200, 2) < F(3, 4)
    assert native_mixing_upper(200, 2) == 2*native_mixing_upper(200)
    assert passive_error_lower(200, 2, 2) < passive_error_lower(200, 2)
    assert passive_error_lower(200, 4) > passive_error_lower(200, 2)
    assert passive_error_lower(300, 2) > passive_error_lower(200, 2)


@pytest.mark.parametrize('mixing', [F(), F(1, 10), F(2, 3)])
def test_independent_decoded_depolarizing_choi_obeys_the_unitary_lower_bound(mixing):
    # An arbitrary unitary encoder, noisy channel, and nonunital decoder.
    x = np.array([[0., 1.], [1., 0.]])
    z = np.diag([1., -1.])
    encoder = expm(-.41j*x)
    target = expm(.27j*z)
    gamma = .4
    decoder = (np.diag([1., np.sqrt(1-gamma)]),
               np.array([[0., np.sqrt(gamma)], [0., 0.]]))
    assert np.allclose(sum(k.conj().T@k for k in decoder), np.eye(2))
    encoded = _apply_first((encoder,), _bell())
    noisy = float(mixing)*encoded+(1-float(mixing))*np.eye(4)/4
    actual = _apply_first(decoder, noisy)
    desired = _apply_first((target,), _bell())
    decoded_replacement = _apply_first(decoder, np.eye(4)/4)
    assert _trace_distance(actual, decoded_replacement) <= float(mixing)+1e-14
    assert _trace_distance(actual, desired) >= float(unitary_error_lower(mixing, 2))-1e-14
    # The maximally entangled projector alone yields the universal 3/4
    # separation from every replacement state, including this nonunital one.
    assert abs(np.trace(desired@decoded_replacement)-.25) < 1e-14


def test_one_prepared_state_can_agree_while_the_all_input_unitary_bound_does_not():
    zero = np.diag([1., 0.])
    replacement_choi = np.kron(zero, np.eye(2)/2)
    replacement_kraus = (np.diag([1., 0.]), np.array([[0., 1.], [0., 0.]]))
    assert np.array_equal(sum(k@zero@k.conj().T for k in replacement_kraus), zero)
    assert _trace_distance(replacement_choi, _bell()) >= float(unitary_error_lower(0, 2))


def test_noiseless_bypass_and_retained_recovery_record_really_escape_the_claim():
    # Keeping the logical state in a bypass memory makes replacement of
    # an unrelated source carrier harmless. That decoder is not admitted.
    state = np.kron(np.eye(2)/2, _bell())
    recovered = np.trace(state.reshape(2, 4, 2, 4), axis1=0, axis2=2)
    assert np.allclose(recovered, _bell())
    # A discarded random-unitary record loses information that is recovered
    # exactly if the decoder receives and uses the record.
    identity = np.eye(2)
    flip = np.array([[0., 1.], [1., 0.]])
    raw = [_apply_first((k/np.sqrt(2),), _bell()) for k in (identity, flip)]
    unread = sum(raw)
    corrected = sum(_apply_first((k,), branch) for k, branch in zip((identity, flip), raw, strict=True))
    assert _trace_distance(unread, _bell()) > .49
    assert np.allclose(corrected, _bell())


@pytest.mark.parametrize('bad', [True, -1, .5, '10'])
def test_invalid_or_inexact_attempts_and_dimensions_are_rejected(bad):
    with pytest.raises(ValueError):
        native_mixing_upper(bad)
    with pytest.raises(ValueError):
        passive_error_lower(10, bad)


def test_invalid_copy_counts_and_mixing_bounds_are_rejected():
    with pytest.raises(ValueError):
        native_mixing_upper(10, 0)
    with pytest.raises(ValueError):
        unitary_error_lower(.1, 2)
    with pytest.raises(ValueError):
        unitary_error_lower(-1, 2)
