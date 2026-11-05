"""Independent operator, complete-output, record and negative controls."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm

from tools import strict_json
from ..field_thermal.thermal import car_restriction, restrict
from ..reference_recovery.recovery import execute_complete
from . import flow, packet


@pytest.fixture(scope='module')
def certificate():
    prior = strict_json.load(packet.ROOT/flow.entropy.PRIOR)
    return flow.exact_bounds(prior)


@pytest.fixture(scope='module')
def prepared():
    # This is the original complete native preparation, including all256
    # column records and every first failure, not a supplied Gibbs array.
    return flow.preparation.prepare()


def test_exact_Pauli_attachment_preserves_boundary_terms_and_all_states(certificate):
    p, terms = flow.generator()
    assert len(terms) == 39 and len(p) == 256
    assert any(label == 'IIII' and value for label, value in terms)
    bare = flow.expansion.plaquette()['coefficients'][0]/32
    assert p*bare-bare*p != flow.sp.zeros(16)
    assert flow.expansion.centered(p-bare) != flow.sp.zeros(16)
    assert certificate['resources']['native_pulses'] == 37568
    assert certificate['resources']['first_failure_Kraus_slots'] > 37568
    assert certificate['complete_half_diamond_error_upper'] < F(1, 10**8)
    assert certificate['ideal_symmetric_product_error'] > 0
    assert certificate['literal_complete_output_error'] > 0


def test_actual_native_complete_operator_matches_positive_sign_and_canonical_clock(certificate):
    run = flow.execute()
    p = np.asarray(flow.generator()[0], complex)
    target = expm(.5j*p)
    error = np.linalg.norm(run.success-target, 2)
    assert error < float(certificate['ideal_symmetric_product_error']+
                         certificate['literal_complete_output_error'])+1e-10
    assert np.linalg.norm(run.success-expm(-.5j*p), 2) > .01
    np.testing.assert_allclose(run.success.conj().T@run.success+run.failure_gram,
                               np.eye(16), atol=1e-10)
    assert run.initial.shape == (16, 16)
    assert len(run.word) == 37568 and len(run.failure_weights) == 37568
    # A complex observable distinguishes the Heisenberg sign convention.
    probe = np.zeros((16, 16), complex)
    probe[0, 3], probe[3, 0] = 1j, -1j
    actual = run.success.conj().T@probe@run.success
    expected = target.conj().T@probe@target
    assert np.linalg.norm(actual-expected, 2) < 1e-9
    assert np.linalg.norm(actual-target@probe@target.conj().T, 2) > 1e-3


def test_complex_untouched_reference_and_complete_failure_grams():
    rng = np.random.default_rng(2731)
    columns = rng.normal(size=(16, 3))+1j*rng.normal(size=(16, 3))
    columns /= np.linalg.norm(columns)
    run = flow.execute()
    fresh = execute_complete(flow.word(), 4, columns)
    np.testing.assert_allclose(fresh.success, run.success@columns, atol=2e-13)
    np.testing.assert_allclose(fresh.failure_gram, columns.conj().T@run.failure_gram@columns,
                               atol=1e-27)
    np.testing.assert_allclose(fresh.success.conj().T@fresh.success+fresh.failure_gram,
                               columns.conj().T@columns, atol=2e-12)
    # Explicit rows include every code complement at selected first failures.
    for index in (0, len(run.word)//2, len(run.word)-1):
        branches = run.failure(index)
        kind, _, _ = run.word[index]
        assert len(branches) == (2 if kind == 'J' else 4)
        for a, b in zip(fresh.failure(index), branches, strict=True):
            np.testing.assert_allclose(a, b@columns, atol=2e-26)


def test_original_preparation_continues_all_records_with_actual_regional_reference(prepared, certificate):
    continued = flow.continue_prepared(prepared)
    diagnostics = flow.prepared_diagnostic(prepared, continued)
    assert diagnostics['all256_classical_records_retained']
    assert diagnostics['success_renormalized'] is False
    assert diagnostics['actual_reference_modular_action_error_up_to_scalar'] < float(
        certificate['complete_half_diamond_error_upper'])+1e-10
    assert diagnostics['regional_reference_stationarity_error'] < 2e-9
    assert diagnostics['complete_trace_error'] < 2e-11
    assert continued.prepared is prepared
    assert prepared.readout.events == 255
    # The same256 columns retain their complete classical record index.
    for record in (0, 73, 255):
        expected = np.kron(np.eye(16), continued.regional_run.success)@prepared.columns[:, record]
        np.testing.assert_allclose(continued.columns[:, record], expected, atol=2e-17)
    assert continued.failure_probability == prepared.failure_probability+continued.actuator_failure_probability
    for event in (0, len(continued.regional_run.word)-1):
        for failure, local in zip(continued.failure(event), continued.regional_run.failure(event), strict=True):
            assert failure.shape[1] == 256
            np.testing.assert_array_equal(failure, np.kron(np.eye(16), local)@prepared.columns)
    # Exterior occupation records have the same ideal marginal. This does
    # not claim that the local action preserves all global correlations.
    maps = car_restriction(8, (4, 5, 6, 7))
    np.testing.assert_allclose(restrict(continued.field_density, maps),
                               restrict(prepared.field_density, maps), atol=2e-12)


def test_inverse_and_group_comparison_use_complete_source_words():
    positive = flow.execute()
    negative = flow.execute(F(-1, 2), 32)
    quarter = flow.execute(F(1, 4), 16)
    np.testing.assert_allclose(negative.success@positive.success, np.eye(16), atol=2e-10)
    np.testing.assert_allclose(quarter.success@quarter.success, positive.success, atol=2e-10)


@pytest.mark.parametrize('damage', ['removed', 'target', 'sign', 'extra'])
def test_uncompiled_routing_and_changed_source_control_fail(damage):
    word = flow.word()
    if damage == 'removed':
        word.pop(0)
    elif damage == 'target':
        word[0]['targets'] = [3]
    elif damage == 'sign':
        # The scalar P8 phase is nonzero and is also an accounted request.
        word[0]['angle'] *= -1
    else:
        word.append(deepcopy(word[0]))
    prior = strict_json.load(packet.ROOT/flow.entropy.PRIOR)
    with pytest.raises((ValueError, ArithmeticError)):
        flow.exact_bounds(prior, emitted=word)


def test_missing_abort_export_cannot_pass(monkeypatch):
    original = flow.compiler.transfer_kraus
    monkeypatch.setattr(flow.compiler, 'transfer_kraus',
                        lambda dimension, coherent=True: original(dimension, coherent)[:-1])
    prior = strict_json.load(packet.ROOT/flow.entropy.PRIOR)
    with pytest.raises(ArithmeticError, match='branch was dropped'):
        flow.exact_bounds(prior)


@pytest.mark.parametrize('parameter,steps', [(.5, 32), (True, 32), (F(3, 2), 32),
                                           (F(1, 2), True), (F(1, 2), 0), (F(1, 2), 65)])
def test_request_domain_fails_closed(parameter, steps):
    with pytest.raises(ValueError):
        flow.word(parameter, steps)
    # Validate before memoization: .5 and Fraction(1,2) compare equal.
    flow.execute(F(1, 2), 32)
    with pytest.raises(ValueError):
        flow.execute(parameter, steps)


def test_zero_parameter_retains_source_requests():
    word = flow.word(0, 1)
    assert len(word) == 1174
    run = execute_complete(word, 4, np.eye(16))
    np.testing.assert_allclose(run.success, np.eye(16), atol=2e-12)
    assert len(run.failure_weights) == 1174


def test_cold_receipt_and_scope_mutations():
    result = packet.report()
    packet.validate(strict_json.load(packet.OUTPUT), result)
    for flag in ('bare_regional_Gibbs_action_substituted',
                 'regional_action_identified_with_global_modular_restriction',
                 'physical_time_or_energy_unit_calibrated',
                 'Einstein_or_interacting_continuum_gate_accepted'):
        changed = deepcopy(result)
        changed['scope'][flag] = True
        with pytest.raises(ValueError, match='differs'):
            packet.validate(changed, result)
    changed = deepcopy(result)
    changed['exact_certificate']['complete_half_diamond_error_upper'] = '0'
    with pytest.raises(ValueError, match='differs'):
        packet.validate(changed, result)
