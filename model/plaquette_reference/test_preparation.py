"""Independent full-domain action, channel, source and modular regressions."""
from fractions import Fraction as F
import math

import mpmath as mp
import numpy as np
import pytest
from scipy.linalg import expm

from ..field_symmetry.symmetry import onsite_frame
from ..full_fock_field.field import FullFock
from ..observer.source_scalar import compiler
from . import preparation as p
from .packet import regional_comparison, resource_report


@pytest.fixture(scope='module')
def prepared():
    return p.prepare()


def test_exact_polynomial_matches_original_complete_action():
    data = p.design()
    original = FullFock(8).hamiltonian(p.GEOMETRY).toarray()
    assert np.array_equal(data.hamiltonian, original)
    x = -(original-np.eye(256))/64
    value, term = np.eye(256, dtype=complex), np.eye(256, dtype=complex)
    for order in range(1, 9):
        term = term@x/order
        value += term
    routed = data.frame.conj().T@(np.array(data.numerator, float)/data.denominator)@data.frame
    assert np.linalg.norm(value-routed, 2) < 1e-15
    assert np.linalg.norm(original@routed-routed@original, 2) < 1e-15
    assert len(data.column_norm_squared) == 256
    assert all(x > 0 for x in data.column_norm_squared)
    assert np.linalg.matrix_rank(routed) == 256


def test_rational_tail_and_faithfulness_without_float_cancellation():
    result = p.error_theorem()
    with mp.workdps(80):
        q = mp.mpf(3)/64
        tail = mp.exp(q)-sum(q**n/math.factorial(n) for n in range(9))
        bound = result['exponential_operator_error']
        assert 0 < tail < mp.mpf(bound.numerator)/bound.denominator
    assert result['regional_modular_operator_error'] < F(3, 10**16)
    assert result['regional_minimum_weight'] > F(1, 20)
    assert result['full_minimum_weight'] > F(1, 400)


@pytest.mark.parametrize('probability', [False, True])
def test_integer_tree_scaling_and_host_rounding_independence(monkeypatch, probability):
    values = [(j % 7)-3 for j in range(128)]
    if probability:
        values = [value*value for value in values]
    expected = p.integer_amplitude_word(values, probability=probability)
    def forbidden(*args, **kwargs):
        raise AssertionError('an exact-data tree must not use host norm/atan2')
    with monkeypatch.context() as context:
        context.setattr(np.linalg, 'norm', forbidden)
        context.setattr(math, 'atan2', forbidden)
        with mp.workdps(5):
            actual = p.integer_amplitude_word([value*(10**100+7) for value in values],
                                              probability=probability)
    assert actual == expected
    assert len(actual) == 3*(6*128+1)


@pytest.mark.parametrize('probability', [False, True])
def test_integer_tree_retains_signed_target_and_all_zero_subtrees(probability):
    from ..certified_reference import certificate
    values = [0]*128
    values[2], values[99] = (9, 16) if probability else (-3, 4)
    word = p.integer_amplitude_word(values, probability=probability)
    proof = certificate.amplitude_certificate(word, values, probability=probability, bits=7)
    assert proof['offset'] == len(word)
    assert proof['checked_conditional_nodes'] == 127
    assert proof['tree_error'] < F(1, 10**13)
    state = np.zeros((128, 1), complex)
    state[0] = 1
    paulis = {'X': np.array([[0, 1], [1, 0]]),
              'Y': np.array([[0, -1j], [1j, 0]]),
              'Z': np.diag([1, -1]), 'I': np.eye(2)}
    paulis['H'] = (paulis['X']+paulis['Z'])/math.sqrt(2)
    for event in word:
        kind, angle = event['kind'], event['angle']
        gate = (np.diag([np.exp(1j*angle)]*3+[1]) if kind == 'J' else
                math.cos(angle)*np.eye(2)+1j*math.sin(angle)*paulis[kind])
        state = compiler._apply(state, gate, event['targets'], 7)
    expected = np.sqrt(values) if probability else np.array(values)
    expected = expected/5
    assert np.linalg.norm(state[:, 0]-expected) < 1e-12
    assert any(event['angle'] == 0 for event in word)


@pytest.mark.parametrize('values,probability', [
    ([0]*128, False), ([1]*127, False), ([1.0]*128, False),
    ([True]*128, False), ([-1]*128, True), ([1]*128, 1),
])
def test_integer_tree_rejects_changed_exact_data_domain(values, probability):
    with pytest.raises(ValueError):
        p.integer_amplitude_word(values, probability=probability)


@pytest.mark.parametrize('site', p.GEOMETRY.sites)
def test_native_onsite_frame_all_four_occupations(site):
    actual = compiler.execute(p.onsite_word(site), 2)
    expected = np.eye(4, dtype=complex)
    expected[np.ix_([1, 2], [1, 2])] = onsite_frame(site)
    assert np.linalg.norm(actual-expected, 2) < 1e-13


def test_inverse_frame_complete_fock_columns():
    actual = compiler.execute(p.inverse_frame_word(), 8)
    assert np.linalg.norm(actual-p.design().frame.conj().T, 2) < 1e-13


def test_source_readout_all_complex_inputs_and_abort_vectors():
    rng = np.random.default_rng(193)
    vector = rng.normal(size=256)+1j*rng.normal(size=256)
    vector /= np.linalg.norm(vector)
    read = p.record_readout(vector)
    assert np.array_equal(read.columns, np.diag(vector))
    assert read.events == 255
    assert len(read.failure_weights) == 255
    assert sum(map(sum, read.failure_weights)) == 0
    for event in (0, 2, 16, 254):
        aborts = read.failure(event)
        assert len(aborts) == 4
        assert all(a.shape == (1, 128) and np.count_nonzero(a) == 0 for a in aborts)
    with pytest.raises(ValueError):
        read.initial.setflags(write=True)
    with pytest.raises(ValueError):
        read.failure(255)


def test_readout_missing_output_fails_closed(monkeypatch):
    original = compiler.transfer_kraus
    monkeypatch.setattr(compiler, 'transfer_kraus', lambda *a, **kw: original(*a, **kw)[:-1])
    with pytest.raises(ArithmeticError, match='lost an outcome'):
        p.record_readout(np.eye(256)[0])


@pytest.mark.parametrize('bad', [True, -1, 256, 1.2, '1'])
def test_conditional_record_domain_is_complete_and_bounded(bad):
    with pytest.raises(ValueError):
        p.conditional_word(bad)
    with pytest.raises(ValueError):
        p.design().column(bad)


@pytest.mark.parametrize('bad', [np.zeros(255), np.ones(256)*np.nan, np.ma.array(np.zeros(256))])
def test_invalid_readout_inputs_fail_closed(bad):
    with pytest.raises(ValueError):
        p.record_readout(bad)


def test_all_executed_records_and_failures_accounted(prepared):
    assert len(prepared.conditional_runs) == 256
    assert prepared.columns.shape == (256, 256)
    assert prepared.complete_output_error < 2e-11
    assert prepared.failure_probability >= 0
    assert abs(np.vdot(prepared.columns, prepared.columns).real+prepared.failure_probability-1) < 2e-11
    norms = np.sum(abs(prepared.columns)**2, axis=0)
    assert min(norms) > 0
    assert np.max(abs(norms-p.design().probabilities())) < 1e-13
    for index in (0, 85, 255):
        run = prepared.conditional_runs[index]
        assert np.count_nonzero(run.initial) == 1
        assert abs(run.initial[index, 0]) > 0
        event = len(run.word)-1
        actual_failure = sum(np.vdot(a, a).real for a in run.failure(event))
        assert abs(actual_failure-run.failure_weights[event].sum()) < 1e-35


def test_unconditioned_field_is_full_rank_gibbs_approximation(prepared):
    h = p.design().hamiltonian
    raw = expm(-np.asarray(h)/32)
    target = raw/np.trace(raw)
    assert np.linalg.norm(prepared.field_density-target, 'nuc') < 1e-11
    assert np.linalg.eigvalsh(prepared.field_density)[0] > 1/400
    # A selected record is a pure column state; only the full marginal is Gibbs.
    branch = prepared.columns[:, 85]
    conditional = np.outer(branch, branch.conj())/np.vdot(branch, branch)
    assert np.linalg.norm(conditional-target, 'nuc') > 1
    assert np.linalg.norm(h@prepared.field_density-prepared.field_density@h, 2) < 1e-13


def test_source_resource_count_and_regional_modular_attachment(prepared):
    resources = resource_report(prepared)
    assert resources['probability_preparation_pulses'] == 5379
    assert resources['conditional_pulses_min'] == 5468
    assert resources['conditional_pulses_max'] == 5476
    assert resources['trajectory_pulses_max'] == 10855
    assert resources['readout_first_failure_Kraus_slots'] == 1020
    result = regional_comparison(prepared)
    assert result['executed_normalized_P8_operator_error'] < 2e-9
    assert result['executed_vs_Gibbs_regional_modular_operator_diagnostic'] < 2e-11
    assert result['combined_P8_and_execution_diagnostic_bound'] < 2e-9
    assert F(result['conditional_admission_combined_P8_bound_exact']) < F(271, 10**11)
    assert F(result['conditional_admission_combined_P8_bound_exact']) > F(27, 10**10)
