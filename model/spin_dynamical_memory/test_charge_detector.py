from fractions import Fraction as F
import numpy as np
import pytest

from model.native_spin_observer import circuits, execute
from . import charge_detector as detector


@pytest.mark.parametrize('bits', [1, 2, 3, 4])
def test_literal_fourier_word_all_reference_columns(bits):
    n = 1 << bits
    word = detector.fourier(tuple(range(bits)))
    actual = execute.complete(circuits.emit(word), bits, np.eye(n)/np.sqrt(n))
    expected = np.exp(2j*np.pi*np.outer(np.arange(n), np.arange(n))/n)/n
    np.testing.assert_allclose(actual.success, expected, atol=3e-14)
    assert actual.failure_probability < 1e-25


def test_complete_native_charge_detector_and_resource_plan():
    result = detector.native_two_site_execution()
    assert result['native_pulses'] == 55
    assert result['native_first_failure_slots'] == 196
    assert result['reference_columns_dephased'] is False
    assert [r['ideal_Choi_probability'] for r in result['complete_class_records']] == ['1/2', '1/4', '0', '1/4']
    assert F(result['complete_trace_error_upper']) < F(1, 10**12)
    op, _ = detector.source_projectors(16)
    plan = detector.precision_plan(op['signs'], 16)
    assert plan['source_pulses'] == 222
    assert F(plan['complete_trace_error_upper']) <= F(1, 10**6)


def test_source_charge_class_measurement_preserves_within_class_coherence():
    op, projectors = detector.source_projectors(16)
    columns = [j for j, q in enumerate(op['charges']) if q == 1]
    u, v = op['walsh'][:, columns[0]]/16, op['walsh'][:, columns[1]]/16
    state = (u+1j*v)/np.sqrt(2)
    coherent = projectors[1]@state/256
    np.testing.assert_allclose(coherent, state, atol=1e-15)
    # Independent fine-string channel retains probabilities but removes this
    # relative phase; it is not the instrument certified above.
    fine_density = (np.outer(u, u)+np.outer(v, v))/2
    assert abs(np.trace(fine_density@fine_density)-.5) < 1e-14
    assert abs(np.vdot(coherent, fine_density@coherent)-.5) < 1e-14


def test_exact_cube_all_charge_classes_and_common_clock():
    op, projectors = detector.source_projectors(16)
    q = op['twice_charge']/2
    values, vectors = np.linalg.eigh(q)
    for r, numerator in enumerate(projectors):
        selected = vectors[:, [i for i, value in enumerate(values) if round(value) % 16 == r]]
        np.testing.assert_allclose(numerator/256, selected@selected.T, atol=2e-14)
    assert [int(np.trace(p))//256 for p in projectors] == [70, 56, 28, 8, 1, 0, 0, 0, 0, 0, 0, 0, 1, 8, 28, 56]


def test_alias_probability_uses_original_weight_and_all_unwanted_classes():
    # A complete finite charge distribution with both +/-1 aliases and
    # unrelated other outcomes. Filtering cannot increase its charge moment.
    charges = [-31, -15, -1, 0, 1, 17, 33]
    weights = [F(1, 10000), F(1, 2000), F(1, 10), F(1, 5), F(1, 5), F(1, 2000), F(1, 10000)]
    moment = sum(q*q*w for q, w in zip(charges, weights))
    for modulus in (4, 8, 16):
        wrong = sum(w for q, w in zip(charges, weights) if q not in (-1, 1) and q % modulus in (1, modulus-1))
        assert wrong <= moment/(modulus-1)**2
        for contraction in (F(1), F(1, 3)):
            assert contraction**2*wrong <= moment/(modulus-1)**2
    row = detector.alias_bound(F(1), 16, F(1, 2))
    assert F(row['absolute_original_weight_wrong_charge_probability_upper']) == F(2, 75)
    assert F(row['conditional_wrong_charge_probability_upper']) == F(4, 75)


def test_same_source_prelimit_moment_threshold():
    for amplitude in (F(1, 7), F(1), F(2000), F(-5)):
        threshold = min(F(1), 1/(2*abs(amplitude)))
        row = detector.prelimit_charge_bound(amplitude, threshold)
        assert F(row['charge_second_moment_upper']) < 7*amplitude*amplitude
        assert row['below_seven_amplitude_squared'] is True


def test_charge_derivative_bound_on_actual_selected_source():
    op = detector.charge.operators()
    _, vectors = np.linalg.eigh(op['H']/op['scale'])
    ground = vectors[:, 0]
    charge = op['twice_charge']/2
    indices = np.arange(256)
    f = np.array([1, -3, 2, 4, -2, 5, 1, -4])/16
    q3 = np.diag(sum(op['signs'][x]*f[x]*(1-2*((indices>>x)&1))/2 for x in range(8)))
    q2 = np.zeros((256, 256), complex)
    for x in range(8):
        q2[indices^(1<<x), indices] += 1j*f[x]*(1-2*((indices>>x)&1))/2
    np.testing.assert_allclose(1j*(charge@q3-q3@charge), q2, atol=1e-15)
    for amplitude in (.125, 1., 3.):
        created = np.sin(amplitude*np.diag(q3))*ground
        actual = np.linalg.norm(charge@created)
        bound = abs(amplitude)*np.linalg.norm(q2@ground)+amplitude**2*np.linalg.norm(q2@q3-q3@q2, 2)/2
        assert actual <= bound+1e-13


@pytest.mark.parametrize('signs,modulus', [((1, 1), 4), ((1, -1), 3), ((1, -1), 2), ((1, 0), 4)])
def test_unqualified_charge_programs_rejected(signs, modulus):
    with pytest.raises(ValueError):
        detector.template(signs, modulus)
