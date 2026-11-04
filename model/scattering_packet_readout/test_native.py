from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm

from . import native


@pytest.fixture(scope='module')
def controls():
    return {interacting: native.execute(interacting) for interacting in (True, False)}


def test_complete_source_echo_orientation_and_free_counterfactual(controls):
    actual, free = controls[True]['report'], controls[False]['report']
    assert actual['joint_quantum_rows'] == free['joint_quantum_rows'] == 32
    assert len(actual['pointer_probabilities']) == len(free['pointer_probabilities']) == 6
    assert .739 < actual['pointer_probabilities'][0] < .741
    assert free['pointer_probabilities'][0] == pytest.approx(.5, abs=1e-10)
    assert actual['complete_vector_error'] < .001
    for value in controls.values():
        report = value['report']
        assert sum(report['pointer_probabilities'])+report['failure_probability'] == pytest.approx(1, abs=1e-10)
        assert len(value['first_failure_events']) == len(value['word'])
        assert report['reference_gram_error'] < 1e-10
        assert report['analytic_scattering_phase_gate_used'] is False
    h0, h = native.hamiltonian(False), native.hamiltonian(True)
    wrong = expm(1j*h0/4)@expm(-1j*h/2)@expm(1j*h0/4)
    assert (1+wrong[3, 3].imag)/2 < .27
    assert abs((1+wrong[3, 3].imag)/2-actual['pointer_probabilities'][0]) > .47


def test_source_controls_retain_arbitrary_full_fock_reference_columns():
    rng = np.random.default_rng(117)
    state = rng.normal(size=(16, 2))+1j*rng.normal(size=(16, 2))
    state /= np.linalg.norm(state)
    result = native.execute(True, F(1, 8), 4, state)
    assert result['report']['reference_columns'] == 2
    assert result['report']['reference_gram_error'] < 1e-10
    assert result['report']['complete_vector_error'] < .002
    assert np.linalg.matrix_rank(state) == 2


def test_scalar_terms_and_controlled_free_identity_are_retained():
    assert any(label == 'IIII' for label, _ in native.source_terms(True))
    matrix = native.hamiltonian(True)
    stripped = matrix-np.trace(matrix)*np.eye(16)/16
    free = expm(-1j*native.hamiltonian(False)/4)
    wrong = free@expm(1j*stripped/2)@free
    assert abs(wrong[3, 3]-native.exact_echo(True)[3, 3]) > .2
    np.testing.assert_allclose(native.exact_echo(False), np.eye(16), atol=4e-16, rtol=0)


def test_literal_source_certificate_checks_every_fixed_angle_and_export(controls):
    for flag, result in controls.items():
        certificate = native.literal_certificate(result['word'], flag)
        assert certificate['native_pulses'] == len(result['word'])
        assert certificate['literal_control_complete_isometry_error_upper'] == F(1, 10**8)
        assert certificate['proper_code_exports_complete'] is True
    wrong = [dict(event) for event in controls[True]['word']]
    wrong[0]['angle'] = 0.
    with pytest.raises(ArithmeticError):
        native.literal_certificate(wrong, True)


def test_native_control_domains_are_explicit():
    for args in ((True, F(0), 8), (True, F(2), 8), (1, F(1, 4), 8)):
        with pytest.raises(ValueError):
            native.execute(*args)
    with pytest.raises(ValueError):
        native.controlled_evolution(True, .25, 8)
