"""Independent source, finite-time and complete-output controls."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest
import sympy as s

from . import comparison as c, packet


@pytest.fixture(scope='module')
def source():
    return c.checked_source()


@pytest.fixture(scope='module')
def spectral(source):
    data, probe = source
    pi = np.array(data['pi'], float).ravel()
    p = np.array(data['P'], float)
    transfer = p*np.sqrt(pi)[None, :]/np.sqrt(pi)[:, None]
    eigenvalues, vectors = np.linalg.eigh(transfer)
    h = (vectors*(-np.log(eigenvalues)))@vectors.T
    return pi, h, np.array(probe['Y'], complex), eigenvalues, vectors


def test_exact_polynomial_derivative_encloses_independent_spectral_result(source, spectral):
    pi, h, y, _, _ = spectral
    actual = np.trace(-1j*(h@np.diag(pi)-np.diag(pi)@h)@y).real
    certificate = c.derivative_interval(*source)
    assert float(F(certificate['lower'])) < actual < float(F(certificate['upper']))
    assert actual > .007 and actual < .008


def test_finite_time_response_and_clock_orientation(source, spectral):
    pi, h, y, values, vectors = spectral
    certificate = c.finite_time_certificate(*source)
    means = []
    for time in (-.1, .1):
        u = (vectors*np.exp(1j*time*np.log(values)))@vectors.T
        means.append(np.trace(y@u@np.diag(pi)@u.conj().T).real)
    assert abs(sum(means)) < 1e-14
    assert means[1] > float(F(certificate['mean_lower'])) > 1/2000
    assert np.linalg.norm(h, 2) < float(c.NORM_UPPER)


def test_reference_is_load_bearing(spectral):
    _, h, y, values, vectors = spectral
    u = (vectors*np.exp(.1j*np.log(values)))@vectors.T
    uniform = np.eye(66)/66
    assert np.linalg.norm(u@uniform@u.conj().T-uniform, 'nuc') < 1e-13
    assert abs(np.trace(y@u@uniform@u.conj().T)) < 1e-13


def test_logarithm_tail_is_not_discarded(source):
    early = c.derivative_interval(*source, degree=1)
    assert F(early['lower']) < 0 < F(early['upper'])
    final = c.derivative_interval(*source)
    assert F(final['lower']) > c.DERIVATIVE_LOWER
    assert F(final['upper'])-F(final['lower']) > 2*F(final['derivative_tail_upper'])


@pytest.mark.parametrize('degree', [True, -1, 0, 1.0])
def test_original_integer_logarithm_degree_required(source, degree):
    with pytest.raises(ValueError):
        c.derivative_interval(*source, degree=degree)


@pytest.mark.parametrize('damage', ['population', 'chart', 'y_sign', 'pair'])
def test_changed_source_or_detector_dictionary_fails(source, damage):
    data, probe = deepcopy(source)
    if damage == 'population':
        data['P'] = data['P'].T
    elif damage == 'chart':
        data['states'] = tuple(reversed(data['states']))
    elif damage == 'y_sign':
        probe['Y'] = -probe['Y']
    else:
        probe['pairs'] = probe['pairs'][:-1]
    with pytest.raises(ArithmeticError):
        c.checked_source(data, probe)


def test_all_raw_effects_reproduce_the_score_without_odd_conditioning(source):
    _, probe = source
    total = s.zeros(66); score = s.zeros(66)
    for row in probe['branches']:
        v = c.observable.branch_operator(probe, row, 'Y')
        effect = s.Rational(row.weight.numerator, row.weight.denominator)*v.H*v
        total += effect
        score += (row.contrast or 0)*effect
    assert len(probe['branches']) == 25
    assert total == s.eye(66)
    assert score == probe['Y']
    assert probe['zero'].rank() == 46


def test_complete_outputs_with_reference_obey_separation(source, spectral):
    _, probe = source
    pi, _, _, values, vectors = spectral
    u = (vectors*np.exp(.1j*np.log(values)))@vectors.T
    reference = np.array([[.5, .2j], [-.2j, .5]])
    stationary = np.kron(np.diag(pi), reference)
    evolved = np.kron(u@np.diag(pi)@u.conj().T, reference)
    a = c.observable.ideal_complete_output(probe, stationary, 'Y')
    b = c.observable.ideal_complete_output(probe, evolved, 'Y')
    assert a.keys() == b.keys() and len(a) == 25
    separation = sum(np.linalg.norm(b[k]-a[k], 'nuc') for k in a)/2
    # Trace the working carrier as in an actual complete-record/reference readout.
    # Explicit diagonal blocks avoid using a partial-trace helper from the proof.
    def reduced(matrix):
        return sum((matrix[2*i:2*i+2, 2*i:2*i+2] for i in range(66)), np.zeros((2, 2), complex))
    record_separation = sum(np.linalg.norm(reduced(b[k])-reduced(a[k]), 'nuc') for k in a)/2
    assert separation >= record_separation-1e-13
    assert record_separation > 1/4000


def test_all721_record_reference_blocks_are_blind_to_diagonal_phases(source):
    data, _ = source
    rng = np.random.default_rng(107)
    psi = rng.normal(size=(66, 3))+1j*rng.normal(size=(66, 3))
    psi /= np.linalg.norm(psi)
    phases = np.exp(1j*rng.normal(size=66))
    for row in data['branches']:
        a = np.zeros_like(psi); b = np.zeros_like(psi)
        for i, j in row.pairs:
            a[j] = psi[i]
            b[j] = phases[i]*psi[i]
        ra = float(row.weight)*a.T@a.conj()
        rb = float(row.weight)*b.T@b.conj()
        assert np.linalg.norm(ra-rb) < 1e-15
    certificate = c.blindness_certificate(data)
    assert certificate['complete_raw_slots'] == 721
    assert certificate['partial_injection_columns'] == 7920


def test_terminal_coherent_readout_is_not_in_the_blindness_domain(source):
    data, probe = source
    i, j, _, _ = next(row for row in probe['pairs']
                     if not data['states'][row[0]][0] and not data['states'][row[0]][5])
    branch = next(row for row in data['branches'] if row.law == 'C'
                  and row.seam == (0, 5) and (i, i) in row.pairs and (j, j) in row.pairs)
    a = np.zeros(66, complex); a[i] = a[j] = 2**-.5
    b = a.copy(); b[i] *= 1j
    y = np.array(probe['Y'], complex)
    assert abs(np.vdot(a, y@a)) < 1e-14
    assert abs(np.vdot(b, y@b)) > .99
    assert branch.weight > 0


@pytest.mark.parametrize('damage', ['source', 'allowance', 'word'])
def test_owning_detector_receipt_cannot_be_weakened(damage):
    receipt = c.strict_json.load(c.observable_packet.OUTPUT)
    if damage == 'source':
        key = next(iter(receipt['source_bindings']))
        receipt['source_bindings'][key] = '0'*64
    elif damage == 'allowance':
        receipt['compiler']['words']['Y']['complete_halfdiamond_upper'] = '0'
    else:
        receipt['compiler']['words']['Y']['word_digest'] = '0'*64
    with pytest.raises(ArithmeticError):
        c.detector_contract(receipt)


def test_preparation_and_compiler_errors_are_charged():
    detector = c.detector_contract(); result = c.operational_certificate(detector)
    assert F(result['complete_output_halftrace_lower']) > F(1, 5000)
    assert F(result['each_complete_output_error_upper']) > c.CLOCK_ERROR
    detector['complete_halfdiamond_upper'] = '1/1000'
    with pytest.raises(ArithmeticError):
        c.operational_certificate(detector)


def test_receipt_is_exact_source_replay():
    assert packet.verify()
    value = packet.produce()
    assert value['scope']['a_transfer_clock_service_is_provided'] is False
    assert value['scope']['clock_specification_readout_is_forbidden'] is False
