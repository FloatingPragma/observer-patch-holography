"""Independent complete-source matrices, noncommuting variations and controls."""
import copy
from fractions import Fraction as F
import math
import subprocess
import sys

import numpy as np
import pytest
from scipy.integrate import quad

from model.field_thermal.thermal import car_restriction, restrict
from model.full_fock_field.field import FieldGeometry, FullFock
from model.prepared_regional_entropy import entropy as e, packet
from model.regional_modular_expansion import expansion
from tools import strict_json


def log_matrix(a):
    values, basis = np.linalg.eigh(a)
    assert min(values) > 0
    return (basis*np.log(values))@basis.conj().T


def entropy(a):
    values = np.linalg.eigvalsh(a)
    assert min(values) > 0
    return -float(np.dot(values, np.log(values)))


def directions():
    out = []
    for i in range(15):
        d = np.zeros((16, 16), complex)
        d[i, i], d[15, 15] = 1, -1
        out.append(d)
    for i in range(16):
        for j in range(i+1, 16):
            for coefficient in (1, 1j):
                d = np.zeros((16, 16), complex)
                d[i, j], d[j, i] = coefficient, coefficient.conjugate()
                out.append(d)
    return out


@pytest.fixture(scope='module')
def source():
    # Independent full-Fock assembly and spectral polynomial evaluation.
    # This checks the analytic bounds on the exact design reference. It is
    # not presented as reexecution of all literal preparation controls.
    h = FullFock(8).hamiltonian(FieldGeometry((2, 2, 1))).toarray()
    eigenvalues, eigenvectors = np.linalg.eigh(h-np.eye(256))
    q = sum(((-eigenvalues/64)**n/math.factorial(n) for n in range(9)),
            np.zeros(256))
    full = (eigenvectors*(q*q/np.dot(q, q)))@eigenvectors.conj().T
    maps = car_restriction(8, (0, 1, 2, 3))
    regional = restrict(full, maps)
    polynomial = np.asarray(expansion.polynomial(expansion.plaquette()), complex)
    prior = strict_json.loads((packet.ROOT/e.PRIOR).read_text())
    bounds = e.source_bounds(prior)
    return full, maps, regional, polynomial, bounds


def test_all_directions_and_parity_subspaces_have_exact_rank():
    basis = directions()
    flat = np.asarray([x.ravel() for x in basis])
    gram = (flat.conj()@flat.T).real
    expected = np.zeros((255, 255))
    expected[:15, :15] = np.eye(15)+np.ones((15, 15))
    expected[15:, 15:] = 2*np.eye(240)
    np.testing.assert_array_equal(gram, expected)
    parity = np.diag([(-1)**i.bit_count() for i in range(16)])
    number = np.diag([i.bit_count() for i in range(16)])
    assert sum(np.array_equal(x@parity, parity@x) for x in basis) == 127
    assert sum(np.array_equal(x@number, number@x) for x in basis) == 69
    result = e.direction_certificate()
    assert result['real_dimension'] == len(basis) == 255
    assert result['parity_preserving_real_dimension'] == 127
    assert result['number_preserving_real_dimension'] == 69
    assert result['full_Gram_determinant'] == 16*2**240
    for d in basis:
        assert np.trace(d) == 0
        assert np.linalg.norm(d, 2) == 1
        assert np.linalg.norm(d, ord='nuc') == 2


def test_same_complete_source_global_lift_and_noncommuting_finite_entropy(source):
    full, maps, sigma, polynomial, bound = source
    k = -log_matrix(sigma)
    assert min(np.linalg.eigvalsh(sigma)) > float(bound['regional_minimum'])
    assert min(np.linalg.eigvalsh(full)) > float(bound['full_minimum'])
    basis = directions()
    # Include imaginary, number-changing, parity-odd and even directions.
    selected = [basis[i] for i in (0, 4, 14, 15, 16, 20, 41, 127, 200, 254)]
    assert any(np.linalg.norm(sigma@d-d@sigma) > 1e-5 for d in selected)
    for d in selected:
        x = d/40
        # Each CAR restriction map selects a complete complement basis state.
        lifted = sum((a.conj().T@x@a for a in maps), np.zeros((256, 256), complex))/16
        np.testing.assert_allclose(restrict(lifted, maps), x, atol=1e-17)
        np.testing.assert_allclose(np.trace(full+lifted), np.trace(full), atol=1e-15)
        assert min(np.linalg.eigvalsh(full+lifted)) > 1/640
        assert min(np.linalg.eigvalsh(full-lifted)) > 1/640
        relative = float(np.trace((sigma+x)@(log_matrix(sigma+x)-log_matrix(sigma))).real)
        change = entropy(sigma+x)-entropy(sigma)
        pairing = float(np.trace(x@k).real)
        assert change == pytest.approx(pairing-relative, abs=3e-15)
        exact = e.remainder_bounds(F(1, 40), minimum=bound['regional_minimum'],
                                   modular_error=bound['modular_polynomial_error'])
        assert float(exact['relative_entropy_lower']) <= relative <= float(exact['relative_entropy_upper'])
        assert abs(change-np.trace(x@polynomial).real) <= float(exact['polynomial_finite_entropy_error_upper'])


def test_hessian_resolvent_matches_finite_differences_without_commutation(source):
    _, _, sigma, _, _ = source
    d = directions()[20]
    assert np.linalg.norm(sigma@d-d@sigma) > 1e-5
    values, frame = np.linalg.eigh(sigma)
    transformed = frame.conj().T@d@frame
    hs = float(np.vdot(d, d).real)
    resolvent, error = quad(lambda t: float(np.sum(
        abs(transformed)**2/((values[:, None]+t)*(values[None, :]+t))).real),
        0, np.inf, epsabs=1e-10)
    assert error < 1e-6
    step = 1e-4
    finite = -(entropy(sigma+step*d)+entropy(sigma-step*d)-2*entropy(sigma))/step**2
    assert finite == pytest.approx(resolvent, rel=2e-6)
    assert hs <= resolvent <= hs/min(values)
    # The entrywise reciprocal shortcut is wrong for a noncommuting direction.
    shortcut = float(np.trace(d@d@np.linalg.inv(sigma)).real)
    assert abs(shortcut-resolvent) > 1e-6


def test_symmetric_estimator_and_cubic_remainder(source):
    _, _, sigma, polynomial, bound = source
    log_sigma = log_matrix(sigma)
    d = directions()[0]
    values, frame = np.linalg.eigh(sigma)
    transformed = frame.conj().T@d@frame
    metric = quad(lambda t: float(np.sum(
        abs(transformed)**2/((values[:, None]+t)*(values[None, :]+t))).real),
        0, np.inf, epsabs=1e-10)[0]
    for amplitude in (F(1, 40), F(1, 400), F(1, 4000)):
        h = float(amplitude)
        a, b = sigma+h*d, sigma-h*d
        exact = e.remainder_bounds(amplitude, minimum=bound['regional_minimum'],
                                   modular_error=bound['modular_polynomial_error'])
        estimate = (entropy(a)-entropy(b))/(2*h)
        assert abs(estimate-np.trace(d@polynomial).real) <= float(exact['symmetric_entropy_derivative_error_upper'])+1e-12
        relative = np.trace(a@(log_matrix(a)-log_sigma)).real
        assert abs(relative-h*h*metric/2) <= float(exact['quadratic_Taylor_remainder_upper'])+1e-14


def test_subnormalization_abort_complement_and_kernel_countercontrol():
    sigma = np.array([[.3, .04j], [-.04j, .2]])
    x = np.array([[.01, .005], [.005, -.01]])
    abort = np.diag([.2, .3])
    before = np.block([[sigma, np.zeros((2, 2))], [np.zeros((2, 2)), abort]])
    after = np.block([[sigma+x, np.zeros((2, 2))], [np.zeros((2, 2)), abort]])
    change = entropy(sigma+x)-entropy(sigma)
    assert entropy(after)-entropy(before) == pytest.approx(change, abs=5e-16)
    relative = np.trace((sigma+x)@(log_matrix(sigma+x)-log_matrix(sigma))).real
    assert change == pytest.approx(-np.trace(x@log_matrix(sigma)).real-relative, abs=5e-16)
    # Normalization would multiply this entropy difference by 1/p=2.
    conditioned = entropy(2*(sigma+x))-entropy(2*sigma)
    assert abs(conditioned-change) > 1e-3
    # A kernel reference admits no two-sided full variation ball.
    kernel = np.diag([1., 0.])
    direction = np.diag([-1., 1.])
    assert min(np.linalg.eigvalsh(kernel-.001*direction)) < 0


@pytest.mark.parametrize('amplitude,kwargs', [
    (.1, {}), (True, {}), (F(-1, 100), {}), (F(1, 20), {}),
    (0, {'minimum': 0}), (0, {'modular_error': -1}),
    (0, {'operator_norm': 0, 'hs_squared': 1}),
])
def test_bounds_reject_nonrational_or_nonpositive_domains(amplitude, kwargs):
    args = {'minimum': F(1, 20), 'modular_error': F(1, 10**9)}
    args.update(kwargs)
    with pytest.raises(ValueError):
        e.remainder_bounds(amplitude, **args)


def test_corrupted_preparation_floor_and_success_normalization_fail_closed():
    prior = strict_json.loads((packet.ROOT/e.PRIOR).read_text())
    for key in ('actual_regional_success_minimum_weight',
                'actual_full_success_minimum_weight',
                'regional_P8_total_operator_and_trace_zero_tangent_bound'):
        damaged = copy.deepcopy(prior)
        damaged['modular_certificate'][key] = '0'
        with pytest.raises(ValueError, match='recomputation'):
            e.source_bounds(damaged)
    damaged = copy.deepcopy(prior)
    damaged['modular_certificate']['regional_log_is_of_unnormalized_success_block'] = False
    with pytest.raises(ValueError, match='normalization'):
        e.source_bounds(damaged)
    damaged = copy.deepcopy(prior)
    path = next(iter(damaged['source_sha256']))
    damaged['source_sha256'][path] = '0'*64
    with pytest.raises(ValueError, match='custody'):
        e.source_bounds(damaged)
    damaged = copy.deepcopy(prior)
    damaged['source_sha256'].pop(path)
    with pytest.raises(ValueError, match='custody coverage'):
        e.source_bounds(damaged)
    damaged = copy.deepcopy(prior)
    damaged['domain']['interaction'] = '0'
    with pytest.raises(ValueError, match='domain'):
        e.source_bounds(damaged)


def test_retained_exact_receipt_and_fresh_process_mutation(tmp_path):
    result = packet.report()
    assert packet.canonical(result) == packet.canonical(strict_json.loads(packet.OUTPUT.read_text()))
    assert result['scope']['record_resolved_lift_or_native_intervention_executed'] is False
    assert result['scope']['entropy_area_stationarity_proved'] is False
    result['uniform_bounds']['relative_entropy_HS_squared_upper_coefficient'] = 1
    wrong = tmp_path/'wrong.json'
    wrong.write_text(packet.canonical(result))
    run = subprocess.run([sys.executable, '-m', 'model.prepared_regional_entropy', '--check', str(wrong)],
                         cwd=packet.ROOT, capture_output=True, text=True)
    assert run.returncode != 0
    assert 'differs from its exact replay' in run.stderr
