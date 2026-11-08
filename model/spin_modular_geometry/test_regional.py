from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm

from model.native_spin_observer import algebra
from model.native_ordered_ground import source
from . import pulse, regional


def independent_source():
    def site(axis, index):
        word = ['I']*8
        word[index] = axis
        return algebra.tensor(''.join(word))
    spin = {(x, a): site(a, x) for x in range(8) for a in 'XYZ'}
    h = np.zeros((256, 256), complex)
    e = np.zeros_like(h)
    edges, odd = source.cube()
    for x, y in edges:
        bond = (sum(spin[x, a]@spin[y, a] for a in 'XYZ')-np.eye(256))/4
        h += bond
        if 0 in (x, y): e += bond/2
    for x in range(8):
        h += (1 if x in odd else -1)*spin[x, 'X']/8
    e -= spin[0, 'X']/8
    return h, e


def restrict(vector):
    matrix = vector.reshape(64, 4).T
    return matrix@matrix.conj().T


def entropy(matrix):
    eigenvalues = np.linalg.eigvalsh(matrix)
    return -np.sum(eigenvalues*np.log(eigenvalues))


def test_complete_source_marginal_and_exact_sylvester_floor():
    h, e = independent_source()
    data = regional.state()
    assert np.array_equal(h, data['energy_integer']/8)
    assert np.array_equal(e, data['density_integer']/8)
    eigenvalues, eigenvectors = np.linalg.eigh(h)
    assert eigenvalues[1]-eigenvalues[0] > 1
    rho = restrict(eigenvectors[:, 0])
    assert np.linalg.norm(rho-np.asarray(data['rho'], float), 2) < float(data['ground_vector_error'])
    assert np.linalg.eigvalsh(rho)[0] > 1/41
    assert len(data['principal_minors']) == 4
    assert all(v > 0 for v in data['principal_minors'])
    assert data['full_source_certificate']['all_256_eigenvalues_certified'] is True


def test_modular_logarithm_against_independent_diagonalization():
    h, _ = independent_source()
    _, vectors = np.linalg.eigh(h)
    rho = restrict(vectors[:, 0])
    values, columns = np.linalg.eigh(rho)
    exact = -(columns*np.log(values))@columns.conj().T
    for order, bits in ((1024, 96), (128, 32)):
        row = regional.logarithm(order, bits)
        assert np.linalg.norm(np.asarray(row['matrix'], float)-exact, 2) < float(row['operator_error_upper'])
    response = regional.source_response()
    assert response['XX_minus_YY']+response['coefficient_error_upper'] < F(-13, 100)
    # A scaled isotropic bond plus arbitrary local fields cannot have this coefficient.
    k = 0.7*(algebra.tensor('XX')+algebra.tensor('YY')+algebra.tensor('ZZ'))
    k += 1.3*algebra.tensor('XI')-2.1*algebra.tensor('IZ')+2*np.eye(4)
    assert abs(np.trace(k@(algebra.tensor('XX')-algebra.tensor('YY')))) < 1e-14


def test_actual_entropy_tangent_and_finite_source_pulse_independently():
    h, e = independent_source()
    _, vectors = np.linalg.eigh(h)
    ground = vectors[:, 0]
    v = 1j*(h@e-e@h)
    source_terms = sum((float(value)*algebra.tensor(label) for label, value in pulse.current_terms().items()),
                       np.zeros((256, 256), complex))
    assert np.array_equal(source_terms, v)
    expectation = lambda lam: entropy(restrict(expm(-1j*lam*v)@ground))
    slope = (expectation(1e-5)-expectation(-1e-5))/(2e-5)
    result = regional.source_response()
    assert abs(slope-float(result['regional_entropy_slope'])) < 1e-10
    assert abs(np.vdot(ground, 1j*(v@h-h@v)@ground)) < 1e-13
    bound = regional.finite_pulse()
    actual = expectation(1e-4)-expectation(0)
    assert actual < float(bound['actual_regional_entropy_change_upper']) < -4e-6
    assert np.linalg.norm(np.asarray(regional.state()['tangent'], float)) > 0


def test_general_complete_pauli_fisher_is_positive_and_has_normalized_tangent():
    rng = np.random.default_rng(1145)
    labels = [a+b for a in 'IXYZ' for b in 'IXYZ']
    rho = np.diag([0.1, 0.2, 0.3, 0.4])
    for _ in range(7):
        raw = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
        delta = raw+raw.conj().T
        delta -= np.trace(delta)*np.eye(4)/4
        fisher, total, rate = 0.0, 0.0, 0.0
        for label in labels:
            p = algebra.tensor(label)
            r, a = np.trace(rho@p).real, np.trace(delta@p).real
            total += (1+r)/32+(1-r)/32
            rate += a/32-a/32
            if label != 'II': fisher += a*a/(16*(1-r*r))
        hs = np.trace(delta@delta).real
        assert total == 1 and rate == 0
        assert 0 < hs/4 <= fisher <= hs/(32*0.1*(1-0.2))
    response = regional.source_response()
    records = response['complete_pauli_response_records']
    assert len(records) == 16
    assert sum(sum(r['probabilities']) for r in records) == 1
    assert sum(sum(r['derivatives']) for r in records) == 0
    assert records[0]['probabilities'] == (F(1, 16), F(0))
    assert response['actual_source_fisher_lower'] > F(21, 1000)
    assert response['actual_source_fisher_upper'] < F(22, 1000)


def test_same_clock_average_identity_and_zero_current_limit():
    h, e = independent_source()
    v = 1j*(h@e-e@h)
    result = regional.current_filter(F(16))
    assert (1j*(h@v-v@h))[0, 1] == F(3, 16)
    for time in (1, 3, 16):
        u = expm(1j*time*h)
        averaged = (u@e@u.conj().T-e)/time
        derivative = (u@v@u.conj().T-v)/time
        assert np.linalg.norm(1j*(h@averaged-averaged@h)-derivative, 2) < 1e-13
        assert np.linalg.norm(averaged, 2) <= 13/(4*time)
        assert np.linalg.norm(derivative, 2) <= 3/time
    assert result['nonzero_spatial_translation_generator_constructed'] is False


@pytest.mark.parametrize('args', [(0, 96), (1024, 8), (True, 96), (1.5, 96)])
def test_logarithm_domain(args):
    with pytest.raises(ValueError): regional.logarithm(*args)


@pytest.mark.parametrize('value', [F(0), F(-1), F(1), 0.0001, True])
def test_finite_pulse_requires_faithful_exact_domain(value):
    with pytest.raises(ValueError): regional.finite_pulse(value)
