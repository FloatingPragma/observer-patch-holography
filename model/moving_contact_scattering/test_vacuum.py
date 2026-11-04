"""Independent CAR matrix controls for the additive Fock obstruction."""
from fractions import Fraction as F
from functools import lru_cache

import numpy as np
import pytest
from scipy import sparse
import sympy as s

from . import vacuum


@lru_cache(maxsize=2)
def jordan_wigner(modes):
    lowering = sparse.csr_matrix([[0., 1.], [0., 0.]])
    parity = sparse.diags([1., -1.], format='csr')
    identity = sparse.eye(2, format='csr')
    result = []
    for target in range(modes):
        operator = sparse.csr_matrix([[1.]])
        for address in reversed(range(modes)):
            factor = lowering if address == target else parity if address < target else identity
            operator = sparse.kron(operator, factor, format='csr')
        result.append(operator.astype(complex))
    return result


def independent_operator(sites, *, include_cross=True):
    annihilators = jordan_wigner(2*sites)
    number = [c.conjugate().T@c for c in annihilators]
    operator = sum((4*n for n in number), sparse.csr_matrix(number[0].shape, dtype=complex))
    operator -= sum((10*number[2*x]@number[2*x+1] for x in range(sites)),
                    sparse.csr_matrix(number[0].shape, dtype=complex))
    links = [(0, 1)] if sites == 2 else [(0, 2), (2, 1), (1, 3), (3, 0)]
    if sites == 4 and not include_cross:
        links = [(0, 2), (1, 3)]
    for x, y in links:
        for spin in (0, 1):
            term = annihilators[2*x+spin].conjugate().T@annihilators[2*y+1-spin]
            operator += -1j*term+1j*term.conjugate().T
    return operator.toarray()


def test_independent_jordan_wigner_source_and_all_particle_lower_bounds():
    matrix = independent_operator(4)
    exact = np.asarray(vacuum.finite_operator(4)).astype(complex)
    np.testing.assert_array_equal(matrix, exact)
    np.testing.assert_array_equal(matrix.conj().T, matrix)
    for number in range(9):
        states = [mask for mask in range(256) if mask.bit_count() == number]
        energies = np.linalg.eigvalsh(matrix[np.ix_(states, states)])
        assert energies[0] >= float(vacuum.fixed_number_lower_bound(number, 10))-1e-12
        assert energies[-1] <= 8*number+1e-12
    # All pairwise CAR identities are checked on the entire 256-dimensional space.
    annihilators = jordan_wigner(8)
    identity = sparse.eye(256, format='csr')
    for i, a in enumerate(annihilators):
        for j, b in enumerate(annihilators):
            anticomm = a@b.conjugate().T+b.conjugate().T@a
            expected = identity if i == j else sparse.csr_matrix((256, 256))
            assert (anticomm-expected).nnz == 0


def test_independent_fermionic_product_has_additive_energy_with_cross_links():
    coefficients = {3: 1, 12: 2, 5: -1j}
    creators = [c.conjugate().T for c in jordan_wigner(8)]
    state = np.zeros(256, dtype=complex)
    for left, x in coefficients.items():
        for right, y in coefficients.items():
            word = ([address for j, address in enumerate((0, 1, 4, 5)) if left & (1 << j)]+
                    [address for j, address in enumerate((2, 3, 6, 7)) if right & (1 << j)])
            component = np.eye(1, 256, 0, dtype=complex).ravel()
            for address in reversed(word):
                component = creators[address]@component
            state += x*y*component
    exact = vacuum.wedge_vector({3: 1, 12: 2, 5: -s.I}, {3: 1, 12: 2, 5: -s.I})
    np.testing.assert_array_equal(state, np.asarray(exact).astype(complex).ravel())
    full = independent_operator(4)
    internal = independent_operator(4, include_cross=False)
    assert np.count_nonzero(full-internal) > 0
    assert np.vdot(state, (full-internal)@state) == 0
    assert np.vdot(state, full@state)/np.vdot(state, state) == pytest.approx(-8/3)
    local = np.zeros(16, complex)
    for mask, value in coefficients.items():
        local[mask] = value
    assert np.vdot(local, independent_operator(2)@local)/np.vdot(local, local) == pytest.approx(-4/3)
    result = vacuum.finite_additivity_certificate()
    assert result['graded_product_energy'] == '-8/3'
    assert result['ungraded_product_energy'] == '-2'


@pytest.mark.parametrize('gamma', [F(1, 7), F(1), F(11, 3)])
def test_negative_compact_packet_budget_at_boundary(gamma):
    result = vacuum.negative_packet_certificate(1/(4096*gamma), gamma)
    # Independently bound the bracket using |P|^2<=3 gamma^2/16.
    assert result['bound_wavepacket_energy_upper'] < -gamma*gamma/8
    assert result['compact_energy_approximation_error_upper'] == gamma*gamma/32
    assert (result['bound_wavepacket_energy_upper']+
            result['compact_energy_approximation_error_upper']) < -gamma*gamma/16
    assert result['compact_support_radius_effectively_computed'] is False


@pytest.mark.parametrize('a,gamma', [(F(1, 4095), 1), (0, 1), (1, 0), (True, 1), (1, -1)])
def test_negative_packet_domain_rejects_invalid_inputs(a, gamma):
    with pytest.raises(ValueError):
        vacuum.negative_packet_certificate(a, gamma)


def test_particle_bound_domain_and_empty_sector():
    assert vacuum.fixed_number_lower_bound(0, 10) == 0
    assert vacuum.fixed_number_lower_bound(3, 10) == -10
    assert vacuum.fixed_number_lower_bound(4, 10) == -20
    for number in (True, F(2), -1, 2.0):
        with pytest.raises(ValueError):
            vacuum.fixed_number_lower_bound(number, 10)
    with pytest.raises(ValueError):
        vacuum.fixed_number_lower_bound(2, -1)


def test_lorentz_quartic_obstruction_survives_affine_dictionary():
    radius = s.symbols('radius', real=True)
    direction = s.Matrix([1, 2, -1])
    linear = s.Matrix([[2, 1, 0], [0, 1, 0], [0, 0, 3]])
    shift = s.Matrix([3, -2, 1])
    momentum = linear*(radius*direction)+shift
    alpha, speed, mass, offset, mu = 3, 5, 7, 11, 13
    energy = alpha*(radius**2*direction.dot(direction)/4-1+2*mu)+offset
    polynomial = s.Poly(energy**2-speed**2*momentum.dot(momentum)-mass**2*speed**4, radius)
    assert polynomial.coeff_monomial(radius**4) == s.Rational(81, 4)
    certificate = vacuum.affine_shell_certificate()
    assert certificate['exact_Lorentz_shell_on_any_open_momentum_set'] is False
    assert certificate['nonrelativistic_low_momentum_approximation_excluded'] is False
