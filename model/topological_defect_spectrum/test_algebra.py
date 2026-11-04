from itertools import product

import pytest
import sympy as s

from . import algebra


def test_full_source_internal_fock_map_against_independent_jordan_wigner():
    gamma, beta, rho, b, chirality = algebra.generators()
    unitary = algebra.internal_fock_basis()
    identity, parity = s.eye(2), s.diag(1, -1)
    lower, y = s.Matrix([[0, 1], [0, 0]]), s.Matrix([[0, -s.I], [s.I, 0]])
    for axis in range(3):
        factors = [parity if j < axis else lower if j == axis else identity for j in range(3)]
        oracle = s.kronecker_product(*factors, identity)
        assert b[axis]*unitary == unitary*oracle
        rho_oracle = s.kronecker_product(*(y if j == axis else identity for j in range(3)), identity)
        assert rho[axis] == rho_oracle
    assert unitary.H*unitary == s.eye(16)
    assert algebra.clifford_certificate()['internal_Fock_multiplicity'] == 2
    for orientation in (-1, 1):
        basis = algebra.zero_basis(orientation)
        assert basis.H*basis == s.eye(2)
        assert chirality*basis == orientation*basis
    assert sum((algebra.occupation_projector(x) for x in product((0, 1), repeat=3)), s.zeros(16)) == s.eye(16)


def test_squared_differential_operator_on_independent_polynomial_columns():
    x = s.symbols('x0:3', real=True)
    v = s.Rational(2)
    phi = (x[0]**2+2*x[1]**2+3*x[2]**2)/2+x[0]*x[1]/3+x[2]**4/4
    gradient = [s.diff(phi, z) for z in x]
    hessian = s.hessian(phi, x)
    gamma, beta, _, b, _ = algebra.generators()
    test = s.Matrix([j+1+(j % 3+1)*x[j % 3]+s.I*(j+2)*x[(j+1) % 3]*x[(j+2) % 3]
                     for j in range(16)])
    def H(vector):
        return sum((-s.I*v*gamma[j]*vector.diff(x[j])+beta[j]*gradient[j]*vector
                    for j in range(3)), s.zeros(16, 1)).applyfunc(s.expand)
    actual = H(H(test))
    expected = s.zeros(16, 1)
    for j in range(3):
        A = v*test.diff(x[j])+gradient[j]*test
        expected += -v*A.diff(x[j])+gradient[j]*A
    expected += 2*v*sum((b[i].H*b[j]*hessian[i, j]*test for i in range(3) for j in range(3)),
                       s.zeros(16, 1))
    assert (actual-expected).applyfunc(s.expand) == s.zeros(16, 1)


@pytest.mark.parametrize('level', range(5))
def test_complete_oscillator_blocks_and_signed_multiplicities(level):
    operator = algebra.oscillator_block(level)
    certificate = algebra.oscillator_multiplicity(level)
    assert operator == operator.H
    assert operator*operator == level*s.eye(operator.rows)
    assert operator.rows == certificate['squared_energy_multiplicity']
    assert s.trace(operator) == 0
    if level:
        assert operator.rows//2 == 4*level*level+2
    else:
        assert operator == s.zeros(2)


def test_smooth_positive_hessian_and_invalid_domains():
    assert algebra.radial_hessian_certificate()['normalizable_zero_mode_dimension'] == 2
    for value in (True, -1, 1.0):
        with pytest.raises(ValueError):
            algebra.oscillator_multiplicity(value)
    for value in (True, 0, 2):
        with pytest.raises(ValueError):
            algebra.zero_basis(value)
    with pytest.raises(ValueError):
        algebra.occupation_projector((0, True, 1))
