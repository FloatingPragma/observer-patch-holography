from itertools import product
import pytest
import sympy as s

from model.weyl_geometric_modular import algebra as a


def test_all_valleys_and_chiralities_and_exact_source_symbol():
    rows = a.source_frames()['rows']
    assert {tuple(row['valley']) for row in rows} == set(product((0, 1), repeat=3))
    assert [row['chirality'] for row in rows].count(1) == 4
    p = (s.Rational(2, 3), -s.Rational(1, 4), 7)
    for valley in a.VALLEYS:
        signs, chi, U = a.frame(valley)
        native = sum((-signs[j]*p[j]*a.pauli()[j] for j in range(3)), s.zeros(2))
        assert U.adjoint()*native*U == chi*sum((p[j]*a.pauli()[j] for j in range(3)), s.zeros(2))


@pytest.mark.parametrize('value', [(0, 1), (0, 1, True), [0, 0, 0], (0, 0, 2)])
def test_bad_valleys_rejected(value):
    with pytest.raises(ValueError): a.frame(value)


def test_stress_polynomial_shell_identities_and_unsymmetrized_control():
    assert a.stress_certificate()['Ward_polynomial_identities'] == 4
    p = (3, 1, 2, 2); q = (5, 3, 4, 0)
    assert a.stress_kernel(0, 1, p, q) == a.stress_kernel(1, 0, p, q)
    assert a.stress_kernel(0, 1, p, q, symmetric=False) != a.stress_kernel(1, 0, p, q, symmetric=False)


def test_gradient_correction_required_not_just_constant_fit():
    assert a.zero(a.weighted_generator()['formal_adjoint_defect'])
    for coefficient in (0, 1, s.Rational(1, 3), s.Rational(-1, 2)):
        assert not a.zero(a.weighted_generator(coefficient)['formal_adjoint_defect'])
    with pytest.raises(ValueError): a.weighted_generator(0.5)


def test_projector_commutator_all_eight_clifford_blocks():
    assert a.projector_certificate()['first_derivative_identities'] == 24
    # A constant wrong projector need not commute even with the principal symbol.
    P = s.diag(1, 0)
    assert not a.zero(a.pauli()[0]*P-P*a.pauli()[0])


def test_full_fock_weighted_charge_generates_the_same_car_motion():
    # Independent finite CAR calculation of 1/2{w,h}, normal ordering and derivation.
    n = 3; dim = 1 << n; c = []
    for j in range(n):
        op = s.zeros(dim)
        for col in range(dim):
            if col & (1 << j):
                op[col ^ (1 << j), col] = (-1)**((col & ((1 << j)-1)).bit_count())
        c.append(op)
    h = s.Matrix([[1, s.I, 0], [-s.I, -2, 1], [0, 1, 3]])
    w = s.diag(2, -1, s.Rational(1, 3)); q = (w*h+h*w)/2
    Q = sum((q[i,j]*c[i].adjoint()*c[j] for i in range(n) for j in range(n)), s.zeros(dim))
    for j in range(n):
        assert Q*c[j]-c[j]*Q == -sum((q[j,k]*c[k] for k in range(n)), s.zeros(dim))
    # A scalar vacuum subtraction changes neither commutator nor normalization.
    Q -= s.Rational(7, 5)*s.eye(dim)
    assert Q == Q.adjoint()
    assert w*h != (w*h).adjoint()
