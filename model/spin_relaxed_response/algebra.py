"""Independent exact finite checks of resolvent relaxation and its certificate."""
from fractions import Fraction as F
import sympy as sp


def exact_checks():
    # Positive excitations of one actual antiferromagnetic bond. The scalar
    # subtraction makes the singlet energy zero and every triplet energy one.
    sx = sp.Matrix([[0, 1], [1, 0]])/2
    sy = sp.Matrix([[0, -sp.I], [sp.I, 0]])/2
    sz = sp.diag(1, -1)/2
    identity = sp.eye(2)
    h = sum((sp.kronecker_product(s, s) for s in (sx, sy, sz)), sp.zeros(4))-sp.eye(4)/4
    ground = sp.Matrix([0, 1, -1, 0])/sp.sqrt(2)
    k = h+sp.eye(4)
    a = sp.kronecker_product(sz, identity)*ground
    eta = sp.Rational(1, 3)
    r = (k+eta*sp.eye(4)).inv()
    vector = r*a
    chi = (a.conjugate().T*vector)[0]
    norm = (vector.conjugate().T*vector)[0]
    energy = (vector.conjugate().T*k*vector)[0]
    assert h*ground == -ground and k*ground == sp.zeros(4, 1)
    assert chi == sp.Rational(3, 16) and norm == sp.Rational(9, 64)
    assert energy == chi-eta*norm == sp.Rational(9, 64)
    tangent = vector/(2*chi)
    assert 2*(a.conjugate().T*tangent)[0] == 1
    assert (tangent.conjugate().T*k*tangent)[0] == 1
    # The exact susceptibility is 1/4, so the constrained energy infimum is 1.
    assert (a.conjugate().T*a)[0] == sp.Rational(1, 4)
    angle = sp.Symbol('s', real=True, nonzero=True)
    probe = sp.kronecker_product(sz, identity)
    sine = 2*sp.sin(angle/2)*probe
    cosine = sp.cos(angle/2)*sp.eye(4)
    assert sp.simplify(sine.conjugate().T*sine+cosine.conjugate().T*cosine) == sp.eye(4)
    bounded_vector = sine*ground/angle
    bounded_chi = sp.simplify((bounded_vector.conjugate().T*r*bounded_vector)[0])
    assert sp.limit(bounded_chi, angle, 0) == chi
    # Complex, deliberately inaccurate trial vectors test the gap-free dual
    # residual identity independently of the optimal-vector construction.
    trial = sp.Matrix([sp.I/7, sp.Rational(1, 5), sp.Rational(-1, 9), sp.Rational(1, 11)])
    residual = a-(k+eta*sp.eye(4))*trial
    dual = 2*sp.re((a.conjugate().T*trial)[0])-(trial.conjugate().T*(k+eta*sp.eye(4))*trial)[0]
    defect = sp.simplify(chi-dual)
    assert sp.simplify(defect-(residual.conjugate().T*r*residual)[0]) == 0
    upper = (residual.conjugate().T*residual)[0]/eta
    assert sp.simplify(defect) >= 0 and sp.simplify(upper-defect) >= 0
    # A two-frequency exact spectral law prevents replacing the relaxed
    # response by bare energy divided by the variance squared.
    weights = (F(1, 3), F(2, 3)); energies = (F(1), F(4))
    mass = sum(weights); first = sum(w*e for w, e in zip(weights, energies))
    inverse = sum(w/e for w, e in zip(weights, energies))
    assert inverse == F(1, 2) and first == 3
    assert F(1, 4)/inverse == F(1, 2) < first/(4*mass*mass) == F(3, 4)
    return {'native_single_bond_ground_exact': True,
            'native_single_bond_regularized_inverse': '3/16',
            'native_single_bond_inverse': '1/4',
            'native_single_bond_unit_tangent_cost': '1',
            'complex_trial_dual_residual_identity_exact': True,
            'complete_bounded_sine_probe_exact': True,
            'bounded_probe_inverse_response_limit_exact': True,
            'two_frequency_relaxed_cost': '1/2',
            'two_frequency_bare_probe_cost': '3/4',
            'spectral_control_is_an_ordered_infinite_ground': False}
