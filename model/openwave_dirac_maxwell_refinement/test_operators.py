from fractions import Fraction as F
import numpy as np
from scipy.linalg import eigh
from .operators import pauli_operator, polynomial_member, root_tail, charged_masses


def test_nonzero_magnetic_commutator_is_in_actual_square():
    p = pauli_operator(4, -3)
    assert np.linalg.norm(p['magnetic_commutator'].toarray()) > .01
    assert np.max(abs((p['D2']-p['scalar_square']-p['magnetic_commutator']).toarray())) == 0
    assert np.max(abs((p['D']-p['D'].conj().T).toarray())) == 0


def test_mass_dependent_root_is_not_old_common_hopping():
    low = polynomial_member(3, -3, F(1), regulator=0)
    high = polynomial_member(3, -3, F(207), regulator=0)
    assert np.linalg.norm(low['kinetic']-high['kinetic']) > .01
    assert np.linalg.norm(low['first_gauge_derivative']) > 100*np.linalg.norm(high['first_gauge_derivative'])
    assert F(low['uniform_root_operator_error_exact']) < F(1, 10**12)
    assert len(charged_masses()) == 3


def test_noncommuting_current_derivatives_match_full_source_probe():
    h = .00005
    d = polynomial_member(2, -3, F(1))
    p = polynomial_member(2, -3, F(1), theta=h)
    m = polynomial_member(2, -3, F(1), theta=-h)
    first = (p['kinetic']-m['kinetic'])/(2*h)
    second = (p['kinetic']+m['kinetic']-2*d['kinetic'])/(h*h)
    assert np.max(abs(first-d['first_gauge_derivative'])) < 1e-9
    assert np.max(abs(second-d['second_gauge_derivative'])) < 1e-8


def test_relative_dipole_has_exact_complete_commutator_and_contact():
    d = polynomial_member(2, -3, F(1))
    h = d['H_minus_mass']; x = np.diag(d['relative_x_dipole_diagonal'])
    j = 1j*(h@x-x@h)
    assert np.max(abs(j-d['relative_dipole_current'])) < 1e-16
    second = 2*x@h@x-x@x@h-h@x@x
    assert np.max(abs(second-d['relative_dipole_second_derivative'])) < 2e-14


def test_actual_Pauli_square_uses_the_Maxwell_symbol_normalization():
    from .operators import ELECTRIC, MAGNETIC, MOMENTUM_SQUARED, SPEED_SQUARED
    assert SPEED_SQUARED == 2*ELECTRIC*MAGNETIC == F(1, 128)
    assert MOMENTUM_SQUARED == F(1, 512)
    p = pauli_operator(3, -3)
    assert np.max(abs((p['D']@p['D']-p['D2']).toarray())) < 1e-17
