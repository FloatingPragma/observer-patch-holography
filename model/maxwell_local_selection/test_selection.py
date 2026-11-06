import pytest
import sympy as s

from . import certificate as c


def test_exhaustive_local_quadratic_classification():
    result = c.classify()
    assert result['local_gauge_quadratic_dimension_exact'] == 20
    assert result['rotation_covariant_local_gauge_dimension_exact'] == 1
    assert result['invariant_linear_seam_dimension_exact'] == 0


def test_random_exact_weighted_faces_and_gauge_mutation():
    D, C, _ = c.carrier()
    weights = [s.Rational(1+i*i, 1+i) for i in range(20)]
    H = C.T*s.diag(*weights)*C
    assert c.is_local_gauge_hessian(H)
    assert not c.is_local_gauge_hessian(H+D*D.T)
    # Orientation erasure changes the gauge object despite identical support.
    unsigned = C.applyfunc(abs)
    assert not c.is_local_gauge_hessian(unsigned.T*unsigned)


def test_all_hypotheses_have_exact_countercontrols():
    result = c.controls()
    assert result['gyroscopic_term_keeps_energy_Faraday_Gauss_but_breaks_time_reversal']
    assert result['H_squared_keeps_gauge_covariance_but_breaks_face_locality']


def test_evolution_keeps_actual_carrier_constraints_and_coefficient():
    D, C, _ = c.carrier()
    G, W, R = c.evolution(s.Rational(7, 3), s.Rational(11, 5))
    E, B = C.T*s.ones(20, 1), C*s.Matrix(range(30))
    # Replace the vanishing all-face electric seed with an actual face cycle.
    E = C.row(0).T
    state = E.col_join(B)
    derivative = G*state
    assert derivative[:30, :] == s.Rational(33, 35)*C.T*B
    assert D.T*derivative[:30, :] == s.zeros(12, 1)
    assert s.ones(1, 20)*derivative[30:, :] == s.zeros(1, 1)
    assert (state.T*W*derivative)[0] == 0
    assert R*G*R == -G


@pytest.mark.parametrize('electric,magnetic', [(0, 1), (1, -1), (1.0, 1)])
def test_no_nonpositive_or_rounded_coefficient_admission(electric, magnetic):
    with pytest.raises(ValueError):
        c.evolution(electric, magnetic)
