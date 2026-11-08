"""Exact source controls and adversarial boundaries for response closure."""
from fractions import Fraction as F
import numpy as np
import pytest
import sympy as s
from . import closure


def complex_control():
    # Complex coupling exercises adjoints in the off-response block.
    k = s.Matrix([[2, s.I, 0], [-s.I, 2, 0], [0, 0, 4]])
    v = s.Matrix([[1, 0], [0, 0], [0, s.I]])
    return k, v


def test_complex_response_defect_equals_leakage_with_positive_and_zero_directions():
    k, v = complex_control()
    row = closure.finite_defect(k, v, F(1))
    assert row['defect'] == s.diag(s.Rational(1, 64), 0)
    assert row['residual'] == s.Matrix([[0, 0], [s.I/8, 0], [0, 0]])
    assert row['compression'] == s.diag(s.Rational(3, 8), s.Rational(1, 5))
    assert row['defect'] == row['residual'].H*row['residual']
    # A nonreducing chosen space still contains a genuinely sharp branch.
    sharp = closure.finite_defect(k, v[:, 1], F(1))
    assert sharp['defect'] == s.zeros(1)


def test_telescoping_identity_for_every_power_and_exact_polynomial():
    k, v = complex_control()
    row = closure.finite_defect(k, v, F(3, 2))
    bounded, a, r = (row[key] for key in ('resolvent', 'compression', 'residual'))
    coefficients = [F(7), F(-2), F(3, 5), F(-1, 3), F(2)]
    for power in range(1, 7):
        rhs = sum((bounded**(power-1-i)*r*a**i for i in range(power)), s.zeros(3, 2))
        assert bounded**power*v-v*a**power == rhs
    coefficient = closure.polynomial_coefficient(coefficients)
    assert coefficient == F(61, 5)
    fs = sum((c*bounded**j for j, c in enumerate(coefficients)), s.zeros(3))
    fa = sum((c*a**j for j, c in enumerate(coefficients)), s.zeros(2))
    difference = np.array((fs*v-v*fa).evalf(), complex)
    defect = np.array(row['defect'].evalf(), complex)
    assert np.linalg.norm(difference, 2) <= float(coefficient)*np.sqrt(np.linalg.norm(defect, 2))


@pytest.mark.parametrize('clock,embedding', [
    ([[1, 1], [0, 1]], [[1], [0]]),
    ([[-1, 0], [0, 1]], [[1], [0]]),
    ([[1, 0], [0, 2]], [[2], [0]]),
    ([[1, 0], [0, 2]], [[1]]),
    ([], []),
])
def test_nonpositive_clock_or_nonisometric_response_rejected(clock, embedding):
    with pytest.raises(ValueError):
        closure.finite_defect(clock, embedding)


def test_floating_point_or_zero_resolution_does_not_masquerade_as_exact_proof():
    with pytest.raises(TypeError):
        closure.finite_defect([[1.0]], [[1]])
    with pytest.raises(TypeError):
        closure.finite_defect([[1]], [[1]], .1)
    with pytest.raises(ValueError):
        closure.finite_defect([[1]], [[1]], 0)


def test_original_score_error_retains_rare_success_denominator():
    w, e1, e2 = F(3, 1000), F(1, 10**7), F(1, 10**6)
    assert closure.scalar_score_error(w, e1, e2) == F(1, 2500)
    bound = closure.scalar_score_error(w, e1, e2)
    # Independent two-point records, including scores at both endpoints.
    for a in (F(0), w/3, w):
        b = a  # A Bernoulli bounded score has identical first/second moments.
        for d1 in (-e1, e1):
            for d2 in (-e2, e2):
                estimate = min(w, max(F(0), a+d1))
                actual_error = abs((b+d2)/w-(estimate/w)**2-(b/w-(a/w)**2))
                assert actual_error <= bound
    assert closure.scalar_score_error(w/10, e1, e2) == 10*bound


@pytest.mark.parametrize('args', [(0, 0, 0), (-1, 0, 0), (2, 0, 0), (1, -1, 0), (1, 0, -1)])
def test_invalid_score_domains_rejected(args):
    with pytest.raises(ValueError):
        closure.scalar_score_error(*args)


def test_native_dimer_uses_full_source_and_retains_all_original_probabilities():
    row = closure.native_dimer_control()
    assert row['source_spectrum'] == ['-9/8', '0', '0', '1/8']
    probabilities = list(map(F, row['complete_record_probabilities']))
    assert probabilities == [F(1, 2), F(7999, 48400), F(16201, 48400)]
    assert sum(probabilities) == 1
    weights = list(map(F, row['original_success_spectral_weights']))
    assert weights == [F(81, 400), F(16, 121)]
    assert sum(weights) == probabilities[-1] < 1
    assert F(row['normalized_rank_one_defect']) == F(5184, 262472401)
    assert row['two_component_response_defect'] == [['0', '0'], ['0', '0']]
    assert row['native_pulse_tape_executed'] is False
    assert row['continuum_Gram_projection_implemented'] is False
    assert row['physical_particle_mass'] is None


def test_dimer_weights_against_independent_double_precision_diagonalization():
    raw, scale = closure.source.matrix(2, ((0, 1),), F(1), F(3, 8))
    energies, vectors = np.linalg.eigh(raw/scale)
    ground = vectors[:, 0]
    x = np.eye(4)[[1, 0, 3, 2]]
    z = np.diag([1, -1, 1, -1])
    k = raw/scale-energies[0]*np.eye(4)
    filt = k@np.linalg.inv(k+np.eye(4)/8)
    success = filt@(x+z)@ground/2
    weights = abs(vectors.conj().T@success)**2
    assert abs(weights[0]) < 1e-28
    assert abs(weights[1]+weights[2]-81/400) < 1e-13
    assert abs(weights[3]-16/121) < 1e-13
    bounded = 1/(1+8*(energies-energies[0]))
    w = np.vdot(success, success).real
    defect = weights@(bounded*bounded)/w-(weights@bounded/w)**2
    assert abs(defect-5184/262472401) < 1e-14


def test_actual_source_clock_mutation_fails_spectral_anchor(monkeypatch):
    original = closure.source.matrix
    def altered(*args, **kwargs):
        raw, scale = original(*args, **kwargs)
        return raw+np.eye(len(raw), dtype=raw.dtype), scale
    monkeypatch.setattr(closure.source, 'matrix', altered)
    with pytest.raises(AssertionError):
        closure.native_dimer_control()


def test_small_defect_cannot_certify_uniform_time_evolution():
    rows = [closure.high_energy_control(n) for n in (1, 10, 100)]
    defects = [F(row['defect']) for row in rows]
    assert defects[2] < defects[1] < defects[0]
    assert all(row['original_probability_outside_response'] == '1' for row in rows)
    assert all(row['uniform_time_control_from_small_resolvent_defect'] is False for row in rows)
    # The obstruction persists at one fixed time on an odd-energy subsequence.
    v = s.Matrix([1, 1])/s.sqrt(2)
    for m in (1, 11, 101):
        u = s.diag(s.exp(-s.I*m*s.pi/2), s.exp(-3*s.I*m*s.pi/2))
        assert s.simplify((v.H*u*v)[0]) == 0


def test_two_pair_circle_preserves_scalar_twists_without_closing_visible_force():
    row = closure.two_pair_circle_control()
    assert row['visible_regular_CCR_pairs'] == 2
    assert row['full_internal_U1_preserves_clock_and_vacuum'] is True
    assert row['two_speeds_squared'] == ['6', '18']
    assert row['hidden_force_coefficient'] == '-6'
    assert row['inherited_single_pair_control']['Q_first_energy_moment'] == '3/2'
    assert row['inherited_single_pair_control']['P_first_energy_coefficient_over_momentum_squared'] == '1/8'
    assert row['projected_memory_proves_scattering'] is False
    assert row['actual_OPH_source_limit_claimed'] is False
    assert row['source_specific_force_estimate_replaced'] is False


def test_countercontrol_replays_inherited_science_not_an_unchecked_copy(monkeypatch):
    original = closure.asymptotic.two_speed_control
    def altered():
        row = original()
        row['order'] = '1/8'
        return row
    monkeypatch.setattr(closure.asymptotic, 'two_speed_control', altered)
    with pytest.raises(AssertionError):
        closure.two_pair_circle_control()


def test_certificate_excludes_automatic_continuum_projection_or_dispersion():
    result = closure.certificate()
    assert result['arbitrary_response_projection_implemented'] is False
    assert result['actual_continuum_response_defect_proved_zero'] is False
    assert result['actual_source_dispersion_identified'] is False


@pytest.mark.parametrize('m,kappa', [(F(1, 4), F(1, 4)), (F(1, 2), F(1, 8)),
                                    (F(11, 100), F(1, 2)), (F(2), F(3))])
def test_general_control_matches_both_measured_coefficients_without_branch_collapse(m, kappa):
    row = closure.general_two_speed_coefficients(m, kappa)
    a = F(row['kinetic_coefficient'])
    b1, b2 = map(F, row['gradient_coefficients'])
    assert m*a/2 == 6*kappa == F(row['Q_first_energy_moment'])
    assert m*(b1+b2)/4 == kappa/2 == F(row['P_first_energy_coefficient_over_momentum_squared'])
    speeds = list(map(F, row['two_speeds_squared']))
    assert speeds == [a*b1, a*b2] and speeds[1] == 3*speeds[0]
    assert F(row['hidden_Delta_coefficient']) == a*(b1-b2)/2 != 0
    assert row['source_coefficients_selected_by_this_control'] is False
    assert row['actual_OPH_source_limit_claimed'] is False
    if F(1, 10) < m <= F(1, 2) and m*m/2 <= kappa <= F(1, 2):
        assert row['within_stated_joint_source_parameter_range'] is True
        assert F(row['sufficient_l1_cone_speed_squared']) < 72**2
        assert F(row['Q_inverse_energy_coefficient_over_inverse_momentum_squared']) <= 4
        assert F(row['P_inverse_energy_moment']) <= F(1, 12)


@pytest.mark.parametrize('args', [(0, 1), (1, 0), (-1, 1), (1, -1)])
def test_nonpositive_general_stiffness_rejected(args):
    with pytest.raises(ValueError):
        closure.general_two_speed_coefficients(*args)
