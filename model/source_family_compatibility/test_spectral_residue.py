"""Independent spectral controls and actual-source range/moment checks."""
from fractions import Fraction as F

import numpy as np
import pytest
import sympy as s

from . import spectral_residue as m


def test_noncommuting_exact_operator_comparison():
    result = m.noncommuting_control()
    assert result['exact_noncommuting_real_imaginary_memory']
    assert result['matrix_lower_and_linear_upper_PSD']
    assert not result['is_actual_many_body_source_realization']


def test_noncommuting_inequalities_against_full_independent_resolvents():
    rng = np.random.default_rng(793)
    for _ in range(32):
        r = rng.normal(size=(3, 2))+1j*rng.normal(size=(3, 2))
        r /= 2*np.linalg.norm(r, 2)
        d = rng.normal(size=(3, 3))+1j*rng.normal(size=(3, 3))
        d = (d+d.conj().T)/2
        eta = float(rng.uniform(.1, 2))
        alpha = float(rng.uniform(.5, 4))
        L = np.block([[np.zeros((2, 2)), alpha*np.eye(2), np.zeros((2, 3))],
                      [alpha*np.eye(2), np.zeros((2, 2)), r.conj().T],
                      [np.zeros((3, 2)), r, d]])
        g = np.linalg.inv(1j*eta*np.eye(7)-L)[:2, :2]
        q = -eta*(g-g.conj().T)/(2j)
        memory = r.conj().T@np.linalg.solve(d@d+eta*eta*np.eye(3), eta*eta*r)
        beta = alpha*alpha+eta*eta
        r2 = np.linalg.norm(r, 2)**2
        lower = (eta*eta*np.eye(2)+memory)@np.linalg.inv(beta*np.eye(2)+memory)
        upper = eta*eta/beta*np.eye(2)+alpha*alpha*(beta+r2)/beta**3*memory
        assert np.linalg.eigvalsh(q-lower).min() > -2e-13
        assert np.linalg.eigvalsh(upper-q).min() > -2e-13


def test_kernel_formula_handles_complex_matrix_weight():
    result = m.kernel_control()
    assert result['full_kernel_dimension'] == result['hidden_kernel_dimension'] == 2
    assert result['kernel_graph_identity']
    r = s.Matrix([[1, s.I], [1, 1], [2, -s.I]])/3
    actual = m._kernel_projection(m.block(4, r, s.diag(0, 0, 2)))[:2, :2]
    wrong = (r.H*r)*(16*s.eye(2)+r.H*r).inv()
    assert m._simple(actual-wrong) != s.zeros(2)
    assert result['physical_mass_assigned'] is None


def test_hidden_kernel_without_source_overlap_is_dark():
    r = s.Matrix([[0], [1]])
    L = m.block(4, r, s.diag(0, 2))
    assert len(L.nullspace()) == 1
    assert m._kernel_projection(L)[0, 0] == 0


def test_six_moments_and_energy_reflection_do_not_determine_zero_weight():
    result = m.matched_moment_controls()
    first, second = result['controls']
    common = [1, 0, 16, 0, 272, 0, 4640]
    assert first['original_moments_0_through_8'][:7] == common
    assert second['original_moments_0_through_8'][:7] == common
    assert first['original_zero_weight'] == 0
    assert second['original_zero_weight'] == F(3, 67)
    assert first['original_moments_0_through_8'][8] == 79168
    assert second['original_moments_0_through_8'][8] == 79216
    assert first['characteristic_coefficients'] == [1, 0, -18, 0, 16]
    assert second['characteristic_coefficients'] == [1, 0, -21, 0, 67, 0]
    assert not result['actual_Hubbard_source_realizations']


def test_finite_broadening_limits_match_both_exact_controls():
    for weights, residue in [([4, 1, 1], 0), ([4, 1, 1, np.sqrt(3)], 3/67)]:
        L = np.diag(weights, 1)+np.diag(weights, -1)
        eig, vec = np.linalg.eigh(L)
        p = np.abs(vec[0])**2
        eta = 1e-5
        q = np.sum(p*eta**2/(eig**2+eta**2))
        assert abs(q-residue) < 2e-10


def test_actual_complete_source_moment_range_and_bloch_gluing():
    result = m.source_locality_control()
    assert result['whole_Fock_dimension'] == 64
    assert result['exact_outside_range_zero_checks'] == 16
    assert result['range_two_nonzero_checked']
    assert result['three_actual_Bloch_boundary_gluings_checked']
    assert not result['uniform_low_frequency_modulus_inferred']


def test_physical_cell_cos_symbol_requires_probe_gluing():
    gammas = m.hole_dispersion.gamma_matrices()
    from itertools import product
    parities = list(product(range(2), repeat=3))
    k = np.array([.23, -.7, .6])
    symbol = sum(np.cos(k[j])*gammas[j] for j in range(3))
    shifted = sum(np.cos(k[j]+(np.pi if j == 0 else 0))*gammas[j] for j in range(3))
    glue = np.diag([(-1)**x[0] for x in parities])
    assert np.linalg.norm(symbol-shifted, 2) > 1
    assert np.linalg.norm(glue@symbol@glue-shifted, 2) < 1e-14


def test_third_derivative_block_identity_before_estimation():
    h = s.Matrix([[1, s.I], [-s.I, -1]])/5
    c = -2*h
    r = s.Matrix([[1, s.I], [2, 1], [0, 1-s.I]])/3
    d = s.Matrix([[0, 1, s.I], [1, 2, 1-s.I], [-s.I, 1+s.I, -1]])
    alpha = s.Integer(4)
    L = m.block(alpha, r, d, h, c)
    original = s.eye(7)[:, :2]
    predicted = s.Matrix.vstack(h**3+2*alpha**2*h+alpha**2*c,
        alpha*(h*h+alpha**2*s.eye(2)+c*(h+c)+r.H*r),
        alpha*(r*(h+c)+d*r))
    assert m._simple(L**3*original-predicted) == s.zeros(7, 2)


def test_actual_subcontact_floor_consumes_source_residual_and_sixth_moment():
    result = m.subcontact_band()
    parent = m.nodal_memory.default_budget(F(1))
    assert result['square_energy_variance_lower'] == 16*parent['neighbourhood_Gram_lower']
    assert result['unit_hopping_sixth_moment_upper'] == 431850305
    assert result['third_absolute_centered_square_energy_moment_upper'] == 431912440
    assert result['nodal_region_S_upper'] == F(1296, 843578984375)
    assert result['original_combined_charge_probability_lower'] == F(2592, 14340842734375)
    assert result['absolute_energy_squared_over_c_squared_upper'] < 16
    assert not result['positive_lower_energy_edge_proved']
    assert not result['physical_ballistic_window_proved']
    assert not result['original_charge_probabilities_conditioned']


def test_negative_deviation_and_window_steps_on_independent_distributions():
    for values, weights in [([0, 16, 36], [F(1, 10), F(4, 5), F(1, 10)]),
                            ([15, 17], [F(1, 2), F(1, 2)]),
                            ([0, 20], [F(1, 5), F(4, 5)])]:
        mean = sum(w*x for w, x in zip(weights, values, strict=True))
        variance = sum(w*(x-mean)**2 for w, x in zip(weights, values, strict=True))
        third = sum(w*abs(x-mean)**3 for w, x in zip(weights, values, strict=True))
        negative = sum(w*max(F(), mean-x) for w, x in zip(weights, values, strict=True))
        assert negative >= variance**2/(2*third)
        floor = variance**2/(2*third)
        delta = floor/4
        event = sum(w for w, x in zip(weights, values, strict=True) if x <= mean-2*delta)
        assert event >= floor/(2*mean)
    # A contact atom has zero variance: its second moment alone gives no
    # strictly subcontact probability.
    assert 16**2-16**2 == 0


def test_inherited_hidden_memory_floor_vanishes_instead_of_selecting_an_atom():
    values = [m.default_memory_floor(F(1, 2**n)) for n in (0, 4, 8)]
    scores = [x['hidden_Lorentzian_score_lower'] for x in values]
    assert scores[0] > scores[1] > scores[2] > 0
    for n, result in zip((0, 4, 8), values, strict=True):
        eta = F(1, 2**n)
        assert result['hidden_Lorentzian_score_lower'] <= eta**2*F(18, 125)**2/5175**2
        assert not result['hidden_score_has_proved_nonzero_limit']


def test_joint_broadening_keeps_constant_kinetic_error():
    for a in (F(1, 2), F(1, 16), F(1, 256)):
        result = m.score_bounds(4, 2*a, 0, 18, 18, a, 3*a)
        assert result['kinetic_comparison_error'] == F(3, 2)
        assert result['original_score_lower'] == 0
        assert result['original_score_upper'] == 1
        assert not result['shrinking_window_weight_assumed']


def test_original_window_conversion_retains_mass_without_conditioning():
    result = m.window_from_score(F(3, 4), 1)
    assert result['original_window_probability_lower'] == F(1, 2)
    assert not result['success_conditioning_used']
    assert m.window_from_score(F(1, 4), 1)['original_window_probability_lower'] == 0
    # Source atom at zero remains part of the admitted original event.
    assert m.window_from_score(1, 2)['original_window_probability_lower'] == 1
    energies, probabilities = [F(), F(1), F(4)], [F(1, 5), F(3, 10), F(1, 2)]
    score = sum(p/(1+e*e) for p, e in zip(probabilities, energies, strict=True))
    assert m.window_from_score(score, 1)['original_window_probability_lower'] <= F(1, 2)


@pytest.mark.parametrize('arguments', [(0, 1, 0, 1, 1), (4, 0, 0, 1, 1),
    (4, 1, 2, 1, 2), (4, 1, 0, 3, 2), (4., 1, 0, 1, 1), (True, 1, 0, 1, 1)])
def test_exact_score_inputs_reject_invalid_premises(arguments):
    with pytest.raises(ValueError):
        m.score_bounds(*arguments)


def test_block_validation_does_not_silently_symmetrize_a_clock():
    with pytest.raises(ValueError):
        m.block(4, [[1], [0]], [[0, 1], [0, 0]])
    with pytest.raises(ValueError):
        m.block(4, [[1]], [[0]], h=[[s.I]])
    with pytest.raises(ValueError):
        m.window_from_score(F(5, 4), 1)
