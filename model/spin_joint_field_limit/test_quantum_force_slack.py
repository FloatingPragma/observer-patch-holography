"""Independent full-spin and variational controls for the force mismatch."""
from fractions import Fraction as F

import pytest
import sympy as s

from model.spin_joint_field_limit import quantum_force_slack as q


@pytest.mark.parametrize('weights', [(1, 2, 3, 4), (0, 1, 1, 8), (2, 0, 0, 7), (0, 9, 1, 0)])
def test_full_density_matrix_correlations(weights):
    # Mixtures of up-up/down-down, two product states, and the singlet.
    corner, updown, downup, singlet = map(s.Rational, weights)
    norm = corner+updown+downup+singlet
    a = corner/(2*norm)
    b = (updown+singlet/2)/norm
    c = (downup+singlet/2)/norm
    z = -singlet/(2*norm)
    rho = s.Matrix([[a, 0, 0, 0], [0, b, z, 0], [0, z, c, 0], [0, 0, 0, a]])
    assert s.trace(rho) == 1
    assert all(value >= 0 for value in rho.eigenvals())
    spin = (s.Matrix([[0, 1], [1, 0]])/2,
            s.Matrix([[0, -s.I], [s.I, 0]])/2, s.diag(1, -1)/2)
    mean = s.trace(rho*s.kronecker_product(spin[2], s.eye(2)))
    cov = [s.trace(rho*s.kronecker_product(axis, axis)) for axis in spin]
    d = -a-z
    if d < 0:
        with pytest.raises(ValueError):
            q.marginal(*(F(x) for x in (a, b, c, z)))
        return
    got = q.marginal(*(F(x) for x in (a, b, c, z, d)))
    assert got['magnetization'] == mean
    assert got['bond_coefficient'] == -cov[2]-cov[0]
    assert got['spin_dot_energy_coefficient'] == -3*sum(cov)
    assert got['force_coefficient_excess'] >= got['excess_lower']


def test_sharp_local_boundary_and_classical_control():
    boundary = q.marginal(0, F(9, 10), F(1, 10), F(-3, 10), F(3, 10))
    assert boundary['force_coefficient_excess'] == F(6, 25)
    assert boundary['force_coefficient_excess'] == boundary['excess_lower']
    assert all(x == 0 for x in boundary['nonnegative_remainder_terms'])
    classical = q.marginal(0, 1, 0, 0)
    assert classical['force_coefficient_excess'] == 0
    assert classical['spin_dot_energy_coefficient'] == F(3, 4)
    with pytest.raises(ValueError):
        q.marginal(0, 1, 0, 0, F(1, 120))


@pytest.mark.parametrize('args', [
    (0, F(9, 10), F(1, 10), F(-1, 2)),  # negative central determinant
    (F(-1, 10), F(3, 5), F(3, 5), 0),
    (0, F(1, 2), F(1, 3), 0),
    (0, F(9, 10), F(1, 10), F(-3, 10), F(31, 100)),
])
def test_rejects_invalid_or_forged_marginal_premise(args):
    with pytest.raises(ValueError):
        q.marginal(*args)


@pytest.mark.parametrize('n', [2, 3, 4, 16, 101])
def test_actual_joint_trajectory_and_variational_margin(n):
    got = q.joint_floor(n)
    assert got['p'] == F(1, n*n)
    assert got['B_over_J']**2 == got['p']**5
    delta = got['exchange_energy_e_lower']/3-F(1, 4)
    assert delta == got['marginal_energy_improvement_lower']
    assert delta >= F(13, 1600)
    assert got['kappa_minus_m_squared_lower'] == delta**2+delta/2
    assert F(6, 25) > got['dimer_trial_order_lower']**2
    assert got['dimer_trial_order_lower'] > 0


@pytest.mark.parametrize('n', [0, 1, -2, True, F(2), 2.0])
def test_joint_trajectory_domain_is_not_relaxed(n):
    with pytest.raises(ValueError):
        q.joint_floor(n)


def test_matrix_source_force_controls_and_no_continuum_promotion():
    got = q.certificate()
    source = got['actual_dimer_force_control']
    assert source['ground_energy'] == '-4/3'
    assert source['Q_residual_ground_score'] == '5/192'
    assert source['P_residual_ground_score'] == '59/300'
    assert source['source_degree'] == 1
    assert not source['finite_dimer_substituted_for_cubic_ground']
    assert got['liminf_kappa_minus_m_squared_lower'] == F(61, 14400)
    for name in ('corrected_kinetic_coefficients_excluded', 'infrared_sharp_component_excluded',
                 'source_to_history_unbounded_force_norm_preservation_assumed',
                 'physical_mass_or_species_assigned', 'additional_native_detector_execution'):
        assert got[name] is False


def test_corrupted_actual_source_is_rejected(monkeypatch):
    actual = q.source.matrix
    def changed(*args, **kwargs):
        matrix, scale = actual(*args, **kwargs)
        matrix = matrix.copy()
        matrix[0, 0] += 1
        return matrix, scale
    monkeypatch.setattr(q.source, 'matrix', changed)
    with pytest.raises(ArithmeticError, match='source frame'):
        q.certificate()


def test_returned_rationals_do_not_mutate_inputs_or_later_certificates():
    value = F(3, 10)
    row = q.marginal(0, F(9, 10), F(1, 10), -value, value)
    snapshot = q.certificate()['liminf_kappa_minus_m_squared_lower']
    row['energy_improvement']._numerator = 999
    assert value == F(3, 10)
    changed = q.certificate()
    changed['liminf_kappa_minus_m_squared_lower']._numerator = 999
    assert q.certificate()['liminf_kappa_minus_m_squared_lower'] == snapshot


def test_exact_rational_inputs_required():
    with pytest.raises(TypeError):
        q.marginal(0, .9, .1, -.3)


def test_owner_binds_marginal_proof_source_trial_and_scalar_Lean():
    from model.spin_joint_field_limit import packet
    paths = packet.bindings()
    assert all(name in paths for name in (
        'model/spin_joint_field_limit/quantum_force_slack.py',
        'model/spin_joint_field_limit/QUANTUM_FORCE_SLACK.md',
        'model/spin_charge_dispersion/dimer.py',
        'model/spin_charge_dispersion/DIMER_BOUND.md',
        'Lean/SourceLaw/SpinMarginalForce.lean'))


@pytest.mark.parametrize('flag', ('corrected_kinetic_coefficients_excluded',
    'infrared_sharp_component_excluded', 'physical_mass_or_species_assigned',
    'source_to_history_unbounded_force_norm_preservation_assumed'))
def test_owner_rejects_forged_force_promotion(monkeypatch, flag):
    from model.spin_joint_field_limit import packet
    import json
    true = {'source_quantum_correction_to_mean_spin_force': packet.serialize(q.certificate())}
    forged = json.loads(json.dumps(true))
    forged['source_quantum_correction_to_mean_spin_force'][flag] = True
    monkeypatch.setattr(packet, 'build', lambda: true)
    with pytest.raises(ValueError, match='scope differs'):
        packet.validate(forged)
