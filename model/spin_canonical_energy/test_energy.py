"""Exact identities, positive spectral margins and non-tautological controls."""
import copy
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm

from . import algebra, bounds, packet


def test_exact_all_axis_spin_and_rotation_identities():
    report = algebra.exact_checks()
    assert report['all_three_double_commutators']
    assert report['exact_single_axis_rotation_energy']
    assert report['mixed_quadratic_rotation_term_zero']
    assert report['transverse_bond_eigenvalues'] == {'-1/2': 1, '0': 2, '1/2': 1}


def test_returned_algebra_report_cannot_mutate_cache():
    report = algebra.exact_checks()
    report['transverse_bond_eigenvalues']['0'] = 100
    assert algebra.exact_checks()['transverse_bond_eigenvalues']['0'] == 2


def test_rotation_energy_against_noncommuting_full_matrices():
    # The singlet supplies C=-I/4. This checks rotations and orientation,
    # not a claim that the infinite ordered ground is a bond singlet.
    sx = np.array([[0, 1], [1, 0]], complex)/2
    sy = np.array([[0, -1j], [1j, 0]], complex)/2
    sz = np.diag([1., -1.])/2
    h = sum(np.kron(s, s) for s in (sx, sy, sz))
    g = np.array([0., 1., -1., 0.])/np.sqrt(2.)
    ax, ay, bx, by = .137, -.214, .091, -.073
    ux, uy = expm(1j*bx*sz)@expm(1j*ax*sy), expm(1j*by*sz)@expm(1j*ay*sy)
    v = np.kron(ux, uy)@g
    actual = np.vdot(v, h@v).real - np.vdot(g, h@g).real
    def rotation(u):
        return np.array([[2*np.trace(u.conj().T@s@u@t).real for t in (sx, sy, sz)]
                         for s in (sx, sy, sz)])
    delta = rotation(ux)-rotation(uy)
    predicted = np.trace(delta.T@delta)/8
    assert abs(actual-predicted) < 2e-15


@pytest.mark.parametrize('channel', ['P', 'Q'])
@pytest.mark.parametrize('radius', [F(1, 100), F(1), F(17, 3)])
def test_window_positive_margin_and_homogeneity(channel, radius):
    w = bounds.window(channel, F(1, 1000), radius)
    assert w['energy_lower'] > 0
    assert w['retained_weight_lower'] > 0
    assert w['total_weight_lower']-bounds.tail_bound(
        w['first_moment_upper'], w['inverse_moment_upper'],
        w['energy_lower'], w['energy_upper']) == w['retained_weight_lower']
    base = bounds.window(channel, F(1, 1000), 1)
    assert w['energy_upper'] == radius*base['energy_upper']
    assert w['retained_weight_lower'] == base['retained_weight_lower']*(radius if channel == 'P' else 1/radius)


def test_actual_conservative_numbers():
    p = bounds.window('P', F(1, 1000), 1)
    assert (p['energy_lower'], p['energy_upper'], p['retained_weight_lower']) == (
        F(1, 8000000), F(24000000), F(1, 48000000))
    assert bounds.stiffness_interval(F(1, 1000)) == (F(1, 2000000), F(1, 2))
    assert bounds.rotation_remainder_coefficient() == F(9, 4)


@pytest.mark.parametrize('bad', [True, .001, 0, F(-1, 3), F(2, 3), '1/1000'])
def test_reject_unadmitted_order_parameters(bad):
    with pytest.raises(ValueError):
        bounds.window('P', bad, 1)


def test_reject_bad_window_inputs():
    for channel, radius in [('S', 1), ('P', False), ('Q', 0), ('P', -.5)]:
        with pytest.raises(ValueError):
            bounds.window(channel, F(1, 1000), radius)


def test_energy_first_moment_is_not_continuous_under_weak_limit():
    previous = bounds.escaping_energy(2)
    current = bounds.escaping_energy(1024)
    assert current['mass'] == 1
    assert current['inverse'] < 1
    assert current['mass_at_escaping_energy'] < previous['mass_at_escaping_energy']
    assert current['first_moment_at_escaping_energy'] == 1
    assert current['first'] > F(199, 100)*current['weak_limit_first']


def test_atom_free_measure_obeys_the_coarse_moment_requirements():
    for m in (F(1, 1000), F(1, 2)):
        control = bounds.continuous_countercontrol(m, F(3, 7))
        assert control['atom_free']
        assert control['obeys_coarse_P_bounds']
        assert control['is_a_Heisenberg_state_construction'] is False


def test_scope_and_current_report_validate():
    report = packet.build()
    packet.validate(report)
    assert report['theorem']['actual_continuum_dynamics_nonstatic_and_gapless']
    assert report['theorem']['continuum_twist_energy_equality_proved'] is False


def test_boolean_scope_laundering_is_rejected():
    report = packet.build()
    report['theorem']['sharp_mass_shell_proved'] = 0
    with pytest.raises(ValueError):
        packet.validate(report)
    report['theorem']['sharp_mass_shell_proved'] = True
    with pytest.raises(ValueError):
        packet.validate(report)


def test_result_and_source_mutations_rejected():
    report = packet.build()
    changed = copy.deepcopy(report)
    changed['spectral_windows'][0]['retained_weight_lower'] = '1'
    with pytest.raises(ValueError):
        packet.validate(changed)
    changed = copy.deepcopy(report)
    changed['source_bindings'].pop(next(iter(changed['source_bindings'])))
    with pytest.raises(ValueError):
        packet.validate(changed)


def test_source_race_rejected(monkeypatch):
    real = packet.bindings
    count = 0
    def racing():
        nonlocal count
        result = real()
        count += 1
        if count > 1:
            result['model/__init__.py'] = 'changed'
        return result
    monkeypatch.setattr(packet, 'bindings', racing)
    with pytest.raises(RuntimeError, match='source changed'):
        packet.build()
