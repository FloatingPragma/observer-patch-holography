"""Actual translations, original outcomes and integer-certified source width."""
from fractions import Fraction as F

import numpy as np
import pytest

from . import directed_bands as d


@pytest.mark.parametrize('delta', [F(1, 4), F(1, 10), F(1, 1000)])
def test_arbitrarily_narrow_caps_keep_original_probability(delta):
    row = d.directed_band(delta, F(2, 7))
    assert row['each_charge_and_direction_original_probability_lower'] == F(1, 774400)
    lo, hi = row['momentum_norm_interval']
    r = row['annular_parameter']
    assert r/2 <= lo < hi <= r
    assert row['cap_radius'] == delta*F(2, 7)
    assert len(row['complete_records']) == 7
    assert row['profile_changes_when_cap_shrinks']
    assert row['bounded_creator_has_no_assumed_momentum_support']
    assert not row['energy_shell_or_plane_wave_atom_inferred']


@pytest.mark.parametrize('delta,rho', [(0, 1), (F(1, 3), 1), (0.1, 1),
                                      (True, 1), (F(1, 4), 0), (F(1, 4), 0.2)])
def test_invalid_or_inexact_caps_rejected(delta, rho):
    with pytest.raises(ValueError):
        d.directed_band(delta, rho)


def test_spatial_phase_precision_preserves_original_weight():
    first = d.translation_score(F(1, 10), 3, F(1, 1000))
    second = d.translation_score(F(1, 20), 3, F(1, 1000))
    assert first['original_interference_score_defect_upper'] == F(9, 200000)
    assert second['original_interference_score_defect_upper']*4 == first['original_interference_score_defect_upper']
    assert not first['success_probability_divided_out']
    assert not first['time_phase_or_dispersion_inferred']
    with pytest.raises(ValueError):
        d.translation_score(1, 1, F(1, 3))


def _projections(t):
    size = len(t)
    return [(np.eye(size)+np.exp(-2j*np.pi*l/3)*t
             +np.exp(-4j*np.pi*l/3)*(t @ t))/3 for l in (0, 1, -1)]


def test_actual_order_three_translation_and_source_charge():
    row = d.finite_operators()
    t, h, q = row['even_translation'], row['H_integer'], row['twice_charge']
    assert np.array_equal(t @ t @ t, np.eye(64, dtype=int))
    assert not np.array_equal(t @ t, np.eye(64, dtype=int))
    assert np.array_equal(t @ h, h @ t)
    assert np.array_equal(t @ q, q @ t)
    assert np.array_equal(row['charge_pairing'] @ q @ row['charge_pairing'], -q)
    projectors = _projections(t)
    assert np.max(np.abs(sum(projectors)-np.eye(64))) < 1e-14
    for i, p in enumerate(projectors):
        for j, r in enumerate(projectors):
            assert np.max(np.abs(p @ r-(p if i == j else 0))) < 1e-14


def test_native_source_and_probes_share_the_same_Marshall_frame():
    row = d.finite_operators()
    identity = np.eye(64)
    words = np.arange(64)
    x = [identity[words ^ (1 << site)] for site in range(6)]
    z = [np.diag(1-2*((words >> site) & 1)) for site in range(6)]
    y = [1j*x[site] @ z[site] for site in range(6)]
    epsilon = [(-1)**site for site in range(6)]
    marshall = np.diag([(-1)**sum((int(w) >> s) & 1 for s in (1, 3, 5)) for w in words])
    physical = sum((x[a] @ x[b]+y[a] @ y[b]+z[a] @ z[b]-identity
                    for a, b in row['edges']), np.zeros((64, 64)))/4
    physical -= sum((epsilon[s]*x[s] for s in range(6)), np.zeros((64, 64)))/64
    assert np.array_equal(marshall @ physical @ marshall, row['H_integer']/row['H_scale'])
    profile = [1, -F(1, 2), -F(1, 2), 1, -F(1, 2), -F(1, 2)]
    for q in (-1, 1):
        creator = sum((epsilon[s]*float(profile[s])*(z[s]-q*1j*y[s])
                       for s in range(6)), np.zeros((64, 64)))/16
        assert np.array_equal(marshall @ creator @ marshall, row['creator_numerators'][q]/32)


def test_full_thirteen_records_are_complete_before_conditioning():
    row = d.finite_operators()
    h = row['H_integer']/row['H_scale']
    values, vectors = np.linalg.eigh(h)
    ground = vectors[:, 0]
    ground_projection = np.outer(ground, ground)
    energies = [np.eye(64)-ground_projection, ground_projection]
    effects, probabilities = [], []
    for q in (-1, 1):
        b = row['creator_numerators'][q]/32
        for p in _projections(row['even_translation']):
            for energy in energies:
                effect = b.T @ p @ energy @ p @ b
                assert np.max(np.abs(effect-effect.conj().T)) < 1e-14
                effects.append(effect)
                probabilities.append(float(np.real(ground @ effect @ ground)))
    failure = row['complement_numerator']/1024
    effects.append(failure)
    probabilities.append(float(ground @ failure @ ground))
    assert len(effects) == 13
    assert np.max(np.abs(sum(effects)-np.eye(64))) < 1e-13
    assert abs(sum(probabilities)-1) < 1e-13
    successes = [p for p in probabilities[:-1] if p > 1e-12]
    assert len(successes) == 4
    assert max(successes)-min(successes) < 1e-14
    assert min(successes) > 1/208
    assert sum(abs(p) < 1e-12 for p in probabilities) == 8
    scaled = (values[1:]-values[0])*4
    assert min(scaled) > 8/3 and max(scaled) < 18


def test_exact_interval_admission_and_nonzero_directional_energy_width():
    row = d.certificate()['finite_periodic_source_control']
    assert row['complete_record_count'] == 13
    assert row['relative_energy_variance_lower'] == F(1, 800)
    assert row['certified_intervals']['scaled_energy_variance'][0] > F(1, 21)
    assert row['ground_vector_error_upper'] == F(1, 10**8)
    assert not row['floating_solver_trusted']
    assert not row['finite_control_proves_continuum_floor']
    assert not row['finite_translation_labels_identified_with_continuum_momentum']
    assert not row['physical_particle_or_mass_assignment']


def test_incomplete_or_falsified_ground_proposals_are_rejected():
    row = d.finite_operators()
    h, scale = row['H_integer'], row['H_scale']
    e, v = d.ground_proposal(h, scale)
    with pytest.raises(ValueError):
        d.verify_ground(h, scale, e[:-1], v[:, :-1])
    broken = v.copy()
    broken[:, 1] = broken[:, 0]
    with pytest.raises(ArithmeticError):
        d.verify_ground(h, scale, e, broken)
    shifted = e.copy()
    shifted[0] -= 1 << 40
    with pytest.raises(ArithmeticError):
        d.verify_ground(h, scale, shifted, v)


def test_source_moment_interval_matches_independent_complex_pvm():
    row = d.finite_operators()
    h = row['H_integer']/row['H_scale']
    energies, basis = np.linalg.eigh(h)
    ground = basis[:, 0]
    image = _projections(row['even_translation'])[1] @ (row['creator_numerators'][1] @ ground/32)
    weights = np.abs(basis.T @ image)**2
    original = float(sum(weights))
    mean = float(sum(weights*(energies-energies[0])*4)/original)
    variance = float(sum(weights*((energies-energies[0])*4-mean)**2)/original)
    checked = d.finite_control()['actual_intervals']
    for key, actual in [('each_original_charge_direction_probability', original),
                        ('scaled_energy_mean', mean), ('scaled_energy_variance', variance)]:
        lo, hi = checked[key]
        assert float(lo) <= actual <= float(hi)


def test_public_certificate_does_not_share_mutable_rationals():
    first = d.certificate()
    before = repr(d.certificate())
    first['finite_periodic_source_control']['certified_intervals']['scaled_energy_mean'][0]._numerator = 0
    assert repr(d.certificate()) == before
