"""Exact root/gap/window admissions and independent full-source comparisons."""
from fractions import Fraction as F
import json

import numpy as np
import pytest
from scipy.linalg import eigh
import sympy as s

from . import residue as r


def array(matrix):
    return np.asarray(matrix.tolist(), dtype=complex)


def test_sturm_counts_restore_repeated_eigenvalue_multiplicity():
    matrix = ((F(0), F(0), F(0)), (F(0), F(0), F(0)), (F(0), F(0), F(1)))
    assert r.root_count(matrix, -1, F(1, 2)) == 2
    assert r.root_count(matrix) == 3
    assert r.root_count(matrix, F(1, 2), 2) == 1
    with pytest.raises(ValueError, match='off the spectrum'):
        r.root_count(matrix, 0, 2)
    with pytest.raises(ValueError, match='increasing'):
        r.root_count(matrix, 2, 1)


def test_full_source_ground_and_rank_two_charged_windows_are_isolated():
    data = r.spectral_windows()
    assert data['ground_rank'] == 1
    h = array(r.inherited.source.exact_control()['H'])
    energy = eigh(h, eigvals_only=True)
    assert sum(energy < -2) == 1
    # The48-bit exact enclosure is narrower than the dense binary64 solver's
    # backward error; this independent check is only a numerical comparison.
    assert float(data['ground_interval'][0])-1e-12 <= energy[0] <= float(data['ground_interval'][1])+1e-12
    for n, row in data['charged'].items():
        indices = [i for i in range(208) if (i//13).bit_count() == n]
        values = eigh(h[np.ix_(indices, indices)], eigvals_only=True)
        assert sum(values < -.6) == 2
        lo, hi = map(float, row['energy_interval'])
        assert np.max(np.abs(values[:2]-(lo+hi)/2)) < 1e-13
        assert row['projector_rank'] == 2


def test_rational_trial_has_a_certified_global_ground_error():
    data = r.ground_trial()
    assert data['occupation'] == 2
    assert data['half_trace_distance_upper'] < F(1065, 10**10)
    assert data['half_trace_distance_upper']**2 >= data['ground_leakage_upper']
    assert r.ground_trial(3)['ground_leakage_upper'] > data['ground_leakage_upper']
    assert not data['source_density_exactly_stationary']
    v = data['vector']; h = s.Matrix(r.spectral_windows()['ground_block_matrix'])
    assert F((v.T*h*v)[0])/data['norm_squared'] == data['mean_energy']
    assert data['mean_energy'] >= r.spectral_windows()['ground_interval'][0]
    with pytest.raises(ValueError):
        r.ground_trial(True)


def test_every_rational_collision_has_complete_effects_and_unitary_dilation():
    data = r.spectral_windows()
    unique = set(matrix for _, n, matrix in r.inherited.rational_blocks() if n == 1)
    lower = data['charged'][1]['energy_interval'][0]
    for matrix in unique:
        collision = r.collision_block(matrix, lower)
        d, t = collision['D'], collision['success_effect']
        eye = s.eye(d.rows)
        assert d == d.T and t == t.T and d*t == t*d
        assert t*(eye+d*d) == eye
        # These exact real identities imply the displayed rational unitary
        # for every source block, avoiding numerical square-root admission.
        assert t*t+(d*t)*(d*t) == t
        assert (eye-t)*(eye-t)+(d*t)*(d*t) == eye-t
        summed = t**r.FILTER_POWER+sum((t**j*(eye-t) for j in range(r.FILTER_POWER)), s.zeros(d.rows))
        assert summed == eye
    # Direct full complex dilation check on a nonscalar Hermitian input.
    collision = r.collision_block(((F(2), F(1)), (F(1), F(3))), F(0))
    a, b = collision['A'], collision['B']
    u = a.row_join(-b.H).col_join(b.row_join(a.H))
    assert (u.H*u).applyfunc(s.simplify) == s.eye(4)


def test_all_failure_times_are_original_outcomes_with_exact_half_sign_weight():
    detector = r.complete_window_detector()
    rows = detector['records']
    assert len(rows) == 40
    assert {row['outcome'] for row in rows} == {'window', 'failure_0', 'failure_1', 'failure_2', 'failure_3'}
    assert sum(row['original_probability'] for row in rows) == 1
    assert all(0 <= row['original_probability'] <= 1 for row in rows)
    for sign in (1, -1):
        assert sum(row['original_probability'] for row in rows if row['sign'] == sign) == F(1, 2)
    assert any(row['original_probability'] > 0 for row in rows if row['outcome'] == 'failure_3')
    assert not detector['conditioned_on_success']


def test_complete_records_match_independent_dense_source_detector():
    h = array(r.inherited.source.exact_control()['H'])
    trial = r.ground_trial()
    vector = np.zeros(208, dtype=complex)
    phases = r.inherited.source.real_computation_chart()[0]
    vector[list(trial['indices'])] = np.array([complex(phases[i])*float(v)
        for i, v in zip(trial['indices'], trial['vector'], strict=True)]) / np.sqrt(float(trial['norm_squared']))
    expected = {}
    for row in r.inherited.instrument():
        n = 2+row['sign']
        indices = [i for i in range(208) if (i//13).bit_count() == n]
        y = (array(row['K'])@vector)[indices]
        target = h[np.ix_(indices, indices)]
        lo = float(r.spectral_windows()['charged'][n]['energy_interval'][0])
        d = (target-lo*np.eye(len(indices)))/float(r.FILTER_WIDTH)
        effect = np.linalg.inv(np.eye(len(indices))+d@d)
        power = np.eye(len(indices), dtype=complex)
        for j in range(r.FILTER_POWER):
            expected[row['sign'], row['mode'], f'failure_{j}'] = float(np.vdot(y, power@(np.eye(len(indices))-effect)@y).real)
            power = power@effect
        expected[row['sign'], row['mode'], 'window'] = float(np.vdot(y, power@y).real)
    for row in r.complete_window_detector()['records']:
        assert abs(float(row['original_probability'])-expected[row['sign'], row['mode'], row['outcome']]) < 2e-12


def test_certified_residue_encloses_independent_exact_ground_spectral_weight():
    h = array(r.inherited.source.exact_control()['H'])
    _, vectors = eigh(h)
    ground = vectors[:, 0]
    result = r.certified_residues()
    actual = {}
    for sign in (1, -1):
        indices = [i for i in range(208) if (i//13).bit_count() == 2+sign]
        _, target = eigh(h[np.ix_(indices, indices)])
        lowest = target[:, :2]
        weight = 0.
        for row in r.inherited.instrument():
            if row['sign'] == sign:
                value = (array(row['K'])@ground)[indices]
                weight += float(np.linalg.norm(lowest.conj().T@value)**2)
        actual[sign] = weight
        low, high = result['channels'][sign]['original_ground_residue_interval']
        assert float(low) < weight < float(high)
        assert low > F(498557, 1000000)
    assert abs(actual[1]-actual[-1]) < 1e-12
    assert result['exact_charge_conjugate_ground_weights_equal']


def test_filter_bounds_follow_from_exact_separated_windows():
    result = r.certified_residues()
    for sign in (1, -1):
        row = result['channels'][sign]
        assert row['outside_window_success_upper'] < F(1, 10**6)
        assert row['inside_window_failure_upper'] < F(1, 10**20)
        assert row['state_transfer_probability_error_upper'] == r.ground_trial()['half_trace_distance_upper']/2
        assert row['original_charge_branch_weight'] == F(1, 2)


def test_finite_clock_window_plan_retains_readout_error_and_original_weights():
    plan = r.finite_clock_window_plan()
    assert plan['alias_probability_upper'] == 0
    assert 0 < plan['filter_leakage_upper'] < F(3, 100)
    assert 0 < plan['complete_clock_error_upper'] < F(1, 10000)
    assert F(48, 100) < plan['original_readout_weight_per_sign_interval'][0] < F(1, 2)
    assert plan['original_readout_weight_per_sign_interval'][1] <= F(1, 2)
    assert plan['maximum_source_time'] > 8000
    assert not plan['full_clock_histogram_executed']
    with pytest.raises(ValueError):
        r.finite_clock_window_plan(radius=F(1, 128))


def test_exact_receipt_pins_domains_and_never_uses_floating_numbers():
    value = r.certificate()
    assert value['source_sha256'] == r.authenticate_sources()
    assert value['source_dimension'] == 208
    assert not value['physical_particle_identification']
    assert value['GeV_mass'] is None
    def no_float(item):
        assert not isinstance(item, float)
        if isinstance(item, dict):
            for child in item.values(): no_float(child)
        if isinstance(item, (tuple, list)):
            for child in item: no_float(child)
    no_float(value)
    assert json.loads(json.dumps(value)) == value
