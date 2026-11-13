from fractions import Fraction as F
import numpy as np
from .clock import (STEP, PUMP, collective, coherence_band, evolve, perron_certificate,
                    complete_native_control, records, single_branch_records)


def test_complete_native_operator_and_battery_work_identity():
    assert complete_native_control(3)['complete_four_state_pair_dimension'] == 64


def test_every_band_matches_the_original_collision_with_all_states():
    cap = 4
    d = collective(cap)
    v = np.arange(1, cap + 2) + 1j * np.arange(cap + 1)
    v = v / np.linalg.norm(v)
    rho = np.outer(v, v.conj())
    U = np.diag(np.exp(-1j * 1.7 * float(STEP) * d['n']))
    expected = rho.copy()
    for _ in range(7):
        expected = U @ sum(K @ expected @ K.conj().T for K in d['K']) @ U.conj().T
    assert np.max(abs(expected - evolve(rho, 7, cap, 1.7))) < 4e-15
    # A differential source-energy perturbation changes the phase; it is not
    # spuriously reported as arbitrary-drive rigidity.
    assert np.max(abs(evolve(rho, 2000, cap, 1.7) - evolve(rho, 2000, cap, 1.9))) > .001


def test_exact_linewidth_bounds_and_phase_invariant_onset_control():
    intervals = []
    for cap in (4, 8, 16):
        c = perron_certificate(cap)
        lo, hi = map(F, c['linewidth_source_units_interval_exact'])
        assert 0 < lo < hi < F(1, 16)
        intervals.append((lo, hi))
        for k in range(2, cap + 1):
            b = perron_certificate(cap, k)
            assert F(b['spectral_radius_interval_exact'][1]) < 1
    assert intervals[1][1] < intervals[0][0]
    assert intervals[2][1] < intervals[1][0]
    d = collective(4)
    vacuum = np.zeros((5, 5)); vacuum[0, 0] = 1
    after = evolve(vacuum, 1, 4)
    assert after[1, 1].real == float(STEP * d['gamma'])
    assert np.trace(d['a'] @ after) == 0


def test_complete_records_include_measurement_and_battery_work():
    r = records(cap=4, burn_steps=1 << 16, flight_steps=1 << 14)
    assert len(r['complete_eight_probe_records']) == 8
    assert r['work_including_all_probe_and_battery_branches_diagnostic']['pumping_from_native_vacuum'] > 0
    assert r['unconditional_phase_mean_before_first_probe_exact'] == '0'


def test_single_particle_gap_is_not_array_total_energy():
    d = collective(1)
    b = perron_certificate(1)
    assert b['angular_frequency_over_native_mass_exact'] == '1'
    assert np.count_nonzero(d['K'][2]) == 0
    vacuum = np.diag([1., 0.])
    after = evolve(vacuum, 32768, 1)
    expected = 1 - (1 - float(STEP * d['gamma'])) ** 32768
    assert abs(after[1, 1].real - expected) < 2e-12
    assert after[0, 1] == 0


def test_zero_pump_retains_the_phase_prepared_by_the_observer():
    k = 65536
    driven = single_branch_records(k, PUMP)
    free = single_branch_records(k, F(0))
    assert driven['prior_pumped_occupation_cancels_exactly'] is True
    a = driven['complete_eight_records'][0]['probability_diagnostic'] - .125
    b = free['complete_eight_records'][0]['probability_diagnostic'] - .125
    assert abs(a / b - (1 - float(STEP * PUMP)) ** (k / 2)) < 2e-14
    assert abs(b) > .05


def test_exact_witness_is_checked_instead_of_trusting_a_solver_frame():
    import pytest
    row = perron_certificate(4)
    bad = list(row['positive_dyadic_witness'])
    bad[0] = str(F(1, 1 << 500))
    with pytest.raises(AssertionError):
        perron_certificate(4, witness=bad)


def test_single_branch_counterfactual_matches_actual_complete_instruments():
    steps = 65536
    sx = np.array([[0., 1.], [1., 0.]])
    sy = np.array([[0., -1j], [1j, 0.]])
    H = np.diag([0., 1.])
    target = single_branch_records(steps)['complete_eight_records']
    for population in (0., .3, 1.):
        rho = np.diag([1 - population, population])
        actual = []; first_energy = between_energy = final_energy = 0.
        for sign in (-1, 1):
            K = (np.eye(2) + sign * sx) / 2
            state = K @ rho @ K
            first_energy += np.trace(state @ H).real
            state = evolve(state, steps, 1)
            between_energy += np.trace(state @ H).real
            for O in (sx, sy):
                for final in (-1, 1):
                    L = (np.eye(2) + final * O) / 2
                    post = .5 * L @ state @ L
                    actual.append(float(np.trace(post).real))
                    final_energy += np.trace(post @ H).real
        assert np.max(abs(np.array(actual) - [x['probability_diagnostic'] for x in target])) < 3e-12
        assert abs(first_energy - .5) < 3e-15
        assert abs(final_energy - .5) < 3e-12
        assert between_energy > first_energy
        # Pumping work during flight is extracted by the complete final probe.
        assert abs((between_energy - first_energy) + (final_energy - between_energy)) < 3e-12
