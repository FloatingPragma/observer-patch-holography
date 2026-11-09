"""Exact source attachments and independent detector/Lehmann checks."""
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import eigh, norm
import sympy as s

from . import records as r
from model.particle_spatial_ground_state import ground


def array(matrix):
    return np.asarray(matrix.tolist(), dtype=complex)


def test_all_original_source_blocks_and_every_occupation_are_retained():
    blocks = r.rational_blocks()
    assert sorted(i for indices, _, _ in blocks for i in indices) == list(range(208))
    assert {n: sum(len(indices) for indices, occupation, _ in blocks if occupation == n)
            for n in range(5)} == {0: 13, 1: 52, 2: 78, 3: 52, 4: 13}
    assert max(map(lambda row: len(row[0]), blocks)) == 12


def test_exact_LDL_test_rejects_zero_and_negative_pivots():
    assert r._shift_is_positive(((F(2), F(1)), (F(1), F(2))), F(0))
    assert not r._shift_is_positive(((F(2), F(1)), (F(1), F(2))), F(1))
    assert not r._shift_is_positive(((F(1), F(2)), (F(2), F(1))), F(0))
    assert not r._shift_is_positive(((F(0),),), F(0))


def test_exact_sector_intervals_enclose_independent_full_matrix_spectra():
    h = array(r.source.exact_control()['H'])
    minima = r.sector_minima()
    for occupation, (lower, upper) in minima.items():
        indices = [i for i in range(208) if (i//13).bit_count() == occupation]
        value = eigh(h[np.ix_(indices, indices)], eigvals_only=True)[0]
        assert float(lower) <= value <= float(upper)
        assert 0 < upper-lower <= F(1, 1 << 32)
    assert minima[1] == minima[3]
    assert r.global_ground_interval() == minima[2]
    assert all(r.charged_thresholds()[n][0] > 1 for n in (1, 3))


def test_compensated_instrument_is_complete_and_retains_Pauli_zero_outcomes():
    rows = r.instrument()
    assert len(rows) == 8
    total = sum((row['K'].H*row['K'] for row in rows), s.zeros(208))
    assert total == s.eye(208)
    number = s.kronecker_product(r.source.field_route.rotor_coefficients()['number'], s.eye(13))
    plus = sum((row['K'].H*row['K'] for row in rows if row['sign'] == 1), s.zeros(208))
    minus = sum((row['K'].H*row['K'] for row in rows if row['sign'] == -1), s.zeros(208))
    assert plus == s.eye(208)-number/4
    assert minus == number/4
    empty_flux_zero, filled_flux_zero = 6, 15*13+6
    assert plus[empty_flux_zero, empty_flux_zero] == 1
    assert minus[empty_flux_zero, empty_flux_zero] == 0
    assert plus[filled_flux_zero, filled_flux_zero] == 0
    assert minus[filled_flux_zero, filled_flux_zero] == 1
    assert all(row['source_charge_change']+row['reference_charge'] == 0 for row in rows)
    assert all(row['reference_parity'] == 1 for row in rows)


def test_charge_conjugation_reverses_source_charge_and_pairs_probe_modes():
    c = r.charge_conjugation()['matrix']
    assert all(r.charge_conjugation()['identities'].values())
    rows = {(row['sign'], row['mode']): row['K'] for row in r.instrument()}
    for mode in range(4):
        transformed = c*rows[1, mode].conjugate()*c.H
        target = rows[-1, 3-mode]
        # The phase depends on the native gamma convention; the complete
        # rank-one instrument is unchanged by this branch phase.
        assert transformed in (s.I*target, -s.I*target)


def test_ground_tail_bound_includes_full_instrument_without_channel_factor():
    energies = np.array([0., .3, 2.])
    probability = np.array([.8, .15, .05])
    # A complete permutation instrument with an explicit downward branch.
    matrix = np.array([[0., 1., 0.], [0., 0., 1.], [1., 0., 0.]])
    weights = matrix*probability[None, :]
    epsilon = F(145, 1000)
    for eta in (F(1, 4), F(1), F(7, 4)):
        negative = float(weights[energies[:, None]-energies[None, :] <= -float(eta)].sum())
        assert negative <= float(r.negative_transfer_bound(epsilon, eta))+1e-15
    assert r.negative_transfer_bound(0, F(1, 4)) == 0
    assert r.negative_transfer_bound(F(1, 100), F(1, 4)) == F(1, 25)
    # A success-conditioned 5%-weight downward event would violate the
    # original-weight bound; the implementation retains its weight0.05.
    assert .05 <= float(r.negative_transfer_bound(epsilon, 1)) < 1


@pytest.mark.parametrize('args', [(-1, 1), (1, 0), (1, -1), (.1, 1), (True, 1)])
def test_negative_transfer_bound_rejects_wrong_domains(args):
    with pytest.raises(ValueError):
        r.negative_transfer_bound(*args)


def test_exact_clock_records_match_independent_full_matrix_computation():
    # A low-order genuine source polynomial exposes detector errors while
    # keeping this independent dense comparison inexpensive.
    prepared = r.source.prepared_density(order=1)
    law = r.complete_filter_records(prepared)
    h = array(r.source.exact_control()['H'])
    phase = np.asarray([complex(p) for p in prepared['computational_basis_phases']])
    q = np.zeros((208, 208), dtype=complex)
    for indices, part in prepared['amplitude_blocks']:
        q[np.ix_(indices, indices)] = array(part)
    q = phase[:, None]*q*phase.conjugate()[None, :]
    z = norm(q, 'fro')**2
    eye = np.eye(208)
    step = (eye-1j*h/512)@np.linalg.inv(eye+1j*h/512)
    u = np.linalg.matrix_power(step, 4)
    probe = {(row['sign'], row['mode']): array(row['K']) for row in r.instrument()}
    assert sum(row['original_probability'] for row in law['records']) == 1
    for row in law['records']:
        x = probe[row['sign'], row['mode']]@q
        y = u@x@u.conj().T
        expected = norm(x+(-1)**row['phase_bit']*y, 'fro')**2/(4*z)
        assert abs(expected-float(row['original_probability'])) < 3e-14
    assert law['complete_original_branch_count'] == 16
    assert law['prepared_exact_clock_output_half_trace_error_upper'] < F(1, 2000)


def test_missing_source_block_and_nonstationary_amplitude_are_rejected():
    prepared = dict(r.source.prepared_density(order=1))
    prepared['amplitude_blocks'] = prepared['amplitude_blocks'][1:]
    with pytest.raises(ValueError, match='all exact source blocks'):
        r.complete_filter_records(prepared)
    prepared = dict(r.source.prepared_density(order=1))
    blocks = list(prepared['amplitude_blocks'])
    for j, (indices, q) in enumerate(blocks):
        if len(indices) > 1:
            bad = s.MutableDenseMatrix(q); bad[0, 0] += 1
            blocks[j] = (indices, s.ImmutableMatrix(bad)); break
    prepared['amplitude_blocks'] = tuple(blocks)
    with pytest.raises(ValueError, match='commute with the same H'):
        r.complete_filter_records(prepared)


def test_cooled_complete_detector_certificate_has_original_half_branch_weights():
    result = r.cooled_control_certificate()
    law = result['complete_clock_apparatus']
    assert law['probability_sum'] == 1
    assert result['mean_ground_excess_interval'][1] < F(157, 1000)
    assert result['complete_negative_transfer_probability_upper'] < F(63, 100)
    for sign in (1, -1):
        intervals = [row['original_probability_interval'] for row in law['records'] if row['sign'] == sign]
        assert sum(lo for lo, _ in intervals) <= F(1, 2) <= sum(hi for _, hi in intervals)
    assert law['Gibbs_exact_clock_output_half_trace_error_upper'] < F(1, 2000)
    assert not law['negative_bin_is_negative_energy_measurement']
    assert result['resolved_readout_plan_from_executed_control']['original_negative_record_probability_upper'] < F(89, 100)


def test_source_Fourier_resolution_plan_keeps_alias_and_clock_errors():
    plan = r.resolved_negative_record_bound(F(1, 100))
    assert plan['alias_probability_upper'] == 0
    assert 0 < plan['filter_leakage_upper'] < F(3, 100)
    assert 0 < plan['complete_Cayley_clock_error_upper'] < F(1, 1000)
    assert plan['original_negative_record_probability_upper'] < F(11, 100)
    assert not plan['executed_full_distribution']


def test_diagnostic_Lemann_weights_cool_without_deleting_negative_transfers():
    result = r.spectral_diagnostics()
    rows = result['cooling']
    assert all(abs(row['all_transition_weight_sum_decimal']-1) < 1e-12 for row in rows)
    assert all(row['negative_transfer_at_least_quarter_weight_decimal'] > 0 for row in rows)
    assert all(row['negative_transfer_at_least_quarter_weight_decimal'] <=
               row['mean_excess_over_threshold_decimal']+1e-12 for row in rows)
    assert rows[-1]['negative_transfer_at_least_quarter_weight_decimal'] < 1e-6
    first_pair = result['ground_lines_above_display_threshold'][:2]
    assert {row['sign'] for row in first_pair} == {-1, 1}
    assert all(.49 < row['original_ground_probe_weight_decimal'] < .5 for row in first_pair)
    assert result['physical_particle_masses'] == []


def test_outward_dyadic_intervals_and_certificate_scope():
    for value in (F(-13, 7), F(1, 3), F(5, 8)):
        lower, upper = r.dyadic_interval(value)
        assert lower <= value <= upper
        assert upper-lower <= F(1, 1 << 40)
    certificate = r.certificate(include_control=False)
    assert not certificate['observed_particle_identification']
    assert not certificate['relativistic_mass_shell_derived']
    assert 'spectral_diagnostics' not in certificate
