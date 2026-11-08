from fractions import Fraction as F
import math

import numpy as np
import pytest
from scipy.linalg import expm

from model.spin_modular_geometry import regional
from . import algebra, response, native


def test_actual_source_weighted_identity_all_axes_and_charge():
    for axis in range(3):
        row = algebra.operators(axis)
        h, b, p = row['H8']/8, row['B16']/16, 1j*row['P_imag128']/128
        assert np.array_equal(p, 1j*(h@b-b@h))
        assert np.array_equal(h@p-p@h, -1j*row['force1024']/1024)
        assert np.any(h@p-p@h)
        assert not np.any(row['twice_charge']@p-p@row['twice_charge'])
        h0 = row['zero_field_H8']/8
        assert np.any(h0@p-p@h0)
        assert np.array_equal(1j*(h0@p-p@h0), row['zero_field_force1024']/1024)
        assert len(native.terms(axis)) == 32


def test_source_ground_defect_and_variance_independently():
    for row in algebra.source_controls():
        op = algebra.operators(row['axis'])
        h, b, p = op['H8']/8, op['B16']/16, 1j*op['P_imag128']/128
        values, vectors = np.linalg.eigh(h)
        ground, energy = vectors[:, 0], values[0]
        relative = h-energy*np.eye(256)
        second = np.vdot(b@ground, b@ground).real
        moment = np.vdot(b@ground, relative@b@ground).real
        d = 1j*(b@p-p@b)+relative
        defect = np.vdot(ground, d@ground).real
        assert abs(np.vdot(ground, b@ground)) < 1e-13
        assert float(row['ground_boost_second_moment_lower']) <= second <= float(row['ground_boost_second_moment_upper'])
        assert float(row['ground_boost_energy_moment_lower']) <= moment <= float(row['ground_boost_energy_moment_upper'])
        assert float(row['ground_boost_defect_lower']) <= defect <= float(row['ground_boost_defect_upper'])
        assert abs(defect+2*moment) < 1e-13
        assert np.trace(d).real > 0


def test_coarse_transport_requires_nonzero_dipole_even_when_energy_ward_holds():
    row = algebra.coarse_transport()
    assert row['axis_block_maps_checked'] == 24
    assert row['maps_with_nonzero_subcell_energy_dipole'] == 12
    op = algebra.operators(0)
    assert np.any(op['B16'])
    # One total-energy cell has representative0 and zero coarse polarization.
    assert not np.array_equal(op['B16'], np.zeros_like(op['B16']))
    assert np.any(op['P_imag128'])


def test_forced_stress_attaches_to_created_charge_channels():
    rows = algebra.forced_stress()
    assert len(rows) == 3
    for row in rows:
        assert row['local_balance_cells'] == 8
        assert row['nonzero_antisymmetric_stress_pairs'] == 28
        assert row['all_eight_coordinate_block_maps_preserve_force'] is True
        for channel in row['created_charge_channel_force']:
            threshold = F(44, 100) if channel['charge'] == -1 else F(13, 100)
            assert channel['unconditioned_energy_flux_acceleration_lower'] > threshold


def test_both_hyperbolic_conjugation_bounds_against_actual_source():
    op = algebra.operators(0)
    h, b, p = op['H8']/8, op['B16']/16, 1j*op['P_imag128']/128
    e0 = np.linalg.eigvalsh(h)[0]
    relative = h-e0*np.eye(256)
    for eta in (F(-1, 16), F(1, 32), F(1, 8)):
        u = expm(1j*float(eta)*b)
        h_error = u@relative@u.conj().T-relative*math.cosh(float(eta))+p*math.sinh(float(eta))
        p_error = u@p@u.conj().T-p*math.cosh(float(eta))+relative*math.sinh(float(eta))
        bounds = response.history_bound(eta, algebra.source_controls()[0]['boost_defect_norm_upper'], [F(1)], [F(1)])
        assert np.linalg.norm(h_error, 2) < float(bounds['clock_generator_error_upper'])
        assert np.linalg.norm(p_error, 2) < float(bounds['candidate_flux_generator_error_upper'])
        # The flux error is first order: its source-ground expectation is nonzero.
        _, vectors = np.linalg.eigh(h)
        assert abs(np.vdot(vectors[:, 0], p_error@vectors[:, 0]).real/float(eta)) > 0.12


def test_positive_speed_normalization_cannot_cancel_ground_defect():
    for speed in (F(1, 2), F(1), F(2), F(100)):
        for row in algebra.speed_normalization(speed):
            assert row['ground_normalized_boost_defect_upper'] < 0
            assert row['physical_speed_selected_by_this_normalization'] is False


def test_actual_source_translation_log_is_not_energy_flux():
    value = algebra.finite_translations()
    assert value['nonidentity_shifts_form_Z2_times_Z2'] is True
    for row in value['even_shift_controls']:
        assert row['same_H_B_charge_and_selected_ground_preserved'] is True
        assert row['principal_log_conserved_by_same_H_B'] is True
        assert len(row['energy_flux_axes_reversed_by_shift']) == 2
        assert row['energy_flux_is_this_conserved_translation_log'] is False


@pytest.mark.parametrize('bad', [True, 3, -1, 0.0])
def test_axis_domain(bad):
    with pytest.raises(ValueError): algebra.operators(bad)


@pytest.mark.parametrize('bad', [F(0), F(-1), 1.0, True])
def test_candidate_speed_domain(bad):
    with pytest.raises(ValueError): algebra.speed_normalization(bad)


def test_rational_hyperbolic_enclosures_and_domains():
    for x in (F(0), F(1, 32), F(-1, 4), F(1)):
        row = response.hyperbolic_bounds(x)
        assert float(row['cosh_minus_one_lower'])-2e-16 <= math.cosh(float(x))-1 <= float(row['cosh_minus_one_upper'])+2e-16
        assert float(row['absolute_sinh_lower'])-2e-16 <= abs(math.sinh(float(x))) <= float(row['absolute_sinh_upper'])+2e-16
    with pytest.raises(ValueError): response.history_bound(F(2), F(1), [], [])
    with pytest.raises(ValueError): response.history_bound(F(1), F(-1), [], [])
    with pytest.raises(ValueError): response.history_bound(F(1), F(1), [0.2], [])
