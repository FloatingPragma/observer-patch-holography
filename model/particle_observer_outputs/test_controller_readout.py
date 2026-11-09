"""Controller observability with original source normalization and records."""
from copy import deepcopy

import pytest
import sympy as s

from model.particle_observer_outputs import controller_readout as c


@pytest.fixture(scope='module')
def controls():
    z = (63+16*s.I)/65
    return {name: c.source_example(r, phase, name) for name, r, phase in (
        ('zero_phase', '1', 1), ('nonzero_phase', '1/2', z), ('zero_pass', '0', z))}


@pytest.fixture(scope='module')
def joined():
    return c.spectral_inverse_join()


def test_literal_moments_and_complete_spectrum(controls):
    assert c.moments.moment_certificate()['literal_three_cycle_replayed']
    for data in controls.values():
        result = data['controller_readout']
        assert c._clean(s.trace(data['literal_Y'])-result['trace_Y']) == 0
        assert c._clean(s.trace(data['literal_Y']**2)-result['trace_Y_squared']) == 0
        assert len(data['all16_original_tracial_flag_probabilities']) == 16
        assert c._clean(sum(data['all16_original_tracial_flag_probabilities'])-1) == 0
        assert result['positive_root_branch_verified']
        assert not result['physical_particle_mass_claimed']
        assert not result['target_calibration_ancestry_erased']


def test_zero_and_exact_admitted_nonzero_phase_controls(controls):
    a = controls['zero_phase']['controller_readout']
    b = controls['nonzero_phase']['controller_readout']
    assert a['r'] == 1 and a['cos_three_theta'] == 1 and a['principal_abs_theta'] == 0
    assert b['r'] == s.Rational(1, 2)
    assert b['cos_three_theta'] == s.Rational(201663, 274625)
    assert b['principal_abs_theta'] == s.acos(s.Rational(201663, 274625))/3
    assert not b['controller_uniquely_selected_by_source_dynamics']
    # The inherited source parser has a narrower exact input interface.
    with pytest.raises(ValueError, match='rational unit-circle'):
        c.attenuation.family('1/2', s.exp(s.I*s.pi/12).expand(complex=True))


def test_zero_pass_has_no_spectral_phase_information(controls):
    result = controls['zero_pass']['controller_readout']
    assert result['sigma_y'] == (s.Rational(5, 108),)*3
    assert result['r'] == 0
    assert result['cos_three_theta'] is None
    assert result['principal_abs_theta'] is None
    assert not result['phase_observable']
    with pytest.raises(ValueError, match='phase-independent second moment'):
        c.inverse_moments(s.Rational(5, 36), s.Rational(1, 100))


def test_unordered_spectrum_does_not_fix_phase_orientation(controls):
    z = (63+16*s.I)/65
    _, reversed_values = c.moments.roots_and_masses(s.re(z), -s.im(z), '1/2')
    original = controls['nonzero_phase']['controller_readout']
    reversed_result = c.controller_from_coefficients(reversed_values)
    assert reversed_result['cos_three_theta'] == original['cos_three_theta']
    assert reversed_result['principal_abs_theta'] == original['principal_abs_theta']
    # A family permutation leaves the estimator unchanged.
    permuted = c.controller_from_coefficients(list(original['sigma_y'])[1:]+list(original['sigma_y'])[:1])
    assert permuted['principal_abs_theta'] == original['principal_abs_theta']
    assert permuted['unordered_family_permutations_unresolved']
    assert '2*pi*k/3' in permuted['phase_equivalence']


def test_declared_scale_and_ratio_normalization_give_same_controller(controls):
    values = controls['nonzero_phase']['full_spectrum']
    scaled = [7*v for v in values]
    declared = c.controller_from_coefficients(scaled, coefficient_scale=7)
    ratio = c.controller_from_ratios(scaled)
    for key in ('r', 'cos_three_theta', 'principal_abs_theta'):
        assert declared[key] == ratio[key]
    assert ratio['common_coefficient_scale_reconstructed'] == 7
    assert ratio['source_spectral_Koide_ratio'] == s.Rational(5, 12)
    assert not ratio['physical_energy_calibration_required']
    assert ratio['positive_root_trace_assumed_for_scale_reconstruction']
    with pytest.raises(ValueError, match='attenuation'):
        c.controller_from_coefficients(scaled)


def test_input_admission_rejects_wrong_moments_and_signed_root_branch():
    with pytest.raises(ValueError, match='attenuation'):
        c.inverse_moments(s.Rational(1, 100), s.Rational(1, 100))
    with pytest.raises(ValueError, match='cos'):
        c.inverse_moments(s.Rational(5, 18), 1)
    _, negative_c_root = c.moments.roots_and_masses(-1, 0, '1')
    with pytest.raises(ValueError, match='positive-C-root'):
        c.controller_from_coefficients(negative_c_root)
    for values in ((1, 2), (1, 2, 0), (1, 2, -1)):
        with pytest.raises(ValueError, match='three positive'):
            c.controller_from_coefficients(values)
    with pytest.raises(ValueError):
        c.controller_from_coefficients((0.1, 0.2, 0.3))
    with pytest.raises(ValueError, match='attenuation'):
        c.controller_from_ratios((1, 1, 1000))


def test_complete_charged_inverse_retains_original_weights_and_degeneracy(joined):
    rows = joined['complete_original_spectral_rows']
    probability = sum(sum(line['unconditional_probe_probability'] for line in row['addition_lines'])
                      + row['removal_unconditional_probe_probability'] for row in rows)
    assert c._clean(probability-1) == 0
    assert sorted(row['multiplicity'] for row in rows if row['internal_channel'] == 'left_doublet') == [4, 8]
    for result in joined['reconstructed_channels'].values():
        assert result['r'] == result['cos_three_theta'] == 1
        assert result['full_family_multiplicity'] == 3
        assert not result['conditioning_or_renormalization_of_original_records']
    assert joined['same_full_CAR_Gram_theorem_consumed']
    assert not joined['large_carrier_simulation_reexecuted']


def test_incomplete_and_invisible_reconstructed_channels_are_rejected(joined):
    rows = joined['complete_original_spectral_rows']
    y = c.attenuation.family('1', 1)['Y']
    gap = c._clean(s.sqrt(1+16*s.trace(y.H*y)))
    reconstructed = [c.bare_effective.channel_reconstruction(gap, row) for row in rows]
    bad = deepcopy(reconstructed)
    next(row for row in bad if row['internal_channel'] == 'left_doublet')['multiplicity'] = 3
    with pytest.raises(ValueError, match='multiplicity'):
        c.controller_from_reconstructed(bad)
    bad = deepcopy(reconstructed)
    next(row for row in bad if row['internal_channel'] == 'left_doublet')['original_addition_probabilities'] = (1, 0)
    with pytest.raises(ValueError, match='strictly positive'):
        c.controller_from_reconstructed(bad)
    with pytest.raises(ValueError, match='three positive'):
        c.controller_from_reconstructed([])
