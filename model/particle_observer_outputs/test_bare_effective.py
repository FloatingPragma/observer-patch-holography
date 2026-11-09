"""Same-source spectral inversion, visibility and exact uncertainty controls."""
from copy import deepcopy
from fractions import Fraction as F

import pytest
import sympy as s

from model.native_flavor_orbit import intervals as iv
from model.particle_observer_outputs import bare_effective as b
from model.particle_observer_outputs import charged_higgs_records as charged


@pytest.fixture(scope='module')
def fixture():
    data = charged.charged_data(s.diag(1, 2, 3), 1, 0)
    return data, charged.response_lines(data)


def test_actual_complete_source_lines_recover_both_family_coefficient_patterns(fixture):
    data, rows = fixture
    recovered = [b.channel_reconstruction(data['neutral']['gap'], row) for row in rows]
    assert b.coefficient_ratios(recovered, 'left_doublet') == (1, 2, 3)
    assert b.coefficient_ratios(recovered, 'charged_singlet') == (1, 2, 3)
    for result, row in zip(recovered, rows, strict=True):
        assert result['reconstructed_sigma_squared'] == row['sigma_squared']
        assert all(p > 0 for p in result['original_addition_probabilities'])
        assert not result['charged_additive_offset_used']
        assert not result['physical_rest_mass_claimed']


def test_actual_charge_offset_cancels_even_with_negative_relative_energies(fixture):
    original, _ = fixture
    shifted = charged.charged_data(s.diag(1, 2, 3), 1, -100)
    a, z = b.reconstruct_source(original), b.reconstruct_source(shifted)
    assert b.coefficient_ratios(a) == b.coefficient_ratios(z)
    assert all(x['charged_line_separation'] == y['charged_line_separation']
               for x, y in zip(a, z, strict=True))
    assert all(line['ground_referenced_energy'] < 0
               for row in charged.response_lines(shifted) for line in row['addition_lines'])


def test_common_clock_scale_changes_coefficients_and_cancels_in_ratios(fixture):
    scale = s.Rational(3, 7)
    scaled = charged.charged_data(scale*s.diag(1, 2, 3), scale, 0)
    a, z = b.reconstruct_source(fixture[0]), b.reconstruct_source(scaled)
    assert b.coefficient_ratios(z) == (1, 2, 3)
    for old, new in zip(a, z, strict=True):
        assert s.simplify(new['reconstructed_sigma']-scale*old['reconstructed_sigma']) == 0
        assert new['original_addition_probabilities'] == old['original_addition_probabilities']


def test_actual_coherent_z_one_keeps_its_repeated_pair():
    recovered = b.reconstruct_source(charged.coherent_charged_data())
    ratios = b.coefficient_ratios(recovered)
    assert ratios[:2] == (1, 1)
    assert s.simplify(ratios[2]-(34+24*s.sqrt(2))) == 0
    assert sorted(row['multiplicity'] for row in recovered if row['internal_channel'] == 'left_doublet') == [4, 8]


def test_zero_residue_is_rejected_rather_than_called_an_observed_line():
    data = charged.charged_data(s.diag(0, 1, 2), 1, 0)
    row = next(row for row in charged.response_lines(data) if row['sigma_squared'] == 0)
    assert row['addition_lines'][1]['unconditional_probe_probability'] == 0
    with pytest.raises(ValueError, match='strictly positive'):
        b.channel_reconstruction(data['neutral']['gap'], row)


def test_incompatible_spectral_data_and_incomplete_channels_fail(fixture):
    _, rows = fixture
    with pytest.raises(ValueError, match='negative squared'):
        b.channel_reconstruction(1, rows[0])
    incomplete = deepcopy(rows[0]); incomplete['addition_lines'] = incomplete['addition_lines'][:1]
    with pytest.raises(ValueError, match='both complete'):
        b.channel_reconstruction(15, incomplete)
    invisible = deepcopy(rows[0]); invisible['addition_lines'][0]['unconditional_probe_probability'] = 0
    with pytest.raises(ValueError, match='strictly positive'):
        b.channel_reconstruction(15, invisible)


def test_rational_spectral_intervals_enclose_coefficients_and_ratios(fixture):
    data, rows = fixture
    boxes = []
    error = iv.I(-F(1, 10**8), F(1, 10**8))
    for target in (1, 4):
        row = next(r for r in rows if r['internal_channel'] == 'left_doublet' and r['sigma_squared'] == target)
        low, high = row['addition_lines']
        lower, upper = (iv.expression(line['ground_referenced_energy'])+error for line in (low, high))
        weights = tuple(iv.expression(line['unconditional_probe_probability']).lo for line in (low, high))
        result = b.interval_reconstruction(iv.I(15)+error, lower, upper, kappa=1,
                                           original_weight_lower_bounds=weights)
        assert result['squared_coefficient'].lo <= target <= result['squared_coefficient'].hi
        assert result['coefficient'].lo <= s.sqrt(target) <= result['coefficient'].hi
        shifted = b.interval_reconstruction(iv.I(15)+error, lower+F(131, 7), upper+F(131, 7),
                                            kappa=1, original_weight_lower_bounds=weights)
        assert shifted['squared_coefficient'] == result['squared_coefficient']
        boxes.append(result['coefficient'])
    ratio = b.ratio_interval(boxes[1], boxes[0])
    assert ratio.lo <= 2 <= ratio.hi


def test_uncertainty_admission_does_not_invent_visibility_or_positivity():
    with pytest.raises(ValueError, match='visible'):
        b.interval_reconstruction((15, 15), (0, 1), (14, 15), kappa=1,
                                  original_weight_lower_bounds=(0, F(1, 10)))
    with pytest.raises(ValueError, match='nonnegative'):
        b.interval_reconstruction((1, 2), (0, 1), (14, 15), kappa=1,
                                  original_weight_lower_bounds=(F(1, 10), F(1, 10)))
    with pytest.raises(ValueError, match='denominator'):
        b.ratio_interval((1, 2), (0, 1))


def test_ambient_variational_bound_does_not_claim_the_physical_vacuum():
    row = b.ambient_ground_control()
    assert row['trial_coupling_squared'] == 63
    assert row['ambient_ground_energy_upper_bound'] < -7
    assert not row['trial_span_is_invariant']
    assert not row['neutral_sector_ground_is_full_ambient_ground']
    assert not row['trial_gauge_singlet_claimed']
    assert not row['physical_Gauss_vacuum_minimum_claimed']
