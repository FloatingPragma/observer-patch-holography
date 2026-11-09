from fractions import Fraction as F

import pytest

from model.particle_common_contract import contract
from . import window


def test_complete_finite_controls_distinguish_zero_and_ultraviolet_weight():
    rows = {row['name']: row for row in window.exact_phase_controls()}
    for name in ('positive_excitation', 'mixed_complete_response', 'small_negative_tail'):
        assert rows[name]['bound']['original_positive_band_weight_lower'] > 0
    for name in ('zero_atom', 'ultraviolet_escape'):
        assert rows[name]['bound']['original_positive_band_weight_lower'] == 0
    assert all(not row['is_selected_spatial_source'] for row in rows.values())


def test_original_measure_scaling_and_error_are_not_success_normalization():
    base = window.positive_window(1, 0, 1, 1, 1, 0, F(1, 4), 16)
    small = window.positive_window(F(1, 100), 0, 1, F(1, 100), F(1, 100), 0, F(1, 4), 16)
    assert small['original_positive_band_weight_lower'] == base['original_positive_band_weight_lower']/100
    uncertain = window.positive_window(1, 1, 1, 1, 1, 0, F(1, 4), 16)
    assert uncertain['original_positive_band_weight_lower'] == 0


def test_negative_transfer_cannot_be_mistaken_for_positive_energy():
    value = window.positive_window(2, 0, 4, 1, 1, 1, F(1, 4), 16)
    assert value['original_positive_band_weight_lower'] == 0


def test_selected_allowances_reuse_actual_cooling_without_inventing_a_contrast():
    member = contract.member(2)
    result = window.selected_requirement(member, 0, 0)
    assert result['member_id'] == member['member_id']
    assert result['negative_tail']['original_negative_spectral_weight_upper'] == '1/16'
    assert result['conditional_window_bound']['original_positive_band_weight_lower'] == '0'
    assert not result['primary_positive_energy_weight_derived']
    assert result['contrast_is_supplied_not_computed_for_primary']


def test_inconsistent_or_reversed_bounds_rejected():
    with pytest.raises(ValueError, match='contradicts'):
        window.positive_window(3, 0, 1, 1, 1, 0, 1, 2)
    with pytest.raises(ValueError, match='exceed'):
        window.positive_window(0, 0, 1, 1, 1, 0, 2, 1)
    with pytest.raises(ValueError):
        window.positive_window(1, -1, 1, 1, 1, 0, 1, 2)
    with pytest.raises(ValueError):
        window.positive_window(True, 0, 1, 1, 1, 0, 1, 2)
