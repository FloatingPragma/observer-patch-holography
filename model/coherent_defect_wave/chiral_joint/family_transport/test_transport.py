"""Actual-current families, physical source pulses and complete contrasts."""
import pytest
import sympy as s

from .model import (preparation, preparation_pulse_column, empty_key,
    validate_key, local_arrival_projector, native_charges)
from .checks import dictionary_checks, preparation_and_detector_checks, full_powers


def test_actual_nonabelian_dictionary_and119_hardware():
    row = dictionary_checks()
    assert row['native_families'] == 1 and row['native_internal_dimensions_sum'] == 15
    assert [r['dimension'] for r in row['native_irrep_rows']] == [3, 6, 1, 2, 3]
    assert all(r['actual119_hardware_complete_transport_instrument'] for r in row['native_irrep_rows'])
    assert row['common_source_native_encoding']['all_sector_native_current_intertwiner_proved']


def test_full_reference_preparations_and_complete_local_detectors():
    row = preparation_and_detector_checks()
    assert not row['static_reference_discarded']
    assert all(r['joint_Gauss_checks_on_preparation_and_first_image'] == 192 for r in row['rows'])
    assert all(r['complete_local_and_exclusive_arrival_derivative_norm_squared'] == '1/8' for r in row['rows'])
    assert all(s.Rational(r['wrong_reference_orientation_residual_norm_squared']) > 0 for r in row['rows'])


def test_all30_full_source_second_powers_and_same_clock_contrast():
    row = full_powers()
    assert len(row['source_mode_powers']) == 30
    assert row['common_source_time'] == '1/1508'
    assert row['common_outward_second_image_norm_bound'] == 377
    low = s.Rational(row['prepared_local_arrival_probability_interval'][0])
    background = s.Rational(row['preparation_off_local_arrival_probability_interval'][1])
    assert low - background == s.Rational(row['strict_local_record_contrast_lower_bound']) > 0
    for name in ('prepared_exclusive', 'prepared_local', 'preparation_off_local'):
        a, b = map(s.Rational, row[name + '_arrival_probability_interval'])
        c, d = map(s.Rational, row[name + '_complement_probability_interval'])
        assert 0 <= a <= b <= 1 and a + d == b + c == 1


def test_reference_and_preparation_domain_rejections():
    with pytest.raises(ValueError):
        preparation(0)
    with pytest.raises(ValueError):
        preparation_pulse_column(empty_key(), -6)
    with pytest.raises(ValueError):
        validate_key((empty_key()[0], 1 << 30))
    for charge in native_charges():
        assert not local_arrival_projector(preparation(charge), charge)


def test_custody_rejects_missing_declared_and_import_inputs(monkeypatch):
    from . import certificate as custody
    with monkeypatch.context() as patch:
        patch.setattr(custody.action, 'SOURCE_FILES', (*custody.action.SOURCE_FILES, 'missing-family-transport-seed.json'))
        with pytest.raises(FileNotFoundError, match='missing declared'):
            custody.input_hashes()
    monkeypatch.setattr(custody, 'import_closure', lambda seeds: set(seeds) | {'missing-family-transport-import.py'})
    with pytest.raises(FileNotFoundError, match='missing family-transport import'):
        custody.input_hashes()
