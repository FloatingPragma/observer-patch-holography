"""Exact reconstruction of the actual fixed-repair coefficient readout."""
import sympy as s

from . import readout


def test_fixed_native_readout_reconstructs_its_excluded_complete_action():
    value = readout.certificate()
    assert value['raw_native_records_recomputed'] == 128
    assert value['native_advance_probability'] == '1/12'
    assert value['CKM_in_common_mass_basis'] == 'I3'
    assert value['Majorana_coefficient'] == 'zero3x3'
    assert value['each_window_disjoint_exact'] == [True, True]
    assert value['expensive_detector_reexecuted_by_this_package'] is False
    for interval, window in zip(value['charged_mass_ratios_after_lightest_exact'], value['retained_comparison_windows'], strict=True):
        assert s.Rational(interval[1]) < s.Rational(window[0])
