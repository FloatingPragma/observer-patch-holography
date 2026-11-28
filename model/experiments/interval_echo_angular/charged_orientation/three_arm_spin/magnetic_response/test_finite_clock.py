"""Replay the all-time current certificate on the original preparations."""

from .finite_clock import run


def test_three_threshold_guides_have_strict_native_record_ceilings():
    cases = run()
    assert [case["cyclic_ordinal_cuts"] for case in cases] == [19, 20, 21]
    for case, blank, written in zip(
        cases,
        (2.152e-6, 2.339e-6, 2.552e-6),
        (2.017e-7, 2.193e-7, 2.392e-7),
        strict=True,
    ):
        assert case["maximum_marked_crossings_per_arm"] == 1
        assert 0 < case["conditional_blank_meter_score_upper"] <= blank
        assert 0 < case["conditional_marked_meter_score_upper"] <= written
