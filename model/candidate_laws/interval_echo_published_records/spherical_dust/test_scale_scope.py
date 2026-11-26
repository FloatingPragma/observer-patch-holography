"""Independent scale-clock receipts and two physical controls outside the survey menu."""
from copy import deepcopy
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from math import pi
from types import SimpleNamespace

import pytest

from . import experiment, survey


def test_dust_scale_from_three_clock_and_circumference_records():
    with localcontext() as ctx:
        ctx.prec = 90
        for denominator in (64, 257, 8192):
            for actual in (D(0), (D(19)/66).sqrt(), D(1)):
                tape = survey.DustEmitter(actual).calibration(denominator)
                recovered, bound = survey.decode_scale(tape)
                assert abs(recovered-actual) <= bound
                assert abs(bound*denominator**10-256) < D("1e-75")

        tape = survey.DustEmitter(D("0.5")).calibration(512)
        altered = deepcopy(tape)
        altered["epochs"][0]["inner_radius"] += 40
        with pytest.raises(ValueError, match="two tagged shells"):
            survey.decode_scale(altered)
        # Preserve the two-shell ratio while breaking the independent timing comparison.
        altered["epochs"][0]["outer_radius"] += 80
        with pytest.raises(ValueError, match="incompatible scale histories"):
            survey.decode_scale(altered)


def test_receiver_uses_its_scale_tape_with_the_source_optical_records():
    denominator = 64
    emitter = survey.DustEmitter(D("0.5"))
    # The receiver has exactly two recording interfaces and no access to h or the stock.
    lines = SimpleNamespace(endpoint=emitter.endpoint, circumferences=emitter.circumferences)
    good = survey.survey_ball(lines, F(1), denominator, emitter.calibration(denominator))
    exported = good.scale_record
    exported["epochs"].clear()
    assert good.scale_record == emitter.calibration(denominator)
    flat_tape = survey.DustEmitter(D(0)).calibration(denominator)
    swapped = survey.survey_ball(lines, F(1), denominator, flat_tape)
    assert abs(good.measured_h-D("0.5")) <= good.scale_error_bound
    assert swapped.measured_h == 0
    assert abs(swapped.area_over_solid_angle-good.area_over_solid_angle)*denominator**4 > 1


def test_interleaved_rest_balls_and_off_axis_centers_cut_no_enrolled_memory():
    with localcontext() as ctx:
        ctx.prec = 85
        actual_h = (D(19)/66).sqrt()
        held, blank = survey.DustEmitter(actual_h), survey.DustEmitter(D(0))
        fixed_pairs = experiment.registered_pairs()
        # The center (1/8, -1/16, 0) is orthogonal to the enrolled (1,2,4) line.
        center_offset_squared = D(5)/256
        expected_deficit = -38*pi/165
        errors = []
        for denominator in (512, 2048):
            target = D(3)/denominator
            ball = survey.survey_ball(held, F(1), denominator, held.calibration(denominator),
                                      target_radius=target)
            reference = survey.survey_ball(blank, F(1), denominator, blank.calibration(denominator),
                                           target_radius=target)
            crossing = survey.comoving_distance(ball.radius, ball.measured_h, D(0))
            assert all(right < crossing or crossing < left for left, right in fixed_pairs.values())
            assert any(right < crossing for _, right in fixed_pairs.values())
            if denominator > 512:
                assert any(crossing < left for left, _ in fixed_pairs.values())
            # Translating the rod/clock apparatus leaves the prepared memory locations fixed.
            assert all(center_offset_squared+q*q > crossing**2
                       for pair in fixed_pairs.values() for q in pair)
            # Whole-pair unitaries conserve entropy. The signed residual is the measured area.
            signed = 4*pi*float((ball.area_over_solid_angle-reference.area_over_solid_angle)/target**4)
            errors.append(abs(signed-expected_deficit))
            assert abs((ball.volume_over_solid_angle-target**3/3)/target**5) < D("1e-10")
        assert errors[1] < errors[0]/8 < 2e-6
