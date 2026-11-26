"""Independent scale-clock receipts and two physical controls outside the survey menu."""
from copy import deepcopy
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from math import pi
from types import SimpleNamespace
import math

import pytest

from . import body, experiment, survey
from .clock_pairs import regions


def test_ordinary_release_gain_and_both_interior_surveys_follow_the_source(monkeypatch):
    mass = F(1, 3)
    frame = (F(1), F(0), F(0), F(0))
    for coupling in (F(1, 4), F(1, 2)):
        monkeypatch.setattr(body, "G", coupling)
        mounted = body.Body(mass, 0)
        duration, radius = 1/128, 2
        measured = 2*(radius-mounted.radius_after_release(radius, duration))*radius**2/(float(mass)*duration**2)
        record = body.release_calibration(mass, 128)
        assert record["gain"] == pytest.approx(measured, abs=1e-10)
        assert record["gain"] == pytest.approx(float(coupling), abs=2e-7)
        expected_h = math.sqrt(2*float(coupling*mass))
        assert float(survey.read_ball(mass, F(1), 128).measured_h) == pytest.approx(expected_h)
        assert regions.geometry(mass, frame, 128, 0)["h"] == pytest.approx(expected_h)


def test_release_meter_uses_the_quantized_physical_radius(monkeypatch):
    # One changed physical record must change the estimate even at sub-binary64 resolution.
    def altered_release(self, start, duration):
        return start - 11*duration**10

    monkeypatch.setattr(body.Body, "radius_after_release", altered_release)
    record = body.release_calibration(F(1, 3), 8192)
    assert record["start"] - record["released"] == 11
    expected_gain = 2*11*F(1, 8192)**8*4*3
    assert record["gain"] == float(expected_gain)


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


def test_unenrolled_moving_clock_misses_every_pair_at_an_enrolled_radius():
    with localcontext() as ctx:
        ctx.prec = 85
        frame = (F(25, 24), F(0), F(0), F(7, 24))
        assert frame[0]**2-sum(x*x for x in frame[1:]) == 1
        pairs = experiment.registered_pairs()
        held_h = (D(19)/66).sqrt()
        want = -38*pi/165*float(frame[0])**2
        for denominator in (512, 2048):
            ell = D(1)/denominator
            held = survey.read_ball(F(19, 33), frame[0], denominator)
            blank = survey.read_ball(F(0), frame[0], denominator)
            stamps, _ = survey.frame_clocks(frame, held_h, ell)
            velocity = experiment.measured_frame(stamps, denominator)
            boundary = experiment._crossing(velocity, held.radius, held.measured_h)
            assert all(right < boundary or boundary < left for left, right in pairs.values())
            # The different moving observer still measures positive dust and the predicted area.
            assert abs(velocity[0]-frame[0]) < F(1, denominator**3)
            actual = 4*pi*float((held.area_over_solid_angle-blank.area_over_solid_angle)/ell**4)
            assert abs(actual-want) < 4e-6
