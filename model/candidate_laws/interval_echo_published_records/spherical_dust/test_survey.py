"""Regression tests for the common dust source, its actual clocks and independent work."""
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
import json
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from . import endpoint_counts, experiment, survey
from .body import Body, release_calibration


def test_received_stamp_reconstruction_and_spacelike_shooting():
    with localcontext() as ctx:
        ctx.prec = 65
        h, s, gamma = D("0.55"), D(1)/32, D(3)
        b = (gamma*gamma-1).sqrt()
        for direction in [(D(1), D(0), D(0)), (D(0), D(0), D(-1))]:
            t, x = survey.ray(s,h,gamma,b,direction)
            w = float(b*direction[2])
            momentum = np.array(direction, dtype=float)
            momentum[2] *= float(gamma)

            def rhs(_, state):
                factor = 1-1.5*float(h)*state[0]
                return [state[1], float(h)*(1+w*w)*factor**(-7/3), *momentum*factor**(-4/3)]

            solution = solve_ivp(rhs,(0,float(s)),[0,w,0,0,0],method="DOP853",rtol=1e-13,atol=1e-15)
            np.testing.assert_allclose([float(t),*map(float,x)],solution.y[[0,2,3,4],-1],atol=2e-12,rtol=0)
            tick = s**10
            stamps = survey.timestamp_event(t,x,h,tick)
            clock, location = survey.locate(stamps,h,tick)
            assert abs(clock-t) <= tick and max(abs(a-b) for a,b in zip(x,location,strict=True)) < 10*tick
            with pytest.raises(ValueError):
                survey.locate((0, 2, 2, 2, 2),h,tick)


def test_each_shrinking_ball_and_held_variants():
    data = experiment.measure()
    experiment.verify(data)
    assert data["registered_pairs"] == 21
    for stock in ("ordinary", "held"):
        previous = {}
        for row in data["ball_surveys"][stock]:
            assert len(row["observers"]) == 10
            for key in ("bulk_error_over_r4","area_error_over_r4","calibrated_balance_over_r4","max_velocity_component_error"):
                largest = max(abs(o[key]) for o in row["observers"])
                limit = previous[key]/10 if key in previous else 0.001
                assert largest < limit
                previous[key] = largest
            assert min(o["membership_margin_over_r"] for o in row["observers"]) > 0.005
            assert max(o["regional_relative_entropy"] for o in row["observers"]) < 1e5/row["inverse_radius"]**8
    held = data["ball_surveys"]["held"][-1]
    assert max(abs(o["calibrated_balance_over_r4"]) for o in held["observers"]) < 1.6e-7
    for trial, floor in (("half_collar_force",1.36),("blind_collar_work",2.73)):
        assert max(abs(o["calibrated_balance_over_r4"]) for o in data["ball_surveys"][trial][-1]["observers"]) > floor
    assert abs(held["rod_einstein_tensor"][0]-19/22) < 2e-8


def test_ordinary_calibration_forecasts_held_clocks_and_release():
    errors = [abs(release_calibration(F(1,3),d)["gain"]-0.25) for d in experiment.SCALES]
    assert errors[2] < errors[1]/15 < errors[0]/225
    data = experiment.measure()
    observed = data["exterior"]["held"]
    assert abs(observed["received_frequency_per_emitted"]-math.sqrt(226/245)) < 2e-13
    assert abs(data["held_forecast"]["held_clock_ratio"]-observed["received_frequency_per_emitted"]) < 8e-8
    releases = observed["signed_radius_acceleration"]
    assert abs(releases[1]+19/528) < abs(releases[0]+19/528)/3
    assert data["exterior"]["half_collar_force"]["received_frequency_per_emitted"] > observed["received_frequency_per_emitted"]+0.008
    with pytest.raises(ValueError):
        Body(1/3,0).shell(0.4,20)


def test_counted_entropy_is_read_independently_of_geometry():
    candidates = experiment.measure()
    counts = endpoint_counts.replay(candidates)
    last = counts["levels"][-1]["comparisons"]
    assert last["held"]["stationarity_upper"] < 0.00016
    assert last["half_collar_force"]["stationarity_lower"] > 1.36
    assert last["blind_collar_work"]["stationarity_lower"] > 2.73
    assert abs(last["collar"]["entropy_stress_tensor"][0]-2/(11*math.pi)) < last["collar"]["component_sampling_radius"]
    # Two preparations with identical mass/area records still have different phase counts.
    saved = json.loads((experiment.HERE/"endpoint_counts.json").read_text())["levels"]
    for level in saved:
        for name in experiment.NAMES:
            level["banks"]["held"][name]["ones"] = level["banks"]["blank"][name]["ones"]
            level["banks"]["held"][name]["zeros"] = level["banks"]["blank"][name]["zeros"]
    altered = endpoint_counts.evaluate(candidates, saved)["levels"][-1]["comparisons"]["held"]
    # Invert each bank in its own measured frame. The finite frame differences
    # multiply only the common blank sampling noise after this substitution.
    assert max(map(abs, altered["entropy_stress_tensor"])) < 1e-12
    assert altered["stationarity_lower"] > 6.5


def test_endpoint_emitter_rejects_corrupted_independent_weighing():
    data = experiment.measure()
    frozen = json.loads((experiment.HERE/"endpoint_counts.json").read_text())["levels"]
    survey_row = data["ball_surveys"]["held"][0]
    scale, inner, outer = survey_row["mass_tallies"][1]
    survey_row["mass_tallies"][1] = (scale, inner, inner+(outer-inner)//2)
    with pytest.raises(ValueError, match="deposited finite emitter differs"):
        endpoint_counts.evaluate(data, frozen)


def test_measured_frame_inverse_rejects_a_repeated_observer():
    observers = experiment.measure()["ball_surveys"]["held"][0]["observers"]
    observers[-1]["frame_phase_stamps"] = observers[-2]["frame_phase_stamps"]
    with pytest.raises(ValueError, match="singular"):
        experiment.frame_inverse(observers, experiment.SCALES[0])
