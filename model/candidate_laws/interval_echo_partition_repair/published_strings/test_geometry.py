import json
import math
from pathlib import Path

import numpy as np
import pytest

from model.candidate_laws.interval_echo_partition_repair.published_strings import geometry as g


def test_stress_from_action_at_fixed_material_labels_and_product_curvature():
    metric = np.diag([-1., 1., 1., 1.])
    stress = np.array(g.energy_tensor(0.25, 0.375))
    for i in range(4):
        for j in range(i, 4):
            change = np.zeros((4, 4))
            change[i, j] = change[j, i] = 1e-5
            forward = g.sheet_action(0.625, metric + change)
            backward = g.sheet_action(0.625, metric - change)
            assert abs(forward-backward - np.sum(stress*change)) < 2e-15
    K = g.curvature(0.25, 0.375)
    G = np.asarray(g.einstein_tensor(K))
    np.testing.assert_allclose(G, g.KAPPA*stress, rtol=1e-15, atol=0)
    charges = [g.project(G, u) for u in g.FRAME]
    for u, charge in zip(g.FRAME, charges, strict=True):
        assert math.isclose(charge, K*float(g.spatial_stretch(u)))
    from matched_einstein_limit import frame
    np.testing.assert_allclose([float(t) for t in frame.decode(charges)],
                               [G[i, j] for i, j in frame.PAIRS], atol=3e-14)


def test_ordinary_holonomy_and_action_predict_every_held_ball_and_controls():
    data = g.survey_report()
    ordinary = data["circle_surveys"]["ordinary"]
    mass, collar = data["densities"]["held"]
    read_gain = ordinary["holonomy"]/ordinary["disk_area"]/mass
    assert abs(read_gain-2*math.pi) < 3e-15
    for row in data["balls"]:
        ell, b = row["ell"], row["b2"]
        expected = -(4*math.pi/15)*read_gain*b*collar
        assert abs(row["stress_forecast"]-expected) < 2e-13
        assert abs(row["held_slope"]-expected) < ell**2*4*math.pi*5*(g.KAPPA*b)**2
        assert abs(row["held_slope"]/expected-1) < 0.016
        assert 0.49 < row["half_strength_slope"]/row["held_slope"] < 0.51
        assert row["zero_strength_slope"] == 0


def test_volume_stops_and_weighted_volume_agree_with_independent_quadrature():
    for ell in g.RADII:
        for u in g.FRAME:
            b = float(g.spatial_stretch(u))
            for K in (0, g.KAPPA/2, 2*g.KAPPA):
                ball = g.measured_ball(K, ell, u)
                target = 4*math.pi*ell**3/3
                assert abs(ball.volume-target) < 3e-15*target
                assert abs(ball.area-g.measured_ball(K, ell, u, order=24).area) < 3e-15*ball.area
                W = g.weighted_volume(K, ball.radius, u)
                # Direct integral of the weight against the derivative of V, read
                # here as independently quadrature-evaluated concentric sphere area.
                nodes, weights = np.polynomial.legendre.leggauss(8)
                r = ball.radius
                direct = sum(w * (r*r-(r*(t+1)/2)**2)/4 *
                             g.product_measures(K, r*(t+1)/2, b, order=14)[1]
                             for t, w in zip(nodes, weights, strict=True))
                assert abs(W-direct) < 3e-14*W
                C = 4*math.pi*ell**4/15
                assert abs(W-C) <= 4*math.pi/20*(K*b)*ell**6


def test_deposited_records_and_units():
    stored = json.loads((Path(__file__).parent/"geometry.json").read_text())
    assert len(stored["balls"]) == 40
    assert stored["records_are_model_integrals"]
    assert stored["flat_longitudinal_axis"] == "z"
    for row in stored["balls"]:
        read = g.measured_ball(g.curvature(*stored["densities"]["held"]), row["ell"],
                               g.FRAME[row["frame"]])
        assert read.area == pytest.approx(row["held_area"], rel=5e-14, abs=0)
    with pytest.raises(ValueError):
        g.curvature(-1, 0)
