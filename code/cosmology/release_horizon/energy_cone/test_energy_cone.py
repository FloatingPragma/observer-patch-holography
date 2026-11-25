"""Independent checks on the saved work readouts, energy words and full input envelope."""

from fractions import Fraction as F
import json
from math import atan2, pi, prod
from pathlib import Path

D = json.loads(Path(__file__).with_name("records.json").read_text())


def covariance(w):
    means = [0, 0, 0]
    for k, n in enumerate(w):
        a = (k & 1) - ((k >> 1) & 1)
        b = ((k >> 2) & 1) - ((k >> 3) & 1)
        for i, v in enumerate((a, b, a * b)):
            means[i] += n * v / sum(w)
    return means[2] - means[0] * means[1]


def test_conditional_density_and_original_denominators():
    for w in D["energy_cone"]["windows"]:
        h = w["thermal"]
        weights = [F(0)] * 16
        for t, p, offers in h["trigger_strata"]:
            # Evolve the two local occupations on each independent contact,
            # starting on the Gibbs populations of the probe and reservoir.
            means = []
            for count in offers:
                a, b = F(1, 5), F(1, 3)
                for _ in range(count):
                    a, b = (3 * a + b) / 4, (a + 3 * b) / 4
                means.append(a)
            for i in range(16):
                weights[i] += F(p) * prod(
                    q if i & (1 << k) else 1 - q for k, q in enumerate(means)
                )
        fired, parked = h["root_interventions"]
        assert fired["root_launch"] == 1 and parked["root_launch"] is None
        assert all(1 <= t <= w["reset_layer"] for t in fired["first_crest"])
        assert all(x is None for x in parked["first_crest"])
        assert parked["contrast_means"] == ["0", "0"]
        assert list(map(F, fired["contrast_means"])) == (
            [F(1, 240), F(1, 512)]
            if w["reset_layer"] == 24
            else [F(1, 7680), F(7, 480)]
        )
        assert sum(weights) == 1 and weights == list(map(F, h["exact_energy_words"]))
        assert covariance(weights) == F(h["exact_covariance"])
        assert covariance(list(map(F, h["independent_trigger_words"]))) == 0
        sampled = list(map(F, h["sampled_energy_words"]))
        assert sum(sampled) == 2**50
        observed = float(covariance(sampled))
        assert observed == h["observed_covariance"]
        assert (
            abs(observed - float(covariance(weights)))
            < h["simultaneous_99_halfwidth"]
            < observed
        )


def test_boundary_injections_reach_and_parked_control():
    c = D["energy_cone"]
    assert c["single_rail_leakage"] == [0.5] * 4
    assert [r["distance"] for r in c["box_axis_returns"]] == [4] * 3
    assert [r["distance"] for r in c["cross_ranges"]] == [3] * 4
    ticks = D["rail_hop_ticks"]
    for w in c["windows"]:
        frames = 25 - w["reset_layer"]
        assert w["clock_strokes"] == [100 * ticks] * frames
        assert w["publications"] == [256] * frames
        assert w["parked_origins"] == c["probes"]
        assert w["marker_readings"] == [
            [60 - 4 * k] * 4 for k in reversed(range(frames + 1))
        ]
        for i, v in enumerate(c["probes"]):
            y = (v // 5) % 5
            y_source = v + (5 if y % 2 == 0 else -5)
            expected = [v, y_source]
            if frames == 2:
                expected.append(y_source + (25 if y_source // 25 % 2 == 0 else -25))
            assert expected[-1] == w["read_origins"][i]
            assert sorted(expected) == [r["start"] for r in w["impulse_ranges"][i]]
            assert all(
                r["return_strokes"] == 2 * ticks * r["distance"]
                for r in w["impulse_ranges"][i]
            )
        assert w["endpoint_radii"] == [frames] * 4
        assert w["separation_minus_radii"] == [3 - 2 * frames] * 4


def test_recover_dilation_and_work_from_quadratures_and_returns():
    w = D["work"]
    ref = w["reference_ticks"]
    for row in w["rows"]:
        widths = [F(t, 2 * ref) for t in row["return_strokes"]]
        assert prod(widths) == F(row["relative_volume"])
        assert len(set(widths)) == 1
        for family, target in (("fixed_action", 1 / widths[0]), ("fixed_duty", 1)):
            bank = row[family]
            z, q = bank["quadratures"]
            measured = 3 * ref * atan2(2 * q - 1, 1 - 2 * z) / (pi * row["cycle_ticks"])
            assert abs(measured - bank["gap_in_reference_units"]) < 3e-12
            assert abs(measured - float(target)) < 3e-12
            assert bank["energy_outcomes"] == [[1, 0], [0, 1]]
    assert w["rows"][4]["fixed_action"] == w["rows"][4]["fixed_duty"]
    for sec in w["secants"]:
        d = sec["half_offset"] / ref
        upper, constant = sec["upper_pressures"]
        assert abs(upper - 1 / ((1 - d * d) * (3 + d * d))) < 4e-11 and constant == 0
    assert w["compression_from_eight_to_one"] == {
        "fixed_action": [0, 0.1, 0.5],
        "fixed_duty": [0, 0, 0],
    }
