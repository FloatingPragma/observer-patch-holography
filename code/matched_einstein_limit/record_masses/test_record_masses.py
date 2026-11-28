"""Local clocks, energy references and identifying the exterior from interval records."""

from copy import deepcopy
import json
from math import log, sqrt, tan

import numpy as np
import pytest

from . import exterior, guide, reference


def receipt(name):
    return json.loads((guide.HERE / (name + ".json")).read_text())


def test_two_readers_kinetic_curvature_and_origin_invariance():
    data = guide.run()
    saved = receipt("guide")
    for label, x in data["guides"].items():
        assert x["tapes"] == saved["guides"][label]["tapes"]
        delta = x["audit"]["analytic_gap"]
        assert x["audit"]["max_gate_phase_residual"] < 2e-12
        assert x["audit"]["max_preparation_error"] < 2e-12
        assert abs(x["audit"]["relative_bias"]) < 0.00012
        for t, r in zip(x["tapes"], x["decoded"], strict=True):
            assert r["rest_half_gap"][0] < delta < r["rest_half_gap"][1]
            for sense in r["senses"].values():
                assert sense["kinetic_mass"][0] < tan(delta) < sense["kinetic_mass"][1]
            shifted = dict(
                t,
                hand=[v + 951 for v in t["hand"]],
                echo_hand=[v + 951 for v in t["echo_hand"]],
            )
            assert guide.decode(shifted) == r
            rescaled = dict(
                t,
                hand=[v * 2 for v in t["hand"]],
                echo_hand=[v * 2 for v in t["echo_hand"]],
            )
            new = guide.decode(rescaled)
            np.testing.assert_allclose(
                np.array(new["rest_half_gap"]) * 2, r["rest_half_gap"]
            )
            np.testing.assert_allclose(
                new["senses"]["1"]["kinetic_over_gap"],
                r["senses"]["1"]["kinetic_over_gap"],
            )


def test_parked_reference_reads_occupied_sector_phase():
    data = reference.run()
    saved = receipt("reference")
    for key, g in data["guides"].items():
        for label, row in g["rows"].items():
            assert row["tape"] == saved["guides"][key]["rows"][label]["tape"]
            contrast, sign = map(int, label.split("/"))
            expected = 4 * contrast + (sign < 0)
            r = row["reading"]
            lo, hi = r["rate_per_empty_clock"]
            assert lo < expected < hi
            assert r["turns"] > 10
        no_arm = np.array(g["no_reference_control"])[:, 0] / reference.SHOTS
        assert np.hypot(no_arm[0] - no_arm[2], no_arm[3] - no_arm[1]) < 2 * np.sqrt(
            np.log(48000 / 0.01) / reference.SHOTS
        )


def test_three_clock_inverse_without_radii_and_held_material_read():
    data = exterior.run()
    saved = receipt("exterior")
    rows = data["sources"]
    for tag, row in rows.items():
        t = row["clock_tape"]
        assert t == saved["sources"][tag]["clock_tape"]
        shift = deepcopy(t)
        for j, link in enumerate(shift["links"]):
            link["emission"] = [i + 17 for i in link["emission"]]
            link["return_tick"] += 17
            link["reception"] = [i - 114 * j for i in link["reception"]]
        assert exterior.decode(shift) == exterior.decode(t)
        if tag == "blank":
            assert not row["inferred"]["resolved"]
            continue
        result = row["inferred"]
        expect = 19 / 132 if tag == "held" else 1 / 12
        assert result["mu_interval"][0] < expect < result["mu_interval"][1]
        for n, (a, b) in enumerate(result["radius_intervals"], 2):
            assert a < n < b
        step = 19 / 11 if tag == "light_alone_collar_blind" else 1
        assert all(abs(v - step) < 1e-7 for v in result["relative_free_fall"].values())
        h = exterior.spans(t)[0]
        for alternative in result["two_station_alternatives"]:
            m = alternative["mu"]
            a, b = alternative["radii"]
            assert a > 1 and abs(m - expect) > 0.01
            z = sqrt((1 - 2 * m / a) / (1 - 2 * m / b))
            dt = (
                2
                * sqrt(1 - 2 * m / a)
                * (b - a + 2 * m * log((b - 2 * m) / (a - 2 * m)))
            )
            assert abs(z - h[0] / h[1]) < 1e-12 and abs(dt - t["tick"] * h[2]) < 1e-12
    assert (
        rows["held"]["word_counts"] == rows["light_alone_collar_blind"]["word_counts"]
    )
    assert rows["ordinary"]["one_end_modular_change_nats"] == 0
    assert abs(data["mu_forecast_from_single_pair"]["held"] - 19 / 132) < 1e-8
    assert abs(data["held_active_mass_in_reference_units"] - 19 / 33) < 1e-8


def test_missing_link_or_failed_guide_echo():
    t = receipt("exterior")["sources"]["ordinary"]["clock_tape"]
    t["links"].pop()
    with pytest.raises(ValueError, match="Two separately"):
        exterior.decode(t)
    t = receipt("guide")["guides"]["27"]["tapes"][0]
    t["echo_hits"][0] += t["echo_hits"][-1]
    t["echo_hits"][-1] = 0
    with pytest.raises(ValueError, match="echo"):
        guide.decode(t)
