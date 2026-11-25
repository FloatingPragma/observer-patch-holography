"""Checks independent of the sampling seed for the common-event measurement."""

from fractions import Fraction as Q
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from background_clock import measure
from horizon import moments, exact_table

DATA = json.loads((HERE / "records.json").read_text())


def test_complete_crest_records_and_density():
    for level, row in enumerate(DATA["boards"]):
        times = row["clock"]
        for j, counts in enumerate(row["counts_per_face"]):
            q = times[j + 1] - times[j]
            assert len(counts) == 20 and counts == [q**3] * 20
            assert row["fixed_enrollment_per_layer"][j] == q
            assert counts[0] == q**4 if q == 1 else counts[0] < q**4
        assert row["parked_census_count"] == 0
        assert row["static_counts_per_face"] == [1] * 20
        assert sum(Q(c["probability"]) for c in row["trigger_mixture"]) == 1
        for c in [row["commanded"], *row["trigger_mixture"]]:
            e = c.get("trigger", 0)
            expected = [s for s in range(e, e + 7, 2) if s <= 36]
            assert c["emission_layers"] == expected
            for v, d in zip(row["receivers"], row["radii"], strict=True):
                obs = c["receivers"][str(v)]
                arrivals = [s + d for s in expected if s + d <= 36]
                assert (
                    obs["arrival_layers"] == arrivals
                    and obs["received"] + obs["missing"] == 4
                )
                assert obs["local_clocks"] == [times[a] for a in arrivals]
                ratios = [
                    Q(times[b] - times[a], times[b - d] - times[a - d])
                    for a, b in zip(arrivals, arrivals[1:])
                ]
                assert ratios == [Q(k) for k in obs["period_ratios"]]
        for bank in row["static"]["receivers"].values():
            assert bank["period_ratios"] == ["1", "1", "1"]
        assert all(v["received"] == 0 for v in row["cut"]["receivers"].values())
        muted = row["muted_receiver"]
        control = row["commanded"]
        assert muted["emitter_clock"] == control["emitter_clock"]
        for v in map(str, row["receivers"]):
            m, c = muted["receivers"][v], control["receivers"][v]
            assert m["arrival_layers"] == c["arrival_layers"]
            if v == str(row["receivers"][-1]):
                assert Q(m["period_ratios"][0]) == Q(c["period_ratios"][0]) - Q(1, 2)
            else:
                assert m == c
    largest = DATA["boards"][-1]
    assert [a["period_ratios"] for a in largest["commanded"]["receivers"].values()] == [
        ["1", "1", "2"],
        ["2", "2", "2"],
        ["4", "4", "4"],
    ]


def test_shared_seam_anchors_and_interventions():
    for board in DATA["boards"]:
        banks = board["comparison"]
        native = banks["native"]["events"]
        for m, case in banks.items():
            assert len(case["attempts"]) == 37
            assert [r[:3] for r in case["attempts"] if r[3:] == [1, 1]] == case[
                "events"
            ]
            assert case["quiet"] + case["invalid"] + len(case["events"]) == 37
            if m == "cut":
                assert case["invalid"] == 4 and not case["events"]
                continue
            (_, x, y), (_, u, v) = case["events"][:2]
            for (_, a, b), residual in zip(
                case["events"][2:], case["calibration"]["residual"], strict=True
            ):
                assert Q((b - y) * (u - x) - (a - x) * (v - y), u - x) == Q(residual)
        assert native[:2] == banks["drift"]["events"][:2]
        assert banks["drift"]["calibration"]["residual"][-1] == "-1"
        assert banks["native"]["calibration"]["residual"] == ["0", "0"]
        assert banks["native"]["calibration"]["rate"] == "1/2"
        # The purposely wider family of all quiet and positive contacts mixes sampling phases.
        assert any(
            Q(r[2]) != Q(r[1], 2) + Q(banks["native"]["calibration"]["offset"])
            for r in banks["native"]["attempts"]
        )


def test_same_clock_thermal_contrast():
    th = DATA["boards"][-1]["thermal"]
    clock = DATA["boards"][-1]["clock"]
    assert th["marker"] - th["reset_marker"] == 8 and not th["shared_postreset_owners"]
    a, b, c, d = map(set, th["pulse_inputs_per_probe"])
    assert not (a | b) & (c | d)
    reference = sum(
        (-1) ** (i + j) * moments(exact_table("release", 25, u, v))[-1]
        for i, u in enumerate([3, 4])
        for j, v in enumerate([3, 4])
    )
    assert reference == Q(th["layer"]["covariance"])
    for mode in ("layer", "microstep"):
        bank = th[mode]
        total = np.zeros(16, dtype=np.int64)
        for row in bank["ages"]:
            t = row["trigger"]
            n = [min(25, t + d) for d in (3, 4, 3, 4)]
            assert row["offers"] == (n if mode == "layer" else [clock[i] for i in n])
            assert row["misses"] == [t + d > 25 for d in (3, 4, 3, 4)]
            words = np.array(row["energy_words"])
            assert np.all(words >= 0)
            total += words
        assert (
            total.tolist() == bank["energy_words"]
            and total.sum() == measure.CONFIG["contrast_originals"]
        )
        words = np.arange(16)
        u = (words & 1) - ((words >> 1) & 1)
        v = ((words >> 2) & 1) - ((words >> 3) & 1)
        prob = total / total.sum()
        est = prob @ (u * v) - (prob @ u) * (prob @ v)
        assert abs(est - bank["observed"]) < 1e-15
        assert (
            abs(est - float(Q(bank["covariance"])))
            < bank["simultaneous_99pct_halfwidth"]
        )


def test_native_replay_preserves_exact_observables():
    fresh = measure.run()
    for old, new in zip(DATA["boards"], fresh["boards"], strict=True):
        assert old["comparison"] == new["comparison"]
        assert old["counts_per_face"] == new["counts_per_face"]
        assert old["trigger_mixture"] == new["trigger_mixture"]
        if "thermal" in old:
            for bank in ("layer", "microstep"):
                assert (
                    old["thermal"][bank]["covariance"]
                    == new["thermal"][bank]["covariance"]
                )
