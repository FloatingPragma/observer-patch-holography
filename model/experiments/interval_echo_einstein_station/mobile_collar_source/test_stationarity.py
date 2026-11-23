"""Cross-check the new time perturbation against independent exposure-four acquisition."""

import hashlib
import json

import numpy as np

from . import stationarity as reader


def test_independent_exposure_four_reference_and_outcomes():
    packet = json.loads((reader.HERE / "stationarity.json").read_text())
    earlier = json.loads((reader.HERE / "results.json").read_text())
    assert (
        packet["protocol_sha256"]
        == hashlib.sha256(
            (reader.HERE / "STATIONARITY_PROTOCOL.md").read_bytes()
        ).hexdigest()
    )
    assert packet["blank_commutator_nonzero_entries"] == 0
    for arm in packet["arms"]:
        formation = next(
            b for b in earlier["experiments"] if b["cells"] == arm["cells"]
        )
        original = next(r for r in formation["arms"] if r["charge"] == arm["charge"])
        original_ci = reader.bin_intervals(
            np.ravel(original["modular"]["calibration_counts"])
        )
        pilot_ci = reader.bin_intervals(arm["contact"]["reference_counts"])
        assert np.all(
            np.maximum(pilot_ci[:, 0], original_ci[:, 0])
            <= np.minimum(pilot_ci[:, 1], original_ci[:, 1])
        )
        tv = [r["reference_tv_lower_99"] for r in arm["contact"]["reads"]]
        assert tv[2] == 0 and min(tv[:2] + tv[3:]) > 0.00305
        assert all(r["reference_tv_lower_99"] == 0 for r in arm["removed"]["reads"])
        assert arm["restored_material_head_energy"] == 16
        for name, exact in arm["exact_contact_entropy_derivatives"].items():
            for slope in arm["contact"]["symmetric_slopes"]:
                lo, hi = slope["entropies"][name]["interval"]
                assert lo > 0 if exact > 0 else hi < 0
            measured = [
                s["entropies"][name]["nats_per_exposure"]
                for s in arm["contact"]["symmetric_slopes"]
            ]
            assert abs(measured[1] - exact) < abs(measured[0] - exact)
            for slope in arm["removed"]["symmetric_slopes"]:
                lo, hi = slope["entropies"][name]["interval"]
                assert lo <= 0 <= hi
            assert arm["uniform_armed_sharp_control"]["entropy"][name]["nats"] == 0


def test_entropy_interval_including_rare_and_unanimous_counts():
    for counts in (
        [reader.SHOTS, 0],
        [1, reader.SHOTS - 1],
        [reader.SHOTS // 2, reader.SHOTS // 2],
    ):
        result = reader.entropy_read(np.array(counts))
        assert 0 <= result["interval"][0] <= result["nats"] <= result["interval"][1]
