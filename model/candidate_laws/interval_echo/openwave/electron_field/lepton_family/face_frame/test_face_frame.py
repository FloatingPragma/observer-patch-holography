import json

import numpy as np

from . import measure, observers


def test_field_spectroscopy_and_native_memory_reversal():
    r = measure.run()
    saved = json.loads((measure.HERE / "results.json").read_text())
    field = r["declared_field"]
    np.testing.assert_allclose(
        field["ordered_ratios"], [1, 206.770316, 3477.472837], rtol=2e-9
    )
    assert len(field["joint_field"]["field_banks"]) == 8
    for result in field["banks"]:
        d = result["long_clock"]
        assert max(m["simultaneous_lower"] for m in d["hankel_minors"]) > 0.77
        assert d["fits"][2]["held_rms"] < 0.002 < 0.5 < d["fits"][1]["held_rms"]
    for actual, prior in zip(
        r["native_written_face"]["forms"],
        saved["native_written_face"]["forms"],
        strict=True,
    ):
        assert actual["global_energy_enclosure"] == prior["global_energy_enclosure"]
        assert actual["rows"][0]["conjugate_spectral_distance"] < 1e-10
        assert actual["rows"][1]["conjugate_spectral_distance"] > 0.70
        out = np.array([row["return_probability_time_one"] for row in actual["rows"]])
        np.testing.assert_allclose(out[1], out[2, [0, 2, 1]], atol=2e-13)
        assert min(min(row["leakage_squared"]) for row in actual["rows"]) > 1


def test_independently_calibrated_observers_predict_charge_and_energy():
    r = observers.run()
    f = r["declared_electroweak_and_colour"]
    stats = f["calibration"]["held"]
    assert len(stats) == 24 and all(x["overlap"] for x in stats)
    assert not any(x["untransported_frame_fits"] for x in stats)
    for row in f["joint_rest_energy_and_weak_read"]:
        assert row["untransported_total_variation"] > 0.12
        for bank in row["two_receiver_record_banks"]:
            assert sum(bank) == r["shots_per_bank"]
    for row in r["native_local_port_frames"]:
        assert (
            abs(row["agreed_probability"] - row["untransported_detector_probability"])
            > 0.02
        )
