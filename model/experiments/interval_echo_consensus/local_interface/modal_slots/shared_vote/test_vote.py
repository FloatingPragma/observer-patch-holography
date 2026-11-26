import json
import numpy as np
from . import measure


def test_native_archive_and_tie_thresholds():
    data = measure.acquire()
    assert json.loads(json.dumps(data)) == json.loads(
        (measure.HERE / "counts.json").read_text()
    )
    assert data["summary"] == {
        "sharp_triples": 2366,
        "minimum_private_edits": {"86": 23, "88": 20, "247": 45},
    }


def test_joint_frequencies_including_held_and_parked_pulses():
    data = json.loads((measure.HERE / "counts.json").read_text())
    eps = np.sqrt(np.log(2 * 48 * 4 / 0.01) / (2 * measure.SHOTS))
    for bank in data["banks"]:
        q = (
            np.sin(bank["pulse_pi_over_three"] * np.pi / 6) ** 2
            if bank["flight"] == "native"
            else 0
        )
        if bank["private_edits"] != "boundary":
            expected = [0, 0, 0, 1]
        elif bank["observers"] == [86, 88]:
            expected = [1 - q, 0, 0, q]
        else:
            expected = [(1 - q) ** 2, q * (1 - q), q * (1 - q), q**2]
        assert sum(bank["rail_counts"]) == sum(bank["joint_counts"]) == measure.SHOTS
        assert (
            max(abs(np.asarray(bank["joint_counts"]) / measure.SHOTS - expected)) < eps
        )
