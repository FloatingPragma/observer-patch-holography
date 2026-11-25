"""Check the receiver against independent standards and the untouched input counts."""

import json
import numpy as np
from . import measure


def test_contact_on_archived_arrivals():
    fresh = measure.run()
    saved = json.loads((measure.HERE / "results.json").read_text())
    assert fresh["summary"] == saved["summary"]
    assert fresh["held"] == saved["held"]
    for sensor, site in fresh["sites"].items():
        assert tuple(site["receiving"]) == measure.SENSORS[sensor]
        assert site["independent_standards"]["stocked_words_and_signed_head"] == [
            [0, 2, 2],
            [2, 0, -2],
        ]
        assert set(site["persistent_flag"]).isdisjoint(
            site["workstation"] + site["receiving"]
        )
    for row in fresh["calibration"]:
        word = row["word"]
        g = row["guard"]
        target = word[0] - 1 if g == word[1] else 0
        p = sum(o["probability"] for o in row["outputs"] if o["regional_charge"] != 0)
        expected = (
            (
                1
                if row["mode"] == "live"
                else 0.5
                if row["mode"] == "half_contrast"
                else 0
            )
            if target != 0
            else 0
        )
        np.testing.assert_allclose(p, expected, atol=2e-13, rtol=0)
        for out in row["outputs"]:
            assert out["donor_change"] == -out["regional_charge"]
            assert out["regional_charge"] in (0, target)
    assert len(fresh["held"]) == 48
    for rb in fresh["held"]:
        counts = rb["input_counts"]
        assert sum(counts) == 2**46
        for setting, signals in rb["matched"].items():
            original = [counts[2], counts[6]]
            assert [j[-1] for j in signals["flag_and_outgoing"]] == original
            assert sum(signals["charges"]) == sum(original)
            for i, (f1, f2, first, second, donor, _) in enumerate(
                signals["flag_and_outgoing"]
            ):
                assert [f1, f2] == ([2, 0] if i == 0 else [0, 2])
                stock = int(setting[0]) if "blocked_donor" in setting else 1
                assert first + second - 2 + donor - stock == 0
        if rb["release_control"] == "two_gates_closed":
            assert counts[2] == counts[6] == 0
    total = fresh["summary"]
    assert total["free/2/live"][:3] == [80942914, 80938876, 0]
    assert total["free/0/live"][:3] == [0, 80942914, 80938876]
    assert all(sum(counts) == 24 * 2**46 for counts in total.values())
