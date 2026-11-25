import json

from . import measure


def test_committed_native_and_c_quadratures_and_held_counts():
    receipt = measure.run()
    assert receipt == json.loads((measure.HERE / "results.json").read_bytes())
    for label, truth in (("E_native", 0), ("C_control", 0.5), ("bypass", 1)):
        bank = receipt["banks"][label]
        assert abs(bank["inferred_real_imag"][0] - truth) <= 2 * receipt["simultaneous_99_percent_radius"]
        assert all(sum(p["counts"]) == measure.SHOTS and len(p["counts"]) == 78 for p in bank["panels"])
