import json

import numpy as np

from ..action_response.four_ticks import parse_interval
from . import bridge, improvement, phase_read


def test_native_bridge_conservation_and_four_tick_sign():
    table = bridge.report()
    certified = json.loads((bridge.HERE / "certificate.json").read_text())
    for q, groups in table["charges"].items():
        mixture = groups["equal_tagged_launch"]
        for key, value in mixture.items():
            np.testing.assert_allclose(
                value,
                (groups["forward_launch"][key] + groups["reverse_launch"][key]) / 2,
                atol=1e-12,
            )
        assert mixture["reference"] in parse_interval(certified[q]["reference"])
        assert mixture["ordinary"] in parse_interval(certified[q]["ordinary"])
        assert mixture["pairing"] / mixture["ordinary"] in parse_interval(
            certified[q]["pairing_over_ordinary"]
        )


def test_phase_records_and_changed_flight():
    data = phase_read.run()
    nominal = json.loads((bridge.HERE / "bridge.json").read_text())["charges"]
    for q, result in data["charges"].items():
        active = result["exact_mean_predictions"]["active"]
        changed = result["exact_mean_predictions"]["zero_contrast_strength"]
        for name, out in result["slopes"].items():
            np.testing.assert_allclose(
                active[name], nominal[q]["equal_tagged_launch"][name], atol=1e-11
            )
            assert out["estimate"] + out["simultaneous_halfwidth"] < 0 < changed[name]


def test_collar_witness_and_conserved_pressure():
    result = improvement.run()
    source = result["source"]
    contrast = parse_interval(source["matched_pairing_response"])
    assert source["matched_ordinary_response_exact"] == 0
    assert contrast.a > 0.20029 and contrast.b < 0.20035
    estimate = parse_interval(result["simulated_counts"]["contrast"])
    assert source["nominal_source"]["pairing"] in estimate
    assert result["finite_chart"]["quadrupole"]["local_energy_at_boost_half"] == "1/3"
