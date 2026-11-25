"""Reproduce the geometric effect, signed Ramsey clocks and free-release falsifier."""

import json
import numpy as np
from . import circuit
from .measure import HERE, SHOTS, run


def test_native_feedback_and_simultaneous_rates():
    fresh = run()
    saved = json.loads((HERE / "readings.json").read_text())
    for data in (fresh, saved):
        assert len(set(sum(data["physical_module"], []))) == 12
        assert (
            data["material"]["reference_port"]
            not in data["material"]["pointer_references"]
        )
        p = np.array(data["static_probabilities"])
        assert abs(p.sum() - 1 / 3) < 1e-12
        assert abs(p.sum() / np.sqrt(p).sum() ** 2 - 2 / 3) < 1e-12
        for panel in data["shared_dwell_banks"]:
            eps = panel["common_exposure"]
            expected = -np.angle(1 + np.expm1(-1j * eps) * p)
            for j, character in enumerate(panel["characters"]):
                assert (
                    character["phase_rate99"][0]
                    < expected[j]
                    < character["phase_rate99"][1]
                )
                assert (
                    character["recovered_effect99"][0]
                    < p[j]
                    < character["recovered_effect99"][1]
                )
                for row in character["observations"]:
                    assert all(sum(count) == SHOTS for count in row["counts"])
            assert panel["recovered_effect_ratio99"][1][0] > 206.7682827
            assert (
                panel["recovered_effect_ratio99"][1][0]
                < 206.770316
                < panel["recovered_effect_ratio99"][1][1]
            )
        a, b = data["shared_dwell_banks"]
        assert a["measured_phase_ratio99"][2][1] < b["measured_phase_ratio99"][2][0]
        assert np.max(abs(np.asarray(data["phase_pi_control_per_tick"]))) < 1e-12
        assert data["unforced_release"]["grams_times_four"] == [
            [5, 2, 2],
            [2, 5, 2],
            [2, 2, 5],
        ]
        assert max(data["unforced_release"]["survival_on_three_ports_at_one"]) < 0.565
    np.testing.assert_allclose(
        circuit.tick(-1 / 128), np.conj(circuit.tick(1 / 128)), rtol=1e-12
    )
    np.testing.assert_allclose(
        saved["weak_limit_ratios"], fresh["weak_limit_ratios"], rtol=1e-10
    )


def test_unitary_final_only_clock():
    from . import coherent

    fresh = coherent.run()
    saved = json.loads((HERE / "coherent.json").read_text())
    for data in (saved, fresh):
        assert len(data["all_physical_phases"]) == 12
        assert (
            sum(abs(np.exp(1j * np.asarray(data["all_physical_phases"])) + 1) < 1e-12)
            == 6
        )
        for panel in data["characters"]:
            assert panel["low_mode_weight"] > 0.999996
            assert (
                panel["observed_rate99"][0]
                < panel["low_radians_per_tick"]
                < panel["observed_rate99"][1]
            )
            assert all(
                sum(counts) == SHOTS
                for row in panel["observations"]
                for counts in row["counts"]
            )
        for interval, rate in zip(
            data["measured_low_ratios"][1:], data["exact_low_ratios"][1:], strict=True
        ):
            assert interval[0] < rate < interval[1]
    np.testing.assert_allclose(
        fresh["exact_low_ratios"], saved["exact_low_ratios"], rtol=1e-9
    )
