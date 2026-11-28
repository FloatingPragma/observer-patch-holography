"""Independent native timing and flight on the actual clock-equipped board."""

import json

import numpy as np

from . import measure as m
from . import neutral_clock as nc


def test_calibrated_neutral_hands_with_arbitrary_local_loads():
    clock = nc.clock_register()
    assert clock["register_levels"] == {"0": [0, 0, 1], "8": [0, 32, 33]}
    assert clock["verified_native_columns"] == 27
    assert clock["all_magnitude_two_memories_in_single_source_columns"] == 12
    np.testing.assert_allclose(
        clock["hand_calibration"]["8"]["ramsey_yes"][:3], [1, 0, 1], atol=2e-14
    )
    np.testing.assert_allclose(
        clock["hand_calibration"]["0"]["ramsey_yes"], 1, atol=2e-14
    )
    assert clock["source_induced_frequency_change"] == [0, 0]


def test_instrumented_native_board_and_simulated_clocks_replay():
    data = nc.run()
    m.fp.check(
        json.loads((m.HERE / "neutral_clock.json").read_text()),
        json.loads(json.dumps(data)),
    )
    material = data["material"]
    assert material["clock_edges_near_far_by_owner"] == [
        (0, 1),
        (2, 3),
        (14, 15),
        (12, 13),
    ]
    assert material["ordinary_raw_increment"] == 255 / 8
    assert material["ordinary_rematched_probability"] == 1 / 15
    for name, sign in [
        ("target", -1),
        ("phase_feedback", 1),
        ("exchange_1_1", 1),
        ("energy_rematched_mass", 1),
    ]:
        stock = material["standards"][name]
        assert abs(stock["initial_material_energy_increment"] - 17 / 8) < 1e-13
        for flight in stock["launches"]:
            assert sign * flight["radial"][2] > 0.006
    sample = data["pointer_counts"]
    for table in sample["records"].values():
        for reference in table:
            for key, counts in reference["counts"].items():
                yes, no = np.array(counts).T
                np.testing.assert_array_equal(yes + no, sample["shots"])
                assert (
                    max(
                        abs(
                            yes / sample["shots"]
                            - np.array(data["clock"]["hand_calibration"]["8"][key])
                        )
                    )
                    < sample["individual_bin_simultaneous_99_half_width"]
                )
