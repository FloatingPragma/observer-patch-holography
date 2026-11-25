import json

import numpy as np

from . import reciprocal


def test_opposite_sources_and_retained_original_starts():
    stored = json.loads((reciprocal.HERE / "reciprocal.json").read_text())
    result = reciprocal.run()
    assert result["terminal_ports"] == [4, 160]
    assert result["largest_Choi_discrepancy"] < 1e-13
    assert len(result["native_lap_controls"]) == 6
    for direction in ("forward", "reciprocal"):
        measure = result[direction]
        assert measure["seed"] == stored[direction]["seed"]
        assert measure["observed_upper99"] < -6e-8
        assert 1.5e-7 < measure["all_launches_code_acceptance"] < 1.6e-7
        assert 2.1e-6 < measure["bare_remote_arrival"] < 2.3e-6
        assert measure["bare_local_return"] > 0.98
        np.testing.assert_allclose(
            measure["observed_PT"], stored[direction]["observed_PT"], atol=3e-10
        )
        for bank in measure["observations"]:
            assert sum(bank["bare_ports_home_remote_else_refused"]) == reciprocal.SHOTS
            assert all(sum(row) == reciprocal.SHOTS for row in bank["settings"])
    assert result["forward"]["seed"] != result["reciprocal"]["seed"]
    # The independently diagonalized symmetry-reduced forward propagator
    # gives the same code aperture, including both local-erasure probabilities.
    other = json.loads((reciprocal.HERE / "results.json").read_text())["autonomous"][
        "32"
    ]
    np.testing.assert_allclose(
        other["original_PT_minimum"],
        result["forward"]["all_launches_PT_minimum"],
        rtol=2e-8,
        atol=1e-14,
    )
    np.testing.assert_allclose(
        other["ideal_total_yield"],
        result["forward"]["all_launches_code_acceptance"],
        rtol=2e-8,
        atol=1e-14,
    )
