import json
from types import SimpleNamespace

import numpy as np
from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    cold_release as cold,
    ordinary_release as ordinary,
)


def test_joint_columns_and_reduction_on_irregular_board():
    ordinary.check_joint_columns(8)
    board = SimpleNamespace(n=4, seams=[(0, 1), (0, 2), (0, 3), (2, 3)])
    h, double, single = ordinary.two_loads(board, 8, [(0, 1)])
    # Reassemble each configuration column directly from the local exchange
    # of its endpoint occupancies, plus the degree-completing private energy.
    for column, (i, j) in enumerate(zip(double, single, strict=True)):
        diag = -6 * 8 * 5
        expected = np.zeros(len(double))
        loads = np.zeros(4, dtype=int)
        loads[[i, j]] = [2, 1]
        degrees = np.bincount(np.asarray(board.seams).ravel(), minlength=4)
        diag += 8 * np.dot(6 - degrees, loads**2)
        for edge in board.seams[1:]:
            a, b = edge
            diag += 8 * (loads[a] - loads[b]) ** 2
            if loads[a] != loads[b]:
                diag += 0.5
                swapped = loads.copy()
                swapped[[a, b]] = swapped[[b, a]]
                row = (swapped[double] == 2) & (swapped[single] == 1)
                assert row.sum() == 1
                expected[row] -= 0.5
        expected[column] += diag
        np.testing.assert_array_equal(h[:, column].toarray().ravel(), expected)


def test_measured_masks_predict_flux_and_counts():
    record = json.loads((cold.HERE / "cold_release.json").read_bytes())
    board = cold.federation(1)
    for name, case in record["rows"].items():
        route, _, stiffness = name.split("_")
        mask, reading = cold.read_map(board, case["tape"])
        cold.apparatus.check(reading, case["modular_read"])
        m = np.bincount(np.array(mask).ravel(), minlength=board.n)
        boundary = sum(m[i] + m[j] for i, j in board.seams if bool(m[i]) != bool(m[j]))
        for q, run in case["probes"].items():
            assert np.isclose(
                run["mean_distance_acceleration"],
                -int(stiffness) * int(q) ** 2 * boundary / board.n,
            )
            np.testing.assert_allclose(
                run["inward_current_derivative"][1:], 0, atol=1e-12
            )
            if int(stiffness):
                assert (
                    min(run["localized"]["source_counts"]) / cold.STARTS
                    >= run["localized"]["all_time_source_probability_lower"]
                )
                diff = ordinary.difference(
                    ordinary.probability_band(run["cold"][0], cold.STARTS),
                    ordinary.probability_band(run["blank"][0], cold.STARTS),
                )
                assert diff[0, 0] > 0
    mixed = json.loads((cold.HERE / "ordinary_release.json").read_bytes())
    x = ordinary.decode(mixed)
    assert [0, 1] in x["opposite_signs"] and [0, 0] in x["resolved_nonadditivity"]
    np.testing.assert_array_equal(
        mixed["arms"]["joint"]["acceleration"],
        mixed["arms"]["ordinary"]["acceleration"],
    )
