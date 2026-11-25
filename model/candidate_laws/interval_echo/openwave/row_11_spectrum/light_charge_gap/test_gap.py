"""Independent fixtures for topology, native heads, energy zero and arrivals."""

import json
from pathlib import Path

import numpy as np
from . import measure


def test_native_edges_and_passive_head_costs():
    for sign, reference in [(-1, -8), (1, 24)]:
        blank, cost = measure.running_coefficients(sign)
        np.testing.assert_array_equal(blank, [[8.5, -0.5], [-0.5, 8.5]])
        assert cost == reference
        for parity in [-1, 1]:
            strength = 32 + (parity == -1)
            block = measure.local_sector(sign, (2, parity))
            np.testing.assert_array_equal(
                block, strength / 2 * np.kron(np.ones((2, 2)), np.eye(2))
            )
            np.testing.assert_array_equal(
                measure.pendant_block(sign, parity),
                -(8 + (parity == -1)) / 2 * np.ones((2, 2)),
            )


def test_parity_resolved_spectrum_and_compensator_controls():
    ports, store, writes, matrices = measure.bench(-1)
    assert len(set(ports + store)) == 6 and writes > 0
    h = matrices[True]
    np.testing.assert_array_equal(np.diag(h), [0.5, 17, 17, 0.5])
    reflection = np.eye(4)[::-1]
    np.testing.assert_array_equal(h @ reflection, reflection @ h)
    expected = sorted(
        e
        for central in [16.5, 17.5]
        for e in np.linalg.eigvalsh([[0.5, -0.5], [-0.5, central]])
    )
    np.testing.assert_allclose(np.linalg.eigvalsh(h), expected, atol=1e-13)
    np.testing.assert_array_equal(matrices[False][:2, 2:], 0)


def test_released_charge_and_checked_receipt():
    data = measure.run()
    receipt = json.loads(Path(measure.__file__).with_name("results.json").read_text())
    for field in ["ports", "store", "exterior_writes", "seed", "starts_per_bank"]:
        assert data[field] == receipt[field]
    for arm, held in zip(data["banks"], receipt["banks"], strict=True):
        assert sum(arm["counts"]) == 2**20 and arm["counts"] == held["counts"]
        np.testing.assert_allclose(
            arm["relative_levels"], held["relative_levels"], atol=1e-12
        )
        assert abs(arm["right_probability"] - held["right_probability"]) < 2e-8
        if not arm["connected"] or arm["fraction"] == 0:
            assert arm["counts"][-1] == 0
        elif arm["fraction"] == 1:
            assert arm["counts"][-1] > 0.95 * 2**20
    for sign, lower in [("-1", 3.5), ("1", 7)]:
        active = data["active_pendants"][sign]
        assert active["native_levels"] == 1536
        assert active["minimum_positive"] >= lower and active["largest_negative"] < -4
        assert (
            abs(
                active["minimum_positive"]
                - receipt["active_pendants"][sign]["minimum_positive"]
            )
            < 1e-12
        )


def test_absolute_meter_resolves_both_low_lines():
    from . import meter

    data = meter.read()
    saved = json.loads(Path(measure.__file__).with_name("meter.json").read_text())
    for row, old in zip(data["rows"], saved["rows"], strict=True):
        assert (
            row["port_counts"] == old["port_counts"]
            and row["odd_bond_counts"] == old["odd_bond_counts"]
        )
        assert sum(row["port_counts"]) == data["starts_per_bank"]
        assert row["simultaneous_99"][0] < row["exact"] < row["simultaneous_99"][1]
    assert data["rows"][0]["simultaneous_99"][1] < data["rows"][1]["simultaneous_99"][0]
