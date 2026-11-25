"""Native counter provenance, finite spectral bound and held original-start trials."""

import json

import networkx as nx
import numpy as np

from ....apparatus import fixture
from .detectors import verify
from .flight import chart, frame
from .geometry import layout
from .measure import HERE, SHORT, LARGE


def test_physical_release_and_ternary_receivers():
    ring, marks, leaves, owners, candidates = layout()
    assert (len(ring), len(owners), len(marks), len(leaves), candidates) == (
        14,
        5,
        2,
        2,
        35,
    )
    f = fixture()
    original = nx.Graph(f["blank"] + f["marked"])
    assert 167 not in original
    banks = verify()
    assert {tuple(x["receiver"]) for x in banks} == {(207, 206), (106, 96)}
    data = json.loads((HERE / "results.json").read_text())
    assert data["local_receiver_meters"] == banks
    assert all(
        row["original_pair_parking"][0] != row["original_pair_parking"][1]
        for row in data["setup"]
    )


def test_static_pair_band_and_orthogonal_side_memories():
    ring, *_ = layout()
    for initial in [(160, 167), (167, 160)]:
        sets = [set(chart((*initial, mask))[0]) for mask in range(4)]
        assert len(set.union(*sets)) == sum(map(len, sets))
    for use_heads in (True, False):
        stats, near, _ = frame(-1, 0, written=use_heads)
        assert stats["upper_dimension"] == 2 * len(ring)
        assert (
            stats["diagonal_exclusion"][0] <= 36 < 38 < stats["diagonal_exclusion"][1]
        )
        assert stats["analytic_all_time_floor"] > (0.68 if use_heads else 0.98)
        assert stats["high_energy_weight"] > 0.997
        assert min(near) > 0.995


def test_held_original_starts_and_controls():
    data = json.loads((HERE / "results.json").read_text())
    held = data["held"]
    assert len(held) == 48 and len(data["grid"]) == 216
    for _, _, _, _, _, adj, full, *_ in data["grid"]:
        assert sum(adj) == sum(full) == SHORT
    for row in held:
        a, b, c = (
            row["data"][k]
            for k in ["free", "independent_marginals", "two_gates_closed"]
        )
        assert sum(a["counts"]) == sum(b["counts"]) == sum(c["counts"]) == row["starts"]
        assert c["coincidences"] == 0
        if row["bench"] == "retired" or row["starts"] == LARGE:
            assert a["original_99"][0] > max(b["original_99"][1], c["original_99"][1])
        assert a["counts"][0] == a["counts"][8] == 0
    # In the original gate configuration, 167 has no free exit. With the
    # reconfigured closed cycle all six opposite pairs are resolved on fresh banks.
    np.testing.assert_array_less(data["original_source_to_neutral_probabilities"], 1e-6)
