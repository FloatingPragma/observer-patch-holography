"""Exact native leading terms and the occupied metric and relabeling certificates."""

import json
import numpy as np
from scipy.sparse.csgraph import shortest_path

from . import measure


def test_all_one_cell_leading_hamiltonian_powers():
    b = measure.monitored.federation(0)
    _, contacts, _, words = measure.monitored.apparatus(b, [0])
    h = (2 * measure.monitored.generator(b, words, contacts)).astype(np.int64)
    hops = (h < 0).astype(np.int64)
    d = shortest_path(hops, directed=False, unweighted=True)
    native = walks = np.eye(len(words), dtype=object)
    for k in range(int(d.max()) + 1):
        assert np.all(native[d > k] == 0)
        assert np.all(walks[d == k] > 0)
        np.testing.assert_array_equal(native[d == k], (-1) ** k * walks[d == k])
        native = h @ native
        walks = hops @ walks


def test_two_inequivalent_rooted_boundaries_and_exact_entropy_degeneracy():
    data = json.loads(measure.RECORD.read_text())
    for row in data["rows"]:
        n = len(row["owners"])
        assert row["unoccupied"]["best_gap"] == -1
        assert (
            sum(row["protected_word_gap_counts"].values())
            == row["monitoring"]["complete_words"]
        )
        assert (
            sum(n for _, n in row["branching_witnesses_by_resident_word_count"])
            == row["monitoring"]["complete_words"]
        )
        valid = row["possible_protected_balls"]
        if n > 1:
            assert not valid and max(map(int, row["protected_word_gap_counts"])) == 0
            continue
        assert row["any_center_equilibrium_probability"] == "1/11"
        assert row["any_center_outward_rejection_equilibrium_probability"] == "3/67"
        assert len({v["center"] for v in valid}) == len(valid) == 6
        for v in valid:
            k = v["cut_occupancy"]
            center, *rest = v["branching_center_incoming_outgoing"]
            assert center == v["center"] and len(set([center, *rest])) == 4
            assert (
                v["volume_ports"] == 10 and v["equilibrium_interior_multiplicity"] == 3
            )
            shell, gates, reference = (
                ([1, 4, 5], [2, 2], [1, 2, 2])
                if k == 1
                else ([1, 3, 5, 1], [3], [1, 1, 3])
            )
            assert (
                v["shell_counts"],
                v["exposed_rim_radii"],
                v["unoccupied_center_to_gates"],
            ) == (shell, gates, reference)
            same = [a for a in valid if a["cut_occupancy"] == k]
            assert {tuple(p["word"]) for p in v["local_rooted_cut_isometries"]} == {
                tuple(s["resident_word"]) for s in same
            }


def test_native_replay():
    assert measure.run() == json.loads(measure.RECORD.read_text())
