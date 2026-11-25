import base64
import itertools
import json
from math import sqrt
import zlib

import networkx as nx
import numpy as np

from . import measure as m


def test_full_marked_echo_and_reference_interactions():
    law = m.candidate.native_law
    ex, contrast = law.local_term_parts()
    h = ex + 8 * contrast
    choice = [
        9 + 9 * i + j for i, (_, d) in enumerate(law.LABELS) if d == 0 for j in range(9)
    ]
    np.testing.assert_allclose(h[:, choice], 0, atol=1e-13)
    np.testing.assert_allclose(h[:9, :9], law.seam_term(8), atol=2e-13)
    np.testing.assert_allclose(h[9:, :9], 0, atol=2e-13)
    # An echo on a filled, blank seam writes either zero-contrast coin,
    # with unchanged endpoint loads; a neighboring write can use them again.
    filled = law.SRC.index((1, 1))
    support = set(np.flatnonzero(abs(law.echo()[:, filled]) > 1e-12))
    assert support == {j for j in choice if (j - 9) % 9 == filled}
    archived = json.loads((m.HERE / "results.json").read_text())
    assert archived["local"] == m.seam_certificate()
    for level, row in enumerate(archived["boards"]):
        board = m.federation(level)
        assert m.path_fixture(board) == row["live_path"]
        assert len(board.seams) - 2 == row["marked_seams"]
        assert sorted(board.degree[row["live_path"]].tolist()) == [5] * 3
        assert row["gateway_graph_is_forest"]
        assert row["maximum_gateway_degree"] == 2
        assert row["maximum_gateway_component"] == [6, 6, 12][level]


def test_bound_and_strict_ground_sign_on_all_small_graph_types():
    """Connected graph atlas with all physical degree-five/six assignments."""
    floor = (19 - sqrt(297)) / 4
    for graph in nx.graph_atlas_g():
        if not 1 <= len(graph) <= 5 or not nx.is_connected(graph):
            continue
        adjacency = nx.to_numpy_array(graph)
        live = adjacency.sum(axis=1)
        lap = np.diag(live) - adjacency
        for private in itertools.product((0, 1), repeat=len(graph)):
            p = np.array(private)
            vacancy, excess = [
                np.linalg.eigvalsh(lap / 2 + np.diag(8 * (live + c * p)))
                for c in (-1, 3)
            ]
            assert np.all((vacancy < 1e-10) | (vacancy >= floor - 1e-10))
            assert np.all((excess < 1e-10) | (excess >= 8 - 1e-10))
            assert max(vacancy) <= 54 + 1e-10 and max(excess) <= 69 + 1e-10
            assert (abs(excess[0] - vacancy[0]) < 1e-10) == (sum(private) == 0)
    leaf = m.local_matrix(m.federation(0), m.path_fixture(m.federation(0)), -1, 8)
    np.testing.assert_allclose(
        np.linalg.eigvalsh(leaf),
        [(19 - sqrt(297)) / 4, 0.5, (19 + sqrt(297)) / 4],
        atol=2e-13,
    )


def unpack(bank):
    raw = zlib.decompress(
        base64.b64decode(bank["counts_u32_le_zlib_b64"], validate=True)
    )
    return np.frombuffer(raw, dtype="<u4").reshape(bank["counts_shape"])


def test_all_4096_position_banks_against_even_and_odd_path_amplitudes():
    archived = json.loads((m.HERE / "results.json").read_text())
    for bank in archived["banks"]:
        strength = 8 if bank["law"] == "native" else 4
        sign = bank["sign"]
        shift = 4 * strength if sign == 1 else 0
        # At the leaf reflection the odd mode occupies the two ends. The even
        # leaf sum and hub form this two-dimensional independent generator.
        small = np.array(
            [[0.5, -1 / sqrt(2)], [-1 / sqrt(2), strength + 1]]
        ) + shift * np.eye(2)
        eig, v = np.linalg.eigh(small)
        even = (np.exp(-1j * np.outer(m.TIMES, eig)) * v[0]) @ v.T / sqrt(2)
        odd = np.exp(-1j * m.TIMES * (0.5 + shift)) / sqrt(2)
        analytic = np.column_stack(
            (
                abs(even[:, 0] + odd) ** 2 / 2,
                abs(even[:, 1]) ** 2,
                abs(even[:, 0] - odd) ** 2 / 2,
                np.zeros(len(m.TIMES)),
            )
        )
        counts = unpack(bank)
        assert counts.shape == (1024, 4) and (counts.sum(axis=1) == m.SHOTS).all()
        assert not counts[:, 3].any()
        np.testing.assert_allclose(counts / m.SHOTS, analytic, atol=0.015)
        fast = m.read_frequencies(counts[:, 1] / m.SHOTS)
        slow = m.read_frequencies((counts[:, 0].astype(float) - counts[:, 2]) / m.SHOTS)
        np.testing.assert_allclose(
            [fast, slow],
            [bank["fitted_hub_gap"], bank["fitted_leaf_imbalance_gap"]],
            atol=2e-8,
        )
        np.testing.assert_allclose(
            [fast, slow], [eig[1] - eig[0], 0.5 + shift - eig[0]], atol=5e-5
        )
