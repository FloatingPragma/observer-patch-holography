"""Independent incidence, full native contact, region, and charge checks."""
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import expm
from scipy.sparse.csgraph import shortest_path

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from ..collar_pairing.measure import compare
from . import clock_balls, measure


def test_bundled_cell_neighbors_are_the_uniform_triangle_duals(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[4] / "code"))
    from source_selection_model.geometry import tower

    for layer in tower(3):
        owners = {}
        dual = nx.Graph()
        for index, face in enumerate(layer["faces"]):
            for edge in combinations(sorted(face), 2):
                if edge in owners:
                    dual.add_edge(index, owners.pop(edge))
                else:
                    owners[edge] = index
        assert not owners
        native = nx.Graph()
        native.add_edges_from((a // 12, b // 12) for a, b in federation(layer["level"]).external)
        assert nx.vf2pp_is_isomorphic(dual, native)


def test_both_probe_transfers_and_changed_double_load_exchange():
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    h = exchange + candidate.KAPPA * contrast
    # The changed law has the same addresses, reversal symmetry and one-load hand.
    # In the blank total-two sector add one more copy of the exchange generator.
    sector = np.array([int(a + b == 2) for a, b in law.SRC])
    changed = h.copy()
    changed[:law.N_SRC, :law.N_SRC] += np.diag(sector) @ exchange[:law.N_SRC, :law.N_SRC]
    for source, target in [((1, 0), (0, 1)), ((2, 0), (0, 2))]:
        launch, arrival = (np.eye(law.DIM)[:, law.SRC.index(x)] for x in (source, target))
        for generator, expected in [(h, arrival), (changed, arrival if source[0] == 1 else launch)]:
            np.testing.assert_allclose(expm(-1j * np.pi * generator) @ launch,
                                       expected, rtol=0, atol=2e-13)


def test_full_ball_records_against_port_graph_and_global_word_law():
    data = clock_balls.run()
    compare(data, json.loads((measure.HERE / "clock_balls.json").read_text()))
    local = data["local_gram"]
    cosine = 1 / np.sqrt(5)
    expected = (1 - cosine) * np.sqrt(1 + 2 * cosine) / 6
    assert local["tetrahedra"] == 180
    np.testing.assert_allclose([local["minimum_volume"], local["maximum_volume"]], expected, atol=2e-14)
    for bank in local["home_readings"]:
        np.testing.assert_array_equal(np.sum(bank["counts_by_launch"], axis=1), [2**20] * 12)
    laws, _, _ = measure.source_reader()
    reverse, fresh, shared = laws
    assert np.isclose(np.sum((shared - fresh) * np.log(fresh / reverse)),
                      data["ideal_source_contrast_nats"], atol=1e-15, rtol=0)
    assert Fraction(data["rational_contrast_lower"]) == Fraction(101, 33696)
    assert data["ideal_source_contrast_nats"] > float(Fraction(101, 33696))
    np.testing.assert_allclose(np.sum((shared - fresh) * (-np.log(fresh))), 0., atol=1e-15)
    for row in data["bundled_boards"]:
        board = federation(row["level"])
        kept = np.flatnonzero(np.arange(board.n) % 12 != 4)
        graph = board.adjacency()[kept][:, kept]
        d = shortest_path(graph, directed=False, indices=0, unweighted=True)[kept % 12 == 0]
        mask = d <= row["radius"]
        assert int(d.max()) == row["farthest_station"]
        assert int(mask.sum()) == row["stations_in_closed_ball"]
        clock = row["boundary_exit_clock"]
        assert clock["cut_seams"] * Fraction(clock["mean_steps_exact"]) == 53 * int(mask.sum())
        samples = clock["trials"]
        mean = clock["sum_steps"] / samples
        se = np.sqrt((clock["sum_squared_steps"] / samples - mean**2) / (samples - 1))
        assert abs(mean - float(Fraction(clock["mean_steps_exact"]))) < 5 * se
        g = np.arange(board.cells) // 4**row["level"]
        hit = mask[:, None] & mask[None, :] & ~np.eye(board.cells, dtype=bool)
        total = board.cells * (board.cells - 1)
        assert Fraction(int(hit.sum()), total) == Fraction(row["global_pair_probability"])
        assert Fraction(int((hit & (g[:, None] == g)).sum()), total) == Fraction(
            row["same_ancestor_pair_probability"])
