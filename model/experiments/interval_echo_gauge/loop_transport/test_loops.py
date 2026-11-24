"""Independent native loop, sector-ground, and observation-denominator checks."""

import hashlib
import json

import numpy as np
from scipy.linalg import expm
from scipy.stats import fisher_exact

from model.executable_model import cell, federation
from model.experiments.material_inertia.neutral_defect import generators, ordered_index
from .measure import HERE, faces, paths


def test_vacancy_blocks_and_dense_history():
    board = cell()
    _, vacancy, _, _, h, _ = generators(board)
    labels = faces(board)
    root, a, b, c = labels
    routes = paths(h, board, labels)
    record = json.loads((HERE / "results.json").read_bytes())
    first = record["boards"][0]
    full = expm(-0.5j * h.toarray())
    for row in first["history"]:
        if row["hold"] != 0.5:
            continue
        state = np.zeros(h.shape[0], complex)
        state[ordered_index(b, root, board.n)] = 1
        for next_root in routes[row["order"]]["route"]:
            state = full @ state
            state[vacancy != next_root] = 0
        assert abs((abs(state) ** 2).sum() / row["success"] - 1) < 1e-10
        for j, port in enumerate((a, b, c)):
            np.testing.assert_allclose(
                abs(state[ordered_index(port, root, board.n)]) ** 2,
                row["probabilities"][6 + j],
                rtol=1e-9,
                atol=1e-20,
            )
    assert routes["UV"]["excess_from_abc"] == [c, a, b]
    assert routes["VU"]["excess_from_abc"] == [b, c, a]


def test_records_and_source_ground():
    d = json.loads((HERE / "results.json").read_bytes())
    assert (
        d["protocol_sha256"]
        == hashlib.sha256((HERE / "DESIGN.md").read_bytes()).hexdigest()
    )
    for src in d["boards"]:
        assert (
            src["single_vacancy_commutator"] == 1
            and src["marked_plane_half_commutator_trace"] == -0.5
        )
        for read in src["history"]:
            counts = np.array(read["counts"])
            p = np.array(read["probabilities"])
            n = d["history_trials"]
            assert len(counts) == 10 and counts.sum() == n
            # Stable binomial scale for nearly certain or nearly impossible records.
            assert np.all(
                abs(counts - n * p) < 8 * np.sqrt(n * p * np.maximum(1 - p, 0)) + 12
            )
        v = src["vacuum"]
        lo, hi = v["ground_bracket"]
        assert hi - lo < 1e-8 and lo - 1e-12 <= src["ground_energy"] <= hi + 1e-12
        assert v["ground_residual"] < 1e-9
        assert v["correlation"][0] <= 1 / src["ports"]
        assert np.all(np.diff(v["correlation"]) <= 0)
        if src["ports"] == 12:
            assert v["mean"] < 1e-27
    s = json.loads((HERE / "placement.json").read_bytes())
    assert (
        s["protocol_sha256"]
        == hashlib.sha256((HERE / "PLACEMENT.md").read_bytes()).hexdigest()
    )
    assert (s["rooted_triangles"], s["eligible_pairs"]) == (1200, 13500)
    assert len(set(s["home"] + s["remote"])) == 6
    for r in s["rows"]:
        a, b = r["counts"]
        p = fisher_exact([[a, s["shots"] - a], [b, s["shots"] - b]]).pvalue
        assert p == r["two_sided_p"] and (p < 0.01 / 6) == r["resolved"]
    assert s["rows"][0]["home_connected"] > d["boards"][-1]["vacuum"]["correlation"][0]


def test_shared_spokes_on_every_federation():
    for level in (0, 1, 2):
        board = federation(level)
        # Count unordered pairs of faces sharing a rooted local spoke.
        counts = np.zeros(board.n, int)
        incidence = {}
        for face in board.faces():
            for a in face:
                for b in face:
                    if a != b:
                        incidence[a, b] = incidence.get((a, b), 0) + 1
        for (root, _), number in incidence.items():
            counts[root] += number * (number - 1) // 2
        np.testing.assert_array_equal(counts, np.full(board.n, 5))
