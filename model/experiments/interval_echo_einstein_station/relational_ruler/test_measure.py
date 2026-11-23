import json
from types import SimpleNamespace

import numpy as np
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import expm_multiply

from model.executable_model import federation
from model.experiments.material_inertia import neutral_defect as defect
from . import measure


def test_clock_census_matches_independent_full_metric():
    board = federation(0)
    routes = measure.Routes(SimpleNamespace(external=np.array(board.external)))
    distance = shortest_path(board.adjacency(), directed=False, unweighted=True)
    home = routes.range_ball(0, 10)
    census = measure.census(routes, home, 3)
    anchors = distance[0] <= 3
    obs = np.sum((distance[anchors, :] > 0) & (distance[anchors, :] <= 3))
    assert obs == census["ordered_pairs"]
    assert sum(anchors) * 3 == census["labels"]


def test_local_flight_matches_full_federation_without_cutoff():
    board = federation(0)
    routes = measure.Routes(SimpleNamespace(external=np.array(board.external)))
    x = 0
    home = routes.range_ball(0, 20)
    y = next(port for port, d in home.items() if d == 4)
    xsites, ysites, _, _, h, _ = defect.generators(board)
    initial = np.zeros(len(xsites))
    initial[defect.ordered_index(x, y, board.n)] = 1
    full = expm_multiply(
        -1j * measure.T * h, initial, traceA=-1j * measure.T * h.diagonal().sum()
    )
    truncated, allowance = measure.sharp_evolution(routes, x, y)
    terms = np.zeros(len(xsites), complex)
    for (a, b), z in truncated.items():
        terms[defect.ordered_index(a, b, board.n)] = z
    assert np.linalg.norm(terms - full * np.exp(14j * measure.T)) < allowance + 1e-12


def test_held_bearings_and_shell_leakage():
    packets = json.loads((measure.HERE / "results.json").read_text())
    shell_results = json.loads((measure.HERE / "shell.json").read_text())
    for p, r in zip(packets["panels"][1:], shell_results["records"][1:], strict=True):
        low, high = p["bearings"]
        assert low["degree"] == high["degree"]
        assert low["vacancy"] != high["vacancy"]
        assert (
            low["closer_interval_99"][1] + low["probability_error_bound"]
            < high["closer_interval_99"][0] - high["probability_error_bound"]
        )
        assert sum(low["counts"]) == sum(high["counts"]) == measure.SHOTS
        for block in r["blocks"]:
            vals = [
                int(k)
                for k, v in block["inward_multiplicities"].items()
                for i in range(v)
            ]
            assert (
                abs(
                    block["fraction_orthogonal"]
                    - np.var(vals) / np.mean(np.asarray(vals) ** 2)
                )
                < 2e-15
            )
    fine = packets["panels"][-1]
    assert abs(fine["dimensions"]["labels"] - 3) < 0.001
    assert abs(fine["dimensions"]["ordered_pairs"] - 4) < 0.02
