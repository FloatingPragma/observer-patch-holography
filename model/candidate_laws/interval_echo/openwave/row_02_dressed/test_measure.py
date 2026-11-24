import json
from itertools import combinations

import numpy as np
from scipy import sparse

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.experiments.interval_echo_encounter.conjugation import defect_generator
from model.candidate_laws.interval_echo.openwave.row_02_dressed import measure


def test_separated_species_identity_on_the_full_cell():
    board = measure.cell()
    holes = PairSector(board)
    states = [
        (a, b, x)
        for a, b in combinations(range(board.n), 2)
        for x in range(board.n)
        if x not in (a, b)
    ]
    ids = [holes.index(a, b) * board.n + x for a, b, x in states]
    location = {s: i for i, s in enumerate(states)}
    for k in (8, -8):
        potential = k * (
            2 * (board.degree[holes.aa] + board.degree[holes.bb])
            - 12
            - 2 * holes.adjacent
        )
        apart = sparse.kronsum(
            defect_generator(board, 1, k, 0),
            holes.exchange() + sparse.diags(potential, dtype=float),
            format="csr",
        )
        difference = (
            measure.restriction(board, states, k) - apart[ids][:, ids]
        ).toarray()
        contact = np.zeros_like(difference)
        for i, (a, b, x) in enumerate(states):
            for v, w in ((a, b), (b, a)):
                if x in board.neighbours[v]:
                    contact[i, i] += 2 * k - 0.5
                    pair = sorted((w, x))
                    contact[i, location[(*pair, v)]] -= 0.5
        np.testing.assert_allclose(difference, contact, rtol=0, atol=1e-12)


def test_rank_and_binding_record():
    j = json.loads((measure.HERE / "results.json").read_text())
    assert sum(s["occurrences"] for s in j["census"]["layouts"]) == 420
    assert all(
        s["ranks"][-1] == s["interior_states"] == 252 for s in j["census"]["layouts"]
    )
    for row in j["upper_band"]:
        assert row["whole_sector_upper_bound"] < row["upper_band_target_lower"]
        assert row["intact_probability_all_exposures"] > 0.57
    for name in ("high_release", "exchange_only_release"):
        for r in j[name]["exposures"]:
            assert sum(r["counts"]) == measure.SHOTS
            assert 0 <= r["full_board_bin_error_upper"] < 0.03
    assert sum(j["high_release"]["exposures"][-1]["probabilities"][:3]) > 0.97
    assert sum(j["exchange_only_release"]["exposures"][-1]["probabilities"][:3]) < 0.57


def test_assembly_difference_from_counts_only():
    j = json.loads((measure.HERE / "results.json").read_text())
    for d, expected in zip(j["strengths"].values(), (31, -33), strict=True):
        rows = d["transfer"]["acquisitions"]
        phases = []
        for row in rows:
            counts = np.array(row["counts_by_phase"])
            np.testing.assert_array_equal(counts.sum(axis=1), measure.SHOTS)
            coeff = (counts[:, 0] - counts[:, 1]) / measure.SHOTS
            phases.append(np.angle(complex(coeff[0] - coeff[1], coeff[3] - coeff[2])))
        assert abs(-(phases[-1] - phases[0]) / rows[-1]["time"] - expected) < 0.2
