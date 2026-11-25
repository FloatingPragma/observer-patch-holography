import json
from itertools import product

import numpy as np
from scipy.linalg import expm
from scipy.sparse import coo_matrix, diags, eye, kron

from . import exchange as ex


def test_true_two_body_generator_and_exchange_sign():
    for closed in (False, True):
        one, native, fermion = ex.matrices(closed)
        # Independent labeled tensor model with the measured diagonal correction.
        sites = list(product(range(ex.SIZE), repeat=2))
        correction = [-17 if a != b and one[a, b] else 0 for a, b in sites]
        lifted = kron(one, eye(ex.SIZE)) + kron(eye(ex.SIZE), one) + diags(correction, dtype=float)
        ij = np.array([ex.SIZE * i + j for i, j in ex.PAIRS])
        ji = np.array([ex.SIZE * j + i for i, j in ex.PAIRS])
        for sign, target in [(1, native), (-1, fermion)]:
            embedding = coo_matrix(
                (
                    np.r_[np.ones(len(ij)), np.full(len(ij), sign)] / np.sqrt(2),
                    (np.r_[ij, ji], np.tile(np.arange(len(ij)), 2)),
                ),
                shape=(ex.SIZE**2, len(ij)),
            ).tocsr()
            np.testing.assert_allclose(
                (embedding.T @ lifted @ embedding).toarray(), target, atol=1e-12
            )
        if not closed:
            np.testing.assert_array_equal(native, fermion)
        else:
            crossing = abs(native - fermion) > 1e-10
            assert crossing.sum() == 2 * (ex.SIZE - 2)
            for a, b in np.argwhere(crossing):
                assert set(ex.PAIRS[a]).symmetric_difference(ex.PAIRS[b]) == {
                    0,
                    ex.SIZE - 1,
                }


def test_holds_and_independent_large_time_propagation():
    data = json.loads((ex.HERE / "exchange.json").read_text())
    one, plus, minus = ex.matrices(True)
    src, _ = ex.prep_two(ex.source(1), ex.source(1, "right"))
    for h, expected in [
        (plus, data["cases"][1]["native_two_a"]),
        (minus, data["cases"][1]["fermion_same_interactions_two_a"]),
    ]:
        # Start at a short true native interval and square the complete propagator.
        shift = h[0, 0]
        one_step = expm(-1j * (h - shift * np.eye(len(h))) / 128)
        for _ in range(34):
            one_step = one_step @ one_step
        p = abs(one_step @ src) ** 2
        assert abs(sum(p) - 1) < 4e-5
        assert abs(sum(p[ex.SECOND < ex.N]) - expected) < 4e-5
    for case in data["cases"]:
        for name in (
            "open_path_joint",
            "independent_rail_joint",
            "next_label_joint",
            "two_loop_by_port_pair",
            "one_loop_by_port",
        ):
            assert sum(case[name]) == ex.SHOTS
        lo, hi = case["two_a_interval"]
        assert lo > 0.9996 - 0.0001 and lo < case["native_two_a"] < hi
        assert case["fermion_same_interactions_two_a"] < 0.48
        lo, hi = case["one_a_interval"]
        assert lo < case["single_native_a"] < hi < 0.70 < case["single_pi_flux_a"]
        assert abs(case["calibrated_single_transfer"] - 0.5) < 0.012


def test_physical_braid_and_reversed_phase_control():
    d = json.loads((ex.HERE / "exchange.json").read_text())
    phase = d["native_phase"]
    ring = ex.hardware()["route"]
    board = ex.paired.federation(0)
    for name, expected in [("one_around_loop", [0]), ("two_around_each_other", [1, 0])]:
        p = phase[name]
        index = dict(zip(p["start"], range(len(p["start"]))))
        wraps = 0
        for a, b in p["hops"]:
            assert (
                ring[b] in board.neighbours[ring[a]] and a in index and b not in index
            )
            wraps += int(a == len(ring) - 1 and b == 0)
            index[b] = index.pop(a)
        assert wraps == 1 and [index[x] for x in p["start"]] == expected
        assert (
            p["xy_bright"][0] == ex.SHOTS
            and abs(p["xy_bright"][1] / ex.SHOTS - 0.5) < 0.001
        )
    assert all(
        r["real_imag_and_pi_shifted_bright"][0] == ex.SHOTS
        and r["real_imag_and_pi_shifted_bright"][2] == 0
        for r in phase["same_orbital_records"]
    )
