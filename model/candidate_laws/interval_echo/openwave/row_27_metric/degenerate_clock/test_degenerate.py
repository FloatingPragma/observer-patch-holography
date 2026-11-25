"""Check the two independent propagators, the exact path witness and the counted tests."""

import json
from fractions import Fraction
from math import factorial

import numpy as np
from scipy.sparse.linalg import expm_multiply

from model.executable_model import federation
from . import measure, replicate


def test_direct_sharp_read_matches_coherent_branch_sum():
    board = federation(0)
    h = measure.free.board_generator(board, 8)
    mask = np.diag([int(j in (measure.A, measure.B)) for j in range(board.n)])
    for changed in (False, True):
        g, start, decode = replicate.device(changed)
        second_gap = 32 if changed else 33
        folded = measure.bins(measure.branch(h + 32 * mask, [1])[1],
                              measure.branch(h + second_gap * mask, [1])[1])
        for index in (0, 1, 3, 5, 7):
            w = measure.WAITS[index]
            state = expm_multiply((-1j*w) * g, start.reshape(-1, 2)).reshape(start.shape)
            direct = replicate.material_counts_law(state, decode)
            np.testing.assert_allclose(direct[:, 0], folded[index], atol=2e-12, rtol=0)
            np.testing.assert_allclose(direct[:, 1], folded[index, [1, 0, 2, 3, 4]], atol=2e-12, rtol=0)


def test_rational_leading_terms_from_full_integer_recurrence():
    board = federation(0)
    a, b = measure.A, measure.B
    data = json.loads((measure.HERE/'results.json').read_text())
    neighbours = [[q for q in board.neighbours[j] if {q, j} != {a, b}] for j in range(board.n)]
    for d, start in data["launch_by_native_distance"].items():
        record = data["short_time_certificates"][d]
        target = record["terminal_excess_port"]
        order = record["probability_power"] // 2
        initial = [int(i == start) for i in range(board.n)]
        arms = [initial, initial]
        for k in range(1, order+1):
            # Integer powers of twice the two full spatial generators. All
            # walks, diagonal pauses and interactions contribute in this check.
            for i, gap in enumerate((32, 33)):
                prev = arms[i]
                arms[i] = [(17*len(ns) + 48*(6-int(board.degree[j])) + 2*gap*int(j in (a,b)))
                           *prev[j] - sum(prev[q] for q in ns) for j, ns in enumerate(neighbours)]
            delta = arms[0][target] - arms[1][target]
            if k < order:
                assert delta == 0
        coeff = Fraction(delta, 2**(order+1) * factorial(order))**2
        assert coeff == Fraction(record["leading_probability_coefficient"]) > 0


def test_sensitivity_and_complete_observer_banks():
    r = json.loads((measure.HERE/"replication.json").read_text())
    for name, entry in r["decisions"].items():
        keep = np.array(r["counts"][name])
        np.testing.assert_array_equal(keep.sum(axis=1), replicate.SHOTS)
        assert not keep[:, 4].any()
        if name.startswith("native_") and name.endswith("_mobile"):
            assert all(entry["resolved_against_blank"])
            assert np.all(keep[:, 2:4].sum(axis=1) > 0)
        else:
            assert not any(entry["opposite_order"])
        assert not any(entry["resolved_hoeffding"])
        for w, (lo, hi) in zip(replicate.WAITS, r["dimer_duration_intervals"][name], strict=True):
            assert 0 < lo <= w <= hi < 8*np.pi
    with np.load(measure.HERE/"counts.npz") as pilot:
        for sign, col in (("plus", 1), ("minus", 0)):
            assert not pilot[f"native_{sign}_free1"][3:, col].any()
            assert pilot[f"native_{sign}_free0"][5, col] > 500
            assert not pilot[f"no_sign_exchange_at_two_{sign}_free0"][:, col].any()
        for k in pilot.files:
            if not k.endswith("_timer"):
                np.testing.assert_array_equal(pilot[k].sum(1), measure.SHOTS)
