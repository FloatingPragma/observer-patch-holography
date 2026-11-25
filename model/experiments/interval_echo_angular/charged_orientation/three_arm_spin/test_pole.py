"""Direct full-chart free evolution and integer low-order transfer checks."""

import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from . import measure as ex

DATA = json.loads(Path(__file__).with_name("results.json").read_text())


def test_native_forward_evolution_without_symmetry_reduction():
    h = ex.native_graph()
    initial = np.zeros((h.shape[0], 2))
    initial[ex.ONE] = ex.BASE
    for time in (16, 32):
        state = expm_multiply(-1j * time * h, initial)
        expected = DATA["autonomous"][str(time)]
        # Read the charge and the receiving owner's heads, keeping every
        # source configuration as a distinct discarded environment value.
        cols = np.stack([state[64 + b + 8 * np.array(ex.ONE)] for b in range(8)])
        operators = np.einsum("ro,bri->boi", ex.BASE, cols)
        p, lmin, _ = ex.pt(operators)
        assert abs(p / 9 - expected["ideal_total_yield"]) < 1e-15
        assert abs(lmin / p - expected["conditional_PT_minimum"]) < 1e-8


def test_exact_earliest_target_transfer_and_s3_order():
    # Exact integer powers on source and target differences. A coordinate
    # at the index of a zero bit-pattern denotes the all-ready reference.
    h = (2 * ex.native_graph()).tocoo()
    links = list(
        zip(h.row.tolist(), h.col.tolist(), np.rint(h.data).astype(int).tolist())
    )
    v = [0] * h.shape[0]
    v[1], v[2] = 1, -1
    moments = []
    for k in range(26):
        moments.append(v[64 + 8] - v[64 + 16])
        if k < 25:
            w = [0] * len(v)
            for row, col, n in links:
                w[row] += n * v[col]
            v = w
    assert moments == [0] * 25 + [-128]
    # Two graph transpositions are distinguishable from a spinor order sign.
    a = np.eye(3)[[1, 0, 2]]
    b = np.eye(3)[[0, 2, 1]]
    order = ex.BASE.T @ (a @ b).T @ (b @ a) @ ex.BASE
    np.testing.assert_allclose((order + order.T) / 2, -np.eye(2) / 2, atol=1e-14)


def test_energy_scope_and_original_launch_records():
    c = DATA["isolated_energy_certificate"]
    e = DATA["selected_energy"]
    assert c["interval"][0] <= e <= c["interval"][1]
    assert c["nearest_other_standard_distance"] > 3.7e-4
    assert c["ideal_witness_upper_bound"] < -2.1e-8
    assert DATA["clock"]["measured_line_weight"] > 0.999999
    assert DATA["material"]["hub"] == 4 and DATA["material"]["tip"] == 160
    assert DATA["filtered"]["original_bell_witness_99"][1] < 0
    assert DATA["autonomous"]["32"]["original_bell_witness_99"][1] < 0
    assert DATA["autonomous"]["16"]["original_bell_witness_99"][0] > 0
    lo, hi = DATA["dephased_source"]["original_bell_witness_99"]
    assert lo < 0 < hi
    for bank in (
        DATA["filtered"],
        DATA["dephased_source"],
        *DATA["autonomous"].values(),
    ):
        for row in bank["counts"]:
            assert sum(row["counts"]) == bank["starts_per_axis"]
