"""Check the exact exterior speed, arrival bounds and source-dependence scope."""

import json

import numpy as np
import pytest
from scipy import linalg, sparse

from . import certificate, measure


@pytest.mark.parametrize("charge", [-1, 1])
def test_native_speed_for_nonuniform_meter(charge):
    board = measure.formation.federation(0)
    route = measure.formation.region(board)
    cell_values = np.random.default_rng(7392).integers(-17, 24, board.cells)
    f = np.repeat(cell_values, 48)
    edges = np.asarray(board.external)
    change = cell_values[edges[:, 1] // 12] - cell_values[edges[:, 0] // 12]
    slot = int(np.argmax(abs(change)))
    x, y = edges[slot]
    state = measure.launch(
        board, [x, y], [1 / np.sqrt(2), 1j * np.sign(change[slot]) / np.sqrt(2)]
    )
    for odd in (False, True):
        h = measure.front.native.generator(board, route, charge, odd)
        fm = sparse.diags(f, dtype=int)
        current = 1j * (h @ fm - fm @ h) / 2
        expected = abs(change[slot]) / 2
        # Check the exact operator squaring on every row as well as the
        # maximizing physical input, leaving all four local coin tags active.
        square = (current @ current).tocsr()
        assert (square - sparse.diags(square.diagonal())).count_nonzero() == 0
        np.testing.assert_allclose(square.diagonal().max(), expected**2)
        np.testing.assert_allclose(current @ state, expected * state, atol=1e-14)


def test_exponential_probability_bound_independently():
    # Dense propagation on a five-cell path with disjoint exterior bridges.
    # Internal energies and mixing differ by cell. A cell-centred imaginary
    # gauge change must leave their anti-Hermitian part zero.
    rng = np.random.default_rng(73134)
    h = linalg.block_diag(*[b + b.T for b in rng.normal(size=(5, 3, 3))])
    for i in range(4):
        h[3 * i + 2, 3 * i + 3] = h[3 * i + 3, 3 * i + 2] = -0.5
    d = np.repeat(np.arange(5), 3)
    for mu in (0.1, 0.65, 1.3):
        tilted = np.exp(mu * (d[:, None] - d)) * h
        np.testing.assert_allclose(
            np.linalg.norm((tilted - tilted.T) / (2j), 2), 0.5 * np.sinh(mu)
        )
    for t in (0.4, 1.0, 2.5):
        u = linalg.expm(-1j * t * h)
        for r in range(1, 5):
            actual = np.linalg.norm(u[d >= r, :3], 2) ** 2
            assert actual < measure.probability_bound(r, t)
    assert measure.probability_bound(1, 3) == 1
    assert measure.probability_bound(1, 0) == 0


def test_full_source_certificate_and_prepared_repetition():
    saved = json.loads((measure.HERE / "certificate.json").read_text())
    assert certificate.run() == saved
    chi, tangents, *_ = measure.source_and_changes()
    for native, name in [("ordinary", "ordinary"), ("unpairing", "unpair")]:
        d = tangents[native]
        b = saved["score"][name]["absolute_log_score_upper"]
        for x, dx in zip(chi, d):
            assert np.linalg.eigvalsh(b * x + dx).min() > 0
            assert np.linalg.eigvalsh(b * x - dx).min() > 0
    # A many-copy parity count distinguishes the ordinary preparation without
    # altering the positive-effect likelihood bound on one prepared head.
    p = chi[1].trace().real
    assert 0.20 < p < 0.21
    assert 0.29 < tangents["ordinary"][1].trace().real < 0.30


def test_terminal_and_clock_records():
    receipt = json.loads((measure.HERE / "results.json").read_text())
    assert [r["cell_diameter"] for r in receipt["hardware"]] == [5, 11, 23]
    data = receipt["tested_station"]
    for q, read in data["charges"].items():
        for label, v in [("active", 0.5), ("half_bridge_hop", 0.25)]:
            row = read["velocity"][label]
            lo, hi = row["speed_interval_99pct"]
            assert lo < v < hi and hi - lo < 0.014
            for pair in row["counts"].values():
                assert sum(pair) == 2**44
        for time, shells in read["shells"].items():
            for distance, row in shells.items():
                assert (
                    0
                    < row["population"]
                    <= measure.probability_bound(int(distance), float(time))
                )
                assert abs(row["matched_log_slopes"]["ordinary"]) < 479
                assert abs(row["matched_log_slopes"]["pairing"]) < 144
