import json
from math import factorial

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from .measure import HERE, SHOTS, band, events, fixture, material


def test_projection_of_physical_seams_and_dense_flight():
    law = candidate.native_law
    board = federation(0)
    dimers, contacts, links = fixture(board, 2)
    h, _ = material(links, 8)
    exchange, contrast = law.local_term_parts()
    physical = exchange + 8 * contrast
    words = law.SRC
    basis = np.eye(law.DIM)[:, : law.N_SRC]
    dimer_cols = basis[:, [words.index((1, 0)), words.index((0, 1))]]
    local = dimer_cols.T @ physical @ dimer_cols
    np.testing.assert_allclose(physical @ dimer_cols, dimer_cols @ local, atol=1e-12)
    donor = (basis[:, words.index((2, 0))] + basis[:, words.index((0, 2))]) / np.sqrt(2)
    collar = (law.echo() @ donor).reshape(-1, len(words))[:, words.index((1, 1))]
    bits = np.arange(2 ** len(dimers))[:, None] >> np.arange(len(dimers)) & 1
    dense = np.zeros(h.shape)
    for i, (a, b) in enumerate(dimers):
        indices = np.arange(len(dense))
        dense[indices, indices ^ (1 << i)] += local[0, 1]
        dense[indices, indices] += local[bits[:, i], bits[:, i]]
        assert board.degree[a] == board.degree[b]
    for u, v in contacts:
        i, j = [next(t for t, pair in enumerate(dimers) if p in pair) for p in (u, v)]
        loads = np.stack(
            (bits[:, i] == dimers[i].index(u), bits[:, j] == dimers[j].index(v)), axis=1
        ).astype(int)
        values = {}
        for pair in set(map(tuple, loads)):
            col = np.kron(collar, np.eye(law.N_SRC)[:, words.index(pair)])
            value = col @ physical @ col
            np.testing.assert_allclose(physical @ col, value * col, atol=1e-12)
            values[pair] = value
        dense += np.diag([values[tuple(p)] for p in loads])
    shift = dense[0, 0] - h[0, 0]
    np.testing.assert_allclose(
        dense - h.toarray(), shift * np.eye(len(dense)), atol=2e-12
    )
    psi = np.zeros(h.shape[0])
    psi[11] = 1.0
    left = expm(-0.25j * dense) @ psi
    right = expm(-0.25j * h.toarray()) @ psi
    np.testing.assert_allclose(left, np.exp(-0.25j * shift) * right, atol=2e-12)


def test_native_short_time_orders():
    for size in (2, 3):
        _, _, links = fixture(federation(0), size)
        for kappa in (8, 0):
            h, walls = material(links, kappa)
            integer = (2 * h).astype(np.int64)
            event = events(size)
            for launch in np.flatnonzero(walls == event["input"]):
                column = np.eye(len(walls), dtype=np.int64)[launch]
                for order in range(1, size + 2):
                    column = integer @ column
                    for name, first in {
                        "along_A": 1,
                        "split_at_center": 1,
                        "turn_B": size + 1,
                        "turn_C": size + 1,
                        "split_in_B_and_C": 3,
                    }.items():
                        entries = column[walls == event[name]]
                        if order < first:
                            assert not entries.any()
                        if order == first:
                            assert sorted(abs(entries)) == [0, factorial(first)]


def test_archive_keeps_every_read_and_changed_law():
    summary = json.loads((HERE / "results.json").read_text())
    with np.load(HERE / "counts.npz") as bank:
        for key, row in summary["runs"].items():
            ref = bank[key + "_blank"]
            data = bank[key + "_incident"]
            for setting in (ref, data):
                assert np.all(setting.sum(axis=2) == SHOTS)
            if "_k0_" in key:
                lo, hi = row["remote_both"]["difference_band"]
                assert lo <= 0 <= hi
    for size in (2, 3):
        for reverse in (0, 1):
            result = summary["runs"][f"b20_L{size}_k8_r{reverse}"]
            assert result["remote_both"]["difference_band"][0] > 0
            assert (
                min(abs(x) for x in result["center_coherence"]["difference_band"])
                > 0.02
            )


def test_scalar_turn_and_straight_predictions():
    s = json.loads((HERE / "scalar.json").read_text())
    for arm in s["settings"].values():
        assert sum(arm["counts"]) == SHOTS
        assert arm["reads"]["along_A"]["compatible"]
        for target in ("turn_B", "turn_C"):
            read = arm["reads"][target]
            assert not read["compatible"]
            assert read["count"] <= 1
            assert read["forecast"][0] > 0.000615


def test_rare_count_intervals_invert_both_binomial_tails():
    from scipy.stats import binom
    from .measure import TAIL

    for k in (1, 1000, 1041, SHOTS // 2, SHOTS - 1):
        lo, hi = band(k)
        assert lo < k / SHOTS < hi
        upper, lower = binom.sf(k - 1, SHOTS, lo), binom.cdf(k, SHOTS, hi)
        assert max(upper, lower) <= TAIL * 1.0001
        if k < SHOTS - 1:
            np.testing.assert_allclose([upper, lower], [TAIL, TAIL], rtol=0.002, atol=0)
