"""Independent checks of pulse words, all archived detector bins, and the scale assignment."""

import gzip
import hashlib
import json
from math import cos, pi, sin
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from . import apparatus as app
from . import measure

RECEIPT = json.loads(gzip.decompress((measure.HERE / "records.json.gz").read_bytes()))


def test_native_calendars_and_full_grating_phases():
    points, axes = app.symbol.chart()
    layout = RECEIPT["layout"]
    board = app.symbol.carrier.federation(0)
    links = {frozenset(edge) for edge in board.seams}
    phase_port = layout["phase"]
    assert board.degree[phase_port] == 5
    assert [head["owner"] for head in layout["heads"]] == [18, 19]

    def transposition(word, first, last):
        slots = list(range(board.n))
        for a, b in word:
            assert frozenset((a, b)) in links
            slots[a], slots[b] = slots[b], slots[a]
        expected = list(range(board.n))
        expected[first], expected[last] = expected[last], expected[first]
        assert slots == expected

    for source, moves in layout["private_routes"].items():
        transposition(moves, int(source), phase_port)
    for head, routes in zip(layout["heads"], layout["observer_routes"], strict=True):
        assert all(board.degree[p] == 5 for p in head["pair"])
        assert all(
            frozenset(pair) in links
            for pair in zip(head["pair"], head["references"], strict=True)
        )
        assert all(
            p // 12 == head["owner"] for p in (*head["pair"], *head["references"])
        )
        assert [tuple(r["source"]) for r in routes] == [app.ROD, *app.EDGES]
        for r in routes:
            for name, first, last in zip(
                ("first", "second"), r["source"], head["pair"], strict=True
            ):
                transposition(r[name], first, last)

    h = app.symbol.local_block(1)
    for label, word in RECEIPT["preparations"].items():
        source, arm, key, *index = label.split("/")
        scale, gain, _ = measure.ARMS[arm]
        if source == "held":
            support, target = range(12), (8 / 3) * scale * (points @ axes[key])
        else:
            k = gain * (pi / 5, pi / 7)[int(key)]
            if source == "rod":
                support, target = app.ROD, np.array([0, k])
            else:
                u, v = app.EDGES[int(index[0])]
                support, target = (
                    range(12),
                    k * scale * (points @ (points[v] - points[u])),
                )
        state = np.zeros(board.n, complex)
        state[word["split"][0][0]] = 1
        for a, b, t in word["split"]:
            assert frozenset((a, b)) in links
            state[[a, b]] = expm(-1j * float(t) * h) @ state[[a, b]]
        for u, t in word["private_dwells"]:
            moves = layout["private_routes"][str(u)]
            for a, b in moves:
                state[a], state[b] = state[b], state[a]
            state[phase_port] *= np.exp(-1j * 4 * app.symbol.KAPPA * float(t))
            for a, b in reversed(moves):
                state[a], state[b] = state[b], state[a]
        np.testing.assert_allclose(
            state[list(support)],
            np.exp(1j * target) / np.sqrt(len(support)),
            atol=3e-9,
            rtol=0,
        )
        assert all(state[p] == 0 for p in set(range(board.n)).difference(support))


def test_every_count_and_independent_energy_predictions():
    shots = RECEIPT["shots"]
    points, axes = app.symbol.chart()
    n_banks = 0

    def check(counts, p_in, cond):
        nonlocal n_banks
        assert len(counts) == 3 and sum(counts) == shots and min(counts) >= 0
        # Frozen one-exit bound also applies to the in-pair indicator.
        for k, p in ((counts[1], p_in * cond), (sum(counts[:2]), p_in)):
            lo, hi = measure.interval(k, shots, False)
            assert lo <= p <= hi
        lo, hi = measure.interval(counts[1], sum(counts[:2]), True)
        assert lo <= cond <= hi
        n_banks += 1

    for arm, rows in RECEIPT["calibration"].items():
        scale, gain, mode = measure.ARMS[arm]
        for k, row in zip((pi / 5, pi / 7), rows, strict=True):
            native = mode == "native"
            for rod in row["reference"]:
                check(rod, 1, sin(gain * k / 2) ** 2 if native else 0.5)
            for head in row["edges"]:
                assert len(head) == 30
                for counts in head:
                    check(
                        counts, 1 / 6, sin(gain * k * scale / 2) ** 2 if native else 0.5
                    )

    for arm, directions in RECEIPT["held"].items():
        q = (8 / 3) * measure.ARMS[arm][0]
        for axis, banks in directions.items():
            phases = q * (points @ axes[axis])
            edge_probs = [(1 - cos(phases[a] - phases[b])) / 2 for a, b in app.EDGES]
            for head in banks:
                for c, p in zip(head, edge_probs, strict=True):
                    check(c, 1 / 6, p)
            lo, hi = RECEIPT["energies"][arm][axis]
            assert lo <= sum(edge_probs) / 6 <= hi
    assert n_banks == 1100
    assert 3 * n_banks <= RECEIPT["budget"]


def test_receiver_vetoes_and_kernel_control():
    assert measure.analysis(RECEIPT["calibration"], RECEIPT["held"]) == (
        RECEIPT["intervals"],
        RECEIPT["energies"],
    )
    for arm, row in RECEIPT["intervals"].items():
        if arm in ("nominal", "dilated", "gain"):
            a = measure.ARMS[arm][0]
            lo, hi = row["joint_intersection"]
            assert row["status"] == "calibrated" and lo < a < hi and hi - lo < 0.04
        else:
            assert row["joint_intersection"] is None
    folded = RECEIPT["intervals"]["folded"]
    assert folded["status"] == "no_common_branch"
    for (lo, hi), value in zip(folded["calibration_hulls"], (1.5, 5.5), strict=True):
        assert lo < value < hi
    assert RECEIPT["intervals"]["blind"]["status"] == "unresolved_reference"
    points, axes = app.symbol.chart()
    unit = points / np.linalg.norm(points, axis=1)[:, None]
    vertex = {
        axis: 5 / 24 * sum(1 - np.cos(4 * (unit @ n))) for axis, n in axes.items()
    }
    energies = RECEIPT["energies"]["nominal"]
    assert energies["face"][0] - energies["port"][1] > 0.08
    assert vertex["face"] - vertex["port"] < -0.38
    for path, digest in RECEIPT["source_sha256"].items():
        assert (
            hashlib.sha256(
                (Path(measure.HERE.parents[2]) / path).read_bytes()
            ).hexdigest()
            == digest
        )
