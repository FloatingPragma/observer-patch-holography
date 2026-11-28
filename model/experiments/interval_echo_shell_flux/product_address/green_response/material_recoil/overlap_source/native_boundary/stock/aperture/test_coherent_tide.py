from collections import defaultdict
from fractions import Fraction as Q
from itertools import product
from types import SimpleNamespace
import json

import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    coherent_tide as co,
)


def test_time_zero_curvature_from_the_full_three_station_hamiltonian():
    board = SimpleNamespace(seams=[(0, 1), (1, 2), (0, 2)])
    ends = (0, 2)
    h = co.lab.energy.pair.stock.kernels()[2]
    # The two blank contacts have no off-block matrix elements. Embed their
    # full amplitudes and the written head before squaring the total matrix.
    assert np.all(h[9:, :9] == 0)
    words = list(product(range(3), repeat=3))
    basis = [(w, m) for w in words for m in range(11)]
    lookup = {k: i for i, k in enumerate(basis)}
    total = np.zeros((len(basis), len(basis)), int)
    for col, (w, m) in enumerate(basis):
        for i, j in board.seams:
            old = 9 * (m if (i, j) == ends else 0) + 3 * w[i] + w[j]
            for row in np.flatnonzero(h[:, old]):
                x = list(w)
                x[i], x[j] = divmod(int(row) % 9, 3)
                dest = (tuple(x), int(row) // 9 if (i, j) == ends else m)
                total[lookup[dest], col] += h[row, old]
    np.testing.assert_array_equal(total, total.T)
    hop = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
    base, use = co.tidal.chart(hop, ends, Q(3, 2))
    d = base + Q(3, 2) * use
    D = d[0, 2]
    axis = (d[0] ** 2 - d[2] ** 2) / (2 * D**2)
    labels = (2, 9)
    dist = {((1, 1, 1), (labels,)): Q(1)}
    rates, _ = co.spectral_contact()
    for t in (1, 2):
        phi = np.array(
            [sum((axis[i] for i in range(3) if w[i] == t), Q(0)) for w, m in basis],
            object,
        )
        pair = np.array(
            [
                sum(
                    (
                        d[i, j] / D
                        for i in range(3)
                        for j in range(i)
                        if w[i] == w[j] == t
                    ),
                    Q(0),
                )
                for w, m in basis
            ],
            object,
        )
        results = []
        for f in (phi, pair):
            launches = []
            for a, b in [(0, 0), (t, 0), (0, t), (t, t)]:
                value = Q(0)
                for y, m in product(range(3), labels):
                    col = lookup[((a, y, b), m)]
                    value += ((f - f[col]) @ total[:, col] ** 2) * Q(
                        1, 2 * 3 * len(labels)
                    )
                launches.append(value)
            results.append(launches)
        v, p = results
        out = co.observe(hop, ends, Q(3, 2), co.closure(dist, t, board, [ends], rates))
        assert v[2] - v[1] == Q(out["normalized_chord_acceleration"])
        assert p[3] - p[2] - p[1] + p[0] == Q(out["normalized_pair_acceleration"])


def test_exterior_witness_on_full_board_native_sharp_columns():
    from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
        exterior_tide as ex,
    )

    _, words, *_ = co.lab.energy.pair.setup()
    board, _, marks, *_ = co.lab.energy.returned(words)
    adj, hops, _ = co.flight.records.thermal_radar.chart(board, marks)
    receipt = json.loads((co.HERE / "exterior_tide.json").read_bytes())
    U, V = receipt["held_four_hop"]["stations"]
    h = co.lab.energy.pair.stock.kernels()[2]
    rng = np.random.default_rng(581211)
    for arm in ("source", "original"):
        cost = Q(receipt["chart_costs"][arm])
        base, used = co.tidal.chart(hops, marks[-1], cost)
        d = base + cost * used
        length = d[U, V]
        phi = (d[U] ** 2 - d[V] ** 2) / (2 * length**2)
        _, chord, joint = ex.fields(hops, marks[-1], cost, adj, np.array([[U, V]]))
        old = np.zeros(board.n, dtype=int)
        old[rng.choice(board.n, 5, replace=False)] = [1, 1, 1, 2, 2]
        for t in (1, 2):
            one = []
            pair = []
            for a, b in ((0, 0), (t, 0), (0, t), (t, t)):
                loads = old.copy()
                loads[[U, V]] = a, b
                amplitudes = defaultdict(int)
                # Assemble each nontrivial material transition of the full
                # Hamiltonian before squaring. Pure-head changes have the
                # identical material read for every outcome and cancel.
                for i, j in board.seams:
                    m = (3, 8, 2)[marks.index((i, j))] if (i, j) in marks else 0
                    col = 9 * m + 3 * loads[i] + loads[j]
                    for row in np.flatnonzero(h[:, col]):
                        if row % 9 == col % 9:
                            continue
                        x = loads.copy()
                        x[i], x[j] = divmod(int(row) % 9, 3)
                        head_change = (i, j, int(row) // 9) if row // 9 != m else ()
                        amplitudes[x.tobytes(), head_change] += int(h[row, col])

                def scores(x):
                    on = np.flatnonzero(x == t)
                    return sum(phi[on], Q(0)), sum(
                        (d[i, j] / length for i in on for j in on if i < j), Q(0)
                    )

                one0, pair0 = scores(loads)
                result = [Q(0), Q(0)]
                for (blob, _), value in amplitudes.items():
                    s, p = scores(np.frombuffer(blob, dtype=int))
                    result[0] += Q(value**2, 2) * (s - one0)
                    result[1] += Q(value**2, 2) * (p - pair0)
                one.append(result[0])
                pair.append(result[1])
            assert one[2] - one[1] == chord[0]
            assert pair[3] - pair[2] - pair[1] + pair[0] == joint[0]
