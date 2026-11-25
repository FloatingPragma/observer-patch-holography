from collections import defaultdict
from itertools import product
import json

import numpy as np
from scipy.linalg import expm

from . import measure as run
from . import spectrum as spec


def times(x, y):
    a, b = x
    c, d = y
    return a * c + 5 * b * d, a * d + b * c


def add(x, y):
    return x[0] + y[0], x[1] + y[1]


def test_full_native_three_defect_operator_on_touching_owners():
    board = run.federation(2)
    charts = spec.compact.gate_profiles(board)
    ends = next(
        (a, b)
        for a, b in board.external
        if 1 in spec.modes(charts[a // 12]) and 1 in spec.modes(charts[b // 12])
    )
    owners = [p // 12 for p in ends]
    neighbour = next(
        y // 12
        for a, b in board.external
        for x, y in ((a, b), (b, a))
        if x // 12 in owners and y // 12 not in owners
    )
    owners.append(neighbour)
    vectors = [
        spec.modes(charts[i])[tag] for i, tag in zip(owners, [1, 0, 0], strict=True)
    ]
    terms = {}
    for i, j, k in product(range(12), repeat=3):
        coeff = (1, 0)
        for pair, x in zip(vectors, (i, j, k), strict=True):
            coeff = times(coeff, tuple(int(row[x]) for row in pair))
        if coeff != (0, 0):
            a, b = sorted([12 * owners[0] + i, 12 * owners[1] + j])
            terms[a, b, 12 * owners[2] + k] = coeff

    # Two vacancies, one excess. Acting with twice the full native operator
    # on each occupied word also checks every zero row reached by one swap.
    def applied(psi):
        hv = defaultdict(lambda: (0, 0))
        for key, amp in psi.items():
            a, b, c = key
            deviations = {a: -1, b: -1, c: 1}
            edges = {tuple(sorted((x, y))) for x in key for y in board.neighbours[x]}
            pot = sum(
                (6 - int(board.degree[x])) * (s * s + 2 * s)
                for x, s in deviations.items()
            )
            unequal = 0
            for x, y in edges:
                u, v = deviations.get(x, 0), deviations.get(y, 0)
                if u == v:
                    continue
                unequal += 1
                pot += (u - v) ** 2

                def sw(i):
                    return y if i == x else x if i == y else i

                aa, bb = sorted([sw(a), sw(b)])
                dest = (aa, bb, sw(c))
                hv[dest] = add(hv[dest], times((-1, 0), amp))
            hv[key] = add(hv[key], times((16 * pot + unequal, 0), amp))
        return {k: v for k, v in hv.items() if v != (0, 0)}

    # Twice the eigenenergy is 2*(35+35+67)−1+sqrt(5).
    assert applied(terms) == {k: times((273, 1), v) for k, v in terms.items()}
    changed = terms.copy()
    key = next(iter(changed))
    changed[key] = times((-1, 0), changed[key])
    assert applied(changed) != {k: times((273, 1), v) for k, v in changed.items()}


def test_sorted_support_storage_spectroscopy_and_free_phase_impulse():
    stored = json.loads((spec.HERE / "results.json").read_text())
    fresh = run.run()
    assert [r["distinct_energies"] for r in fresh["levels"]] == [10, 10478, 698858]
    light = fresh["levels"][-1]["comparisons"]["ordinary_compact"]["measured_muon"]
    assert light["twice_energy"] == [14514, -18]
    obs = fresh["complete_storage_scan"]
    lo, hi = obs["simultaneous_interval_99"]
    assert lo < light["energy"] < hi < float(fresh["targets"]["measured_muon"]) * 35
    assert hi - lo < 0.002
    assert (
        abs(
            obs["formation_in_dimer_units"]
            - stored["complete_storage_scan"]["formation_in_dimer_units"]
        )
        < 0.0004
    )
    assert sum(v["owners"] for v in obs["owner_groups"].values()) == 143
    for row in obs["owner_groups"].values():
        assert np.all(
            np.array(row["pooled_counts"]).sum(axis=2)
            == obs["shots_per_owner"] * row["owners"]
        )
    for row in fresh["spectroscopy"]["0"]["bright_dark_else_by_wait"].values():
        assert row[0] == [run.SHOTS, 0, 0] and row[2] == [0, run.SHOTS, 0]
    for tag in (-1, 1):
        row = fresh["spectroscopy"][str(tag)]
        assert (
            abs(row["signed_gap_from_time_one"] - row["predicted_difference"]) < 0.012
        )
    law = run.store.LAW
    swap = np.eye(9)[:, [law.SRC.index((b, a)) for a, b in law.SRC]]
    np.testing.assert_allclose(
        expm(-1j * np.pi * run.store.anchor.free_seam(8)), swap, atol=5e-13
    )
