"""Matched-filling single-charge spectra and position readout on actual ports."""

import json
from itertools import product
from pathlib import Path

import numpy as np

from model.executable_model import federation
from model.experiments.interval_echo_angular.charged_orientation.live_trimer.measure import (
    LAW,
    blocks,
)
from ..conjugate_flight.measure import edge, seam_data
from ..marked_bands import material


def local_sector(sign, head):
    """Subtract the head on two filled ports before selecting a coin basis."""
    matrix = blocks()[head]
    filled = LAW.SRC.index((1, 1))
    host = matrix[filled::9, filled::9]
    relative = matrix - np.kron(host, np.eye(9, dtype=int))
    words = [LAW.SRC.index(w) for w in [(1 + sign, 1), (1, 1 + sign)]]
    take = [9 * coin + i for coin in range(len(host)) for i in words]
    assert not relative[np.ix_(np.setdiff1d(range(len(matrix)), take), take)].any()
    return relative[np.ix_(take, take)] / 8


def running_coefficients(sign):
    """Native edge costs, including the physical cut, and ordinary private cost."""
    blank = local_sector(sign, "blank")
    dark = local_sector(sign, (0, 1))
    assert not dark.any()
    _, private = seam_data(8)
    cost = 8 * ((1 + sign) ** 2 - 1)
    np.testing.assert_allclose(
        8 * (private[3 * (1 + sign)] - private[3]), cost, atol=1e-13
    )
    return blank, cost


def pendant_block(sign, parity):
    matrix = local_sector(sign, (1, parity))
    active = np.flatnonzero(np.diag(matrix))
    rest = np.setdiff1d(range(4), active)
    assert len(active) == 2 and active[0] % 2 != active[1] % 2
    assert not matrix[rest].any()
    term = matrix[np.ix_(active, active)]
    np.testing.assert_array_equal(term, -np.ones((2, 2)) * (8 + (parity == -1)) / 2)
    return term


def periodic_check():
    ports, marks = material.geometry()
    board = federation(2)
    assert all(board.degree[q] == 5 for _, q in marks)
    assert len(ports) == 3 * len(marks) and len(marks) == 48
    report = {}
    for sign in [-1, 1]:
        blank, private = running_coefficients(sign)
        values = []
        for parity, hub_private, end_private in product([-1, 1], [0, 1], [0, 1]):
            side = pendant_block(sign, parity)
            for wave in range(len(marks)):
                phase = np.exp(2j * np.pi * wave / len(marks))
                h = np.diag(
                    [0, hub_private * private, 0, end_private * private]
                ).astype(complex)
                for a, b, twist in [(0, 1, 1), (1, 2, 1), (2, 0, phase)]:
                    h[a, a] += blank[0, 0]
                    h[b, b] += blank[1, 1]
                    h[a, b] += twist * blank[0, 1]
                    h[b, a] = h[a, b].conjugate()
                h[np.ix_([1, 3], [1, 3])] += side
                values.extend(np.linalg.eigvalsh(h))
        report[str(sign)] = dict(
            native_levels=len(values),
            minimum_positive=min(v for v in values if v > 0),
            largest_negative=max(v for v in values if v < 0),
        )
    return report


def bench(sign):
    board = federation(0)
    a, b = board.external[0]

    def tip(gateway):
        return min(p for p in board.neighbours[gateway] if board.degree[p] == 5)

    ports = [tip(a), a, b, tip(b)]
    assert ports[0] // 12 == a // 12 != b // 12 == ports[3] // 12
    assert list(board.degree[ports]) == [5, 6, 6, 5]
    domain = set(ports)
    store = next(
        (i, j)
        for i, j in board.seams
        if all(
            p not in domain and not domain.intersection(board.neighbours[p])
            for p in (i, j)
        )
    )
    live = {edge(i, j) for i, j in zip(ports, ports[1:])} | {edge(*store)}
    occupied = {ports[0]: 1 + sign, store[0]: 1 - sign}
    writes = []
    for i, j in board.seams:
        if not (domain | set(store)).intersection((i, j)) or edge(i, j) in live:
            continue
        # Temporarily move each charged endpoint away on its own live arm.
        swaps = []
        for start, relay in [(ports[0], ports[1]), store]:
            if start in (i, j):
                assert relay not in (i, j)
                occupied[relay] = occupied.pop(start)
                swaps.append((start, relay))
        assert occupied.get(i, 1) == occupied.get(j, 1) == 1
        writes.append(edge(i, j))
        for start, relay in reversed(swaps):
            occupied[start] = occupied.pop(relay)
    assert occupied == {ports[0]: 1 + sign, store[0]: 1 - sign}
    blank, private = running_coefficients(sign)
    matrices = {}
    for connected in [True, False]:
        h = np.diag([private, 0, 0, private]).astype(float)
        for i in range(3):
            if connected or i != 1:
                ix = [i, i + 1]
                h[np.ix_(ix, ix)] += blank
        matrices[connected] = h
    return list(map(int, ports)), list(map(int, store)), len(writes), matrices


def run():
    ports, store, writes, vac = bench(-1)
    eig = np.linalg.eigvalsh(vac[True])
    time = np.pi / (eig[1] - eig[0])
    rng = np.random.default_rng(193874)
    banks = []
    for sign, connected in product([-1, 1], [True, False]):
        h = bench(sign)[3][connected]
        e, v = np.linalg.eigh(h)
        for fraction in [0, 0.5, 1]:
            state = v @ (np.exp(-1j * e * time * fraction) * v[0].conj())
            prob = np.abs(state) ** 2
            count = rng.multinomial(2**20, prob / prob.sum())
            banks.append(
                dict(
                    sign=sign,
                    connected=connected,
                    fraction=fraction,
                    counts=count.tolist(),
                    right_probability=float(prob[-1]),
                    relative_levels=e.tolist(),
                )
            )
    return dict(
        ports=ports,
        store=store,
        exterior_writes=writes,
        exposure=time,
        starts_per_bank=2**20,
        seed=193874,
        banks=banks,
        active_pendants=periodic_check(),
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    receipt = json.dumps(run(), indent=2) + "\n"
    if parser.parse_args().write:
        Path(__file__).with_name("results.json").write_text(receipt)
    else:
        print(receipt, end="")
