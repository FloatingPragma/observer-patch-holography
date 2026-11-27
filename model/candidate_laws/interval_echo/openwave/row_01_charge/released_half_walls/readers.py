"""Routed sharp word receivers and local-potential interferometers."""

from functools import cache
from math import pi

import numpy as np
from scipy.linalg import expm

from model.experiments.interval_echo_consensus.seam_word_meter.finite_program.machine import (
    Bench,
    READ_CODES,
    template_matching,
)
from . import ring


def stations(level, n):
    board = ring.federation(level)
    ports = ring.geometry(level, n)
    other = [o for o in range(board.cells) if o not in {p // 12 for p in ports}]
    matching = template_matching(set(range(12)))
    sites = [(12 * other[0] + a, 12 * other[0] + b) for a, b in matching[:4]]
    heads = [12 * o for o in other[1:3]]
    refs = [
        tuple(sorted(i for i in board.neighbours[p] if i // 12 == p // 12)[:2])
        for p in heads
    ]
    bench = Bench(
        board=board,
        words=[sites[:1]],
        flags=sites[1:3],
        writer=sites[3],
        heads=heads,
        references=refs,
    )
    bench.start([0])
    book = []
    for p in ports:
        rows = []
        for load in range(3):
            bench.loads[p] = load
            prior = bench.loads.copy()
            rows.append([bench.observe(p, at) for at in (0, 1)])
            assert rows[-1] == [READ_CODES[load]] * 2 and bench.loads == prior
        book.append(rows)
    return dict(ports=ports, heads=heads, references=refs, native_port_codes=book)


def phase(n, theta):
    """Two positive clock holds on the actual native contrast columns."""
    _, contrast = ring.columns()
    law = ring.candidate.native_law
    clocks = np.mod([-theta / 4, theta / 4], 2 * pi) / ring.K
    eigen = [
        contrast[law.SRC.index((int(n), r)), law.SRC.index((int(n), r))] for r in (0, 2)
    ]
    return np.exp(-1j * ring.K * np.dot(clocks, eigen))


def phase_word(loads, coordinates, theta):
    return np.prod([phase(n, theta * w) for n, w in zip(loads, coordinates) if w])


def signed_motion(n, sign):
    words, s = ring.basis(n)
    if sign < 0:
        words = ((1 << n) - 1) ^ words
        s = 1 - s
    bonds = s + np.roll(s, -1, axis=1) - 1
    rows = {
        int(w): tuple(np.flatnonzero(z == sign))
        for w, z in zip(words, bonds)
        if np.count_nonzero(z) == 2
    }
    approved = []
    for a in sorted(rows):
        for j in range(n):
            jj = (j + 1) % n
            b = a ^ (1 << j) ^ (1 << jj)
            if b not in rows or b <= a:
                continue
            shared = set(rows[a]) & set(rows[b])
            if len(shared) != 1:
                continue
            fixed = shared.pop()
            (x,) = set(rows[a]) - {fixed}
            (y,) = set(rows[b]) - {fixed}
            if (y - x) % n == 2:
                first, last, origin = a, b, int(x)
            elif (x - y) % n == 2:
                first, last, origin = b, a, int(y)
            else:
                continue
            if (
                min((fixed - x) % n, (x - fixed) % n, (fixed - y) % n, (y - fixed) % n)
                < 3
            ):
                continue
            approved.append((a, b, first, last, origin, int(fixed)))
    return sorted(approved)[:2]


@cache
def exchange_plate(amount):
    law = ring.candidate.native_law
    h = law.seam_term(ring.K)
    indices = [law.SRC.index(w) for w in ((0, amount), (amount, 0))]
    gate = expm(-1j * np.pi / 2 * h)
    plate = gate[np.ix_(indices, indices)]
    np.testing.assert_allclose(plate.T.conj() @ plate, np.eye(2), atol=1e-12)
    return plate


def arm_probabilities(a, b, theta, coords):
    changed = np.flatnonzero(a != b)
    assert len(changed) == 2
    site = changed[0]
    delta = int(b[site] - a[site])
    gate = exchange_plate(abs(delta))
    z = phase_word(b, coords, theta) * phase_word(a, coords, theta).conjugate()
    # The quarter-quadrature reference is the same calibrated local service
    # on one of the two differing sites. The bright receiver is output word b.
    quarter = np.pi / (2 * delta)
    rotation = phase(b[site], quarter) * phase(a[site], quarter).conjugate()
    source = np.array([1, z]) / np.sqrt(2)
    return np.array(
        [(abs(gate @ (source * [1, offset])) ** 2)[::-1] for offset in (rotation, 1)]
    )


def interferometry(n, rng, shots):
    detector = {}
    for amount in (1, 2):
        for sign in (1, -1):
            scans = []
            for _, __, a, b, start, fixed in signed_motion(n, sign):
                coords = np.zeros(n, int)
                coords[(start + np.arange(4)) % n] = np.arange(4)
                arm = np.array(
                    [[amount * ((w >> i) & 1) for i in range(n)] for w in (a, b)]
                )
                ordinary = np.zeros((2, n), int)
                ordinary[0, (start + 1) % n] = amount
                ordinary[1, (start + 2) % n] = amount
                plots = {}
                for name, (source, angles) in {
                    "wall": (arm, ring.ANGLES),
                    "standard": (ordinary, ring.ANGLES),
                    "parked": (arm, np.zeros_like(ring.ANGLES)),
                }.items():
                    laws = np.array(
                        [arm_probabilities(*source, theta, coords) for theta in angles]
                    )
                    np.testing.assert_allclose(laws.sum(axis=-1), 1, atol=1e-12)
                    counts = rng.binomial(shots, np.clip(laws[:, :, 0], 0, 1))
                    plots[name] = {
                        "bright": counts.tolist(),
                        "bright_probability": laws[:, :, 0].tolist(),
                    }
                scans.append(
                    dict(
                        initial=a,
                        displaced=b,
                        first_bond=start,
                        final_bond=(start + 2) % n,
                        unchanged_bond=fixed,
                        coordinates=coords.tolist(),
                        ports=plots,
                    )
                )
            assert len(scans) == 2
            detector[f"{amount}_{sign}"] = scans
    return detector


def receiver_phases(n, amount, sign, probs, gate, tapered, flat, rng, shots):
    # Each commuting word effect can be read jointly with the terminal
    # voltage probes. The two calibrated quadratures have separate banks.
    tables = {}
    for selector, charge in (("tapered", tapered / 2), ("flat", flat)):
        # Integrate all words with identical physical phases and gate.
        rows = np.unique(np.column_stack((gate, charge)), axis=0)
        grouped = [
            (
                int(accept),
                left,
                right,
                probs[
                    (gate == accept) & (charge[:, 0] == left) & (charge[:, 1] == right)
                ].sum(),
            )
            for accept, left, right in rows
        ]
        fields = []
        for theta in (0, pi / 2, pi, 2 * pi):
            panels = []
            for axes in ((0, 0), (0, 1), (1, 0), (1, 1)):
                bins = np.zeros((2, 2, 2))
                for accept, left, right, w in grouped:
                    bright = [
                        (1 + (np.cos if not axis else np.sin)(theta * q)) / 2
                        for axis, q in zip(axes, (left, right))
                    ]
                    polar = [np.array([p, 1 - p]) for p in bright]
                    bins[accept] += w * np.outer(*polar)
                p = bins.ravel().clip(0, 1)
                counts = rng.multinomial(shots, p / p.sum())
                panels.append(
                    {
                        "axes": axes,
                        "probabilities": p.tolist(),
                        "counts": counts.tolist(),
                    }
                )
            fields.append(panels)
        tables[selector] = fields
    return tables
