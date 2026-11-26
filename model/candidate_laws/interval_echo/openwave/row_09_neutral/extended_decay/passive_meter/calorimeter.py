"""Measure stored exchange energy at an independently prepared passive dimer."""

import argparse
import json
from collections import deque
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from model.candidate_laws.interval_echo import candidate
from ..material import edge, route
from ..measure import rebuild
from ...joined_decay.measure import core_pulse
from ...passive_decay.material import columns, commission
from ...passive_decay.measure import SHOTS, hamiltonians, physical
from ...passive_decay.readers import X, Z, charge, circuit, rotate

HERE = Path(__file__).resolve().parent
WAITS = (0, 4, 8, 16, 32, 64)
SEED = 26112609
# Two observables for each of 48 flight panels and two quadratures per clock hold.
RADIUS = float(np.sqrt(np.log(2 * 128 / 0.01) / (2 * SHOTS)))


def apparatus():
    board, old, heads, _, _ = route((3,), crossing=True)
    free = set(range(board.n)) - set(sum(old, ())) - set(sum(heads.keys(), ()))
    endpoint = old[-1]
    candidates = (
        (emitter, p, q)
        for emitter in endpoint
        for p in sorted(set(board.neighbours[emitter]) & free)
        for q in sorted(
            set(board.neighbours[p]) & set(board.neighbours[emitter]) & free
        )
        if board.degree[p] == board.degree[q] == 5
    )
    emitter, p, q = next(candidates)
    dimers = old + ((p, q),)
    other = next(k for k in endpoint if k != emitter)
    pin = min(
        k
        for k in board.neighbours[other]
        if k not in set(sum(dimers, ())) and edge(other, k) not in heads
    )
    heads = {**heads, edge(emitter, p): 33, edge(emitter, q): 32, edge(other, pin): 32}
    distance = {k: 0 for k in sum(dimers[:2], ())}
    queue = deque(distance)
    while queue:
        s = queue.popleft()
        for k in board.neighbours[s]:
            if k not in distance:
                distance[k] = distance[s] + 1
                queue.append(k)
    return board, dimers, heads, emitter, [distance[p], distance[q]]


def generators(board, dimers, heads, emitter, original):
    upstream = edge(dimers[-3][0], dimers[-2][0])
    capture = edge(emitter, dimers[-1][0])
    labels = np.arange(64)
    bits = (labels[:, None] >> np.arange(5, -1, -1)) & 1
    a, b, c, d, e, f = bits.T
    contact = -(e if emitter == dimers[-2][0] else 1 - e) * f
    kinetic = -np.kron(np.eye(32), X) / 2
    native = np.kron(original, np.eye(2)) + kinetic + np.diag(contact)
    source_terms = 32 * a + 32 * b - 32 * a * b
    comparison = native - np.diag(np.diag(native)) + np.diag(source_terms)
    settings = (
        ("native", heads, native),
        ("upstream_cut", {**heads, upstream: 32}, native + np.diag(d * e)),
        ("detector_cut", {**heads, capture: 32}, native - np.diag(contact)),
        ("no_reversal", dict.fromkeys(heads, 32), comparison),
    )
    for _, deposits, h in settings:
        real_diagonal = physical(board, dimers, deposits)
        np.testing.assert_array_equal(real_diagonal - real_diagonal[0], np.diag(h))
        expected = h - np.diag(np.diag(h))
        for i in range(6):
            expected[labels, labels ^ (1 << i)] += 0.5
        np.testing.assert_array_equal(expected, np.zeros_like(h))
    return settings


def native_meter():
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    native = exchange + candidate.KAPPA * contrast
    indices = [law.SRC.index(word) for word in ((0, 1), (1, 0))]
    embedding = np.eye(law.DIM)[:, indices]
    small = embedding.T @ native @ embedding
    np.testing.assert_allclose(native @ embedding, embedding @ small, atol=2e-13)
    return small - np.eye(2) * small[0, 0]


def gap_calibration(rng, meter):
    eigenvalues, eigenvectors = eigh(meter)
    phases = (1, 1j)
    hardware_frame = rotate(np.pi / 4) @ Z
    bras = [np.array([1, phase]) @ hardware_frame.T / np.sqrt(2) for phase in phases]
    sharp_stock = np.array([1, 0])
    initial = eigenvectors.T @ sharp_stock
    # Sharp initial orientation has equal weights on both exchange lines.
    t_last, angle_last = 0.0, 0.0
    receipt = []
    for t in 4.0 ** np.arange(7) / 64:
        evolved = eigenvectors @ (np.exp(-1j * t * eigenvalues) * initial)
        outputs = [abs(bra @ evolved) ** 2 for bra in bras]
        counts = [int(rng.binomial(SHOTS, p)) for p in outputs]
        c, s = 2 * np.array(counts) / SHOTS - 1
        angle = np.arctan2(s, c)
        if t_last:
            angle += 2 * np.pi * round((angle_last * t / t_last - angle) / (2 * np.pi))
        receipt.append(dict(wait=float(t), cos_sin_counts=counts))
        angle_last, t_last = angle, t
    return dict(
        holds=receipt,
        gap=float(angle_last / t_last),
        half_width=float(np.arcsin(2 * np.sqrt(2) * RADIUS) / t_last),
    )


def acquire():
    board, dimers, heads, emitter, distances = apparatus()
    old_heads = route((3,), crossing=True)[2]
    old = rebuild(board, dimers[:-1], old_heads, False, True)
    _, guide = eigh(old)
    core_eval, core = eigh(hamiltonians(4)[0])
    core_read = circuit(core, 32)
    pulse = np.kron(circuit(core_pulse(core), 32), np.eye(8))
    # Columns in increasing exchange energy are |+> and |->.
    meter_read = rotate(np.pi / 4) @ Z
    np.testing.assert_allclose(
        meter_read.T @ (-X / 2) @ meter_read, np.diag([-0.5, 0.5]), atol=1e-12
    )
    read = np.kron(core_read, np.kron(np.eye(8), meter_read))
    settings = generators(board, dimers, heads, emitter, old)
    rng = np.random.default_rng(SEED)
    calibration = gap_calibration(rng, native_meter())
    cold = guide[:, 0]
    stocks = [
        (label, np.kron(state, meter_read[:, 0]))
        for label, state in (("pulse", pulse @ cold), ("idle", cold))
    ]
    rows = []
    for name, deposits, h in settings:
        energies, v = eigh(h)
        for arm, psi in stocks:
            incoming = v.T @ psi
            for t in WAITS:
                after = v @ (np.exp(-1j * t * energies) * incoming)
                spectrum = abs(read.T @ after).reshape(4, 8, 2) ** 2
                p = spectrum.sum(axis=1)
                counts = rng.multinomial(SHOTS, (p / p.sum()).ravel()).reshape(4, 2)
                rows.append(
                    dict(
                        apparatus=name,
                        source=arm,
                        wait=t,
                        core_by_meter=counts.tolist(),
                    )
                )
    return dict(
        seed=SEED,
        shots_per_panel=SHOTS,
        familywise_99_half_width=RADIUS,
        boards_cells=int(board.cells),
        dimers=dimers,
        source_distances=distances,
        mount=commission(board, dimers, heads),
        native_and_changed_column_errors=[
            round(columns(), 13),
            round(columns(True), 13),
        ],
        retained_lower_core_levels=[1, 2],
        core_gaps=(core_eval - core_eval[0]).tolist(),
        meter_clock=calibration,
        meter_charge=charge(board, (*dimers[:2], dimers[-1])),
        rows=rows,
    )


def decode(data, setting, source, wait):
    [row] = [
        r
        for r in data["rows"]
        if (r["apparatus"], r["source"], r["wait"]) == (setting, source, wait)
    ]
    counts = np.asarray(row["core_by_meter"])
    assert counts.shape == (4, 2) and counts.sum() == data["shots_per_panel"]
    return int(counts[:, 1].sum()), int(counts[1:3, 1].sum())


def summarize(data):
    for name in ("native", "upstream_cut", "detector_cut", "no_reversal"):
        for arm in ("pulse", "idle"):
            print(name, arm, [decode(data, name, arm, t) for t in WAITS])
    print(
        "calorimeter gap and familywise half-width",
        data["meter_clock"]["gap"],
        data["meter_clock"]["half_width"],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=Path)
    args = parser.parse_args()
    if args.acquire:
        data = acquire()
        args.acquire.write_text(json.dumps(data, indent=1) + "\n")
    else:
        data = json.loads((HERE / "records.json").read_text())
    summarize(data)
