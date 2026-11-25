"""Extract and test the finite Euclidean slab from addressed Interval Echo records."""

import argparse
import hashlib
import json
from fractions import Fraction
from itertools import product
from math import comb
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
from model.candidate_laws.interval_echo import observations
from model.experiments.interval_echo_consensus.refinement_cooling.dimer_band import (
    measure as dimer,
)
from model.executable_model import federation


def source():
    law = observations.native_law
    e, v = law.local_term_parts()
    h = e + observations.KAPPA * v
    cols = [law.SRC.index(w) for w in ((2, 0), (0, 2))]
    h2 = h[np.ix_(cols, cols)]
    assert np.linalg.norm(np.delete(h[:, cols], cols, axis=0)) < 1e-12
    dimer.pointer()  # native echo/projector/undo reads the occupied rail

    def sharp(theta):
        prop = expm(-1j * theta * h)[:, cols]
        assert np.linalg.norm(np.delete(prop, cols, axis=0)) < 1e-12
        return abs(prop[cols]) ** 2

    kernel = sharp(np.pi / 3)
    ground, upper = np.linalg.eigvalsh(h2)
    survival = np.linalg.eigvalsh(kernel)[0]
    dt = -np.log(survival) / (upper - ground)
    transfer = expm(-dt * (h2 - ground * np.eye(2)))
    np.testing.assert_allclose(kernel, transfer, atol=1e-13, rtol=0)
    # A fresh fair enable followed by a fair half exchange if enabled.
    np.testing.assert_allclose((np.eye(2) + sharp(np.pi / 2)) / 2, transfer, atol=1e-13)
    changed_h = observations.KAPPA * v
    changed_columns = expm(-1j * np.pi / 3 * changed_h)[:, cols]
    np.testing.assert_allclose(
        abs(changed_columns) ** 2, np.eye(law.DIM)[:, cols], atol=1e-13
    )
    changed_gap = np.ptp(np.linalg.eigvalsh(changed_h[np.ix_(cols, cols)]))
    np.testing.assert_allclose(changed_gap, 0, atol=1e-13)
    return kernel, dict(
        hamiltonian=h2.tolist(),
        exchange_off=dict(
            gap=float(changed_gap),
            sharp_transition=(abs(changed_columns[cols]) ** 2).tolist(),
        ),
        eigenvalues=[ground, upper],
        sharp_transition=kernel.tolist(),
        slab_duration=dt,
        transfer_residual=float(abs(transfer - kernel).max()),
        uninterrupted_two=float(sharp(2 * np.pi / 3)[1, 0]),
        monitored_two=float((kernel @ kernel)[1, 0]),
        full_exchange_eigenvalues=np.linalg.eigvalsh(sharp(np.pi)).tolist(),
        parked=sharp(0).tolist(),
    )


def refinement():
    # Exact product pushforward, allowing all hidden inputs to vary.
    words = list(product(range(2), repeat=4))
    p = Fraction(1, 4)
    worst = Fraction(0)
    parity = set()
    for start in words:
        rows = []
        for stop in words:
            diff = sum(a != b for a, b in zip(start, stop))
            rows.append(p**diff * (1 - p) ** (4 - diff))
        assert sum(rows) == 1
        for keep in range(4):
            flipped = sum(w for w, b in zip(rows, words) if b[keep] != start[keep])
            worst = max(worst, abs(flipped - p))
        parity.add(sum(w for w, b in zip(rows, words) if sum(b) % 2 != sum(start) % 2))
    assert parity == {Fraction(15, 32)}
    return dict(
        inputs=len(words),
        retained_child_choices=4,
        max_flip_probability_error=str(worst),
        same_duration_parity_flip=str(parity.pop()),
    )


def make_log(owners, initial, flips, ids):
    state = list(map(int, initial))
    writers = [None] * len(state)
    events = []
    for owner, flip, name in zip(owners, flips, ids, strict=True):
        before = state[owner]
        after = before ^ int(flip)
        events.append(
            dict(
                id=name,
                register=int(owner),
                read=before,
                writer=writers[owner],
                write=after,
            )
        )
        state[owner] = after
        writers[owner] = name
    return events


def authenticate(events, initial):
    state = list(map(int, initial))
    writers = [None] * len(state)
    height = {}
    for event in events:
        name, owner = event["id"], event["register"]
        if (
            not (0 <= owner < len(state))
            or name in height
            or event["write"] not in (0, 1)
        ):
            raise ValueError("Invalid owner, identifier or rail value")
        parent = event["writer"]
        if parent != writers[owner] or event["read"] != state[owner]:
            raise ValueError("Stale writer, wrong register or discontinuous value")
        height[name] = 0 if parent is None else height[parent] + 1
        state[owner], writers[owner] = event["write"], name
    return height


def worst_word_distance(owners, appointments):
    """Exact total variation from a sharp word, maximal by convexity and symmetry."""
    scale = 2**appointments
    delta = sum(
        comb(owners, k)
        * abs((scale - 1) ** k * (scale + 1) ** (owners - k) - scale**owners)
        for k in range(owners + 1)
    )
    return Fraction(delta, 2 * (2 * scale) ** owners)


def one_percent_mixing(owners):
    t = 0
    previous = distance = worst_word_distance(owners, t)
    while distance > Fraction(1, 100):
        t += 1
        previous, distance = distance, worst_word_distance(owners, t)
    return dict(
        appointments=t, previous_distance=float(previous), distance=float(distance)
    )


def capture(kernel):
    lookup = dimer.layouts()
    rng = np.random.default_rng(1158201)
    boards = []
    for level in range(3):
        board = federation(level)
        flats = []
        for owner in range(board.cells):
            off = 12 * owner
            blocked = tuple(np.flatnonzero(board.degree[off : off + 12] == 6))
            flats.extend(off + i for i in lookup[blocked])
        assert len(set(flats)) == 4 * board.cells
        # Sample the full monitored joint trajectory on each sharp input.
        # Counts per owner avoid pretending that refinement reuses a coarse launch.
        panels = []
        for start in (0, 1):
            word = np.zeros((board.cells, 256), np.uint8)
            previous = np.full_like(word, start)
            for _ in range(3):
                after = (rng.random(previous.shape) < kernel[1, previous]).astype(
                    np.uint8
                )
                word = 2 * word + after
                previous = after
            panels.append(np.stack([(word == w).sum(axis=1) for w in range(8)], axis=1))
        boards.append(
            dict(
                level=level,
                owners=board.cells,
                worst_word_one_percent=one_percent_mixing(board.cells),
                shots_per_start_per_owner=256,
                address_sha256=hashlib.sha256(
                    np.asarray(flats, dtype="<i4").tobytes()
                ).hexdigest(),
                counts_by_word=np.sum(panels, axis=1).tolist(),
                counts=np.array(panels),
            )
        )
    return boards


def run(out):
    kernel, local = source()
    boards = capture(kernel)
    init = [0, 1, 0, 1]
    schedule = [0, 1, 2, 3, 2, 0, 0, 3, 3, 1, 0, 0]
    flips = np.random.default_rng(82115).random(len(schedule)) < kernel[1, 0]
    events = make_log(schedule, init, flips, [f"r{x}" for x in range(len(schedule))])
    heights = authenticate(events, init)
    independent_order = sorted(events, key=lambda e: e["register"])
    assert authenticate(independent_order, init) == heights
    out.mkdir(parents=True, exist_ok=True)
    counts = out / "counts.npz"
    np.savez_compressed(counts, **{f"L{b['level']}": b.pop("counts") for b in boards})
    receipt = dict(
        source=local,
        refinement=refinement(),
        boards=boards,
        log=dict(initial=init, events=events, heights=heights),
        counts_sha256=hashlib.sha256(counts.read_bytes()).hexdigest(),
        inputs={
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "PROTOCOL.md",
                Path(__file__),
                Path(dimer.__file__),
                Path(observations.native_law.__file__),
                Path(observations.__file__),
            )
        },
    )
    (out / "results.json").write_text(
        json.dumps(receipt, indent=2, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=HERE)
    run(parser.parse_args().out)
