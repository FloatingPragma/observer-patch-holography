"""Cooling measured free-contact modes on a four-child address refinement."""

import argparse
import json
from itertools import combinations
from math import pi
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo.candidate import KAPPA, native_law
from model.executable_model import CELL_SEAMS, federation

HERE = Path(__file__).resolve().parent
STEPS = (0, 1, 3, 6)
SHOTS_PER_RAIL = 1024
FLIGHT_TIMES = (0, pi / 3, pi / 2, pi, 3 * pi / 2)


def local_paths():
    links = {frozenset(edge) for edge in CELL_SEAMS}
    return [
        (a, b, c, d)
        for a, b in sorted(tuple(edge) for edge in CELL_SEAMS)
        for c in range(12)
        for d in range(12)
        if len({a, b, c, d}) == 4
        and frozenset((b, c)) in links
        and frozenset((c, d)) in links
    ]


def layouts():
    paths = local_paths()
    return {
        blocked: next(p for p in paths if set(p).isdisjoint(blocked))
        for blocked in combinations(range(12), 3)
    }


def gates(mode):
    exchange, contrast = native_law.local_term_parts()
    h = KAPPA * contrast + (0 if mode == "exchange_off" else exchange)
    ports = np.eye(3, dtype=int) * 2
    blocks = []
    for a, b in [(0, 1), (1, 2)]:
        ix = [native_law.SRC.index((w[a], w[b])) for w in ports]
        if np.max(np.abs(np.delete(h[:, ix], ix, axis=0))) > 1e-13:
            raise ValueError("The three material addresses are not invariant")
        blocks.append(h[np.ix_(ix, ix)])
    phase = np.diag(np.exp(-1j * pi / (8 * KAPPA) * KAPPA * ports[:, 1] ** 2))
    decoder = expm(-1j * pi / 2 * blocks[0]) @ phase
    leak = expm(-1j * (0 if mode == "parked_filter" else pi / 2) * blocks[1])
    return blocks[0], decoder, leak


def pointer():
    written = native_law.echo()
    bright = np.tile(
        [pair == (1, 1) for pair in native_law.SRC],
        native_law.DIM // len(native_law.SRC),
    )
    for value in (0, 2):
        i = native_law.SRC.index((value, 0))
        branch = written[:, i] * (bright if value == 2 else ~bright)
        np.testing.assert_allclose(
            written.conj().T @ branch, np.eye(native_law.DIM)[:, i], atol=1e-13, rtol=0
        )


def outcomes(decoder, leak, start, steps):
    state = decoder[:, start].copy()
    probs = []
    for _ in range(steps):
        state = leak @ state
        probs.append(abs(state[2]) ** 2)
        state[2] = 0  # unnormalized dark branch, so every original start is retained
    restored = decoder.conj().T @ state
    terminal = abs(decoder @ restored) ** 2
    np.testing.assert_allclose(terminal[2], 0, atol=1e-13)
    p = np.r_[probs, terminal[:2]]
    np.testing.assert_allclose(p.sum(), 1, atol=1e-13, rtol=0)
    return p / p.sum()


def run():
    rng = np.random.default_rng(850082)
    pointer()
    lookup = layouts()
    all_boards = []
    for level in range(3):
        board = federation(level)
        layouts_here = []
        for owner in range(board.cells):
            gateways = tuple(
                np.flatnonzero(board.degree[12 * owner : 12 * owner + 12] == 6)
            )
            layouts_here.append(lookup[gateways])
        panels = []
        for mode in ("native", "parked_filter", "exchange_off"):
            h, decoder, leak = gates(mode)
            flight = []
            for t in FLIGHT_TIMES:
                p = float(abs(expm(-1j * t * h)[1, 0]) ** 2)
                flight.append(int(rng.binomial(65536, np.clip(p, 0, 1))))
            banks = []
            for steps in STEPS:
                # One independent cohort starts on each sharp rail of every owner.
                counts = np.stack(
                    [
                        rng.multinomial(
                            SHOTS_PER_RAIL,
                            outcomes(decoder, leak, j, steps),
                            size=board.cells,
                        )
                        for j in (0, 1)
                    ],
                    axis=1,
                )
                banks.append(dict(steps=steps, counts=counts.tolist()))
            panels.append(
                dict(
                    mode=mode,
                    energy_levels=np.linalg.eigvalsh(h[:2, :2]).tolist(),
                    flight_transfers=flight,
                    banks=banks,
                )
            )
        all_boards.append(dict(owners=board.cells, layouts=layouts_here, panels=panels))
    return dict(
        shots_per_rail=SHOTS_PER_RAIL,
        flight_shots=65536,
        times=FLIGHT_TIMES,
        boards=all_boards,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    path = HERE / "records.json"
    encoded = json.dumps(result, separators=(",", ":")) + "\n"
    if args.write:
        path.write_text(encoded)
    elif encoded != path.read_text():
        raise SystemExit("Record mismatch")
    for board in result["boards"]:
        for panel in board["panels"]:
            last = np.array(panel["banks"][-1]["counts"]).sum(axis=(0, 1))
            print(board["owners"], panel["mode"], last[-2:].tolist(), last.sum())
