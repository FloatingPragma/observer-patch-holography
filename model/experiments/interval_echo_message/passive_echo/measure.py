"""Complete material-clock records for a passive out-and-return mark and a separate one-way bit."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from ..contact_deadline import measure as pilot
from ..ranged_deadlines.measure import COMMON, OTHER, ROUTE

HERE = Path(__file__).parent
GATEWAY, HAND = 1, 2
TICKS = (1, 2, 3, 4)
PAIRS = tuple(zip(ROUTE[1:], OTHER, strict=True))
PORTS = (GATEWAY, *(p for pair in PAIRS for p in pair))
SHOTS = 2**30


def histories(h, x, y, size):
    """One forward column per sharp release, with a real local tick projection after each slot."""
    hand_rows = np.flatnonzero(x == HAND)
    local = y[hand_rows]
    state = np.zeros((len(x), len(PORTS)), complex)
    state[hand_rows[np.searchsorted(local, PORTS)], np.arange(len(PORTS))] = 1
    if not np.array_equal(local, np.delete(np.arange(size), HAND)):
        raise ValueError("clock subspace ordering")
    gen = (-1j * h).tocsr()
    trace = gen.diagonal().sum()
    out = []
    for k in range(2 * max(TICKS) + 1):
        live = np.zeros((size, state.shape[1]), complex)
        live[local] = state[hand_rows]
        out.append(live)
        if k == 2 * max(TICKS):
            break
        if k == max(TICKS):
            state = state[:, :1].copy()
        state = expm_multiply(gen, state, traceA=trace)
        state[x != HAND] = 0
    return out


def dist(psi):
    mass = abs(psi) ** 2
    return np.r_[mass[:12], mass[12:].sum()]


def mass(psi):
    return np.vdot(psi, psi).real


def sample(labels, rates, rng):
    rates = np.asarray(rates, dtype=float)
    if np.min(rates) < -2e-13 or abs(rates.sum() - 1) > 1e-11:
        raise ValueError("incomplete history law")
    # Round far below one sampling count; retain all labelled bins, including zeros.
    p = np.maximum(0, np.round(rates, 15))
    counts = rng.multinomial(SHOTS, p / p.sum())
    return dict(zip(labels, map(int, counts), strict=True))


def signal(history, source, n, rng):
    j = PORTS.index(source)
    previous = 1.0
    rates = []
    for k in range(1, n + 1):
        current = mass(history[k][:, j])
        rates.append(previous - current)
        previous = current
    rates.extend(dist(history[n][:, j]))
    labels = [f"timeout/{k}" for k in range(1, n + 1)]
    labels.extend(f"local/{p}" for p in range(12))
    labels.append("quiet")
    return sample(labels, rates, rng)


def ranging(history, source, n, rng):
    j = PORTS.index(source)
    labels, rates = [], []
    previous = 1.0
    for k in range(1, n + 1):
        current = mass(history[k][:, 0])
        labels.append(f"outward/timeout/{k}")
        rates.append(previous - current)
        previous = current
    visited_amplitude = history[n][source, 0]
    for visited in (True, False):
        tag = "marked" if visited else "unmarked"
        # The remote local sharp read splits the successful clock branch into
        # U^n|home> on its marked port and the coherent complement. Linearity
        # propagates each subnormalised branch through the second n tick tests.
        for k in range(n + 1):
            psi = visited_amplitude * history[k][:, j]
            if not visited:
                psi = history[n + k][:, 0] - psi
            current = mass(psi)
            if k:
                labels.append(f"{tag}/timeout/{k}")
                rates.append(previous - current)
            previous = current
        for addr, p in enumerate(dist(psi)):
            labels.append(f"{tag}/" + ("quiet" if addr == 12 else f"local/{addr}"))
            rates.append(p)
    return sample(labels, rates, rng)


def run():
    pilot.native_gates()
    board, x, y, _, kinetic, contrast = pilot.flight()
    for c, (z, o) in zip(COMMON, PAIRS, strict=True):
        if (
            c // 12 != z // 12
            or c // 12 != o // 12
            or not {z, o} <= set(board.neighbours[c])
        ):
            raise ValueError("invalid sender-local exchange")
    rng = np.random.default_rng(731905)
    output = {"shots": SHOTS, "seed": 731905, "port_pairs": PAIRS, "banks": {}}
    for law, strength in (("native", 1), ("half_exchange", 0.5)):
        print("propagating", law, flush=True)
        history = histories(strength * kinetic + contrast, x, y, board.n)
        panel = {}
        for near, far in PAIRS:
            panel[str(near)] = {}
            for n in TICKS:
                panel[str(near)][str(n)] = {
                    "one_way": [signal(history, p, n, rng) for p in (near, far)],
                    "echo": [ranging(history, p, n, rng) for p in (near, far)],
                    "sham": [signal(history, near, n, rng) for _ in range(2)],
                }
        output["banks"][law] = panel
    return output


if __name__ == "__main__":
    arg = argparse.ArgumentParser(description=__doc__)
    arg.add_argument("--record", type=Path)
    args = arg.parse_args()
    result = json.dumps(run(), separators=(",", ":")) + "\n"
    if args.record:
        with args.record.open("x") as f:
            f.write(result)
    elif result != (HERE / "records.json").read_text():
        raise ValueError("reacquired records differ")
    print("complete")
