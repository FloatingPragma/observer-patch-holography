"""Replay the held protocol on the historical ports of the entire twelve-port cell."""

import argparse
import json
from itertools import combinations_with_replacement
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from . import measure as m

WORDS = list(combinations_with_replacement(range(12), 2))
X, W = WORDS.index((2, 3)), WORDS.index((0, 1))
ROWS = {word: i for i, word in enumerate(WORDS)}
ROUTE = ((5, 4), (4, 2), (7, 6), (6, 3))


def lift(block, edge):
    if tuple(sorted(edge)) not in m.carrier.CELL_SEAMS:
        raise ValueError("Nonlocal contact")
    a, b = edge
    out = np.zeros((len(WORDS), len(WORDS)), dtype=complex)
    for i, word in enumerate(WORDS):
        key = (word.count(a), word.count(b))
        spectators = [p for p in word if p not in edge]
        for j, val in enumerate(block[:, m.LAW.SRC.index(key)]):
            if val == 0:
                continue
            left, right = m.LAW.SRC[j]
            if left + right != sum(key):
                raise ValueError("Conservation violated")
            dest = tuple(sorted(spectators + [a] * left + [b] * right))
            out[ROWS[dest], i] += val
    return out


def phase(port, angle):
    return np.diag([np.exp(1j * word.count(port) * angle) for word in WORDS])


def preparation():
    local = m.LAW.seam_term(m.candidate.KAPPA)

    def flight(e, t):
        return lift(expm(-1j * t * local), e)

    route = np.eye(len(WORDS))
    for e in ROUTE:
        route = flight(e, np.pi) @ route
    local_source = np.eye(len(WORDS))
    for e in ((0, 5), (1, 7)):
        local_source = phase(e[0], -np.pi / 2) @ flight(e, np.pi / 2) @ local_source
    refocus = phase(0, np.pi) @ flight((0, 1), np.pi / (4 * (m.candidate.KAPPA + 0.5)))
    local_source = flight((1, 7), np.pi / 2) @ refocus @ refocus @ local_source
    return route @ local_source @ route.conj().T


def run(commissioning):
    source = preparation()
    state = source[:, X]
    expected = np.zeros(len(WORDS))
    expected[[W, X]] = 0.5
    np.testing.assert_allclose(abs(state) ** 2, expected, atol=2e-14)
    np.testing.assert_allclose(state[W], state[X], atol=2e-14)
    instruments = {}
    for name, u in (("native", m.LAW.echo()), ("total", m.total_sensitive())):
        blocks = u[:, : m.LAW.N_SRC].reshape(-1, m.LAW.N_SRC, m.LAW.N_SRC)
        instruments[name] = [lift(block, (0, 1)) for block in blocks if np.any(block)]
    rng = np.random.default_rng(914731)
    budget = m.rate_radius(len(WORDS) + len(m.WORDS) + 2)
    runs = {}
    occupation = np.array([word.count(0) + word.count(1) for word in WORDS])
    for mode in m.MODES:
        if mode == "bypass":
            branches = [state]
            rejected = 0
        else:
            ks = instruments["total" if mode == "total" else "native"]
            allowed = {"only_empty": (0,), "only_full": (2,), "which": (0, 1, 2)}.get(
                mode
            )
            groups = (
                [state]
                if allowed is None
                else [state * (occupation == n) for n in allowed]
            )
            rejected = (
                0
                if allowed is None
                else sum(abs(state[~np.isin(occupation, allowed)]) ** 2)
            )
            branches = [k @ s for s in groups for k in ks]
        cal = commissioning["calibration"][mode]
        offset = 0.25 if mode.startswith("only_") else 0.5
        angles = []
        for angle in m.PHASES:
            read = source.conj().T @ phase(0, angle)
            probabilities = np.r_[sum(abs(read @ x) ** 2 for x in branches), rejected]
            np.testing.assert_allclose(probabilities.sum(), 1, atol=4e-14)
            probabilities = np.round(probabilities, 13)
            counts = rng.multinomial(m.SHOTS, probabilities / probabilities.sum())
            c, s = np.cos(angle), np.sin(angle)
            pred = offset + cal["cosine"] * c + cal["sine"] * s
            angles.append(
                dict(
                    phase=float(angle),
                    counts=counts.tolist(),
                    bright_forecast=float(pred),
                    residual=float(abs(counts[X] / m.SHOTS - pred)),
                    joint_error=float(budget * (1 + abs(c) + abs(s))),
                )
            )
        runs[mode] = angles
    return json.loads(
        json.dumps(
            dict(
                words=WORDS,
                shots=m.SHOTS,
                route=ROUTE,
                bright_word=WORDS[X],
                combined_rate_radius=budget,
                banks=runs,
            )
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    obs = run(json.loads((m.HERE / "results.json").read_text()))
    if args.output:
        args.output.write_text(json.dumps(obs, indent=2) + "\n")
    elif obs != json.loads((m.HERE / "historical.json").read_text()):
        raise ValueError("Historical-address replay differs")
    for mode in m.MODES:
        print(mode, [r["counts"][X] for r in obs["banks"][mode]])
