"""Determine the gap from independent phase counts and transfer it to a held daughter assay."""

import argparse
import json

import numpy as np
from scipy.sparse import csr_matrix, eye
from scipy.sparse.linalg import expm_multiply

from . import measure as m

WIDTH = float(np.sqrt(2 * np.log(2 * 500 / 0.01) / m.N))


def read_clock(counts):
    low, high = 0.0, 256.0
    for j, (nx, ny) in enumerate(counts):
        t = 2**j / 128
        x, y = 2 * np.array([nx, ny]) / m.N - 1
        magnitude = np.hypot(x, y)
        radius = np.arcsin(np.sqrt(2) * WIDTH / magnitude)
        phase = np.arctan2(y, x)
        least = int(np.ceil((t * low - phase - radius) / (2 * np.pi)))
        last = int(np.floor((t * high - phase + radius) / (2 * np.pi)))
        m.check(least == last, "phase is unresolved on this frequency arm")
        angle = phase + 2 * np.pi * least
        low = max(low, (angle - radius) / t)
        high = min(high, (angle + radius) / t)
        m.check(low <= high, "phase observations are incompatible")
    return [float(low), float(high)]


def frequency(gap, rng):
    ts = 2.0 ** np.arange(28) / 128
    probabilities = (1 + np.stack([np.cos(gap * ts), np.sin(gap * ts)], axis=-1)) / 2
    counts = rng.binomial(m.N, np.round(probabilities, 14))
    return dict(quadrature_counts=counts.tolist(), gap_interval=read_clock(counts))


def terminal(wave, pairs, aims, rng):
    probs = []
    for state, pair, aim in zip(wave, pairs, aims, strict=True):
        _, words, embed, e, v, _, levels = pair
        d = (e - levels[0] - aim) * m.DT
        size = m.CLOCK_BINS[-1]
        coeff = (
            np.exp(-0.5j * (size - 1) * d)
            * np.sinc(size * d / (2 * np.pi))
            / np.sinc(d / (2 * np.pi))
        )
        # The observer holds a second, ground-prepared copy in the reference arm.
        selected = v @ (coeff[:, None] * (v.T @ embed @ state))
        p = (abs(selected) ** 2).sum(axis=1)
        probs.extend([*p, max(0.0, float(np.vdot(state, state).real - p.sum()))])
    count = m.sample([*probs, 1 - float(np.vdot(wave, wave).real)], rng)
    words_count = count[:-1].reshape(2, -1)[:, :-1].sum(axis=0)
    accepted = int(words_count.sum())
    charge = words[:, -2:].sum(axis=1) - 2
    return dict(
        complete_counts=count.tolist(),
        line_counts=accepted,
        conditional_rate=m.interval(accepted, m.N - int(count[-1])),
        outer_signed_counts=[int(words_count[charge == q].sum()) for q in (-2, 0, 2)],
    )


def run():
    rng = np.random.default_rng(9061202)
    boundary, initial, order = m.reference()
    h = csr_matrix(m.walls.operators(8, False, 3))
    wave = expm_multiply(
        -32j * (h - eye(len(initial)) * np.mean(h.diagonal())), initial
    ).reshape(16, 8, 8, 2)[order]
    rows = {}
    for mode in m.SETTINGS:
        pairs = [
            m.wing(
                8, 0, int(boundary[0] if phase == 0 else 2 - boundary[0]), False, mode
            )
            for phase in range(2)
        ]
        clocks = [frequency(p[-1][1] - p[-1][0], rng) for p in pairs]
        aims = [np.mean(c["gap_interval"]) for c in clocks]
        rows[mode] = dict(
            clocks=clocks,
            held=[terminal(wave[:, :, :, i], pairs, aims, rng) for i in (0, 1)],
        )
    return dict(
        shots=m.N, initial_frequency_bracket=[0, 256], seed=9061202, settings=rows
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=m.Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as f:
            json.dump(result, f, indent=1, allow_nan=False)
            f.write("\n")
    else:
        m.compare(
            result, json.loads((m.HERE / "clocks.json").read_text()), "daughter clock"
        )
        print(
            "Replayed the count-calibrated soft and hard clocks and the held material words."
        )
