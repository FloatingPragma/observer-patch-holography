"""Independent, simulated binary energy records after an isolated bridge phase kick."""

import argparse
import json

import numpy as np
from scipy.sparse.linalg import expm_multiply

from . import bridge


def sums_for_charge(board, route, rings, q):
    owners = tuple(p // 12 for p in route[1:3])
    ports = [12 * owner + i for owner in owners for i in range(12)]
    assert sum(a // 12 in owners and a // 12 == b // 12 for a, b in board.seams) == 60
    idx = bridge.beam.indices(ports)
    src = bridge.beam.indices(list(rings[0]) + list(rings[1]))
    b = bridge.bridges(board, [route[1:3]])[idx, :][:, idx].toarray() / 2
    level, basis = np.linalg.eigh(b)
    phase = 2 * np.pi * np.arange(19) / 19
    freq = np.arange(-9, 10)
    weights = (np.exp(-1j * np.outer(phase, freq)) @ (1j * freq)).real / (2 * 19)
    after_kick = {mode: [] for mode in ("active", "zero_contrast_strength")}
    for odd in (False, True):
        cells = bridge.operators(board, route, q, odd)
        diff = bridge.masked(cells, [owners[1]]) - bridge.masked(cells, [owners[0]])
        d = diff[idx, :][:, idx].toarray() / 2
        for mode in after_kick:
            if mode == "active":
                ht = -2j * bridge.native.generator(board, route, q, odd)
            else:
                # Frozen source and active phase/energy meters. Only the free
                # flight has its contrast-squared strength removed.
                ht = -4j * bridge.native.iso.ranging.generator(board, route, q, 0, odd)
            launch = np.zeros((ht.shape[0], len(src)))
            launch[src] = np.eye(len(src))
            w = expm_multiply(ht, launch, traceA=ht.diagonal().sum())[idx]
            samples = []
            for theta in phase:
                u = (basis * np.exp(-1j * theta * level)) @ basis.T
                v = u @ w
                effects = [arm.conj().T @ d @ arm for arm in np.hsplit(v, 10)]
                samples.append(np.mean(effects, axis=0))
            after_kick[mode].append(samples)
    return weights, after_kick


def run():
    seed = 19262738
    rng = np.random.default_rng(seed)
    n = 2**60
    rows = {}
    states = bridge.action.preparations()
    tangent = bridge.tangents()
    board = bridge.native.iso.pairing.formation.federation(1)
    route = bridge.native.iso.pairing.formation.region(board)
    rings = bridge.beam.screens(board, route)
    angles = {
        ("ordinary" if row["source"] == "ordinary" else "pairing"): row[
            "angle_per_unit"
        ]
        for row in json.loads(
            (bridge.HERE.parent / "signed_transfers.json").read_text()
        )["tests"]
    }
    # Select one of 60 internal seam terms and 24 individual degree terms.
    # Each restricted term has spectral radius <=64, including sector scalars
    # removed at the two marked seams by a reversal-label dependent offset.
    magnitude = (60 + 24) * 64
    term_error = magnitude * np.sqrt(2 * np.log(2 * 19 * len(states) * 2 / 0.01) / n)
    for q in (-1, 1):
        weights, masks = sums_for_charge(board, route, rings, q)
        settings = {}
        raw = {}
        for mode, effects in masks.items():
            measured = np.einsum("t,btij->bij", weights, np.array(effects))
            raw[mode] = bridge.evaluate(measured, tangent)
            if mode == "active":
                for name, state in states.items():
                    means = np.einsum("btij,bji->t", np.array(effects), state).real
                    probability = (1 + means / magnitude) / 2
                    assert np.all((0 < probability) & (probability < 1))
                    counts = rng.binomial(n, probability)
                    estimate = magnitude * (2 * (counts / n) - 1) @ weights
                    settings[name] = dict(
                        plus_counts=counts.tolist(), estimated_current=float(estimate)
                    )
        d = {k: v["estimated_current"] for k, v in settings.items()}
        slopes = {}
        for name, angle in angles.items():
            estimate = d[name] - d["equilibrium"] + angle * (d["plus"] - d["minus"]) / 2
            halfwidth = term_error * np.linalg.norm(weights, ord=1) * (2 + abs(angle))
            assert abs(estimate - raw["active"][name]) < halfwidth
            slopes[name] = dict(
                estimate=estimate, simultaneous_halfwidth=float(halfwidth)
            )
        rows[str(q)] = dict(
            phase_settings=settings, slopes=slopes, exact_mean_predictions=raw
        )
    return dict(
        seed=seed,
        angles="2 pi m / 19, m=0..18",
        trials_per_setting=n,
        simultaneous_coverage=0.99,
        binary_readout_magnitude=magnitude,
        charges=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (bridge.HERE / "phase_read.json").write_text(data)
    else:
        print(data, end="")
