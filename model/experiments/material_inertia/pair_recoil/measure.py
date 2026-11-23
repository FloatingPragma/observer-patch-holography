"""Prepare a coherent native bound pair on opposite sides of a real bridge."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import jv

from model.executable_model import federation
from .. import neutral_defect as defect

HERE = Path(__file__).resolve().parent
SEED = 770409
PLACEMENTS = (0, 10, 20)
FLIGHTS = (1.0, 4.0, 12.0)
PRIMARY = (1 / 512, 1 / 1024)
CROSS = (1 / 8, 1 / 16, 1 / 32)
EPSILON = 2**-18


def propagation(h, sources, times):
    """Unfiltered Chebyshev expansion on the certified shifted spectrum [-16,44]."""
    center, radius = 14.0, 30.0
    # Even at the largest argument this tail is below the precision of the read.
    limit = math.ceil(
        radius * max(times) + 20 * (radius * max(times) + 1) ** (1 / 3) + 30
    )
    order = np.arange(limit + 1)
    times = np.asarray(times)
    coeff = 2 * (-1j) ** order[:, None] * jv(order[:, None], radius * times)
    coeff[0] /= 2
    series = [coeff[0, i] * sources for i in range(len(times))]
    previous = sources
    current = (h @ previous - center * previous) / radius
    for j in range(1, limit + 1):
        for i, c in enumerate(coeff[j]):
            if abs(c) > 1e-19:
                series[i] += c * current
        if j < limit:
            previous, current = (
                current,
                2 / radius * (h @ current - center * current) - previous,
            )
    for i, value in enumerate(series):
        value *= np.exp(1j * (37 - center) * times[i])
        np.testing.assert_allclose(
            value.conj().T @ value, sources.T @ sources, rtol=0, atol=3e-12
        )
    return dict(zip(times, series, strict=True))


def source_columns(board, bridge, h):
    idx, pairs = [], []
    for port in bridge:
        regular = min(
            v
            for v in board.neighbours[port]
            if v // 12 == port // 12 and board.degree[v] == 5
        )
        pairs.append([regular, port])
        idx.append(defect.ordered_index(regular, port, board.n))
    sharp = np.zeros((h.shape[0], 2))
    sharp[idx, [0, 1]] = 1.0
    lifted, bound = defect.high_component(h, sharp, 420)
    gram = lifted.T @ lifted
    values, basis = np.linalg.eigh(gram)
    prepared = lifted @ ((basis * (1 / np.sqrt(values))) @ basis.T)
    return pairs, sharp, prepared, gram, bound


def read_ramsey(probability, shots, seed):
    clicks = int(np.random.default_rng(seed).binomial(shots, probability))
    return dict(probability=float(probability), click=clicks, quiet=shots - clicks)


def spectral(h, prepared, flights, index):
    compressed = prepared.T @ (h @ prepared)
    primary, cross = [], []
    for i, t in enumerate(PRIMARY):
        records = []
        for j, sign in enumerate((1, -1)):
            initial = (prepared[:, 0] + sign * prepared[:, 1]) / math.sqrt(2)
            final = (flights[t][:, 2] + sign * flights[t][:, 3]) / math.sqrt(2)
            probability = (1 - np.imag(initial @ final)) / 2
            records.append(read_ramsey(probability, 2**49, [SEED, index, 1, i, j]))
        observed = (records[0]["click"] - records[1]["click"]) / (t * 2**49)
        se = (
            math.sqrt(
                sum(r["probability"] * (1 - r["probability"]) for r in records)
                / (2**49)
            )
            / t
        )
        primary.append(
            dict(time=t, counts=records, sampled_J=observed, standard_error=se)
        )
    for i, t in enumerate(CROSS):
        z = prepared[:, 1] @ flights[t][:, 2]
        records = []
        scale = math.sqrt(1 - EPSILON**2)
        for j, sign in enumerate((1, -1)):
            prob = abs(scale * z + sign * 1j * EPSILON) ** 2 / 2
            records.append(read_ramsey(prob, 2**59, [SEED, index, 4, i, j]))
        denom = -2 * t * scale * EPSILON
        measured = (records[0]["click"] - records[1]["click"]) / (2**59 * denom)
        se = math.sqrt(
            sum(r["probability"] * (1 - r["probability"]) for r in records) / 2**59
        ) / abs(denom)
        cross.append(
            dict(
                time=t,
                counts=records,
                computed_quadrature_J=float(-z.imag / t),
                sampled_J=float(measured),
                standard_error=se,
            )
        )
    return dict(
        matrix=compressed.tolist(), sine_comparison=primary, weak_reference=cross
    )


def detection(state, bins, signs, seed):
    prob = abs(state) ** 2
    categories = 3 * bins + signs + 1
    partition = np.bincount(categories, weights=prob, minlength=9).reshape(3, 3)
    counts = np.random.default_rng(seed).multinomial(
        2**36, (partition / partition.sum()).ravel()
    )
    return dict(probabilities=partition.tolist(), counts=counts.reshape(3, 3).tolist())


def bridge_study(board, bridge, h, x, y, high, index):
    pairs, sharp, prepared, gram, error = source_columns(board, bridge, h)
    flights = propagation(
        h, np.column_stack((sharp, prepared)), PRIMARY + CROSS + FLIGHTS
    )
    spectroscopy = spectral(h, prepared, flights, index)
    right = bridge[1] // 12
    left = bridge[0] // 12
    bins = np.full(len(x), 2, dtype=int)
    bins[high & (x // 12 == left) & (y // 12 == left)] = 0
    bins[high & (x // 12 == right) & (y // 12 == right)] = 1
    signs = (x // 12 == right).astype(int) - (y // 12 == right).astype(int)
    readings = []
    from scipy.linalg import expm

    for i, t in enumerate((0.0,) + FLIGHTS):
        data = np.column_stack((sharp, prepared)) if not t else flights[t]
        for j, name in enumerate(("sharp", "filtered")):
            cols = data[:, 2 * j : 2 * j + 2]
            records = [
                detection(
                    (cols[:, 0] + sign * 1j * cols[:, 1]) / math.sqrt(2),
                    bins,
                    signs,
                    [SEED, index, 2, i, j, k],
                )
                for k, sign in enumerate((1, -1))
            ]
            rates = np.array([r["probabilities"] for r in records]).sum(axis=-1)
            frequency = np.array([r["counts"] for r in records]).sum(axis=-1) / 2**36
            values = np.array([-1.0, 1.0, 0.0])
            signal = (rates[0] - rates[1]) @ values
            measured = (frequency[0] - frequency[1]) @ values
            se = np.sqrt(sum(r @ values**2 - (r @ values) ** 2 for r in rates) / 2**36)
            contrast = dict(
                exact=float(signal), sampled=float(measured), standard_error=float(se)
            )
            if name == "filtered":
                two = expm(-1j * t * np.array(spectroscopy["matrix"]))
                predicted = []
                for sign in (1, -1):
                    p = (
                        abs(
                            prepared @ (two @ np.array([1.0, 1j * sign])) / math.sqrt(2)
                        )
                        ** 2
                    )
                    predicted.append(np.bincount(bins, weights=p, minlength=3) @ values)
                contrast["two_site_prediction"] = float(predicted[0] - predicted[1])
            readings.append(
                dict(time=t, preparation=name, records=records, imbalance=contrast)
            )
    return dict(
        bridge=list(bridge),
        source_pairs=pairs,
        filter_gram=gram.tolist(),
        filter_error=error,
        energy=spectroscopy,
        flights=readings,
    )


def run():
    board = federation(0)
    x, y, adj, kinetic, h, cost = defect.generators(board)
    data = [
        bridge_study(board, board.external[j], h, x, y, cost == 4 * defect.KAPPA, j)
        for j in PLACEMENTS
    ]
    return dict(
        protocols={
            n: hashlib.sha256((HERE / n).read_bytes()).hexdigest()
            for n in ("PROTOCOL.md", "CROSS_READ_PROTOCOL.md")
        },
        registered_cells=board.cells,
        source_configurations=len(x),
        placements=data,
        ranges=dict(
            position=["left_cell_bound", "right_cell_bound", "elsewhere_or_unbound"],
            right_cell_relative_load=[-1, 0, 1],
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = run()
    if args.output:
        with args.output.open("x") as file:
            json.dump(data, file, indent=2, allow_nan=False)
            file.write("\n")
    else:
        old = json.loads((HERE / "results.json").read_text())

        def check(a, b):
            if isinstance(a, dict):
                assert a.keys() == b.keys()
                for k in a:
                    check(a[k], b[k])
            elif isinstance(a, list):
                for u, v in zip(a, b, strict=True):
                    check(u, v)
            elif isinstance(a, float):
                np.testing.assert_allclose(a, b, rtol=2e-8, atol=2e-13)
            else:
                assert a == b

        check(data, old)
    for row in data["placements"]:
        print(
            row["bridge"],
            row["energy"]["matrix"][0][1],
            row["energy"]["weak_reference"][1]["sampled_J"],
            [
                (r["time"], r["imbalance"])
                for r in row["flights"]
                if r["preparation"] == "filtered"
            ],
        )
