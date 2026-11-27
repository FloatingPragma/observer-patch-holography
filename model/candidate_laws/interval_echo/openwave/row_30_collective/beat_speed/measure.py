"""Measure velocity from moving adjacent-tone fringes in a passive polar material."""

import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from model.experiments.material_inertia.pair_spectrum.neutral_guide.polar_ring import (
    measure as polar,
    source,
    stock,
)
from ..polar_recoil.measure import angle_error, line_error, unwrap


HERE = Path(__file__).resolve().parent
SHOTS = 2**22
WAITS = [0, *(2**j for j in range(14))]


def character(n, k):
    return np.exp(2j * np.pi * k * np.arange(n) / n) / np.sqrt(n)


def bank(rng, one):
    probs = abs(one) ** 2
    complete = np.r_[probs, 1 - probs.sum()]
    if complete.min() < -1e-10:
        raise ValueError("Spatial probabilities exceed one")
    complete = complete.clip(0).round(13)
    return rng.multinomial(SHOTS, complete / complete.sum()).tolist()


def launch(a, b, phase=1):
    amplitudes = (a + phase * b) / np.sqrt(2)
    sites = 1 << np.arange(len(a))
    return source.prepared(amplitudes)[sites]


def acquire():
    rng = np.random.default_rng(9302601)
    data = []
    for case in stock.charts():
        if case["cells"] != 80 or not case["stock"]:
            continue
        ring = [list(dimer) for dimer in case["stock"]["ring"]]
        n = len(ring)
        uniform = character(n, 0)
        shape = {"0": [bank(rng, source.prepared(uniform)[1 << np.arange(n)])]}
        for k in range(1, n):
            shape[str(k)] = [
                bank(rng, launch(uniform, character(n, k), phase)) for phase in (1, 1j)
            ]
        flights = {}
        for contrast in (8, 4):
            _, blocks = polar.spectrum(n, contrast)
            returns, main = [], []
            for _, energies, _, origin, principal in blocks:
                weights = n * abs(origin) ** 2
                returns.append(np.exp(-1j * np.outer(WAITS, energies)) @ weights)
                main.append(float(weights[principal]))
            rows = []
            for direction in (1, -1):
                for j in range(n // 2):
                    a, b = direction * j % n, direction * (j + 1) % n
                    ca, cb = character(n, a), character(n, b)
                    start = launch(ca, cb)
                    release = (
                        returns[a][:, None] * ca + returns[b][:, None] * cb
                    ) / np.sqrt(2)
                    np.testing.assert_allclose(release[0], start, atol=2e-13)
                    rows.append(
                        dict(
                            carriers=[a, b], banks=[bank(rng, word) for word in release]
                        )
                    )
            flights[str(contrast)] = dict(line_weights=main, rows=rows)
        data.append(dict(ring=ring, sources=shape, flights=flights))
    return dict(shots=SHOTS, waits=WAITS, materials=data)


def moments(banks, n, appointments):
    counts = np.asarray(banks)
    if (
        counts.shape != (appointments, n + 1)
        or counts.dtype.kind not in "iu"
        or np.any(counts < 0)
        or np.any(counts.sum(axis=1) != SHOTS)
    ):
        raise ValueError("Incomplete terminal count bank")
    return np.fft.fft(counts[:, :-1] / SHOTS, axis=1)


def identify(records, n):
    if set(records) != set(map(str, range(n))):
        raise ValueError("Missing source character meter")
    if np.max(abs(moments(records["0"], n, 1)[0, 1:])) >= 0.01:
        raise ValueError("Spatial modulation in uniform control")
    observed = {0: 0}
    for k in range(1, n):
        cosine, minus_sine = moments(records[str(k)], n, 2)
        strength = abs((cosine - 1j * minus_sine)[1:])
        if np.sort(strength)[-1] - np.sort(strength)[-2] <= 0.1:
            raise ValueError("Unresolved signed spatial character")
        observed[k] = int(np.argmax(strength)) + 1
        if observed[k] != k:
            raise ValueError("Source preparation disagrees with spatial receiver")
    return observed


def read(data):
    if data["shots"] != SHOTS or data["waits"] != WAITS:
        raise ValueError("Unrecognized observation schedule")
    summary = {}
    if sorted(len(r["ring"]) for r in data["materials"]) != [10, 12]:
        raise ValueError("Both polar rings are required")
    for material in data["materials"]:
        n = len(material["ring"])
        wave = identify(material["sources"], n)
        if set(material["flights"]) != {"8", "4"}:
            raise ValueError("Both material settings are required")
        for contrast, flight in material["flights"].items():
            rays = {}
            for item in flight["rows"]:
                a, b = item["carriers"]
                step = (wave[b] - wave[a]) % n
                if step not in (1, n - 1):
                    raise ValueError("Non-adjacent measured wavelengths")
                sign = 1 if step == 1 else -1
                j = (sign * wave[a]) % n
                if (sign, j) in rays:
                    raise ValueError("Duplicate or reversed carrier pair")
                z = moments(item["banks"], n, len(WAITS))
                if np.argmax(abs(z[0, 1 : n // 2 + 1])) + 1 != min(step, n - step):
                    raise ValueError(
                        "Measured beat period differs from source separation"
                    )
                spread = line_error(flight["line_weights"][a]) + line_error(
                    flight["line_weights"][b]
                )
                obs = [
                    (v.conjugate(), spread + angle_error(v, SHOTS)) for v in z[1:, step]
                ]
                gap = np.asarray(unwrap((-1, 1), WAITS[1:], obs))
                rays[sign, j] = gap / (2 * np.pi / n)
            if set(rays) != {(s, j) for s in (1, -1) for j in range(n // 2)}:
                raise ValueError("Missing outward momentum interval")
            for sign in (1, -1):
                bands = [rays[sign, j] for j in range(n // 2)]
                summary[f"{n}/{contrast}/{sign}"] = dict(
                    speed=[v.tolist() for v in bands],
                    near_rest_increase_at_least=float(bands[1][0] - bands[0][1]),
                    near_edge_decrease_at_least=float(bands[-2][0] - bands[-1][1]),
                )
    return summary


def load():
    with gzip.open(HERE / "counts.json.gz", "rt") as source_file:
        return json.load(source_file)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument(
        "--output", type=Path, help="Write new counts to this unused gzip path"
    )
    args = cli.parse_args()
    records = acquire()
    result = read(records)
    if args.output:
        with gzip.open(args.output, "xt") as out:
            json.dump(records, out, separators=(",", ":"), allow_nan=False)
    elif records != load():
        raise ValueError("Acquisition differs from archived counts")
    print(json.dumps(result, indent=2, allow_nan=False))
