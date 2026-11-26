"""Measure formation and unordered joint arrivals after two local island flips."""

import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np

from . import measure
from .. import measure as ring

TIMES = (0, 2, 8, 16, 32, 64, 128)
SHOTS = 2**22
LABELS = 16384
PERIOD = 512


def bank(probs, success, rng):
    probs = np.asarray(probs, dtype=float)
    rest = 1 - probs.sum()
    assert probs.min() >= -1e-12 and rest >= -1e-12
    all_probs = np.r_[
        1 - success, success * np.maximum(probs, 0), success * max(0.0, rest)
    ]
    return rng.multinomial(SHOTS, all_probs / all_probs.sum()).tolist()


def window(energy, center):
    bins = np.arange(
        np.ceil((center - 1) * LABELS / PERIOD),
        np.floor((center + 1) * LABELS / PERIOD) + 1,
    )
    phase = 2 * np.pi * (energy[:, None] / PERIOD - bins / LABELS)
    return (abs(measure.aperture(phase, LABELS)) ** 2).sum(axis=1)


def panel(base):
    n, kappa = base["dimer_count"], base["kappa"]
    rng = np.random.default_rng([970428, n, int(kappa)])
    _, blocks = ring.spectrum(n, kappa)
    cold, success = measure.sharp_cold(blocks[0])
    frame, e, v, _, _ = blocks[0]
    source = frame @ (v @ cold)
    words = np.arange(len(source))
    ds = (1, 2, n // 2)
    prepares = [source[words ^ 1], *(source[words ^ 1 ^ (1 << d)] for d in ds)]
    inputs = np.column_stack(prepares)
    waves = [np.zeros_like(inputs) for t in TIMES]
    centers = [2 * ring.read_clock(base["clocks"][0]["plus_counts"]), 6 * kappa]
    freq = np.zeros((len(prepares), 2))
    for tr, e, v, _, _ in blocks:
        assert min(e) > -PERIOD / 2 and max(e) < PERIOD / 2
        amp = v.conj().T @ (tr.conj().T @ inputs)
        freq += (abs(amp) ** 2).T @ np.column_stack([window(e, c) for c in centers])
        for t, wave in zip(TIMES, waves, strict=True):
            wave += tr @ (v @ (np.exp(-1j * e * t)[:, None] * amp))
    energy = [bank(f, success, rng) for f in freq[1:]]
    addr = np.array(list(combinations(range(n), 2)))
    pair_words = (1 << addr[:, 0]) + (1 << addr[:, 1])
    rows = []
    for t, wave in zip(TIMES, waves, strict=True):
        np.testing.assert_allclose((abs(wave) ** 2).sum(axis=0), 1, atol=1e-12)
        # Independent detector tapes for each of the two comparison copies.
        singleton = abs(wave[1 << np.arange(n), 0]) ** 2
        rows.append(
            dict(
                time=t,
                reference_left=bank(singleton, success, rng),
                reference_right=bank(singleton, success, rng),
                pairs=[
                    bank(abs(wave[pair_words, j + 1]) ** 2, success, rng)
                    for j in range(len(ds))
                ],
            )
        )
    return dict(
        n=n,
        kappa=kappa,
        separations=ds,
        frequency_centers=centers,
        energy=energy,
        word_order=addr.tolist(),
        spatial=rows,
    )


def conditional(row):
    return np.asarray(row[1:]) / (sum(row) - row[0])


def compare(p):
    i, j = np.array(p["word_order"]).T
    tv = []
    for row in p["spatial"]:
        left = conditional(row["reference_left"])[:-1]
        errors = []
        for d, observed in zip(p["separations"], row["pairs"], strict=True):
            right = np.roll(conditional(row["reference_right"])[:-1], d)
            pred = left[i] * right[j] + left[j] * right[i]
            errors.append(
                float(
                    0.5
                    * np.abs(conditional(observed) - np.r_[pred, 1 - pred.sum()]).sum()
                )
            )
        tv.append(errors)
    return np.array(tv)


def collect():
    return dict(
        starts=SHOTS,
        seed=970428,
        labels=LABELS,
        period=PERIOD,
        panels=[
            panel(p)
            for p in measure.read(measure.HERE.parent / "records.json")["panels"]
            if p["cells"] == 80
        ],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    dest = parser.parse_args().output
    data = collect()
    if dest:
        dest.write_text(json.dumps(data, indent=1, allow_nan=False) + "\n")
    else:
        assert json.loads(json.dumps(data)) == measure.read(measure.HERE / "pairs.json")
    for p in data["panels"]:
        print(
            p["n"],
            p["kappa"],
            [conditional(v)[:2].round(5).tolist() for v in p["energy"]],
            compare(p)[:, [0, -1]].round(4).tolist(),
        )
