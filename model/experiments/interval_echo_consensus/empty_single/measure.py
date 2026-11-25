"""The C/E coherence observable on its declared twelve-port source and reader."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from ..empty_full import historical as cell

HERE = Path(__file__).resolve().parent
LAW = cell.m.LAW
Z = cell.WORDS.index((1, 2))
PHASES = np.pi * np.array([0, 0.5, 1, 1.5, 1 / 3, -1 / 3, 0.75])
SHOTS = 2**18


def require(ok, message):
    if not ok:
        raise ValueError(message)


def channels():
    e = LAW.echo()[:, :9].reshape(-1, 9, 9)
    c = np.zeros((2, 3, 9, 9))
    for i, (a, b) in enumerate(LAW.SRC):
        for coin in (0, 1):
            out = LAW.native_event(a, b, coin)
            current = out[0] - a
            c[coin, 1 + current, LAW.SRC.index(out), i] = 1 / np.sqrt(2)
    c = c.reshape(-1, 9, 9)
    require(np.allclose((c**2).sum(axis=0), (e**2).sum(axis=0), atol=1e-14, rtol=0),
            "control has different sharp transition probabilities")
    instruments = {"E_native": e, "C_control": c, "bypass": np.eye(9)[None]}
    for blocks in instruments.values():
        require(np.max(abs(sum(k.conj().T @ k for k in blocks) - np.eye(9))) < 1e-14,
                "incomplete local record family")
    return instruments


def run():
    h = LAW.seam_term(cell.m.candidate.KAPPA)
    route = cell.lift(expm(-1j * np.pi * h), (3, 8))
    mixer = cell.lift(expm(-0.5j * np.pi * h), (1, 8))
    source = cell.phase(1, -np.pi / 2) @ route @ mixer @ route
    state = source[:, cell.X]
    target = np.zeros_like(state)
    target[[cell.X, Z]] = 1 / np.sqrt(2)
    require(abs(abs(np.vdot(target, state)) - 1) < 1e-13, "wrong coherent source")
    rng = np.random.Generator(np.random.PCG64(138507))
    instruments = channels()
    family = len(cell.WORDS) * len(PHASES) * len(instruments)
    radius = np.sqrt(np.log(2 * family / 0.01) / (2 * SHOTS))
    banks = {}
    for label, blocks in instruments.items():
        operators = [cell.lift(k, (0, 1)) for k in blocks]
        empty, single = LAW.SRC.index((0, 0)), LAW.SRC.index((0, 1))
        g = np.vdot(blocks[:, single, single], blocks[:, empty, empty]).real
        offset = (np.linalg.norm(blocks[:, empty, empty])**2 +
                  np.linalg.norm(blocks[:, single, single])**2) / 4
        expected = {"E_native": (0, 3/8), "C_control": (0.5, 3/8), "bypass": (1, 0.5)}
        require(np.max(abs(np.array([g, offset]) - expected[label])) < 2e-14,
                "full native record disagrees with the sharp-key Gram entry")
        panel, frequencies = [], []
        for theta in PHASES:
            prepared = cell.phase(1, theta) @ state
            public = sum(abs(source.conj().T @ k @ prepared)**2 for k in operators)
            require(abs(public.sum() - 1) < 4e-13, "lost terminal mass")
            theory = offset + g * np.cos(theta) / 2
            require(abs(public[cell.X] - theory) < 2e-13, "native bright return differs")
            counts = rng.multinomial(SHOTS, public / public.sum())
            require(np.max(abs(counts / SHOTS - public)) <= radius, "held counts outside tolerance")
            frequencies.append(counts[cell.X] / SHOTS)
            panel.append({"phase_in_pi": round(theta / np.pi, 12), "counts": counts.tolist()})
        a, b, c, d = frequencies[:4]
        center, real, imag = (a + b + c + d) / 4, a - c, b - d
        for j, theta in enumerate(PHASES[4:], start=4):
            pred = center + (real * np.cos(theta) + imag * np.sin(theta)) / 2
            bound = radius * (2 + abs(np.cos(theta)) + abs(np.sin(theta)))
            require(abs(pred - frequencies[j]) < bound, "quadrature forecast failed")
            panel[j]["bright_forecast"] = round(pred, 9)
        banks[label] = {"gram_multiplier": round(g, 12), "offset": round(offset, 12),
                        "inferred_real_imag": [round(real, 9), round(imag, 9)], "panels": panel}
    return {"seam": [0, 1], "word_order": [list(w) for w in cell.WORDS],
            "bright_word": [2, 3], "source_partner": [1, 2], "shots_per_setting": SHOTS,
            "seed_pcg64": 138507, "original_starts": SHOTS * len(instruments) * len(PHASES),
            "simultaneous_99_percent_radius": round(float(radius), 12), "banks": banks}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    text = json.dumps(run(), indent=1) + "\n"
    path = HERE / "results.json"
    if args.write:
        path.write_text(text)
    require(text == path.read_text(), "held receipt differs")
    print("Verified all 78 terminal words under the native event, record control and bypass.")
