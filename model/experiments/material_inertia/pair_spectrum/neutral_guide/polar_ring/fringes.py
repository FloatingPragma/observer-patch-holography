"""Measure dispersed one-island fringes and host-relative energy on sharp input stocks."""

import argparse
import json
from pathlib import Path

import numpy as np

from .measure import HOLDS, SHOTS, spectrum
from .source import prepared
from .stock import charts

TIMES = [0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
OFFSET = 117 / 4
SEED = 191919


def counts(z, rng):
    probabilities = np.clip([(1 + z.real) / 2, (1 - z.imag) / 2], 0, 1)
    return rng.binomial(SHOTS, probabilities).tolist()


def amplitude(pair):
    return complex(2 * pair[0] / SHOTS - 1, 1 - 2 * pair[1] / SHOTS)


def unwrap(times, signals):
    estimate = None
    for t, z in zip(times, signals, strict=True):
        phase = np.angle(z)
        if estimate is not None:
            phase += 2 * np.pi * round((estimate * t - phase) / (2 * np.pi))
        estimate = float(phase / t)
    return estimate


def recoil_from_words(table, n):
    sites = 1 << np.arange(n)
    harmonics = np.fft.fft(table[:, sites] / SHOTS, axis=1)
    k = 1 + int(np.argmax(abs(harmonics[0, 1 : n // 2 + 1])))
    signal = harmonics[:, k]
    return dict(
        mode=k,
        initial_amplitude=float(abs(signal[0])),
        min_amplitude=float(min(abs(signal))),
        max_nonconstant=float(abs(harmonics[:, 1:]).max()),
        signed_recoil=unwrap(TIMES[1:], signal[0] * signal[1:].conj()),
    )


def run_one(n, kappa, cells, store):
    ground, blocks = spectrum(n, kappa)
    addr = np.arange(n)
    sharp = np.eye(1 << n, 1)[:, 0]
    waves = np.array(
        [prepared(np.exp(2j * np.pi * k * addr / n) / np.sqrt(n)) for k in (0, 1, -1)]
    )
    starts = np.stack(
        [
            waves[0],
            (waves[0] + waves[1]) / np.sqrt(2),
            (waves[0] + waves[2]) / np.sqrt(2),
        ],
        axis=1,
    )
    # Verify the same Fourier-sum stocks by a single amplitude-distribution circuit.
    for s, sign in ((1, 1), (2, -1)):
        profile = (1 + np.exp(2j * np.pi * sign * addr / n)) / np.sqrt(2 * n)
        np.testing.assert_allclose(starts[:, s], prepared(profile), atol=3e-14)
    rng = np.random.default_rng([SEED, n, int(kappa), cells])
    rng_clock = np.random.default_rng([SEED + 1, n, int(kappa), cells])
    panels = [np.zeros_like(starts) for t in TIMES]
    sharp_return = np.zeros(len(HOLDS), complex)
    wave_returns = np.zeros((3, len(HOLDS)), complex)
    for trans, e, v, local, line in blocks:
        weights = abs(v.conj().T @ (trans.conj().T @ sharp)) ** 2
        ew = np.exp(-1j * np.outer(HOLDS, e + OFFSET))
        sharp_return += ew @ weights
        for i, wave in enumerate(waves):
            w = abs(v.conj().T @ (trans.conj().T @ wave)) ** 2
            wave_returns[i] += ew @ w
        coeff = v.conj().T @ (trans.conj().T @ starts)
        for j, t in enumerate(TIMES):
            panels[j] += trans @ (v @ (np.exp(-1j * e * t)[:, None] * coeff))
    host_counts = [counts(z, rng_clock) for z in sharp_return]
    clocks = []
    for signal in wave_returns:
        bank = [counts(z, rng_clock) for z in signal]
        er = unwrap(
            HOLDS,
            [
                amplitude(h) * amplitude(r).conjugate()
                for h, r in zip(host_counts, bank, strict=True)
            ],
        )
        clocks.append(dict(plus_counts=bank, host_gap=er))
    probs = abs(np.array(panels).transpose(2, 0, 1)) ** 2
    np.testing.assert_allclose(probs.sum(axis=2), 1, atol=5e-13)
    reads = []
    for i, series in enumerate(probs):
        table = np.array([rng.multinomial(SHOTS, p / p.sum()) for p in series])
        key = f"b{cells}_n{n}_g{int(kappa)}_s{i}"
        store[key] = table
        reads.append(
            dict(
                source=["uniform", "forward", "reverse"][i],
                counts_key=key,
                **recoil_from_words(table, n),
            )
        )
    return dict(
        cells=cells,
        n=n,
        kappa=kappa,
        ground_line_weight=float(abs(ground[0]) ** 2),
        sharp_lower_return=host_counts,
        clock_order=[0, 1, -1],
        clocks=clocks,
        fringes=reads,
    )


def run(output):
    store = {}
    panels = []
    for c in charts():
        if c["stock"]:
            for kappa in (8, 4):
                panels.append(
                    run_one(len(c["stock"]["ring"]), kappa, c["cells"], store)
                )
    np.savez_compressed(output / "fringe_counts.npz", **store)
    data = dict(
        shots=SHOTS,
        seed=SEED,
        common_bypass_offset=OFFSET,
        fringe_times=TIMES,
        clock_holds=HOLDS,
        panels=panels,
    )
    (output / "fringes.json").write_text(
        json.dumps(data, indent=1, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args().output)
