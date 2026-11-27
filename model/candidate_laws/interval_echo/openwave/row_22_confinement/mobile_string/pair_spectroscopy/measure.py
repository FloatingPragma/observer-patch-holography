"""Native static two-wall spectroscopy, full clock and material exit banks."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, eigvalsh

from .. import material

HERE = Path(__file__).resolve().parent
SHOTS = 2**24
CLOCK = 8192
SPAN = 512
FILTER = 256
TIME = (0, 2, 16, 32, 128)


def populations(n):
    u = material.bits(n)
    end = np.column_stack((np.zeros(len(u), int), u, np.zeros(len(u), int)))
    signed = np.diff(end, axis=1)
    rising = (signed == 1).sum(axis=1)
    upper = u.sum(axis=1)
    return upper, rising


def words(n):
    site = n // 2 - 1
    return dict(
        host=0,
        one=1 << site,
        two=3 << site,
        four=15 << (site - 1),
        separated=(1 << 1) | (1 << (n - 2)),
    )


def aperture(energies, target):
    x = (energies - target) / SPAN
    return (
        np.exp(-1j * np.pi * (FILTER - 1) * x) * np.sinc(FILTER * x) / np.sinc(x)
    ) ** 4


def draw(p, rng):
    if p.min() < -1e-9 or abs(p.sum() - 1) > 1e-8:
        raise ArithmeticError("incomplete material exits")
    p = np.maximum(p, 0.0)
    return rng.multinomial(SHOTS, p / p.sum())


def raw_clock(e, v, sources):
    """FFT of the finite Fejer kernel gives all delay-pointer probabilities."""
    j = np.arange(CLOCK)
    weights = abs(v[list(sources.values())]) ** 2
    chi = np.empty((len(weights), CLOCK), complex)
    for i in range(0, CLOCK, 128):
        grid = 2j * np.pi * np.outer(e, j[i : i + 128]) / SPAN
        chi[:, i : i + 128] = weights @ np.exp(-grid)
    spectrum = np.fft.ifft(chi * (1 - j / CLOCK), axis=1).real * 2 - 1 / CLOCK
    return spectrum


def capture(c, v, prior, first):
    yes = abs(v @ c) ** 2
    return np.r_[1 - first, first - prior, prior - yes.sum(), yes]


def conversion(v, first, target, pulse, q, rng):
    start = first * v[q]
    passed = float(np.vdot(start, start).real)
    actual = v @ start
    if pulse is not None:
        actual = actual[np.arange(len(actual)) ^ (1 << pulse)]
    after = target * (v.T @ actual)
    return draw(capture(after, v, passed, 1.0), rng)


def measure(n, mode, original, rng, banks):
    g = 0.0 if mode == "degree_off" else float(original[2])
    h = 0.0 if mode == "exchange_off" else float(original[1] / 2)
    # Two frozen lower boundary dimers contribute a common scalar. The internal
    # sites of core(n+2,0,...) have exactly the physical signed end contacts.
    static = material.core(n + 2, 0, g, h)
    if h:
        e, v = eigh(static, driver="evd")
    else:
        order = np.argsort(np.diag(static), kind="stable")
        e, v = static.diagonal()[order], np.eye(len(static))[:, order]
    up, pairs = populations(n)
    name = f"{mode}_{n}"
    src = words(n)
    clocks = raw_clock(e, v, src)
    for key, p in zip(src, clocks, strict=True):
        banks[f"{name}_clock_{key}"] = draw(p, rng)

    band_result = {}
    for length, target in ((1, 40), (2, 48)):
        band = (e >= target - 4) & (e < target + 4)
        active = v[:, band]
        match = (up == length) & (pairs == 1)
        if active.shape[1]:
            gram = active[match].T @ active[match]
            enclosure = [float(e[band].min() - e[0]), float(e[band].max() - e[0])]
            survival = float(max(0.0, eigvalsh(gram)[0]))
        else:
            enclosure, survival = [], None
        band_result[str(length)] = dict(
            states=int(band.sum()),
            vacuum_relative_range=enclosure,
            single_pair_length_weight_lower=survival,
        )
    targets = {target: aperture(e, target) for target in (40, 48)}
    for length, target in ((1, 40), (2, 48)):
        seed = src["one" if length == 1 else "two"]
        incoming = targets[target] * v[seed]
        passed = float(np.vdot(incoming, incoming).real)
        for t in TIME:
            c = incoming * np.exp(-1j * e * t)
            banks[f"{name}_flight_{length}_{t}"] = draw(capture(c, v, passed, 1.0), rng)
        for ctrl, site in [("near", n // 2), ("far", n - 2), ("parked", None)]:
            banks[f"{name}_convert_{length}_{ctrl}"] = conversion(
                v, targets[target], targets[88 - target], site, seed, rng
            )

    far = (pairs == 1) & (up >= 4)
    separate = pairs >= 2

    def constrained(mask):
        block = static[np.ix_(mask, mask)]
        if h:
            return float(
                eigh(block, subset_by_index=[0, 0], eigvals_only=True)[0] - e[0]
            )
        return float(block.diagonal().min() - e[0])

    return dict(
        n=n,
        mode=mode,
        ground=float(e[0]),
        first_gap=float(e[1] - e[0]),
        two_separated_lowest_products=float(2 * (e[1] - e[0])),
        far_pair_conditioned_min=constrained(far),
        two_pair_conditioned_min=constrained(separate),
        band=band_result,
        input_words=src,
    )


def acquisition():
    rng = np.random.default_rng(181831)
    calibrated = material.calibration_energies()
    np.testing.assert_allclose(calibrated, [32, 1, 8], atol=1e-12)
    # Physical deployment and its checked source states are inherited from the
    # 16-dimer stock; the cut here holds dimer 0 and n+1 in the lower state.
    route = material.route(0)
    banks = {}
    panels = [
        measure(n, mode, calibrated, rng, banks)
        for n in (8, 10)
        for mode in ("native", "degree_off", "exchange_off")
    ]
    summary = dict(
        seed=181831,
        shots=SHOTS,
        clock_bins=CLOCK,
        energy_period=SPAN,
        clock_steps_per_filter=FILTER,
        physical_ports=route,
        calibration=calibrated.tolist(),
        panels=panels,
        terminal_order=[
            "input refusal (unused)",
            "first temporal refusal",
            "second temporal refusal",
            "all orientation words as nonnegative integers, least bit first",
        ],
    )
    return summary, banks


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    path = parser.parse_args().output
    data, banks = acquisition()
    path.with_suffix(".json").write_text(
        json.dumps(data, indent=1, allow_nan=False) + "\n"
    )
    np.savez_compressed(path.with_suffix(".npz"), **banks)
