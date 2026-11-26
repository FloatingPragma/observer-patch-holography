"""Bilateral elapsed-phase comparison through the native equilibrium heads."""
import argparse
import cmath
import json
import math
from pathlib import Path

import numpy as np

from ...full_algebra import native_precision as native
from ...full_algebra.reuse import clock_history
from ..measure import density_center

HERE = Path(__file__).resolve().parent
SHOTS = 2**38
GRID = 2**40
SEED = 805318
PANELS = ((0.5, 6), (1.0, 11))
# Two charges, directions, panels and apertures (prompt and delayed), two
# epochs, two components, local and remote. Distinct cohorts for every score.
SCORES = 2**8
SCORE_ERROR = math.sqrt(2 * math.log(2 * SCORES / 0.01) / SHOTS) + 7e-7


def arrivals():
    """First-hit probabilities in half-tick bins from the native integer propagator."""
    states, _, _, _ = density_center()
    board = native.iso.pairing.formation.federation(1)
    marks = native.iso.pairing.formation.region(board)
    ports = (marks[1], marks[2])
    assert ports == (1, 38)
    addresses = np.array([4 * port + i for port in ports for i in range(4)])
    weights = native.prior.weights()
    result = {}
    for charge in (-1, 1):
        probabilities = np.zeros((2, 11))
        for parity, state in enumerate(states):
            h = native.generator(board, marks, charge, bool(parity))
            panels = native.kernel_panels(h, addresses, weights)
            for sender in (0, 1):
                amplitudes = clock_history.first_hits(panels, sender, 1 - sender)[1:]
                for j, (real, imag) in enumerate(amplitudes):
                    a = (np.asarray(real, dtype=float) + 1j * np.asarray(imag, dtype=float)) / native.prior.FRAC
                    probabilities[sender, j] += np.trace(a @ state @ a.conj().T).real
        if np.min(probabilities) < 0 or np.max(probabilities.sum(axis=1)) > 1:
            raise ArithmeticError("invalid first-hit distribution")
        result[str(charge)] = probabilities.tolist()
    return result


def clock_bank(rng, epochs, probabilities):
    """Terminal sign for each quadrature, including a zero score on every miss."""
    phases = [cmath.exp(1j * t) for t in epochs]
    rows = []
    for quadrature in ("real", "imag"):
        cells = []
        for p, phase in zip(probabilities, phases, strict=True):
            score = getattr(phase, quadrature)
            plus = round(p * (1 + score) * GRID / 2)
            minus = round(p * (1 - score) * GRID / 2)
            cells.extend([plus, minus])
        cells.append(GRID - sum(cells))
        if min(cells) < 0:
            raise ArithmeticError("invalid clock-sign sampler")
        counts = rng.multinomial(SHOTS, np.array(cells) / GRID)
        rows.append(dict(by_hit=counts[:-1].reshape(-1, 2).tolist(), missed=int(counts[-1])))
    return rows


def acquire(distributions):
    rng = np.random.default_rng(SEED)
    panels = []
    for charge, direction_bins in distributions.items():
        for period, ticks in PANELS:
            for delay in (0, 2):
                heads = []
                for direction, bins in enumerate(direction_bins):
                    prob = bins[:ticks]
                    local, remote = [], []
                    for t in (0, period):
                        local.append(clock_bank(rng, [t], [1]))
                        epochs = [(1 + delay) * t + 0.5 * n for n in range(1, ticks + 1)]
                        remote.append(clock_bank(rng, epochs, prob))
                    heads.append(dict(sender=(1, 38)[direction], local=local, remote=remote))
                panels.append(dict(charge=int(charge), period=period, ticks=ticks,
                                   late_wait_multiple=delay, directions=heads))
    return panels


def phase_and_error(bank):
    scores = []
    for row in bank:
        bins = row["by_hit"]
        if not bins or any(len(pair) != 2 or min(pair) < 0 for pair in bins):
            raise ValueError("an incomplete Ramsey bank")
        total = sum(plus + minus for plus, minus in bins)
        if row["missed"] < 0 or row["missed"] + total != SHOTS:
            raise ValueError("an incomplete Ramsey bank")
        scores.append(sum(plus - minus for plus, minus in bins) / SHOTS)
    z = complex(*scores)
    radius = math.sqrt(2) * SCORE_ERROR
    if abs(z) <= radius:
        raise ValueError("no resolved Ramsey phase")
    return z, math.asin(radius / abs(z))


def interval(bank):
    (base, a), (later, b) = (phase_and_error(row) for row in bank)
    phase = cmath.phase(later * base.conjugate())
    if not 0 < phase - a - b or not phase + a + b < math.pi:
        raise ValueError("elapsed phase not within the registered unwrapped branch")
    return (phase - a - b, phase + a + b)


def infer(panel):
    # The local-to-received elapsed interval is the received-to-sent frequency.
    low = high = 1
    for arm in panel["directions"]:
        x = interval(arm["local"])
        y = interval(arm["remote"])
        low *= x[0] / y[1]
        high *= x[1] / y[0]
    return dict(product=[low, high],
                radial_beta=[(1 - high) / (1 + high), (1 - low) / (1 + low)],
                radial_rapidity=[-math.log(high) / 2, -math.log(low) / 2])


def replay():
    record = json.loads((HERE / "results.json").read_text())
    native_bins = arrivals()
    for charge, bins in native_bins.items():
        np.testing.assert_allclose(bins, record["first_hit"][charge], atol=2e-13, rtol=0)
    for panel in record["panels"]:
        for field, bounds in infer(panel).items():
            np.testing.assert_allclose(bounds, panel["inferred"][field], atol=1e-13, rtol=0)
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    if parser.parse_args().write:
        bins = arrivals()
        panels = acquire(bins)
        for panel in panels:
            panel["inferred"] = infer(panel)
        output = dict(seed=SEED, shots_per_quadrature=SHOTS, first_hit=bins, panels=panels)
        (HERE / "results.json").write_text(json.dumps(output, indent=2) + "\n")
    else:
        replay()
        print("Verified native transfer probabilities and both signed clock comparisons.")
