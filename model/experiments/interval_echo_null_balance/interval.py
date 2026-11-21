"""Finite counts in the patch energy meter over a held free-flight interval."""

import argparse
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh

from model.executable_model import cell
from .experiment import ANGLES, HERE, KAPPA, WEIGHTS, hardware, portable, states

SHOTS = 2**19
SEED = 114702
START, END = 1.125, 1.375


class Counter:
    def __init__(self, sector, meter):
        masks = np.zeros(sector.P, dtype=int)
        for index in (sector.aa, sector.bb):
            inside = index < 12
            masks[inside] |= 1 << index[inside]
        self.words, self.word_index = np.unique(masks, return_inverse=True)
        self.output_words = []
        self.clicks = []
        self.meter = meter
        # A single diagonal port-word reading counts every seam mismatch at once.
        self.scores = KAPPA * np.array([
            sum(((word >> a) ^ (word >> b)) & 1 for a, b in cell().seams)
            for word in self.words
        ])
        self.rng = np.random.default_rng(SEED)

    def capture(self, state):
        occupancy = np.bincount(self.word_index, weights=abs(state)**2, minlength=len(self.words))
        words = self.rng.multinomial(SHOTS, occupancy / occupancy.sum())
        x, y = state[self.meter.left], state[self.meter.right]
        probability = (.5 * abs(x - y)**2).sum(axis=1)
        clicks = self.rng.binomial(SHOTS, probability)
        self.output_words.append(words)
        self.clicks.append(clicks)
        weights = words / SHOTS
        p = clicks / SHOTS
        mean = weights @ self.scores
        variance = weights @ (self.scores - mean)**2 + np.sum(p * (1-p))
        true_energy = occupancy @ self.scores + np.sum(probability)
        if abs(true_energy - self.meter.read(state)) > 1e-10:
            raise ArithmeticError("the counted words and seam settings lost energy")
        return float(mean + np.sum(p)), float(variance / (SHOTS-1))


def appointments(order):
    points, weights = leggauss(order)
    return (START + END + (END - START) * points) / 2, (END - START) * weights / 2


def measure():
    sector, _, _, meter, cell_h, bonds = hardware()
    energy, vectors = eigh((sum(cell_h) + sum(b.h for b in bonds)).toarray())
    source = next(s for p, sign, s in states(sector) if p == 9 and sign == 1)
    initial = vectors.T @ source
    def flight(t):
        return vectors @ (np.exp(-1j * t * energy) * initial)

    counter = Counter(sector, meter)
    start, end = counter.capture(flight(START)), counter.capture(flight(END))
    delta = np.array([end[0] - start[0], end[1] + start[1]])

    def integrate(order, sample):
        total = np.zeros((2, 2))
        times, weights = appointments(order)
        for t, weight in zip(times, weights, strict=True):
            state = flight(t)
            for j, bond in enumerate(bonds):
                for angle, derivative in zip(ANGLES, WEIGHTS, strict=True):
                    wave = bond.expose(state, angle)
                    mean, variance = counter.capture(wave) if sample else (meter.read(wave), 0.)
                    gain = weight * derivative
                    total[j] += [gain * mean, gain**2 * variance]
        return total

    actual = integrate(6, True)
    low, high = integrate(6, False)[:, 0], integrate(12, False)[:, 0]
    difference = meter.read(flight(END)) - meter.read(flight(START))
    def describe(value, variance):
        se = np.sqrt(variance)
        return dict(value=float(value), standard_error=float(se), excludes_zero_at_4se=bool(abs(value) > 4 * se))

    report = dict(interval=[START, END], quadrature_points=appointments(6)[0].tolist(),
                  current_port=9, current_sign=1, preparations_per_setting=SHOTS, seed=SEED,
                  energy_readouts=len(counter.clicks), total_preparations=len(counter.clicks)*31*SHOTS,
                  expected_change=difference, expected_gate_integrals=high.tolist(),
                  largest_change_doubling_quadrature=float(max(abs(low-high))),
                  expected_residual_6=float(difference - sum(low)),
                  endpoint_change=describe(*delta), gateways=[describe(*x) for x in actual],
                  endpoint_minus_both=describe(delta[0]-actual[:,0].sum(), delta[1]+actual[:,1].sum()),
                  endpoint_minus_first=describe(delta[0]-actual[0,0], delta[1]+actual[0,1]))
    archives = dict(local_word_labels=counter.words, diagonal_word_counts=np.array(counter.output_words),
                    antisymmetric_click_counts=np.array(counter.clicks))
    return portable(report), archives


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    report, archives = measure()
    payload = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "interval.json").write_text(payload)
        np.savez_compressed(HERE / "counts.npz", **archives)
    else:
        if json.loads((HERE / "interval.json").read_text()) != report:
            raise ValueError("finite-sample estimates differ from the saved receipt")
        with np.load(HERE / "counts.npz") as saved:
            for name, value in archives.items():
                np.testing.assert_array_equal(value, saved[name])
    print(payload)


if __name__ == "__main__":
    main()
