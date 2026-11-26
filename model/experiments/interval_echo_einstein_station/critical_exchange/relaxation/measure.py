"""Bound local return after an arbitrary charge-preserving two-station disturbance."""

import argparse
import json
from fractions import Fraction
from math import factorial, isqrt
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo.openwave.row_30.experiment import select_route
from model.executable_model.carrier import federation
from ..measure import ring_modes

HERE = Path(__file__).resolve().parent
TIMES = (0, 8, 16, 32, 64, 80, 96)
END = 96


def square_root_upper(q):
    scale = 10**15
    return Fraction(isqrt(q.numerator * scale**2 // q.denominator) + 1, scale)


def bessel_abs_upper(order, t):
    """Exact rational alternating-series enclosure at the specified integer."""
    half = Fraction(t, 2)
    term = half**order / factorial(order)
    total = term
    m = 0
    while True:
        m += 1
        term *= -half**2 / (m * (m + order))
        if m * (m + order) > half**2 and abs(term) < Fraction(1, 10**30):
            return max(abs(total), abs(total + term))
        total += term


def window(n, start):
    # Off-diagonal endpoints differ by at most one station.
    first = n - 1
    half_end = Fraction(END, 2)
    tail = 2 * half_end**first / factorial(first) / (1 - half_end / (first + 1))
    envelope = square_root_upper(sum(bessel_abs_upper(k, start)**2 for k in (0, 1)))
    bound = 2 * (envelope + 2 * tail)
    # Round the displayed upper bound outwards at one part in a billion.
    units = (bound.numerator * 10**9 + bound.denominator - 1) // bound.denominator
    return {
        "start": start, "end": END,
        "universal_bond_error_upper": units / 10**9,
        "propagator_winding_error_upper": float(tail),
    }


def run():
    route, _ = select_route(federation(2))
    n = len(route)
    k, basis, filled = ring_modes(n)
    base = basis[:, filled] @ basis[0, filled].conj()
    reference_energy = -base[1].real
    h = np.array([[0, -0.5], [-0.5, 0]])
    reference = basis[:2, filled] @ basis[:2, filled].conj().T
    rng = np.random.default_rng(253296)
    shots = 2**18

    def terminal_counts(c):
        double = np.linalg.det(c).real
        one = np.trace(c).real - 2 * double
        coherence = c[0, 1].real
        return rng.multinomial(shots, [1 - one - double, one / 2 + coherence, one / 2 - coherence, double]).tolist()

    times = []
    counts = {"stationary_reference": terminal_counts(reference)}
    for t in TIMES:
        phase = np.exp(1j * t * np.cos(k))
        column = basis[:2] @ (phase * basis[0].conj())
        centred = basis[:2] @ (phase * (filled.astype(float) - 0.5) * basis[0].conj())
        coefficient = np.vdot(centred, h @ column)
        pi_change = -4 * coefficient.real
        largest = 2 * (abs(coefficient.real) + abs(coefficient))
        times.append({
            "time": t, "reference_energy": reference_energy,
            "phase_flip_energy": reference_energy + pi_change,
            "phase_flip_change": pi_change,
            "maximum_change_over_one_station_channels": largest,
        })
        if t in (0, END):
            changed = reference - 2 * (np.outer(column, centred.conj()) + np.outer(centred, column.conj()))
            counts[f"phase_flip_time_{t}"] = terminal_counts(changed)
    return {
        "stations": n, "disturbance_and_receiver_ports": route[:2],
        "one_station_phase_flip_port": route[0],
        "half_filling": int(filled.sum()),
        "stationary_uniform_filling_energy": 0,
        "measurements": times, "certified_windows": [window(n, t) for t in (16, 32, 64)],
        "synthetic_terminal_counts": counts, "shots_per_setting": shots,
        "count_order": ["empty", "symmetric", "antisymmetric", "double"],
        "seed": 253296,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = run()
    for read in data["measurements"]:
        for name in ("reference_energy", "phase_flip_energy", "phase_flip_change", "maximum_change_over_one_station_channels"):
            read[name] = round(float(read[name]), 12)
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    elif data != json.loads(path.read_text()):
        raise ArithmeticError("Local return probabilities or certified envelopes differ")
    print("Verified the signed bond readout and all pre-return error bounds.")


if __name__ == "__main__":
    main()
