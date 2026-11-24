"""Run actual pulsed native contacts at and beside the critical duty."""

import argparse
import json

import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.executable_model.carrier import federation
from . import experiment as ex

DELTA = 1 / 8
CYCLES = 2 ** np.arange(11)


def rotate_pairs(matrix, left, angle):
    right = left + 1
    a, b = matrix[left].copy(), matrix[right].copy()
    matrix[left] = np.cos(angle) * a + np.sin(angle) * b
    matrix[right] = -np.sin(angle) * a + np.cos(angle) * b


def parity_phase(n, j, h, f):
    # The open wall chain has n+1 Jordan-Wigner modes, two Majoranas apiece.
    # Similarity by staggered Majorana signs can reverse either layer.
    width = 2 * (n + 1)
    transform = np.eye(width)
    rotate_pairs(transform, np.arange(1, width - 1, 2), 2 * h * DELTA)
    rotate_pairs(transform, np.arange(0, width, 2), 2 * j * f * DELTA)
    phases = np.angle(np.linalg.eigvals(transform))
    # The two signs have equal magnitude, including at an unresolved quiet edge.
    return float(min(abs(phases)))


def run():
    rng = np.random.default_rng(201125)
    static_rng = np.random.default_rng(201126)
    board = federation(2)
    _, path, _ = ex.guide(board)
    j, h = ex.couplings(candidate.KAPPA)
    # Local gap: symmetric/antisymmetric exchange pair; unequal/equal Zeno contact.
    cal = {
        name: ex.phase_reader(value, rng)
        for name, value in (("two_J", 2 * j), ("two_h", 2 * h))
    }
    ratio = h / j
    schedules = []
    for factor in (0.5, 1.0, 2.0):
        f = factor * ratio
        period = DELTA * (1 + f)
        rows = []
        for n in range(2, (len(path) - 2) // 2 + 1):
            sign = 1 if n % 2 == 0 else -1
            phase = sign * parity_phase(n, j, h, f)
            read = ex.phase_reader(phase, rng, CYCLES)
            rows.append(
                dict(
                    n=n,
                    phase=phase,
                    energy=phase / period,
                    counts=read,
                    measured_energy=read["energy"] / period,
                    simultaneous_error=read["simultaneous_error"] / period,
                )
            )
        schedules.append(dict(factor=factor, duty=f, actual_period=period, rows=rows))
    static_odd = []
    for n in range(3, (len(path) - 2) // 2 + 1, 2):
        val = -ex.splitting(n, j, h)
        static_odd.append(dict(n=n, exact=val, read=ex.phase_reader(val, static_rng)))
    return dict(
        calibration=cal,
        calibrated_duty=cal["two_h"]["energy"] / cal["two_J"]["energy"],
        selected_duty=ratio,
        cycles=CYCLES.tolist(),
        schedules=schedules,
        registered_static_odd=static_odd,
        shots_per_quadrature=ex.SHOTS,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = ex.portable(run())
    path = ex.HERE / "schedule.json"
    if args.write:
        path.write_text(json.dumps(data, indent=1, allow_nan=False) + "\n")
    elif json.loads(path.read_text()) != data:
        raise ValueError("Pulse-schedule replay differs from the receipt")
    print(
        "Verified signed boundary response for all scheduled and odd-length controls."
    )


if __name__ == "__main__":
    main()
