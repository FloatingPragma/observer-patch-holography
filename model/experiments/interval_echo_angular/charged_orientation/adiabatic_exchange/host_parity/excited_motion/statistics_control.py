"""Replay the commissioned Y pump with ordinary or antisymmetric vacancy hopping."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from . import pump


def hopping(count, antisymmetric):
    h, _, _ = pump.source_geometry(count)
    h = h.copy()
    if antisymmetric:
        # Exterior-power signs on ordered vacancy ports. All diagonal seam energies survive.
        words = pump.stocks(count)
        for i, j in zip(*np.nonzero(np.triu(h, 1)), strict=True):
            a, b = (p for p in range(4) if words[i][p] != words[j][p])
            spectators_between = words[i][a + 1:b].count(0)
            h[i, j] *= (-1) ** spectators_between
            h[j, i] = h[i, j]
    return h


def ordered_product(steps):
    """Last pulse acts on the left. Pairwise reduction retains chronological order."""
    while len(steps) > 1:
        n = len(steps) // 2 * 2
        pairs = steps[1:n:2] @ steps[:n:2]
        steps = np.concatenate((pairs, steps[n:])) if n < len(steps) else pairs
    return steps[0]


def read_exits(forward, reference, source):
    a, b = forward @ source, reference @ source
    dark = float(np.linalg.norm(a - b) ** 2 / 4)
    bright = float(np.linalg.norm(a + b) ** 2 / 4)
    if not np.isclose(dark + bright, 1, rtol=0, atol=8e-8):
        raise ArithmeticError("Balanced exits failed normalization")
    return dark


def run_case(detuning, dt, appointments, duration, antisymmetric):
    _, pairs, central = pump.source_geometry(2)
    flight = expm(-1j * dt * hopping(2, antisymmetric))
    count, rest = divmod(duration, appointments * dt)
    if rest or appointments % 12 or count < 1:
        raise ValueError("Use twelve equal address windows and whole microsteps")
    source = np.eye(6)[:, pump.index_pair(pump.PAIRS[0])]
    forward = np.eye(6, dtype=complex)
    addresses = [pump.index_pair(pair) for pair in pump.PAIRS]
    home_reads = [1.0]
    crossing_reads = []
    frame_size = appointments // 12
    for frame in range(1, 13):
        for first in range((frame - 1) * frame_size, frame * frame_size, 512):
            ids = np.arange(first, min(first + 512, frame * frame_size))
            angles = (ids + 0.5) * (2 * np.pi / appointments)
            aperture = 2 * (1 + np.cos(angles[:, None] - pump.PHI))
            diagonal = detuning * central - aperture @ pairs.T
            kick = np.exp(-0.5j * dt * diagonal)
            steps = kick[:, :, None] * flight * kick[:, None, :]
            forward = ordered_product(np.linalg.matrix_power(steps, int(count))) @ forward
        if frame == 6:
            # Each split microstep is symmetric. Transpose exactly reverses the
            # pulse order, retaining every dynamical phase on the reference arm.
            reference = forward.T @ forward
        probs = abs(forward @ source) ** 2
        if frame % 4 == 0:
            home_reads.append(float(probs[addresses[(frame // 4) % 3]]))
        if frame % 4 == 2:
            leg = frame // 4
            crossing_reads.append(float(probs[addresses[leg]] + probs[addresses[(leg + 1) % 3]]))
    return {
        "detuning": detuning, "microstep": dt, "appointments": appointments, "duration": duration,
        "statistics": "antisymmetric" if antisymmetric else "native",
        "dark_once": read_exits(forward, reference, source),
        "dark_twice": read_exits(forward @ forward, reference @ reference, source),
        "pinned_probabilities": home_reads, "crossing_probabilities": crossing_reads,
    }


def direct_route(antisymmetric):
    """Use the parent's actual addressed exchange and retracing routes on the sharp stock."""
    amplitudes = []
    for route in ((1, 3, 2, 1, 3, 2), (1, 3, 2, 2, 3, 1)):
        holes = {1, 2}
        amp = 1
        for leaf in route:
            if (0 in holes) != (leaf in holes):
                if antisymmetric:
                    amp *= (-1) ** sum(0 < p < leaf for p in holes)
                holes ^= {0, leaf}
        if holes != {1, 2}:
            raise ArithmeticError("Direct route did not restore the original sharp pair")
        amplitudes.append(amp)
    return (amplitudes[0] - amplitudes[1]) ** 2 / 4


def run():
    certificate = pump.route_certificate()
    return {
        "ports": certificate["ports"], "initial_stock": certificate["sharp_stock"],
        "direct_sharp_dark": {"native": direct_route(False), "antisymmetric": direct_route(True)},
        "pumps": [
            run_case(g, dt, n, 393216, anti)
            for dt, n in ((1 / 16, 49152), (1 / 32, 98304))
            for g in (16, -16) for anti in (False, True)
        ],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    dest = parser.parse_args().output
    result = json.dumps(run(), indent=2) + "\n"
    if dest is None:
        print(result, end="")
    else:
        dest.write_text(result)
