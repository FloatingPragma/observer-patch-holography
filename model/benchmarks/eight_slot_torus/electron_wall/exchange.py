"""Two noninteracting exchange completions of precisely the same one-body walk."""
import argparse
import json

import numpy as np

from .walk import links, prepare, step, wall_mode, wave


def site_pairs(first, second):
    """Unordered position-pair records. Internal sockets are traced at the detector."""
    if first.shape != second.shape or first.ndim != 7:
        raise ValueError("both exchange probes must share the wall carrier")
    if not np.isfinite(first).all() or not np.isfinite(second).all():
        raise ValueError("exchange probes must have finite amplitudes")
    sites = int(np.prod(first.shape[1:4]))
    first = np.moveaxis(first, 0, 3).reshape(sites, -1)
    second = np.moveaxis(second, 0, 3).reshape(sites, -1)
    overlap = np.sum(first.conj() * second, axis=1)
    a = np.sum(abs(first) ** 2, axis=1)
    b = np.sum(abs(second) ** 2, axis=1)
    if max(abs(a.sum() - 1), abs(b.sum() - 1), abs(overlap.sum())) > 1e-10:
        raise ValueError("exchange probes must be normalized and orthogonal")
    direct = (np.outer(a, b) + np.outer(b, a)) / 2
    cross = np.outer(overlap, overlap.conj()).real
    detectors = {}
    for name, sign in (("symmetric", 1), ("tagged", 0), ("antisymmetric", -1)):
        ordered = direct + sign * cross
        unordered = np.diag(np.diag(ordered)) + np.triu(ordered + ordered.T, 1)
        if unordered.min() < -1e-12 or abs(unordered.sum() - 1) > 1e-12:
            raise ArithmeticError("invalid pair probabilities")
        detectors[name] = unordered
    return detectors


def histories(period=13, length=3, waits=(0, 1, 4, 16)):
    shape = (period, 1, 1)
    conn = links(shape)
    modes = []
    for sign in (-1, 1):
        p = sign * 2 * np.pi / period
        vector, _ = wall_mode(length, p, 0)
        modes.append(prepare(wave(shape, p, vector)))
    records = {}
    for n in range(max(waits) + 1):
        if n in waits:
            records[n] = site_pairs(*modes)
        if n < max(waits):
            modes = [step(state, conn, -1) for state in modes]
    return records


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--details", action="store_true", help="emit every unordered detector pair")
    args = parser.parse_args()
    print(json.dumps({n: {kind: matrix.tolist() if args.details else float(np.trace(matrix))
                         for kind, matrix in readings.items()}
                      for n, readings in histories().items()}, indent=2))
