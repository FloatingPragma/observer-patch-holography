"""Subsequent implementation check for the supplied sharp and coherent sources."""

import json
from pathlib import Path

import numpy as np
from . import apparatus as lab
from .measure import TARGETS


def build(state, exchange):
    target = state.ravel()
    out = np.zeros_like(target)
    out[0] = 1.0
    holds = 0
    # Split each remaining amplitude on the next physical neighbor. Resolve
    # very small tails by their remaining norms rather than accumulating roundoff.
    residual = np.sqrt(np.cumsum(abs(target[::-1]) ** 2))[::-1]
    for j in range(len(out) - 1):
        if residual[j] < 1e-14:
            break
        theta = np.arctan2(residual[j + 1], abs(target[j]))
        t = 2 * theta / exchange
        common = np.exp(-1j * (32 + exchange / 2) * t)
        a = out[j] * common
        out[j], out[j + 1] = a * np.cos(theta), 1j * a * np.sin(theta)
        if abs(target[j]) > 1e-14:
            dwell = (-np.angle(target[j] / out[j]) % (2 * np.pi)) / 32
            out[j] *= np.exp(-32j * dwell)
        holds += 1
    if abs(target[-1]) > 1e-14:
        dwell = (-np.angle(target[-1] / out[-1]) % (2 * np.pi)) / 32
        out[-1] *= np.exp(-32j * dwell)
    return holds, float(np.max(abs(target - out)))


def report():
    rows = []
    for lev in (1, 2):
        route = lab.mount(lev)
        board = lab.federation(lev)
        private = next(x for x in route if board.degree[x] == 5)
        n = len(route) // 2
        for ex in (1.0, 0.8):
            layers = lab.gates(ex, 8.0)
            error = 0.0
            number = 0
            for target in TARGETS:
                h = round(target * n / (2 * np.pi))
                for plane in (False, True):
                    a, b = lab.arms(n, h, layers, ex, 8.0, plane=plane)
                    for psi in [a, b] + [
                        (a + q * b) / np.sqrt(2) for q in (1, 1j, -1, -1j)
                    ]:
                        holds, e = build(psi, ex)
                        number += 1
                        error = max(e, error)
            rows.append(
                dict(
                    level=lev,
                    exchange=ex,
                    private_phase_port=private,
                    states=number,
                    maximum_adjacent_holds=len(route) - 1,
                    maximum_column_error=error,
                )
            )
    return rows


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, required=True)
    p.parse_args().output.write_text(
        json.dumps(report(), indent=1, allow_nan=False) + "\n"
    )
