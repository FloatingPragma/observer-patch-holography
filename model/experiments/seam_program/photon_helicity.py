"""Resolve axial rotation characters of the full eight-slot one-load walk."""

import argparse
import json
from pathlib import Path

import numpy as np

from .gate import real_step


def quarter_turn(axis):
    """Right-handed quarter turn of the six spatial slots, fixing the private slots."""
    first, second = (axis + 1) % 3, (axis + 2) % 3
    target = np.arange(8)
    for sign in (0, 1):
        target[2 * first + sign] = 2 * second + sign
        target[2 * second + sign] = 2 * first + (1 - sign)
    return np.eye(8)[:, target]


def bloch_columns(axis, order, n):
    """Send each socket through the literal spatial update and read its Fourier return."""
    shape = [1] * 4
    shape[axis] = n
    wave = (np.exp(2j * np.pi * order * np.arange(n) / n) / np.sqrt(n)).reshape(shape)
    columns = [real_step(wave[..., None] * slot) for slot in np.eye(8)]
    return np.stack([wave.ravel().conj() @ column.reshape(n, 8) for column in columns], axis=1)


def run(n=32):
    if n < 4 or n % 2:
        raise ValueError("Expected an even axial period of at least four")
    identity = np.eye(8)
    rows = []
    for axis in range(3):
        rotation = quarter_turn(axis)
        transverse = [
            sum(np.linalg.matrix_power(rotation, j) / phase**j for j in range(4)) / 4
            for phase in (1j, -1j)
        ]
        for order in range(1, n // 2):
            walk = bloch_columns(axis, order, n)
            values, vectors = np.linalg.eig(walk)
            mobile = (abs(values - 1) > 1e-8) & (abs(values + 1) > 1e-8)
            moving = vectors[:, mobile]
            rows.append(dict(
                axis=axis, order=order, moving_rank=int(mobile.sum()),
                moving_phases=sorted(float(v) for v in np.angle(values[mobile])),
                moving_quarter_turn_error=float(np.linalg.norm((rotation - identity) @ moving)),
                transverse_ranks=[int(np.linalg.matrix_rank(p)) for p in transverse],
                transverse_unit_step_errors=[float(np.linalg.norm((walk - identity) @ p)) for p in transverse],
            ))
    return dict(period=n, private_order=0, rows=rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.parse_args().output.write_text(json.dumps(run(), indent=2, allow_nan=False) + "\n")
