"""Reconstruct a common state from three sharp reads and predict untrained unsharp events."""

import argparse
import json
from pathlib import Path
import numpy as np

from . import circuit as c
from . import measure as m

HERE = Path(__file__).parent
PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]])
DIRECTIONS = np.array([[3 / 5, 0, 4 / 5], [0, 1 / np.sqrt(2), -1 / np.sqrt(2)]])
SHOTS = 2**12 * np.array([111, 68])


def isometry(effects, rotate):
    blocks = []
    for j, e in enumerate(effects):
        v, u = np.linalg.eigh(e)
        if min(v) < -1e-14:
            raise ValueError("Nonpositive effect")
        block = (u * np.sqrt(v.clip(0))) @ u.conj().T
        if rotate:
            angle = (j + 1) * np.pi / 6
            block = (
                np.array(
                    [[np.cos(angle), np.sin(angle)], [-np.sin(angle), np.cos(angle)]]
                )
                @ block
            )
        blocks.append(block)
    return np.vstack(blocks)


def programs():
    unit = np.eye(2, dtype=complex)
    a, b = [(unit + np.tensordot(n, PAULI, axes=1)) / 4 for n in DIRECTIONS]
    return {
        "X": c.word_rotation(0, np.pi / 4, np.pi),
        "Y": c.word_rotation(0, np.pi / 4, -np.pi / 2),
        "Z": [],
        "fine": c.compile_reader(isometry((a, b, unit - a - b), False)),
        "direct": c.compile_reader(isometry((a + b, unit - a - b, 0 * unit), True)),
    }


def readback(banks):
    if set(banks) != {"X", "Y", "Z", "fine", "direct"}:
        raise ValueError("Incomplete acquisition")
    pooled = {}
    for key, bank in banks.items():
        n = np.array(bank)
        if (
            n.shape != (2, 7)
            or n.dtype.kind not in "iu"
            or n.min() < 0
            or not np.array_equal(n.sum(axis=1), SHOTS)
        ):
            raise ValueError("Invalid full source cohorts")
        pooled[key] = n.sum(axis=0) / SHOTS.sum()
    radius = float(np.sqrt(np.log(2 * 256 / 0.01) / (2 * SHOTS.sum())))
    # Other words and failures score zero in these signed tests.
    r = np.array([pooled[axis][0] - pooled[axis][1] for axis in ("X", "Y", "Z")])
    extra = 2 * radius
    norm = float(np.linalg.norm(r))
    observed = [
        pooled["fine"][:2].sum(),
        pooled["fine"][2:4].sum(),
        pooled["fine"][:4].sum(),
        pooled["direct"][:2].sum(),
    ]
    dirs = [*DIRECTIONS, DIRECTIONS.sum(axis=0), DIRECTIONS.sum(axis=0)]
    forecasts = []
    for index, (n, p) in enumerate(zip(dirs, observed, strict=True)):
        predicted = (1 + int(index >= 2) + r @ n) / 4
        allowance = float(radius + extra * abs(n).sum() / 4)
        forecasts.append(
            {
                "predicted": round(float(predicted), 9),
                "observed": round(float(p), 9),
                "allowance": round(allowance, 9),
                "passes": bool(abs(p - predicted) <= allowance),
            }
        )
    return {
        "estimated_bloch": np.round(r, 9).tolist(),
        "coordinate_radius": round(extra, 9),
        "positive_at_99": bool(norm + np.sqrt(3) * extra < 1),
        "lower_eigenvalue_99": round((1 - norm - np.sqrt(3) * extra) / 2, 9),
        "forecasts": forecasts,
    }


def run():
    rng = np.random.default_rng(113031209)
    result = {}
    for k in (c.K, 0):
        sharp = np.eye(729, dtype=complex)[:, c.INDEX[:2]]
        result[str(k)] = {}
        for source, v in [
            ("original", sharp),
            ("rotated", c.play(c.word_rotation(0, np.pi / 5, np.pi / 7), sharp, k)),
        ]:
            bank = {}
            for key, pulses in programs().items():
                out = c.play(pulses, v, k)
                bank[key] = [
                    m.counts(rng, n, c.read(out[:, j])).tolist()
                    for j, n in enumerate(SHOTS)
                ]
            result[str(k)][source] = {"banks": bank, "state_readback": readback(bank)}
    return {
        "seed": 113031209,
        "separate_sharp_input_counts": SHOTS.tolist(),
        "law": result,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    row = run()
    path = HERE / "tomography.json"
    if args.write:
        path.write_text(json.dumps(row, indent=2) + "\n")
    elif row != json.loads(path.read_bytes()):
        raise ValueError("Changed state reconstruction")
    print(json.dumps(row["law"]["8.0"]["rotated"]["state_readback"]))
