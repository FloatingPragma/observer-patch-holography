"""Evaluate the preregistered comparisons from recorded head, wall and clock counts."""

import numpy as np
from scipy.stats import beta

TAIL = 0.01 / (2 * 1024)


def interval(count, n):
    return [
        float(beta.ppf(TAIL, count, n - count + 1)) if count else 0.0,
        float(beta.isf(TAIL, count + 1, n - count)) if count < n else 1.0,
    ]


def subtract(a, b):
    return [a[0] - b[1], a[1] - b[0]]


def nonzero(ci):
    return ci[0] > 0 or ci[1] < 0


def evaluate(banks, shots, source_bonds, times):
    grid = {(r["k"], r["instrument"], r["t"], r["source_bond"]): r for r in banks}
    instruments = ("magnitude", "sign", "coin", "full")
    observations = {}
    differences = {}
    clock_error = 0.0
    conditional = []
    for k in (8, 0):
        observations[str(k)] = {}
        differences[str(k)] = {}
        for meter in instruments:
            deltas = []
            conditional_rows = []
            for t in times:
                pair = [grid[k, meter, t, s] for s in source_bonds]
                complete = []
                for r in pair:
                    v = np.reshape(r["counts"], (-1, 2))
                    assert v.sum() == shots
                    complete.append(interval(int(v[:, 0].sum()), shots))
                    hand = np.array(r["clock"])
                    box = 2 * np.array([interval(int(n), shots) for n in hand]) - 1
                    z = 2 * hand / shots - 1
                    angle = float(np.arctan2(z[1], z[0]))
                    # Appointment order fixes the winding because each spacing is 2 < π.
                    angle += 2 * np.pi * round((t - angle) / (2 * np.pi))
                    radius = np.linalg.norm(
                        np.maximum(abs(box[:, 0] - z), abs(box[:, 1] - z))
                    )
                    delta = np.arcsin(radius / np.linalg.norm(z))
                    clock_error = max(clock_error, float(abs(angle - t) + delta))
                deltas.append(subtract(*complete))
                for i, record in enumerate(pair[0]["record"]):
                    sign_intervals = []
                    for r in pair:
                        yes, no = r["counts"][2 * i : 2 * i + 2]
                        sign_intervals.append(interval(yes, yes + no))
                    conditional_rows.append(
                        dict(t=t, record=record, interval=subtract(*sign_intervals))
                    )
            observations[str(k)][meter] = deltas
            conditional.append(dict(k=k, instrument=meter, rows=conditional_rows))
        obs = observations[str(k)]
        for a, b, scale in [
            ("sign", "magnitude", 1),
            ("full", "magnitude", 1),
            ("coin", "magnitude", 0.5),
            ("full", "sign", 0.5),
        ]:
            key = f"{a} minus {scale:g} {b}"
            differences[str(k)][key] = [
                subtract(left, [scale * v for v in right])
                for left, right in zip(obs[a], obs[b], strict=True)
            ]
    return dict(
        clock_max_error=clock_error,
        source_response=observations,
        model_contrasts=differences,
        conditional=conditional,
    )
