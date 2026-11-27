"""Measure an operational rank bound and the held three-port cancellation."""

from math import ceil, log, pi, sin, sqrt

import numpy as np

from . import measure as m
from ..apparatus import geometry


def ideal_geometry():
    rings, anti, _ = geometry()
    bare = []
    for cycle in rings:
        q = np.zeros((12, 12))
        for a, i in enumerate(cycle):
            for b, j in enumerate(cycle):
                q[i, j] = sin(2 * pi * (a - b) / 5) / 5
        bare.append(q)
    bare = np.array(bare)
    np.testing.assert_allclose(bare.sum(axis=0), 0, atol=1e-15)
    five = (bare + bare[list(anti)]) / 2
    i, j = m.COORDS.T
    rows = [0.25 * q[i, j] for q in bare]
    for b in m.CENTERS:
        comm = five[0] @ five[b] - five[b] @ five[0]
        rows.append(comm[i, j] / 16)
    return np.column_stack(rows)


def read_counts(data, key, shape, shots):
    if key not in data:
        raise ValueError(f"Missing {key} bank")
    arr = np.asarray(data[key])
    if (
        arr.shape != shape
        or arr.dtype.kind not in "iu"
        or (arr < 0).any()
        or (arr > shots).any()
    ):
        raise ValueError(f"Malformed {key} records")
    if not (arr.sum(axis=-1) == shots).all():
        raise ValueError(f"Incomplete {key} bank")
    return arr.astype(np.int64)


def order_contrast(arr):
    signed = (arr[..., 6] - arr[..., 7]) / m.PAIR_SHOTS
    alternating = signed[0] - signed[1]
    cross = (
        alternating[0, 0] - alternating[0, 1] - alternating[1, 0] + alternating[1, 1]
    )
    return cross.mean(axis=0) / 4


def halfwidth(n, family, budget):
    # One bounded [-1,1] read per start, n is reciprocal sum(w_i²/N_i).
    return sqrt(2 * log(2 * family / budget) / n)


def decode(data):
    for key, expected in [
        ("coords", m.COORDS),
        ("cycle", m.CYCLE),
        ("shots", [m.SINGLE_SHOTS, m.PAIR_SHOTS]),
        ("seed", m.SEED),
        ("path_lengths", [len(p) - 1 for p in m.lab.layout()["contacts"]]),
    ]:
        if key not in data or not np.array_equal(data[key], expected):
            raise ValueError(f"Wrong {key} enrollment")
    singles = read_counts(data, "single", (2, 12, 2, 66, 9), m.SINGLE_SHOTS)
    orders = read_counts(data, "double", (2, 2, 2, 2, 4, 66, 9), m.PAIR_SHOTS)
    held = read_counts(data, "held", (3, 2, 2, 2, 2, 4, 3, 9), m.PAIR_SHOTS)
    standards = read_counts(data, "standards", (2, 4, 2), m.SINGLE_SHOTS)
    if not np.all(standards[:, :, 0] == m.SINGLE_SHOTS * np.array([0, 1, 0, 1])):
        raise ValueError("Polarity or gain standard failed")
    if singles[..., 3:6].any():
        raise ValueError("Unscheduled second encounter recorded")
    pred = ideal_geometry()
    scalar_tests = (
        4096  # 1584 marked columns plus the raw order coordinates and held contrasts.
    )
    score = (singles[..., 6] - singles[..., 7]) / m.SINGLE_SHOTS
    one = (score[:, :, 1] - score[:, :, 0]) / 2
    # Per conditional binomial, Hoeffding confines the ideal outcome to at
    # most 2r+1 central integer bins except for mass eps. Each central bin
    # and each tail interval changes by at most one midpoint CDF grid tick.
    # Couple successive conditional draws and sum their total-variation costs.
    exposure_banks = sum(
        int(np.prod(v.shape[:-1])) for v in (singles, orders, held, standards)
    )
    eps = 2.0**-80
    radius = ceil(sqrt(m.SINGLE_SHOTS * log(2 / eps) / 2))
    quantile_budget = 8 * exposure_banks * ((2 * radius + 5) / 2**52 + 2 * eps)
    quant = 9 / 2**26  # signed expected event read, including the largest-bin residual.
    single_width = (
        halfwidth(2 * m.SINGLE_SHOTS, scalar_tests, 0.005 - quantile_budget) + quant
    )
    order_width = (
        halfwidth(8 * m.PAIR_SHOTS, scalar_tests, 0.005 - quantile_budget) + 2 * quant
    )
    max_single = float(abs(one - pred[:, :12].T).max())
    if max_single > single_width:
        raise ValueError("Single-source calibration failed")
    cols = [order_contrast(a) for a in orders]
    obs = np.column_stack([one.mean(axis=0).T, *cols])
    max_order = float(abs(obs[:, 12:] - pred[:, 12:]).max())
    if max_order > order_width:
        raise ValueError("Measured order response misses forecast")
    # Singular vectors come exclusively from the pre-enrolled geometrical forecast.
    left, _, right = np.linalg.svd(pred, full_matrices=False)
    left = left[:, :13]
    right = right[:13].T
    small = left.T @ obs @ right
    least = float(np.linalg.svd(small, compute_uv=False)[-1])
    weighted_n = min(4 * m.SINGLE_SHOTS, 8 * m.PAIR_SHOTS)
    spectral_radius = 13 * halfwidth(weighted_n, 13**2, 0.005) + 2 * quant * sqrt(
        pred.size
    )
    if least <= spectral_radius:
        raise ValueError("Thirteenth response direction unresolved")
    idx = [m.COORDS.tolist().index(x) for x in m.CYCLE.tolist()]
    weights = -1 / np.sin(np.array([1, 2, 1]) * pi / 5)
    weights /= abs(weights).sum()
    np.testing.assert_allclose(weights @ pred[idx, :12], 0, atol=1e-15)
    predicted = float(weights @ pred[idx, 12])
    observed = float(weights @ obs[idx, 12])
    bound = (
        halfwidth(
            8 * m.PAIR_SHOTS / np.dot(weights, weights),
            scalar_tests,
            0.005 - quantile_budget,
        )
        + 2 * quant
    )
    held_results = []
    for experiment, (rate, time) in enumerate(m.CONTROLS):
        scale = sin(rate * time) ** 4
        values = [float(weights @ order_contrast(h)) for h in held[experiment]]
        if max(abs(values[0] - scale * predicted), abs(values[1])) > bound:
            raise ValueError("Held cycle quadrature or response scale failed")
        held_results.append(
            {
                "rate": rate,
                "dwell": time,
                "forecast": predicted * scale,
                "real": values[0],
                "quarter_phase": values[1],
            }
        )
    if abs(observed - predicted) > bound:
        raise ValueError("Primary cycle missed")
    total = sum(int(v.sum()) for v in (singles, orders, held, standards))
    failures = sum(int(v[..., :6].sum()) for v in (singles, orders, held))
    return dict(
        total_simulated_starts=total,
        count_banks=exposure_banks,
        failed_encounters=failures,
        smallest_projected_singular=least,
        spectral_error=spectral_radius,
        lower_99=least - spectral_radius,
        maximum_single_error=max_single,
        joint_single_radius=single_width,
        maximum_order_error=max_order,
        joint_order_radius=order_width,
        cycle=observed,
        cycle_forecast=predicted,
        cycle_joint_radius=bound,
        held=held_results,
    )


def main():
    """Score archived terminal counts without regenerating the encounters."""
    import json

    with np.load(m.HERE / "counts.npz", allow_pickle=False) as archive:
        print(json.dumps(decode(archive), indent=2))


if __name__ == "__main__":
    main()
