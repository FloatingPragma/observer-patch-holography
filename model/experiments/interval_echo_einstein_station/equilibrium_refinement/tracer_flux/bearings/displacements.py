"""Locate the calibrated marker stencil among sharp, readable two-resident positions."""

import json
from itertools import combinations
from math import log, sqrt

import numpy as np

from .. import measure as tracer
from .measure import bearing_frame, exchange_certificate


def shell_records(edges, outside):
    axes, gram = bearing_frame(edges)
    eye = np.eye(12, dtype=int)
    opposite = np.zeros_like(eye)
    for a, b in axes:
        opposite[a, b] = opposite[b, a] = 1
    f = eye - opposite
    s = edges.astype(int) - (1 - eye - opposite - edges)
    words = list(combinations((i for i in range(12) if i != outside), 2))
    # Word vectors give cumulative signed load relative to a balanced opposite pair.
    z = np.array([2 * (eye[a] + eye[b]) for a, b in words])
    qf, qs = z @ f @ z.T, z @ s @ z.T
    ff, ss = qf.diagonal(), qs.diagonal()
    centered = np.flatnonzero((ff == 0) & (ss == 0))
    assert len(centered) == 5
    d0 = ff[:, None] + ff[None, :] - 2 * qf
    d1 = ss[:, None] + ss[None, :] - 2 * qs
    dsq = d0 + d1 / np.sqrt(5)
    nonzero = (d0 != 0) | (d1 != 0)
    min_d2 = dsq[nonzero].min()
    tables = {}
    for side in (-1, 1):
        shell = np.flatnonzero((ff == 8) & (ss == side * 8))
        pairs = [
            (a, b)
            for a, b in combinations(shell, 2)
            if qf[a, b] == -8 and qs[a, b] == -side * 8
        ]
        assert len(pairs) == 10 and len(shell) == 25
        reps = [a for a, _ in pairs]
        sq = ((z[reps] @ gram[:, outside]) ** 2) / (8 * (1 + side / np.sqrt(5)))
        tilted = sq > 0.01
        x = (1 - side / np.sqrt(5)) / 2
        np.testing.assert_allclose(sq, np.where(tilted, x, 0), rtol=0, atol=3e-13)
        assert tilted.sum() == 5
        # The board has a fivefold rotation preserving the outside registration.
        # Verify the isotropic identity independently in the twelve receiver labels.
        mix = 1 / (3 * x)
        weights = np.where(tilted, mix / 5, (1 - mix) / 5)
        rl2 = 8 * (1 + side / np.sqrt(5))
        projection = z[reps] @ gram
        np.testing.assert_allclose(
            (projection.T * weights) @ projection / rl2, gram / 3, atol=3e-14, rtol=0
        )
        # In Q(sqrt(5)) the per-axis weights times 30 are 5+side*sqrt(5)
        # on the tilted orbit and 1-side*sqrt(5) on the equatorial orbit.
        y, zz = z[reps] @ f, z[reps] @ s
        a, b = np.where(tilted, 5, 1), np.where(tilted, side, -side)
        const = 5 * y[:, :, None] * y[:, None, :] + zz[:, :, None] * zz[:, None, :]
        root = y[:, :, None] * zz[:, None, :] + zz[:, :, None] * y[:, None, :]
        np.testing.assert_array_equal(
            np.einsum("i,ijk->jk", a, const) + 5 * np.einsum("i,ijk->jk", b, root),
            400 * f + 80 * side * s,
        )
        np.testing.assert_array_equal(
            np.einsum("i,ijk->jk", a, root) + np.einsum("i,ijk->jk", b, const),
            80 * (s + side * f),
        )
        # A squared axial cosine exceeds 1/3 on the tilted orbit only for the inner shell.
        for idx, i in enumerate(reps):
            y0, z0 = int(y[idx, outside]), int(zz[idx, outside])
            if not tilted[idx]:
                assert y0 == 0 and z0 == 0
            else:
                # r^2 * (1-side/sqrt(5))/2 is 16/5 on both shells.
                assert 5 * y0 * y0 + z0 * z0 == 16 and y0 * z0 == 0
        tables[str(side)] = dict(
            radius_squared=float(rl2),
            endpoints=len(shell),
            antipodal_word_pairs=[[words[a], words[b]] for a, b in pairs],
            tilted_total_weight=float(mix),
            max_paired_squared_axial_cosine=float(x),
            positive_isotropic_weights=bool(np.all(weights > 0)),
        )
    # For the inner radius every half-aperture surrounding a marker is a single
    # Gram position. The centered larger aperture adds the entire smaller shell.
    for a, b in zip(d0[nonzero] - 2, d1[nonzero] + 2, strict=True):
        # a*sqrt(5)+b>0 exactly, for the squared displacement minus r_inner^2/4.
        assert (
            (a >= 0 and b >= 0 and (a or b))
            or (a > 0 and b < 0 and 5 * a * a > b * b)
            or (a < 0 and b > 0 and b * b > 5 * a * a)
        )
    count_c = len(centered) + tables["-1"]["endpoints"]
    count_h = len(centered)
    # Common factor 1/(55 times the large flat reference volume) cancels.
    h, b = 8 * count_h, 8
    curvature_numerator = 3 * (18 * b + 22 * h - 40 * count_c)
    curvature_denominator = 4 * h - count_c
    # Independent copies of the uniform monitored source for every point survey.
    # There are twenty equal-probability small apertures; their calibration weights
    # are nonnegative with unit sum. One simultaneous error bound covers all 22 banks.
    shots = 2**32
    rng = np.random.default_rng(280571)
    counts = rng.binomial(shots, np.array([count_c, count_h] + [1] * 20) / len(words))
    orbit = tables["-1"]["tilted_total_weight"]
    endpoints = tables["-1"]["antipodal_word_pairs"]
    is_tilted = np.array(
        [
            abs((z[words.index(tuple(w))] @ gram)[outside]) > 0.1
            for pair in endpoints
            for w in pair
        ]
    )
    weights = np.where(is_tilted, orbit / 10, (1 - orbit) / 10)
    c_meas, h_meas = counts[0] / shots, 8 * counts[1] / shots
    b_meas = 8 * float(weights @ counts[2:]) / shots
    upperr = sqrt(log(2 * len(counts) / 0.01) / (2 * shots))
    num = 3 * (18 * b_meas + 22 * h_meas - 40 * c_meas)
    den = 4 * h_meas - c_meas
    limits = [
        x / y
        for x in (num - 1080 * upperr, num + 1080 * upperr)
        for y in (den - 33 * upperr, den + 33 * upperr)
    ]
    assert den > 33 * upperr
    return dict(
        shots_per_marker_bank=shots,
        marker_counts=counts.tolist(),
        observed_scalar_times_radius_squared=num / den,
        simultaneous_99_count_interval=[min(limits), max(limits)],
        words=len(words),
        shared_origin_records=[words[i] for i in centered],
        shells=tables,
        minimum_nonzero_distance_squared=float(min_d2),
        small_center_count=count_h,
        large_center_count=count_c,
        displaced_small_count=1,
        apparent_scalar_times_radius_squared=[
            curvature_numerator,
            curvature_denominator,
        ],
    )


def run():
    twice, basis, _, outside = tracer.source(8)
    tags = np.array([p for p, _ in basis])
    _, edges = exchange_certificate(twice, tags)
    return shell_records(edges, outside)


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
