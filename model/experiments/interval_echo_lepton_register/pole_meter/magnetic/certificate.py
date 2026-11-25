"""Enclose the principal weights of the exact magnetic band on the compact projector."""

from fractions import Fraction
import numpy as np
from model.experiments.interval_echo_angular.charged_orientation.flux_doublet import (
    bands,
)
from model.experiments.interval_echo_angular.charged_orientation.flux_doublet import (
    measure as magnetic,
)


def enclose(columns, pnum, pden):
    size = 12
    grid = bands.GRID
    # Real representation of the rounded magnetic Hamiltonian, as in band_counts.
    link = magnetic.links(1)
    cuts = bands.band_counts(link, magnetic.CELL.seams)
    assert cuts["inertia"][0] == {"cut": 41, "counts_below": 2}
    roots = bands.integer_roots()
    real = np.eye(size, dtype=np.int64) * (85 * grid)
    imag = np.zeros_like(real)
    for a, b in magnetic.CELL.seams:
        c, s = roots[int(link[a, b]) % 20]
        real[b, a] = real[a, b] = -c
        imag[b, a], imag[a, b] = -s, s
    h = np.block([[real, -imag], [imag, real]])
    x = np.rint(
        grid * np.block([[columns.real, -columns.imag], [columns.imag, columns.real]])
    ).astype(np.int64)
    energy = round(40.373373 * 2 * grid)
    # Use Python integers for the exact comparison, independently of BLAS.
    x = x.astype(object)
    identity = np.eye(4, dtype=object)
    gram = x.T @ x
    err_gram = Fraction(
        int(max(np.sum(abs(gram - grid**2 * identity), axis=1))), grid**2
    )
    residual = (h.astype(object) - energy * np.eye(24, dtype=object)) @ x
    err_res = Fraction(10 * int(np.max(abs(residual))), 2 * grid**2)
    # sqrt(24*4)<10. ||X||<2 and 1/smin(X)<2 at the certified Gram error.
    assert err_gram < Fraction(1, 10**5)
    distance = 41 - Fraction(energy, 2 * grid)
    tilt = 2 * (err_res + 2 * Fraction(250, grid)) / distance
    projector_error = tilt + err_gram
    num = np.block([[pnum, 0 * pnum], [0 * pnum, pnum]]).astype(object)
    compressed = x.T @ num @ x
    thresholds = [1318, 1320, 48940, 48943]
    inertia = []
    for cut in thresholds:
        d = 10**6 * compressed - cut * pden * grid**2 * identity
        inertia.append(bands.inertia(d))
    assert inertia == [0, 2, 2, 4]
    # PXX*P and X*PX have the same nonzero eigenvalues.
    # Replacing the orthonormal projector by XX* costs at most err_gram.
    intervals = [
        [
            float(Fraction(lo, 10**6) - projector_error),
            float(Fraction(hi, 10**6) + projector_error),
        ]
        for lo, hi in zip(thresholds[::2], thresholds[1::2], strict=True)
    ]
    assert intervals[0][0] > 0.001 and intervals[1][1] < 0.05
    return dict(
        grid=grid,
        certified_inertias=inertia,
        rounded_frame_cuts_per_million=thresholds,
        frame_gram_error=str(err_gram),
        residual_bound=str(err_res),
        projector_distance_bound=str(projector_error),
        exact_weight_intervals=intervals,
    )
