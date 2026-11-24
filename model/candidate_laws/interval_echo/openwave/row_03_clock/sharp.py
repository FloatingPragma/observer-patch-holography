"""Return records of an unmonitored vacancy and its certified persistent beat."""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo.openwave.row_04_compact import (
    anchored_vacancy as anchor,
)
from model.candidate_laws.interval_echo.openwave.row_04_compact import (
    experiment as compact,
)


def unmonitored(level, shots, rng):
    board = anchor.federation(level)
    local, axes = compact.geometry()
    home = next(
        i
        for i, tips in enumerate(compact.gateway_table(board))
        if any(set(axis) <= set(tips) for axis in axes)
    )
    h = anchor.defect_generator(board, -1, anchor.K, 0)
    site = 12 * home + 4
    sharp = np.zeros(board.n)
    sharp[site] = 1
    # In the interior H = (69 I - A)/2. Both integer identities below are checked
    # on the entire federation, so they include any amplitude leaving home.
    u = np.zeros(board.n, dtype=int)
    u[12 * home + np.array([4, 5, 6, 7])] = [1, 1, -1, -1]
    b = (69 * sparse.eye(board.n) - 2 * h).tocsr()
    np.testing.assert_array_equal(b.data, np.rint(b.data))
    b = b.astype(np.int64)
    v = b @ u
    np.testing.assert_array_equal(b @ v, 5 * u)
    overlaps = []
    for sign in (-1, 1):
        eigenvector = np.sqrt(5) * u + sign * v
        overlaps.append(
            float(eigenvector[site] ** 2 / np.dot(eigenvector, eigenvector))
        )
    centered = h - 32 * sparse.eye(board.n, format="csr")
    waves = expm_multiply(-1j * centered, sharp, start=0, stop=16, num=65)
    returns = np.abs(waves[:, site]) ** 2
    click_counts = [int(rng.binomial(shots, round(p, 12))) for p in returns]
    return {
        "cells": board.cells,
        "source_cell": home,
        "source_port": site % local.n,
        "squared_overlap_each_compact_mode": overlaps,
        "return_cosine_amplitude_lower_bound": 2 * overlaps[0] * overlaps[1],
        "inner_products_uu_uv_vv": [int(x) for x in (u @ u, u @ v, v @ v)],
        "sharp_components_u_v": [int(u[site]), int(v[site])],
        "clock_appointments": "0 to 16 exchange exposures in equal increments of 1/4",
        "return_counts": click_counts,
        "other_counts": [shots - n for n in click_counts],
    }
