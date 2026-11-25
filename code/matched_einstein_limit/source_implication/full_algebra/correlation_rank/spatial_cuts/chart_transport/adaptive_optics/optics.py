"""Bounded native optical ports and directional counts for two independently calibrated receivers."""

import numpy as np
from mpmath import iv
from scipy.stats import beta

from matched_einstein_limit.source_implication.full_geometry.station_balls.reciprocal_screen import (
    measure as camera,
)
from matched_einstein_limit.source_implication.directional_stress.traction.certificate import (
    four_tick_columns,
)
from matched_einstein_limit.source_implication.full_algebra import (
    native_precision as native,
)


def matrices(blocks):
    return [
        iv.matrix(
            [[iv.mpc(float(z.real), float(z.imag)) for z in row] for row in block]
        )
        for block in blocks
    ]


def adjoint(a):
    return iv.matrix(
        [
            [iv.mpc(a[j, i].real, -a[j, i].imag) for j in range(a.rows)]
            for i in range(a.cols)
        ]
    )


def kernels():
    board = camera.camera.front.optics.pairing.formation.federation(1)
    route = camera.camera.front.optics.pairing.formation.region(board)
    rings = camera.camera.screens(board, route)
    pos = camera.camera.indices(np.r_[rings[0], rings[1]])
    weights = native.prior.weights()[7]
    for charge in (-1, 1):
        exact = []
        largest = iv.mpf(0)
        for odd in (False, True):
            h = native.generator(board, route, charge, odd)
            raw = (
                four_tick_columns(h, pos, pos, weights)
                .reshape(2, 10, 4, 10, 4)
                .transpose(0, 1, 3, 2, 4)
            )
            blocks = {}
            for i in range(10):
                for j in range(10):
                    if (i < 5) == (j < 5):
                        continue
                    re, im = raw[:, i, j]
                    m = iv.matrix(
                        [
                            [
                                iv.mpc(int(a), int(b)) / native.prior.FRAC
                                for a, b in zip(row, col)
                            ]
                            for row, col in zip(re, im)
                        ]
                    )
                    blocks[i, j] = m, adjoint(m)
                    norm = (
                        iv.sqrt(
                            iv.mpf(
                                sum(
                                    int(a) ** 2 + int(b) ** 2
                                    for a, b in zip(re.flat, im.flat)
                                )
                            )
                        )
                        / native.prior.FRAC
                    )
                    largest = max(largest, (norm + iv.mpf("2e-19")).b)
            exact.append(blocks)
        yield charge, exact, largest**2, camera.propagate(board, route, rings, charge)


def directional(kernels, source, effects):
    x = matrices(source)
    terminal = [matrices(e) for e in effects]
    retained = np.zeros((10, 10))
    for i, j in kernels[0]:
        s = iv.mpf(0)
        for pi in range(2):
            k, t = kernels[pi][i, j]
            op = terminal[int(i >= 5)][pi] * k * x[pi] * t
            s += sum(op[u, u] for u in range(4)).real
        assert s.b - s.a < iv.mpf("1e-28")
        retained[i, j] = float(s.mid)
    # For exact arithmetic on these numerical model stocks the true kernel
    # induces <3e-19 error per probability, before final float conversion.
    return retained[5:, :5], retained[:5, 5:]


def counts(f, r, rng):
    # Bonferroni coverage for 25 pairs, two directions, both probes and both
    # frozen/adapted banks at one source (200 scalar binomial marginals).
    obs = []
    limits = []
    for arm in (f, r):
        totals = [
            rng.multinomial(camera.SHOTS, np.r_[col, 1 - col.sum()]) for col in arm.T
        ]
        n = np.array(totals)[:, :5].T
        alpha = 0.01 / (2 * 200)
        lo = np.nan_to_num(beta.ppf(alpha, n, camera.SHOTS - n + 1))
        hi = np.nan_to_num(beta.ppf(1 - alpha, n + 1, camera.SHOTS - n), nan=1.0)
        obs.append(np.asarray(totals).tolist())
        limits.append((lo, hi))
    (fl, fh), (rl, rh) = limits
    return dict(
        forward=obs[0],
        reverse=obs[1],
        simultaneous99_max_pair_mismatch_lower=float(
            max(0, np.max(fl - rh.T), np.max(rl.T - fh))
        ),
    )
