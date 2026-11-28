"""Reported-word contrasts and regional thermal instruments for the same register."""

from fractions import Fraction as Q

import numpy as np
from scipy.sparse import eye, diags
from scipy.sparse.linalg import splu

from . import closed_balls as balls
from . import sector_regional
from . import sector_retention as experiment
from . import sector_thermal as thermal


def sources(divisors=(1, 8, 64)):
    words = thermal.apparatus()[0]
    original = experiment.peer_population.fibers.source.records()[1]["original"]
    _, arms = experiment.peer_population.histogram(16)
    contrast = arms["agreed"] - original
    full, initial, _, service, _ = experiment.acquisition.returns.retention.setup()
    ix = {tuple(row): i for i, row in enumerate(full)}
    subset = np.array([ix[tuple(row)] for row in words])
    assert set(subset) == set(np.flatnonzero(initial))
    mats = [M for M, _ in service]
    S = (mats[2] @ mats[1] @ mats[0]).T.tocsr()
    den = experiment.ReturnClock.den
    assert np.all(S.sum(0) == den) and np.all(S.sum(1) == den)
    qfull = np.zeros(len(full))
    qfull[subset] = contrast.astype(float)
    out = {0: (qfull[subset], 1e-14)}
    for div in divisors:
        report = 1 / Q(div)
        numer, denom = report.numerator, report.denominator
        survive = np.full(len(full), denom)
        survive[subset] -= numer
        killed = eye(len(full)) - S @ diags(survive / (denom * den))
        x = splu(killed.tocsc()).solve(qfull, trans="T")
        # Certify the exit measure from an integer residual. Each remaining
        # unit of unsigned injection generates at most one accepted word.
        # The existing return-time certificate and geometric reporting imply
        # absorption with probability one for each positive reporting rate.
        grid = 2**39
        while (float(np.max(abs(x))) * grid + 1) * den * denom >= 2**61:
            grid //= 2
        v = np.rint(grid * x).astype(np.int64)
        assert int(abs(v).max()) * den * denom < 2**61
        arrival = S.T @ v
        resid = v * den * denom - arrival * survive
        scale = grid * den * denom
        start = {int(i): s for i, s in zip(subset, contrast, strict=True)}
        bound = sum(
            abs(Q(int(r), scale) - start.get(i, 0)) for i, r in enumerate(resid)
        )
        d = arrival[subset] * numer / scale
        eps = float(bound) + 1e-14
        assert abs(d.sum()) <= eps and eps < 1e-7 and np.min(1 / len(words) + d) > 0
        out[div] = (d, eps)
    return out, contrast, subset, full, mats


def meters(words, banks, requested):
    rho = thermal.bath()
    measured = sector_regional.calibrate(banks)
    sharp = {}
    coherent = {}
    intervals = {}
    for mask in sorted(requested):
        ids = balls.mask_ids(words, mask)
        ext = balls.mask_ids(words, 127 ^ mask)
        i, j = np.nonzero(ext[:, None] == ext[None, :])
        size = int(ids.max()) + 1
        region = (ids, i, j, ids[i] * size + ids[j], size)
        local = {}
        for tag, base in rho.items():
            K, price = thermal.operators(base, region)
            prob = np.bincount(ids, weights=np.diag(base))
            np.testing.assert_allclose(price, -np.log(prob), atol=1e-11)
            lohi = np.vstack([np.bincount(ids, weights=q) for q in measured[tag].T]).T
            assert np.all((lohi[:, 0] <= prob) & (prob <= lohi[:, 1]))
            local[tag] = (
                -np.log(prob[ids]),
                K.diagonal()[ids],
                -np.log(lohi[ids, ::-1]),
            )
        r, c = local["retained"], local["exchanged"]
        sharp[mask] = r[0] - c[0]
        coherent[mask] = r[1] - c[1]
        intervals[mask] = r[2] - c[2][:, ::-1]
        # The sharp regional likelihood ratio is the retained-reference
        # conditional mean of the full likelihood given the visible word.
        ratio = np.diag(rho["exchanged"]) / np.diag(rho["retained"])
        pooled = np.bincount(
            ids, weights=np.diag(rho["retained"]) * ratio
        ) / np.bincount(ids, weights=np.diag(rho["retained"]))
        np.testing.assert_allclose(np.exp(sharp[mask]), pooled[ids], rtol=1e-12)
    np.testing.assert_allclose(sharp[127], coherent[127], atol=1e-10)
    return sharp, coherent, intervals


def group_mean(ids, n, score):
    return (np.bincount(ids, weights=score, minlength=len(n)) / n)[ids]
