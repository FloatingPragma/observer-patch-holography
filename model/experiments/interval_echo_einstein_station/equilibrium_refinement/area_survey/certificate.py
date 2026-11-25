"""Independent integer residuals for the exact 80-cell, all-start return-record kernel."""

import json
from fractions import Fraction as F
from math import isqrt

import numpy as np
from mpmath import iv
from scipy.linalg import eigh

from . import measure

UNIT = 2**53
WAVE = 2**66
METER = 2**55
INVERSE = 2**34


def integers(v, den):
    return np.vectorize(lambda x: int(round(float(x) * den)), otypes=[object])(v)


def root_upper(n):
    return isqrt(int(n)) + 1


def squared_norm(m):
    return sum(int(x) ** 2 for x in m.flat)


def transition_bound(h):
    """Exact rational residuals of a floating proposal; only a single host/clock."""
    vals, vec = eigh(h)
    v = integers(vec, UNIT)
    l = integers(vals, UNIT)
    identity = np.eye(len(v), dtype=int).astype(object)
    defect = v.T @ v - UNIT**2 * identity
    orth = F(root_upper(squared_norm(defect)), UNIT**2)
    h2 = (2 * h).astype(int)
    assert np.array_equal(h2, 2 * h)
    hv = np.array(
        [sum(h2[i, j] * v[j] for j in np.flatnonzero(h2[i])) for i in range(len(v))],
        dtype=object,
    )
    res = hv * UNIT - 2 * v * l[None, :]
    eig = F(root_upper(squared_norm(res)), 2 * UNIT**2)
    iv.dps = 45
    parts = []
    for fn in (iv.cos, iv.sin):
        trig = []
        for eig_int in l:
            angle = -iv.mpf(int(eig_int)) / (2 * UNIT)
            interval = fn(angle)
            t = int(float(interval.mid) * WAVE)
            assert abs(interval - iv.mpf(t) / WAVE).b < iv.mpf("2e-16")
            trig.append(t)
        q = (v * np.asarray(trig, dtype=object)) @ v.T
        parts.append((q + UNIT**2 // 2) // UNIT**2)
    n = len(v)
    # U-V exp(-itL)V^T = U(I-VV^T)+(UV-V exp(-itL))V^T.
    # Duhamel bounds the second term by t*||HV-VL|| times ||V||.
    wave_error = (
        orth + eig * (1 + orth) / 2 + (1 + orth) * F(2, 10**16) * 2 + F(n, WAVE)
    )
    row_error = 2 * wave_error + wave_error**2
    assert row_error < F(2, 10**11)
    return (
        parts[0] ** 2 + parts[1] ** 2,
        row_error,
        dict(
            gram_error_upper=float(orth),
            eigen_residual_frobenius_upper=float(eig),
            transition_row_l1_upper=float(row_error),
        ),
    )


def response_rows(s, pnum, row_error):
    n = len(s.words)
    c = measure.MEAN_TIMEOUT
    den = WAVE**2
    rows = {}
    max_error = F()
    max_res = F()
    for a, b in measure.EDGES:
        outcomes = measure.discounted_pair(
            s.transition, s.occupied[:, a], s.occupied[:, b], c
        )
        for ray, (enrolled, _, score, home_score) in zip(
            ((a, b), (b, a)), outcomes, strict=True
        ):
            away = s.occupied[:, ray[1]] & ~s.occupied[:, ray[0]]
            trip = integers(score, METER)
            home = integers(home_score, METER)
            assert all(0 <= x <= METER for x in np.r_[trip, home])
            assert np.array_equal(trip[away], home[away]) and all(
                home[enrolled] == METER
            )
            total = F()
            for x, mask in ((home, enrolled), (trip, away)):
                # x solves x=z P x off the absorbing target; a boundary error
                # in the second solve is transmitted with gain at most one.
                d = c * den * x[~mask] - (c - 1) * (pnum[~mask] @ x)
                r = F(max(map(abs, d)), c * den * METER)
                max_res = max(r, max_res)
                total += c * (r + row_error)
            max_error = max(max_error, total)
            rows[ray] = trip * enrolled
    included = (*measure.EDGES, (1, 0), (2, 0), (3, 0))
    # H_ab-H_ba=n_a-n_b, so these nine span every ordered launch rate.
    e2 = (2 * s.h.diagonal()).astype(int)
    known = [np.ones(n, dtype=int) * METER, e2.astype(object) * (METER // 512)]
    known.extend(
        s.starts[:, measure.RAYS.index(ray)].astype(object) * METER for ray in included
    )
    known.extend(rows[ray] for ray in measure.RAYS)
    known.append((s.cut_count == 2).astype(object) * METER)
    return np.asarray(known, dtype=object), max_error, max_res


def run():
    s = measure.survey(1)
    pnum, row_error, diag = transition_bound(s.h)
    a, prob_error, solve_error = response_rows(s, pnum, row_error)
    n = a.shape[0]
    middle = np.asarray(a, dtype=float) / METER
    b = integers(np.linalg.pinv(middle), INVERSE)
    eye = np.eye(n, dtype=int).astype(object) * METER * INVERSE
    raw = a @ b - eye
    residual = F(max(abs(raw).sum(axis=1)), METER * INVERSE)
    uncertainty = F(sum(abs(b).flat), INVERSE) * prob_error
    q = residual + uncertainty
    assert n == 24 and q < F(1, 10)
    # Exact right inverse B(AB)^−1. p* ± tD is positive at all words.
    norm = F(max(abs(b).sum(axis=1)), INVERSE) / (1 - q)
    step = 1 / (len(s.words) * norm)
    return dict(
        host=s.board_cells,
        stations=s.rim,
        **diag,
        solve_defect_upper=float(solve_error),
        trip_probability_error_upper=float(prob_error),
        rows=n,
        right_inverse_defect_upper=float(q),
        direction_infinity_norm_upper=float(norm),
        positive_interval_halfwidth_lower=float(step),
        target="unit increase in probability of two occupied rim stations",
        fixed_linear_readouts=[
            "normalization",
            "native mean energy",
            "all 12 launch probabilities",
            "all 12 successful return probabilities including the independent timeout",
        ],
        entropy_slope_nats="(4 log(4) - 3 log(6))/3",
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
