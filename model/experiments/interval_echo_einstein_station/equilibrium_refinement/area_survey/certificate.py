"""Independent integer residuals for the exact 80-cell, all-start return-record kernel."""

import json
from fractions import Fraction as F
from math import isqrt
from itertools import permutations

import numpy as np
from mpmath import iv
from scipy.linalg import eigh

from . import measure

UNIT = 2**53
WAVE = 2**66
METER = 2**55
INVERSE = 2**34


def up(x):
    return float(np.nextafter(float(x), np.inf))


def integers(v, den):
    return np.vectorize(lambda x: int(round(float(x) * den)), otypes=[object])(v)


def root_upper(n):
    return isqrt(int(n)) + 1


def squared_norm(m):
    return sum(int(x) ** 2 for x in m.flat)


def transition_bound(h):
    """Exact rational residuals of a floating proposal; only a single host/clock."""
    # A scalar shift leaves every transition probability unchanged. Centering
    # the integer energies reduces cancellation in the proposed eigenpairs.
    shift = int(np.rint(np.trace(h) / len(h)))
    h = h - shift * np.eye(len(h))
    vals, vec = eigh(h, driver="evd")
    v = integers(vec, UNIT)
    values = integers(vals, UNIT)
    identity = np.eye(len(v), dtype=int).astype(object)
    defect = v.T @ v - UNIT**2 * identity
    orth = F(root_upper(squared_norm(defect)), UNIT**2)
    h2 = (2 * h).astype(int)
    assert np.array_equal(h2, 2 * h)
    hv = np.array(
        [sum(h2[i, j] * v[j] for j in np.flatnonzero(h2[i])) for i in range(len(v))],
        dtype=object,
    )
    res = hv * UNIT - 2 * v * values[None, :]
    eig = F(root_upper(max(abs(res.T @ res).sum(axis=1))), 2 * UNIT**2)
    iv.dps = 45
    parts = []
    for fn in (iv.cos, iv.sin):
        trig = []
        for eig_int in values:
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
            gram_error_upper=up(orth),
            scalar_energy_shift=shift,
            eigen_residual_operator_upper=up(eig),
            transition_row_l1_upper=up(row_error),
        ),
    )


def response_rows(s, pnum, row_error):
    n = len(s.words)
    c = measure.MEAN_TIMEOUT
    den = WAVE**2
    rows = {}
    max_error = F()
    max_res = F()
    interiors = sorted({v for w in s.words for v in w} - set(s.rim))
    stations = [*s.rim, *interiors]
    occupied = np.array([[v in w for v in stations] for w in s.words])
    spans = {}
    edges = [
        *measure.EDGES,
        *((i, j) for i in range(4, len(stations)) for j in range(4)),
    ]
    for a, b in edges:
        outcomes = measure.discounted_pair(
            s.transition, occupied[:, a], occupied[:, b], c
        )
        leg = []
        for ray, (enrolled, _, score, home_score) in zip(
            ((a, b), (b, a)), outcomes, strict=True
        ):
            away = occupied[:, ray[1]] & ~occupied[:, ray[0]]
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
            if max(ray) < 4:
                rows[ray] = trip * enrolled
            leg.append(F(int(trip[enrolled].sum()), int(enrolled.sum()) * METER))
        d = (1 - sum(leg) / 2) * F(c, 4)
        spans[stations[a], stations[b]] = spans[stations[b], stations[a]] = d
    included = (*measure.EDGES, (1, 0), (2, 0), (3, 0))
    # H_ab-H_ba=n_a-n_b, so these nine span every ordered launch rate.
    e2 = (2 * s.h.diagonal()).astype(int)
    known = [np.ones(n, dtype=int) * METER, e2.astype(object) * (METER // 512)]
    known.extend(
        s.starts[:, measure.RAYS.index(ray)].astype(object) * METER for ray in included
    )
    known.extend(rows[ray] for ray in measure.RAYS)
    known.append((s.cut_count == 2).astype(object) * METER)
    return (
        np.asarray(known, dtype=object),
        max_error,
        max_res,
        rim_test(s.rim, interiors, spans, max_error * F(c, 4)),
    )


def rim_test(rim, interiors, spans, error):
    def ivrat(x):
        return iv.mpf(x.numerator) / x.denominator

    def d2(a, b):
        if a == b:
            return iv.mpf(0)
        return (ivrat(spans[a, b]) + iv.mpf([-1, 1]) * ivrat(error)) ** 2

    def determinant(g):
        n = len(g)
        total = iv.mpf(0)
        for q in permutations(range(n)):
            term = (-1) ** sum(q[i] > q[j] for i in range(n) for j in range(i + 1, n))
            for i in range(n):
                term *= g[i][q[i]]
            total += term
        return total

    test = [*rim, interiors[0]]
    gram = [
        [(d2(test[0], a) + d2(test[0], b) - d2(a, b)) / 2 for b in test[1:]]
        for a in test[1:]
    ]
    height = determinant(gram) / determinant([r[:3] for r in gram[:3]])
    spread = (
        min(
            max(spans[p, b] for b in rim) - min(spans[p, b] for b in rim)
            for p in interiors
        )
        - 2 * error
    )
    assert height.a > 36000 and spread > F(7139, 100), (height, float(spread))
    return dict(
        held_station=interiors[0],
        flat_3d_extra_squared_length_lower=36000,
        native_interior_centers_checked=len(interiors),
        minimum_rim_radius_spread_lower=71.39,
    )


def run():
    s = measure.survey(1)
    pnum, row_error, diag = transition_bound(s.h)
    a, prob_error, solve_error, rim = response_rows(s, pnum, row_error)
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
    b_norm = F(max(abs(b).sum(axis=1)), INVERSE)
    # Only the last right inverse column is used in this trajectory.
    norm = F(max(abs(b[:, -1])), INVERSE) + b_norm * q / (1 - q)
    step = 1 / (len(s.words) * norm)
    return dict(
        host=s.board_cells,
        stations=s.rim,
        **diag,
        solve_defect_upper=up(solve_error),
        trip_probability_error_upper=up(prob_error),
        rows=n,
        held_geometry=rim,
        right_inverse_defect_upper=up(q),
        direction_infinity_norm_upper=up(norm),
        positive_interval_halfwidth_lower=float(np.nextafter(float(step), -np.inf)),
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
