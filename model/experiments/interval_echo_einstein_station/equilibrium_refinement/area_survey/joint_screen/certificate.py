"""Exact all-start old-clock and simultaneous two-resident coincidence readout audit."""

import json
from fractions import Fraction as F

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from mpmath import iv

from ..certificate import up
from . import measure, wave

M = 2**80
B_UNIT = 2**48
C = measure.clock.MEAN_TIMEOUT


def quantize(x, unit):
    return np.vectorize(lambda y: int(round(float(y) * unit)), otypes=[object])(x)


def estimates(s, im, pnum, perr):
    n = len(s.words)
    den = wave.UNIT**2
    nominal = np.asarray(pnum, dtype=float) / den
    z = 1 - 1 / C
    old = {}
    old_err = F()
    joint_err = F()

    def refine(base, mask, boundary, f):
        free = ~mask
        mesh = np.eye(sum(free)) - z * nominal[np.ix_(free, free)]
        solver = cho_factor(mesh, check_finite=False)
        x = quantize(base, M)
        x[mask] = boundary
        for k in range(3):
            target = pnum[free] @ ((C - 1) * x + f)
            diff = target - C * den * x[free]
            resid = F(max(map(abs, diff.flat)), C * den * M)
            if k < 2:
                dx = cho_solve(
                    solver,
                    np.asarray(diff, dtype=float) / (C * den * M),
                    check_finite=False,
                )
                x[free] += quantize(dx, M)
        assert all(0 <= a <= M for a in x.flat)
        # Each inhomogeneous target (f+(C−1)x)/C lies in [0,1].
        return x, C * (resid + perr)

    for a, b in measure.EDGES:
        values = measure.clock.discounted_pair(
            s.transition, s.occupied[:, a], s.occupied[:, b], C
        )
        for ray, (start, _, trip, ret) in zip(((a, b), (b, a)), values, strict=True):
            away = s.occupied[:, ray[1]] & ~s.occupied[:, ray[0]]
            base, er = refine(ret, start, M, np.zeros(n, dtype=object))
            complete, et = refine(trip, away, base[away], np.zeros(n, dtype=object))
            old[ray] = start * complete
            old_err = max(old_err, er + et)
    y, joint_err = refine(
        im.terminal,
        np.zeros(n, dtype=bool),
        np.empty((0, 10), dtype=object),
        im.support.astype(int).astype(object) * M,
    )
    # Weighted rows include the source herald in the retained outcome.
    joint = np.vstack(
        [(s.occupied[:, i, None] * y[:, e]).T for i, e in enumerate(measure.PAIRS)]
    )
    return old, old_err, y, joint, joint_err


def exact_read(x, den, error=F()):
    interval = iv.mpf(int(x)) / int(den)
    return interval + iv.mpf([-1, 1]) * iv.mpf(error.numerator) / error.denominator


def interval_images(s, old, e_old, terminal, e_joint, tickets):
    n = len(s.words)
    den = int(sum(tickets))
    w = np.array(tickets, dtype=object)

    def lengths(weights):
        edges = []
        for a, b in measure.EDGES:
            means = []
            for ray in ((a, b), (b, a)):
                h = int(weights @ s.starts[:, measure.clock.RAYS.index(ray)])
                means.append(1024 * (1 - exact_read(weights @ old[ray], h * M, e_old)))
            edges.append(sum(means) / 2)
        return edges

    frozen = lengths(np.ones(n, dtype=object))
    live = lengths(w)

    def root3(q):
        a, b, c = q
        return iv.sqrt(3 * (a * b + a * c + b * c))

    vals = np.array([iv.mpf(0) for _ in range(8)], dtype=object)
    for i, pairs in enumerate(measure.PAIRS):
        face = []
        for edges in (frozen, live):
            a, b, c = [edges[j] ** 2 for j in pairs]
            face.append(iv.sqrt(4 * a * b - (a + b - c) ** 2) / 4)
        wt = w * s.occupied[:, i]
        hs = sum(wt)
        h = iv.mpf(int(hs)) / den
        row = [exact_read(v, den * M, e_joint * F(int(hs), den)) for v in wt @ terminal]
        joint = [row[j] for j in pairs]
        singles = [t / (2 * h) for t in row[6:]]
        sham = [
            2 * singles[a] * singles[b]
            for a, b in measure.combinations(measure.FACES[i], 2)
        ]
        sj = root3(joint)
        ss = root3(sham)
        vals += np.array(
            [
                area * v
                for area in face
                for v in (sj / sum(joint), sj, ss / sum(sham), h * h * ss)
            ],
            dtype=object,
        )
    return vals


def rank_proof(s, old, e_old, joint, e_joint, tickets):
    rays = (*measure.EDGES, (1, 0), (2, 0), (3, 0))
    n = len(s.words)
    known = [
        np.full(n, M, dtype=object),
        np.rint(2 * s.h.diagonal()).astype(int).astype(object) * (M // 512),
    ]
    known.extend(
        s.starts[:, measure.clock.RAYS.index(e)].astype(object) * M for e in rays
    )
    known.extend(old[r] for r in measure.clock.RAYS)
    known.extend(joint * 256)
    known.append((s.cut_count == 2).astype(object) * M)
    a = np.array(known, dtype=object)
    size = len(a)
    b = quantize(np.linalg.pinv(np.asarray(a, dtype=float) / M, rcond=1e-13), B_UNIT)
    eye = np.eye(size, dtype=int).astype(object) * M * B_UNIT
    defect = F(max(abs(a @ b - eye).sum(1)), M * B_UNIT)
    bsum = F(sum(abs(b).flat), B_UNIT)
    q = defect + bsum * max(e_old, 256 * e_joint)
    assert size == 36 and q < F(1, 10), float(q)
    # Stronger bound on the particular last column than on the inverse norm.
    # Its residual involves only the proposed source direction, not all columns.
    last = b[:, -1]
    rr = a @ last - eye[:, -1]
    dlast = F(max(abs(rr)), M * B_UNIT) + F(sum(abs(last)), B_UNIT) * max(
        e_old, 256 * e_joint
    )
    binf = F(max(abs(b).sum(1)), B_UNIT)
    norm = F(max(abs(last)), B_UNIT) + binf * dlast / (1 - q)
    # Move the disclosed dyadic stock to the exact joint fiber. Use each
    # registered row's launch support when propagating its kernel error.
    den = int(tickets.sum()) * n
    dw = tickets.astype(object) * n - int(tickets.sum())
    v = a @ dw
    v[-1] += den * M // 65536
    widths = [F(0)] * 11
    widths.extend(e_old * F(int(abs(dw) @ s.starts[:, i]), den) for i in range(12))
    for i in range(4):
        widths.extend([256 * e_joint * F(int(abs(dw) @ s.occupied[:, i]), den)] * 3)
    widths.append(F(0))
    assert len(widths) == size
    vnorm = F(int(max(abs(v))), den * M) + max(widths)
    correction = F(int(max(abs(b @ v))), B_UNIT * den * M)
    correction += max(abs(b) @ np.asarray(widths, dtype=object)) / B_UNIT
    correction += binf * q * vnorm / (1 - q)
    assert F(int(min(tickets)), int(tickets.sum())) - correction > F(1, 1000)

    def bulk(weights):
        probs = {}
        for word, w in zip(s.words, weights, strict=True):
            interior = tuple(v for v in word if v not in s.rim)
            probs[interior] = probs.get(interior, 0) + int(w)
        denominator = int(sum(weights))
        masses = [iv.mpf(w) / denominator for w in probs.values()]
        return -sum(t * iv.log(t) for t in masses)

    entropy_gap = bulk([1] * n) - bulk(tickets)
    # Every positive interior probability stays above .001 on the joining
    # segment; the entropy gradient's absolute coordinates are below seven.
    exact_gap = (
        entropy_gap - 7 * n * iv.mpf(correction.numerator) / correction.denominator
    )
    assert exact_gap.a > iv.mpf("0.0046"), exact_gap
    return dict(
        rank=size,
        right_inverse_defect_upper=up(q),
        finite_direction_infinity_upper=up(norm),
        reachable_halfwidth_lower=float(np.nextafter(float(1 / (n * norm)), -np.inf)),
        correction_to_exact_joint_fiber_linf_upper=up(correction),
        certified_finite_entropy_loss_lower_nats=float(
            np.nextafter(float(exact_gap.a), -np.inf)
        ),
    )


def interior_control(s):
    groups, tickets = measure.interior_reallocation(s)
    dw = tickets - 2
    assert np.count_nonzero(s.occupied.T @ abs(dw)) == 0
    energy2 = np.rint(2 * s.h.diagonal()).astype(int)
    assert dw.sum() == dw @ energy2 == 0
    # Each changed interior word has a unique interior label and no exterior resident.
    pairs = sum(count // 2 for count in groups.values())
    assert pairs == sum(dw == -1) == sum(dw == 1) == 94
    a, b = iv.mpf(3) / 2, iv.mpf(1) / 2
    drop = pairs * (a * iv.log(a) + b * iv.log(b)) / len(tickets)
    return dict(
        energy_doubled_counts=groups,
        disjoint_equal_energy_pairs=pairs,
        tickets_per_276_words=int(tickets.sum()),
        exact_bulk_entropy_loss_nats=[
            float(np.nextafter(float(drop.a), -np.inf)),
            float(np.nextafter(float(drop.b), np.inf)),
        ],
    )


def run(kernel=None):
    iv.dps = 42
    s = measure.clock.survey(1)
    im = measure.instrument(s)
    pnum, perr = wave.kernel(s.h) if kernel is None else kernel
    old, e_old, terminal, joint, e_joint = estimates(s, im, pnum, perr)
    nominal_err = F(max(abs(quantize(im.terminal, M) - terminal).flat), M) + e_joint
    assert nominal_err < F(1, 10**10)
    nominal_old = quantize(s.starts - s.scores / 1024, M).T
    clock_err = (
        F(
            max(abs(nominal_old - np.array([old[r] for r in measure.clock.RAYS])).flat),
            M,
        )
        + e_old
    )
    assert clock_err < F(1, 10**10)
    tickets = measure.dyadic_stock(s, im)
    proof = rank_proof(s, old, e_old, joint, e_joint, tickets)
    t = json.loads((measure.clock.HERE / "counts.json").read_text())["held_lottery"][
        "source_word_tickets"
    ]
    a = interval_images(s, old, e_old, terminal, e_joint, [1] * len(s.words))
    b = interval_images(s, old, e_old, terminal, e_joint, t)
    intervals = [
        [
            float(np.nextafter(float(v.a), -np.inf)),
            float(np.nextafter(float(v.b), np.inf)),
        ]
        for v in b - a
    ]
    n = len(tickets)
    total = int(sum(tickets))
    dw = tickets.astype(object) * n - total
    denom = n * total
    l1 = F(sum(abs(dw)), denom)
    old_max = (
        F(max(abs(np.array([old[r] for r in measure.clock.RAYS]) @ dw)), denom * M)
        + l1 * e_old
    )
    pair_max = F(max(abs(joint @ dw)), denom * M) + l1 * e_joint
    gain = F(int((s.cut_count == 2) @ dw), denom)
    return dict(
        **proof,
        interior_complete_record_control=interior_control(s),
        numerical_event_probability_error_upper=up(nominal_err),
        numerical_old_clock_error_upper=up(clock_err),
        new_ticket_rim_two_probability_difference=str(gain),
        new_ticket_largest_old_success_difference_upper=up(old_max),
        new_ticket_largest_joint_pixel_difference_upper=up(pair_max),
        new_ticket_max_old_herald_difference_upper=up(
            F(int(max(abs(s.starts.T @ dw))), denom)
        ),
        native_276_word_transition_row_error_upper=up(perr),
        old_return_probability_error_upper=up(e_old),
        terminal_event_probability_error_upper=up(e_joint),
        old_clock_lottery_image_changes=intervals,
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
