"""Exact head lotteries with equal strain jet and equal outward number currents.

The three input rows use Z=diag(+I,-I) and native exchange time at pump age 4.
Heads, both marginal axes and the classical 2:1 rail preparation are those of
the formed vacancy source. All receipts refer to differences of two lotteries.

The energy meters use the symmetric half-bond cut. Under H → H + cI their
outward currents shift by c times the matched outward number currents. The
exact pair therefore separates these meters for every common origin c.
"""

from functools import lru_cache
from itertools import product
import json

import numpy as np
from scipy import sparse
from scipy.linalg import qr
import sympy as sp

from . import observability
from .. import live_radar, live_source
from ..source import certificate as cert

AXES = tuple(a + b for a, b in product("XYZ", repeat=2))
Q = sp.Rational


def outward(center, width):
    return cert.enclosed(center, width)


def table(rows):
    pairs = [[tuple(map(Q, row[name])) for name in AXES] for row in rows]
    return tuple(
        sp.Matrix([[f(a, b) for a, b in row] for row in pairs])
        for f in (lambda a, b: (a + b) / 2, lambda a, b: (b - a) / 2)
    )


@lru_cache(maxsize=1)
def material():
    (formed, product_head), ferr, _ = live_source.formation_center()
    return observability.native_bank(formed, product_head, ferr)


@lru_cache(maxsize=1)
def ward_rows():
    """Check native operator identities and use the 82-bit inclusive pump."""
    src = live_source.source
    board = src.formation.federation(1)
    route = src.formation.region(board)
    h, _, b2 = src.system(board, route, -1)
    eye = sparse.eye(b2.shape[0], dtype=np.int64, format="csr")
    z = sparse.block_diag((eye, -eye), format="csr")
    c = h @ z - z @ h
    # i[40H,Z]=-40Y, and -[40H,[40H,Z]]= -1600Z+360 X (2b).
    assert (
        c - sparse.bmat([[None, 40 * eye], [-40 * eye, None]], format="csr")
    ).count_nonzero() == 0
    xb2 = sparse.bmat([[None, b2], [b2, None]], format="csr")
    assert ((c @ h - h @ c) + 1600 * z - 360 * xb2).count_nonzero() == 0
    hi, lo, colerr, grams = live_source.pump_center(True)
    indices = live_radar.station_rows(board, route[1:3]).ravel()
    raw = (hi[indices].astype(object) * (1 << 36) + lo[indices]).astype(object)
    re, im = np.split(raw, 2, axis=1)
    local = (sp.Matrix(re) + sp.I * sp.Matrix(im)) / (1 << 82)
    x = sp.Matrix([[0, 1], [1, 0]])
    b = sp.Matrix([[17, -1], [-1, 17]]) / 2
    op = sp.kronecker_product(b, x, sp.eye(4))
    xb = (local.H * op * local).applyfunc(sp.expand)
    v = []
    for axis in observability.AXES:
        out = live_source.rational_read(grams, axis, (2, 1))
        a = sp.re(sum(out[k, k] - out[4 + k, 4 + k] for k in range(4)))
        dy = 2 * sp.im(sum(out[k, 4 + k] for k in range(4)))  # minus <Y>
        cov = sp.re(sp.trace(xb * sp.kronecker_product(sp.diag(2, 1) / 3, axis)))
        v.append([a, dy, -a + Q(9, 20) * cov])
    # All eight input columns together have operator error <=3 times the
    # single-column bound; every fixed-marginal unit Pauli/4 has trace norm 1.
    guard = 6 * colerr + 9 * colerr**2
    return sp.Matrix(v).T, sp.Matrix([[g * guard] * 9 for g in (1, 1, Q(101, 20))])


def intersection(m, w):
    p, e = ward_rows()
    m, w = sp.Matrix(m), sp.Matrix(w)
    if (m.shape, w.shape) != ((3, 9), (3, 9)) or min(w) < 0:
        raise ValueError(
            "supply three native unscaled jet rows and their nonnegative widths"
        )
    lo = sp.Matrix(3, 9, lambda i, j: max(p[i, j] - e[i, j], m[i, j] - w[i, j]))
    hi = sp.Matrix(3, 9, lambda i, j: min(p[i, j] + e[i, j], m[i, j] + w[i, j]))
    if min(hi - lo) < 0:
        raise ValueError(
            "supplied Z jet and independent native Ward enclosures disagree"
        )
    return (hi + lo) / 2, (hi - lo) / 2


def round_center(m, width, bits=78):
    mid = m.applyfunc(lambda x: Q(round(x * 2**bits), 2**bits))
    total = width + (m - mid).applyfunc(abs)
    return mid, total.applyfunc(lambda e: Q(sp.ceiling(e * 2**bits), 2**bits))


def abs_m(m):
    return m.applyfunc(abs)


def inf_norm(m):
    return max(sum(abs(x) for x in m.row(i)) for i in range(m.rows))


def solve_weights(constraints, cerr, bank, berr, reference_index):
    rows = constraints * bank
    scales = sp.diag(*[1 / max(abs(x) for x in rows.row(i)) for i in range(rows.rows)])
    m = scales * constraints
    em = scales * cerr
    rows = m * bank
    errors = em * abs_m(bank) + (abs_m(m) + em) * berr
    a = np.array(rows, dtype=float)
    _, _, vh = np.linalg.svd(a, full_matrices=True)
    null = vh[rows.rows :]
    proposal = null.T @ null[:, reference_index]
    proposal /= np.linalg.norm(proposal)
    w = observability.dyadic_vector(proposal)
    cols = list(map(int, qr(a, pivoting=True)[2][: rows.rows]))
    inv = rows[:, cols].inv()
    adjust = inv * (rows * w)
    for j, x in zip(cols, adjust):
        w[j] -= x
    assert rows * w == sp.zeros(rows.rows, 1)
    d = bank * w
    uncertainty = abs_m(inv) * errors[:, cols]
    defect = inf_norm(uncertainty)
    assert defect < Q(1, 10)
    noise = abs_m(inv) * (em * abs_m(d) + (abs_m(m) + em) * berr * abs_m(w))
    maxerr = max(noise) / (1 - defect)
    e = sp.zeros(9, 1)
    for j, v in zip(cols, noise + uncertainty * sp.ones(rows.rows, 1) * maxerr):
        e[j] = v
    return w, e, cols, defect


def witness(mid, widths):
    """Use native unscaled z,z',z'' (3×9 rational centers and absolute radii)."""
    measures = json.loads((observability.HERE / "current.json").read_text())["regions"]
    regional = json.loads((observability.HERE / "regional.json").read_text())["regions"]
    assert [r["radius"] for r in measures] == [r["radius"] for r in regional] == [1, 2]
    assert all(a["ports"] == b["ports"] for a, b in zip(measures, regional))
    front, fw = intersection(mid, widths)
    number, nw = table([r["derivatives"]["number"] for r in measures])
    targets, tw = table(
        [r["derivatives"]["energy"] for r in measures]
        + [r["regional_entropy_slopes"] for r in regional]
    )
    C, E = round_center(front.col_join(number), fw.col_join(nw))
    names, original, error = material()
    bank, berr = round_center(original, sp.ones(9) * error)
    proposal_name = "L1R2"
    w, we, columns, defect = solve_weights(C, E, bank, berr, names.index(proposal_name))
    # Positive weights select the first lottery, negative weights the second.
    # The unused probability in each prepares the unrotated paired reference.
    scale = max(sum(max(0, s * x + e) for x, e in zip(w, we)) for s in (1, -1)) + Q(
        1, 1000
    )
    # A dyadic upper scale makes all four frozen signed weights independently
    # reproducible as small numerator/denominator pairs.
    scale = Q(sp.ceiling(scale * 2**36), 2**36)
    w, we = w / scale, we / scale
    d = bank * w
    de = berr * abs_m(w) + (abs_m(bank) + berr) * we
    diffs = targets * d
    errs = tw * abs_m(d) + (abs_m(targets) + tw) * de
    # Subtract only certified zero rows. This second enclosure removes the
    # large part of a held covector parallel to the five matching conditions.
    dual = (targets * bank[:, columns]) * (C * bank[:, columns]).inv()
    sharp = targets - dual * C
    sharp_width = tw + abs_m(dual) * E
    assert sharp * d == diffs
    alternate = sharp_width * abs_m(d) + (abs_m(sharp) + sharp_width) * de
    errs = sp.Matrix([min(a, b) for a, b in zip(errs, alternate)])
    quantities = (
        "outward_energy_radius_1",
        "outward_energy_radius_2",
        "regional_entropy_slope_radius_1",
        "regional_entropy_slope_radius_2",
    )
    gaps = {name: outward(v, e) for name, v, e in zip(quantities, diffs, errs)}
    assert all(bounds[0] > 0 for bounds in gaps.values()), gaps

    # Exact probabilities are defined by the true five-row inverse on the
    # selected columns and by the frozen free weights on the other four.
    def lotteries(w, we):
        out = {}
        for side, sign in (("plus", 1), ("minus", -1)):
            lows = [max(0, sign * v - e) for v, e in zip(w, we)]
            highs = [max(0, sign * v + e) for v, e in zip(w, we)]
            assert sum(highs) <= 1
            out[side] = {
                n: [max(0.0, outward(x, 0)[0]), outward(y, 0)[1]]
                for n, x, y in zip(names, lows, highs)
            }
            out[side]["reference"] = [
                max(0.0, outward(1 - sum(highs), 0)[0]),
                outward(1 - sum(lows), 0)[1],
            ]
        return out

    # Finite dyadic weights constitute an independently replayable approximation.
    # Bounding its physical responses directly avoids treating near-zero number
    # current as exact invariance under unbounded energy-origin shifts.
    size = 1 << 52
    weights = {
        s: [int(sp.floor(max(0, sign * v) * size)) for v in w]
        for s, sign in (("plus", 1), ("minus", -1))
    }
    r = (sp.Matrix(weights["plus"]) - sp.Matrix(weights["minus"])) / size
    nominal = bank * r
    nomerr = berr * abs_m(r)
    fixed = abs_m(C * nominal) + E * abs_m(nominal) + (abs_m(C) + E) * nomerr
    gain = targets * nominal
    gerr = tw * abs_m(nominal) + (abs_m(targets) + tw) * nomerr
    dyadic_gaps = {n: outward(v, e) for n, v, e in zip(quantities, gain, gerr)}
    assert min(b[0] for b in dyadic_gaps.values()) > 0
    # Regional support is fixed by the rail and the inside/vacuum sectors.
    # For endpoints mixed into the faithful reference at common fraction t,
    # the entropy error about its modular slope is <= t^2/((1-t) lambda).
    # The two concave Taylor remainders lie in the same one-sided interval.
    floors = [
        min(map(Q, region["minimum_reference_eigenvalue_lower"])) for region in regional
    ]
    slopes = [diffs[i + 2] - errs[i + 2] for i in range(len(regional))]
    # Choose one dyadic fraction for both nested regions and both lotteries.
    t = sp.S.Half
    while t > min(floor * slope / 4 for floor, slope in zip(floors, slopes)):
        t /= 2
    finite = {}
    for i, (region, floor, slope) in enumerate(zip(regional, floors, slopes)):
        radius = t * t / ((1 - t) * floor)
        low = t * slope - radius
        assert low > 0
        finite[str(region["radius"])] = dict(
            common_admixture=float(t),
            finite_entropy_gap_nats=outward(t * diffs[i + 2], t * errs[i + 2] + radius),
        )
    labels = ["z", "dz", "ddz", "outward_number_radius_1", "outward_number_radius_2"]
    return dict(
        fixed_readings=labels,
        basis=list(AXES),
        native_arms=names,
        proposal_arm=proposal_name,
        correction_arms=[names[i] for i in columns],
        fixed_free_weights={
            name: [int(w[i].p), int(w[i].q)]
            for i, name in enumerate(names)
            if i not in columns
        },
        rank=5,
        inverse_defect_upper=outward(defect, 0)[1],
        ward_intersection_row_width_upper=[
            outward(max(fw.row(i)), 0)[1] for i in range(3)
        ],
        lottery_probability_intervals=lotteries(w, we),
        signed_responses=gaps,
        finite_regional_entropy=finite,
        dyadic_control=dict(
            denominator=size,
            positive=weights["plus"],
            negative=weights["minus"],
            positive_reference=size - sum(weights["plus"]),
            negative_reference=size - sum(weights["minus"]),
            fixed_response_absolute_upper={
                n: outward(x, 0)[1] for n, x in zip(labels, fixed)
            },
            signed_responses=dyadic_gaps,
        ),
    )
