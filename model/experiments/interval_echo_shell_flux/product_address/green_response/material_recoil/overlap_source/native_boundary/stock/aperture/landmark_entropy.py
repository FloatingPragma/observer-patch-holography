"""Certify observable frequencies at the entropy-stationary modular repair."""

from fractions import Fraction as Q
import math
import numpy as np
from scipy.optimize import least_squares
from scipy.special import softmax
from sympy import Matrix, diag

from . import collective_source as source


def exp_bounds(x):
    # Exact Taylor remainder at |t|<=1 is at most e/(41!) < 3/(41!).
    parts = 1
    while abs(x) > parts:
        parts *= 2
    t = x / parts
    total = term = Q(1)
    for k in range(1, 41):
        term *= t / k
        total += term
    err = Q(3, math.factorial(41))
    grid = 2**80
    lo = Q(math.floor((total - err) * grid), grid)
    hi = Q(math.ceil((total + err) * grid), grid)
    assert lo > 0
    while parts > 1:
        lo = Q(math.floor(lo**2 * grid), grid)
        hi = Q(math.ceil(hi**2 * grid), grid)
        parts //= 2
    return lo, hi


def as_fraction(x):
    return Q(int(x.p), int(x.q))


def solve(freeze_selective, progress=Q(1)):
    energy = source.tidal.lab.energy
    _, words, cells, _, codes, h, ids, n = energy.pair.setup()
    _, probs, _ = source.records()
    original = probs["original"]
    partial = np.isin(words[:, -1], source.PARTIAL)
    free = np.ones(n, bool)
    if freeze_selective:
        free[ids[partial]] = False
        assert np.array_equal(free[ids], ~partial)
    rows = np.vstack(
        (
            [cells == i for i in range(4)],
            words[:, :4].T,
            h.diagonal()[codes],
            energy.returned(words)[-1],
        )
    )
    target = rows @ original
    target[:4] += progress * (rows @ (probs["source"] - original))[:4]
    group_rows = energy.pair.group_matrix(rows, ids, n)
    sizes = np.bincount(ids)[free]
    total = Q(int(sizes.sum()), len(words))
    rhs = Matrix((target - group_rows[:, ~free].sum(1) * Q(1, len(words))).tolist())
    A = Matrix(group_rows[:, free])
    weight = diag(*[Q(1, int(g)) for g in sizes])
    averages = A * weight
    centered = averages - averages[:, 0] * Matrix.ones(1, A.cols)
    use = list(centered.T.rref()[1])
    scale = diag(*[1 / max(abs(x) for x in centered.row(i)) for i in use])
    features = scale * centered[use, :]
    desired = scale * (rhs[use, :] / total - averages[use, 0])
    V = np.array(features, float)
    y = np.array(desired, float).ravel()

    def law(theta):
        return softmax(np.log(sizes) + theta @ V)

    def jac(theta):
        w = law(theta)
        vbar = V @ w
        return (V * w) @ V.T - np.outer(vbar, vbar)

    result = least_squares(
        lambda t: V @ law(t) - y,
        np.zeros(len(use)),
        jac=jac,
        gtol=1e-14,
        ftol=1e-14,
        xtol=1e-14,
    )
    assert result.success
    theta = Matrix([Q(round(t * 10**10), 10**10) for t in result.x])
    exps = [exp_bounds(as_fraction(x)) for x in theta.T * features]
    denominator = [
        sum(int(n) * row[j] for n, row in zip(sizes, exps, strict=True)) for j in (0, 1)
    ]
    # Intervals for the conditional Gibbs mass of each free exchange component.
    q = [
        (int(n) * lo / denominator[1], int(n) * hi / denominator[0])
        for n, (lo, hi) in zip(sizes, exps, strict=True)
    ]
    start = Matrix(
        [
            Q(round((lo + hi) * total / (2 * int(n)) * 2**52), 2**52)
            for n, (lo, hi) in zip(sizes, q, strict=True)
        ]
    )
    basis = list(A.T.rref()[1])
    C = A[basis, :]
    witness = start + weight * C.T * (C * weight * C.T).inv() * (
        rhs[basis, :] - C * start
    )
    assert A * witness == rhs and min(witness) > 0
    group_p = np.full(n, Q(1, len(words)), object)
    group_p[free] = [as_fraction(v) for v in witness]
    p = group_p[ids]
    assert np.all(rows @ p == target) and sum(p) == 1
    # Fixed features give D(p*||q)<=D(p||q)<=chi²(p,q) for the
    # entropy maximizer p*. Pinsker and the triangle inequality bound
    # total variation between the feasible rational sampler and p*.
    chi = Q(0)
    for count, val, (lo, hi) in zip(sizes, group_p[free], q, strict=True):
        prob = int(count) * val / total
        chi += max(abs(prob - lo), abs(prob - hi)) ** 2 / lo
    bound = Q(1, 2**27)
    assert 2 * chi * total**2 < bound**2
    blank = words[:, -1] == 0
    da = (p - original) @ blank
    db = (p - original) @ partial

    def band(x):
        return [float(x - bound), float(x + bound)]

    return (
        dict(
            multipliers=list(map(str, theta)),
            group_weights=list(map(str, group_p)),
            total_variation_upper=str(bound),
            blank_increase=band(da),
            partial_increase=[0.0, 0.0] if freeze_selective else band(db),
            signed_variance=str(p @ (words[:, 1] - 1) ** 2),
            uniform_proposal_acceptance=str(1 / (len(p) * max(p))),
            ordinary_energy_readings=list(map(str, rows[8:] @ p)),
        ),
        p,
        (q, sizes, free, ids, total),
    )


def certify(freeze_selective):
    report, p, _ = solve(freeze_selective)
    return report, p
