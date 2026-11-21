"""Exact visible subspace and first response of each station to one written collar."""

from fractions import Fraction
from math import factorial, gcd

import numpy as np

from model.candidate_laws.scalar_channel.carrier import require
from model.experiments.cage import construction
from .joint_source import radar_source


def exact_chart(spec):
    car = construction.cell()
    c = car["C"]
    seams, wedges = c.shape
    levels = np.ones(seams, dtype=int)
    levels[spec["endpoints"]] += 1
    # 2 S^-1 H S for S = diag(I,sqrt(2) I). The induced metric is S*S.
    chart = np.block([[2 * np.diag(levels), 2 * c], [c.T, -2 * np.eye(wedges, dtype=int)]])
    weights = np.r_[np.ones(seams, dtype=int), np.full(wedges, 2)]
    root = np.sqrt(weights)
    actual = radar_source(spec, 0)
    require(np.max(abs(.5 * chart * root[:, None] / root - actual)) < 1e-13,
            "integer chart differs from the pulse generator")
    require(np.array_equal(weights[:, None] * chart, chart.T * weights), "wrong integer metric")
    return chart, weights


def marked_space(chart, weights, mark, stations):
    """Integer orthogonal Krylov basis; exact termination certifies maximality."""
    rows = [[(j, int(row[j])) for j in np.flatnonzero(row)] for row in chart]
    def advance(v):
        return [sum(c * v[j] for j, c in row) for row in rows]
    def dot(u, v):
        return sum(int(w) * a * b for w, a, b in zip(weights, u, v))

    basis = []
    v = [int(i == mark) for i in range(len(chart))]
    while any(v):
        require(len(basis) < len(chart), "nonterminating orthogonal recursion")
        basis.append(v)
        step = advance(v)
        # Project the new vector against the last two basis vectors. Symmetry makes
        # all earlier projections zero. Check exact orthogonality before accepting it.
        for q in basis[-2:]:
            norm, prod = dot(q, q), dot(q, step)
            residual = [norm * a - prod * b for a, b in zip(step, q)]
            scale = gcd(*residual) or 1
            step = [a // scale for a in residual]
        require(all(dot(q, step) == 0 for q in basis), "lost exact Krylov orthogonality")
        v = step
    moments = {}
    v = basis[0]
    for k in range(len(basis)):
        for a in stations:
            if a not in moments and v[a]:
                # Mark and enrolled stations are seam basis vectors, where S is identity.
                moments[a] = (k, Fraction(v[a], 2**k))
        v = advance(v)
    return len(basis), moments


def run(spec):
    chart, weights = exact_chart(spec)
    stations = spec["endpoints"]
    rank, moments = marked_space(chart, weights, spec["mark"], stations)
    single = []
    propagation = []
    for a in stations:
        if a not in moments:
            single.append(None)
            continue
        k, c = moments[a]
        single.append({"order_to_mark": k, "matrix_element": str(c),
                       "ramsey_order": 2*k + 1,
                       "minus_imaginary_coefficient": str(c*c / ((2*k+1)*factorial(k)**2))})
    for a in stations:
        row = []
        for b in stations:
            if a not in moments or b not in moments:
                row.append(None)
                continue
            ka, ca = moments[a]
            kb, cb = moments[b]
            n = ka + kb + 1
            row.append({"order": n, "coefficient_of_minus_it_power": str(ca*cb / factorial(n))})
        propagation.append(row)
    return {"dimension": len(chart), "mark_cyclic_dimension": rank,
            "unaffected_at_every_time_dimension": len(chart) - rank,
            "marked_collar_coherences_by_emitter": single, "station_propagation_amplitudes": propagation}
