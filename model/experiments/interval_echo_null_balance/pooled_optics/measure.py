"""Exact charge-pooled coefficients and signed optical controls for a written-field cell."""

from fractions import Fraction
from itertools import combinations
import json

import numpy as np
from scipy import sparse
from scipy.optimize import linprog
from scipy.sparse.linalg import expm_multiply

from model.experiments.interval_echo_einstein_station.field_stress_curvature import (
    measure as field,
)


def polynomial_effects(patch, words, i):
    """Coefficients at noninitial endpoints from three powers of the integer generator 2H."""
    n = len(words)
    per_sign = []
    for sign in (-1, 1):
        patch.compensator(i, sign)
        pos, flight = field.charge_sector(patch, words, i, sign)
        h = (2 * flight).astype(np.int64)
        assert (h - 2 * flight).nnz == 0
        # Every power and entrywise product below is an integer. This row-norm bound
        # controls their intermediate products in the sparse multiplication.
        row = int(abs(h).sum(axis=1).max())
        assert 14 * row**4 * (128 * n) < 2**62
        one = h[:, :n]
        two = h @ one
        three = h @ two
        terms = []
        for j in range(12):
            if j == i:
                terms.append(None)
                continue
            a, b, c = (m[pos == j] for m in (one, two, three))
            terms.append(
                (a.T @ a, a.T @ b - b.T @ a, 3 * (b.T @ b) - 2 * (a.T @ c + c.T @ a))
            )
        per_sign.append(terms)
    edges = set(patch.edges)
    quartic = []
    for j in range(12):
        if i == j:
            quartic.append(None)
            continue
        a, b, c = [
            left + right
            for left, right in zip(per_sign[0][j], per_sign[1][j], strict=True)
        ]
        target = (
            64 * sparse.eye(n, dtype=np.int64, format="csr")
            if tuple(sorted((i, j))) in edges
            else 0
        )
        assert (a - target).nnz == 0 and b.nnz == 0
        quartic.append(c)
    return quartic


def sources(patch, words, head, ring, ops):
    """Unnormalized Gaussian-integer sources. All norms squared are 2*n."""
    orientation = np.where(words & patch.bit[head], -1, 1)
    twice_ring = (2 * ops[ring]).astype(np.int64)
    acted = twice_ring @ orientation
    assert not np.any(acted % 2)
    hit = acted // 2
    # v=(I-(1+i)h_ring) Z_e 1. The two states are 1 +/- i v.
    return [(1 + sign * hit, sign * (orientation - hit)) for sign in (-1, 1)]


def apply_coefficient(matrix, state):
    a, b = state
    return int(a @ (matrix @ a) + b @ (matrix @ b))


def currents(ops, state):
    a, b = state
    real = np.stack([(2 * h).astype(np.int64) @ a for h in ops])
    imag = np.stack([(2 * h).astype(np.int64) @ b for h in ops])
    n = int(a @ a + b @ b)
    return [
        [Fraction(int(value), 2 * n) for value in row]
        for row in -(real @ imag.T - imag @ real.T)
    ]


def raw_quartics(rows, states):
    values = np.zeros((2, 12, 12), dtype=np.int64)
    for i, els in enumerate(rows):
        for j, mat in enumerate(els):
            if i != j:
                for k, state in enumerate(states):
                    values[k, i, j] = apply_coefficient(mat, state)
    return values


def overlap4(patch, i, j, values):
    """Fourth-order coefficient of endpoint overlap, ignoring source-independent calibration."""
    edges = set(patch.edges)
    total = values[i, j] + values[j, i]
    for k in range(12):
        if k in (i, j):
            continue
        ai = tuple(sorted((i, k))) in edges
        aj = tuple(sorted((j, k))) in edges
        if ai == aj:
            total += min(values[i, k], values[j, k])
        elif ai:
            total += values[j, k]
        else:
            total += values[i, k]
    return total


def mixture_kinks(patch, ray, values, denom):
    """Exact switches of the smaller arrival probability in the equal-order bins."""
    i, j = ray
    edges = set(patch.edges)
    switches = []
    for k in range(12):
        if k in ray or (
            (tuple(sorted((i, k))) in edges) != (tuple(sorted((j, k))) in edges)
        ):
            continue
        low, high = (int(v[i, k] - v[j, k]) for v in values)
        if low * high < 0:
            switches.append(
                dict(
                    port=k,
                    positive_weight=str(Fraction(-low, high - low)),
                    downward_slope_jump=str(Fraction(abs(high - low), denom)),
                )
            )
    return switches


def residual_cost(distance, first, second, scale):
    """Solve the finite residual transport LP after removing its common unit cost."""
    difference = first - second
    sources = np.flatnonzero(difference > 0)
    sinks = np.flatnonzero(difference < 0)
    costs = distance[np.ix_(sources, sinks)] - 1
    supplies = np.kron(np.eye(len(sources)), np.ones((1, len(sinks))))
    demands = np.tile(np.eye(len(sinks)), (1, len(sources)))
    # Divide out the t^4 bin scale and the calibration t^2 cost scale to retain
    # small source-dependent transport routes. One redundant final demand is omitted.
    constraints = np.r_[supplies, demands[:-1]]
    marginals = np.r_[difference[sources], -difference[sinks[:-1]]]
    price_scale = np.max(abs(costs))
    answer = linprog(
        costs.ravel() / price_scale,
        A_eq=constraints,
        b_eq=marginals / scale,
        method="highs-ds",
    )
    if not answer.success:
        raise ArithmeticError(answer.message)
    return answer.fun * scale * price_scale


def finite_replay(patch, words, states, ray):
    """Direct unitary Born laws at all endpoints of the two selected injections."""
    times = (1 / 32, 1 / 64, 1 / 128, 1 / 256)
    kernels = np.zeros((len(times), 2, 2, 12))
    columns = np.array([a + 1j * b for a, b in states]).T / np.sqrt(2 * len(words))
    vacuum = np.zeros((len(times), 12, 12))
    for i in range(12):
        for sign in (-1, 1):
            patch.compensator(i, sign)
            ports, h = field.charge_sector(patch, words, i, sign)
            for step, time in enumerate(times):
                generator = (-1j * time * 2 / 3) * h
                quiet = np.zeros(len(ports), dtype=complex)
                quiet[: len(words)] = 1 / np.sqrt(len(words))
                evolved = expm_multiply(
                    generator, quiet, traceA=generator.diagonal().sum()
                )
                vacuum[step, i] += (
                    np.bincount(ports, abs(evolved) ** 2, minlength=12) / 2
                )
            if i not in ray:
                continue
            launch = ray.index(i)
            start = np.zeros((len(ports), 2), dtype=complex)
            start[: len(words)] = columns
            for step, time in enumerate(times):
                generator = (-1j * time) * h
                wave = expm_multiply(
                    generator, start, traceA=generator.diagonal().sum()
                )
                for end in range(12):
                    kernels[step, :, launch, end] += (
                        np.sum(abs(wave[ports == end]) ** 2, axis=0) / 2
                    )
    np.testing.assert_allclose(kernels.sum(axis=-1), 1, atol=5e-14, rtol=0)
    records = []
    np.testing.assert_allclose(vacuum.sum(axis=-1), 1, atol=5e-14, rtol=0)
    for time, kernel, quiet in zip(times, kernels, vacuum, strict=True):
        overlap = np.minimum(kernel[:, 0], kernel[:, 1]).sum(axis=-1)
        mixture = kernel.mean(axis=0)
        middle = np.minimum(mixture[0], mixture[1]).sum()
        d = np.abs(quiet[:, None, :] - quiet[None, :, :]).sum(axis=2) / 2
        correction = np.array(
            [
                residual_cost(d, first, second, time**4)
                for first, second in [*kernel, mixture]
            ]
        )
        curvature_delta = np.diff(overlap)[0] - (correction[1] - correction[0])
        curvature_gap = (
            middle - overlap.mean() - (correction[2] - correction[:2].mean())
        )
        records.append(
            dict(
                dwell=time,
                vacuum_dwell=time * 2 / 3,
                calibrated_curvature_difference_over_t4=float(
                    curvature_delta / (d[ray] * time**4)
                ),
                calibrated_curvature_midpoint_gap_over_t4=float(
                    curvature_gap / (d[ray] * time**4)
                ),
                difference_over_t4=float(np.diff(overlap)[0] / time**4),
                midpoint_gap_over_t4=float((middle - overlap.mean()) / time**4),
            )
        )
    return records


def run():
    patch = field.search.Patch()
    words = np.array(patch.field_component())
    ops = field.components(patch, words)
    head = (4, 11)
    ring = list(patch.faces).index((2, 4, 11))
    states = sources(patch, words, head, ring, ops)
    norm = 2 * len(words)
    assert all(a @ a + b @ b == norm for a, b in states)
    first, other = [currents(ops, s) for s in states]
    assert first == other
    current = max(
        (
            (abs(first[a][b]), a, b, first[a][b])
            for a, b in combinations(range(len(ops)), 2)
        )
    )
    rows = [polynomial_effects(patch, words, i) for i in range(12)]
    values = raw_quartics(rows, states)
    denom = 384 * norm
    result = []
    for i, j in patch.edges:
        neighbors = [
            {b if a == k else a for a, b in patch.edges if k in (a, b)} for k in (i, j)
        ]
        assert len(neighbors[0] & neighbors[1]) == 2
        low, high = [Fraction(int(overlap4(patch, i, j, k)), denom) for k in values]
        middle = Fraction(int(overlap4(patch, i, j, values.sum(axis=0))), 2 * denom)
        result.append(
            dict(
                ray=[i, j],
                low=str(low),
                high=str(high),
                middle=str(middle),
                difference=str(high - low),
                midpoint_gap=str(middle - (low + high) / 2),
            )
        )
    selected = (2, 4)
    # Mix each sign with a common uniform lottery over all legal sharp head words.
    # Since denom=384*2*n, the maximally mixed source contributes 2 Tr(C_num).
    uniform = np.zeros((12, 12), dtype=np.int64)
    for i, effects in enumerate(rows):
        for j, quartic in enumerate(effects):
            if i != j:
                uniform[i, j] = 2 * quartic.diagonal().sum()
    faithful = 15 * values + uniform
    faithful_denom = 16 * denom
    i, j = selected
    ends = [Fraction(int(overlap4(patch, i, j, k)), faithful_denom) for k in faithful]
    midpoint = Fraction(
        int(overlap4(patch, i, j, faithful.sum(axis=0))), 2 * faithful_denom
    )
    switches = mixture_kinks(patch, selected, faithful, faithful_denom)
    (switch,) = switches
    switch_p = Fraction(switch["positive_weight"])
    slope_drop = Fraction(switch["downward_slope_jump"])
    difference = ends[1] - ends[0]
    faithful_row = dict(
        common_uniform_weight="1/16",
        minimum_eigenvalue=str(Fraction(1, 16 * len(words))),
        arm_difference=str(difference),
        initial_slope=str(difference + (1 - switch_p) * slope_drop),
        final_slope=str(difference - switch_p * slope_drop),
        midpoint_gap=str(midpoint - sum(ends) / 2),
        switches=switches,
    )
    # Suppress the addressed ring pulse on the same two head preparations.
    orientation = np.where(words & patch.bit[head], -1, 1)
    control = [(np.ones_like(orientation), sign * orientation) for sign in (-1, 1)]
    controls = raw_quartics(rows, control)
    assert np.array_equal(controls[0], controls[1])
    return dict(
        selected_ray=selected,
        full_rank_lottery=faithful_row,
        selected_switches=mixture_kinks(patch, selected, values, denom),
        direct_unitary_overlap=finite_replay(patch, words, states, selected),
        suppressed_ring_all_quartic_contrasts_zero=bool(np.array_equal(*controls)),
        words=len(words),
        head=head,
        ring=patch.faces[ring],
        common_current=dict(
            face_a=patch.faces[current[1]],
            face_b=patch.faces[current[2]],
            value=str(current[3]),
        ),
        all_rays=result,
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
