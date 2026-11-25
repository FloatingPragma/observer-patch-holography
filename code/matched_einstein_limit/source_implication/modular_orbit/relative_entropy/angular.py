"""A time-stationary quadrupole transfer and its complete regional entropy remainder."""

from itertools import pairwise

import numpy as np
from mpmath import iv
from scipy.linalg import expm
from scipy.sparse.csgraph import shortest_path

from model.experiments.interval_echo_angular import modular
from model.experiments.interval_echo_einstein_station import cap

from ..angular_energy.certificate import filter_bounds, reference_matrix
from ..angular_energy.stationary import waiting_protocol


def probability_interval(distances):
    coefficients = filter_bounds()
    low, high = [reference_matrix(side, coefficients) for side in (False, True)]
    dot = [iv.mpf(1), 1 / iv.sqrt(5), -1 / iv.sqrt(5), iv.mpf(-1)]
    summed = iv.mpf(0)
    for i, j in np.ndindex(distances.shape):
        summed += dot[int(distances[i, j])] * iv.mpf(
            [int(low[i + 1, j + 1]), int(high[i + 1, j + 1])]
        )
    return summed / (12 * iv.mpf([int(np.trace(low)), int(np.trace(high))]))


def local_cost(a, b):
    """Sum nonnegative spectral terms for the trace Bregman divergence."""
    values, ua = np.linalg.eigh(a)
    reference, ub = np.linalg.eigh(b)
    assert min(values) > 0 and min(reference) > 0
    t = (values[:, None] - reference[None, :]) / reference
    scalar = (1 + t) * np.log1p(t) - t
    near = abs(t) < 0.0002
    small = t[near]
    scalar[near] = small**2 * (0.5 + small * (-1 / 6 + small * (1 / 12 - small / 20)))
    return float(np.sum(abs(ua.conj().T @ ub) ** 2 * reference * scalar))


def mode_meter(vectors):
    """Route three orthogonal real modes to three distinct terminal occupation words."""
    graph = cap.cell().neighbours

    def extend(path):
        if len(path) == len(graph):
            return path
        for j in sorted(graph[path[-1]]):
            if j not in path:
                result = extend(path + [j])
                if result:
                    return result
        return []

    path = extend([0])
    assert len(path) == 12 and all(y in graph[x] for x, y in pairwise(path))
    data = vectors[path].copy()
    u = np.eye(12)
    gate_error = 0.0
    appointments = 0
    single = cap.NativeFilter().binary[1:3, 1:3]
    for column in range(3):
        for j in range(11, column, -1):
            x, y = data[j - 1 : j + 1, column]
            if abs(y) < 1e-15:
                continue
            phi = np.arctan2(y, x)
            t = 2 * phi % (4 * np.pi)
            # The phase difference on the two ports converts i hopping to
            # a real rotation. Two local phases compensate the common edge
            # phase relative to all unaddressed positions.
            d = np.diag([1.0, 1j])
            native = np.exp(8.5j * t) * d @ expm(-1j * t * single) @ d.conj().T
            wanted = np.array([[np.cos(phi), np.sin(phi)], [-np.sin(phi), np.cos(phi)]])
            gate_error = max(gate_error, float(np.linalg.norm(native - wanted, ord=2)))
            u[j - 1 : j + 1] = wanted @ u[j - 1 : j + 1]
            data[j - 1 : j + 1] = wanted @ data[j - 1 : j + 1]
            appointments += 1
    expected = np.eye(12)[:, :3]
    assert np.max(abs(np.abs(u @ vectors[path]) - expected)) < 1e-12
    return dict(
        path=list(map(int, path)),
        rotations=appointments,
        worst_rotation_error=gate_error,
    )


def receipt():
    iv.dps = 35
    distance = shortest_path(
        cap.cell().adjacency(), directed=False, unweighted=True
    ).astype(int)
    bound = probability_interval(distance)
    parent = modular.common_reference()[0]
    native = cap.native_energy(
        [cap.seam_indices(e) for e in cap.cell().seams], np.ones(30)
    )[0]
    reference, _ = waiting_protocol(parent, native)
    tri = cap.chord_coordinates() / 2
    axial = tri @ (2 * tri[0])
    n = 2 * tri[0]
    _, axes = np.linalg.eigh(np.eye(3) - np.outer(n, n))
    orientations = np.column_stack((axial, tri @ axes[:, -2:]))
    transfer = np.outer(axial, axial) - (tri @ tri.T - np.outer(axial, axial)) / 2
    p = float(np.trace(tri.T @ reference[1:13, 1:13] @ tri) / 3)
    assert bound.a < p < bound.b
    np.testing.assert_allclose(
        reference[1:13, 1:13] @ tri, p * tri, rtol=1e-12, atol=1e-18
    )
    rng = np.random.default_rng(9150349)
    counts = []
    for x in (-0.5, 0, 0.5):
        probs = [p * (1 + x), p * (2 - x), 1 - 3 * p]
        counts.append(
            dict(
                fraction=x,
                trials=2**42,
                axial_transverse_other=rng.multinomial(2**42, probs).tolist(),
            )
        )
    levels = []
    t = 0.5
    f = (1 + t) * np.log1p(t) + (2 - t) * np.log1p(-t / 2)
    for radius in range(4):
        data = []
        for pole in range(12):
            subset = np.flatnonzero(distance[pole] <= radius)
            reduced = cap.CapRestriction(cap.WORDS, subset).reduce(reference)
            k = len(subset)
            b = reduced[1 : 1 + k, 1 : 1 + k]
            delta = p * transfer[np.ix_(subset, subset)]
            blank = reduced[:1, :1]
            dv = -np.trace(delta).real
            w, axes = np.linalg.eigh(b)
            derivative = -float(
                np.trace(delta @ ((axes * np.log(w)) @ axes.T)).real
                + dv * np.log(blank[0, 0])
            )
            cost = local_cost(b + t * delta, b) + local_cost(blank + t * dv, blank)
            data.append(
                dict(
                    pole=pole,
                    linear_entropy_slope=derivative,
                    relative_entropy=cost,
                    finite_entropy_difference=t * derivative - cost,
                )
            )
        levels.append(
            dict(
                hop_radius=radius,
                sites=k,
                full_information_fraction_home=data[0]["relative_entropy"] / (p * f),
                home=data[0],
                max_cost=max(v["relative_entropy"] for v in data),
            )
        )
    assert -0.231148 < levels[1]["home"]["linear_entropy_slope"] / p < -0.230975
    return dict(
        triplet_eigenweight_interval=str(bound),
        triplet_eigenweight=p,
        complete_information_at_half_transfer=p * f,
        mode_meter=mode_meter(orientations),
        sampled_classification=counts,
        nested_regions=levels,
    )
