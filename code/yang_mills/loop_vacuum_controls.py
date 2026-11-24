"""Haar SU(2) commutator spectrum and certified magnetic response on a cubic box."""

import argparse
from fractions import Fraction as F
from itertools import product
import json
from math import exp, log, sqrt
from pathlib import Path

import numpy as np
from scipy.special import roots_chebyu, roots_legendre


def cubic_faces(shape):
    faces = []
    for x in product(*(range(d) for d in shape)):
        for i, j in ((0, 1), (0, 2), (1, 2)):
            if x[i] + 1 >= shape[i] or x[j] + 1 >= shape[j]:
                continue

            def shift(k):
                return tuple(v + (a == k) for a, v in enumerate(x))

            faces.append(frozenset([(x, i), (x, j), (shift(i), j), (shift(j), i)]))
    return faces


def multiply(a, b):
    return np.r_[
        a[0] * b[0] - a[1:] @ b[1:],
        a[0] * b[1:] + b[0] * a[1:] + np.cross(a[1:], b[1:]),
    ]


def quadrature():
    # Independent group-word multiplication with the Haar radial and relative-axis measures.
    nodes, weights = roots_chebyu(6)
    cs, wc = roots_legendre(5)
    values = np.zeros(4)
    for i, j, k in product(range(6), range(6), range(5)):
        x, y, c = nodes[i], nodes[j], cs[k]
        u = np.array([x, sqrt(1 - x * x), 0.0, 0.0])
        v = np.array([y, c * sqrt(1 - y * y), sqrt((1 - y * y) * (1 - c * c)), 0.0])
        inverse = np.array([1, -1, -1, -1])
        uv = multiply(u, v)
        w = multiply(multiply(uv, u * inverse), v * inverse)[0]
        values += (
            weights[i]
            * weights[j]
            * wc[k]
            * (2 / np.pi) ** 2
            / 2
            * np.array([w, w * w, (x * x - 0.25) * (w - 0.25), uv[0] * x * y])
        )
    np.testing.assert_allclose(
        values, [1 / 4, 1 / 3, 1 / 16, 1 / 16], atol=2e-14, rtol=0
    )
    return values.tolist()


def response(count, coupling):
    # Casimir gap >= 3 on a simple bipartite cubic graph, and ||sum W|| <= count.
    radius = F(3, 4 * count)
    ratio = coupling / radius
    assert ratio < 1
    coefficient = F(7, 216)
    leading = coefficient * coupling**2
    tail = 2 * ratio**4 / (1 - ratio**2)
    lo, hi = leading - tail, leading + tail
    assert lo > 0
    perturbation = count * coupling
    gap = 3 - perturbation
    ground_distance = perturbation / (F(3, 2) - perturbation)
    variance_floor = F(1, 4) - 2 * ground_distance - hi**2
    variance_ceiling = min(F(1), F(1, 4) + 2 * ground_distance)
    first_moment_ceiling = F(6, 4) * (F(3, 4) + 2 * ground_distance)
    assert 0 < variance_floor <= variance_ceiling
    correlation = []
    for t in (0, 0.25, 0.5, 1, 2):
        lower = float(variance_floor) * exp(
            -t * float(first_moment_ceiling / variance_floor)
        )
        upper = float(variance_ceiling) * exp(-t * float(gap))
        assert 0 < lower <= upper
        correlation.append(dict(time=t, lower=lower, upper=upper))
    return dict(
        active_faces=count,
        coupling=str(coupling),
        mixed_derivative=str(coefficient),
        remainder_bound=str(tail),
        rectangle_mean_interval=[str(lo), str(hi)],
        electric_unit_gap_floor=str(gap),
        rectangle_variance_floor=str(variance_floor),
        rectangle_first_moment_ceiling=str(first_moment_ceiling),
        connected_rectangle_correlations=correlation,
    )


def arbitrary_coupling(j):
    # Product-Haar Harnack comparison for the full interacting vacuum at electric time 4.
    n = 2**j
    sites = n * n + 1
    edges = 3 * (sites - 1) * sites**2
    faces = 3 * sites * (sites - 1) ** 2
    a = F(1, n)
    g2 = F(4, 1 + j)
    alpha = g2 / (2 * a)
    b = 4 / g2**2
    log_density = -16 * float(b) * faces + 2 * edges * log(3 / 5)
    return dict(
        spacing=str(a),
        physical_side=n,
        links=edges,
        faces=faces,
        g_squared=str(g2),
        electric_clock=str(alpha),
        magnetic_ratio=str(b),
        log_density_floor=log_density,
        log_physical_gap_floor=log(3 * float(alpha)) + log_density,
        log_plaquette_variance_floor=log_density - log(4),
        log_physical_gap_ceiling=log(4 * float(alpha)) - log_density,
    )


def run():
    faces = cubic_faces((3, 2, 2))
    p, q = (
        faces[0],
        next(
            face for face in faces if ((1, 0, 0), 0) in face and ((1, 0, 0), 1) in face
        ),
    )
    perimeter = p ^ q
    assert len(faces) == 11 and len(perimeter) == 6
    tilings = [
        [i, j]
        for i in range(len(faces))
        for j in range(i)
        if faces[i] ^ faces[j] == perimeter
    ]
    assert len(tilings) == 1
    assert all(
        sum(sum(vertex[:axis]) for vertex, axis in face) % 2 == 1 for face in faces
    )
    shape = []
    for n in (1, 2, 4, 8, 16):
        # Radius, spacing and box length are fixed prospectively by this schedule.
        a = F(1, n)
        sites = 1 + n * n
        face_count = 3 * sites * (sites - 1) ** 2
        b = F(1, 512 * face_count)
        g2 = 2 / sqrt(float(b))
        mass = (3 - 1 / 512) * g2 / (2 * float(a))
        shape.append(
            dict(
                spacing=str(a),
                physical_side=n,
                faces=face_count,
                magnetic_ratio=str(b),
                g_squared=g2,
                physical_gap_floor=mass,
                log_connected_upper_at_unit_time=-mass,
            )
        )
    return dict(
        haar_quadrature=quadrature(),
        commutator_mean="1/4",
        commutator_variance="13/48",
        electric_weights={"one_spin_one_loop": "1/8", "two_spin_one_loops": "7/48"},
        electric_energies_per_loop_link=[2, 4],
        electric_square_correlation=[
            dict(time=t, value=exp(-8 * t) / 8 + 7 * exp(-16 * t) / 48)
            for t in (0, 0.25, 0.5, 1, 2)
        ],
        rectangle_tilings=tilings,
        cubic_extent=[3, 2, 2],
        magnetic=[response(2, F(1, 128)), response(11, F(1, 4096))],
        switch_off_mean_for_two_face_source=[0, 0],
        cofinal_strong_coupling=shape,
        running_coupling_comparison=[arbitrary_coupling(j) for j in range(5)],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.output.open("x") as f:
        json.dump(run(), f, indent=2, allow_nan=False)
        f.write("\n")
