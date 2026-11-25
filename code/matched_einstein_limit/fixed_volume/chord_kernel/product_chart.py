"""Separate sectional curvatures on a curved two-surface times a flat line."""

import math

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq


def product_exit(radius, r, unit, launch, tick, k):
    u, v = np.linalg.norm(unit[:2]), np.linalg.norm(launch[:2])
    mu = np.dot(unit[:2], launch[:2]) / (u * v) if u * v else 0.0
    sine = math.sin if k > 0 else math.sinh
    inverse = math.asin if k > 0 else math.asinh

    def q(distance):
        return 2 * sine(distance / 2) ** 2

    initial = q(r * u)

    def exterior(n):
        t = n * tick
        projected = (
            initial + (1 - k * initial) * q(t * v) + sine(r * u) * sine(t * v) * mu
        )
        surface_distance = 2 * inverse(math.sqrt(max(projected, 0) / 2))
        return surface_distance**2 + (r * unit[2] + t * launch[2]) ** 2 >= radius**2

    low, high = 0, math.ceil(3 * radius / tick)
    assert not exterior(low) and exterior(high)
    while high - low > 1:
        mid = (low + high) // 2
        if exterior(mid):
            high = mid
        else:
            low = mid
    return high


def survey(radius, tick, curvature):
    # Curved factor has unit positive or negative sectional curvature.
    s = np.sin if curvature > 0 else np.sinh
    u, weights = leggauss(32)
    u, weights = (u + 1) / 2, weights / 2
    transverse = np.sqrt(1 - u * u)

    def volume(r):
        a = r * transverse
        # The slice perpendicular to the independent line is a 2D geodesic disk.
        return 8 * np.pi * r * np.dot(weights, s(a / 2) ** 2)

    target = 4 * np.pi * radius**3 / 3
    aperture = brentq(
        lambda r: volume(r) - target, 0.5 * radius, 1.5 * radius, xtol=1e-15
    )
    angle = aperture * transverse
    area = 4 * np.pi * aperture**2 * np.dot(weights, s(angle) / angle)
    area_g = -(area - 4 * np.pi * radius**2) / (4 * np.pi * radius**4 / 15)
    r = radius / 2
    tilted = np.array([1.0, 2.0, 3.0])
    tilted = np.eye(3) - 2 * np.outer(tilted, tilted) / (tilted @ tilted)
    frames = []
    for name, basis in (("aligned", np.eye(3)), ("tilted", tilted)):
        readings = []
        for i, j in ((0, 1), (0, 2), (1, 2)):
            normal = np.cross(basis[i], basis[j])
            exits = [
                product_exit(radius, r, basis[i], direction * basis[j], tick, curvature)
                for direction in (1, -1)
            ]
            d2 = math.prod(exits) * tick**2
            flat = radius**2 - r**2
            readings.append(
                dict(
                    plane=[i, j],
                    exit_ticks=exits,
                    expected_sectional=curvature * normal[2] ** 2,
                    inferred_sectional=3 * (d2 - flat) / (r * r * flat),
                )
            )
        frames.append(
            dict(
                frame=name,
                planes=readings,
                clock_contraction=sum(p["inferred_sectional"] for p in readings),
            )
        )
    return dict(
        curved_factor=curvature,
        radius=radius,
        tick=tick,
        launch_ticks=round(r / tick),
        equal_volume_aperture=aperture,
        equal_volume_area_contraction=area_g,
        frames=frames,
    )
