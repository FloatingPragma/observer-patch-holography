"""Time reflected geodesics in an ultrastatic curved surface times a line."""
from itertools import combinations, product
import math
import numpy as np


def distance(a, b, curvature):
    if curvature == 0:
        return float(np.linalg.norm(a - b))
    p, q = np.linalg.norm(a[:2]), np.linalg.norm(b[:2])
    cosine = np.dot(a[:2], b[:2]) / (p*q) if p*q else 1.
    # Half-angle cosine law in the curved factor, stable for small separations.
    if curvature > 0:
        half = math.sin((p-q)/2)**2 + math.sin(p)*math.sin(q)*(1-cosine)/2
        arc = 2 * math.asin(math.sqrt(max(0., half)))
    else:
        half = math.sinh((p-q)/2)**2 + math.sinh(p)*math.sinh(q)*(1-cosine)/2
        arc = 2 * math.asinh(math.sqrt(max(0., half)))
    return math.hypot(arc, a[2] - b[2])


def interval_product(a, b):
    values = [x*y for x, y in product(a, b)]
    return min(values), max(values)


def decode(radial_ticks, panels, tick):
    unit = tick / 2
    radii = np.array(radial_ticks) * unit
    noise = tick / 4
    rows = []
    for (i, j), clocks in zip(combinations(range(3), 2), panels, strict=True):
        # Round-trip clock rounding contributes at most a quarter tick to length.
        r, s = radii[i], radii[j]
        lengths = np.array(clocks) * unit
        x = r*r + s*s - np.mean(lengths**2)
        low = (r-noise)**2 + (s-noise)**2 - np.mean((lengths+noise)**2)
        high = (r+noise)**2 + (s+noise)**2 - np.mean(np.maximum(0., lengths-noise)**2)
        denom = ((r-noise)**2*(s-noise)**2, (r+noise)**2*(s+noise)**2)
        band = interval_product((3*low, 3*high), (1/denom[1], 1/denom[0]))
        rows.append(dict(axes=[i, j], chord_ticks=clocks, contraction=3*x/(r*r*s*s),
                         clock_interval=band))
    return dict(radial_ticks=radial_ticks, pairs=rows,
                scalar=2*sum(p["contraction"] for p in rows),
                scalar_clock_interval=[2*sum(p["clock_interval"][j] for p in rows) for j in [0,1]])


def run():
    q = np.array([1., 2., 2.]) / 3
    axes = np.eye(3) - 2*np.outer(q,q)
    zero = np.zeros(3)
    reports = []
    for ell in (.25, .125, .0625, .03125):
        tick = ell**5 / 4096
        for kind, k in (("product", -1), ("product", 0), ("product", 1), ("bent_flat", 0)):
            directions = axes if kind == "product" else np.eye(3)
            radii = (ell, 2*ell, 2*ell)
            ends = []
            for i, v in enumerate(directions):
                pair = []
                for sign in (1, -1):
                    direction = sign*v
                    if kind == "bent_flat" and i < 2:
                        direction = direction*math.cos(ell) + directions[2]*math.sin(ell)
                    pair.append(radii[i]*direction)
                ends.append(pair)

            def clock(a, b):
                return round(2*distance(a, b, k)/tick)

            radial = [clock(zero, pair[0]) for pair in ends]
            panels = [[clock(a, b) for a, b in product(ends[i], ends[j])]
                      for i, j in combinations(range(3), 2)]
            reading = decode(radial, panels, tick)
            diameters = [clock(*pair) for pair in ends]
            reports.append(dict(source=kind, sectional_curvature=k, short_radius=ell,
                clock_tick=tick, opposite_ticks=diameters,
                axial_shortfall=[(2*r-d)*tick/2 for r,d in zip(radial, diameters, strict=True)],
                expected_scalar=2*k, plane_limits=[k*axes[j,2]**2 for j in (2,1,0)],
                **reading))
    return reports
