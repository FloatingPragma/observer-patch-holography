"""Test unrestricted radial transport of a fixed inventory at ordinary Newtonian G."""

import argparse
from collections import defaultdict
from fractions import Fraction as Q
from hashlib import sha256
import json

import numpy as np

from . import kernel_budget as budget
from . import well_certificate as dec

HERE = budget.HERE


def enclosed(r, h, a, b):
    if h == 0:
        return max(0, min(b, r) ** 2 - a * a) / (b * b - a * a) if r >= a else 0
    full = max(0, min(b, max(0, r - h)) ** 2 - a * a)
    lo = max(a, abs(r - h))
    hi = min(b, r + h)
    partial = (
        0
        if hi <= lo
        else ((hi - lo) * r * r - ((hi - h) ** 3 - (lo - h) ** 3) / 3) / (2 * h)
    )
    return (full + partial) / (b * b - a * a)


def cdf_matrix(radius, offsets, chart):
    edges, weights = chart
    out = np.zeros((len(radius), len(offsets)))
    for i, r in enumerate(radius):
        a, b = edges[:-1], edges[1:]
        h = offsets[1:, None]
        full = np.maximum(0, np.minimum(b, np.maximum(0, r - h)) ** 2 - a * a)
        lo = np.maximum(a, abs(r - h))
        hi = np.minimum(b, r + h)
        width = np.maximum(0, hi - lo)
        midpoint = (lo + hi) / 2
        partial = width * (r * r - (midpoint - h) ** 2 - width**2 / 12) / (2 * h)
        out[i, 1:] = ((full + partial) / (b * b - a * a)) @ weights
        out[i, 0] = np.clip((r * r - a * a) / (b * b - a * a), 0, 1) @ weights
    if np.min(out) < -1e-8 or np.max(out) > 1 + 1e-8:
        raise ArithmeticError("Invalid enclosed fraction")
    return np.clip(out, 0, 1)


def mass_intervals(rows, mass):
    g = dec.decimal(budget.local.profiles.galaxies.G) * 10**9
    return [
        (r, low * r / (g * dec.decimal(mass)), high * r / (g * dec.decimal(mass)))
        for r, low, high in budget.exact.intervals(rows)
    ]


def cubic_events(radius, edges, weights):
    """Coefficients of h times enclosed fraction, as global polynomials in h."""
    events = defaultdict(lambda: [Q(0)] * 4)
    for a, b, w in zip(edges[:-1], edges[1:], weights, strict=True):
        knots = sorted(
            {Q(0), radius, abs(radius - a), abs(radius - b), radius + a, radius + b}
        )
        previous = [Q(0)] * 4
        for lo, hi in zip(knots, knots[1:]):
            step = (hi - lo) / 3
            pts = [lo + step * i for i in range(4)]
            f = [w * h * enclosed(radius, h, a, b) for h in pts]
            # Newton interpolation at four equal spacings, in global h.
            a3 = (f[3] - 3 * f[2] + 3 * f[1] - f[0]) / (6 * step**3)
            a2 = (f[2] - 2 * f[1] + f[0]) / (2 * step**2) - a3 * sum(pts[:3])
            a1 = (
                (f[1] - f[0]) / step
                - a3 * (pts[0] ** 2 + pts[0] * pts[1] + pts[1] ** 2)
                - a2 * (pts[0] + pts[1])
            )
            a0 = f[0] - lo * (a1 + lo * (a2 + lo * a3))
            current = [a0, a1, a2, a3]
            for i in range(4):
                events[lo][i] += current[i] - previous[i]
            previous = current
        for i in range(4):
            events[knots[-1]][i] -= previous[i]
    return events


def max_ratio(c, a, b):
    """Bernstein controls for a cubic divided by positive h on one interval."""

    def value(t):
        return c[0] + t * (c[1] + t * (c[2] + t * c[3]))

    def slope(t):
        return c[1] + t * (2 * c[2] + 3 * t * c[3])

    ta, tb = value(a), value(b)
    step = (b - a) / 3
    control = (ta, ta + step * slope(a), tb - step * slope(b), tb)
    assert a > 0 or ta == 0
    return max(x / (a + step * i) for i, x in enumerate(control) if a + step * i > 0)


def certify(rows, chart, mass, ray):
    weights = [Q(x) for x in ray]
    intervals = mass_intervals(rows, mass)
    target = sum(
        y * (lo if y > 0 else hi)
        for y, (_, lo, hi) in zip(weights, intervals, strict=True)
    )
    assert target > 0
    radii, seed = [[dec.decimal(x) for x in col] for col in chart]
    combined = defaultdict(lambda: [Q(0)] * 4)
    for y, (r, _, _) in zip(weights, intervals, strict=True):
        if not y:
            continue
        for h, poly in cubic_events(r, radii, seed).items():
            for i in range(4):
                combined[h][i] += y * poly[i]
    knots = sorted(combined)
    current = [Q(0)] * 4
    peak = Q(0)
    for a, b in zip(knots, knots[1:]):
        current = [v + delta for v, delta in zip(current, combined[a], strict=True)]
        # Subdivide to keep interval relaxation tight; every possible offset is covered.
        points = [a + (b - a) * k / 4 for k in range(5)]
        peak = max(
            peak, *(max_ratio(current, x, y) for x, y in zip(points, points[1:]))
        )
    assert all(x + y == 0 for x, y in zip(current, combined[knots[-1]], strict=True))
    return dict(
        impossible_at_any_mass=(peak == 0),
        ratio_lower=float(target / peak) if peak else None,
        covered_intervals=4 * (len(knots) - 1),
    )


def evaluate(rows, chart, mass, caps):
    intervals = mass_intervals(rows, mass)
    lo = np.array([float(lower) for _, lower, _ in intervals])
    hi = np.array([float(u) for _, _, u in intervals])
    # This uses only positivity and monotonicity of spherical enclosed mass.
    for j, (_, _, upper) in enumerate(intervals):
        if upper < 0:
            return dict(impossible_at_any_mass=True, negative_upper_ring=j)
        for i, (_, lower, _) in enumerate(intervals[: j + 1]):
            if lower > upper:
                return dict(impossible_at_any_mass=True, nonmonotone_pair=[i, j])
    r = np.array(rows["radius_kpc"])
    knots = np.unique(
        np.r_[0, r, abs(r[:, None] - chart[0]).ravel(), (r[:, None] + chart[0]).ravel()]
    )
    grid = np.unique(
        np.concatenate(
            [knots]
            + [knots[:-1] + fraction * np.diff(knots) for fraction in (0.25, 0.5, 0.75)]
        )
    )
    matrix = cdf_matrix(r, grid, chart)
    eps = np.minimum(1e-7, (hi - lo) / 100)
    fit = budget.local_inverse(matrix, lo + eps, hi - eps, False)
    proof = certify(rows, chart, mass, fit["dual_ray"])
    witness = None
    if fit["feasible"]:
        xs = grid[fit["knot_indices"]]
        ds = [dec.decimal(x) for x in xs]
        amps = [dec.decimal(x) for x in fit["weights"]]
        edges, seed = [[dec.decimal(x) for x in col] for col in chart]
        margin = []
        for R, lower, upper in intervals:
            val = sum(
                amt
                * sum(
                    w * enclosed(R, h, a, b)
                    for a, b, w in zip(edges[:-1], edges[1:], seed, strict=True)
                )
                for amt, h in zip(amps, ds, strict=True)
            )
            margin.append(min(val - lower, upper - val))
        assert min(margin) > 0
        witness = dict(
            offset_kpc=xs.tolist(),
            mass_ratios=fit["weights"],
            total_ratio=float(sum(amps)),
            verified_rings=len(margin),
            min_margin_ratio=float(min(margin)),
        )
        assert (
            not proof["impossible_at_any_mass"]
            and proof["ratio_lower"] <= witness["total_ratio"]
        )
    verdict = {}
    for k, cap in caps.items():
        if proof["impossible_at_any_mass"] or proof["ratio_lower"] > cap:
            verdict[k] = "excluded"
        elif witness is not None and witness["total_ratio"] <= cap:
            verdict[k] = "admitted"
        else:
            verdict[k] = "unresolved"
    return dict(**proof, dual_ray=fit["dual_ray"], witness=witness, verdict=verdict)


def run():
    inputs, charts, _, masses = budget.load_inputs()
    mixture = HERE / "mixture.json"
    caps = dict(
        trail_32=budget.local.profiles.CAP,
        cosmic_95=json.loads(mixture.read_bytes())["highest_galaxy_cap_at_95"],
    )
    out = {}
    for n, rows in inputs.items():
        panel = evaluate(rows, charts[n], masses[n], caps)
        if panel["impossible_at_any_mass"]:
            panel["verdict"] = dict.fromkeys(caps, "excluded")
        out[n] = panel
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                mixture,
                HERE / "source_profile.json",
                HERE / "KERNEL_BUDGET.md",
                budget.local.profiles.maps.HERE / "leroy_table7.dat",
                budget.local.profiles.maps.HERE / "leroy_sample.dat",
                budget.local.profiles.galaxies.sparc.DATA / "table1.dat",
            )
        },
        G_kpc_kms2_per_solar_mass=budget.local.profiles.galaxies.G,
        caps=caps,
        individual=out,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    out = run()
    path = HERE / "newton_transport.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), out)
    for k, v in out["individual"].items():
        print(
            k,
            v["verdict"],
            v.get("ratio_lower"),
            (v.get("witness") or {}).get("total_ratio"),
        )
