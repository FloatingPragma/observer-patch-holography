"""Check the measured spherical branch on every enrolled SPARC galaxy."""

import argparse
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path

import numpy as np

from cosmology.rar_deep_regime.brightness_pairs import midpoint_halo as sky

HERE = Path(__file__).resolve().parent


@lru_cache(maxsize=24000)
def h_bounds(x):
    """Outward bounds for [1-(1-x)_+^(3/2)]/x, continuous at zero."""
    if x >= 1:
        return 1 / x, 1 / x
    den = 2**52
    a = 1 - x
    n = isqrt((a.numerator * den**2) // a.denominator)
    # h(s)=s+1/(1+s) increases with s=sqrt(1-x).
    u, v = F(n, den), F(n + 1, den)
    return u + 1 / (1 + u), v + 1 / (1 + v)


def shape(k, a, b):
    if k == 1:
        return F(1), F(1)
    # h decreases in its argument. Bound numerator and denominator
    # separately on the whole edge interval, without assuming a ratio slope.
    return k * h_bounds(k * k * b)[0] / h_bounds(a)[1], k * h_bounds(k * k * a)[
        1
    ] / h_bounds(b)[0]


def h_array(x):
    s = np.sqrt(np.maximum(0, 1 - x))
    return np.where(x >= 1, 1 / np.maximum(1, x), s + 1 / (1 + s))


def geometry(rings, bounds, inner):
    lo, hi = bounds
    lo = F(sky.sample.upper(lo) if inner else sky.sample.lower(lo))
    hi = F(sky.sample.lower(hi) if inner else sky.sample.upper(hi))
    bary = [([-1, 0, 0], -lo)]
    for j, mx in ((1, F(1)), (2, F(7, 5))):
        v = [F(1, 10), 0, 0]
        v[j] = -1
        bary.append((v, F(0)))
        v = [-mx, 0, 0]
        v[j] = 1
        bary.append((v, F(0)))
    r0 = rings[0][0]
    specs = [(v, rhs, None, 0) for v, rhs in bary]
    for r, slow, fast, parts, _ in rings:
        specs.append((parts, fast, r / r0, 1))
        specs.append(([-v for v in parts], -slow, r / r0, -1))
    return specs, [hi, hi, F(7, 5) * hi], rings[0][2]


def nodes(rings, n):
    r0 = rings[0][0]
    knots = sorted({F(0), *((r0 / r) ** 2 for r, *_ in rings)})
    return [
        a + (b - a) * i / n
        for a, b in zip(knots[:-1], knots[1:], strict=True)
        for i in range(n)
    ] + [knots[-1]]


def optimize(specs, cap, dark, grid):
    g = np.array(grid, float)
    rows = []
    scales = []
    for base, rhs, k, sign in specs:
        m = [a * b for a, b in zip(base, cap, strict=True)]
        scale = max(F(1), abs(rhs), *map(abs, m), dark * k if k is not None else F(0))
        cols = (
            sign * float(dark * k / scale) * h_array(float(k * k) * g) / h_array(g)
            if k is not None
            else np.zeros(len(g))
        )
        rows.append((np.r_[[float(a / scale) for a in m], cols], float(rhs / scale)))
        scales.append(scale)
    return sky.solve(rows), scales


def witness(res, specs, cap, dark, grid):
    params = [F(float(v)).limit_denominator(10**11) for v in res.x[:-1]]
    bary = [a * b for a, b in zip(params[:3], cap, strict=True)]
    source = [(t, p * dark) for t, p in zip(grid, params[3:], strict=True) if p]
    if min(params) < 0 or max(params) > 1:
        return None
    for base, rhs, k, sign in specs:
        val = sum(a * b for a, b in zip(bary, base, strict=True))
        if k is not None:
            val += sum(
                sign * w * shape(k, t, t)[1 if sign == 1 else 0] for t, w in source
            )
        if val > rhs:
            return None
    return {
        "baryonic_parameters": list(map(str, bary)),
        "edge_coordinate_and_squared_speed_at_first_ring": [
            list(map(str, pair)) for pair in source
        ],
    }


def separate(res, scales, specs, cap, dark, grid):
    weights = [
        F(float(-y)).limit_denominator(10**9) / u
        for y, u in zip(res.ineqlin.marginals, scales, strict=True)
    ]
    assert min(weights) >= 0
    base = [F(0)] * 3
    right = F(0)
    terms = {}
    for w, (row, rhs, k, sign) in zip(weights, specs, strict=True):
        if not w:
            continue
        right += w * rhs
        for j in range(3):
            base[j] += w * row[j]
        if sign:
            terms[k] = terms.get(k, F(0)) + w * sign * dark
    slack = sum(min(F(0), c * v) for c, v in zip(cap, base, strict=True)) - right
    if slack <= 0:
        return None
    # Baryonic squared speed is positive at the innermost ring throughout
    # the allowed box. The total normalized screen weight is at most one.
    queue = list(zip(grid[:-1], grid[1:], strict=True))
    least = slack
    steps = 0
    while queue:
        a, b = queue.pop()
        lower = sum(w * shape(k, a, b)[0 if w >= 0 else 1] for k, w in terms.items())
        margin = slack + min(F(0), lower)
        steps += 1
        if margin > 0:
            least = min(least, margin)
        elif steps > 24000 or b - a < F(1, 2**48):
            return None
        else:
            m = (a + b) / 2
            queue.extend([(a, m), (m, b)])
    return {
        "ring_and_baryon_multipliers": [
            [j, str(w)] for j, w in enumerate(weights) if w
        ],
        "certified_gap_lower_bound": str(least),
        "intervals_examined": steps,
    }


def test_galaxy(rings, bounds):
    for n in (2, 6, 18):
        grid = nodes(rings, n)
        # Every constructive screen has a strictly finite positive edge.
        inner = geometry(rings, bounds, True)
        res, _ = optimize(*inner, grid[1:])
        assert res.success
        fit = witness(res, *inner, grid[1:])
        if fit is not None:
            return {"outcome": "feasible", "subdivisions": n, **fit}
        outer = geometry(rings, bounds, False)
        dual, scales = optimize(*outer, grid)
        assert dual.success
        if dual.x[-1] < -1e-10:
            cert = separate(dual, scales, *outer, grid)
            if cert is not None:
                return {"outcome": "excluded", "subdivisions": n, **cert}
    return {"outcome": "unresolved"}


def one_screen(rings, bounds):
    inner = geometry(rings, bounds, True)
    specs, cap, dark = geometry(rings, bounds, False)
    checked = set()

    def fit(t):
        if t in checked or not t:
            return None
        checked.add(t)
        candidate, _ = optimize(*inner, [t])
        assert candidate.success
        return witness(candidate, *inner, [t])

    for t in reversed(nodes(rings, 1)[1:]):
        answer = fit(t)
        if answer is not None:
            return {"outcome": "feasible", **answer}
    boxes = []
    queue = [(F(0), F(1))]
    while queue:
        a, b = queue.pop()
        value = fit((a + b) / 2)
        if value is not None:
            return {"outcome": "feasible", **value}
        rows = []
        for base, rhs, k, sign in specs:
            coeff = sign * dark * shape(k, a, b)[0 if sign > 0 else 1] if sign else 0
            cols = [*(d * v for d, v in zip(base, cap, strict=True)), F(coeff)]
            scale = max(F(1), abs(rhs), *map(abs, cols))
            rows.append(([v / scale for v in cols], rhs / scale))
        res = sky.solve(rows)
        assert res.success
        if res.x[-1] < -1e-11:
            y = [F(float(-v)).limit_denominator(10**10) for v in res.ineqlin.marginals]
            assert min(y) >= 0
            terms = [(j, w) for j, w in enumerate(y) if w]
            left = [sum(w * rows[j][0][i] for j, w in terms) for i in range(4)]
            gap = sum(min(F(0), x) for x in left) - sum(
                w * rows[j][1] for j, w in terms
            )
            if gap > 0:
                boxes.append(
                    {
                        "ends": [str(a), str(b)],
                        "gap": str(gap),
                        "multipliers": [[j, str(w)] for j, w in terms],
                    }
                )
                continue
        if b - a < F(1, 2**44) or len(checked) > 5000:
            return {"outcome": "unresolved", "unresolved_edge": [str(a), str(b)]}
        mid = (a + b) / 2
        queue.extend([(a, mid), (mid, b)])
    return {"outcome": "excluded", "cover": boxes}


def chord_survey():
    from . import midpoint_bulk as bulk

    vertices, faces = bulk.base_chart()
    price = json.loads((HERE / "resolved_shells.json").read_bytes())
    kappa = bulk.shells.iv.mpf(price["conditional_entropy_slope"])
    word = [F(price["word_event_probabilities"][k]) for k in ("fresh", "shared")]
    radii = np.array([0.25, 0.5, 0.75, 1.0])
    cuts = np.sqrt(1 - radii**2)
    runs = []
    rng = np.random.default_rng(28122879)
    for level in range(5):
        x, w, mesh = bulk.sites(vertices, faces)
        collision = w @ w
        sums = np.zeros(len(radii))
        for i in range(0, len(x), 64):
            chord = np.sqrt(np.maximum(0, (1 - x[i : i + 64] @ x.T) / 2))
            chord[np.arange(len(chord)), np.arange(i, i + len(chord))] = 0
            weights = w[i : i + 64, None] * w * chord
            sums += [np.sum(weights * (chord >= cut - 2e-12)) for cut in cuts]
        ids = rng.choice(len(x), (2, 2**18), p=w)
        ids = ids[:, ids[0] != ids[1]]
        sep = np.linalg.norm(x[ids[0]] - x[ids[1]], axis=1) / 2
        coin = rng.random(len(sep))
        rows = []
        for radius, cut, amount in zip(
            radii, cuts, sums / (1 - collision), strict=True
        ):
            theory = 2 * (1 - cut**3) / 3
            lo = 2 * (1 - min(1, cut + mesh) ** 3) / 3 - mesh - collision
            hi = 2 * (1 - max(0, cut - mesh) ** 3) / 3 + mesh + collision
            assert lo <= amount <= hi
            rows.append(
                {
                    "radius_fraction": float(radius),
                    "population_active_fraction": float(amount),
                    "continuum": float(theory),
                    "refinement_envelope": [max(0.0, lo), min(1.0, hi)],
                    "address_gate_hits": int(
                        np.sum((sep >= cut - 2e-12) & (coin < sep))
                    ),
                    "word_read": bulk.native_read(F(float(amount)), kappa, word, rng),
                }
            )
        runs.append(
            {
                "stations": len(x),
                "mesh_displacement": mesh,
                "collision_probability": collision,
                "address_starts": 2**18,
                "distinct_address_starts": len(sep),
                "radii": rows,
            }
        )
        if level < 4:
            vertices, faces = bulk.geometry._subdivide(vertices, faces)
    return runs


def run():
    sky.sample.iv.dps = 36
    source = sky.HERE / "outer_requirement.json"
    selected = json.loads(source.read_bytes())["galaxies"]
    table1 = sky.sample.measure.DATA / "table1.dat"
    table2 = sky.sample.measure.DATA / "table2.dat"
    meta = {t[:11].strip(): t for t in table1.read_text().splitlines()}
    raw = table2.read_text().splitlines()
    answers = {}
    for name, row in selected.items():
        rings = sky.rings_in(raw, row["rings"])
        answers[name] = test_galaxy(rings, sky.sample.slope_box(meta[name]))
        answers[name]["one_screen"] = one_screen(
            rings, sky.sample.slope_box(meta[name])
        )
        answers[name]["catalogue_lines"] = [r[-1] for r in rings]
    files = [
        source,
        table1,
        table2,
        HERE / "CURVED_ORBITS.md",
        HERE / "curved_shells.json",
        HERE / "resolved_shells.json",
    ]
    from os.path import relpath

    return {
        "chord_survey": chord_survey(),
        "galaxies": answers,
        "groups": {
            k: [g for g, v in answers.items() if v["outcome"] == k]
            for k in ("feasible", "excluded", "unresolved")
        },
        "inputs": {relpath(p, HERE): sha256(p.read_bytes()).hexdigest() for p in files},
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    z = run()
    path = HERE / "curved_orbits.json"
    if args.write:
        path.write_text(json.dumps(z, indent=1) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), z)
    print({g: len(v) for g, v in z["groups"].items()}, z["groups"]["excluded"])
