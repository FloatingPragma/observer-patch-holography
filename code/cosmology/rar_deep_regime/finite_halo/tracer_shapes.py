"""Check exact certificates for spatially local and transported record stock."""

import argparse
from dataclasses import dataclass
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import ceil, factorial, floor
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
CATALOGUE = HERE.parent / "data"
INPUTS = {
    "table1.dat": CATALOGUE / "table1.dat",
    "table2.dat": CATALOGUE / "table2.dat",
    "publication.json": ROOT / "model/experiments/interval_echo_clock/publication.json",
    "exterior_reader.json": ROOT
    / "model/candidate_laws/interval_echo_published_records/spherical_dust/exterior_reader.json",
}
RECEIPT = HERE / "orbit_certificates.json"
CASES = ((2, 0), (5, 0), (5, 5), (5, 10))
MODELS = ("ordinary", "counted_collar", "local", "isothermal_annulus")


@dataclass(frozen=True)
class Ring:
    line: int
    radius: Q
    speed: Q
    error: Q
    template: tuple


def counted_mass():
    publication = json.loads(INPUTS["publication.json"].read_text())
    (specimen,) = (
        s
        for s in publication["specimens"]
        if (s["cells"], s["width"], s["seed"]) == (80, 3, 31)
    )
    words = list(zip(specimen["moving_public"], specimen["home_public"], strict=True))
    ordinary = Q(sum(t for t, _ in words), len(words))
    joint = Q(sum(t * h for t, h in words), len(words))
    record = json.loads(INPUTS["exterior_reader.json"].read_text())
    inner, outer = map(Q, record["areal_radii"])
    tick = Q(record["proper_tick"])

    def read_mu(row):
        sent, received = row["emitted_inner_ticks"], row["received_outer_ticks"]
        ratio = Q(sent[1] - sent[0], received[1] - received[0])
        clock_mu = (1 - ratio**2) / (2 * (1 / inner - ratio**2 / outer))
        radii = row["released_clocks"]["1"]
        acc = {
            s: 2 * (radii[s]["radius_ticks"] * tick - inner) / Q(s) ** 2
            for s in ("1/16", "1/32")
        }
        return (clock_mu, -(inner**2) * (4 * acc["1/32"] - acc["1/16"]) / 3)

    observed = {
        name: read_mu(record["readouts"][name])
        for name in ("ordinary", "held", "light_alone_collar_blind")
    }
    gains = [mu / ordinary for mu in observed["ordinary"]]
    for gain, mu in zip(gains, observed["held"], strict=True):
        assert abs(mu - gain * (ordinary + joint)) < Q(1, 10**9)
    split = observed["light_alone_collar_blind"]
    assert abs(split[1] - split[0]) > Q(6, 100)
    print(
        f"Ordinary, joint stock {ordinary}, {joint}; clock and release gains {float(gains[0]):.10f}, {float(gains[1]):.10f}"
    )
    return (ordinary + joint) / ordinary


def sin_bounds(angle):
    if angle in (0, 90):
        return (Q(angle, 90),) * 2
    # Archimedean rational bounds for pi followed by alternating Taylor bounds.
    # SPARC angle precision is 0.1 degree, safely below 90 for all interior values.
    low, high = [angle * q / 180 for q in (Q(333, 106), Q(355, 113))]

    def series(x, terms):
        return sum(
            (-1) ** k * x ** (2 * k + 1) / factorial(2 * k + 1) for k in range(terms)
        )

    return series(low, 6), series(high, 7)


def geometry(row):
    d, err, inc, di = (Q(row[a:b]) for a, b in ((15, 21), (22, 27), (30, 34), (35, 39)))
    lower_i, upper_i = sin_bounds(inc)
    lower = (
        max(0, 1 - 2 * err / d) * (sin_bounds(max(0, inc - 2 * di))[0] / upper_i) ** 2
    )
    upper = (1 + 2 * err / d) * (sin_bounds(min(90, inc + 2 * di))[1] / lower_i) ** 2
    scale = 10**8
    return Q(floor(lower * scale), scale), Q(ceil(upper * scale), scale)


def read_catalogue():
    disks = {}
    geometries = {}
    for row in INPUTS["table1.dat"].read_text().splitlines():
        if int(row[112:115]) in (1, 2) and Q(row[30:34]) >= 30:
            disks[row[:11].strip()] = Q(row[71:76])
            geometries[row[:11].strip()] = geometry(row)
    galaxies = {}
    for line, row in enumerate(INPUTS["table2.dat"].read_text().splitlines(), 1):
        name = row[:11].strip()
        radius = Q(row[19:25])
        if name not in disks or radius < 3 * disks[name]:
            continue
        speeds = [Q(row[a:b]) for a, b in ((39, 45), (46, 52), (53, 59))]
        gas, disk, bulge = [v * abs(v) for v in speeds]
        if (
            min(
                gas + d * disk + b * bulge
                for d in (Q("0.1"), Q(1))
                for b in (Q("0.1"), Q("1.4"))
            )
            <= 0
        ):
            continue
        ring = Ring(line, radius, Q(row[26:32]), Q(row[33:38]), (gas, disk, bulge))
        galaxies.setdefault(name, []).append(ring)
    return galaxies, geometries


def constraints(rings, model, sigma, extra, geo, stock_factor):
    """Return all rational half-spaces A*x <= b, with public row labels."""
    rows = {}
    for r in rings:
        vector = list(r.template)
        if model == "isothermal_annulus":
            vector += [Q(1), 1 / r.radius]
        if model in ("ordinary", "counted_collar"):
            vector += [Q(0)]
        half = sigma * r.error + extra
        rows[f"{r.line}/upper"] = (vector, (r.speed + half) ** 2)
        rows[f"{r.line}/lower"] = ([-x for x in vector], -(max(0, r.speed - half) ** 2))
    # Tracer amplitudes and annular mass slope have the attractive sign.
    # Grant the annulus an unrestricted central integration constant.
    for i, kind in enumerate(
        ("gas", "disk", "bulge", "annular slope")[
            : 4 if model == "isothermal_annulus" else 3
        ]
    ):
        v = [Q(0)] * len(vector)
        v[i] = -1
        rows[kind] = (v, Q(0))
    if model == "isothermal_annulus":
        inner_radius = min(r.radius for r in rings)
        rows["enclosed mass at inner edge"] = ([0, 0, 0, -1, -1 / inner_radius], Q(0))
    if model in ("ordinary", "counted_collar"):
        factor = stock_factor if model == "counted_collar" else Q(1)
        for i, (name, lo, hi) in enumerate(
            (
                ("gas", Q("0.5"), Q(2)),
                ("disk", Q("0.1"), Q(1)),
                ("bulge", Q("0.1"), Q("1.4")),
            )
        ):
            a = [Q(0)] * 4
            a[i] = 1
            a[3] = -hi * factor
            rows[name + "/ceiling"] = (a, Q(0))
            a = [Q(0)] * 4
            a[i] = -1
            a[3] = lo * factor
            rows[name + "/floor"] = (a, Q(0))
        rows["geometry/ceiling"] = ([0, 0, 0, 1], geo[1])
        rows["geometry/floor"] = ([0, 0, 0, -1], -geo[0])
    return rows


def dot(a, b):
    return sum((x * y for x, y in zip(a, b, strict=True)), Q(0))


def outer_force_requirement(rings):
    # Both tracers fade on this pair. The unused bulge is exactly zero.
    by_line = {r.line: r for r in rings}
    inner, outer = by_line[642], by_line[656]
    assert inner.radius == Q("14.49") and outer.radius == Q("26.08")
    positive = [
        (o, i) for o, i in zip(outer.template, inner.template, strict=True) if i > 0
    ]
    assert all(
        i > 0 or o == i == 0
        for o, i in zip(outer.template, inner.template, strict=True)
    )
    ratio = max(o / i for o, i in positive)
    inner_edge = min(r.radius for r in rings)
    # For v_halo^2 = H + C/r, nonnegative enclosed mass at the first
    # ring gives C >= -H*r_min. The C coefficient in the contrast is
    # negative, so this lower bound gives the largest possible contrast.
    radial_coefficient = 1 / outer.radius - ratio / inner.radius
    assert radial_coefficient < 0
    slope_multiplier = 1 - ratio - inner_edge * radial_coefficient
    for extra in (0, 5, 10):
        high = inner.speed + 5 * inner.error + extra
        low = max(0, outer.speed - 5 * outer.error - extra)
        # An additional inward source at the inner ring can only lower
        # the luminous-template contribution available at the outer one.
        fraction = max(0, 1 - ratio * high**2 / low**2)
        slope = max(0, low**2 - ratio * high**2) / slope_multiplier
        print(
            f"NGC 1003, 5 sigma + {extra} km/s: at least {float(fraction):.6%} of outer radial acceleration beyond local templates; annular H >= {float(slope):.4f} (km/s)^2"
        )


def validate(rows, certificate):
    if "coefficients" in certificate:
        coeff = list(map(Q, certificate["coefficients"]))
        assert all(dot(a, coeff) <= b for a, b in rows.values())
        return "feasible"
    combo = {k: Q(v) for k, v in certificate["contradiction"].items()}
    assert combo and all(y > 0 for y in combo.values())
    cols = zip(*(rows[k][0] for k in combo), strict=True)
    assert all(dot(combo.values(), col) == 0 for col in cols)
    # Positive weighted sum of valid inequalities would say 0 <= negative.
    assert dot(combo.values(), (rows[k][1] for k in combo)) < 0
    return "excluded"


def exact_dual(rows, keys):
    # Solve A_support^T*y=0, sum(y)=1. A linear-program solver selects
    # support only. This verifier recovers the weights over the rationals.
    eqs = [
        list(map(Q, col)) + [Q(0)]
        for col in zip(*(rows[k][0] for k in keys), strict=True)
    ]
    eqs.append([Q(1)] * (len(keys) + 1))
    for c in range(len(keys)):
        pivot = next(i for i in range(c, len(eqs)) if eqs[i][c])
        eqs[c], eqs[pivot] = eqs[pivot], eqs[c]
        scale = eqs[c][c]
        eqs[c] = [x / scale for x in eqs[c]]
        for i in range(len(eqs)):
            if i != c:
                factor = eqs[i][c]
                eqs[i] = [x - factor * y for x, y in zip(eqs[i], eqs[c], strict=True)]
    assert all(not row[-1] for row in eqs[len(keys) :])
    return [eqs[i][-1] for i in range(len(keys))]


def discover(rows):
    # Optional regeneration. The normal replay uses only the Python stdlib.
    import numpy as np
    from scipy.optimize import linprog

    names = list(rows)
    a = np.array([rows[k][0] for k in names], dtype=float)
    b = np.array([rows[k][1] for k in names], dtype=float)
    n, d = a.shape
    scales = np.maximum(1, np.maximum(abs(a).max(axis=1), abs(b)))
    a, b = a / scales[:, None], b / scales
    column_scales = np.max(abs(a), axis=0)
    column_scales[column_scales == 0] = 1
    a = a / column_scales
    tolerances = {
        "dual_feasibility_tolerance": 1e-10,
        "primal_feasibility_tolerance": 1e-10,
    }
    interior = linprog(
        [0] * d + [-1],
        A_ub=np.c_[a, np.ones(n)],
        b_ub=b,
        bounds=[(None, None)] * d + [(None, 1)],
        method="highs",
        options=tolerances,
    )
    assert interior.success, interior.message
    if interior.x[-1] >= -1e-10:
        for limit in (10000, 10**8, 10**14):
            candidate = [
                Q(float(x)).limit_denominator(limit)
                for x in interior.x[:-1] / column_scales
            ]
            if all(dot(aq, candidate) <= bq for aq, bq in rows.values()):
                return {"coefficients": list(map(str, candidate))}
    dual = linprog(
        b,
        A_eq=np.vstack([a.T, np.ones(n)]),
        b_eq=[0] * d + [1],
        bounds=(0, None),
        method="highs",
        options=tolerances,
    )
    assert dual.success and dual.fun < 0, "No exact primal or dual support found"
    support = [
        name for name, weight in zip(names, dual.x, strict=True) if weight > 1e-11
    ]
    weights = exact_dual(rows, support)
    return {"contradiction": dict(zip(support, map(str, weights), strict=True))}


def run(regenerate):
    galaxies, geometries = read_catalogue()
    stock_factor = counted_mass()
    hashes = {
        name: sha256(path.read_bytes()).hexdigest() for name, path in INPUTS.items()
    }
    archive = (
        {"inputs": hashes, "cases": {}}
        if regenerate
        else json.loads(RECEIPT.read_text())
    )
    assert hashes == archive["inputs"]
    print(f"{len(galaxies)} galaxies, {sum(map(len, galaxies.values()))} outer rings")
    outer_force_requirement(galaxies["NGC1003"])
    expected = set()
    for model in MODELS:
        for sigma, extra in CASES:
            label = f"{model} / {sigma} sigma / {extra} km/s"
            expected.add(label)
            certificates = (
                archive["cases"].setdefault(label, {})
                if regenerate
                else archive["cases"][label]
            )
            excluded = []
            for name, rings in galaxies.items():
                rows = constraints(
                    rings, model, sigma, extra, geometries[name], stock_factor
                )
                if regenerate:
                    certificates[name] = discover(rows)
                if validate(rows, certificates[name]) == "excluded":
                    excluded.append(name)
            assert set(certificates) == set(galaxies)
            print(
                f"{label}: {len(excluded)} excluded, {len(galaxies) - len(excluded)} feasible"
            )
    assert expected == set(archive["cases"])
    if regenerate:
        RECEIPT.write_text(json.dumps(archive, sort_keys=True, indent=1) + "\n")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run without -O to enable the scientific checks.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--regenerate",
        action="store_true",
        help="rediscover all exact witnesses using SciPy",
    )
    run(parser.parse_args().regenerate)
