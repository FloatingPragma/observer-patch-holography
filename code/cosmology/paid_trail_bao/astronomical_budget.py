"""Certify the cold-source budget from radial BAO and galaxy ages, with a free ruler."""

import argparse
from dataclasses import dataclass
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction as F
import hashlib
import json
from math import erf, sqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE
HASHES = {
    "mean": "9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585",
    "cov": "252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509",
}
# Loubser, arXiv:2511.02730, equations 4–6: z, H, statistical and systematic errors in km/s/Mpc.
CLOCKS = [
    (".46", "88.48", ".57", "12.32"),
    (".67", "119.45", "6.39", "16.64"),
    (".83", "108.28", "10.07", "15.08"),
]
BARYON_MAX = F(".02383")
NEUTRINO_MAX = F(".003")
RADIATION_MAX = F(".00005")
K = Decimal("2.774")
QUANTUM = Decimal("0.0000001")


def edges(mean, variance):
    with localcontext() as ctx:
        ctx.prec = 55
        width = K * variance.sqrt()
        low = (mean - width).quantize(QUANTUM, rounding=ROUND_FLOOR)
        high = (mean + width).quantize(QUANTUM, rounding=ROUND_CEILING)
    if low <= 0:
        raise ValueError(
            "Squaring or inverting the data needs a positive lower endpoint"
        )
    return F(low), F(high)


@dataclass(frozen=True)
class Inequality:
    # coefficients dot [omega_m, omega_k, omega_v, omega_r, (c/100r_d)^2] <= limit
    coefficients: tuple
    limit: F


def constraints():
    files = {}
    for suffix, digest in HASHES.items():
        raw = (DATA / f"desi_gaussian_bao_ALL_GCcomb_{suffix}.txt").read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError(f"Published DESI {suffix} bytes changed")
        files[suffix] = [
            line.split()
            for line in raw.decode().splitlines()
            if line.strip() and not line.startswith("#")
        ]
    if (
        len(files["mean"]) != 13
        or any(len(row) != 13 for row in files["cov"])
        or len(files["cov"]) != 13
    ):
        raise ValueError("Unexpected DESI dimensionality")
    rows = {}
    for i, (z, value, kind) in enumerate(files["mean"]):
        if kind != "DH_over_rs":
            continue
        low, high = edges(Decimal(value), Decimal(files["cov"][i][i]))
        x = 1 + F(z)
        prefix = f"BAO {float(z):g}"
        rows[prefix + " upper H"] = Inequality((x**3, x**2, 1, x**4, -1 / low**2), 0)
        rows[prefix + " lower H"] = Inequality(
            (-(x**3), -(x**2), -1, -(x**4), 1 / high**2), 0
        )
    if len(rows) != 12:
        raise ValueError("Expected six radial BAO coordinates")
    for z, mean, stat, syst in CLOCKS:
        low, high = edges(Decimal(mean), Decimal(stat) ** 2 + Decimal(syst) ** 2)
        x = 1 + F(z)
        prefix = f"CC {float(z):g}"
        rows[prefix + " upper H"] = Inequality(
            (x**3, x**2, 1, x**4, 0), (high / 100) ** 2
        )
        rows[prefix + " lower H"] = Inequality(
            (-(x**3), -(x**2), -1, -(x**4), 0), -((low / 100) ** 2)
        )
    for idx, name in [
        (0, "matter"),
        (2, "vacuum"),
        (3, "radiation"),
        (4, "ruler scale"),
    ]:
        rows[name + " nonnegative"] = Inequality(
            tuple(-int(i == idx) for i in range(5)), 0
        )
    rows["radiation cap"] = Inequality((0, 0, 0, 1, 0), RADIATION_MAX)
    return rows


def solve(matrix, values):
    """Exact elimination for the five equalities supporting each extremum."""
    a = [[F(x) for x in row] + [F(b)] for row, b in zip(matrix, values, strict=True)]
    n = len(a)
    for col in range(n):
        pivot = next(i for i in range(col, n) if a[i][col])
        a[pivot], a[col] = a[col], a[pivot]
        a[col] = [v / a[col][col] for v in a[col]]
        for i in range(n):
            if i != col:
                a[i] = [u - a[i][col] * v for u, v in zip(a[i], a[col], strict=True)]
    return [row[-1] for row in a]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b, strict=True))


def certificate(label, flat=False, cap=None):
    rows = constraints()
    active = [
        "BAO 0.934 upper H",
        "BAO 2.33 lower H",
        "CC 0.67 lower H",
        "radiation cap",
    ]
    objective = (1, 0, 0, 0, 0)
    if flat:
        rows["flat"] = Inequality(
            (0, 1, 0, 0, 0), 0
        )  # equality; dual multiplier unrestricted
        active.append("flat")
    elif cap is None:
        active.append("vacuum nonnegative")
    else:
        objective = (0, 1, 0, 0, 0)
        rows["mass cap"] = Inequality(
            (1, 0, 0, 0, 0), (1 + cap) * BARYON_MAX + NEUTRINO_MAX
        )
        active.append("mass cap")
    matrix = [rows[key].coefficients for key in active]
    limits = [rows[key].limit for key in active]
    point = solve(matrix, limits)
    dual = solve(list(zip(*matrix)), objective)
    # A feasible point and matching lower bound prove global optimality.
    if any(dot(row.coefficients, point) > row.limit for row in rows.values()) or (
        flat and point[1]
    ):
        raise ValueError(f"Infeasible sky or calibration in {label}")
    if any(v > 0 and key != "flat" for v, key in zip(dual, active, strict=True)):
        raise ValueError(f"Invalid lower-bound multiplier in {label}")
    if [dot(dual, col) for col in zip(*matrix)] != list(objective) or dot(
        dual, limits
    ) != dot(objective, point):
        raise ValueError(f"Primal and dual bounds differ for {label}")
    if sum(point[:4]) <= 0 or point[-1] <= 0:
        raise ValueError(f"Nonphysical present rate or ruler in {label}")

    def rounded(value):
        return round(float(value), 10)

    return {
        "quantity_minimized": "omega_m" if cap is None else "omega_k",
        "minimum": rounded(dot(objective, point)),
        "source_fraction_floor" if cap is None else "source_fraction_cap": rounded(
            (point[0] - NEUTRINO_MAX) / BARYON_MAX - 1 if cap is None else cap
        ),
        "extremizing_parameters": dict(
            zip(
                ("omega_m", "omega_k", "omega_vac", "omega_rad", "ruler_scale"),
                map(rounded, point),
                strict=True,
            )
        ),
        "ruler_Mpc_at_extremum": rounded(2997.92458 / sqrt(float(point[-1]))),
        "inequality_multipliers": {
            name: rounded(v) for name, v in zip(active, dual, strict=True)
        },
    }


def check():
    return {
        "inputs": "Six DESI DR2 radial distances and three Loubser DESI DR1 stellar-age rates",
        "joint_gaussian_marginal_box_coverage_lower": round(
            1 - 9 * (1 - erf(float(K) / sqrt(2))), 7
        ),
        "ruler": "unrestricted positive constant",
        "flat": certificate("flat", flat=True),
        "free_curvature": certificate("free curvature"),
        "unit_ticket_maximum": certificate("unit ticket", cap=F(1)),
        "observed_ticket_fraction": certificate("observed tickets", cap=F(8, 11)),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--print",
        action="store_true",
        help="also print the complete reproduced receipt",
    )
    args = parser.parse_args()
    result = check()
    if result != json.loads((HERE / "astronomical_budget.json").read_text()):
        raise SystemExit("The reconstructed budget differs from its recorded receipt")
    print(
        json.dumps(result, indent=2)
        if args.print
        else "Four free-ruler extrema verified by exact arithmetic."
    )
