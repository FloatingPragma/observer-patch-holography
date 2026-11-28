"""Resolved and unresolved ordinary tracers compared with the same outer orbits."""

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json

from mpmath import iv

from . import thermal_work as orbits
from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    occupied_equilibrium as experiment,
)

HERE = orbits.HERE
LAB = experiment.HERE


def empirical():
    # Rebuild the data intervals from the deposited speeds and mass model.
    inputs, *_ = orbits.budget.load_inputs()
    gal = inputs["NGC2841"]
    lo, hi = orbits.budget.bounds(gal)
    r = gal["radius_kpc"]
    idx = [min(range(len(r)), key=lambda k: abs(r[k] - 10.25 * j)) for j in range(2, 7)]
    rows = [
        dict(kpc=x, extra_squared_kms=[float(lower), float(upper)])
        for x, lower, upper in zip(r, lo, hi, strict=True)
    ]
    work = [
        orbits.integral(rows[a : b + 1]) for a, b in zip(idx[:-1], idx[1:], strict=True)
    ]
    return [float(r[j]) for j in idx], work


def extract(row, key):
    return [iv.mpf([int(a), int(b)]) / 10 ** row["decimal_places"] for a, b in row[key]]


def coefficients(data):
    rows = [data[tag] for tag in ("vacancy", "excess")]
    raw = [extract(s, "native_probability_bounds") for s in rows]
    slope = [extract(s, "mixed_probability_bounds") for s in rows]
    coeff = []
    for u in range(2, 6):
        v = u + 1

        def cross(i, j):
            return slope[i][v] * raw[j][u] - slope[i][u] * raw[j][v]

        coeff.append((cross(0, 0), cross(0, 1) + cross(1, 0), cross(1, 1)))
    return raw, slope, coeff


def work_interval(f, source, raw, slope, coeff):
    # Bernstein products cancel the shell-independent normalization before
    # enclosing a whole range of tracer fractions.
    vacancy = 1 - f
    p = [vacancy * a + f * b for a, b in zip(*raw, strict=True)]
    m = [vacancy * a + f * b for a, b in zip(*slope, strict=True)]
    out = []
    for u, (a, b, c) in enumerate(coeff, 2):
        num = vacancy**2 * a + f * vacancy * b + f**2 * c
        denominator = p[u + 1] * (p[u] + source * m[u])
        if denominator.a <= 0:
            return None
        term = source * num / denominator
        if term.a <= -1:
            return None
        out.append(-iv.log1p(term) / iv.log(iv.mpf(8) / 7))
    return out


def audit(data, delta, work):
    iv.dps = 65
    source = iv.mpf(delta)
    fixed = coefficients(data)
    obs = [iv.mpf(row) for row in work]
    counts = Counter()
    leaves = []
    unresolved = []

    def visit(left, right, depth):
        model = work_interval(iv.mpf([left, right]) / 2**depth, source, *fixed)
        if model is not None:
            for i, j in combinations(range(len(obs)), 2):
                determinant = obs[i] * model[j] - obs[j] * model[i]
                if determinant.a > 0 or determinant.b < 0:
                    label = f"{i + 2}_{i + 3} vs {j + 2}_{j + 3}"
                    counts[label] += 1
                    leaves.append(
                        [
                            left,
                            depth,
                            label,
                            experiment.flight.lab.interval_output(determinant),
                        ]
                    )
                    return
        if depth == 60:
            unresolved.append([left, depth])
            return
        visit(2 * left, left + right, depth + 1)
        visit(left + right, 2 * right, depth + 1)

    visit(0, 1, 0)
    coverage = sum((Fraction(1, 2**k) for _, k, *_ in leaves), Fraction(0))
    assert coverage + sum((Fraction(1, 2**k) for _, k in unresolved), Fraction(0)) == 1
    return dict(
        enclosed_fraction=str(coverage),
        excluded_every_mixture=not unresolved,
        disjoint_dyadic_cells=leaves,
        unresolved=unresolved,
        comparisons=dict(counts),
    )


def run():
    lab = json.loads((LAB / "occupied_equilibrium.json").read_bytes())["results"]
    old = json.loads((LAB / "pair_contact.json").read_bytes())
    prior = json.loads((HERE / "thermal_work.json").read_bytes())
    delta = experiment.flight.lab.probe.log_factor(old["source"]).tolist()
    radii, work = empirical()
    comparisons = audit(lab, delta, work)
    upper = max(row["absolute_outer_to_inner"][1] for row in lab.values())
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "OCCUPIED_WORK.md",
                HERE / "thermal_work.json",
                LAB / "occupied_equilibrium.json",
                LAB / "pair_contact.json",
            )
        },
        galaxy="NGC2841",
        radii_kpc=radii,
        required_work_km2_s2=work,
        source_displacement=delta,
        observational_work_ratio=prior["required_work_ratio"],
        positive_sum_ratio_upper=upper,
        minimum_ratio_discrepancy=prior["required_work_ratio"][0] / upper,
        unresolved_tracer_test=comparisons,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    opt = cli.parse_args()
    output = run()
    path = HERE / "occupied_work.json"
    if opt.write:
        path.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    else:
        experiment.flight.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_bytes()), output
        )
    out = output["unresolved_tracer_test"]
    print(
        out["excluded_every_mixture"],
        len(out["disjoint_dyadic_cells"]),
        out["comparisons"],
    )
