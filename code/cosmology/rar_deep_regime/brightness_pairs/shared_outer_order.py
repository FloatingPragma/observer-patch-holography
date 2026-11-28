"""Test a single monotone force enhancement across the catalogue, including common nuisance parameters."""

import argparse
from bisect import bisect_right
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, vstack

from . import outer_requirement as sample

HERE = sample.HERE


def split(name):
    return "held" if sha256(name.encode()).digest()[0] % 2 else "calibration"


def constraints(groups, data, mode):
    table = sample.measure.DATA
    metas = {z[:11].strip(): z for z in (table / "table1.dat").read_text().splitlines()}
    rings = (table / "table2.dat").read_text().splitlines()
    names = [n for n in data if split(n) in groups]
    entries = [(n, r) for n in names for r in data[n]["rings"]]
    upper = []
    rows = []

    def var(top):
        upper.append(F(top))
        return len(upper) - 1

    def add(terms, rhs=0):
        scale = max(F(1), abs(rhs), *(abs(w) for w in terms.values()))
        rows.append(({i: F(v) / scale for i, v in terms.items()}, F(rhs) / scale))

    common = {}
    for name in names:
        hl, hh = sample.slope_box(metas[name])
        lo = max(F(0), F(sample.lower(hl)))
        hi = F(sample.upper(hh))
        h, d, b = var(hi), var(hi), var(F("1.4") * hi)
        add({h: -1}, -lo)
        for index, L, H in ((d, F("0.1"), 1), (b, F("0.1"), F("1.4"))):
            add({index: -1, h: L})
            add({index: 1, h: -H})
        common[name] = (h, d, b)
    ordering = []
    labels = {}
    for name, r in entries:
        t = rings[r["line"] - 1]
        v, e = F(t[26:32]), F(t[33:38])
        parts = [F(t[a:b]) for a, b in ((39, 45), (46, 52), (53, 59))]
        force = dict(zip(common[name], [p * abs(p) for p in parts], strict=True))
        y = var(1)
        # (v-e)^2*y <= h*B <= (v+e)^2*y with e equal to two quoted errors.
        add({**{i: -z for i, z in force.items()}, y: max(0, v - 2 * e) ** 2})
        add({**force, y: -((v + 2 * e) ** 2)})
        brightness = F(t[60:67]) + F(t[68:76])
        if brightness > 0:
            ordering.append((brightness, y))
        labels[y] = (name, r["line"])
    ordering.sort()
    if mode == "strict":
        for (old, i), (new, j) in zip(ordering, ordering[1:]):
            add({i: 1, j: -1})
            if old == new:
                add({i: -1, j: 1})
    elif mode == "factor_two":
        levels = [v for v, _ in ordering]
        maxima = []
        for j, (v, i) in enumerate(ordering):
            m = var(1)
            add({i: 1, m: -1})
            if j:
                add({maxima[-1]: 1, m: -1})
            last = bisect_right(levels, v / 2) - 1
            if last >= 0:
                add({maxima[last]: 1, i: -1})
            maxima.append(m)
    else:
        raise ValueError(mode)
    for j, u in enumerate(upper):
        add({j: 1}, u)
    return rows, upper, labels


def program(groups, data, mode):
    rows, U, labels = constraints(groups, data, mode)
    ij = [(i, j, float(v)) for i, (row, rhs) in enumerate(rows) for j, v in row.items()]
    i, j, v = zip(*ij, strict=True)
    a = coo_matrix((v, (i, j)), shape=(len(rows), len(U))).tocsr()
    rhs = np.array([float(z) for _, z in rows])
    fit = linprog(np.zeros(len(U)), A_ub=a, b_ub=rhs, bounds=(0, None), method="highs")
    summary = {
        "variables": len(U),
        "constraints": len(rows),
        "rings": len(labels),
        "solver_status": int(fit.status),
    }
    if fit.success:
        summary["max_scaled_primal_residual"] = float(
            max(np.max(a @ fit.x - rhs), np.max(-fit.x))
        )
        return summary
    if fit.status != 2:
        raise ArithmeticError(fit.message)
    # A subprobability combination of rows gives an approximate separating
    # hyperplane. Check the separation using exact rational original rows
    # and the allowed nonnegative variable box, even with a dual residual.
    dual = linprog(
        rhs,
        A_ub=vstack([-a.T, np.ones((1, len(rows)))]).tocsr(),
        b_ub=np.r_[np.zeros(len(U)), 1],
        bounds=(0, None),
        method="highs",
    )
    if not dual.success:
        raise ArithmeticError(dual.message)
    residual = [F(0)] * len(U)
    ceiling = F(0)
    used = []
    for ix, weight in enumerate(dual.x):
        if weight <= 0:
            continue
        y = F(float(weight))
        co, b = rows[ix]
        used.append([ix, y.numerator, y.denominator])
        ceiling += y * b
        for col, n in co.items():
            residual[col] += y * n
    floor = sum(min(z, 0) * u for z, u in zip(residual, U, strict=True))
    assert floor > ceiling
    ncommon = 3 * sum(split(n) in groups for n in data)
    support = sorted(
        {
            labels[j]
            for index, *_ in used
            for j in rows[index][0]
            if j in labels and min(rows[index][0]) < ncommon
        }
    )
    summary.update(
        support_catalogue_rings=support,
        certificate_rows=used,
        exact_gap=[(floor - ceiling).numerator, (floor - ceiling).denominator],
    )
    return summary


def run():
    sample.iv.dps = 36
    base = HERE / "outer_requirement.json"
    d = json.loads(base.read_bytes())["galaxies"]
    sets = ["calibration", "held"]
    rows = [
        (name, r) for name, a in d.items() for r in a["rings"] if r["brightness"] > 0
    ]
    witnesses = {
        k: {"count": 0, "worst": None} for k in ("calibration", "held", "transfer")
    }
    for left, right in combinations(rows, 2):
        if left[0] == right[0]:
            continue
        bright, faint = sorted(
            (left, right), key=lambda x: x[1]["brightness"], reverse=True
        )
        if bright[1]["brightness"] < 2 * faint[1]["brightness"]:
            continue
        lo = bright[1]["total_over_baryonic"][0]
        hi = faint[1]["total_over_baryonic"][1]
        if hi is None or lo <= hi:
            continue
        a, b = (split(n) for n, _ in (bright, faint))
        out = witnesses[a if a == b else "transfer"]
        out["count"] += 1
        if out["worst"] is None or lo / hi > out["worst"]["separation_factor"]:
            out["worst"] = {
                "bright": [bright[0], bright[1]],
                "faint": [faint[0], faint[1]],
                "separation_factor": lo / hi,
            }
    groups_out = {}
    for groups in ([sets[0]], [sets[1]], sets):
        modes = {mode: program(groups, d, mode) for mode in ("factor_two", "strict")}
        for mode, entry in modes.items():
            if "support_catalogue_rings" in entry:
                support = set(entry["support_catalogue_rings"])
                reduced = {
                    n: {
                        "rings": [r for r in d[n]["rings"] if (n, r["line"]) in support]
                    }
                    for n in d
                    if any(g == n for g, _ in support)
                }
                entry["isolated_witness"] = program(groups, reduced, mode)
                entry["support_catalogue_rings"] = list(
                    map(list, entry["support_catalogue_rings"])
                )
                entry["isolated_witness"]["support_catalogue_rings"] = list(
                    map(
                        list,
                        entry["isolated_witness"].get("support_catalogue_rings", []),
                    )
                )
        groups_out["+".join(groups)] = modes
    return {
        "witnesses": witnesses,
        "groups": groups_out,
        "inputs": {
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                base,
                HERE / "SHARED_OUTER_ORDER.md",
                sample.measure.DATA / "table1.dat",
                sample.measure.DATA / "table2.dat",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    d = run()
    path = HERE / "shared_outer_order.json"
    if cli.parse_args().write:
        path.write_text(json.dumps(d, indent=1, allow_nan=False) + "\n")
    else:
        assert json.loads(path.read_bytes()) == d
    print(d["witnesses"])
    print(
        {
            k: {m: v["solver_status"] for m, v in a.items()}
            for k, a in d["groups"].items()
        }
    )
