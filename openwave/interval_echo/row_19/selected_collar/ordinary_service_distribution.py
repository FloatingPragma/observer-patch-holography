"""Test every nonnegative allocation of the ordinary source over twelve service histories."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json

import numpy as np
from scipy.optimize import linprog

from . import density_material as mat

d, lab, reg = mat.density, mat.lab, mat.reg
HERE = mat.HERE
SEED = 28112112
POINTERS = ("direct", "source")


def separate(columns, target, names):
    # For x >= 0, necessary interval inequalities are L x <= target_hi and
    # -U x <= -target_lo. Keep the original binary endpoints as exact rationals.
    a, b, labels = [], [], []
    for name in names:
        low, high = target[name]
        for side, multiplier, endpoint, rhs in (
            ("lower", 1, 0, high),
            ("upper", -1, 1, -low),
        ):
            a.append([multiplier * F(col[name][endpoint]) for col in columns.values()])
            b.append(F(rhs))
            labels.append([name, side])
    units = [
        max(abs(v) for v in [*row, bound]) or F(1)
        for row, bound in zip(a, b, strict=True)
    ]
    A = np.array(
        [[float(v / unit) for v in row] for row, unit in zip(a, units, strict=True)]
    )
    B = np.array([float(v / unit) for v, unit in zip(b, units, strict=True)])
    n = len(labels)
    # Seek a strict dual margin on every unbounded ordinary weight. The
    # floating solution is only a proposal for the rational check below.
    program = np.column_stack(
        (np.vstack((-A.T, B, np.ones(n))), np.r_[np.ones(A.shape[1] + 1), 0])
    )
    dual = linprog(
        np.r_[np.zeros(n), -1],
        A_ub=program,
        b_ub=np.r_[np.zeros(A.shape[1] + 1), 1],
        bounds=(0, None),
        method="highs",
    )
    if not dual.success:
        raise RuntimeError(dual.message)
    weights = [
        F(float(max(0, y))) / unit for y, unit in zip(dual.x[:-1], units, strict=True)
    ]
    totals = [
        sum((y * v for y, v in zip(weights, col, strict=True)), F(0))
        for col in zip(*a, strict=True)
    ]
    rhs = sum((w * v for w, v in zip(weights, b, strict=True)), F(0))
    if min(totals) >= 0 and rhs < 0:
        return {
            "status": "excluded",
            "coefficients": [
                [*z, w.numerator, w.denominator]
                for z, w in zip(labels, weights, strict=True)
                if w
            ],
            "minimum_ordinary_coefficient": float(min(totals)),
            "negative_target_bound": float(rhs),
        }
    fit = linprog(
        np.zeros(A.shape[1]), A_ub=A, b_ub=B, bounds=(0, None), method="highs"
    )
    if fit.success:
        # Only report compatibility if a nonnegative rational allocation
        # satisfies every outward interval use, including the original scale.
        x = [F(float(max(0, v))) for v in fit.x]
        # A solver optimum lies on the boundary. Seek an interior point.
        if any(
            sum(v * w for v, w in zip(row, x, strict=True)) > bound
            for row, bound in zip(a, b, strict=True)
        ):
            inner = linprog(
                np.r_[np.zeros(A.shape[1]), -1],
                A_ub=np.column_stack((A, np.ones(len(B)))),
                b_ub=B,
                bounds=[(0, None)] * A.shape[1] + [(0, 1)],
                method="highs",
            )
            if inner.success:
                x = [F(float(max(0, v))) for v in inner.x[:-1]]
        if all(
            sum(v * w for v, w in zip(row, x, strict=True)) <= bound
            for row, bound in zip(a, b, strict=True)
        ):
            return {
                "status": "compatible_relaxation",
                "allocation": [
                    [k, w.numerator, w.denominator]
                    for k, w in zip(columns, x, strict=True)
                    if w
                ],
            }
    return {"status": "unresolved"}


def run():
    mat.iv.dps = 35
    rng = np.random.default_rng(SEED)
    _, inputs, *_ = mat.c.trials()
    ports, bins = lab.six.survey.return_bins()
    effects, _, _ = mat.beam.flight(mat.METERS, include_edges=False)
    prices = {f"{meter}/{sign}": np.asarray(E) for (sign, meter), E in effects.items()}
    names = list(prices)
    # The timeout constant cancels between normalized ensembles.
    clocks = np.einsum("rtbij,t->rbij", bins, reg.SCORE[:-1] - 3)
    clock_names = []
    for route, price in zip(lab.six.survey.ROUTES, clocks, strict=True):
        key = f"clock/{route[0]}/{route[1]}"
        prices[key] = price
        clock_names.append(key)
    panels, observed = {t: {} for t in d.WINDOWS}, {}
    for sender in ports:
        for sign in (-1, 1):
            label = f"{sender}/{sign}"
            runs = {}
            for arm in ("reference", "ordinary"):
                _, _, runs[arm] = d.history(inputs[arm], (sender,), rng, signs=(sign,))
            for tau in d.WINDOWS:
                tape = {
                    arm: lab.pairs(menu[tau]["dense"].sum(0), rng, False)
                    for arm, menu in runs.items()
                }
                observed[f"{tau}/{label}"] = tape
                read = {
                    k: reg.src.observed(raw, lab.SHOTS, False)
                    for k, raw in tape.items()
                }
                panels[tau][label] = {
                    name: lab.six.bounds(
                        d.DOSE
                        * mat.forecast(price, read["ordinary"], read["reference"])
                    )
                    for name, price in prices.items()
                }
    prior = {
        name: json.loads((HERE / (name + ".json")).read_bytes())
        for name in (
            "density_material_active",
            "density_material_impulse",
            "regional_density",
        )
    }
    windows = {}
    for tau in d.WINDOWS:
        key = str(tau)
        source = prior[
            "density_material_impulse" if tau == 2 else "density_material_active"
        ]["windows"][key]
        meters = source["exposures"]["76" if tau == 2 else "72"]
        target = {k: {} for k in POINTERS}
        ordinary = {}
        for meter, row in meters.items():
            for sign, p in row["probes"].items():
                name = f"{meter}/{sign}"
                target["direct"][name] = p["contrasts"]["sparse"]
                target["source"][name] = p["source_forecasts"]["sparse"]
                ordinary[name] = p["contrasts"]["ordinary"]
        times = prior["regional_density"]["windows"][key]["radar_counts"]
        means, _ = reg.timing_slopes(times, lab.SHOTS, (("sparse", "dense", 1),))
        for j, name in enumerate(clock_names):
            for t in target.values():
                t[name] = lab.six.bounds(means["sparse"][j] - means["dense"][j])
            ordinary[name] = lab.six.bounds(means["ordinary"][j] - means["dense"][j])
        uniform = {
            name: lab.six.bounds(
                sum(mat.iv.mpf(z[name]) for z in panels[tau].values())
                / len(panels[tau])
            )
            for name in prices
        }
        failed = [
            name
            for name in prices
            if max(ordinary[name][0], uniform[name][0])
            > min(ordinary[name][1], uniform[name][1])
        ]
        groups = {
            "impulse": [n for n in names if n.startswith("shift/")],
            "marked": [n for n in names if n.startswith("marked_force/")],
            "material": names,
            "clocks": clock_names,
            "joint": [*names, *clock_names],
        }
        windows[key] = {
            "ordinary_columns": panels[tau],
            "uniform_ordinary": uniform,
            "held_ordinary": ordinary,
            "ordinary_controls_disagree": failed,
            "target": target,
            "groups": groups,
            "comparisons": {
                p: {
                    group: separate(panels[tau], t, rows)
                    for group, rows in groups.items()
                }
                for p, t in target.items()
            },
        }
    return {
        "seed": SEED,
        "active_head_starts": lab.SHOTS,
        "ordinary_head_counts": observed,
        "windows": windows,
        "inputs": {
            n: sha256((HERE / n).read_bytes()).hexdigest()
            for n in (
                "ORDINARY_SERVICE_DISTRIBUTION.md",
                *(k + ".json" for k in prior),
                "finite_pairing_dose.json",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    r = run()
    path = HERE / "ordinary_service_distribution.json"
    if args.write:
        path.write_text(json.dumps(r, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), r)
    for k, w in r["windows"].items():
        print(
            k,
            w["ordinary_controls_disagree"],
            {
                p: {g: z["status"] for g, z in rows.items()}
                for p, rows in w["comparisons"].items()
            },
        )
