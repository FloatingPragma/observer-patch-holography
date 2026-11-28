"""Ordinary-only interval calibration and equal-source tests on a shared reference."""

import argparse
from hashlib import sha256
import json

import numpy as np
from scipy.linalg import qr
import sympy as sp

from .six_meter_pointer import HERE, bounds, components, positive


def single_column(bands):
    return components([[v] for v in bands])


def calibrate(source, instruments, alternative, labels):
    names = list(source)
    columns = np.array([[(a + b) / 2 for a, b in source[k]] for k in names]).T
    scales = np.abs(columns).max(axis=1)
    if np.any(scales == 0):
        return {"status": "unresolved zero midpoint coordinate"}
    d = len(columns)
    _, _, pivots = qr(columns / scales[:, None], pivoting=True, mode="economic")
    basis = [names[i] for i in pivots[:d]]
    A, E = components([[source[k][i] for k in basis] for i in range(d)])
    # Numerical proposals have no evidential role until checked as dyadic rationals.
    try:
        proposal = np.linalg.inv(np.array(A, dtype=float))
    except np.linalg.LinAlgError:
        return {"basis": basis, "status": "unresolved singular midpoint"}
    K = sp.Matrix(proposal).applyfunc(sp.Rational)
    residual = sp.eye(d) - K * A
    H = positive(residual) + positive(K) * E
    estimate = np.linalg.solve(np.eye(d) - np.array(H, dtype=float), np.ones(d))
    v = sp.Matrix([sp.Rational(float(z)) for z in estimate])
    if min(v) <= 0:
        return {
            "basis": basis,
            "status": "unresolved inverse on the measured intervals",
        }
    alpha = max(z / w for z, w in zip(H * v, v, strict=True))
    if alpha >= 1:
        return {
            "basis": basis,
            "status": "unresolved inverse on the measured intervals",
        }
    coeff, uncertain = components([[instruments[k][n] for k in basis] for n in labels])
    price = positive(coeff * K)
    remainder = positive(coeff * residual) + uncertain

    def forecast(bands, observed):
        y, ey = single_column(bands)
        x = K * y
        b = positive(K) * ey + H * positive(x)
        rad = max(z / w for z, w in zip(b, v, strict=True)) / (1 - alpha) * v
        worst = positive(x) + rad
        err = price * (ey + E * worst) + remainder * worst
        mid = coeff * x
        readings = {n: bounds(mid[i], err[i]) for i, n in enumerate(labels)}
        differences = {
            n: [
                float(np.nextafter(observed[n][0] - hi, -np.inf)),
                float(np.nextafter(observed[n][1] - lo, np.inf)),
            ]
            for n, (lo, hi) in readings.items()
        }
        dose = sp.Rational(1)
        while 36 * dose * max(1, max(worst)) > 1:
            dose /= 2
        return {
            "coefficients": {k: bounds(x[j], rad[j]) for j, k in enumerate(basis)},
            "common_physical_dose": [int(dose.p), int(dose.q)],
            "prediction": readings,
            "observed_minus_prediction": differences,
            "excluded": [n for n, (lo, hi) in differences.items() if lo > 0 or hi < 0],
        }

    return {
        "status": "verified inverse",
        "basis": basis,
        "weighted_inverse_residual_upper": bounds(alpha, 0)[1],
        "ordinary": {
            k: forecast(source[k], instruments[k]) for k in names if k not in basis
        },
        "overlap": {k: forecast(y, obs) for k, (y, obs) in alternative.items()},
    }


def run():
    p = HERE / "distributed_source.json"
    data = json.loads(p.read_bytes())
    res = data["result"]
    receivers = [str(k) for k in data["ports"]] + ["complement"]
    labels = res["groups"]["shape"] + res["groups"]["material"]
    rows = {}
    for fields in (
        ["energy"],
        ["energy", "connected"],
        ["energy", "connected", "marginal"],
    ):

        def extract(z):
            return [z[n][f] for f in fields for n in receivers]

        controls = {k: extract(z) for k, z in res["ordinary_regional_sources"].items()}
        observed = res["sources"]["connected"]["target"]
        alt = {"measured": (extract(res["held_sources"]), observed)}
        if fields == ["energy"]:
            for k, arm in res["sources"].items():
                alt[k] = ([arm["target"]["site/" + n] for n in receivers], observed)
        rows["+".join(fields)] = calibrate(controls, res["columns"], alt, labels)
    return {
        "labels": labels,
        "calibrations": rows,
        "inputs": {
            f.name: sha256(f.read_bytes()).hexdigest()
            for f in (p, HERE / "REGIONAL_SOURCE_CALIBRATION.md")
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    r = run()
    path = HERE / "regional_source_calibration.json"
    if args.write:
        path.write_text(json.dumps(r, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), r)
    for k, v in r["calibrations"].items():
        print(
            k,
            v["status"],
            sum(bool(row["excluded"]) for row in v.get("ordinary", {}).values()),
            {s: len(row["excluded"]) for s, row in v.get("overlap", {}).items()},
        )
