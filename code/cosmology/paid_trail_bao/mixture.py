"""Profile a universal native stock against DESI+BBN and SPARC outer speeds."""

import argparse
from fractions import Fraction
from functools import cache
import hashlib
import json

import camb
import numpy as np
from scipy.optimize import brentq, least_squares, minimize_scalar

from . import measure as bao
from cosmology.rar_deep_regime.individual_trail_budget import measure as orbits

BBN, SIGMA = 0.02218, 0.00055
DELTA = 3.84146


def inventory():
    caps = orbits.native_caps()
    boards = [bao.prepare(b) for b in json.loads(bao.RECEIPT.read_bytes())["boards"]]
    paid = [Fraction(b["ratio_exact"]) for b in boards]
    offsets = {
        Fraction(caps[str(b["live_seams"])]["cap_exact"]) - q
        for b, q in zip(boards, paid, strict=True)
    }
    if len(offsets) != 1:
        raise ValueError("No common offset for these trail lengths")
    return paid, offsets.pop()


def profile_fit(data, q):
    def residuals(x):
        h, baryon100 = x
        baryon = baryon100 / 100
        _, v = bao.physical_background(data, baryon, q, h)
        return np.r_[bao.whiten(data.full, data.obs - v), (baryon - BBN) / SIGMA]

    fit = least_squares(
        residuals,
        [68.5, 100 * BBN],
        bounds=([35, 1.2], [115, 3.4]),
        jac="3-point",
        diff_step=5e-5,
        xtol=2e-8,
        ftol=1e-10,
    )
    if not fit.success or any(fit.active_mask):
        raise RuntimeError(fit.message)
    return dict(
        q=q,
        chi2=float(fit.fun @ fit.fun),
        H0=fit.x[0],
        omega_b=fit.x[1] / 100,
        BBN_chi2=fit.fun[-1] ** 2,
    )


def run():
    paid, offset = inventory()
    low, high = float(min(paid)), float(max(paid))
    data = bao.DistanceData()

    @cache
    def profile(q):
        return profile_fit(data, q)

    def score(q):
        return profile(float(q))["chi2"]

    best = minimize_scalar(
        score, method="bounded", bounds=(low, high), options={"xatol": 2e-5}
    )
    if not best.success:
        raise RuntimeError(best.message)
    q0 = float(best.x)
    grid = sorted({*np.linspace(low, high, 17), q0})
    scores = [score(q) for q in grid]
    turn = grid.index(q0)
    assert all(np.diff(scores[: turn + 1]) < 0) and all(np.diff(scores[turn:]) > 0)
    lo, hi = [
        brentq(lambda q: score(q) - best.fun - DELTA, a, b, xtol=1e-5)
        for a, b in ((low, q0), (q0, high))
    ]
    galaxies = orbits.run()
    by_box = {}
    for name in orbits.BOXES:
        required = {
            g["galaxy"]: g["bounds"][name]["q_required_min"]
            for g in galaxies["galaxies"]
        }
        worst = max(required, key=required.get)
        threshold = required[worst] - float(offset)
        constrained = None if threshold > high else profile(max(q0, threshold))
        by_box[name] = dict(
            largest_required_galaxy=worst,
            largest_required=required[worst],
            paid_q_required=threshold,
            missed_at_upper_95=sorted(
                k for k, v in required.items() if v > hi + float(offset)
            ),
            missed_at_best=sum(v > q0 + float(offset) for v in required.values()),
            least_distance_penalty=None
            if constrained is None
            else constrained["chi2"] - best.fun,
            constrained_fit=constrained,
        )
    return dict(
        protocol_sha256=hashlib.sha256(
            (bao.HERE / "MIXTURE_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        source_native_sha256=galaxies["source_native_sha256"],
        source_sparc_sha256=galaxies["source_sparc_sha256"],
        CAMB=camb.__version__,
        desi_sha256={
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(bao.HERE.glob(bao.PREFIX + "*.txt"))
        },
        native_paid=[str(q) for q in paid],
        upper_inventory_minus_paid=str(offset),
        best=profile(q0),
        profile_95=[profile(lo), profile(hi)],
        delta=DELTA,
        highest_galaxy_cap_at_95=hi + float(offset),
        profile_grid=[profile(q) for q in grid],
        galaxy_count=len(galaxies["galaxies"]),
        galaxy_bounds=by_box,
        pressure_sensitivity=galaxies["NGC3741_controls"][
            "removed_pressure_correction"
        ],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    out = json.loads(json.dumps(run(), allow_nan=False))
    path = bao.HERE / "mixture.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), out)
    print(out["best"], out["highest_galaxy_cap_at_95"])
