"""Test nonnegative spherical mass at every eligible SPARC orbit in one galaxy."""

import argparse
import json
from hashlib import sha256
from math import radians, sin

import numpy as np
from scipy.optimize import linprog

from . import mixture
from cosmology.rar_deep_regime.individual_trail_budget import measure as o

TOL = 1e-8


def shape(rows, meta, box, cap):
    k, ml_min, ml_max, gas_min, gas_max = box
    inc = meta["inclination_deg"]
    err_i = k * meta["e_inclination_deg"]
    err_d = k * meta["e_dist_mpc"] / meta["dist_mpc"]

    def angular(angle):
        return (sin(radians(angle)) / sin(radians(inc))) ** 2

    tlo = max(0.0, 1 - err_d) * angular(max(0.0, inc - err_i))
    thi = (1 + err_d) * angular(min(90.0, inc + err_i))
    luminosity = meta["L36_gsun"]
    light_error = k * meta["e_L36_gsun"] / luminosity
    slo = ml_min * max(0.0, 1 - light_error)
    shi = ml_max * (1 + light_error)
    gas = o.sparc.GAS_HELIUM * meta["MHI_gsun"]
    radius = np.array([r["rad_kpc"] for r in rows])
    assert np.all(np.diff(radius) > 0)
    v = np.array([r["vobs_kms"] for r in rows])
    sigma = np.array([r["e_vobs_kms"] for r in rows])
    lo = (np.maximum(0.0, v - k * sigma) / v) ** 2
    hi = ((v + k * sigma) / v) ** 2
    disk = np.array([o.signed_square(r["vdisk_kms"]) for r in rows]) / v**2
    gas_field = np.array([o.signed_square(r["vgas_kms"]) for r in rows]) / v**2
    kernel = 1e9 * o.G / (radius * v**2)

    # Unknowns are common (t,x,y), uniform widening in units of each v²,
    # and nonnegative enclosed masses in scaled units of 10^9 solar masses.
    n = len(rows)
    a = np.zeros((4 + 2 * n + n - 1 + int(cap is not None), 4 + n))
    b = np.zeros(len(a))
    a[:4, :3] = [[slo, -1, 0], [-shi, 1, 0], [gas_min, 0, -1], [-gas_max, 0, 1]]
    for j in range(n):
        for sign, ix, bound in ((1, 4 + j, hi[j]), (-1, 4 + n + j, -lo[j])):
            a[ix, 1:4] = (sign * disk[j], sign * gas_field[j], -1)
            a[ix, 4 + j] = sign * kernel[j]
            b[ix] = bound
        if j + 1 < n:
            a[4 + 2 * n + j, 4 + j : 6 + j] = (1, -1)
    if cap is not None:
        a[-1, [1, 2, -1]] = (-cap * luminosity, -cap * gas, 1)
    objective = np.zeros(4 + n)
    objective[3] = 1
    fit = linprog(
        objective,
        A_ub=a,
        b_ub=b,
        bounds=[(tlo, thi)] + [(0, None)] * (3 + n),
        method="highs",
        options={
            "primal_feasibility_tolerance": 1e-9,
            "dual_feasibility_tolerance": 1e-9,
        },
    )
    if not fit.success:
        raise RuntimeError(fit.message)
    t, x, y, delta = fit.x[:4]
    if t <= 0:
        raise ValueError("Degenerate zero distance/inclination in the fit")
    # The least admissible nondecreasing profile is the running maximum
    # of the lower bounds. This independently checks the fitted witness.
    floor = np.maximum.accumulate(
        np.maximum(0, (lo - delta - x * disk - y * gas_field) / kernel)
    )
    model = x * disk + y * gas_field + kernel * floor
    assert np.all(model <= hi + delta + TOL) and np.all(model >= lo - delta - TOL)
    if cap is not None:
        assert floor[-1] <= cap * (x * luminosity + y * gas) + TOL
    # A dual lower bound certifies any positive required widening. Check
    # all multipliers against stationarity, feasibility and equality of values.
    lam = fit.ineqlin.marginals
    lower, upper = fit.lower.marginals, fit.upper.marginals
    np.testing.assert_allclose(objective, a.T @ lam + lower + upper, atol=TOL, rtol=0)
    assert max(lam) <= TOL and min(lower) >= -TOL and max(upper) <= TOL
    assert np.max(np.abs(upper[1:])) < TOL
    dual_value = float(b @ lam + tlo * lower[0] + thi * upper[0])
    assert abs(dual_value - delta) < TOL
    return dict(
        widening=delta,
        dual_widening=dual_value,
        feasible=bool(delta < TOL),
        t=t,
        scaled_stellar=x,
        scaled_gas=y,
        scaled_enclosed_mass=floor.tolist(),
    )


def run():
    mixture_path = mixture.bao.HERE / "mixture.json"
    context = json.loads(mixture_path.read_bytes())
    paid, extra = mixture.inventory()
    caps = {
        "cosmic_95": context["highest_galaxy_cap_at_95"],
        "maximum_native": float(max(paid) + extra),
        "unlimited": None,
    }
    meta, groups = o.catalogue()
    kept, _ = o.select(meta, groups)
    out = {}
    for box, settings in o.BOXES.items():
        cases = {}
        for name, rows in kept.items():
            curves = {
                tag: shape(rows, meta[name], settings, cap) for tag, cap in caps.items()
            }
            outer = shape(rows[-1:], meta[name], settings, caps["cosmic_95"])
            cases[name] = dict(
                radii_kpc=[r["rad_kpc"] for r in rows],
                curves=curves,
                outer_only_feasible=outer["feasible"],
            )
        missed = {
            tag: sorted(k for k, c in cases.items() if not c["curves"][tag]["feasible"])
            for tag in caps
        }
        out[box] = dict(
            missed=missed,
            outer_only_missed=sorted(
                k for k, c in cases.items() if not c["outer_only_feasible"]
            ),
            galaxies=cases,
        )
    # Diagnose the failed two-error case on the same data. The cap is
    # monotone, so all systems feasible at the smaller cap keep their fits.
    followup = {}
    for name in out["two_errors"]["missed"]["cosmic_95"]:
        rows = kept[name]
        box = o.BOXES["two_errors"]
        a, b = caps["cosmic_95"], caps["maximum_native"]
        if not out["two_errors"]["galaxies"][name]["curves"]["maximum_native"][
            "feasible"
        ]:
            continue
        while b - a > 0.00001:
            mid = (a + b) / 2
            if shape(rows, meta[name], box, mid)["feasible"]:
                b = mid
            else:
                a = mid
        bg = mixture.profile_fit(mixture.bao.DistanceData(), b - float(extra))
        followup[name] = dict(
            minimum_cap_bracket=[a, b],
            distance_fit_at_feasible_cap=bg,
            delta_chi2=bg["chi2"] - context["best"]["chi2"],
            leaving_one_ring_out={
                str(r["rad_kpc"]): shape(
                    rows[:j] + rows[j + 1 :], meta[name], box, caps["cosmic_95"]
                )["widening"]
                for j, r in enumerate(rows)
                if len(rows) > 1
            },
        )
    return dict(
        protocol_sha256=sha256(
            (mixture.bao.HERE / "RINGS_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        mixture_sha256=sha256(mixture_path.read_bytes()).hexdigest(),
        caps=caps,
        source_sparc_sha256=context["source_sparc_sha256"],
        eligible_galaxies=len(kept),
        eligible_rings=sum(map(len, kept.values())),
        boxes=out,
        two_error_followup=followup,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    path = mixture.bao.HERE / "rings.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), result)
    for box, detail in result["boxes"].items():
        print(
            box,
            {k: len(v) for k, v in detail["missed"].items()},
            len(detail["outer_only_missed"]),
        )
