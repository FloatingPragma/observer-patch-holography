"""Test the spherical midpoint source with a free normalization and edge on each SPARC galaxy."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import ceil

import numpy as np
from scipy.optimize import linprog

from . import outer_requirement as sample

HERE = sample.HERE
PROTOCOL = (
    HERE.parents[3] / "openwave/interval_echo/row_19/selected_collar/MIDPOINT_BULK.md"
)


def rings_in(raw, entries, extra_speed_error=F(0)):
    rings = []
    for entry in entries:
        line = entry["line"]
        t = raw[line - 1]
        v, e = F(t[26:32]), F(t[33:38])
        signed = [
            v * abs(v) for v in (F(t[a:b]) for a, b in ((39, 45), (46, 52), (53, 59)))
        ]
        assert (
            min(signed[0] + x * signed[1] + y * signed[2] for x, y in sample.RECT) > 0
        )
        rings.append(
            (
                F(t[19:25]),
                max(0, v - 2 * e - extra_speed_error) ** 2,
                (v + 2 * e + extra_speed_error) ** 2,
                signed,
                line,
            )
        )
    return sorted(rings)


def problem(rings, lo, hi, cut, compact_cap=None):
    # Speeds at catalogue distances and inclinations. A and D absorb the
    # common geometric conversions. Their ratio gives the squared edge in
    # catalogue kpc. Every positive physical radius and gain is allowed.
    radii = sorted({r[0] for r in rings})
    rmin, rmax = radii[0], radii[-1]
    fast = max(r[2] for r in rings)
    mixture = cut is None
    source_caps = [fast / rmin] * len(radii) if mixture else [fast / rmin, fast * rmax]
    bary_caps = [hi, hi, F(7, 5) * hi] if compact_cap is None else [compact_cap]
    cap = [*bary_caps, *source_caps]
    rows = []

    def add(row, rhs):
        data = [F(val) * u for val, u in zip(row, cap, strict=True)]
        scale = max(F(1), abs(rhs), *(abs(v) for v in data))
        rows.append(([v / scale for v in data], F(rhs) / scale))

    zeros = [0] * len(source_caps)
    if compact_cap is None:
        add([-1, 0, 0, *zeros], -lo)
        for j, mx in ((1, F(1)), (2, F(7, 5))):
            row = [F(1, 10), 0, 0, *zeros]
            row[j] = -1
            add(row, 0)
            row = [-mx, 0, 0, *zeros]
            row[j] = 1
            add(row, 0)
    for r, slow, fast, parts, _ in rings:
        if mixture:
            row = [*parts, *(min(r, L * L / r) for L in radii)]
        else:
            inner = cut > 0 and r <= radii[cut - 1]
            row = [*parts, r if inner else 0, 0 if inner else 1 / r]
        add(row, fast)
        add([-a for a in row], -slow)
    if not mixture and cut:
        add([*([0] * len(bary_caps)), radii[cut - 1] ** 2, -1], 0)
    if not mixture and cut < len(radii):
        add([*([0] * len(bary_caps)), -(radii[cut] ** 2), 1], 0)
    return rows, cap


def solve(rows):
    a = np.array([z for z, _ in rows], float)
    b = np.array([z for _, z in rows], float)
    # Maximize the least normalized slack in a finite physical box.
    return linprog(
        [0] * a.shape[1] + [-1],
        A_ub=np.c_[a, np.ones(len(a))],
        b_ub=b,
        bounds=[(0, 1)] * a.shape[1] + [(None, None)],
        method="highs",
    )


def certify(rings, bounds, cut, compact_cap=None):
    lo, hi = bounds
    outer, upper = problem(
        rings, max(F(0), F(sample.lower(lo))), F(sample.upper(hi)), cut, compact_cap
    )
    trial = solve(outer)
    assert trial.success
    if trial.x[-1] < -1e-10:
        y = [F(float(-v)).limit_denominator(10**10) for v in trial.ineqlin.marginals]
        assert min(y) >= 0
        left = [
            sum(w * row[i] for w, (row, _) in zip(y, outer, strict=True))
            for i in range(len(upper))
        ]
        right = sum(w * b for w, (_, b) in zip(y, outer, strict=True))
        gap = sum(min(v, 0) for v in left) - right
        assert gap > 0
        return {
            "possible": False,
            "multipliers": [[i, str(w)] for i, w in enumerate(y) if w],
            "exact_gap": str(gap),
        }
    inner, cap = problem(
        rings, max(F(0), F(sample.upper(lo))), F(sample.lower(hi)), cut, compact_cap
    )
    res = solve(inner)
    assert res.success
    x = [F(float(v)).limit_denominator(10**11) for v in res.x[:-1]]
    assert all(0 <= v <= 1 for v in x)
    for row, b in inner:
        assert sum(w * v for w, v in zip(row, x, strict=True)) <= b, (cut, res.x[-1])
    return {
        "possible": True,
        "parameters": [str(a * b) for a, b in zip(cap, x, strict=True)],
    }


LENS_INPUT = HERE.parent / "compensation/trail_profile/inputs.json"


def lensing():
    # Published velocities are already spherically deprojected ensemble
    # measurements. A common profile across a stack is an extra hypothesis.
    src = json.loads(LENS_INPUT.read_bytes())
    out = {}
    for j, logmass in enumerate(("10.5", "10.8", "11.1", "11.5")):
        # Round the doubled maximum bin mass upward to whole solar masses.
        # Convert the catalogue G from Mpc to kpc in the speed equation.
        cap = F("0.000004300917270") * ceil(
            sample.upper(2 * 10 ** sample.iv.mpf(logmass))
        )
        output = {}
        for factor in (1, 2):
            rings = []
            for k, entry in enumerate(src["rows"]):
                r = F(str(entry["radius_kpc"]))
                v, stat, sys = [F(str(a)) for a in entry["velocity_stat_sys"][j]]
                err = factor * (stat + sys)
                rings.append((r, max(0, v - err) ** 2, (v + err) ** 2, [1 / r], k))
            bounds = [sample.iv.mpf(1)] * 2
            cuts = [certify(rings, bounds, k, cap) for k in range(len(rings) + 1)]
            output[str(factor)] = {
                "compact_GM_cap": str(cap),
                "cuts": cuts,
                "all_inside": cuts[-1]["possible"],
                "any_edge": any(z["possible"] for z in cuts),
                "mixture": certify(rings, bounds, None, cap),
            }
        out[logmass] = output
    return out


COLD_INPUT = HERE.parents[1] / "cold_dilution/results.json"
WORD_INPUT = PROTOCOL.parent / "three_regions.json"


def dilution():
    # In one fixed comoving volume the total gravitating energy of the
    # expanding screens scales as their summed word price divided by L.
    measured = json.loads(WORD_INPUT.read_bytes())["signed_stress_sweep"][-1][
        "difference"
    ]
    empirical = json.loads(COLD_INPUT.read_bytes())["constant_exponent_test"]
    n_lo, n_hi = map(lambda v: F(str(v)), empirical["dilution_exponent_99_7pct"])
    required = [4 - n_hi, 4 - n_lo]
    gain = sample.iv.mpf(measured) * 15 / (32 * sample.iv.pi)
    return {
        "measured_pair_price_nats": measured,
        "mass_times_radius_per_pair_in_natural_units": [
            float(sample.lower(gain)),
            float(sample.upper(gain)),
        ],
        "observational_dilution_power": [str(n_lo), str(n_hi)],
        "required_power_of_total_pair_price_in_comoving_volume": list(
            map(str, required)
        ),
        "required_pair_price_factor_per_scale_doubling": [
            16 * v for v in empirical["density_ratio_for_doubling_99_7pct"]
        ],
        "fixed_price_expanding_screen_density_ratio": "1/16",
        "price_proportional_to_radius_density_ratio": "1/8",
        "bound_screens_with_conserved_mass_density_ratio": "1/8",
    }


def run():
    sample.iv.dps = 36
    src = HERE / "outer_requirement.json"
    selected = json.loads(src.read_bytes())["galaxies"]
    table1 = sample.measure.DATA / "table1.dat"
    table2 = sample.measure.DATA / "table2.dat"
    meta = {t[:11].strip(): t for t in table1.read_text().splitlines()}
    raw = table2.read_text().splitlines()
    results = {}
    for name, record in selected.items():
        rings = rings_in(raw, record["rings"])
        radii = sorted({r[0] for r in rings})
        assert not record["nonpositive_template_possible_lines"]
        bounds = sample.slope_box(meta[name])
        cases = [certify(rings, bounds, k) for k in range(len(radii) + 1)]
        results[name] = {
            "lines": [r[-1] for r in rings],
            "radii_kpc": list(map(str, radii)),
            "by_number_of_inside_radii": cases,
            "any_edge": any(v["possible"] for v in cases),
            "nonnegative_mixture": certify(rings, bounds, None),
            "all_inside": cases[-1]["possible"],
        }
        mixture = results[name]["nonnegative_mixture"]
        if not mixture["possible"]:
            support = sorted(
                {(j - 5) // 2 for j, _ in mixture["multipliers"] if j >= 5}
            )
            limited = [rings[i] for i in support]
            small = certify(limited, bounds, None)
            assert not small["possible"]
            results[name]["mixture_witness"] = {
                "lines": [row[-1] for row in limited],
                "proof": small,
            }
            left, right = F(0), F(1)

            def extra(dv):
                return certify(rings_in(raw, record["rings"], dv), bounds, None)

            while not extra(right)["possible"]:
                right *= 2
            while right - left > F(1, 64):
                mid = (left + right) / 2
                if extra(mid)["possible"]:
                    right = mid
                else:
                    left = mid
            results[name]["additional_speed_halfwidth_kms"] = {
                "excluded_at": str(left),
                "feasible_at": str(right),
                "low_certificate": extra(left),
                "high_certificate": extra(right),
            }
    summary = {
        "galaxies": len(results),
        "rings": sum(len(v["lines"]) for v in results.values()),
        "interior_only_excluded": [
            g for g, v in results.items() if not v["all_inside"]
        ],
        "free_edge_excluded": [g for g, v in results.items() if not v["any_edge"]],
        "arbitrary_centered_mixture_excluded": [
            g for g, v in results.items() if not v["nonnegative_mixture"]["possible"]
        ],
    }
    return {
        "summary": summary,
        "galaxies": results,
        "deprojected_lens_stacks": lensing(),
        "cosmological_dilution": dilution(),
        "inputs": {
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (src, table1, table2, PROTOCOL, LENS_INPUT, COLD_INPUT, WORD_INPUT)
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    data = run()
    path = HERE / "midpoint_halo.json"
    if args.write:
        path.write_text(json.dumps(data, indent=1, allow_nan=False) + "\n")
    else:
        assert json.loads(path.read_bytes()) == data
    print(data["summary"])
