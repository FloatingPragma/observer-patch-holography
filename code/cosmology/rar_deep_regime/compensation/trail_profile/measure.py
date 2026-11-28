"""Amplitude- and hop-length-independent constraints on two aperture lens readings."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from scipy.stats import norm

from .. import measure as lens
from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    measure as native,
)

HERE = Path(__file__).parent


def acquire(raw):
    if sha256(raw).hexdigest() != lens.MRT_HASH:
        raise ValueError("Different Mistele table")
    rows = []
    for line in raw.decode("ascii").splitlines():
        if line[:24].strip() != "All":
            continue
        radius = float(line[25:32])
        if radius > 300:
            continue
        data = np.array(list(map(float, line[33:].split()))).reshape(4, 3)
        rows.append(dict(radius_kpc=radius, velocity_stat_sys=data.tolist()))
    return dict(sha256=lens.MRT_HASH, rows=rows)


def envelope(r1, r2, board, layout, inner_fraction):
    data = board[layout]
    mid = [Fraction(float(x)) for x in data["midpoints"]]
    # m/d enters a radius when d<=radius/m. Use exact boundaries (including
    # equality), an open-interval representative, and both constant tails.
    r, R = Fraction(str(r1)), Fraction(str(r2))
    boundaries = sorted({x / m for m in mid for x in (r, R) if m})
    probes = [
        boundaries[0] / 2,
        2 * boundaries[-1],
        *boundaries,
        *[(x + y) / 2 for x, y in zip(boundaries, boundaries[1:])],
    ]
    g, E = 8.5, board["initial"]
    # General endpoint degree allows a difference of 1/2 in blank energy.
    deficit = 8.5 * len(mid) + 54 - data["board_energy"]
    upper = {
        k: (-np.inf, None, 0, 0) for k in ("at_release", "all_free_times", "evacuated")
    }
    for d in probes:
        inner = [m * d <= r for m in mid]
        outer = [m * d <= R for m in mid]
        i, j = sum(inner), sum(outer)

        def ratio(a, b):
            return (E + b) / (E * inner_fraction + a) * float(r / R)

        # At formation the last head overlaps the only load and loses 12.75.
        closed = ratio(g * i - 12.75 * inner[-1], g * j - 12.75 * outer[-1])
        # 0 <= D_i <= D_j. The ratio grows with their common subtraction
        # when the numerator >= denominator before the common loss.
        lost = min(deficit, 18 * i)
        bound = ratio(g * i - lost, g * j - lost)
        for label, v in (
            ("at_release", closed),
            ("all_free_times", bound),
            ("evacuated", ratio(g * i, g * j)),
        ):
            if v > upper[label][0]:
                upper[label] = (v, float(d), i, j)
    return {
        k: dict(
            speed_ratio_ceiling=np.sqrt(v),
            extremizing_test_hop_kpc=d,
            inner_heads=i,
            outer_heads=j,
        )
        for k, (v, d, i, j) in upper.items()
    }


def observations(inner, outer, ratio):
    # Necessary V_out^2 - ratio^2 V_in^2 <=0. The worst possible variance
    # given two marginal standard deviations is the square of their sum.
    u, su, du = np.array(inner["velocity_stat_sys"]).T
    v, sv, dv = np.array(outer["velocity_stat_sys"]).T
    q = ratio**2
    contrast = v**2 - q * u**2
    syst = 2 * (v * dv + q * u * du)
    sigma = 2 * (v * sv + q * u * su)
    return dict(
        observed_speed_ratios=(v / u).tolist(),
        quadratic_excess_after_full_systematics=(contrast - syst).tolist(),
        minimum_z_under_any_radial_correlation=((contrast - syst) / sigma).tolist(),
        above_bonferroni_four_one_sided_95=(
            (contrast - syst) > norm.ppf(1 - 0.05 / 4) * sigma
        ).tolist(),
    )


def run():
    rows = json.loads((HERE / "inputs.json").read_text())["rows"]
    boards = json.loads((native.HERE / "results.json").read_text())["boards"]
    outer = rows[-1]
    main = next(row for row in rows if row["radius_kpc"] >= 50)
    results = []
    for inner in rows:
        if inner["radius_kpc"] > 150:
            continue
        s = outer["radius_kpc"] / inner["radius_kpc"]
        for board in boards:
            for layout in ("geodesic", "alternative"):
                for f in (1.0, 0.75, 0.5):
                    bounds = envelope(
                        inner["radius_kpc"], outer["radius_kpc"], board, layout, f
                    )
                    coherent = bounds["all_free_times"]["speed_ratio_ceiling"]
                    n = len(board[layout]["midpoints"])
                    q = 8.5 * n / board["initial"]
                    results.append(
                        dict(
                            inner_kpc=inner["radius_kpc"],
                            outer_kpc=outer["radius_kpc"],
                            primary=inner == main,
                            owners=board["owners"],
                            layout=layout,
                            inner_baryon_fraction=f,
                            uniform_continuous_trail_ceiling=np.sqrt(
                                (1 + q) / (s * f + q)
                            ),
                            limits=bounds,
                            comparison=observations(inner, outer, coherent),
                        )
                    )
    return dict(
        primary_radii_kpc=[main["radius_kpc"], outer["radius_kpc"]], cases=results
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    if args.fetch:
        with urlopen(lens.URL, timeout=30) as f:
            raw = f.read()
        (HERE / "inputs.json").write_text(json.dumps(acquire(raw), indent=2) + "\n")
    result = run()
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        native.check(json.loads(path.read_text()), result)
