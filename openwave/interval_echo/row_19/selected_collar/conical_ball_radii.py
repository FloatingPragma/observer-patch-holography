"""Infer the small-radius scaling of equal-volume ball areas from observed clock rates."""

import argparse
from hashlib import sha256
import json

import numpy as np

from . import regional_density as regions

HERE, lab, reg, iv = regions.HERE, regions.lab, regions.reg, regions.iv


def run():
    iv.dps = 35
    path = HERE / "six_meter_clocks.json"
    records = json.loads(path.read_bytes())
    source = HERE / "six_meter_pointer.json"
    assert sha256(source.read_bytes()).hexdigest() == records["inputs"][source.name]
    num, den = json.loads(source.read_bytes())["common_ratio"]
    beta = iv.mpf(num) / den
    dose = regions.density.DOSE
    means, _ = reg.timing_slopes(
        records["clock_counts"], records["starts"], (("held", "dense", 1),)
    )
    means = {k: np.array(v, dtype=object) for k, v in means.items()}
    baseline = (1 - dose) * means["reference"] + dose * means["dense"]
    ordinary = np.array([beta * v for v in dose * (means["ordinary"] - means["dense"])])
    difference = dose * (means["held"] - means["dense"]) - ordinary
    rates = regions.ball_area_rates(baseline + ordinary, difference, tangent=True)
    rows = []
    for (tri, vertex), rate in rates.items():
        # One measured angular derivative enters all three radii. Preserve
        # this dependence when cancelling the nuisance terms.
        rows.append(
            {
                "triangle": list(tri),
                "vertex": vertex,
                "fixed_volume_area_rates": [
                    lab.six.bounds(rate / (4**j)) for j in range(3)
                ],
                "half_minus_full_area_over_radius_fourth": lab.six.bounds(48 * rate),
                "affine_source_contrast": lab.six.bounds(3 * rate / 8),
                "minimum_maximum_residual_lower": lab.six.bounds(abs(rate) / 16)[0],
            }
        )
    return {
        "radii": [0.5, 0.25, 0.125],
        "ordinary_fraction": [num, den],
        "balls": rows,
        "inputs": {
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (path, source, HERE / "CONICAL_BALL_RADII.md")
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    path = HERE / "conical_ball_radii.json"
    if args.write:
        path.write_text(json.dumps(result, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), result)
    rows = result["balls"]
    print(
        "resolved radius contrasts",
        sum(r["minimum_maximum_residual_lower"] > 0 for r in rows),
        "/",
        len(rows),
    )
    print("largest lower bound", max(r["minimum_maximum_residual_lower"] for r in rows))
