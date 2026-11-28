"""Compare the finite equilibrium position work with resolved outer galaxy speeds."""

import argparse
from hashlib import sha256
import json
from math import log

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    thermal_pair as lab,
)
from . import kernel_budget as budget

HERE = budget.HERE


def integral(rows):
    # v_extra^2 d log(r) = d Phi under the stipulated circular, central reduction.
    return [
        sum(
            log(v["kpc"] / u["kpc"])
            * (v["extra_squared_kms"][i] + u["extra_squared_kms"][i])
            / 2
            for u, v in zip(rows[:-1], rows[1:], strict=True)
        )
        for i in range(2)
    ]


def positive_ratio(num, den):
    assert min(den) > 0 and min(num) > 0
    return [num[0] / den[1], num[1] / den[0]]


def run():
    source = lab.HERE / "thermal_pair.json"
    work = json.loads(source.read_bytes())["certified_work"]["enrolled"]
    inputs, *_ = budget.load_inputs()
    galaxy = inputs["NGC2841"]
    low, high = budget.bounds(galaxy)
    rows = [
        dict(
            kpc=r,
            observed_kms=v,
            error_kms=e,
            baryons_kms=b,
            extra_squared_kms=[float(a), float(z)],
        )
        for r, v, e, b, a, z in zip(
            galaxy["radius_kpc"],
            galaxy["observed_kms"],
            galaxy["error_kms"],
            galaxy["baryons_kms"],
            low,
            high,
            strict=True,
        )
        if 20 <= r <= 64
    ]
    i, j, k = [
        min(range(len(rows)), key=lambda n: abs(rows[n]["kpc"] - target))
        for target in (20.5, 41.0, 61.5)
    ]
    inner, outer = integral(rows[i : j + 1]), integral(rows[j : k + 1])
    radial = positive_ratio(outer, inner)
    reverse = {key: [-b, -a] for key, (a, b) in work.items()}
    shape = positive_ratio(reverse["4_6"], reverse["2_4"])
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                source,
                HERE / "THERMAL_WORK.md",
                HERE / "THERMAL_SHAPE.md",
                HERE / "source_profile.json",
            )
        },
        galaxy="NGC2841",
        orbits=rows,
        all_work_km2_s2=integral(rows),
        matched_radii_kpc=[rows[n]["kpc"] for n in (i, j, k)],
        inner_work_km2_s2=inner,
        outer_work_km2_s2=outer,
        required_work_ratio=radial,
        model_work=work,
        model_absolute_ratio=shape,
        minimum_shape_discrepancy=radial[0] / shape[1],
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    path = HERE / "thermal_work.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_bytes()), result
        )
    print({k: v for k, v in result.items() if k not in ("orbits", "sha256")})
