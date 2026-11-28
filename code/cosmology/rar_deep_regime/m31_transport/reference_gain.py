"""Can changing the finite equilibrium reference repair the frozen M31 aperture ratio?"""

import argparse
import json
from pathlib import Path

import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.closed_region import (
    reservoir,
)

HERE = Path(__file__).parent


def run():
    source = json.loads((HERE / "aperture.json").read_text())["surveys"]
    tested = reservoir.run()["contact_bank"]
    valid = [r for r in tested if not r["singular"]]
    # Outer versus inner means of any positive conversion from this family
    # cannot exceed its maximum divided by its minimum, for any radial mixing.
    factors = [1 + np.array(r["fractional_error_billionths"]) * 1e-9 for r in valid]
    gain_ceiling = float(max(x[1] for x in factors) / min(x[0] for x in factors))
    common_scale = [
        (16 * r["beta_equal_to_even"] + r["beta_even_to_odd"]) / 17 for r in valid
    ]
    result = {}
    for survey, panel in source.items():
        cases = {}
        for name, c in panel["source_cases"].items():
            ratio = c["upper_mass_ratio_all_lifetimes"]
            observed = c["required_mass_ratio_nominal_one_two_error_boxes"]
            cases[name] = dict(
                bounds_on_mean_outer_to_inner_gain=[
                    max(0, r - 1) / (ratio - 1) for r in observed
                ],
                maximum_aperture_mass_ratio_with_reference_variation=1
                + gain_ceiling * (ratio - 1),
            )
        result[survey] = dict(
            endpoint_radii_kpc=panel["endpoint_radii_kpc"],
            cases=cases,
            least_two_error_gain_over_tested_sources=min(
                v["bounds_on_mean_outer_to_inner_gain"][2] for v in cases.values()
            ),
        )
    return dict(
        tested_faithful_contacts=[r["contacts"] for r in valid],
        excluded_singular=[r["contacts"] for r in tested if r["singular"]],
        maximum_gain_ratio_with_local_ordinary_calibration=gain_ceiling,
        maximum_gain_ratio_with_one_overall_calibration=max(common_scale)
        / min(common_scale),
        surveys=result,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    path = HERE / "reference_gain.json"
    data = run()
    if args.write:
        path.write_text(json.dumps(data, indent=2) + "\n")
    else:
        from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.closed_region.measure import (
            check,
        )

        check(json.loads(path.read_text()), data)
