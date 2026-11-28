"""Required recovery contrast if lost, flagged seams carry all inferred dark mass."""

import argparse
from hashlib import sha256
import json

import numpy as np

from . import measure as source, tidal

HERE = source.HERE
COSMOLOGY = HERE.parents[3] / "code" / "cosmology" / "record_abundance"
RECEIPT = HERE / "retention_budget.json"


def run():
    files = [
        source.RECEIPT,
        COSMOLOGY / "mass_ratio.json",
        COSMOLOGY / "clustering.json",
    ]
    measured, budget, clusters = [json.loads(file.read_text()) for file in files]
    gates = {}
    for key, arm in measured["gates"].items():
        g = float(key)
        m = source.decode(arm["counts"], arm["pauli_coefficients"], g)["flagged"]
        expected = (1 - g) * measured["gates"]["0.0"]["expected"]["flagged"]
        np.testing.assert_allclose(arm["expected"]["flagged"], expected, atol=1e-15)
        gates[key] = dict(
            predicted=expected, read=m["estimate"], half_width=m["halfwidth"]
        )
    p = budget["p"]
    cosmic = budget["external_2018_cold_to_baryons"]

    def required(stock):
        ratio = (
            stock["dynamical_msun"] - stock["catalogued_stars_and_gas_msun"]
        ) / stock["catalogued_stars_and_gas_msun"]
        return ratio * (1 - p) / cosmic

    gain = required(budget["ngc3741"])
    cluster = clusters["observational_comparison"]
    cluster_max = (
        cluster["dark_density_cap"]
        / cluster["r500_nonbaryonic_density_over_critical_at_upper_baryon_fraction"]
    )
    return dict(
        inputs_sha256={f.name: sha256(f.read_bytes()).hexdigest() for f in files},
        measured_flagged_scores=gates,
        occupation=p,
        externally_anchored_cosmic_ratio=cosmic,
        galaxy_min_failed_fraction_gain=gain,
        galaxy_max_cosmic_failed_fraction=1 / gain,
        galaxy_gain_at_distance_endpoints=[
            required(row) for row in budget["distance_sensitivity"]
        ],
        cluster_max_cosmic_failed_fraction_fixed_sites=cluster_max,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    d = run()
    if args.write:
        with RECEIPT.open("x") as out:
            json.dump(d, out, indent=2, allow_nan=False)
            out.write("\n")
    else:
        tidal.verify(d, json.loads(RECEIPT.read_text()))
    print(
        d["galaxy_min_failed_fraction_gain"],
        d["galaxy_max_cosmic_failed_fraction"],
        d["cluster_max_cosmic_failed_fraction_fixed_sites"],
    )
