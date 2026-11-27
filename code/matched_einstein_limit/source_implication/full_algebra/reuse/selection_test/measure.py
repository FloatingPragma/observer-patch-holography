"""Compare the volume-matched uniform reuse equilibrium and the area-entropy optimum."""

import argparse
import json
from fractions import Fraction
from math import log, sqrt
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from ... import station
from ....material_clocks.measure import density_center
from ..retained_equilibrium.certificate import channels

HERE = Path(__file__).resolve().parent
DATA = HERE.parents[1]
CYCLES = 256
SHOTS = 2**14
SEED = 8026


def box(value):
    return tuple(Fraction(s.strip()) for s in value.strip("[]").split(","))


def record():
    def read(filename):
        return json.loads((DATA / filename).read_text(), parse_float=Fraction)

    theorem = read("certificate.json")
    reuse = read("reuse_equilibrium.json")
    optimizer_zz = box(read("pulse_bath.json")["matched_measurement"]["optimizer_joint_ZZ"])
    reuse_zz = box(reuse["joint_ZZ"])
    gain = box(theorem["calibrated_gain"])
    target = box(theorem["calibrated_volume"])
    overlap = box(reuse["target_volume"])
    if max(target[0], overlap[0]) > min(target[1], overlap[1]):
        raise ArithmeticError("The calibration and reuse volume enclosures disagree")
    co = theorem["global_survey"]["feasible_objective_gap_coefficient"]
    separation = optimizer_zz[0] - reuse_zz[1]
    population_floor = co * separation**2
    if population_floor <= Fraction("0.015834"):
        raise ArithmeticError("Certified selection deficit is too small")

    # Trace over the pre-deadline record; average both charges and all directions.
    bath = channels().sum(axis=3).mean(axis=(0, 2))
    sharp = np.diag([1, 0, 0, 0]).ravel()
    sectors = np.stack([(np.linalg.matrix_power(m, CYCLES) @ sharp).reshape(4, 4) for m in bath])
    sectors = (sectors + sectors.conj().transpose(0, 2, 1)) / 2
    sectors = sectors / sectors.trace(axis1=1, axis2=2)[:, None, None]
    geometry = station.radar_menus()[6]

    def state(odd):
        return sectors * np.array([1 - odd, odd])[:, None, None]

    def survey(x):
        return station.geometry(station.pair(geometry, x))

    nominal_volume = float(sum(target) / 2)
    fitted = brentq(lambda w: survey(state(w))[1] - nominal_volume, 0.134, 0.1346, xtol=1e-14)
    uniform = state(fitted)
    optimum, prepared_odd, _, _ = density_center()
    zz = np.diag([1, -1, -1, 1])
    means = [float(np.einsum("bij,ji->", x, zz).real) for x in (optimum, uniform)]
    for val, certified in zip(means, (optimizer_zz, reuse_zz), strict=True):
        if not float(certified[0]) <= val <= float(certified[1]):
            raise ArithmeticError("Independent native computation left the witness enclosure")
    areas_volumes = [survey(x) for x in (optimum, uniform)]
    eta = float(sum(gain) / 2)
    area_change, volume_change = areas_volumes[0] - areas_volumes[1]
    entropy_change = station.entropy(optimum) - station.entropy(uniform)

    rates = (1 + np.array(means)) / 2
    counts = np.random.default_rng(SEED).binomial(SHOTS, rates)
    radius = sqrt(log(4 / 0.01) / (2 * SHOTS))
    # Approximation-to-exact bias is covered by the complete certified widths.
    bias = sum(float(hi - lo) for lo, hi in (optimizer_zz, reuse_zz))
    difference_lower = 2 * ((counts[0] - counts[1]) / SHOTS - 2 * radius) - bias
    return {
        "shots_per_cohort": SHOTS, "seed": SEED, "uncontrolled_cycles": CYCLES,
        "product_plus_counts_optimum_then_reuse": counts.tolist(),
        "numerical_odd_fraction": fitted, "variational_odd_fraction_center": float(prepared_odd),
        "numerical_ZZ_means": means, "numerical_area_and_volume": [x.tolist() for x in areas_volumes],
        "numerical_entropy_difference": entropy_change,
        "numerical_area_difference": area_change,
        "numerical_volume_difference": volume_change,
        "numerical_objective_deficit": entropy_change + eta * area_change,
        "exact_population_mean_separation_lower": str(separation),
        "exact_population_deficit_lower_nats": str(population_floor),
        "sampling_confidence": 0.99, "rate_radius": radius,
        "certified_mean_bias_allowance": bias,
        "counts_objective_deficit_lower_nats": float(co) * max(0, difference_lower)**2,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = record()
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    else:
        saved = json.loads(path.read_text())
        for name, value in data.items():
            if isinstance(value, str) or name == "product_plus_counts_optimum_then_reuse":
                if value != saved[name]:
                    raise ArithmeticError(f"Selection witness differs for {name}")
            else:
                np.testing.assert_allclose(value, saved[name], rtol=1e-9, atol=2e-12)
    print("Verified the joint-product bound for the uncontrolled radar fixed point.")


if __name__ == "__main__":
    main()
