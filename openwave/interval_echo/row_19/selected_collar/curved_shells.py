"""Recover spherical source masses from two independently acquired aperture entropies."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import pi, sqrt

import networkx as nx
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from . import midpoint_bulk as bulk

HERE = bulk.HERE
SEED = 28122864


def nonlinear_area(r, centers, energy, step):
    # Use areal radius and time-symmetric Einstein constraint. Exact
    # shell jumps are in gravitational mass, not in coordinate density.
    shell_r = np.sqrt(1 - centers**2)
    mass = centers * energy

    def volume(radius, eps):
        edges = [0.0, *(x for x in shell_r if 0 < x < radius), radius]
        total = 0.0
        for lower, upper in zip(edges[:-1], edges[1:], strict=True):
            cumulative = sum(mass[shell_r <= (lower + upper) / 2])

            def element(x):
                f = 1 - x * x - 2 * eps * cumulative / x
                return 4 * pi * x * x / sqrt(f)

            total += quad(element, lower, upper, epsabs=1e-13, epsrel=3e-13)[0]
        return total

    standard = volume(r, 0)
    out = []
    for eps in (step, -step):
        bracket = min((1 - r) / 3, r / 10)
        match = brentq(
            lambda x: volume(x, eps) - standard, r - bracket, r + bracket, xtol=2e-15
        )
        out.append(
            {
                "coupling_times_mass_scale": eps,
                "matched_areal_radius": match,
                "volume_difference": volume(match, eps) - standard,
                "area": 4 * pi * match**2,
            }
        )
    return out, (out[0]["area"] - out[1]["area"]) / (2 * step)


def clock_aperture(prob, laws, rng):
    """Read the full native clock/energy alphabet after an address-only gate."""
    n = 2**51
    iv = bulk.shells.iv
    counts, conditional = [], []

    # At most 3*82*10 conditional or admission probabilities on the ten
    # apertures. Conditional on the accepted starts, the words are iid.
    def uncertainty(size):
        return iv.mpf([-1, 1]) * iv.sqrt(iv.log(2 * 3000 / iv.mpf(".01")) / (2 * size))

    for law in laws:
        admitted = int(rng.binomial(n, float(prob.mid)))
        row = rng.multinomial(admitted, law / 2592).tolist()
        counts.append([*row, n - admitted])
        conditional.append([iv.mpf(k) / admitted + uncertainty(admitted) for k in row])
    admitted = sum(counts[0][:-1])
    gate = iv.mpf(admitted) / n + uncertainty(n)
    price = -sum((s - f) * iv.log(81 * b) for b, f, s in zip(*conditional, strict=True))
    return {"counts": counts, "charge": bulk.shells.bounds(gate * price)}


def run():
    iv = bulk.shells.iv
    iv.dps = 35
    rng = np.random.default_rng(SEED)
    board = bulk.federation(0)
    graphs, points = bulk.geometry._atlas([board])
    hop = nx.floyd_warshall_numpy(graphs[0], nodelist=range(board.cells)).astype(int)
    pairs = np.triu_indices(board.cells, 1)
    separation = hop[pairs]
    dots = [iv.mpf(v) / 3 for v in (3, iv.sqrt(5), 1, -1, -iv.sqrt(5), -3)]
    centers = [iv.sqrt((1 - dots[k]) / 2) for k in (5, 4, 3, 2, 1)]
    shell = np.array([5 - k for k in separation])
    counts = np.bincount(shell)
    angular = np.sqrt(
        np.maximum(
            0, 1 - ((points[0][pairs[0]] + points[0][pairs[1]]) ** 2).sum(axis=1) / 4
        )
    )
    np.testing.assert_allclose(
        angular, np.array([float(c.mid) for c in centers])[shell], atol=2e-15
    )
    src = bulk.shells.dist.read("resolved_shells")
    price = iv.mpf(src["conditional_entropy_slope"])
    word = [F(src["word_event_probabilities"][a]) for a in ("fresh", "shared")]

    def read(p):
        return bulk.native_read(F(float(p.mid)), price, word, rng)

    from . import paired_ranging as survey

    detector = json.loads((HERE / "pair_clock_information.json").read_bytes())[
        "local_reader"
    ]
    symbols = detector["calibrations"][detector["primary"]][
        "sharp_records_in_calibration_order"
    ]
    alphabet = [tuple(z) for z in detector["alphabet"]]
    labels = [alphabet.index(tuple(a + b)) for a in symbols for b in symbols]
    laws = [
        np.bincount(labels, weights=w.ravel()).astype(int)
        for w in survey.standards(board)[0][:3]
    ]
    clock_rng = np.random.default_rng(SEED + 700)
    full = read(iv.mpf(1))
    # Express both measured masses as fractions of the unfiltered full
    # price divided by 2π L. External unit conversion cancels in these reads.
    norm = iv.mpf(full["charge"])
    regimes = {}
    for label in ("complete", "empty_center"):
        w = [
            iv.mpf(int(n)) / len(shell) if label == "complete" or j else iv.mpf(0)
            for j, n in enumerate(counts)
        ]
        outputs = []
        for i, (c, d) in enumerate(
            zip(centers, [*centers[1:], iv.mpf(0)], strict=True)
        ):
            hi = (2 * c + d) / 3
            lo = (c + 2 * d) / 3
            aperture = []
            y = []
            clock_y = []
            for cut in (hi, lo):
                sn = iv.sqrt(1 - cut**2)
                prob = sum(
                    e * (v - cut) / sn
                    for e, v in zip(w[: i + 1], centers[: i + 1], strict=True)
                )
                obs = read(prob)
                aperture.append({"cosine": bulk.shells.bounds(cut), "word": obs})
                y.append(sn * iv.mpf(obs["charge"]) / norm)
                if label == "complete":
                    local = clock_aperture(prob, laws, clock_rng)
                    aperture[-1]["clock_and_energy"] = local
                    clock_y.append(sn * iv.mpf(local["charge"]) / norm)
            mass = (hi * y[1] - lo * y[0]) / (hi - lo)
            proper = (y[1] - y[0]) / (hi - lo)
            exact = sum(
                a * b for a, b in zip(w[: i + 1], centers[: i + 1], strict=True)
            )
            # Independent direct address and word read: roll a fresh tag
            # under the measured shell cosine before receiving the word.
            addr = read(exact)
            independent = iv.mpf(addr["charge"]) / norm
            target = (hi + lo) / 2
            s = iv.sqrt(1 - target**2)
            forecast = -8 * iv.pi * ((y[1] + y[0]) / 2) / s
            ideal = -8 * iv.pi * (exact - target * sum(w[: i + 1])) / s
            result = {
                "cuts": aperture,
                "normalized_mass_read": bulk.shells.bounds(mass),
                "proper_mass_read": bulk.shells.bounds(proper),
                "direct_address_read": addr,
                "direct_normalized_mass": bulk.shells.bounds(independent),
                "population_mass": bulk.shells.bounds(exact),
                "squared_areal_radius": bulk.shells.bounds(s**2),
                "area_slope_per_G_mass_unit": bulk.shells.bounds(forecast),
                "population_area_slope": bulk.shells.bounds(ideal),
                "inward_acceleration_per_G_mass_unit": bulk.shells.bounds(mass / s**2),
            }
            if clock_y:
                result["local_detector_normalized_mass"] = bulk.shells.bounds(
                    (hi * clock_y[1] - lo * clock_y[0]) / (hi - lo)
                )
            if label == "empty_center":
                curves = []
                for step in (2**-14, 2**-16, 2**-18):
                    test, slope = nonlinear_area(
                        float(s.mid),
                        np.array([float(z.mid) for z in centers[1:]]),
                        np.array([float(z.mid) for z in w[1:]]),
                        step,
                    )
                    curves.append(
                        {
                            "step": step,
                            "finite_volume_checks": test,
                            "centered_area_slope": slope,
                        }
                    )
                result["nonlinear_spatial_constraint"] = curves
            outputs.append(result)
        regimes[label] = outputs
    cosmic_file = "../../../../code/cosmology/cold_dilution/results.json"
    observed = json.loads((HERE / cosmic_file).read_bytes())["constant_exponent_test"][
        "density_ratio_for_doubling_99_7pct"
    ]
    f1 = iv.mpf(regimes["complete"][-1]["normalized_mass_read"])
    density = iv.mpf(observed)
    # The weighted mean cosine ratio f2/f1 multiplies K2/(16 K1).
    cosmic = {
        "observed_density_ratio": observed,
        "initial_cosine_interval": bulk.shells.bounds(f1),
        "largest_conserved_price_density_ratio": bulk.shells.bounds(1 / (16 * f1.a)),
        "minimum_price_growth_measured_profile": bulk.shells.bounds(
            16 * density.a * f1.a
        ),
        "minimum_price_growth_isotropic": bulk.shells.bounds(16 * density.a * 2 / 3),
    }
    return {
        "cosmic_doubling": cosmic,
        "seed": SEED,
        "full_price_read": full,
        "shell_multiplicities": counts.tolist(),
        "cosine_shell_locations": [bulk.shells.bounds(c) for c in centers],
        "regimes": regimes,
        "inputs": {
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in (
                "CURVED_SHELLS.md",
                "resolved_shells.json",
                "paired_ranging.json",
                "pair_clock_information.json",
                cosmic_file,
            )
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    z = run()
    target = HERE / "curved_shells.json"
    if args.write:
        target.write_text(json.dumps(z, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(target.read_bytes()), z)
    for tag, g in z["regimes"].items():
        print(tag, [r["normalized_mass_read"] for r in g])
