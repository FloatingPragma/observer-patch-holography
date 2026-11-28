"""Test the disk-calibrated cut-response amplitude against a disjoint shear stack."""

import argparse
from hashlib import sha256
import json

import numpy as np
from scipy.special import gammainc, gammaincinv

from . import fluctuation_cuts as cuts
from cosmology.rar_deep_regime import stellar_lensing as lens

HERE = cuts.HERE


def envelope(p):
    # Concave majorant matching the quartic with endpoint tangents at 1/3 and 2/3.
    return np.where(
        p < 1 / 3,
        4 * p / 27,
        np.where(p > 2 / 3, 4 * (1 - p) / 27, p**2 * (1 - p) ** 2),
    )


def run():
    physics = HERE / "fluctuation_cuts.json"
    results = json.loads(physics.read_bytes())
    raw = lens.HERE / "data/slacs_weak_lenses.dat"
    fields = np.loadtxt(raw, dtype=str, skiprows=2)
    z, sigma, re, eff, ch, ec, sp, es = fields[:, 1:].astype(float).T
    assert len(set(fields[:, 0])) == 22
    weight = eff**2 / sum(eff**2)
    shearfile = lens.RECEIPT
    bins = json.loads(shearfile.read_bytes())["outer_shear_same_22_lenses"][
        "bins_in_h_Msun_pc2"
    ]
    h = 0.7
    x, y = np.array([r["centre_hinv_kpc"] / h for r in bins])
    div = 1e6 * h
    factor = 4 * cuts.profiles.galaxies.G * div / (y - x)

    def contrast(key):
        a, b = (r[key] for r in bins)
        center = factor * (y * y * b - x * x * a)
        width = factor * (
            x * x * bins[0]["quoted_error"] + y * y * bins[1]["quoted_error"]
        )
        return dict(
            central=float(center),
            one_error=[float(center - width), float(center + width)],
            two_errors=[float(center - 2 * width), float(center + 2 * width)],
        )

    peakfile = HERE / "record_profile_maps.json"
    peak = json.loads(peakfile.read_bytes())["peak_msun_pc2"]["mapped_baryons"]
    book_unit = (
        results["possible_addresses"]
        * cuts.rp.galaxy.WIDTH**2
        * cuts.rp.maps.AREA_UNITS
        * peak
    )
    # Integrate the observed stellar sizes over the same physical area per
    # slot. The first slot covers the core and never samples a singular center.
    dl = lens.C_KM_S / 70 * 1000 * np.array([lens.comoving(x) for x in z]) / (1 + z)
    re_kpc = re * dl / lens.ARCSEC_PER_RADIAN
    edges = cuts.rp.galaxy.WIDTH * np.sqrt(
        np.arange(results["possible_addresses"] + 1) / np.pi
    )
    cumulative = gammainc(
        8, gammaincinv(8, 0.5) * (edges[None, :] / re_kpc[:, None]) ** 0.25
    )
    fractions = np.diff(cumulative, axis=1)
    assert fractions.min() >= 0
    projected = {}
    for cap in cuts.rp.CAPACITY:
        for imf, mass in (("Chabrier", 10**ch), ("Salpeter", 10**sp)):
            f = (
                mass[:, None]
                * fractions
                / (book_unit * cap / results["possible_addresses"])
            )
            bad = f.max(axis=1) > 1
            projected[(cap, imf)] = (
                fields[bad, 0].tolist(),
                f.mean(1) / 4,
                6 * (448 / 2110) * (f**2 * (1 - f) ** 2).mean(1) / 4,
            )
    report = {}
    for key, row in results["panels"].items():
        tracer, cap, L, label = key.split("/")
        if tracer != "mapped_baryons" or label not in (
            "single",
            "middle_original",
            "middle_repaired",
            "middle_increment",
        ):
            continue
        n = 0 if L == "infinite" else int(L)
        gain = row["all_disk_least_squares_gain"]
        capacity = float(cap) * book_unit
        pmin, pmax = 10 ** (ch - ec) / capacity, 10 ** (sp + es) / capacity
        assert min(pmin) > 0 and max(pmax) <= 1
        if label == "single":
            report[key] = dict(
                chabrier=float(gain * (weight @ cuts.rate(10**ch / (4 * capacity), n))),
                salpeter=float(gain * (weight @ cuts.rate(10**sp / (4 * capacity), n))),
                mass_box=[
                    float(gain * (weight @ cuts.rate(p / (4), n))) for p in (pmin, pmax)
                ],
            )
        else:
            ystar = np.clip(0.5, pmin, pmax)
            eligible_cap = 6 * (448 / 2110) * envelope(ystar)
            before, after = (eligible_cap * f / 4 for f in (3 / 28, 5 / 28))
            if label == "middle_increment":
                # The finite response E[1-(1-h)^T] is concave and increasing.
                # Hence R(b)-R(a) <= (1-a/b)R(b). The long-ring bound is exact.
                top = (
                    (1 - 3 / 5) * cuts.rate(after, n)
                    if n
                    else np.sqrt(after) - np.sqrt(before)
                )
            else:
                top = cuts.rate(before if label == "middle_original" else after, n)
            report[key] = dict(
                mass_and_placement_box=[0.0, float(gain * (weight @ top))]
            )
        rates = {}
        for imf in ("Chabrier", "Salpeter"):
            rejected, ordinary, signal = projected[(int(cap), imf)]
            if rejected:
                rates[imf] = dict(over_capacity=rejected)
                continue
            a, b = (cuts.rate(signal * t, n) for t in (3 / 28, 5 / 28))
            predicted = {
                "single": cuts.rate(ordinary, n),
                "middle_original": a,
                "middle_repaired": b,
                "middle_increment": b - a,
            }[label]
            rates[imf] = dict(
                predicted_squared_speed=float(gain * (weight @ predicted))
            )
        report[key]["photometric_projection"] = rates
    return dict(
        missed_light_fraction_range=[
            float(min(1 - cumulative[:, -1])),
            float(max(1 - cumulative[:, -1])),
        ],
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (physics, raw, shearfile, peakfile, HERE / "FLUCTUATION_CUTS.md")
        },
        mass_capacity_scale_msun=book_unit,
        radii_physical_kpc=[float(x), float(y)],
        velocity_squared_unit="(km/s)^2",
        tangential=contrast("tangential"),
        cross=contrast("cross"),
        disk_calibrated=report,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    result = run()
    path = HERE / "fluctuation_lensing.json"
    if cli.parse_args().write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        cuts.rp.galaxy.check(json.loads(path.read_bytes()), result)
    print(result["tangential"], result["cross"])
