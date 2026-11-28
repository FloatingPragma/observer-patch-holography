"""Compare observer chord tides with the SPARC in-plane circular acceleration."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json

import numpy as np

from . import moving_budget as budget
from . import release_rotation as rotation

HERE = budget.HERE


def catalogue_tides(meta, rows):
    r, v, e = [
        np.array([z[k] for z in rows]) for k in ("rad_kpc", "vobs_kms", "e_vobs_kms")
    ]
    columns = rotation.signed_columns(rows)
    nominal = columns @ [1.0, 0.5, 0.7]
    upper = np.maximum(
        columns * np.array([0.7, 0.3, 0.4]), columns * np.array([1.5, 0.8, 1.0])
    ).sum(1)
    inc = meta["inclination_deg"]
    k = (1 + 2 * meta["e_dist_mpc"] / meta["dist_mpc"]) * (
        np.sin(np.radians(min(90, inc + 2 * meta["e_inclination_deg"])))
        / np.sin(np.radians(inc))
    ) ** 2
    # The catalogue speed uses its nominal distance and inclination. Dividing
    # both ends by a common positive projection cannot affect relative width.
    endpoints = np.stack((np.maximum(v - 2 * e, 0) ** 2, (v + 2 * e) ** 2))
    return {
        k: (endpoints - bare) / r**2
        for k, bare in [("nominal", nominal), ("favorable", upper * k)]
    }


def precision_test(native):
    meta, groups = rotation.galaxies.catalogue()
    data = {n: catalogue_tides(meta[n], r) for n, r in groups.items()}
    ring_total = sum(z["nominal"].shape[1] for z in data.values())
    out = {}
    for order, p in native["schedules"].items():
        for label, v in p["clock_freedom"].items():
            for clock in ("arbitrary_intervals", "equal_intervals"):
                eps = float(Q(v[clock]["error"]))
                panels = {}
                for box in ("nominal", "favorable"):
                    excluded = {}
                    for n, tab in data.items():
                        lo, hi = tab[box]
                        # Necessary and sufficient for fitting one positive
                        # common tidal band in the closed clock cone. Count
                        # only orbits that require strictly inward excess.
                        miss = (lo > 0) & ((hi - lo) < eps * (hi + lo))
                        excluded[n] = int(miss.sum())
                    panels[box] = dict(
                        excluded_rings=sum(excluded.values()),
                        galaxies_with_exclusions=sum(v > 0 for v in excluded.values()),
                        counts=excluded,
                    )
                out[f"{order}/{label}/{clock}"] = panels
    # A specific public radius gives the units and sign of the comparison.
    n = "DDO154"
    example = dict(
        galaxy=n,
        radius_kpc=groups[n][-1]["rad_kpc"],
        inward_tide_interval_kms2_per_kpc2=data[n]["nominal"][:, -1].tolist(),
    )
    return dict(rings=ring_total, example=example, comparisons=out)


def evaluate_maps(factor, maps, inputs):
    G = budget.inventory.rp.galaxy.transfer.budget.local.profiles.galaxies.G
    panels = {}
    for tracer, menu in maps.items():
        for cap, panel in menu.items():
            all_counts = {}
            invalid = 0
            for name, stock in panel["closed"].items():
                tab = inputs[name]
                r, v, e, b = [
                    np.array(tab[key])
                    for key in (
                        "radius_kpc",
                        "observed_kms",
                        "error_kms",
                        "baryons_kms",
                    )
                ]
                # A spherical enclosed excess gives kappa=GM/r³ on the chord.
                # Add to the actual disk's baryonic circular acceleration.
                speed2 = b * b + G * float(factor) * np.asarray(stock) / r
                ok = (speed2 >= np.maximum(v - 2 * e, 0) ** 2) & (
                    speed2 <= (v + 2 * e) ** 2
                )
                all_counts[name] = [int(ok.sum()), len(ok)]
                invalid += int((speed2 < 0).sum())
            panels[f"{tracer}/{cap}"] = dict(
                within_two=sum(v for v, n in all_counts.values()),
                complete=sum(v == n for v, n in all_counts.values()),
                negative_speed_squared=invalid,
                counts=all_counts,
            )
    return panels


def radial_allowance(maps, inputs, measured):
    """Grant arbitrary annulus-by-annulus gains within *all* observed prices."""
    G = budget.inventory.rp.galaxy.transfer.budget.local.profiles.galaxies.G
    low_gain, high_gain = min(0, min(measured)), max(0, max(measured))
    panels = {}
    for tracer, menu in maps.items():
        for cap, bank in menu.items():
            disks = {}
            for name, profile in bank["closed"].items():
                r, v, e, b = [
                    np.asarray(inputs[name][k])
                    for k in ("radius_kpc", "observed_kms", "error_kms", "baryons_kms")
                ]
                lower = r * (np.maximum(0, v - 2 * e) ** 2 - b * b) / G
                upper = r * ((v + 2 * e) ** 2 - b * b) / G
                stock = np.r_[0, profile]
                increments = np.diff(stock)
                assert np.all(increments >= 0)
                lo = hi = 0.0
                first = None
                for j, amount in enumerate(increments):
                    lo = max(lower[j], lo + amount * low_gain)
                    hi = min(upper[j], hi + amount * high_gain)
                    if lo > hi:
                        first = float(r[j])
                        break
                # Pairwise necessary bound: what extra mass must fit after
                # each earlier admissible upper edge, including the origin?
                need = 0.0
                witness = None
                ceiling = np.r_[0, upper]
                for k in range(1, len(stock)):
                    for j in range(k):
                        mass = lower[k - 1] - ceiling[j]
                        if mass <= 0:
                            continue
                        amount = stock[k] - stock[j]
                        rate = float(mass / amount) if amount else float("inf")
                        if rate > need:
                            need = rate
                            witness = [
                                0.0 if j == 0 else float(r[j - 1]),
                                float(r[k - 1]),
                            ]
                envelope = (lower <= high_gain * stock[1:]) & (
                    upper >= low_gain * stock[1:]
                )
                disks[name] = dict(
                    first_infeasible_kpc=first,
                    individually_compatible=int(envelope.sum()),
                    necessary_peak_gain=need if np.isfinite(need) else "unbounded",
                    witness_radii_kpc=witness,
                )
            panels[f"{tracer}/{cap}"] = dict(
                disks=disks,
                complete=sum(v["first_infeasible_kpc"] is None for v in disks.values()),
                compatible_rings=sum(
                    v["individually_compatible"] for v in disks.values()
                ),
            )
    return dict(allowed_gain=[low_gain, high_gain], by_map=panels)


def exterior_profiles(research, maps, inputs, meta, single_mass):
    radii = np.array(research["registered_radii"], float)
    G = budget.inventory.rp.galaxy.transfer.budget.local.profiles.galaxies.G
    measured = {}
    for src, donor in [("source", "single"), ("positive", "reference_plus_single")]:
        signal = []
        ordinary = []
        for radius in radii:
            row = research["shells"][str(int(radius))]
            data = row["comparisons"][src]["distant_equal_range"]["contrasts"]
            signal.append(-float(Q(data[src]["projected"]["mean"])))
            ordinary.append(-float(Q(data[donor]["projected"]["mean"])))
        # Report the inferred *signed* enclosed profile at each measured
        # radius. A common Newtonian ordinary source would give constant r³K.
        raw = radii**3 * np.array(signal)
        ref = radii**3 * np.array(ordinary)
        tests = {}
        for j in np.flatnonzero(ref > 0):
            for tracer, bank in maps.items():
                for cap, panel in bank.items():
                    counts = {}
                    invalid = 0
                    for name, stock in panel["closed"].items():
                        tab = inputs[name]
                        rr, v, e, b = [
                            np.asarray(tab[k])
                            for k in (
                                "radius_kpc",
                                "observed_kms",
                                "error_kms",
                                "baryons_kms",
                            )
                        ]
                        scale = rr / meta[name]["disk_scale_kpc"]
                        good = (scale >= radii[0]) & (scale <= radii[-1])
                        # One hop per measured stellar disk scale length;
                        # interpolate the measured tidal field, not the orbit.
                        price = (
                            single_mass
                            * scale**3
                            * np.interp(scale, radii, signal)
                            / ref[j]
                        )
                        v2 = (
                            b * b
                            + G * np.asarray(stock) * np.asarray(price, float) / rr
                        )
                        ok = (
                            good
                            & (np.maximum(0, v - 2 * e) ** 2 <= v2)
                            & (v2 <= (v + 2 * e) ** 2)
                        )
                        invalid += int((good & (v2 < 0)).sum())
                        counts[name] = [int(ok.sum()), int(good.sum()), len(rr)]
                    tests[f"{int(radii[j])}/{tracer}/{cap}"] = dict(
                        within_two=sum(x[0] for x in counts.values()),
                        in_domain=sum(x[1] for x in counts.values()),
                        complete_in_domain=sum(
                            0 < n == q for q, n, N in counts.values()
                        ),
                        complete_catalogue_curves=sum(
                            q == N for q, n, N in counts.values()
                        ),
                        negative_speed_squared=invalid,
                        counts=counts,
                    )
        measured[src] = dict(
            radii=radii.tolist(),
            signed_enclosed_profile=raw.tolist(),
            ordinary_enclosed_profile=ref.tolist(),
            forecasts=tests,
        )
    return measured


def annular_transfer(radius, hop, knots, shape, nr, nt):
    # Each independent imaging annulus is uniform in area and azimuth.
    # delta is dimensionless distance from an off-center source element.
    leg, weights = np.polynomial.legendre.leggauss(nr)
    lo = np.r_[0, radius[:-1]]
    R = np.sqrt(lo[:, None] ** 2 + (radius**2 - lo**2)[:, None] * (leg + 1) / 2)[
        :, :, None
    ]
    cos = np.cos(2 * np.pi * (np.arange(nt) + 0.5) / nt)
    mat = []
    for x in radius:
        radial = x - R * cos
        delta = np.sqrt(np.maximum(0, x * x + R**2 - 2 * x * R * cos)) / hop
        kernel = np.interp(delta, knots, shape, left=0, right=0)
        mat.append((kernel * radial).mean(axis=2) @ (weights / 2))
    return np.asarray(mat)


def distributed_disks(research, maps, inputs, meta, single_mass):
    G = budget.inventory.rp.galaxy.transfer.budget.local.profiles.galaxies.G
    knots = np.array(list(map(float, research["shells"])))
    out = {}
    for source, donor in [("source", "single"), ("positive", "reference_plus_single")]:

        def profile(which):
            return np.array(
                [
                    -float(
                        Q(
                            v["comparisons"][source]["distant_equal_range"][
                                "contrasts"
                            ][which]["projected"]["mean"]
                        )
                    )
                    for v in research["shells"].values()
                ]
            )

        shape = profile(source)
        refs = knots**3 * profile(donor)
        assert shape[0] == shape[-1] == 0
        menu = {}
        free = {}
        for n, tab in inputs.items():
            r, v, e, b = [
                np.asarray(tab[k])
                for k in ("radius_kpc", "observed_kms", "error_kms", "baryons_kms")
            ]
            scale = meta[n]["disk_scale_kpc"]
            coarse, fine = [
                annular_transfer(r, scale, knots, shape, nr, na)
                for nr, na in [(12, 96), (24, 192)]
            ]
            for tracer, caps in maps.items():
                for cap, panel in caps.items():
                    mass = np.diff(np.r_[0, panel["closed"][n]])
                    assert np.all(mass >= 0)
                    conv = fine @ mass
                    bands = budget.passage_images.intervals(
                        G * r * conv / scale**3,
                        np.maximum(0, v - 2 * e) ** 2 - b * b,
                        (v + 2 * e) ** 2 - b * b,
                    )
                    free.setdefault(f"{tracer}/{cap}", {})[n] = [
                        (max(0.0, lo), hi) for lo, hi in bands
                    ]
                    change = abs(conv - coarse @ mass)
                    for j in np.flatnonzero(refs > 0):
                        price = G * single_mass * r / (scale**3 * refs[j])
                        v2 = b * b + price * conv
                        lower = np.maximum(0, v - 2 * e) ** 2
                        upper = (v + 2 * e) ** 2
                        ok = (v2 >= lower) & (v2 <= upper)
                        gross = price * change
                        tag = f"{int(knots[j])}/{tracer}/{cap}"
                        stats = menu.setdefault(
                            tag,
                            dict(
                                within_two=0,
                                whole_curves=0,
                                negative_speed_squared=0,
                                changed_classification_on_refinement=0,
                                galaxies={},
                            ),
                        )
                        stats["within_two"] += int(ok.sum())
                        stats["whole_curves"] += int(ok.all())
                        stats["negative_speed_squared"] += int((v2 < 0).sum())
                        vcoarse = b * b + price * (coarse @ mass)
                        stats["changed_classification_on_refinement"] += int(
                            np.count_nonzero(
                                ok != ((vcoarse >= lower) & (vcoarse <= upper))
                            )
                        )
                        stats["galaxies"][n] = dict(
                            within_two=int(ok.sum()),
                            rings=len(r),
                            max_refinement_kms2=float(gross.max()),
                        )
                        if j == 5 and tracer == "mapped_baryons" and cap == "1":
                            stats["galaxies"][n]["predicted_signed_v2_kms2"] = (
                                v2.tolist()
                            )
        shape_only = {}
        for tag, disks in free.items():
            all_bands = [band for rows in disks.values() for band in rows]
            rates = {0.0} | {
                u for a, b in all_bands if a <= b for u in (a, b) if np.isfinite(u)
            }
            best = max(sum(a <= u <= b for a, b in all_bands) for u in rates)
            whole = sum(
                max(a for a, b in rows) <= min(b for a, b in rows)
                for rows in disks.values()
            )
            shape_only[tag] = dict(
                rings_one_common_positive_scale=int(best),
                whole_curves_with_individual_scales=int(whole),
            )
        out[source] = dict(
            shape_only=shape_only,
            shells=knots.tolist(),
            native_inward_tide=shape.tolist(),
            ordinary_r_cubed_times_tide=refs.tolist(),
            forecasts=menu,
        )
    return out


def run():
    lab = budget.LAB / "tidal_release.json"
    motion = budget.LAB / "moving_records.json"
    field = HERE / "closed_inventory_maps.json"
    native = json.loads(lab.read_bytes())
    maps = json.loads(field.read_bytes())
    inputs, _, meta, _ = budget.inventory.rp.galaxy.transfer.budget.load_inputs()
    exterior = budget.LAB / "exterior_tide.json"
    energy = json.loads(motion.read_bytes())["currents"]["sharp_ordinary_standards"]
    mass = [Q(energy[k]["energy"]) / Q(energy["2"]["energy"]) for k in ("1", "2")]
    results = {}
    for order, p in native["schedules"].items():
        common = p["common_ruler"]
        for signal, donors in [
            ("source", ("single", "ordinary")),
            ("positive", ("reference_plus_single", "reference_plus_double")),
        ]:
            for clock in ("equal_sweep", "word_clock"):
                for i, probe in enumerate(("single", "double")):
                    source = Q(common[signal][clock][i])
                    for donor, m in zip(donors, mass, strict=True):
                        ordinary = Q(common[donor][clock][i])
                        tag = f"{order}/{signal}/{clock}/{probe}/{donor}"
                        row = dict(
                            source_chord=str(source), ordinary_chord=str(ordinary)
                        )
                        if ordinary < 0:
                            scale = source / ordinary * m
                            row.update(
                                extra_double_masses_per_patch=str(scale),
                                maps=evaluate_maps(scale, maps, inputs),
                            )
                        else:
                            row["calibration"] = (
                                "zero ordinary response"
                                if ordinary == 0
                                else "outward ordinary response"
                            )
                        results[tag] = row
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                lab,
                exterior,
                motion,
                field,
                HERE / "source_profile.json",
                rotation.galaxies.sparc.DATA / "table1.dat",
                rotation.galaxies.sparc.DATA / "table2.dat",
            )
        },
        error_band_test=precision_test(native),
        imaged_rings=sum(len(v["radius_kpc"]) for v in inputs.values()),
        single_to_double_assumed_mass_from_measured_energy=str(mass[0]),
        conditional_maps=results,
        distributed_disks=distributed_disks(
            json.loads(exterior.read_bytes()), maps, inputs, meta, float(mass[0])
        ),
        exterior_profile=exterior_profiles(
            json.loads(exterior.read_bytes()), maps, inputs, meta, float(mass[0])
        ),
        arbitrary_annular_response=radial_allowance(
            maps,
            inputs,
            [
                float(Q(x["extra_double_masses_per_patch"]))
                for x in results.values()
                if "maps" in x
            ],
        ),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    raw = json.dumps(result, indent=2, allow_nan=False) + "\n"
    path = HERE / "tidal_rotation.json"
    if args.write:
        path.write_text(raw)
    else:
        assert path.read_text() == raw
    print(result["error_band_test"]["example"])
