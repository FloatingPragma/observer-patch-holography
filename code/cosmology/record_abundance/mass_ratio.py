"""Mass budgets for local mismatching-seam sources on any three-regular board."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import json
from math import pi, sqrt
import re

from . import clustering
from . import stationary as stock

HERE = stock.measure.HERE
RECEIPT = HERE / "mass_ratio.json"
BEGUM = (
    "0711.1588v1",
    "N3741_final_resubmit.tex",
    "8343c2215f5bf582df8fb17ed6b3a33b606da8452352618811cf630cfb91a5fe",
)
EARLIER = (
    "astro-ph/0502307v1",
    "NGC3741.tex",
    "c5617f8592a497fc3a8db7d8eaa33059056f85c76da2fc4d7e7787d784e352ce",
)
PLANCK18 = (
    "1807.06209v4",
    "ms.tex",
    "fe0833d2d0fda17427a709c9559dc7b3a988b955f47bb7de7a02c900ffdc21cf",
)
PLANCK13 = (
    "1303.5076v3",
    "abstract.tex",
    "8a2b1b01cd4179b416b426adc5e9628dac361c9015d825a9018e68f4350abacb",
)
SOURCES = (BEGUM, EARLIER, PLANCK18, PLANCK13, clustering.GO)
# Total r500 fractions and quoted errors in the authors' Planck reanalysis.
CLUSTERS = {
    "0122": (0.120, 0.013),
    "1651": (0.151, 0.013),
    "2401": (0.121, 0.014),
    "2721": (0.150, 0.021),
    "2811": (0.145, 0.012),
    "2955": (0.103, 0.011),
    "2984": (0.159, 0.016),
    "3112": (0.163, 0.010),
    "3693": (0.141, 0.014),
    "4010": (0.150, 0.012),
    "S0084": (0.118, 0.011),
    "S0296": (0.101, 0.014),
}
GAS, STAR, DYN = 1.69e8, 1.38e7, 4.03e9
DISTANCE, DISTANCE_ERROR = 3.03, 0.33


def verify_sources():
    with ThreadPoolExecutor(max_workers=len(SOURCES)) as pool:
        dwarf, earlier, recent, first, groups = pool.map(clustering.source, SOURCES)
    table = dwarf.split(r"\caption{Derived parameters for NGC~3741}", 1)[1].split(
        r"\end{tabular}", 1
    )[0]
    rows = re.findall(r"^(.+?)&\s*([\d.]+)\s*\$\\times~10\^(\d)\$~M", table, re.M)
    masses = [(label, float(f"{coeff}e{power}")) for label, coeff, power in rows]
    assert [mass for _, mass in masses] == [GAS, STAR, DYN]
    assert all(
        k in label
        for (label, _), k in zip(
            masses, ("1.3M_{HI}", "Stellar", "dynamical"), strict=True
        )
    )
    assert "a distance of 3.03 Mpc" in dwarf
    assert r"3.0\pm0.4~\times10^9" in earlier and r"1.6\pm0.4 \times{10}^{8}" in earlier
    assert r"dynamic mass to the baryonic mass ($\sim 8$)" in earlier
    for species, last in (("b", r"0.02242\pm 0.00014"), ("c", r"0.11933\pm 0.00091")):
        assert any(
            line.startswith(r"\Omega_{\mathrm{" + species + r"}} h^2&")
            and line.endswith(last + r"\cr")
            for line in recent.splitlines()
        )
    assert r"\Omega_{\rm b}h^2 = 0.02205 \pm 0.00028" in first
    assert r"\Omega_{\rm c}h^2 = 0.1199 \pm 0.0027" in first
    table = groups.split(
        r"\tablecaption{Derived Mass Fractions for Planck Cosmology", 1
    )[1].split(r"\enddata", 1)[0]
    rows = re.findall(
        r"^Abell (\w+)\s*&[^&\n]+&[^&\n]+&\s*\$([\d.]+)\\pm([\d.]+)\$", table, re.M
    )
    assert {name: (float(f), float(e)) for name, f, e in rows} == CLUSTERS
    assert table.count("---") == 6  # The remaining three groups only have gas masses.
    assert r"f_{bary}=0.144\pm0.005" in groups and r"f_{bary}=0.136\pm0.005" in groups


def dwarf_budget(total_cap, distance=DISTANCE):
    scale = distance / DISTANCE
    dyn = DYN * scale
    known = (GAS + STAR) * scale**2
    missing = max(0.0, dyn / total_cap - known)
    dark_cap = total_cap - 1
    return {
        "distance_mpc": distance,
        "dynamical_msun": dyn,
        "catalogued_stars_and_gas_msun": known,
        "required_total_baryons_msun": dyn / total_cap,
        "minimum_extra_interior_baryons_msun": missing,
        "if_all_ionized_extra_hydrogen_msun": missing / 1.3,
        "nonbaryonic_to_catalogued_baryons": dyn / known - 1,
        "maximum_fraction_of_inferred_dark_mass_from_known_baryons": dark_cap
        * known
        / (dyn - known),
        "minimum_one_link_exterior_baryons_if_no_hidden_interior_msun": 2
        * missing
        * total_cap
        / dark_cap,
        "maximum_velocity_fraction_at_fixed_radius_with_known_baryons": sqrt(
            total_cap * known / dyn
        ),
    }


def aperture(edges, occupied, region):
    """Measure source-owned and symmetric cuts and their signed boundary correction."""
    owned = symmetric_twice = external_active = internal_active = boundary = 0
    outside_shell = set()
    for u, v in edges:
        inside_u, inside_v = u in region, v in region
        if inside_u != inside_v:
            boundary += 1
            outside_shell.add(v if inside_u else u)
        if (u in occupied) == (v in occupied):
            continue
        owner_inside = inside_u if u in occupied else inside_v
        owned += owner_inside
        symmetric_twice += int(inside_u) + int(inside_v)
        if inside_u != inside_v:
            internal_active += owner_inside
            external_active += not owner_inside
    return {
        "ordinary": len(occupied & region),
        "owned": owned,
        "symmetric_twice": symmetric_twice,
        "active_outside_crossings": external_active,
        "active_inside_crossings": internal_active,
        "boundary_edges": boundary,
        "external_occupied_neighbors": len(occupied & outside_shell),
    }


def warm_gas(result):
    """Follow-up for case-B, 10,000-K gas in a sphere with the adopted radius."""
    radius_cm = result["outer_radius_upper_scale_kpc"] * 3.08567758e21
    volume = 4 * pi * radius_cm**3 / 3
    protons = (
        result["ngc3741"]["if_all_ionized_extra_hydrogen_msun"]
        * 1.98847e33
        / 1.6726219e-24
    )
    minimum_luminosity = 3.54e-25 * protons**2 / volume
    # L/(4*pi*d^2), averaged over pi*(R/d)^2 steradians.
    mean_flux = minimum_luminosity / (4 * pi**2 * radius_cm**2 * 206265**2)
    reported = 4.3e38
    return {
        "minimum_luminosity_erg_per_s": minimum_luminosity,
        "reported_galaxy_luminosity_erg_per_s": reported,
        "minimum_mean_surface_brightness_erg_s_cm2_arcsec2": mean_flux,
        "survey_typical_surface_brightness_limit_erg_s_cm2_arcsec2": 4e-18,
        "warm_gas_upper_msun_if_reported_light_is_complete": 1.3
        * sqrt(reported * volume / 3.54e-25)
        * 1.6726219e-24
        / 1.98847e33,
    }


def run():
    p = stock.screen.DENSITY
    gain = 1 / (1 - p)
    cosmic_ratio = 0.11933 / 0.02242
    cap = cosmic_ratio * gain
    # Reference at z ~ 0. The primary puts the measured endpoint inside 7.6 kpc.
    radius_kpc = 7.6
    reference_baryons = 0.02242 * 2.77536627e11 * 4 * pi / 3 * (radius_kpc / 1000) ** 3
    dwarf = dwarf_budget(1 + cap)
    # Use the cold/baryon ratio of the Planck 2013 release for the cluster appendix.
    old_ratio = 0.1199 / 0.02205
    floor = 1 / (1 + old_ratio * gain)
    groups = {
        name: {
            "baryon_fraction": f,
            "quoted_error": err,
            "minimum_extra_baryons_over_total_at_central_mass": max(0.0, floor - f),
        }
        for name, (f, err) in CLUSTERS.items()
    }
    return {
        "protocol_sha256": sha256((HERE / "MASS_RATIO.md").read_bytes()).hexdigest(),
        "primary_tex_sha256": {s[0]: s[2] for s in SOURCES},
        "p": p,
        "stationary_cuts_per_ordinary": 3 * (1 - p),
        "relative_enhancement_cap": gain,
        "external_2018_cold_to_baryons": cosmic_ratio,
        "dark_to_baryons_cap": cap,
        "minimum_baryon_fraction": 1 / (1 + cap),
        "ngc3741": dwarf,
        "distance_sensitivity": [
            dwarf_budget(1 + cap, d)
            for d in (DISTANCE - DISTANCE_ERROR, DISTANCE + DISTANCE_ERROR)
        ],
        "required_distance_if_no_unseen_mass_mpc": DISTANCE
        * DYN
        / ((1 + cap) * (GAS + STAR)),
        "favorable_2018_marginal_error_corner_missing_baryons_msun": dwarf_budget(
            1 + gain * (0.11933 + 0.00091) / (0.02242 - 0.00014)
        )["minimum_extra_interior_baryons_msun"],
        "outer_radius_upper_scale_kpc": radius_kpc,
        "mean_missing_hplus_column_lower_if_all_ionized_cm_minus2": dwarf[
            "if_all_ionized_extra_hydrogen_msun"
        ]
        * 1.98847e33
        / (pi * (radius_kpc * 3.08567758e21) ** 2 * 1.6726219e-24),
        "uniform_reference_baryons_msun_inside_upper_scale": reference_baryons,
        "maximum_background_subtraction_correction_msun": (cap - cosmic_ratio)
        * reference_baryons,
        "planck2013_clusters": {
            "external_cold_to_baryons": old_ratio,
            "minimum_baryon_fraction": floor,
            "systems": groups,
            "below_at_center": [name for name, (f, _) in CLUSTERS.items() if f < floor],
            "below_at_upper_two_quoted_errors": [
                name for name, (f, e) in CLUSTERS.items() if f + 2 * e < floor
            ],
            "published_high_mass_weighted_mean": 0.144,
            "published_error_on_mean": 0.005,
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fetch",
        action="store_true",
        help="verify mass and abundance inputs against versioned primary sources",
    )
    parser.add_argument("--write", action="store_true")
    parser.add_argument(
        "--warm-gas",
        action="store_true",
        help="also print the conditional recombination-light budget",
    )
    args = parser.parse_args()
    if args.fetch:
        verify_sources()
    measured = run()
    if args.write:
        if not args.fetch:
            parser.error("--write requires --fetch")
        with RECEIPT.open("x") as out:
            json.dump(measured, out, indent=2, allow_nan=False)
            out.write("\n")
    else:
        stock.check(json.loads(RECEIPT.read_text()), measured)
    print(json.dumps(measured["ngc3741"], indent=2))
    if args.warm_gas:
        print(json.dumps(warm_gas(measured), indent=2))
