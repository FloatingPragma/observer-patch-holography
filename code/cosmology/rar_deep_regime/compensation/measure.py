"""Finite closed-redistribution mass and volume budgets for isolated galaxy lenses."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
from urllib.request import urlopen

import numpy as np

from cosmology.record_abundance import clustering, stationary

HERE = Path(__file__).parent
URL = "https://content.cld.iop.org/journals/2041-8205/969/1/L3/revision1/apjlad54b0t1_mrt.txt"
MRT_HASH = "2d2f446b70fbc74d8fa6fe66b51f69ac7615dcec1214a4f75bb62f900c947c71"
GAVAZZI = (
    "astro-ph/0701589v2",
    "slacs_wl.tex",
    "d01f4d18efdf71d2284486acf857957e74dadb284d54c3ebb244d440dc7e9b00",
)
G = 4.300917270e-9  # Mpc (km/s)^2 / solar mass
MSUN_KG = 1.98847e30
C_M_S = 299792458
OMEGA_B_H2 = 0.02242
OMEGA_C_H2 = 0.11933
LOG_BARYON_UPPER = np.array([10.5, 10.8, 11.1, 11.5])
RHO100 = 3 * 100**2 / (8 * np.pi * G)  # solar masses per Mpc^3


def selected_rows(raw):
    lines = [
        line for line in raw.decode("ascii").splitlines() if line[:24].strip() == "All"
    ]
    assert len(lines) == 15
    radii = [float(line[25:32]) for line in lines]
    assert all(b > a for a, b in zip(radii, radii[1:]))
    enrolled = {}
    for label, cut in (("near", 150), ("outer", 300)):
        i = max(i for i, r in enumerate(radii) if r <= cut)
        components = [float(x) for x in lines[i][33:].split()]
        assert len(components) == 12
        enrolled[label] = {
            "radius_kpc": radii[i],
            "v_stat_syst_kms": np.reshape(components, (4, 3)).tolist(),
        }
    return enrolled


def acquire():
    with urlopen(URL, timeout=30) as f:
        raw = f.read()
    if sha256(raw).hexdigest() != MRT_HASH:
        raise ValueError("Published machine-readable table has changed")
    fields = selected_rows(raw)
    tex = clustering.source(GAVAZZI)
    rows = [
        list(map(float, re.findall(r"-?\d+(?:\.\d+)?", line)[:4]))
        for line in tex.splitlines()
        if re.match(r"^\s*(164\.5|287\.5)\s*&", line)
    ]
    assert len(rows) == 2 and all(len(row) == 4 for row in rows)
    return {
        "mistele_table1_sha256": MRT_HASH,
        "gavazzi_table2_tex_sha256": GAVAZZI[-1],
        "mistele_all": fields,
        "gavazzi_last_two_rows_h_units": rows,
    }


def donor_radius(aperture_mpc, delta_mass, max_density):
    """Minimum common outer radius for a donor with a bounded mass density."""
    return np.cbrt(aperture_mpc**3 + 3 * delta_mass / (4 * np.pi * max_density))


def projected_ball(radius, size, density):
    """Surface and mean interior surface density of a homogeneous sphere."""
    dz = np.sqrt(np.maximum(0, size**2 - radius**2))
    sigma = 2 * density * dz
    mean = 4 * density * (size**3 - dz**3) / (3 * radius**2)
    return sigma, mean


def compensated_shear_mass(radius, inner, outer, deposit):
    """pi R² DeltaSigma for a uniform positive core and depleted donor shell."""
    donor = 3 * deposit / (4 * np.pi * (outer**3 - inner**3))
    core = 3 * deposit / (4 * np.pi * inner**3)
    if radius >= outer:
        return 0.0
    sigma_core, mean_core = projected_ball(radius, inner, core + donor)
    sigma_donor, mean_donor = projected_ball(radius, outer, donor)
    return np.pi * radius**2 * (mean_core - sigma_core - mean_donor + sigma_donor)


def mass_budget(row):
    R = row["radius_kpc"] / 1000
    v, stat, syst = np.array(row["v_stat_syst_kms"]).T
    baryon_cap = 10**LOG_BARYON_UPPER
    cases = {}
    # Authors propagate uncertainties from acceleration linearly into V.
    # Work on V^2 = R*g; the full systematic is taken in the same adverse
    # direction. This is conservative relative to squaring V-k*s-d.
    for label, k, system, b in (
        ("nominal", 0, 0, 1),
        ("one_and_systematic", 1, 1, 1),
        ("two_systematic_double_baryons", 2, 1, 2),
        ("bonferroni_four_normal_95_double_baryons", 2.2414027276, 1, 2),
    ):
        v2 = np.maximum(0, v * v - 2 * v * (k * stat + system * syst))
        total = R * v2 / G
        cap = b * baryon_cap
        excess = np.maximum(0, total - cap)
        radii = {}
        for species, omega in (("ordinary", OMEGA_B_H2), ("cold", OMEGA_C_H2)):
            for z in (0, 0.5):
                radii[f"{species}_z{z:g}"] = donor_radius(
                    R, excess, omega * RHO100 * (1 + z) ** 3
                ).tolist()
        rho_b = OMEGA_B_H2 * RHO100 * 1.5**3
        cases[label] = {
            "total_enclosed_mass_msun": total.tolist(),
            "subtracted_baryon_cap_msun": cap.tolist(),
            "min_net_excess_and_exterior_debit_msun": excess.tolist(),
            "min_exterior_debit_over_allowed_galaxy_baryons": (excess / cap).tolist(),
            "min_rest_energy_equivalent_joule": (excess * MSUN_KG * C_M_S**2).tolist(),
            "minimum_common_compensation_radius_mpc": radii,
            "required_density_over_mean_baryons_z0.5_for_1Mpc_compensation": (
                excess / (4 * np.pi / 3 * (1 - R**3) * rho_b)
            ).tolist(),
        }
    return {"radius_kpc": row["radius_kpc"], "cases": cases}


def direct_shear(inputs):
    # Ordered radial bins: previous bin center conservatively bounds every
    # contributing pair radius in the last bin, without importing its weights.
    previous, last = inputs["gavazzi_last_two_rows_h_units"]
    h = 0.7
    lower_radius_pc = previous[0] / h * 1000
    data = HERE.parent / "data/slacs_weak_lenses.dat"
    lenses = np.loadtxt(data, dtype=str, skiprows=2)
    assert len(lenses) == 22
    efficiency, logmass, error = lenses[:, [4, 7, 8]].astype(float).T
    budget = 2 * 10 ** (logmass + error)
    caps = {}
    for name, mass in (
        ("catalog_mean_source_efficiency", budget),
        ("arbitrary_source_efficiencies_up_to_one", budget / efficiency),
    ):
        cap = max(mass) / (np.pi * lower_radius_pc**2 * h)
        caps[name] = {
            "compact_ceiling_h_msun_pc2": cap,
            "measured_minus_ceiling_over_error": (last[1] - cap) / last[3],
        }
    return {
        "largest_pair_radius_lower_limit_kpc": lower_radius_pc / 1000,
        "reported_outer_bin_center_kpc": last[0] / h,
        "tangential_cross_error_h_msun_pc2": last[1:],
        "mass_catalog_sha256": sha256(data.read_bytes()).hexdigest(),
        "max_ordinary_budget_msun": max(budget),
        "minimum_compact_mass_if_no_extended_contribution_two_error_msun": max(
            0, last[1] - 2 * last[3]
        )
        * np.pi
        * lower_radius_pc**2
        * h,
        "caps": caps,
    }


def run(inputs):
    apertures = {
        label: mass_budget(row) for label, row in inputs["mistele_all"].items()
    }
    anchor = apertures["outer"]
    setting = anchor["cases"]["two_systematic_double_baryons"]
    inner = anchor["radius_kpc"] / 1000
    deposit = setting["min_net_excess_and_exterior_debit_msun"][-1]
    outer = setting["minimum_common_compensation_radius_mpc"]["ordinary_z0.5"][-1]
    radii = [inner / 2, inner, 2 * inner, outer / 2, outer, 1.2 * outer]
    return {
        "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        "rho_ordinary_z0_msun_mpc3": OMEGA_B_H2 * RHO100,
        "rho_cold_z0_msun_mpc3": OMEGA_C_H2 * RHO100,
        "apertures": apertures,
        "saturating_signed_shell_control": {
            "inner_outer_mpc": [inner, outer],
            "deposit_msun": deposit,
            "radii_mpc": radii,
            "pi_R2_DeltaSigma_over_deposit": [
                compensated_shear_mass(r, inner, outer, deposit) / deposit
                for r in radii
            ],
        },
        "independent_direct_shear": direct_shear(inputs),
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    if args.write and not args.fetch:
        p.error("writing measurements requires verified primary tables")
    ipath = HERE / "inputs.json"
    opath = HERE / "results.json"
    if args.fetch:
        inputs = acquire()
        if args.write:
            with ipath.open("x") as handle:
                json.dump(inputs, handle, indent=2)
                handle.write("\n")
        else:
            stationary.check(json.loads(ipath.read_text()), inputs)
    out = run(json.loads(ipath.read_text()))
    if args.write:
        with opath.open("x") as handle:
            json.dump(out, handle, indent=2, allow_nan=False)
            handle.write("\n")
    else:
        stationary.check(json.loads(opath.read_text()), out)
    for label, aperture in out["apertures"].items():
        print(
            label,
            aperture["radius_kpc"],
            aperture["cases"]["two_systematic_double_baryons"][
                "min_net_excess_and_exterior_debit_msun"
            ],
        )
