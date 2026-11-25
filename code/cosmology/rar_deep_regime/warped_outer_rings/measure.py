"""Resolved exterior H I rings of NGC 2685 against the held SPARC calibration."""

import argparse
import json
from math import ceil, floor, pi, radians, sin, sqrt
from pathlib import Path

from ..udg_transfer.measure import GM_SUN, KPC

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
G_KPC = GM_SUN / KPC / 1e6
ARCSEC = pi / 648000
# Convert 10^19 neutral H atoms/cm² to neutral-hydrogen solar masses/pc².
SURFACE = 1e19 * 1.6735575e-24 * (KPC / 10) ** 2 / 1.98847e33


def read(path):
    return json.loads(path.read_text())


def gas_inside(profile, theta, distance, errors):
    """Absolute H I mass with helium inside 0.9 theta in the tilted ring planes.

    Radii are angular, and the local face-on column is interpolated linearly.
    Integrate theta times linear column exactly between knots. The conservative
    envelope lowers every knot by its supplied error times a common multiplier.
    """
    limit = 0.9 * theta
    if limit > profile[-1]["radius_arcsec"]:
        raise ValueError("unobserved gas radius")
    moment = 0.0
    for left, right in zip(profile, profile[1:]):
        x, y = left["radius_arcsec"], right["radius_arcsec"]
        if x >= limit:
            break
        end = min(y, limit)
        a, b = (
            max(0, r["NHI_1e19_cm2"] - errors * r["error_NHI"]) for r in (left, right)
        )
        slope = (b - a) / (y - x)
        u = end - x
        moment += a * (x * u + u * u / 2) + slope * (x * u * u / 2 + u**3 / 3)
    return 1.4 * 2 * pi * (distance * 1e6 * ARCSEC) ** 2 * SURFACE * moment


def enclosed(profile, ring, host):
    theta = ring["radius_arcsec"]
    star = 0.6 * host["stellar_unit_msun"]
    star_in = 0.6 * (
        host["stellar_inside_unit_msun"] - 0.01 * host["stellar_unit_msun"]
    )
    max_mass = star + 10**0.1 * host["helium_hi_msun"]
    table = {}
    for k in (0, 1, 2):
        gas = gas_inside(profile, theta, host["distance_mpc"], k)
        needed = (0.90 * max_mass - star_in) / gas
        table[str(k)] = {
            "cylindrical_gas_msun": gas,
            "fraction_retaining_95pct": (star_in + 0.95 * gas) / max_mass,
            "minimum_retention_for_90pct": needed,
            "max_rms_height_for_90pct_arcsec": theta * sqrt(0.19) * sqrt(1 - needed)
            if needed < 1
            else None,
            "unlocated_extra_mass_allowance_msun": (star_in + 0.95 * gas) / 0.9
            - max_mass,
        }
    return table


def speed(gm, r, a0):
    excess = sqrt(gm * a0 * KPC / 1e6)
    return sqrt(gm / r + excess)


def coherent_phase(rows, gm, spacing, a0):
    """All rings share a single phase on the registered exterior cut lattice."""
    units = [r["radius_kpc"] / spacing for r in rows]
    breaks = sorted({0.0, 1.0, *[1.0 - (u % 1) for u in units]})
    plateau = sqrt(gm * a0 * KPC / 1e6)
    scored = []
    for a, b in zip(breaks, breaks[1:]):
        phase = (a + b) / 2
        maxima = max(
            abs(
                sqrt(gm / r["radius_kpc"] + plateau * floor(u + phase) / u)
                - r["speed_kms"]
            )
            / r["error_speed"]
            for r, u in zip(rows, units, strict=True)
        )
        scored.append((maxima, b - a))
    return dict(
        minimax_error_multiplier_range=[
            min(x for x, w in scored),
            max(x for x, w in scored),
        ],
        uniform_phase_fraction_all_within={
            str(k): sum(w for x, w in scored if x <= k) for k in (1, 2)
        },
    )


def edge_interval(rows, mass, a0, k):
    """Intersect all absolute-speed boxes for a sharp common physical edge."""
    gm = G_KPC * mass
    plateau = sqrt(gm * a0 * KPC / 1e6)
    low = 0.0
    high = None
    overshoot = False
    negative = False
    for row in rows:
        r, v, e = (row[t] for t in ("radius_kpc", "speed_kms", "error_speed"))
        zlo = (max(0, v - k * e) ** 2 - gm / r) / plateau
        zhi = ((v + k * e) ** 2 - gm / r) / plateau
        overshoot |= zlo > 1
        negative |= zhi < 0
        low = max(low, r * zlo)
        if zhi < 1:
            high = r * zhi if high is None else min(high, r * zhi)
    return dict(
        feasible=not (overshoot or negative or (high is not None and high < low)),
        lower_kpc=low,
        upper_kpc=high,
        above_unlimited_plateau=overshoot,
    )


def dwarf_comparison(ring, host, cases):
    # Universal f = M_excess(<r) / M_b, x=r/sqrt(M_b); dM must be positive.
    # One extreme host in each box minimizes required f. x does not depend
    # on adopted distance: radii scale with d, photometric masses with d².
    out = {}
    for label, k, stellar, gas_scale, tilt in (
        ("nominal", 0, 0.8, 1.0, False),
        ("one_error", 1, 1.0, 10**0.1, False),
        ("two_errors", 2, 1.0, 10**0.1, False),
        ("two_errors_plus_tilt", 2, 1.0, 10**0.1, True),
    ):
        m0 = stellar * host["stellar_unit_msun"] + gas_scale * host["helium_hi_msun"]
        d = 1 + k * host["distance_error_mpc"] / host["distance_mpc"]
        v = max(0, ring["speed_kms"] - k * ring["error_speed"])
        if tilt:
            i = ring["inclination_deg"]
            v *= sin(radians(i)) / sin(
                radians(min(90, i + k * ring["error_inclination"]))
            )
        r = ring["radius_kpc"]
        f = r * v * v / (d * G_KPC * m0) - 1
        x = r / sqrt(m0)
        largest_x = r / sqrt(
            (0.8 if k == 0 else 0.6) * host["stellar_unit_msun"]
            + (1 if k == 0 else 10**-0.1) * host["helium_hi_msun"]
        )
        out[label] = dict(
            required_minimum_excess_ratio=f,
            host_x_kpc_per_sqrt_msun=x,
            largest_host_x_in_box=largest_x,
            against={},
        )
        for name, dwarf in cases.items():
            dx = dwarf["radius_kpc"] / sqrt(dwarf["mass_solar"])
            cap = dwarf["max_halo_to_baryonic_total_mass"]
            out[label]["against"][name] = dict(
                dwarf_x_kpc_per_sqrt_msun=dx,
                dwarf_max_excess_ratio=cap,
                radial_order_preserved=largest_x < dx,
                minimum_gap=f - cap,
            )
    return out


def run():
    catalog = read(ROOT / "outer_orbit_transfer" / "inputs.json")["galaxies"]
    host = next(r for r in catalog if r["name"] == "NGC2685")
    archive = read(ROOT / "cut_spacing" / "results.json")
    if any(
        k.replace(" ", "").upper() == "NGC2685"
        for k in archive["held_log_density_by_galaxy"]
    ):
        raise ValueError("test galaxy overlaps the training catalogue")
    spaced = read(ROOT / "cut_spacing" / "disk_domain.json")
    lens = read(ROOT / "receipts" / "lensing_profile.json")
    a0 = read(ROOT / "finite_halo" / "results.json")["sparc_a0_m_s2"]
    dwarf = read(ROOT / "finite_halo" / "transfer.json")["cases"]
    profile = read(HERE / "rings.json")["rows"]
    mass = 0.8 * host["stellar_unit_msun"] + host["helium_hi_msun"]
    gm = G_KPC * mass
    sample = []
    for ring in profile:
        theta = ring["radius_arcsec"]
        if theta < 320:
            continue  # stellar enclosure certified at 320 arcseconds
        r = host["radius_kpc"] * theta / 320
        inventory = enclosed(profile, ring, host)
        g = GM_SUN * mass / (KPC * r) ** 2
        primary = (
            inventory["1"]["fraction_retaining_95pct"] >= 0.9
            and r >= 5 * host["re_kpc"]
            and g <= 1.2e-11
        )
        if not primary:
            continue
        pred = speed(gm, r, a0)
        u = r / spaced["spacing_kpc"]
        space_a = sqrt(gm * spaced["a0_m_s2"] * KPC / 1e6)
        sample.append(
            dict(
                **ring,
                radius_kpc=r,
                gbar_m_s2=g,
                enclosure=inventory,
                forecast_kms=pred,
                residual_over_quoted_error=(ring["speed_kms"] - pred)
                / ring["error_speed"],
                finite_spacing_possible_kms=[
                    sqrt(gm / r + space_a * n / u) for n in (floor(u), ceil(u))
                ],
                required_excess_ratio=r * ring["speed_kms"] ** 2 / gm - 1,
            )
        )
    if not sample:
        raise ValueError("no exterior ring passed the photometric cuts")
    # The only common fit is labelled as a descriptive scale control. On this
    # photometric geometry v_i(D)=v_i(D0)sqrt(D/D0); require one D for all.
    scale_error = max(
        (x["speed_kms"] / x["forecast_kms"] - y["speed_kms"] / y["forecast_kms"])
        / (x["error_speed"] / x["forecast_kms"] + y["error_speed"] / y["forecast_kms"])
        for x in sample
        for y in sample
    )
    z = max(
        (r["speed_kms"] - scale_error * r["error_speed"]) / r["forecast_kms"]
        for r in sample
    )
    outer = sample[-1]
    checks = {}
    for n in (0, 1, 2):
        # Broad independent upper mass, SPARC normalization and distance;
        # this is a necessary lower bound, not a joint best-fit interval.
        m0 = host["stellar_unit_msun"] + 10**0.1 * host["helium_hi_msun"]
        gd = G_KPC * m0
        d = 1 + n * host["distance_error_mpc"] / host["distance_mpc"]
        h = sqrt(gd * lens["sparc_cluster_bootstrap_95pct_m_s2"][1] * KPC / 1e6)
        checks[str(n)] = max(
            0.0,
            max(
                (
                    s["radius_kpc"] * max(0, s["speed_kms"] - n * s["error_speed"]) ** 2
                    - d * gd
                )
                / h
                for s in sample
            ),
        )
    return dict(
        source_table="Jozsa et al. Table C.1, arXiv:0810.3025v1",
        host="NGC2685",
        distance_mpc=host["distance_mpc"],
        mass_msun=mass,
        reference_a0_m_s2=a0,
        finite_spacing_calibration=dict(
            spacing_kpc=spaced["spacing_kpc"], a0_m_s2=spaced["a0_m_s2"]
        ),
        prior_enclosure_lower=read(ROOT / "outer_orbit_transfer" / "results.json")[
            "enrollment"
        ]["NGC2685"]["enclosed_fraction_lower"],
        sample=sample,
        common_distance_control=dict(
            minimax_error_multiplier=scale_error,
            distance_mpc=z * z * host["distance_mpc"],
        ),
        shared_cut_phase=coherent_phase(
            sample, gm, spaced["spacing_kpc"], spaced["a0_m_s2"]
        ),
        sharp_edge_nominal={str(k): edge_interval(sample, mass, a0, k) for k in (1, 2)},
        broad_necessary_edge_lower_kpc_by_speed_distance_error=checks,
        cross_galaxy_at_outer_ring=dwarf_comparison(outer, host, dwarf),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    text = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if p.parse_args().write:
        (HERE / "results.json").write_text(text)
    else:
        print(text, end="")
