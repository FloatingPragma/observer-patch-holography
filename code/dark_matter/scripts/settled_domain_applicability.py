#!/usr/bin/env python3
"""Evaluate the settled-domain applicability predicate on three regimes.

The predicate scopes the static settled response nu_OPH(x) to closed settled
cuts. Clause P-A: the activation variable is the enclosed-flux ratio of the
cut, so divergence-free external fields carry zero opportunity intensity.
Clause P-B: the enclosed-source Jacobi rate exceeds the environmental tidal
Jacobi rate (Gamma_rec = Gamma_J is the declared conditional clock). Clause
P-C: the cut structure is older than ln(1/epsilon)/Gamma_J^int.

Gates checked here with declared numbers and zero fit freedom:
1. Cassini / Solar System: EXEMPT (zero interior anomaly quadrupole).
2. SPARC settled galaxies: RETAINED (predicate reduces to eq:ophrar).
3. Wide binaries: stance EMITTED (anomalous, full boost, no external-field
   suppression), plus the satellite tidal-clock rows.
"""

from __future__ import annotations

import argparse
import json
import math
from typing import Any

from d6_capacity_calculator import C, G, compute as compute_d6
from z6_shared_edge_reserve import P_PIXEL

M_SUN_KG = 1.988_47e30
AU_M = 1.495_978_707e11
KPC_M = 3.085_677_581_491_3673e19
PC_M = KPC_M / 1000.0
SEC_PER_YR = 365.25 * 24.0 * 3600.0
SEC_PER_GYR = SEC_PER_YR * 1.0e9

DEFAULT_N_SCR = 3.31e122
AGE_OF_UNIVERSE_GYR = 13.8

# Cassini quadrupole audit numbers (paper, Cassini Applicability Stress).
Q2_UNIVERSAL_Z6_S2 = 3.62018e-26
CASSINI_Q2_S2 = 1.6e-27
CASSINI_SIGMA_S2 = 1.8e-27
SUPPRESSION_FLOOR = 0.856

# Solar System / Galaxy inputs (public values).
M_SUN = 1.0
R_SATURN_AU = 9.5826
R0_KPC = 8.2
M_MW_ENC_R0_MSUN = 1.0e11
G_EXT_SOLAR_M_S2 = 1.8e-10

# SPARC representative row (paper Jacobi-clock table) and stored diagnostic.
OUTER_DISK_MASS_MSUN = 6.0e10
OUTER_DISK_RADIUS_KPC = 30.0
M31_MASS_MSUN = 1.5e12
MW_M31_SEPARATION_KPC = 780.0
STORED_RAR_ALLPOINT_RMS_DEX = 0.132834

# Wide-binary sample points.
WB_MASS_MSUN = 1.5
WB_SEPARATION_KAU = 10.0
WB_DEEP_GN_M_S2 = 1.0e-11
AQUAL_EFE_BOOST_CEILING = 1.43

# Satellite rows (representative public values).
SATELLITES = [
    {"name": "Crater II", "mb_msun": 3.0e5, "r_kpc": 1.07, "d_kpc": 117.0,
     "m_host_msun": 8.0e11},
    {"name": "Fornax", "mb_msun": 4.0e7, "r_kpc": 0.71, "d_kpc": 140.0,
     "m_host_msun": 9.0e11},
]


def declared_scales(n_scr: float) -> dict[str, float]:
    a0 = compute_d6(n_scr)["d12_benchmarks_not_theorem_outputs"][
        "a0_modular_anomaly_benchmark_m_s2"
    ]
    lam = 1.0 - P_PIXEL / 24.0
    return {
        "a0_oph_m_s2": a0,
        "lambda_collar": lam,
        "a_eff_m_s2": a0 / (lam * lam),
        "c_over_a0_gyr": C / a0 / SEC_PER_GYR,
        "c_over_a_eff_gyr": C / (a0 / (lam * lam)) / SEC_PER_GYR,
    }


def nu_oph(x: float, lam: float) -> float:
    m = lam * math.sqrt(x)
    if m > 700.0:
        return 1.0
    return 1.0 / (1.0 - math.exp(-m))


def nu_minus_one(x: float, lam: float) -> float:
    m = lam * math.sqrt(x)
    if m > 700.0:
        return 0.0
    e = math.exp(-m)
    return e / (1.0 - e)


def tau_j_s(mass_msun: float, radius_m: float) -> float:
    return math.sqrt(radius_m**3 / (G * mass_msun * M_SUN_KG))


def gate_cassini(scales: dict[str, float]) -> dict[str, Any]:
    a0, lam = scales["a0_oph_m_s2"], scales["lambda_collar"]
    r = R_SATURN_AU * AU_M
    g_sun = G * M_SUN * M_SUN_KG / r**2
    x = g_sun / a0
    m = lam * math.sqrt(x)
    interior_fraction = nu_minus_one(x, lam)  # underflows to exactly 0.0
    interior_bound = 1.0e-320 if interior_fraction == 0.0 else interior_fraction
    r_m = math.sqrt(G * M_SUN * M_SUN_KG / a0)
    jacobi_radius_m = (
        R0_KPC * KPC_M * (M_SUN / (3.0 * M_MW_ENC_R0_MSUN)) ** (1.0 / 3.0)
    )
    # P-A: the Galactic field is external to heliocentric cuts, so the solar
    # anomaly is spherically symmetric; interior force and tide vanish. The
    # residual is bounded by the anomalous part of the Galactic tide.
    gal_tide = G * M_MW_ENC_R0_MSUN * M_SUN_KG / (R0_KPC * KPC_M) ** 3
    q2_pred_bound = max(gal_tide, interior_bound * Q2_UNIVERSAL_Z6_S2)
    pull_sigma = abs(CASSINI_Q2_S2 - 0.0) / CASSINI_SIGMA_S2
    suppression = 1.0 - q2_pred_bound / Q2_UNIVERSAL_Z6_S2
    return {
        "regime": "solar_system_saturn",
        "g_sun_m_s2": g_sun,
        "x": x,
        "lambda_sqrt_x": m,
        "nu_minus_one_bound": interior_bound,
        "r_M_sun_au": r_m / AU_M,
        "r_M_sun_pc": r_m / PC_M,
        "sun_jacobi_radius_pc": jacobi_radius_m / PC_M,
        "transition_zone_inside_jacobi_radius": r_m < jacobi_radius_m,
        "galactic_tide_s2": gal_tide,
        "q2_predicted_bound_s2": q2_pred_bound,
        "q2_cassini_s2": CASSINI_Q2_S2,
        "q2_cassini_sigma_s2": CASSINI_SIGMA_S2,
        "pull_sigma": pull_sigma,
        "suppression_vs_universal": suppression,
        "suppression_floor": SUPPRESSION_FLOOR,
        "verdict_exempt": pull_sigma < 2.0 and suppression >= SUPPRESSION_FLOOR,
    }


def gate_sparc(scales: dict[str, float]) -> dict[str, Any]:
    a0, lam, a_eff = (
        scales["a0_oph_m_s2"],
        scales["lambda_collar"],
        scales["a_eff_m_s2"],
    )
    # P-A on spherical cuts is eq:ophrar. Limit checks, not refits:
    x_deep = 1.0e-8
    g_deep = x_deep * a0
    deep_rel_err = abs(
        nu_oph(x_deep, lam) * g_deep / math.sqrt(a_eff * g_deep) - 1.0
    )
    newtonian_dev = nu_minus_one(1.0e6, lam)
    tau_int = tau_j_s(OUTER_DISK_MASS_MSUN, OUTER_DISK_RADIUS_KPC * KPC_M)
    tau_ext = tau_j_s(M31_MASS_MSUN, MW_M31_SEPARATION_KPC * KPC_M)
    age_bound_gyr = math.log(100.0) * tau_int / SEC_PER_GYR
    return {
        "regime": "sparc_settled_galaxies",
        "deep_limit_rel_err_at_x_1e-8": deep_rel_err,
        "newtonian_limit_deviation_at_x_1e6": newtonian_dev,
        "tau_j_internal_outer_disk_gyr": tau_int / SEC_PER_GYR,
        "tau_j_external_mw_m31_gyr": tau_ext / SEC_PER_GYR,
        "clock_ratio_ext_over_int": tau_ext / tau_int,
        "age_clause_bound_gyr_at_eps_1e-2": age_bound_gyr,
        "stored_rar_allpoint_rms_dex": STORED_RAR_ALLPOINT_RMS_DEX,
        "verdict_retained": (
            deep_rel_err < 1.0e-3
            and newtonian_dev < 1.0e-6
            and tau_ext > tau_int
            and age_bound_gyr < AGE_OF_UNIVERSE_GYR
        ),
    }


def gate_wide_binaries(scales: dict[str, float]) -> dict[str, Any]:
    a0, lam = scales["a0_oph_m_s2"], scales["lambda_collar"]
    s = WB_SEPARATION_KAU * 1000.0 * AU_M
    g_n = G * WB_MASS_MSUN * M_SUN_KG / s**2
    boost_sample = math.sqrt(nu_oph(g_n / a0, lam))
    boost_deep = math.sqrt(nu_oph(WB_DEEP_GN_M_S2 / a0, lam))
    tau_int = tau_j_s(WB_MASS_MSUN, s)
    tau_ext = tau_j_s(M_MW_ENC_R0_MSUN, R0_KPC * KPC_M)
    satellites = []
    for sat in SATELLITES:
        t_int = tau_j_s(sat["mb_msun"], sat["r_kpc"] * KPC_M)
        t_ext = tau_j_s(sat["m_host_msun"], sat["d_kpc"] * KPC_M)
        satellites.append(
            {
                "name": sat["name"],
                "tau_j_internal_myr": t_int / SEC_PER_YR / 1.0e6,
                "tau_j_external_myr": t_ext / SEC_PER_YR / 1.0e6,
                "settled": t_int < t_ext,
                "emitted": "full boost" if t_int < t_ext else
                "slaved, suppressed below full boost",
            }
        )
    return {
        "regime": "wide_binaries_and_satellites",
        "g_n_sample_m_s2": g_n,
        "boost_at_sample": boost_sample,
        "boost_at_g_1e-11": boost_deep,
        "external_field_solar_m_s2": G_EXT_SOLAR_M_S2,
        "external_field_over_a_eff": G_EXT_SOLAR_M_S2 / scales["a_eff_m_s2"],
        "external_field_suppression": "none (P-A: zero enclosed flux)",
        "tau_j_internal_kyr": tau_int / SEC_PER_YR / 1.0e3,
        "tau_j_external_myr": tau_ext / SEC_PER_YR / 1.0e6,
        "clock_ratio_ext_over_int": tau_ext / tau_int,
        "aqual_efe_boost_ceiling": AQUAL_EFE_BOOST_CEILING,
        "boost_exceeds_efe_ceiling_in_deep_regime": (
            boost_deep > AQUAL_EFE_BOOST_CEILING
        ),
        "stance": "anomalous_full_boost_no_external_field_effect",
        "stance_sign": "+1 (Chae side; killed by a confirmed Newtonian verdict)",
        "satellites": satellites,
        "verdict_emitted": True,
    }


def route_one_record(scales: dict[str, float]) -> dict[str, Any]:
    cycles = AGE_OF_UNIVERSE_GYR / scales["c_over_a_eff_gyr"]
    return {
        "route": "repair saturation at the capacity rate",
        "c_over_a0_gyr": scales["c_over_a0_gyr"],
        "c_over_a_eff_gyr": scales["c_over_a_eff_gyr"],
        "max_repair_cycles_since_big_bang": cycles,
        "discriminating": cycles >= 1.0,
        "note": (
            "a0 = (15/8pi^2) c^2/r_dS is declared (benchmark grade), so "
            "c/a0 = 5.264 t_Lambda is forced; a capacity-rate cycle count is "
            "below one for every system and separates nothing."
        ),
    }


def compute(args: argparse.Namespace) -> dict[str, Any]:
    scales = declared_scales(args.n_scr)
    cassini = gate_cassini(scales)
    sparc = gate_sparc(scales)
    binaries = gate_wide_binaries(scales)
    return {
        "status": {
            "category": "settled-domain applicability predicate, gate check",
            "paper_grade": False,
            "named_premises": [
                "PR1 cut support: activation variable is the enclosed-flux "
                "ratio on closed settled cuts",
                "PR2 flux-recovery closure on settled cuts",
                "PR3 conditional clock Gamma_rec = Gamma_J",
                "PR4 record age >= ln(1/eps)/Gamma_J^int",
            ],
            "fit_freedom": "none; all scales declared, all masses public",
        },
        "declared_scales": scales,
        "route_one_record": route_one_record(scales),
        "gates": {
            "cassini": cassini,
            "sparc": sparc,
            "wide_binaries": binaries,
        },
        "all_gates_pass": bool(
            cassini["verdict_exempt"]
            and sparc["verdict_retained"]
            and binaries["verdict_emitted"]
        ),
    }


def print_markdown(payload: dict[str, Any]) -> None:
    scales = payload["declared_scales"]
    cas = payload["gates"]["cassini"]
    spc = payload["gates"]["sparc"]
    wb = payload["gates"]["wide_binaries"]
    print("# Settled-Domain Applicability Predicate: Gate Check")
    print()
    print(
        f"Declared scales: `a0 = {scales['a0_oph_m_s2']:.9e} m/s^2`, "
        f"`lambda_collar = {scales['lambda_collar']:.9f}`, "
        f"`a_eff = {scales['a_eff_m_s2']:.9e} m/s^2`."
    )
    print()
    print("| Gate | Key numbers | Verdict |")
    print("| --- | --- | --- |")
    print(
        f"| Cassini | lambda sqrt(x) = {cas['lambda_sqrt_x']:.1f}, "
        f"pull {cas['pull_sigma']:.2f} sigma, suppression "
        f"{cas['suppression_vs_universal']:.6f} | "
        f"{'EXEMPT' if cas['verdict_exempt'] else 'FAIL'} |"
    )
    print(
        f"| SPARC | clock ratio {spc['clock_ratio_ext_over_int']:.1f}, "
        f"RAR RMS {spc['stored_rar_allpoint_rms_dex']} dex | "
        f"{'RETAINED' if spc['verdict_retained'] else 'FAIL'} |"
    )
    print(
        f"| Wide binaries | boost {wb['boost_at_sample']:.3f} at "
        f"{wb['g_n_sample_m_s2']:.3e} m/s^2, {wb['boost_at_g_1e-11']:.3f} at "
        f"1e-11 | EMITTED ({wb['stance']}) |"
    )
    print()
    for sat in wb["satellites"]:
        print(
            f"- {sat['name']}: tau_int {sat['tau_j_internal_myr']:.0f} Myr, "
            f"tau_ext {sat['tau_j_external_myr']:.0f} Myr, {sat['emitted']}."
        )
    print()
    print(f"All gates pass: {payload['all_gates_pass']}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-scr", type=float, default=DEFAULT_N_SCR)
    parser.add_argument("--json", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    payload = compute(args)
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print_markdown(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
