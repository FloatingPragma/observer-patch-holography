#!/usr/bin/env python3
"""Run the executable empirical diagnostics for the OPH dark sector."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import dark_cluster_boltzmann_prelikelihood
import dark_cmb_bao_growth_s8_likelihood
import dark_homogeneous_state_selection
import dark_parent_collar_grid
import dark_sector_measurement_ledger


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PARENT_GRID = ROOT / "outputs" / "dark_parent_collar_grid_diagnostic.json"
DEFAULT_OUT_JSON = ROOT / "outputs" / "dark_empirical_scorecard.json"
DEFAULT_OUT_MD = ROOT / "dark_empirical_implementation_status.md"


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def build_parent_grid(args: argparse.Namespace) -> dict[str, Any]:
    parent_args = argparse.Namespace(
        samples_json=args.samples_json,
        diagnostic_scalar_load=True,
        mu_eq=args.mu_eq,
        B_A=args.B_A,
        n_max=args.n_max,
        hold=args.hold,
        a_grid=args.a_grid,
        k_grid=args.k_grid,
        out=None,
        json=False,
    )
    parent = dark_parent_collar_grid.compute(parent_args)
    write_json(args.parent_grid_out, parent)
    return parent


def build_homogeneous_state(args: argparse.Namespace) -> dict[str, Any]:
    return dark_homogeneous_state_selection.compute(
        argparse.Namespace(
            n_scr=args.n_scr,
            H0_km_s_Mpc=args.H0,
            ombh2=args.ombh2_homogeneous,
            sum_mnu_eV=args.sum_mnu_eV,
            omega_r=args.omega_r,
            json=False,
        )
    )


def build_sparc(args: argparse.Namespace) -> dict[str, Any]:
    return dark_sector_measurement_ledger.compute(
        argparse.Namespace(
            n_scr=args.n_scr,
            rar_all=args.rar_all,
            rar_bins=args.rar_bins,
            mass_models=args.mass_models,
            galaxy_sample=args.galaxy_sample,
            fit_min_a0=args.fit_min_a0,
            fit_max_a0=args.fit_max_a0,
            upsilon_disk=args.upsilon_disk,
            upsilon_bulge=args.upsilon_bulge,
            disk_prior=args.disk_prior,
            bulge_prior=args.bulge_prior,
            ml_prior_sigma_dex=args.ml_prior_sigma_dex,
            intrinsic_scatter_km_s=args.intrinsic_scatter_km_s,
            disk_min=args.disk_min,
            disk_max=args.disk_max,
            bulge_min=args.bulge_min,
            bulge_max=args.bulge_max,
            coordinate_iterations=args.coordinate_iterations,
            max_quality=args.max_quality,
            min_inclination_deg=args.min_inclination_deg,
            profile_distance=args.profile_distance,
            distance_min_mpc=args.distance_min_mpc,
            distance_sigma_window=args.distance_sigma_window,
            profile_inclination=args.profile_inclination,
            inclination_min_deg=args.inclination_min_deg,
            inclination_max_deg=args.inclination_max_deg,
            inclination_sigma_window=args.inclination_sigma_window,
            profile_gas=args.profile_gas,
            gas_prior_sigma_dex=args.gas_prior_sigma_dex,
            gas_min=args.gas_min,
            gas_max=args.gas_max,
        )
    )


def build_cluster(
    args: argparse.Namespace,
    homogeneous: dict[str, Any],
) -> dict[str, Any]:
    fractions = homogeneous["density_fractions"]
    return dark_cluster_boltzmann_prelikelihood.compute(
        argparse.Namespace(
            mu_eq=args.mu_eq,
            n_max=args.n_max,
            hold=args.hold,
            tau_values_gyr=args.tau_values_gyr,
            separation_kpc=args.separation_kpc,
            time_since_passage_gyr=args.time_since_passage_gyr,
            observed_offset_kpc=args.observed_offset_kpc,
            offset_sigma_kpc=args.offset_sigma_kpc,
            omega_b=fractions["Omega_b"],
            omega_anomaly=fractions["Omega_A"],
            omega_lambda=fractions["Omega_Lambda_OPH"],
            omega_r=fractions["Omega_r"],
            h0_km_s_mpc=args.H0,
            redshifts=args.redshifts,
            parent_collar_json=str(args.parent_grid_out),
            json=False,
        )
    )


def build_cmb(args: argparse.Namespace) -> dict[str, Any]:
    try:
        return dark_cmb_bao_growth_s8_likelihood.compute(
            argparse.Namespace(
                parent_grid=str(args.parent_grid_out),
                H0=args.H0,
                ombh2=args.ombh2_cmb,
                sum_mnu_eV=args.sum_mnu_eV,
                tau=args.tau_reio,
                As=args.As,
                ns=args.ns,
                lmax=args.lmax,
                kmax=args.kmax,
                weak_lensing_s8=args.weak_lensing_s8,
                weak_lensing_s8_sigma=args.weak_lensing_s8_sigma,
                weak_lensing_source=args.weak_lensing_source,
                json=False,
            )
        )
    except SystemExit as exc:
        return {
            "status": {
                "category": "compressed dark-sector CMB BAO growth S8 likelihood",
                "ready": False,
                "reason": str(exc),
            },
            "rows": [],
            "total_chi2_diagonal": None,
        }


def scenario(payload: dict[str, Any], lane: str, name: str) -> dict[str, Any]:
    return payload[lane][name]


def render_markdown(payload: dict[str, Any]) -> str:
    sparc = payload["sparc"]
    cmb = payload["cmb_bao_growth"]
    cluster = payload["cluster_boltzmann"]
    parent = payload["parent_grid"]
    homogeneous = payload["homogeneous_state"]
    inputs = payload["inputs"]

    unit = sparc["rar"]["unit"]
    z6 = sparc["rar"]["z6_poisson"]
    best = sparc["rar"]["best"]
    systematic_unit = scenario(sparc, "systematic_likelihood", "unit")["stats"]
    systematic_z6 = scenario(sparc, "systematic_likelihood", "z6_poisson")["stats"]
    systematic_emp = scenario(sparc, "systematic_likelihood", "empirical")["stats"]
    cluster_rows = cluster["cluster"]["rows"]
    best_cluster = min(
        (row for row in cluster_rows if "offset_chi2" in row),
        key=lambda row: row["offset_chi2"],
        default=None,
    )

    lines = [
        "# Dark Empirical Implementation Claim Boundary",
        "",
        "This file is generated by `dark_empirical_scorecard.py`.",
        "It records executable diagnostics at pre-likelihood grade.",
        "",
        "## Inputs",
        "",
        "| Input | Value |",
        "| --- | ---: |",
        f"| N_scr | `{inputs['N_scr']:.9e}` |",
        f"| parent grid category | `{parent['status']['category']}` |",
        f"| parent paper grade | `{parent['status']['paper_grade']}` |",
        f"| rho_A/rho_b supplied to parent grid | `{parent['rho_A_over_rho_b']:.9f}` |",
        f"| homogeneous Omega_A | `{homogeneous['density_fractions']['Omega_A']:.9f}` |",
        f"| homogeneous Omega_Lambda | `{homogeneous['density_fractions']['Omega_Lambda_OPH']:.9f}` |",
        f"| repair gamma_rec | `{parent['repair_matrix']['gap']['gamma_rec']:.9f}` |",
        "",
        "## SPARC Galaxy Tests",
        "",
        "| Quantity | OPH unit | Z6/Poisson | Empirical or best comparison |",
        "| --- | ---: | ---: | ---: |",
        (
            f"| acceleration scale a0 m/s^2 | `{unit['a0_m_s2']:.9e}` | "
            f"`{z6['a0_eff_m_s2']:.9e}` | `{best['a0_m_s2']:.9e}` best same-function RAR |"
        ),
        (
            f"| RAR all-point RMS dex | `{unit['all_data']['rms_residual_dex']:.6f}` | "
            f"`{z6['rar_all_rms_dex']:.6f}` | `{best['all_data']['rms_residual_dex']:.6f}` |"
        ),
        (
            f"| RAR binned RMS dex | `{unit['binned_data']['rms_residual_dex_weighted_by_N']:.6f}` | "
            f"`{z6['rar_binned_weighted_rms_dex']:.6f}` | "
            f"`{best['binned_data']['rms_residual_dex_weighted_by_N']:.6f}` |"
        ),
        (
            f"| systematic SPARC chi2/pt | `{systematic_unit['chi2_per_point']:.6f}` | "
            f"`{systematic_z6['chi2_per_point']:.6f}` | "
            f"`{systematic_emp['chi2_per_point']:.6f}` empirical reference |"
        ),
        (
            f"| systematic SPARC RMS km/s | `{systematic_unit['rms_residual_km_s']:.6f}` | "
            f"`{systematic_z6['rms_residual_km_s']:.6f}` | "
            f"`{systematic_emp['rms_residual_km_s']:.6f}` empirical reference |"
        ),
        "",
        "Interpretation: the unit branch gives the best nuisance-profiled SPARC scaffold under the stated priors. The Z6/Poisson branch gives the closest acceleration normalization and the best binned RAR residual, while its coefficient remains conditional on the finite-thickness reserve theorem.",
        "",
        "## Compressed CMB, BAO, Growth, And S8",
        "",
    ]

    if not cmb["status"].get("ready", False):
        lines.extend(
            [
                f"Ready: `{cmb['status'].get('ready')}`",
                f"Reason: `{cmb['status'].get('reason', 'unspecified')}`",
                "",
            ]
        )
    else:
        lines.extend(
            [
                f"CAMB version: `{cmb['boltzmann']['camb_version']}`",
                "",
                "| Quantity | Prediction | Target | z | chi2 |",
                "| --- | ---: | ---: | ---: | ---: |",
            ]
        )
        for row in cmb["rows"]:
            lines.append(
                f"| {row['name']} | `{row['prediction']:.9g}` | "
                f"`{row['central']:.9g} +/- {row['sigma']:.3g}` | "
                f"`{row['z']:+.3f}` | `{row['chi2']:.3f}` |"
            )
        lines.extend(
            [
                "",
                (
                    f"Diagonal compressed chi2: `{cmb['total_chi2_diagonal']:.3f}` "
                    f"for `{cmb['dof_diagonal']}` rows."
                ),
                "",
            ]
        )

    lines.extend(
        [
            "## Cluster Timing Gate",
            "",
            (
                f"Example target: `{inputs['observed_offset_kpc']:.6g} +/- "
                f"{inputs['offset_sigma_kpc']:.6g} kpc` after "
                f"`{inputs['time_since_passage_gyr']:.6g} Gyr`, with an initial "
                f"`{inputs['separation_kpc']:.6g} kpc` anomaly separation."
            ),
            "",
            "| tau_rec Gyr | t_commit Myr | retained fraction | model offset kpc | z | chi2 |",
            "| ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in cluster_rows:
        lines.append(
            f"| `{row['tau_rec_Gyr']:.6g}` | `{1000.0 * row['t_commit_Gyr_from_K']:.6g}` | "
            f"`{row['retention_fraction']:.6g}` | `{row['model_offset_kpc']:.6g}` | "
            f"`{row.get('offset_z', 0.0):+.3f}` | `{row.get('offset_chi2', 0.0):.3f}` |"
        )
    if best_cluster is not None:
        lines.extend(
            [
                "",
                (
                    f"Smallest example chi2 in this grid: `tau_rec = "
                    f"{best_cluster['tau_rec_Gyr']:.6g} Gyr`, "
                    f"model offset `{best_cluster['model_offset_kpc']:.6g} kpc`."
                ),
            ]
        )

    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "- SPARC numbers use public SPARC tables. The systematic run profiles stellar mass-to-light ratios, distance, inclination, and gas scale. It omits the full hierarchical treatment.",
            "- The CMB, BAO, growth, and S8 rows use a diagnostic scalar-load parent grid passed to CAMB as a cold pressureless component. They remain diagnostic unless a finite-collar parent evaluator and custom Boltzmann module emit the required receipts.",
            "- The cluster row is an offset-timing gate. Publication use requires real lensing maps, gas maps, merger ages, and covariance.",
            "",
            "## External Measurement Sources",
            "",
            "- SPARC and RAR: https://astroweb.cwru.edu/SPARC/ and https://arxiv.org/abs/1609.05917",
            "- Planck 2018 cosmological parameters: https://arxiv.org/abs/1807.06209",
            "- DESI BAO comparison targets in the local target file: https://arxiv.org/abs/2404.03002 and https://arxiv.org/abs/2503.14738",
            "- DES Y3 plus KiDS-1000 compressed S8 row: https://arxiv.org/abs/2305.17173",
            "",
        ]
    )
    return "\n".join(lines)


def compute(args: argparse.Namespace) -> dict[str, Any]:
    parent = build_parent_grid(args)
    homogeneous = build_homogeneous_state(args)
    sparc = build_sparc(args)
    cluster = build_cluster(args, homogeneous)
    cmb = build_cmb(args)
    return {
        "status": {
            "category": "OPH dark empirical implementation scorecard",
            "publication_grade": False,
            "outputs": {
                "parent_grid": str(args.parent_grid_out),
                "json": str(args.out_json),
                "markdown": str(args.out_md),
            },
        },
        "inputs": {
            "N_scr": args.n_scr,
            "mu_eq": args.mu_eq,
            "B_A": args.B_A,
            "H0": args.H0,
            "ombh2_cmb": args.ombh2_cmb,
            "ombh2_homogeneous": args.ombh2_homogeneous,
            "sum_mnu_eV": args.sum_mnu_eV,
            "separation_kpc": args.separation_kpc,
            "time_since_passage_gyr": args.time_since_passage_gyr,
            "observed_offset_kpc": args.observed_offset_kpc,
            "offset_sigma_kpc": args.offset_sigma_kpc,
        },
        "parent_grid": parent,
        "homogeneous_state": homogeneous,
        "sparc": sparc,
        "cluster_boltzmann": cluster,
        "cmb_bao_growth": cmb,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-scr", type=float, default=3.31e122)
    parser.add_argument("--mu-eq", type=float, default=5.363470440729118)
    parser.add_argument("--B-A", dest="B_A", type=float, default=1.0)
    parser.add_argument("--n-max", type=int, default=40)
    parser.add_argument("--hold", type=float, default=0.25)
    parser.add_argument("--a-grid", default=dark_parent_collar_grid.DEFAULT_A_GRID)
    parser.add_argument("--k-grid", default=dark_parent_collar_grid.DEFAULT_K_GRID)
    parser.add_argument("--samples-json", default=None)
    parser.add_argument("--parent-grid-out", type=Path, default=DEFAULT_PARENT_GRID)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--out-md", type=Path, default=DEFAULT_OUT_MD)

    parser.add_argument("--H0", type=float, default=67.4)
    parser.add_argument("--ombh2-cmb", type=float, default=0.0224)
    parser.add_argument("--ombh2-homogeneous", type=float, default=0.02237)
    parser.add_argument(
        "--sum-mnu-eV",
        type=float,
        default=dark_cmb_bao_growth_s8_likelihood.OPH_SUM_MNU_EV,
    )
    parser.add_argument("--omega-r", type=float, default=9.17e-5)
    parser.add_argument("--tau-reio", type=float, default=0.0544)
    parser.add_argument("--As", type=float, default=2.1e-9)
    parser.add_argument("--ns", type=float, default=0.965)
    parser.add_argument("--lmax", type=int, default=1600)
    parser.add_argument("--kmax", type=float, default=2.0)
    parser.add_argument("--weak-lensing-s8", type=float, default=0.790)
    parser.add_argument("--weak-lensing-s8-sigma", type=float, default=0.016)
    parser.add_argument(
        "--weak-lensing-source",
        default="DES_Y3_KiDS_1000_joint_cosmic_shear_compressed",
    )

    parser.add_argument("--tau-values-gyr", default="0.3,0.5,1.0")
    parser.add_argument("--separation-kpc", type=float, default=200.0)
    parser.add_argument("--time-since-passage-gyr", type=float, default=0.2)
    parser.add_argument("--observed-offset-kpc", type=float, default=150.0)
    parser.add_argument("--offset-sigma-kpc", type=float, default=50.0)
    parser.add_argument("--redshifts", default="0,10,1100")

    parser.add_argument("--rar-all", type=Path, default=dark_sector_measurement_ledger.sparc_rar_compare.DEFAULT_RAR_ALL)
    parser.add_argument("--rar-bins", type=Path, default=dark_sector_measurement_ledger.sparc_rar_compare.DEFAULT_RAR_BINS)
    parser.add_argument("--mass-models", type=Path, default=dark_sector_measurement_ledger.sparc_rotation_curve_compare.DEFAULT_MASS_MODELS)
    parser.add_argument("--galaxy-sample", type=Path, default=dark_sector_measurement_ledger.sparc_systematic_likelihood.DEFAULT_GALAXY_SAMPLE)
    parser.add_argument("--fit-min-a0", type=float, default=0.2e-10)
    parser.add_argument("--fit-max-a0", type=float, default=3.0e-10)
    parser.add_argument("--upsilon-disk", type=float, default=0.5)
    parser.add_argument("--upsilon-bulge", type=float, default=0.7)
    parser.add_argument("--disk-prior", type=float, default=0.5)
    parser.add_argument("--bulge-prior", type=float, default=0.7)
    parser.add_argument("--ml-prior-sigma-dex", type=float, default=0.15)
    parser.add_argument("--intrinsic-scatter-km-s", type=float, default=8.0)
    parser.add_argument("--disk-min", type=float, default=0.05)
    parser.add_argument("--disk-max", type=float, default=1.5)
    parser.add_argument("--bulge-min", type=float, default=0.05)
    parser.add_argument("--bulge-max", type=float, default=2.0)
    parser.add_argument("--coordinate-iterations", type=int, default=3)
    parser.add_argument("--max-quality", type=int, default=2)
    parser.add_argument("--min-inclination-deg", type=float, default=30.0)
    parser.add_argument("--profile-distance", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--distance-min-mpc", type=float, default=0.1)
    parser.add_argument("--distance-sigma-window", type=float, default=3.0)
    parser.add_argument("--profile-inclination", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--inclination-min-deg", type=float, default=5.0)
    parser.add_argument("--inclination-max-deg", type=float, default=89.5)
    parser.add_argument("--inclination-sigma-window", type=float, default=3.0)
    parser.add_argument("--profile-gas", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--gas-prior-sigma-dex", type=float, default=0.08)
    parser.add_argument("--gas-min", type=float, default=0.5)
    parser.add_argument("--gas-max", type=float, default=1.5)
    parser.add_argument("--quiet", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    payload = compute(args)
    markdown = render_markdown(payload)
    write_json(args.out_json, payload)
    write_text(args.out_md, markdown)
    if not args.quiet:
        print(markdown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
