#!/usr/bin/env python3
"""Source-bound retrospective calibration of the empirical full SPARC RAR.

This module fits the standard one-parameter interpolation

    g_obs = g_bar / (1 - exp(-sqrt(g_bar / a0)))

to the committed CDS/VizieR SPARC snapshot.  It is deliberately separate from
the OPH conditional deep law in ``rar_deep_regime_diagnostic.py``.  No OPH
candidate value enters the objective, and no source derivation selects the
fitted constant.  Because both the interpolation and the SPARC scale were
known before this analysis, the receipt is a same-data calibration diagnostic,
not an independent prediction test or OPH evidence.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy import optimize


SCHEMA = "oph.cosmology.sparc_full_rar_calibration.v1"
HERE = Path(__file__).resolve().parent
RER_ROOT = HERE.parents[2]
DATA = HERE / "data"
RECEIPT_PATH = HERE / "receipts" / "sparc_full_rar_calibration.json"

SOURCE_PINS = {
    "data/PROVENANCE.md": "8805147aa9879f51ccc5280b87f71392908aa310cac0d54b021e6f70900adbbb",
    "data/ReadMe": "cbc3830f34c0379135f494719e254af0a76f75fee9ebbfe78d2e655caaaf92d5",
    "data/table1.dat": "86771a7f9f7af768ba7932dc57c6a9d81de16a287d98f4cd49bd2b4fb2a42c5c",
    "data/table2.dat": "7d027e515441c6b4ebbf6aadee0327e6ad81156c4e8b151af2f6d62cb44c3962",
}

KPC_M = 3.0856775814913673e19
KM_M = 1.0e3
A0_DISPLAY_UNIT_M_PER_S2 = 1.0e-10
BULGE_MASS_TO_LIGHT = 0.7
DISK_MASS_TO_LIGHT_VALUES = (0.4, 0.5, 0.6)
BOOTSTRAP_REPLICATES = 200
BOOTSTRAP_SEED = 1
LOG10_A0_BOUNDS = (-11.0, -9.0)


class SourceIntegrityError(RuntimeError):
    """Raised when the committed SPARC snapshot fails authentication."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise SourceIntegrityError(message)


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_json_bytes(payload: Any) -> bytes:
    return json.dumps(
        payload,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")


def _rer_relative(path: Path) -> str:
    return path.resolve().relative_to(RER_ROOT).as_posix()


def _verify_pinned_file(path: Path, expected_sha256: str) -> dict[str, Any]:
    _require(path.is_file(), f"source_missing:{path}")
    actual_sha256 = _sha256_file(path)
    _require(actual_sha256 == expected_sha256, f"source_sha256_mismatch:{path.name}")
    return {
        "path": _rer_relative(path),
        "bytes": path.stat().st_size,
        "sha256": actual_sha256,
        "code_pin_match": True,
    }


def verify_source_snapshot() -> dict[str, Any]:
    """Authenticate every committed source/provenance file against code pins."""

    files = [
        _verify_pinned_file(HERE / relative_path, expected_sha256)
        for relative_path, expected_sha256 in SOURCE_PINS.items()
    ]
    return {
        "canonical_source_binding_receipt": True,
        "catalogue": "SPARC, CDS J/AJ/152/157",
        "official_snapshot_url": "https://cdsarc.cds.unistra.fr/ftp/J/AJ/152/157/",
        "retrieved_utc_date": "2026-08-21",
        "files": files,
        "file_receipts_sha256": _sha256_bytes(_canonical_json_bytes(files)),
    }


def _float(field: str) -> float:
    field = field.strip()
    return float(field) if field else math.nan


def read_table1(path: Path = DATA / "table1.dat") -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        rows.append(
            {
                "name": line[0:11].strip(),
                "inclination_deg": _float(line[30:34]),
                "quality_flag": int(line[112:115]),
            }
        )
    _require(len(rows) == 175, "unexpected_table1_row_count")
    _require(len({row["name"] for row in rows}) == len(rows), "duplicate_table1_name")
    return rows


def read_table2(path: Path = DATA / "table2.dat") -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        rows.append(
            {
                "name": line[0:11].strip(),
                "radius_kpc": _float(line[19:25]),
                "velocity_km_per_s": _float(line[26:32]),
                "velocity_error_km_per_s": _float(line[33:38]),
                "gas_velocity_km_per_s": _float(line[39:45]),
                "disk_velocity_km_per_s": _float(line[46:52]),
                "bulge_velocity_km_per_s": _float(line[53:59]),
            }
        )
    _require(len(rows) == 3391, "unexpected_table2_row_count")
    return rows


def passes_galaxy_selection(inclination_deg: float, quality_flag: int) -> bool:
    """Apply the declared catalogue cut: Q != 3 and i >= 30 degrees."""

    return quality_flag != 3 and inclination_deg >= 30.0


def point_selection_mask(
    radius_kpc: np.ndarray,
    velocity_km_per_s: np.ndarray,
    velocity_error_km_per_s: np.ndarray,
) -> np.ndarray:
    """Apply the declared strict point cut: R>0, V>0, and dV/V<0.1."""

    return (
        (radius_kpc > 0.0)
        & (velocity_km_per_s > 0.0)
        & (
            velocity_error_km_per_s
            / np.maximum(velocity_km_per_s, 1.0e-9)
            < 0.1
        )
    )


def passes_point_selection(
    radius_kpc: float,
    velocity_km_per_s: float,
    velocity_error_km_per_s: float,
) -> bool:
    return bool(
        point_selection_mask(
            np.asarray(radius_kpc),
            np.asarray(velocity_km_per_s),
            np.asarray(velocity_error_km_per_s),
        )
    )


def collect(
    table1: list[dict[str, Any]],
    table2: list[dict[str, Any]],
    *,
    disk_mass_to_light: float,
    bulge_mass_to_light: float = BULGE_MASS_TO_LIGHT,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, int]]:
    metadata = {row["name"]: row for row in table1}
    curve_names = {row["name"] for row in table2}
    _require(curve_names == set(metadata), "catalogue_curve_name_set_mismatch")

    grouped_rows: dict[str, list[dict[str, Any]]] = {name: [] for name in metadata}
    for row in table2:
        grouped_rows[row["name"]].append(row)

    gbar_parts: list[np.ndarray] = []
    gobs_parts: list[np.ndarray] = []
    relative_error_parts: list[np.ndarray] = []
    galaxy_labels: list[str] = []
    contributing_galaxies: set[str] = set()
    points_after_galaxy_cuts = 0
    strict_boundary_points = 0
    strict_boundary_galaxies: set[str] = set()
    retained_points_at_inclination_30 = 0

    for name in sorted(metadata):
        meta = metadata[name]
        inclination_deg = float(meta["inclination_deg"])
        quality_flag = int(meta["quality_flag"])
        rows = grouped_rows[name]
        if not passes_galaxy_selection(inclination_deg, quality_flag):
            continue

        points_after_galaxy_cuts += len(rows)
        radius = np.asarray([row["radius_kpc"] for row in rows], dtype=float)
        velocity = np.asarray([row["velocity_km_per_s"] for row in rows], dtype=float)
        error = np.asarray(
            [row["velocity_error_km_per_s"] for row in rows], dtype=float
        )
        gas_velocity = np.asarray(
            [row["gas_velocity_km_per_s"] for row in rows], dtype=float
        )
        disk_velocity = np.asarray(
            [row["disk_velocity_km_per_s"] for row in rows], dtype=float
        )
        bulge_velocity = np.asarray(
            [row["bulge_velocity_km_per_s"] for row in rows], dtype=float
        )

        relative_error = error / np.maximum(velocity, 1.0e-9)
        point_mask = point_selection_mask(radius, velocity, error)
        boundary_mask = (
            (radius > 0.0)
            & (velocity > 0.0)
            & np.isclose(relative_error, 0.1, rtol=0.0, atol=1.0e-15)
            & ~point_mask
        )
        if bool(boundary_mask.any()):
            strict_boundary_galaxies.add(name)
            strict_boundary_points += int(boundary_mask.sum())

        mask = point_mask
        radius_m = radius * KPC_M
        gbar = (
            gas_velocity * np.abs(gas_velocity)
            + disk_mass_to_light * disk_velocity * np.abs(disk_velocity)
            + bulge_mass_to_light * bulge_velocity * np.abs(bulge_velocity)
        ) * KM_M**2 / radius_m
        gobs = velocity**2 * KM_M**2 / radius_m
        mask &= gbar > 0.0

        retained_count = int(mask.sum())
        if retained_count:
            contributing_galaxies.add(name)
        if inclination_deg == 30.0:
            retained_points_at_inclination_30 += retained_count
        gbar_parts.append(gbar[mask])
        gobs_parts.append(gobs[mask])
        relative_error_parts.append(relative_error[mask])
        galaxy_labels.extend([name] * retained_count)

    _require(bool(gbar_parts), "no_sparc_points_collected")
    gbar_all = np.concatenate(gbar_parts)
    gobs_all = np.concatenate(gobs_parts)
    relative_error_all = np.concatenate(relative_error_parts)
    galaxies_all = np.asarray(galaxy_labels)

    quality_excluded = {
        row["name"] for row in table1 if int(row["quality_flag"]) == 3
    }
    inclination_excluded = {
        row["name"] for row in table1 if float(row["inclination_deg"]) < 30.0
    }
    passing_catalogue_cuts = {
        row["name"]
        for row in table1
        if passes_galaxy_selection(
            float(row["inclination_deg"]), int(row["quality_flag"])
        )
    }
    inclination_equal_30 = {
        row["name"] for row in table1 if float(row["inclination_deg"]) == 30.0
    }
    census = {
        "catalogue_galaxies": len(table1),
        "rotation_curve_points": len(table2),
        "galaxies_excluded_quality_q_equals_3": len(quality_excluded),
        "galaxies_excluded_inclination_lt_30": len(inclination_excluded),
        "galaxies_excluded_by_both_catalogue_cuts": len(
            quality_excluded & inclination_excluded
        ),
        "galaxies_passing_catalogue_cuts": len(passing_catalogue_cuts),
        "galaxies_with_retained_points": len(contributing_galaxies),
        "galaxies_passing_catalogue_cuts_with_zero_retained_points": len(
            passing_catalogue_cuts - contributing_galaxies
        ),
        "catalogue_galaxies_at_inclination_eq_30": len(inclination_equal_30),
        "selected_galaxies_at_inclination_eq_30": len(
            inclination_equal_30 & passing_catalogue_cuts
        ),
        "contributing_galaxies_at_inclination_eq_30": len(
            inclination_equal_30 & contributing_galaxies
        ),
        "retained_points_at_inclination_eq_30": retained_points_at_inclination_30,
        "points_after_catalogue_galaxy_cuts": points_after_galaxy_cuts,
        "points_rejected_at_strict_relative_error_eq_0_1": strict_boundary_points,
        "galaxies_lost_at_strict_relative_error_eq_0_1": len(
            strict_boundary_galaxies - contributing_galaxies
        ),
        "retained_points": len(gbar_all),
    }
    return gbar_all, gobs_all, galaxies_all, relative_error_all, census


def rar(gbar: np.ndarray, a0: float) -> np.ndarray:
    return gbar / (1.0 - np.exp(-np.sqrt(gbar / a0)))


def fit(
    gbar: np.ndarray,
    gobs: np.ndarray,
    weights: np.ndarray | None = None,
) -> tuple[float, float, int]:
    _require(len(gbar) >= 2 and len(gbar) == len(gobs), "invalid_fit_sample")
    if weights is None:
        weights = np.ones(len(gbar))
    _require(
        len(weights) == len(gbar)
        and bool(np.all(np.isfinite(weights)))
        and bool(np.all(weights > 0.0)),
        "invalid_fit_weights",
    )

    def objective(log10_a0: float) -> float:
        residual = np.log10(gobs) - np.log10(rar(gbar, 10.0**log10_a0))
        return float(np.sum(weights * residual**2))

    result = optimize.minimize_scalar(
        objective,
        bounds=LOG10_A0_BOUNDS,
        method="bounded",
    )
    _require(bool(result.success), "a0_optimization_failed")
    a0 = float(10.0**result.x)
    residual = np.log10(gobs) - np.log10(rar(gbar, a0))
    return a0, float(np.std(residual, ddof=0)), len(gbar)


def weighting_sensitivity(
    gbar: np.ndarray,
    gobs: np.ndarray,
    galaxies: np.ndarray,
    relative_velocity_error: np.ndarray,
    primary_a0: float,
) -> dict[str, Any]:
    """Compute two named estimator sensitivities, neither a full likelihood."""

    names, counts = np.unique(galaxies, return_counts=True)
    counts_by_name = dict(zip(names.tolist(), counts.tolist(), strict=True))
    equal_galaxy_weights = np.asarray(
        [1.0 / counts_by_name[name] for name in galaxies],
        dtype=float,
    )
    _require(
        bool(np.all(relative_velocity_error > 0.0)),
        "nonpositive_relative_velocity_error",
    )
    # The constant 2/ln(10) converting dV/V to first-order sigma(log10 gobs)
    # cancels from the optimum, so inverse fractional-velocity variance is
    # sufficient. This remains deliberately naive because every gbar and
    # correlated/nuisance uncertainty is omitted.
    velocity_error_only_weights = 1.0 / relative_velocity_error**2

    rows: dict[str, Any] = {}
    for name, weights, warning in (
        (
            "equal_parent_galaxy_total_weight",
            equal_galaxy_weights,
            "Each contributing galaxy receives equal total weight; this is an estimator sensitivity, not a calibrated likelihood.",
        ),
        (
            "inverse_velocity_fractional_variance_only",
            velocity_error_only_weights,
            "Uses only dV/V in the gobs channel and omits gbar, distance, inclination, stellar-population, intrinsic-scatter, and correlation uncertainties; it is not a recommended likelihood.",
        ),
    ):
        a0 = fit(gbar, gobs, weights=weights)[0]
        rows[name] = {
            "a0_si_m_per_s2": a0,
            "a0_in_1e_minus_10_m_per_s2": a0 / A0_DISPLAY_UNIT_M_PER_S2,
            "relative_to_primary_point_weighted_percent": 100.0
            * (a0 - primary_a0)
            / primary_a0,
            "eligible_as_inference_or_oph_evidence": False,
            "warning": warning,
        }
    return rows


def fit_with_galaxy_bootstrap(
    gbar: np.ndarray,
    gobs: np.ndarray,
    galaxies: np.ndarray,
) -> dict[str, Any]:
    a0, scatter, point_count = fit(gbar, gobs)
    galaxy_names = np.unique(galaxies)
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    bootstrap_a0: list[float] = []
    for _ in range(BOOTSTRAP_REPLICATES):
        sampled_names = rng.choice(
            galaxy_names,
            size=len(galaxy_names),
            replace=True,
        )
        indices = np.concatenate(
            [np.flatnonzero(galaxies == name) for name in sampled_names]
        )
        bootstrap_a0.append(fit(gbar[indices], gobs[indices])[0])
    bootstrap_std = float(np.std(bootstrap_a0, ddof=0))
    return {
        "a0_si_m_per_s2": a0,
        "a0_in_1e_minus_10_m_per_s2": a0 / A0_DISPLAY_UNIT_M_PER_S2,
        "galaxy_bootstrap_std_si_m_per_s2": bootstrap_std,
        "galaxy_bootstrap_std_in_1e_minus_10_m_per_s2": (
            bootstrap_std / A0_DISPLAY_UNIT_M_PER_S2
        ),
        "scatter_dex": scatter,
        "n_points": point_count,
        "n_galaxies": len(galaxy_names),
    }


def build_report() -> dict[str, Any]:
    source_binding = verify_source_snapshot()
    _require(
        bool(source_binding["canonical_source_binding_receipt"]),
        "canonical_source_binding_not_earned",
    )
    table1 = read_table1()
    table2 = read_table2()

    fits: dict[str, dict[str, Any]] = {}
    sensitivity_rows: dict[str, dict[str, Any]] = {}
    canonical_census: dict[str, int] | None = None
    for disk_mass_to_light in DISK_MASS_TO_LIGHT_VALUES:
        gbar, gobs, galaxies, relative_velocity_error, census = collect(
            table1,
            table2,
            disk_mass_to_light=disk_mass_to_light,
        )
        if canonical_census is None:
            canonical_census = census
        else:
            _require(census == canonical_census, "mass_to_light_census_drift")
        fit_key = f"disk_mass_to_light_{disk_mass_to_light:.1f}"
        fits[fit_key] = (
            fit_with_galaxy_bootstrap(gbar, gobs, galaxies)
        )
        sensitivity_rows[fit_key] = weighting_sensitivity(
            gbar,
            gobs,
            galaxies,
            relative_velocity_error,
            float(fits[fit_key]["a0_si_m_per_s2"]),
        )
    _require(canonical_census is not None, "missing_selection_census")

    return {
        "schema": SCHEMA,
        "scope": "retrospective_same_data_empirical_full_rar_calibration",
        "receipts": {
            "canonical_public_source_binding_receipt": True,
            "deterministic_calibration_reproduction_receipt": True,
            "public_data_attached": True,
            "independent_held_out_test_receipt": False,
            "oph_specific_empirical_evidence_receipt": False,
            "oph_validation_or_falsification_receipt": False,
        },
        "interpretation": {
            "classification": "retrospective_same_data_calibration_diagnostic",
            "registered_before_data_exposure": False,
            "independent_held_out_data": False,
            "oph_candidate_values_consumed_by_fit": False,
            "source_derived_a0": False,
            "oph_specific_empirical_evidence": False,
            "verdict": "NOT_ELIGIBLE_AS_OPH_EVIDENCE",
            "warnings": [
                "The empirical interpolation and the SPARC acceleration scale were known before this calibration.",
                "The fit is an unweighted point-level log-residual diagnostic, not a complete measurement likelihood.",
                "Distance, inclination, gas, and stellar mass-to-light nuisance uncertainties are not marginalized.",
                "The galaxy bootstrap is conditional cluster-resampling dispersion, not total statistical or systematic uncertainty.",
                "The disk mass-to-light scan with a fixed bulge value is a sensitivity scan, not propagated systematic uncertainty.",
            ],
        },
        "source_binding": source_binding,
        "producer": {
            "script_path": _rer_relative(Path(__file__)),
            "script_sha256": _sha256_file(Path(__file__).resolve()),
        },
        "analysis_contract": {
            "rar_function": "g_obs = g_bar / (1 - exp(-sqrt(g_bar/a0)))",
            "objective": "unweighted_sum_squared_log10_gobs_residuals",
            "galaxy_selection": {
                "quality_flag": "Q != 3",
                "inclination_deg": "i >= 30",
            },
            "point_selection": {
                "radius_kpc": "R > 0",
                "velocity": "V > 0",
                "relative_velocity_error": "dV/V < 0.1",
                "baryonic_acceleration": "g_bar > 0",
            },
            "bulge_mass_to_light": BULGE_MASS_TO_LIGHT,
            "disk_mass_to_light_values": list(DISK_MASS_TO_LIGHT_VALUES),
            "bootstrap_unit": "galaxy",
            "bootstrap_replicates": BOOTSTRAP_REPLICATES,
            "bootstrap_seed": BOOTSTRAP_SEED,
            "standard_deviation_ddof": 0,
            "display_acceleration_unit_m_per_s2": A0_DISPLAY_UNIT_M_PER_S2,
        },
        "sample_census": canonical_census,
        "fits": fits,
        "weighting_sensitivity": {
            "status": "DETERMINISTIC_ESTIMATOR_SENSITIVITY__NOT_INFERENCE",
            "rows": sensitivity_rows,
        },
        "reference_literature": {
            "citation": "McGaugh, Lelli and Schombert (2016), arXiv:1609.05917",
            "reported_a0_si_m_per_s2": 1.20e-10,
            "reported_random_uncertainty_si_m_per_s2": 0.02e-10,
            "reported_systematic_uncertainty_si_m_per_s2": 0.24e-10,
            "reported_scatter_dex": 0.13,
            "reported_points": 2693,
            "reported_galaxies": 153,
            "census_reconciliation_status": "THREE_POINT_DIFFERENCE_OPEN__153_GALAXY_SELECTION_REPRODUCED",
            "warning": "The Q != 3 and i >= 30 selection reproduces the published 153-galaxy parent sample. After the strict point cut, 147 galaxies contribute 2696 points, three more than the published 2693; only that point-level difference remains to reconcile.",
        },
    }


def render_report(report: dict[str, Any]) -> bytes:
    return (json.dumps(report, indent=2, sort_keys=True) + "\n").encode("utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=RECEIPT_PATH)
    parser.add_argument(
        "--check",
        action="store_true",
        help="require the canonical receipt to be byte-identical to recomputation",
    )
    args = parser.parse_args(argv)
    _require(
        args.output.resolve() == RECEIPT_PATH.resolve(),
        "canonical_receipt_output_path_mismatch",
    )

    output_bytes = render_report(build_report())
    if args.check:
        if not RECEIPT_PATH.is_file() or RECEIPT_PATH.read_bytes() != output_bytes:
            raise SystemExit("stale_sparc_full_rar_calibration_receipt")
        print(f"PASS: {RECEIPT_PATH} is current and source-bound")
        return 0

    RECEIPT_PATH.parent.mkdir(parents=True, exist_ok=True)
    RECEIPT_PATH.write_bytes(output_bytes)
    print(output_bytes.decode("utf-8"), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
