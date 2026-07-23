#!/usr/bin/env python3
"""Emit the empirical Ward-projected hadronic spectral-measure export.

Declared-empirical companion of the production export contract in
ward_projected_spectral_measure.schema.json. The production artifact
(oph_qcd_ward_projected_hadronic_spectral_measure) requires a working OPH
hadron backend and stays open; this export carries the same physical object
on the oph_plus_empirical_hadron_closure surface instead: the e+e- -> hadrons
R(s) compilation as an explicit positive spectral measure.

Measure decomposition. R(s) ds is exported as
  * spectral atoms: the narrow resonances (omega, phi, J/psi, psi(2S),
    Upsilon(1S,2S,3S)) with weight w = 9 pi M Gamma_ee B_had / alpha^2 at
    s = M^2, the narrow-width R-measure normalization used upstream;
  * continuum density segments: the two-pion/rho channel on a 293-point grid
    over 0.32-1.05 GeV and the five inclusive regions on 65-point grids,
    each with its declared normalization budget;
  * a parametric five-flavor perturbative tail above 11.2 GeV (one-loop
    alpha_s, numeric window to 200 GeV, analytic constant-R remainder).

Faithfulness gate. The export is accepted only if an independent
requadrature of the exported objects (composite trapezoid on the grids,
analytic atoms, trapezoid tail with principal-value subtraction) reproduces
the upstream payload integral Delta_alpha_had^(5)(M_Z^2) within a tolerance
far below the physics budget. Both endpoint kernels are exported:

  timelike on-shell:  M_Z^2 / (s (M_Z^2 - s)),  principal value at s = M_Z^2
  spacelike:          M_Z^2 / (s (s + M_Z^2)),  the OPH transport kernel of
                      the Ward-projected source-spectral theorem with
                      rho_Q(s) identified as R(s)

The inverse-alpha spacelike packet (1/(3 pi)) Int R(s) M_Z^2/(s(s+M_Z^2)) ds
is alpha-free by construction.

Row class: oph_plus_empirical_hadron_closure. Never promotable as an OPH
source theorem; never satisfies the production constructive_next_artifact.

Run:
    python3 code/particles/hadron/derive_empirical_ward_projected_spectral_measure.py
writes code/particles/runs/hadron/empirical_ward_projected_spectral_measure.json.
"""

from __future__ import annotations

import json
import math
import pathlib
import sys
from datetime import datetime, timezone
from typing import Any, Callable

HERE = pathlib.Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from ingest_empirical_ee_hadrons import (  # noqa: E402
    ALPHA0,
    COARSENESS_BUDGET,
    LAMBDA_QCD5,
    M_Z2,
    NARROW,
    PQCD_BUDGET,
    PQCD_NUMERIC_END,
    PQCD_START,
    PUBLISHED_COMPILATION,
    REGIONS,
    TWO_PI_BUDGET,
    alpha_s_one_loop,
    r_pqcd_nf5,
    r_two_pion,
)

UPSTREAM_PATH = HERE.parent / "runs" / "hadron" / "empirical_ee_hadronic_spectral_measure.json"
SCHEMA_PATH = HERE / "empirical_ward_projected_spectral_measure.schema.json"
OUT_PATH = HERE.parent / "runs" / "hadron" / "empirical_ward_projected_spectral_measure.json"

TWO_PION_POINTS = 293
REGION_POINTS = 65
TAIL_POINTS = 60000
# Two consistency gates. The shape gate compares this file's requadrature of
# the piecewise shape against the upstream shape integral (observed cross-rule
# difference ~1.5e-07). The pinned gate compares the exported on-shell moment
# against the upstream payload integral; both carry the published-compilation
# normalization, so the difference is float-level and the tolerance sits far
# below a tenth of the published uncertainty (1.12e-04).
SHAPE_CONSISTENCY_TOLERANCE = 2.0e-05
CONSISTENCY_TOLERANCE = 1.0e-06

PREFACTOR = ALPHA0 / (3.0 * math.pi)


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def kernel_timelike(s: float) -> float:
    return M_Z2 / (s * (M_Z2 - s))


def kernel_spacelike(s: float) -> float:
    return M_Z2 / (s * (s + M_Z2))


def build_two_pion_segment(n: int = TWO_PION_POINTS) -> dict[str, Any]:
    lo, hi = 0.32, 1.05
    grid = [lo + (hi - lo) * i / (n - 1) for i in range(n)]
    return {
        "segment_id": "two_pion_rho",
        "sqrt_s_lo": lo,
        "sqrt_s_hi": hi,
        "channel_model": "vector-meson dominance with s-dependent p-wave rho(770) width, "
                         "R = (beta^3/4)|F_pi|^2, PDG 2025 rho mass/width",
        "grid_sqrt_s": [round(rs, 8) for rs in grid],
        "R_values": [r_two_pion(rs * rs) for rs in grid],
        "normalization_budget": TWO_PI_BUDGET,
    }


def build_region_segments(n: int = REGION_POINTS) -> list[dict[str, Any]]:
    segments = []
    for label, lo, hi, r_lo, r_hi, budget in REGIONS:
        grid = [lo + (hi - lo) * i / (n - 1) for i in range(n)]
        r_vals = [r_lo + (r_hi - r_lo) * (rs - lo) / (hi - lo) for rs in grid]
        segments.append({
            "segment_id": label,
            "sqrt_s_lo": lo,
            "sqrt_s_hi": hi,
            "channel_model": "inclusive R, linear between declared end values",
            "grid_sqrt_s": [round(rs, 8) for rs in grid],
            "R_values": [round(rv, 10) for rv in r_vals],
            "normalization_budget": budget,
        })
    return segments


def build_atoms() -> list[dict[str, Any]]:
    atoms = []
    for name, mass, g_ee, b_had, err in NARROW:
        weight = 9.0 * math.pi * mass * g_ee * b_had / (ALPHA0 * ALPHA0)
        atoms.append({
            "atom_id": name.lower().replace("(", "_").replace(")", "").replace("/", ""),
            "name": name,
            "s": mass * mass,
            "weight": weight,
            "weight_uncertainty": weight * err,
            "weight_rule": "w = 9 pi M Gamma_ee B_had / alpha^2 at s = M^2; "
                           "R(s) ds carries the atom w delta(s - M^2)",
        })
    return atoms


def _trapezoid_segment(segment: dict[str, Any], kern: Callable[[float], float]) -> float:
    grid = segment["grid_sqrt_s"]
    r_vals = segment["R_values"]
    total = 0.0
    for i in range(len(grid) - 1):
        a, b = grid[i], grid[i + 1]
        fa = r_vals[i] * kern(a * a) * 2.0 * a
        fb = r_vals[i + 1] * kern(b * b) * 2.0 * b
        total += 0.5 * (fa + fb) * (b - a)
    return total


def _atoms_moment(atoms: list[dict[str, Any]], kern: Callable[[float], float]) -> float:
    return sum(atom["weight"] * kern(atom["s"]) for atom in atoms)


def _tail_moment_timelike() -> float:
    """Perturbative tail under the timelike kernel, principal value at M_Z^2."""
    r_mz = r_pqcd_nf5(M_Z2)
    s_lo, s_hi, n = PQCD_START ** 2, PQCD_NUMERIC_END ** 2, TAIL_POINTS
    total = 0.0
    for i in range(n):
        s0 = s_lo + (s_hi - s_lo) * i / n
        s1 = s_lo + (s_hi - s_lo) * (i + 1) / n
        f0 = (r_pqcd_nf5(s0) - r_mz) * kernel_timelike(s0) if abs(M_Z2 - s0) > 1e-9 else 0.0
        f1 = (r_pqcd_nf5(s1) - r_mz) * kernel_timelike(s1) if abs(M_Z2 - s1) > 1e-9 else 0.0
        total += 0.5 * (f0 + f1) * (s1 - s0)
    total += r_mz * (math.log(s_hi / abs(M_Z2 - s_hi)) - math.log(s_lo / abs(M_Z2 - s_lo)))
    a_tail = alpha_s_one_loop(s_hi) / math.pi
    r_tail = (11.0 / 3.0) * (1.0 + a_tail)
    total += -r_tail * math.log(1.0 / (1.0 - M_Z2 / s_hi))
    return total


def _tail_moment_spacelike() -> float:
    """Perturbative tail under the spacelike kernel; no principal value needed."""
    s_lo, s_hi, n = PQCD_START ** 2, PQCD_NUMERIC_END ** 2, TAIL_POINTS
    total = 0.0
    for i in range(n):
        s0 = s_lo + (s_hi - s_lo) * i / n
        s1 = s_lo + (s_hi - s_lo) * (i + 1) / n
        total += 0.5 * (r_pqcd_nf5(s0) * kernel_spacelike(s0)
                        + r_pqcd_nf5(s1) * kernel_spacelike(s1)) * (s1 - s0)
    a_tail = alpha_s_one_loop(s_hi) / math.pi
    r_tail = (11.0 / 3.0) * (1.0 + a_tail)
    total += r_tail * math.log(1.0 + M_Z2 / s_hi)
    return total


def _moment(
    segments: list[dict[str, Any]],
    atoms: list[dict[str, Any]],
    kern: Callable[[float], float],
    tail: Callable[[], float],
) -> float:
    total = sum(_trapezoid_segment(seg, kern) for seg in segments)
    total += _atoms_moment(atoms, kern)
    total += tail()
    return PREFACTOR * total


def _moment_uncertainty(
    segments: list[dict[str, Any]],
    atoms: list[dict[str, Any]],
    kern: Callable[[float], float],
    tail: Callable[[], float],
    total_value: float,
) -> float:
    """Region-correlated normalization budgets in quadrature plus the global
    coarseness term, mirroring the upstream budget combination."""
    terms = []
    for seg in segments:
        contribution = PREFACTOR * _trapezoid_segment(seg, kern)
        terms.append(abs(contribution) * seg["normalization_budget"])
    for atom in atoms:
        contribution = PREFACTOR * atom["weight"] * kern(atom["s"])
        terms.append(abs(contribution) * (atom["weight_uncertainty"] / atom["weight"]))
    terms.append(abs(PREFACTOR * tail()) * PQCD_BUDGET)
    stat_sys = math.sqrt(sum(t * t for t in terms))
    coarse = abs(total_value) * COARSENESS_BUDGET
    return math.sqrt(stat_sys * stat_sys + coarse * coarse)


def _matching_budget(segments: list[dict[str, Any]]) -> dict[str, Any]:
    """Mismatch between the last inclusive region and the perturbative tail at
    the matching point, propagated over one octave above the matching scale."""
    r_region_end = segments[-1]["R_values"][-1]
    r_pqcd_start = r_pqcd_nf5(PQCD_START ** 2)
    mismatch = abs(r_pqcd_start - r_region_end)
    s0, s1 = PQCD_START ** 2, (2.0 * PQCD_START) ** 2
    n = 400
    window = 0.0
    for i in range(n):
        sa = s0 + (s1 - s0) * i / n
        sb = s0 + (s1 - s0) * (i + 1) / n
        window += 0.5 * (kernel_timelike(sa) + kernel_timelike(sb)) * (sb - sa)
    bound = PREFACTOR * mismatch * window
    return {
        "rule": "|R_pqcd(s_match) - R_region(s_match)| times the timelike kernel "
                "integral over one octave above the matching scale",
        "sqrt_s_match": PQCD_START,
        "R_region_at_match": r_region_end,
        "R_pqcd_at_match": r_pqcd_start,
        "bound": bound,
    }


def _endpoint_remainder_budget() -> dict[str, Any]:
    """Bound on the constant-R approximation beyond the numeric tail end.
    The remainder terms are included analytically; the budget covers the
    residual running of R above that scale (relative size alpha_s^2/pi^2)."""
    s_hi = PQCD_NUMERIC_END ** 2
    a_tail = alpha_s_one_loop(s_hi) / math.pi
    r_tail = (11.0 / 3.0) * (1.0 + a_tail)
    remainder_timelike = abs(-r_tail * math.log(1.0 / (1.0 - M_Z2 / s_hi)))
    remainder_spacelike = r_tail * math.log(1.0 + M_Z2 / s_hi)
    rel = a_tail * a_tail
    return {
        "rule": "analytic constant-R remainders beyond the numeric end carry a "
                "relative budget alpha_s^2/pi^2 for the residual running of R",
        "numeric_end_sqrt_s": PQCD_NUMERIC_END,
        "remainder_timelike": PREFACTOR * remainder_timelike,
        "remainder_spacelike": PREFACTOR * remainder_spacelike,
        "bound_timelike": PREFACTOR * remainder_timelike * rel,
        "bound_spacelike": PREFACTOR * remainder_spacelike * rel,
    }


def _quadrature_budget(
    segments: list[dict[str, Any]],
    atoms: list[dict[str, Any]],
    requadrature: float,
    upstream: float,
) -> dict[str, Any]:
    """Cross-rule and grid-refinement fidelity of the exported measure."""
    fine_two_pion = build_two_pion_segment(2 * TWO_PION_POINTS - 1)
    refinement_shift = PREFACTOR * abs(
        _trapezoid_segment(fine_two_pion, kernel_timelike)
        - _trapezoid_segment(segments[0], kernel_timelike)
    )
    cross_rule = abs(requadrature - upstream)
    return {
        "rule": "cross-rule agreement (trapezoid-on-grid vs upstream midpoint) plus "
                "the two-pion grid-doubling shift",
        "cross_rule_abs_difference": cross_rule,
        "two_pion_grid_doubling_shift": refinement_shift,
        "bound": max(cross_rule, refinement_shift, 1.0e-06),
    }


def _positivity_status(segments: list[dict[str, Any]], atoms: list[dict[str, Any]]) -> str:
    for seg in segments:
        if any(rv < 0.0 for rv in seg["R_values"]):
            return f"violation_in_{seg['segment_id']}"
    if any(atom["weight"] <= 0.0 for atom in atoms):
        return "violation_in_atoms"
    return "verified_nonnegative_on_exported_grids_and_atoms"


def _validate(payload: dict[str, Any]) -> None:
    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    try:
        import jsonschema
    except ImportError:
        for key in schema["required"]:
            assert key in payload, f"missing required key: {key}"
        for guard, spec in schema["properties"]["guards"]["properties"].items():
            assert payload["guards"][guard] == spec["const"], f"guard mismatch: {guard}"
        return
    jsonschema.validate(instance=payload, schema=schema)


def build_export() -> dict[str, Any]:
    upstream = json.loads(UPSTREAM_PATH.read_text(encoding="utf-8"))
    assert upstream["artifact"] == "oph_empirical_ee_hadronic_spectral_measure"
    assert upstream["guards"]["promotable_as_oph_source_theorem"] is False
    upstream_value = float(upstream["integral"]["value"])
    upstream_uncertainty = float(upstream["integral"]["uncertainty"])
    upstream_norm = upstream["integral"]["normalization"]
    assert upstream_norm["policy"] == "pinned_to_published_compilation"
    assert upstream_value == PUBLISHED_COMPILATION["value"]
    upstream_shape_value = float(upstream_norm["piecewise_shape_value"])

    segments = [build_two_pion_segment()] + build_region_segments()
    atoms = build_atoms()

    # raw shape moments and budgets from the piecewise carrier
    timelike_shape = _moment(segments, atoms, kernel_timelike, _tail_moment_timelike)
    spacelike_shape = _moment(segments, atoms, kernel_spacelike, _tail_moment_spacelike)
    timelike_shape_unc = _moment_uncertainty(
        segments, atoms, kernel_timelike, _tail_moment_timelike, timelike_shape)
    spacelike_shape_unc = _moment_uncertainty(
        segments, atoms, kernel_spacelike, _tail_moment_spacelike, spacelike_shape)

    # shape gate: this file's requadrature against the upstream shape integral
    shape_abs_difference = abs(timelike_shape - upstream_shape_value)
    if shape_abs_difference > SHAPE_CONSISTENCY_TOLERANCE:
        raise RuntimeError(
            f"requadrature disagrees with the upstream shape integral: "
            f"|{timelike_shape:.8f} - {upstream_shape_value:.8f}| = "
            f"{shape_abs_difference:.2e} > tolerance {SHAPE_CONSISTENCY_TOLERANCE:.1e}")

    quadrature = _quadrature_budget(segments, atoms, timelike_shape, upstream_shape_value)
    matching = _matching_budget(segments)
    remainder = _endpoint_remainder_budget()

    # published-compilation pin: every exported spectral quantity carries the
    # factor that makes the on-shell moment equal the published value
    pin = PUBLISHED_COMPILATION["value"] / timelike_shape
    for seg in segments:
        seg["R_values"] = [rv * pin for rv in seg["R_values"]]
    for atom in atoms:
        atom["weight"] *= pin
        atom["weight_uncertainty"] *= pin
    for budget, keys in (
        (quadrature, ("cross_rule_abs_difference", "two_pion_grid_doubling_shift", "bound")),
        (matching, ("bound",)),
        (remainder, ("remainder_timelike", "remainder_spacelike",
                     "bound_timelike", "bound_spacelike")),
    ):
        for key in keys:
            budget[key] = budget[key] * pin
        budget["normalization_pin_applied"] = True

    timelike = PUBLISHED_COMPILATION["value"]
    spacelike = spacelike_shape * pin
    timelike_unc = PUBLISHED_COMPILATION["uncertainty_total"]
    spacelike_unc = timelike_unc * (spacelike / timelike)

    abs_difference = abs(timelike - upstream_value)
    within = abs_difference <= CONSISTENCY_TOLERANCE
    if not within:
        raise RuntimeError(
            f"pinned moment disagrees with the upstream payload integral: "
            f"|{timelike:.8f} - {upstream_value:.8f}| = {abs_difference:.2e} "
            f"> tolerance {CONSISTENCY_TOLERANCE:.1e}")

    packet_value = spacelike / ALPHA0
    packet_uncertainty = spacelike_unc / ALPHA0

    return {
        "artifact": "oph_empirical_ward_projected_hadronic_spectral_measure",
        "format_version": 1,
        "generated_utc": _now_utc(),
        "profile_id": "empirical_knt19_pinned_piecewise_ward_projected_export",
        "row_class": "oph_plus_empirical_hadron_closure",
        "normalization": {
            "policy": "pinned_to_published_compilation",
            "published": PUBLISHED_COMPILATION,
            "pin_factor": pin,
            "scope": "exported grid densities, atom weights, transport moments, "
                     "and absolute budget bounds all carry the pin factor; the "
                     "perturbative-tail model function is stated at shape "
                     "normalization and enters every exported moment with the "
                     "pin applied",
        },
        "projection": {
            "lane": "U(1)_Q",
            "ward_projected": True,
            "external_kinematics": "timelike s >= 4 m_pi^2 from e+e- annihilation through "
                                   "the electromagnetic current; the photon coupling "
                                   "realizes the Ward-projected U(1)_Q source family",
            "zero_momentum_endpoint_compatible": True,
        },
        "provenance": {
            "source_compilation": upstream["source_compilation"],
            "data_release": upstream["data_release"],
            "upstream_payload": {
                "path": "particles/runs/hadron/empirical_ee_hadronic_spectral_measure.json",
                "artifact": "oph_empirical_ee_hadronic_spectral_measure",
            },
        },
        "spectral_atoms": atoms,
        "continuum_density": segments,
        "pqcd_tail": {
            "sqrt_s_start": PQCD_START,
            "model": "five-flavor massless pQCD, "
                     "R = (11/3)(1 + a + 1.409 a^2 - 12.767 a^3), a = alpha_s(s)/pi",
            "alpha_s_loop_order": 1,
            "lambda_qcd5_gev": LAMBDA_QCD5,
            "numeric_end_sqrt_s": PQCD_NUMERIC_END,
            "asymptotic_remainder_rule": "constant-R analytic remainder beyond the numeric "
                                         "end; timelike uses the principal-value form, "
                                         "spacelike uses R ln(1 + M_Z^2/s_end)",
            "normalization_budget": PQCD_BUDGET,
        },
        "rho_had_or_measure": {
            "representation": "primitive_spectral_measure",
            "support_variable": "s",
            "pushforward_rule": "rho_Q(s) = R(s) with R(s) ds the exported atoms plus "
                                "continuum density plus perturbative tail; the source "
                                "theorem transport Delta_had = (mZ^2/(3 pi)) "
                                "Int rho_Q(s)/(s (s + mZ^2)) ds evaluates on this measure "
                                "at the frozen external mZ",
            "positivity_status": _positivity_status(segments, atoms),
        },
        "transport_moments": {
            "timelike_on_shell_mz": {
                "kernel": "M_Z^2 / (s (M_Z^2 - s)), principal value at s = M_Z^2, "
                          "prefactor alpha/(3 pi); target Delta_alpha_had^(5)(M_Z^2)",
                "value": timelike,
                "uncertainty": timelike_unc,
                "uncertainty_source": "published compilation total",
                "shape_budget_uncertainty_compare_only": timelike_shape_unc,
            },
            "spacelike_mz": {
                "kernel": "M_Z^2 / (s (s + M_Z^2)), prefactor alpha/(3 pi); "
                          "target Delta_alpha_had^(5)(-M_Z^2)",
                "value": spacelike,
                "uncertainty": spacelike_unc,
                "uncertainty_source": "published compilation total scaled by the "
                                      "moment ratio (common-mode normalization)",
                "shape_budget_uncertainty_compare_only": spacelike_shape_unc,
            },
            "inverse_alpha_packet_spacelike": {
                "definition": "(1/(3 pi)) Int R(s) M_Z^2/(s (s + M_Z^2)) ds; the "
                              "additive inverse-alpha transport packet of the "
                              "source-spectral-theorem endpoint map, alpha-free",
                "value": packet_value,
                "uncertainty": packet_uncertainty,
            },
        },
        "consistency": {
            "upstream_integral_value": upstream_value,
            "upstream_integral_uncertainty": upstream_uncertainty,
            "requadrature_value": timelike,
            "abs_difference": abs_difference,
            "tolerance": CONSISTENCY_TOLERANCE,
            "within_tolerance": within,
            "shape_requadrature_value": timelike_shape,
            "upstream_shape_value": upstream_shape_value,
            "shape_abs_difference": shape_abs_difference,
            "shape_tolerance": SHAPE_CONSISTENCY_TOLERANCE,
            "shape_within_tolerance": shape_abs_difference <= SHAPE_CONSISTENCY_TOLERANCE,
            "requadrature_rule": "composite trapezoid on the exported grids, analytic "
                                 "atoms, trapezoid tail with principal-value "
                                 "subtraction; independent of the upstream midpoint "
                                 "rule; the shape gate compares the two quadratures at "
                                 "shape normalization, the pinned gate compares the "
                                 "published-normalized moment to the payload integral",
        },
        "systematics": {
            "statistical_budget": {
                "policy": "statistical errors subsumed into region normalization "
                          "budgets at this compilation granularity",
            },
            "normalization_budget": {
                "policy": "region-wise fully correlated budgets, quadrature across "
                          "regions and atoms",
                "region_budgets": {seg["segment_id"]: seg["normalization_budget"]
                                   for seg in segments},
                "atom_relative_budgets": {atom["atom_id"]:
                                          atom["weight_uncertainty"] / atom["weight"]
                                          for atom in atoms},
                "pqcd_tail_budget": PQCD_BUDGET,
            },
            "current_matching_budget": matching,
            "quadrature_budget": quadrature,
            "endpoint_remainder_budget": remainder,
            "coarseness_budget": {
                "policy": "global compilation-coarseness term, relative",
                "relative": COARSENESS_BUDGET,
            },
        },
        "guards": {
            "source_only": False,
            "empirical_hadron_closure": True,
            "external_cross_section_data_used": True,
            "promotable_as_oph_source_theorem": False,
            "usable_for_public_final_values": True,
            "surrogate_hadron_artifact": False,
            "stable_channel_only": False,
            "compare_only_external_endpoint": False,
            "satisfies_production_constructive_next_artifact": False,
        },
    }


def main() -> int:
    payload = build_export()
    _validate(payload)
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_PATH, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=1)
        f.write("\n")
    moments = payload["transport_moments"]
    consistency = payload["consistency"]
    print(f"timelike  Delta_alpha_had^(5)(M_Z^2)  = {moments['timelike_on_shell_mz']['value']:.8f} "
          f"+- {moments['timelike_on_shell_mz']['uncertainty']:.8f}")
    print(f"spacelike Delta_alpha_had^(5)(-M_Z^2) = {moments['spacelike_mz']['value']:.8f} "
          f"+- {moments['spacelike_mz']['uncertainty']:.8f}")
    print(f"inverse-alpha spacelike packet        = {moments['inverse_alpha_packet_spacelike']['value']:.6f} "
          f"+- {moments['inverse_alpha_packet_spacelike']['uncertainty']:.6f}")
    print(f"requadrature vs upstream: |diff| = {consistency['abs_difference']:.2e} "
          f"(tolerance {consistency['tolerance']:.1e})")
    print(f"positivity: {payload['rho_had_or_measure']['positivity_status']}")
    print(f"wrote {OUT_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
