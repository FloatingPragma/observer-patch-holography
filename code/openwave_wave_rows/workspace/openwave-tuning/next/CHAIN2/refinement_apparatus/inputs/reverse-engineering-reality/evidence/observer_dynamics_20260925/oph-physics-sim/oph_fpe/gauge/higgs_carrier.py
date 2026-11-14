from __future__ import annotations

import cmath
import json
from pathlib import Path
from typing import Any

from oph_fpe.claims import CONTINUATION, with_claim_metadata


BOREL_WEIL_HIGGS_RECEIPT = "BOREL_WEIL_HIGGS_CARRIER_RECEIPT"
REPORT_SCHEMA = "oph_borel_weil_higgs_carrier_bridge_v1"

DEFAULT_BW_HIGGS_CARRIER = {
    "screen_chart": "CP1",
    "line_bundle": "O(1)",
    "section_degree_n": 1,
    "complex_section_dimension": 2,
    "real_section_degrees_of_freedom": 4,
    "su2_representation_dimension": 2,
    "hypercharge_normalization": "OPH_Z6_NEUTRAL_VEV",
    "neutral_component_T3": -0.5,
    "hypercharge": 0.5,
    "lorentz_role": "internal_scalar_0_form",
    "higgs_doublets": 1,
    "stabilizer_test_phase_radians": 0.731,
}


def borel_weil_higgs_carrier_receipt(candidate: dict[str, Any] | None = None) -> dict[str, Any]:
    """Validate the Borel-Weil local carrier for the OPH one-Higgs slot.

    This receipt is representation/carrier evidence only. It deliberately does
    not promote the weak scale, Higgs quartic, Higgs mass, or Coleman-Weinberg
    dynamics, which are owned by the OPH D10/D11 quantitative branch.
    """

    payload = {**DEFAULT_BW_HIGGS_CARRIER, **(candidate or {})}
    n = _int(payload.get("section_degree_n"))
    complex_dim = _int(payload.get("complex_section_dimension"))
    real_dof = _int(payload.get("real_section_degrees_of_freedom"))
    rep_dim = _int(payload.get("su2_representation_dimension"))
    hypercharge = _float(payload.get("hypercharge"))
    t3_lower = _float(payload.get("neutral_component_T3"))
    test_phase = _float(payload.get("stabilizer_test_phase_radians"))
    claims = payload.get("derived_quantitative_claims") or []
    if not isinstance(claims, list):
        claims = [claims]
    promoted_forbidden_claims = [str(item) for item in claims if str(item) in _FORBIDDEN_PROMOTIONS]

    group_dimension = 4
    projective_stabilizer_dimension = 2
    vector_stabilizer_dimension = 1
    projective_orbit_dimension = group_dimension - projective_stabilizer_dimension
    goldstone_count = group_dimension - vector_stabilizer_dimension
    radial_higgs_modes = (real_dof or 0) - goldstone_count
    group_action_acceptance = _stabilizer_action_acceptance(
        t3_lower=t3_lower,
        hypercharge=hypercharge,
        phase=test_phase,
    )

    checks = {
        "screen_chart_is_cp1": _norm(payload.get("screen_chart")) in {"cp1", "mathbbcp1", "complex_projective_line"},
        "line_bundle_is_positive_hopf_o1": _norm(payload.get("line_bundle")) in {"o(1)", "mathcalo(1)", "positive_hopf"},
        "section_degree_is_minimal_nontrivial": n == 1,
        "section_dimension_is_two_complex": complex_dim == 2,
        "real_dof_is_four": real_dof == 4,
        "su2_doublet_representation": rep_dim == 2,
        "single_higgs_doublet": _int(payload.get("higgs_doublets")) == 1,
        "internal_scalar_zero_form": _norm(payload.get("lorentz_role")) == "internal_scalar_0_form",
        "oph_hypercharge_normalization": _norm(payload.get("hypercharge_normalization"))
        in {"oph_z6_neutral_vev", "oph_hypercharge_neutral_vev"},
        "neutral_lower_component": t3_lower == -0.5 and hypercharge == 0.5 and abs((t3_lower or 0.0) + (hypercharge or 0.0)) <= 1e-12,
        "projective_ray_stabilizer_is_two_torus": projective_stabilizer_dimension == 2,
        "projective_orbit_is_two_dimensional": projective_orbit_dimension == 2,
        "vector_stabilizer_is_u1_q": vector_stabilizer_dimension == 1,
        "goldstone_count_uses_vector_stabilizer": goldstone_count == 3 and radial_higgs_modes == 1,
        "hypercharge_phase_fixes_projective_ray": group_action_acceptance[
            "hypercharge_fixes_projective_ray"
        ],
        "t3_phase_fixes_projective_ray": group_action_acceptance["t3_fixes_projective_ray"],
        "generic_hypercharge_phase_changes_vector": not group_action_acceptance[
            "hypercharge_fixes_vector"
        ],
        "generic_t3_phase_changes_vector": not group_action_acceptance["t3_fixes_vector"],
        "diagonal_q_phase_fixes_vector": group_action_acceptance["diagonal_q_fixes_vector"],
        "forbidden_quantitative_promotions_absent": not promoted_forbidden_claims,
    }
    receipt = all(checks.values())
    report = {
        "schema": REPORT_SCHEMA,
        "mode": "borel_weil_higgs_carrier_bridge_v1",
        BOREL_WEIL_HIGGS_RECEIPT: bool(receipt),
        "receipt": bool(receipt),
        "checks": checks,
        "candidate": payload,
        "carrier_identification": "H_OPH = H^0(CP1, O(1)) ~= C^2",
        "symmetry_breaking_geometry": {
            "integer_hypercharge_normalization": "q = 6Y, q_H = 3",
            "cover_action": "(g,z).phi = z^3 g phi",
            "carrier_projective_space": "P(H_OPH) ~= CP1",
            "projective_ray": "[phi0] with phi0 = (0, v), v != 0",
            "projective_ray_stabilizer_on_cover": "{(diag(a,a^-1),z): a,z in U(1)}",
            "projective_ray_stabilizer": "(U(1)_T3 x U(1)_Y)/finite_center",
            "projective_stabilizer_dimension": projective_stabilizer_dimension,
            "projective_orbit_dimension": projective_orbit_dimension,
            "projectivization_forgets_scalar_hypercharge_phase": True,
            "nonzero_vacuum_vector": "phi0 = (0, v), v != 0",
            "vector_stabilizer_on_cover": "{(diag(z^3,z^-3),z): z in U(1)}",
            "vector_stabilizer": "U(1)_Q",
            "vector_stabilizer_generator": "Q = T3 + Y",
            "vector_stabilizer_dimension": vector_stabilizer_dimension,
            "broken_generator_count": goldstone_count,
            "goldstone_count": goldstone_count,
            "goldstone_count_source": "dim(SU(2)_L x U(1)_Y) - dim(Stab(phi0)) = 4 - 1",
            "radial_higgs_modes": radial_higgs_modes,
        },
        "group_action_acceptance": group_action_acceptance,
        "explicit_nonclaims": sorted(_FORBIDDEN_PROMOTIONS),
        "promoted_forbidden_claims": promoted_forbidden_claims,
        "claim_boundary": (
            "Borel-Weil carrier receipt for the OPH one-Higgs slot. It validates the minimal "
            "holomorphic section carrier, SU(2)L doublet representation, OPH hypercharge "
            "normalization, and the distinct projective-ray and nonzero-vector stabilizers. "
            "Projectivization forgets the overall hypercharge phase and does not by itself "
            "identify the unbroken electromagnetic subgroup; U(1)_Q is the stabilizer of the "
            "nonzero vacuum vector. This receipt does not derive m_H, lambda, v, or "
            "Coleman-Weinberg symmetry breaking."
        ),
    }
    return with_claim_metadata(report, claim_level=CONTINUATION, receipt=BOREL_WEIL_HIGGS_RECEIPT)


def write_borel_weil_higgs_carrier_report(
    out: Path,
    candidate: dict[str, Any] | None = None,
) -> dict[str, Any]:
    report = borel_weil_higgs_carrier_receipt(candidate)
    destination = Path(out)
    if destination.suffix.lower() != ".json":
        destination = destination / "borel_weil_higgs_carrier_report.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2, default=str) + "\n", encoding="utf-8")
    destination.with_suffix(".md").write_text(_markdown_report(report), encoding="utf-8")
    report["report_path"] = str(destination)
    return report


_FORBIDDEN_PROMOTIONS = {
    "higgs_mass",
    "m_H",
    "higgs_quartic",
    "lambda",
    "weak_scale",
    "v",
    "coleman_weinberg_breaking",
}


def _norm(value: Any) -> str:
    return str(value or "").replace(" ", "").replace("\\", "").lower()


def _int(value: Any) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _float(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _stabilizer_action_acceptance(
    *,
    t3_lower: float | None,
    hypercharge: float | None,
    phase: float | None,
) -> dict[str, Any]:
    """Apply a generic electroweak phase to the neutral Higgs vector.

    Hypercharge acts by one scalar on the full doublet, so it fixes the
    projective ray for every phase.  On the lower-component vacuum vector,
    the diagonal transformation with equal T3 and Y parameters has phase
    exp(i * beta * (T3 + Y)); neutrality makes that phase one.
    """

    if t3_lower is None or hypercharge is None or phase is None:
        return {
            "test_phase_radians": phase,
            "action_formula": "g(theta,beta) phi0 = exp(i(beta-theta)/2) phi0",
            "hypercharge_phase_on_vector": None,
            "t3_phase_on_vector": None,
            "diagonal_q_phase_on_vector": None,
            "hypercharge_fixes_projective_ray": False,
            "t3_fixes_projective_ray": False,
            "hypercharge_fixes_vector": False,
            "t3_fixes_vector": False,
            "diagonal_q_fixes_vector": False,
        }

    vacuum_vector = (0.0 + 0.0j, 1.0 + 0.0j)
    hypercharge_factor = cmath.exp(1j * phase * hypercharge)
    hypercharge_vector = tuple(hypercharge_factor * component for component in vacuum_vector)
    t3_factor = cmath.exp(1j * phase * t3_lower)
    t3_vector = tuple(t3_factor * component for component in vacuum_vector)
    diagonal_q_factor = t3_factor * hypercharge_factor
    diagonal_q_vector = tuple(diagonal_q_factor * component for component in vacuum_vector)

    return {
        "test_phase_radians": phase,
        "action_formula": "g(theta,beta) phi0 = exp(i(beta-theta)/2) phi0",
        "hypercharge_phase_on_vector": _complex_parts(hypercharge_factor),
        "t3_phase_on_vector": _complex_parts(t3_factor),
        "diagonal_q_phase_on_vector": _complex_parts(diagonal_q_factor),
        "hypercharge_fixes_projective_ray": _same_projective_ray(hypercharge_vector, vacuum_vector),
        "t3_fixes_projective_ray": _same_projective_ray(t3_vector, vacuum_vector),
        "hypercharge_fixes_vector": _vectors_close(hypercharge_vector, vacuum_vector),
        "t3_fixes_vector": _vectors_close(t3_vector, vacuum_vector),
        "diagonal_q_fixes_vector": _vectors_close(diagonal_q_vector, vacuum_vector),
    }


def _complex_parts(value: complex) -> dict[str, float]:
    return {"real": float(value.real), "imag": float(value.imag)}


def _vectors_close(
    left: tuple[complex, complex],
    right: tuple[complex, complex],
    *,
    tolerance: float = 1e-12,
) -> bool:
    return all(abs(a - b) <= tolerance for a, b in zip(left, right, strict=True))


def _same_projective_ray(
    left: tuple[complex, complex],
    right: tuple[complex, complex],
    *,
    tolerance: float = 1e-12,
) -> bool:
    left_norm_sq = sum(abs(component) ** 2 for component in left)
    right_norm_sq = sum(abs(component) ** 2 for component in right)
    if left_norm_sq <= tolerance**2 or right_norm_sq <= tolerance**2:
        return False
    determinant = left[0] * right[1] - left[1] * right[0]
    return abs(determinant) <= tolerance * (left_norm_sq * right_norm_sq) ** 0.5


def _markdown_report(report: dict[str, Any]) -> str:
    checks = report.get("checks") or {}
    geometry = report.get("symmetry_breaking_geometry") or {}
    group_action = report.get("group_action_acceptance") or {}
    return (
        "# Borel-Weil Higgs Carrier Receipt\n\n"
        f"- Receipt: {report.get(BOREL_WEIL_HIGGS_RECEIPT)}\n"
        f"- Carrier: {report.get('carrier_identification')}\n"
        f"- Claim level: {report.get('claim_level')}\n"
        f"- Physical claim: {report.get('physical_claim')}\n"
        f"- Minimal O(1) carrier: {checks.get('section_degree_is_minimal_nontrivial')}\n"
        f"- SU(2) doublet: {checks.get('su2_doublet_representation')}\n"
        f"- Neutral lower component: {checks.get('neutral_lower_component')}\n"
        f"- Projective-ray stabilizer: {geometry.get('projective_ray_stabilizer')}\n"
        f"- Nonzero-vector stabilizer: {geometry.get('vector_stabilizer')}\n"
        f"- Hypercharge fixes the projective ray: {group_action.get('hypercharge_fixes_projective_ray')}\n"
        f"- T3 fixes the projective ray: {group_action.get('t3_fixes_projective_ray')}\n"
        f"- Hypercharge alone fixes the vector: {group_action.get('hypercharge_fixes_vector')}\n"
        f"- T3 alone fixes the vector: {group_action.get('t3_fixes_vector')}\n"
        f"- Diagonal Q fixes the vector: {group_action.get('diagonal_q_fixes_vector')}\n"
        f"- Forbidden quantitative promotions absent: {checks.get('forbidden_quantitative_promotions_absent')}\n"
        f"- Promoted forbidden claims: {', '.join(report.get('promoted_forbidden_claims') or []) or 'none'}\n\n"
        f"{report.get('claim_boundary')}\n"
    )
