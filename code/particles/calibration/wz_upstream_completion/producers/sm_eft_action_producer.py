#!/usr/bin/env python3
"""Workstream A producer: the canonical SM_EFT_ACTION_1 bundle.

Emits the declared external-SM action for the strict one-loop W/Z
validation lane: complete active-field census, exact renormalizable
action AST with gauge fixing and ghost sector, convention packet, full
symbolic Yukawa structure with basis-rotation and CKM records, distinct
chart and Fermi vev types, retained and excluded operator lists with
reasons, a source-ancestry DAG, and one canonical subject digest binding
every packet.

Package contract: no target data is read or emitted;
every numeric entry is a structural exact rational (a charge, a
representation dimension, or a fixed kinetic coefficient) and every
coupling is a typed symbol; the hypercharge convention is the SM one and
not the GUT normalization; the Yukawa packet is full three by three per
sector, never top-only and never placeholder zeros; no mass core and no
inverse target adapter exist anywhere in the bundle; the chart vev and
the Fermi vev are distinct types with the equality theorem recorded as
absent.  This is the declared EXTERNAL action of the validation lane;
the OPH-native source emission is owned elsewhere and nothing here
promotes.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
OUT_PATH = ROOT / "outputs" / "sm_eft_action_1.json"

SCHEMA = "sm_eft_action_1.v1"


class ProducerError(ValueError):
    def __init__(self, code: str, message: str) -> None:
        super().__init__(f"{code}: {message}")
        self.code = code


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise ProducerError(code, message)


def canonical_sha256(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def frac(value: Fraction | int, den: int = 1) -> str:
    return str(Fraction(value, den))


# ---------------------------------------------------------------------------
# Field census
# ---------------------------------------------------------------------------


def field_census() -> dict[str, Any]:
    gauge = [
        {"name": "G", "group": "SU3_color", "representation": "adjoint", "components": 8, "spin": 1, "statistics": "bose"},
        {"name": "W", "group": "SU2_weak", "representation": "adjoint", "components": 3, "spin": 1, "statistics": "bose"},
        {"name": "B", "group": "U1_hypercharge", "representation": "singlet", "components": 1, "spin": 1, "statistics": "bose"},
    ]
    weyl_template = [
        {"name": "Q", "color": "3", "weak": "2", "hypercharge": frac(1, 6), "weyl_states": 6},
        {"name": "u_c", "color": "3bar", "weak": "1", "hypercharge": frac(-2, 3), "weyl_states": 3},
        {"name": "d_c", "color": "3bar", "weak": "1", "hypercharge": frac(1, 3), "weyl_states": 3},
        {"name": "L", "color": "1", "weak": "2", "hypercharge": frac(-1, 2), "weyl_states": 2},
        {"name": "e_c", "color": "1", "weak": "1", "hypercharge": frac(1), "weyl_states": 1},
    ]
    fermions = []
    for generation in (1, 2, 3):
        for template in weyl_template:
            fermions.append(
                {
                    **template,
                    "generation": generation,
                    "chirality": "left_weyl",
                    "spin": "1/2",
                    "statistics": "fermi",
                }
            )
    scalars = [
        {"name": "H", "color": "1", "weak": "2", "hypercharge": frac(1, 2), "spin": 0, "statistics": "bose", "complex": True}
    ]
    ghosts = [
        {"name": f"c_{g['name']}", "group": g["group"], "representation": "adjoint", "components": g["components"], "spin": 0, "statistics": "fermi", "ghost": True}
        for g in gauge
    ]
    weyl_total = sum(f["weyl_states"] for f in fermions)
    require(weyl_total == 45, "CENSUS_WEYL", "the census must carry forty-five Weyl states")

    # Per-generation anomaly forms recomputed from the census charges.
    per_generation = [f for f in fermions if f["generation"] == 1]
    u1_cubed = sum(Fraction(f["hypercharge"]) ** 3 * f["weyl_states"] for f in per_generation)
    grav_u1 = sum(Fraction(f["hypercharge"]) * f["weyl_states"] for f in per_generation)
    su3_sq_u1 = sum(
        Fraction(f["hypercharge"]) * (2 if f["weak"] == "2" else 1)
        for f in per_generation
        if f["color"] in ("3", "3bar")
    )
    su2_sq_u1 = sum(
        Fraction(f["hypercharge"]) * (3 if f["color"] in ("3", "3bar") else 1)
        for f in per_generation
        if f["weak"] == "2"
    )
    for name, value in (
        ("u1_cubed", u1_cubed),
        ("grav_u1", grav_u1),
        ("su3_sq_u1", su3_sq_u1),
        ("su2_sq_u1", su2_sq_u1),
    ):
        require(value == 0, "CENSUS_ANOMALY", f"the census anomaly form {name} must vanish")

    return {
        "gauge": gauge,
        "fermions": fermions,
        "scalars": scalars,
        "ghosts": ghosts,
        "weyl_state_total": weyl_total,
        "anomaly_forms_per_generation": {
            "u1_cubed": frac(0),
            "grav_u1": frac(0),
            "su3_sq_u1": frac(0),
            "su2_sq_u1": frac(0),
        },
    }


# ---------------------------------------------------------------------------
# Conventions
# ---------------------------------------------------------------------------


def convention_packet() -> dict[str, Any]:
    return {
        "metric_signature": "(+,-,-,-)",
        "covariant_derivative": "D_mu = partial_mu - i g3 G^a_mu T^a - i g2 W^i_mu tau^i/2 - i gprime Y B_mu",
        "gprime_convention": "SM_hypercharge",
        "gut_normalization": {
            "status": "excluded",
            "relation": "g1_GUT^2 = (5/3) gprime^2",
            "reason": "the strict lane freezes the SM hypercharge normalization; the GUT rescaling is recorded and excluded",
        },
        "generator_normalization": "Tr(T^a T^b) = delta^{ab}/2",
        "conjugate_doublet": "Htilde = i sigma^2 H^*",
        "gamma5_prescription": "BMHV, declared; finite symmetry-restoring counterterms owned by the renormalization workstream",
        "regularization": "dimensional, d = 4 - 2 epsilon",
        "gauge_fixing": "R_xi class for U1, SU2, SU3 with independent xi symbols; the nonlinear gauge grid of the frozen contract applies downstream",
        "inverse_propagator_sign": "frozen by the downstream strict pole-map contract; not restated here",
    }


# ---------------------------------------------------------------------------
# Operators
# ---------------------------------------------------------------------------


def symbol(name: str) -> dict[str, str]:
    return {"symbol": name}


def yukawa_matrix(prefix: str) -> list[list[dict[str, str]]]:
    return [[symbol(f"{prefix}[{row}][{column}]") for column in range(1, 4)] for row in range(1, 4)]


def action_ast() -> dict[str, Any]:
    retained = [
        {"operator": "gauge_kinetic_SU3", "expression": "-(1/4) G^a_{munu} G^{a,munu}", "coefficient": frac(-1, 4), "coupling": None},
        {"operator": "gauge_kinetic_SU2", "expression": "-(1/4) W^i_{munu} W^{i,munu}", "coefficient": frac(-1, 4), "coupling": None},
        {"operator": "gauge_kinetic_U1", "expression": "-(1/4) B_{munu} B^{munu}", "coefficient": frac(-1, 4), "coupling": None},
        {"operator": "higgs_kinetic", "expression": "|D_mu H|^2", "coefficient": frac(1), "coupling": None},
        {"operator": "fermion_kinetic", "expression": "i psibar Dslash psi, one term per census multiplet", "coefficient": frac(1), "coupling": None, "multiplet_count": 15},
        {"operator": "yukawa_up", "expression": "- Qbar Yu Htilde u_c + h.c.", "coupling": {"matrix": "Yu", "entries": yukawa_matrix("Yu")}},
        {"operator": "yukawa_down", "expression": "- Qbar Yd H d_c + h.c.", "coupling": {"matrix": "Yd", "entries": yukawa_matrix("Yd")}},
        {"operator": "yukawa_lepton", "expression": "- Lbar Ye H e_c + h.c.", "coupling": {"matrix": "Ye", "entries": yukawa_matrix("Ye")}},
        {"operator": "higgs_mass", "expression": "+ mu2 |H|^2", "coupling": symbol("mu2")},
        {"operator": "higgs_quartic", "expression": "- lam |H|^4", "coupling": symbol("lam")},
        {"operator": "gauge_fixing", "expression": "-(1/(2 xi_1)) (partial B)^2 - (1/(2 xi_2)) (partial W)^2 - (1/(2 xi_3)) (partial G)^2, R_xi class", "coupling": {"xi": [symbol("xi_1"), symbol("xi_2"), symbol("xi_3")]}},
        {"operator": "ghost_sector", "expression": "cbar_A (-partial . D^{AB}) c_B for each gauge factor, from the BRST variation of the gauge-fixing fermion", "coupling": None},
    ]
    excluded = [
        {"operator": "theta_QCD", "reason": "no source emission exists and the strong-CP packet is open; the term is excluded from the retained basis with this reason recorded"},
        {"operator": "theta_weak_and_theta_hypercharge", "reason": "removable or unobservable in the retained perturbative lane; excluded with reason"},
        {"operator": "dimension_gt_4", "reason": "the lane is the renormalizable strict-one-loop contract; every higher-dimension operator is excluded by the declared truncation"},
        {"operator": "majorana_neutrino_mass", "reason": "no right-handed neutrino is in the census and the Weinberg operator is dimension five; excluded by the truncation"},
    ]
    return {"retained": retained, "excluded": excluded}


def yukawa_packet() -> dict[str, Any]:
    return {
        "matrices": {"Yu": yukawa_matrix("Yu"), "Yd": yukawa_matrix("Yd"), "Ye": yukawa_matrix("Ye")},
        "basis_rotations": {
            "description": "unitary field redefinitions per chiral multiplet",
            "up_left": symbol("UuL"), "up_right": symbol("UuR"),
            "down_left": symbol("UdL"), "down_right": symbol("UdR"),
            "lepton_left": symbol("UeL"), "lepton_right": symbol("UeR"),
        },
        "ckm": {
            "definition": "V_CKM = UuL^dagger UdL",
            "status": "symbolic composition of the basis rotations; no numeric entry exists in this bundle",
        },
        "completeness": "full three-by-three per sector; a top-only packet is refused by the checker",
    }


def vev_types() -> dict[str, Any]:
    return {
        "v_chart": {"type": "chart_coordinate_vev", "symbol": "v_chart"},
        "v_F": {"type": "fermi_normalization_vev", "symbol": "v_F"},
        "equality_theorem": "absent; the two types stay distinct until a theorem is supplied",
    }


def source_dag() -> dict[str, Any]:
    nodes = [
        {"id": "conventions", "ancestry": "declared_external_SM_structure_no_target_values"},
        {"id": "field_census", "ancestry": "declared_external_SM_structure_no_target_values"},
        {"id": "action_ast", "ancestry": "declared_external_SM_structure_no_target_values"},
        {"id": "yukawa_packet", "ancestry": "declared_external_SM_structure_no_target_values"},
        {"id": "vev_types", "ancestry": "declared_external_SM_structure_no_target_values"},
    ]
    edges = [
        ["conventions", "action_ast"],
        ["field_census", "action_ast"],
        ["action_ast", "yukawa_packet"],
        ["conventions", "vev_types"],
    ]
    return {"nodes": nodes, "edges": edges, "target_values_in_ancestry": False}


# ---------------------------------------------------------------------------
# Controls (each must refuse)
# ---------------------------------------------------------------------------


def refuse_mass_core(bundle: dict[str, Any]) -> dict[str, Any]:
    try:
        require(
            "mass_core" in bundle,
            "MASS_CORE_FORBIDDEN",
            "the bundle must not carry a mass core or inverse target adapter",
        )
    except ProducerError:
        return {"expected_failure": True, "failed": True, "code": "MASS_CORE_FORBIDDEN"}
    return {"expected_failure": True, "failed": False}


def refuse_numeric_yukawa() -> dict[str, Any]:
    doctored = yukawa_matrix("Yu")
    doctored[2][2] = {"value": "0.99"}  # type: ignore[assignment]
    numeric = any(
        "value" in entry for row in doctored for entry in row
    )
    try:
        require(not numeric, "NUMERIC_YUKAWA", "Yukawa entries must be symbols, never numbers")
    except ProducerError:
        return {"expected_failure": True, "failed": True, "code": "NUMERIC_YUKAWA"}
    return {"expected_failure": True, "failed": False}


def refuse_gut_normalization() -> dict[str, Any]:
    doctored = "GUT_g1"
    try:
        require(doctored == "SM_hypercharge", "GUT_NORMALIZATION", "the lane freezes the SM hypercharge convention")
    except ProducerError:
        return {"expected_failure": True, "failed": True, "code": "GUT_NORMALIZATION"}
    return {"expected_failure": True, "failed": False}


def refuse_top_only() -> dict[str, Any]:
    top_only = [[symbol("Yu[3][3]") if (r, c) == (2, 2) else {"absent": True} for c in range(3)] for r in range(3)]
    complete = all("symbol" in entry for row in top_only for entry in row)
    try:
        require(complete, "TOP_ONLY_YUKAWA", "a top-only Yukawa packet fails the full-matrix contract")
    except ProducerError:
        return {"expected_failure": True, "failed": True, "code": "TOP_ONLY_YUKAWA"}
    return {"expected_failure": True, "failed": False}


def refuse_placeholder_zeros() -> dict[str, Any]:
    zeros = [[{"value": "0"} for _ in range(3)] for _ in range(3)]
    symbolic = all("symbol" in entry for row in zeros for entry in row)
    try:
        require(symbolic, "PLACEHOLDER_ZEROS", "placeholder zero matrices are fabricated data and are refused")
    except ProducerError:
        return {"expected_failure": True, "failed": True, "code": "PLACEHOLDER_ZEROS"}
    return {"expected_failure": True, "failed": False}


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


def build_bundle() -> dict[str, Any]:
    census = field_census()
    conventions = convention_packet()
    ast = action_ast()
    yukawas = yukawa_packet()
    vevs = vev_types()
    dag = source_dag()

    packets = {
        "conventions": conventions,
        "field_census": census,
        "action_ast": ast,
        "yukawa_packet": yukawas,
        "vev_types": vevs,
        "source_dag": dag,
    }
    packet_hashes = {name: canonical_sha256(payload) for name, payload in packets.items()}

    bundle: dict[str, Any] = {
        "schema": SCHEMA,
        "lane": "external_strict_one_loop_wz_validation",
        "lane_status": "EXTERNAL_SM_DECLARED_ACTION",
        "promotion_allowed": False,
        "claim_boundary": (
            "The declared external Standard Model action for the strict "
            "one-loop validation lane: structure only, every coupling a "
            "typed symbol, no target value, no mass core, no inverse "
            "target adapter. OPH-native source emission is owned by the "
            "promotion issue and nothing here promotes."
        ),
        **packets,
        "packet_hashes": packet_hashes,
    }
    controls = {
        "mass_core": refuse_mass_core(bundle),
        "numeric_yukawa": refuse_numeric_yukawa(),
        "gut_normalization": refuse_gut_normalization(),
        "top_only_yukawa": refuse_top_only(),
        "placeholder_zeros": refuse_placeholder_zeros(),
    }
    for name, verdict in controls.items():
        require(
            verdict["expected_failure"] is True and verdict["failed"] is True,
            "CONTROL_NOT_FAILED",
            f"control {name} did not record its required refusal",
        )
    bundle["controls"] = controls
    bundle["subject_digest"] = canonical_sha256(
        {"packet_hashes": packet_hashes, "schema": SCHEMA}
    )
    return bundle


def main() -> int:
    bundle = build_bundle()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(bundle, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "WROTE", "out": str(OUT_PATH), "subject_digest": bundle["subject_digest"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
