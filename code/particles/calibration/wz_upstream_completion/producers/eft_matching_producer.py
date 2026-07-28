#!/usr/bin/env python3
"""Workstream B producer: the EFT_MATCHING_1 interval packet.

Emits the renormalization-group matching packet for the strict one-loop
W/Z validation lane, bound to the Workstream A subject digest.  The
scale axis carries one pure Standard Model interval; the packet lists
the active fields with representations, spin, and statistics read from
the action bundle census, and derives every one-loop gauge beta
coefficient FROM that census through exact group invariants:

    b = -(11/3) C2(adjoint) + (2/3) sum_Weyl T(R) + (1/3) sum_cx_scalar T(R)

with T(fundamental) = 1/2, C2(SU3) = 3, C2(SU2) = 2, and the abelian
index Y^2 in the SM hypercharge normalization.  The derived vector is
(b1, b2, b3) = (41/6, -19/6, -7).  The GUT-normalized value 41/10 and
every MSSM coefficient are refused by controls.  Yukawa and quartic
one-loop beta functions enter as declared imported structure with their
census-derived group factors recorded separately from the imported
kernel shapes, so nothing numeric is fabricated.

Schemes are frozen (MSbar for the running, the frozen FJ coordinate for
the output vector), the DRbar-to-MSbar finite shifts are recorded as
imported finite maps with their known one-loop coefficients labeled as
imports, decoupling maps are recorded empty for the pure interval with
the reason stated, the interval Jacobian is the symbolic one-loop
linearization, and the truncation remainder is a declared symbol, never
a fabricated number.  The output names one canonical SM_MSbar_FJ(Q)
vector as a typed symbol list.  Nothing promotes.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
ACTION_PATH = ROOT / "outputs" / "sm_eft_action_1.json"
OUT_PATH = ROOT / "outputs" / "eft_matching_1.json"

SCHEMA = "eft_matching_1.v1"


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


def frac(value: Fraction) -> str:
    return str(value)


# ---------------------------------------------------------------------------
# Census-derived group invariants
# ---------------------------------------------------------------------------


DYNKIN_FUNDAMENTAL = Fraction(1, 2)
CASIMIR_ADJOINT = {"SU3_color": Fraction(3), "SU2_weak": Fraction(2)}


def load_action_bundle() -> dict[str, Any]:
    require(ACTION_PATH.is_file(), "ACTION_MISSING", "run the Workstream A producer first")
    return json.loads(ACTION_PATH.read_text(encoding="utf-8"))


def color_dynkin(label: str) -> Fraction:
    if label in ("3", "3bar"):
        return DYNKIN_FUNDAMENTAL
    require(label == "1", "REP_LABEL", f"unknown color representation {label}")
    return Fraction(0)


def weak_dynkin(label: str) -> Fraction:
    if label == "2":
        return DYNKIN_FUNDAMENTAL
    require(label == "1", "REP_LABEL", f"unknown weak representation {label}")
    return Fraction(0)


def color_dimension(label: str) -> int:
    return 3 if label in ("3", "3bar") else 1


def weak_dimension(label: str) -> int:
    return 2 if label == "2" else 1


def gauge_beta_coefficients(census: dict[str, Any]) -> dict[str, Any]:
    fermions = census["fermions"]
    scalars = census["scalars"]

    # SU(3): Weyl Dynkin sum over color multiplets (weak dimension counts copies).
    su3_weyl = sum(
        color_dynkin(f["color"]) * weak_dimension(f["weak"]) for f in fermions
    )
    su3_scalar = sum(
        color_dynkin(s["color"]) * weak_dimension(s["weak"]) for s in scalars
    )
    b3 = -Fraction(11, 3) * CASIMIR_ADJOINT["SU3_color"] + Fraction(2, 3) * su3_weyl + Fraction(1, 3) * su3_scalar

    # SU(2): Weyl Dynkin sum over weak doublets (color dimension counts copies).
    su2_weyl = sum(
        weak_dynkin(f["weak"]) * color_dimension(f["color"]) for f in fermions
    )
    su2_scalar = sum(
        weak_dynkin(s["weak"]) * color_dimension(s["color"]) for s in scalars
    )
    b2 = -Fraction(11, 3) * CASIMIR_ADJOINT["SU2_weak"] + Fraction(2, 3) * su2_weyl + Fraction(1, 3) * su2_scalar

    # U(1) in the SM hypercharge normalization: index Y^2 times multiplicity.
    u1_weyl = sum(
        Fraction(f["hypercharge"]) ** 2 * color_dimension(f["color"]) * weak_dimension(f["weak"])
        for f in fermions
    )
    u1_scalar = sum(
        Fraction(s["hypercharge"]) ** 2 * color_dimension(s["color"]) * weak_dimension(s["weak"])
        for s in scalars
    )
    b1 = Fraction(2, 3) * u1_weyl + Fraction(1, 3) * u1_scalar

    require(b3 == Fraction(-7), "BETA_SU3", f"census-derived b3 must be -7, got {b3}")
    require(b2 == Fraction(-19, 6), "BETA_SU2", f"census-derived b2 must be -19/6, got {b2}")
    require(b1 == Fraction(41, 6), "BETA_U1", f"census-derived b1 must be 41/6 in the SM normalization, got {b1}")

    return {
        "convention": "d g_i / d ln Q = b_i g_i^3 / (16 pi^2), SM hypercharge normalization for b1",
        "derivation": {
            "formula": "b = -(11/3) C2(adjoint) + (2/3) sum_Weyl T(R) + (1/3) sum_complex_scalar T(R)",
            "su3_weyl_dynkin_sum": frac(su3_weyl),
            "su3_scalar_dynkin_sum": frac(su3_scalar),
            "su2_weyl_dynkin_sum": frac(su2_weyl),
            "su2_scalar_dynkin_sum": frac(su2_scalar),
            "u1_weyl_index_sum": frac(u1_weyl),
            "u1_scalar_index_sum": frac(u1_scalar),
        },
        "coefficients": {"b1": frac(b1), "b2": frac(b2), "b3": frac(b3)},
        "gut_normalized_b1": {
            "value": frac(Fraction(41, 10)),
            "status": "excluded",
            "reason": "the lane freezes the SM hypercharge normalization; the GUT value is recorded for refusal only",
        },
    }


def yukawa_quartic_structure(census: dict[str, Any]) -> dict[str, Any]:
    """Census-derived group factors for the imported one-loop Yukawa and
    quartic kernels; the kernel shapes are declared imports and the
    numeric factors here are exact census consequences."""

    nc = 3
    return {
        "imported_kernel_shapes": [
            "one-loop Yukawa beta: (Nc + 3/2) y^3 - like terms minus gauge Casimir subtractions",
            "one-loop quartic beta: quartic, Yukawa-trace, and gauge-quartic terms",
        ],
        "census_group_factors": {
            "color_number": nc,
            "yukawa_trace_multiplicity": "Tr over the full three-by-three sector matrices of the action bundle",
            "gauge_casimirs": {
                "C2_fundamental_SU3": frac(Fraction(4, 3)),
                "C2_fundamental_SU2": frac(Fraction(3, 4)),
            },
        },
        "status": "imported_structure_with_census_factors; no numeric coupling value exists in this packet",
    }


# ---------------------------------------------------------------------------
# Interval, schemes, maps
# ---------------------------------------------------------------------------


def interval_packet(census: dict[str, Any], betas: dict[str, Any]) -> dict[str, Any]:
    active = {
        "gauge": [g["name"] for g in census["gauge"]],
        "fermion_multiplets": len(census["fermions"]),
        "scalars": [s["name"] for s in census["scalars"]],
        "ghosts": [g["name"] for g in census["ghosts"]],
    }
    return {
        "intervals": [
            {
                "id": "pure_sm_interval",
                "range": {"lower": {"symbol": "Q_low"}, "upper": {"symbol": "Q_high"}},
                "active_fields": active,
                "operator_basis": "the retained renormalizable basis of the action bundle",
                "beta_coefficients": betas["coefficients"],
                "decoupling_maps": {
                    "maps": [],
                    "reason": "no threshold crosses a pure Standard Model interval; the empty list is a recorded statement, not an omission",
                },
            }
        ],
        "scheme_freeze": {
            "running": "MSbar",
            "output_coordinate": "the frozen FJ coordinate of the downstream strict pole-map contract",
            "drbar_to_msbar_finite_maps": {
                "status": "imported_finite_maps",
                "gauge_shift": "g_DRbar^2 = g_MSbar^2 (1 + C2(adjoint) g^2/(24 pi^2) x epsilon-scheme term), recorded as imported structure with symbolic coupling",
                "reason": "the finite scheme shifts are standard imported structure; no numeric value is fabricated",
            },
        },
        "jacobian": {
            "definition": "the one-loop linearization d g_i(Q2) / d g_j(Q1) of the interval flow",
            "form": "diagonal at one loop in the gauge sector: delta_ij (1 + 3 b_i g_i^2 t / (16 pi^2)) with t = ln(Q2/Q1), symbolic",
        },
        "remainder": {
            "definition": "the two-loop truncation term of the interval flow",
            "form": {"symbol": "R2[g, t]"},
            "status": "declared deterministic remainder symbol; it is carried, never estimated without inputs",
        },
        "output_vector": {
            "name": "SM_MSbar_FJ(Q)",
            "components": [
                {"symbol": "g1(Q)"}, {"symbol": "g2(Q)"}, {"symbol": "g3(Q)"},
                {"symbol": "Yu(Q)"}, {"symbol": "Yd(Q)"}, {"symbol": "Ye(Q)"},
                {"symbol": "lam(Q)"}, {"symbol": "mu2(Q)"},
                {"symbol": "v_chart(Q)"},
            ],
            "note": "one canonical typed vector; every component is a symbol and the chart vev keeps its type from the action bundle",
        },
    }


# ---------------------------------------------------------------------------
# Controls (each must refuse)
# ---------------------------------------------------------------------------


def refuse_mssm_beta() -> dict[str, Any]:
    mssm_b3 = Fraction(-3)
    try:
        require(mssm_b3 == Fraction(-7), "MSSM_BETA", "an MSSM coefficient in a pure-SM interval is refused")
    except ProducerError:
        return {"expected_failure": True, "failed": True, "code": "MSSM_BETA"}
    return {"expected_failure": True, "failed": False}


def refuse_gut_b1() -> dict[str, Any]:
    gut = Fraction(41, 10)
    try:
        require(gut == Fraction(41, 6), "GUT_B1", "the GUT-normalized b1 is refused in the SM-normalization lane")
    except ProducerError:
        return {"expected_failure": True, "failed": True, "code": "GUT_B1"}
    return {"expected_failure": True, "failed": False}


def refuse_scalar_as_weyl(census: dict[str, Any]) -> dict[str, Any]:
    """Counting the Higgs with the Weyl weight must change b2 and be refused."""

    su2_weyl = sum(
        weak_dynkin(f["weak"]) * color_dimension(f["color"]) for f in census["fermions"]
    )
    wrong_b2 = -Fraction(11, 3) * CASIMIR_ADJOINT["SU2_weak"] + Fraction(2, 3) * (su2_weyl + DYNKIN_FUNDAMENTAL)
    try:
        require(wrong_b2 == Fraction(-19, 6), "SCALAR_WEIGHT", "a scalar counted at the Weyl weight is refused")
    except ProducerError:
        return {"expected_failure": True, "failed": True, "code": "SCALAR_WEIGHT"}
    return {"expected_failure": True, "failed": False}


def refuse_missing_generation(census: dict[str, Any]) -> dict[str, Any]:
    reduced = [f for f in census["fermions"] if f["generation"] != 3]
    u1 = Fraction(2, 3) * sum(
        Fraction(f["hypercharge"]) ** 2 * color_dimension(f["color"]) * weak_dimension(f["weak"])
        for f in reduced
    ) + Fraction(1, 6)
    try:
        require(u1 == Fraction(41, 6), "GENERATION_COUNT", "a two-generation census is refused")
    except ProducerError:
        return {"expected_failure": True, "failed": True, "code": "GENERATION_COUNT"}
    return {"expected_failure": True, "failed": False}


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


def build_bundle() -> dict[str, Any]:
    action = load_action_bundle()
    census = action["field_census"]
    betas = gauge_beta_coefficients(census)
    yukawa = yukawa_quartic_structure(census)
    intervals = interval_packet(census, betas)

    controls = {
        "mssm_beta": refuse_mssm_beta(),
        "gut_b1": refuse_gut_b1(),
        "scalar_as_weyl": refuse_scalar_as_weyl(census),
        "missing_generation": refuse_missing_generation(census),
    }
    for name, verdict in controls.items():
        require(
            verdict["expected_failure"] is True and verdict["failed"] is True,
            "CONTROL_NOT_FAILED",
            f"control {name} did not record its required refusal",
        )

    bundle: dict[str, Any] = {
        "schema": SCHEMA,
        "lane": "external_strict_one_loop_wz_validation",
        "promotion_allowed": False,
        "claim_boundary": (
            "The interval EFT/RG packet of the validation lane: one pure "
            "Standard Model interval with census-derived one-loop gauge "
            "beta coefficients, frozen schemes, imported finite maps "
            "recorded as imports, a symbolic interval Jacobian, a declared "
            "truncation remainder, and one canonical typed output vector. "
            "No numeric coupling value exists anywhere in this packet, and "
            "nothing promotes."
        ),
        "action_subject_digest": action["subject_digest"],
        "gauge_betas": betas,
        "yukawa_quartic_structure": yukawa,
        "matching": intervals,
        "controls": controls,
    }
    bundle["subject_digest"] = canonical_sha256(
        {
            "schema": SCHEMA,
            "action_subject_digest": action["subject_digest"],
            "gauge_betas": betas,
            "matching": intervals,
        }
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
