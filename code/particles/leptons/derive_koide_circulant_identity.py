#!/usr/bin/env python3
"""Emit the exact C3-circulant Koide identity and its claim boundary.

For a positive Hermitian circulant square-root-mass carrier

    C = a I + b R + conjugate(b) R^2,

the eigenvalues are ``a + 2 |b| cos(delta + 2 pi k / 3)``.  In a chamber
where all three eigenvalues are nonnegative, direct trigonometric summation
gives

    Q = 1/3 + (2/3) (|b|/a)^2.

Thus ``Q = 2/3`` if and only if ``|b|/a = 1/sqrt(2)``.  The identity does
not determine the phase and therefore does not determine the two independent
mass ratios.  Outside the positive chamber, physical square roots are
absolute eigenvalues and the signed-eigenvalue formula is not physical Q.

The artifact also evaluates the historically target-informed MCPR coordinate
and a compare-only PDG central-mass coordinate.  Neither comparison promotes a
source-only charged-lepton prediction.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from decimal import Decimal, localcontext
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
MCPR = (
    ROOT
    / "particles"
    / "runs"
    / "leptons"
    / "charged_mcpr_completion_conditional.json"
)
FINITE_GNS = (
    ROOT
    / "particles"
    / "runs"
    / "leptons"
    / "charged_koide_orientation_isometry.json"
)
FACE_CARRIER = (
    ROOT
    / "particles"
    / "runs"
    / "leptons"
    / "charged_icosahedral_face_carrier_frontier.json"
)
PDG_REFERENCES = ROOT / "particles" / "data" / "particle_reference_values.json"
DEFAULT_OUT = (
    ROOT
    / "particles"
    / "runs"
    / "leptons"
    / "koide_circulant_identity.json"
)
EXPECTED_PDG_2026_MASSES_MEV = {
    "electron": "0.51099895069",
    "muon": "105.6583755",
    "tau": "1776.93",
}
EXPECTED_PDG_2026_Q = Decimal(
    "0.6666644634026367029382786248892424584933042461962970542137359728464935"
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def decimal_text(value: Decimal) -> str:
    return format(value, "f")


def koide_from_modulus_ratio(ratio: Decimal) -> Decimal:
    return (Decimal(1) + Decimal(2) * ratio * ratio) / Decimal(3)


def modulus_ratio_from_koide(q_value: Decimal) -> Decimal:
    with localcontext() as context:
        context.prec = 70
        radicand = (Decimal(3) * q_value - Decimal(1)) / Decimal(2)
        if radicand < 0:
            raise ValueError("Q must be at least 1/3")
        return radicand.sqrt()


def roots(a_value: float, abs_b: float, phase: float) -> tuple[float, ...]:
    return tuple(
        a_value
        + 2.0
        * abs_b
        * math.cos(phase + 2.0 * math.pi * index / 3.0)
        for index in range(3)
    )


def physical_koide(root_values: tuple[float, ...]) -> float:
    masses = [value * value for value in root_values]
    return sum(masses) / sum(abs(value) for value in root_values) ** 2


def measured_central_q(references: dict[str, Any]) -> tuple[Decimal, dict[str, str]]:
    entries = references["entries"]
    masses = {
        name: Decimal(str(entries[name]["raw_value"]))
        for name in ("electron", "muon", "tau")
    }
    with localcontext() as context:
        context.prec = 70
        denominator = sum(value.sqrt() for value in masses.values()) ** 2
        q_value = sum(masses.values()) / denominator
    return q_value, {name: decimal_text(value) for name, value in masses.items()}


def build_artifact() -> dict[str, Any]:
    mcpr = json.loads(MCPR.read_text(encoding="utf-8"))
    finite_gns = json.loads(FINITE_GNS.read_text(encoding="utf-8"))
    face = json.loads(FACE_CARRIER.read_text(encoding="utf-8"))
    references = json.loads(PDG_REFERENCES.read_text(encoding="utf-8"))

    with localcontext() as context:
        context.prec = 70
        one = Decimal(1)
        two = Decimal(2)
        three = Decimal(3)
        exact_balance = one / two.sqrt()
        a_value = one
        rho = Decimal(mcpr["regular_C3_shape"]["rho"])
        abs_b = rho / two
        modulus_ratio = abs_b / a_value
        q_from_ratio = koide_from_modulus_ratio(modulus_ratio)
        q_mcpr = Decimal(mcpr["regular_C3_shape"]["sqrt_mass_invariant"])
        measured_q, measured_masses = measured_central_q(references)
        measured_ratio = modulus_ratio_from_koide(measured_q)
        two_thirds = two / three
        measured_mu_over_e = (
            Decimal(measured_masses["muon"]) / Decimal(measured_masses["electron"])
        )
        mcpr_mu_over_e = Decimal(
            mcpr["conditional_prediction"]["ratios"]["m_mu_over_m_e"]
        )
        modulus_ratio_minus_exact_balance = modulus_ratio - exact_balance
        relative_modulus_defect_ppm = (
            modulus_ratio / exact_balance - one
        ) * Decimal(1_000_000)
        q_mcpr_minus_two_thirds = q_mcpr - two_thirds
        measured_q_minus_two_thirds = measured_q - two_thirds
        q_mcpr_minus_measured = q_mcpr - measured_q
        mcpr_mu_over_e_minus_measured = mcpr_mu_over_e - measured_mu_over_e

    phase = float(mcpr["regular_C3_shape"]["delta"])
    current_roots = roots(float(a_value), float(abs_b), phase)
    positive_phase_samples = (0.0, 0.1, phase)
    outside_phase = 0.4
    phase_samples = [
        {
            "phase": format(sample, ".17g"),
            "roots": [format(value, ".17g") for value in roots(1.0, 2.0**-0.5, sample)],
            "all_roots_nonnegative": min(roots(1.0, 2.0**-0.5, sample))
            >= -1.0e-14,
            "physical_Q": format(
                physical_koide(roots(1.0, 2.0**-0.5, sample)),
                ".17g",
            ),
        }
        for sample in positive_phase_samples
    ]
    outside_roots = roots(1.0, 2.0**-0.5, outside_phase)

    checks = {
        "face_carrier_has_C3_stabilizer": (
            face.get("face_orbit_theorem", {}).get("face_stabilizer_order") == 3
        ),
        "face_carrier_physical_attachment_open": (
            face.get("public_charged_mass_promotion_allowed") is False
        ),
        "finite_GNS_theorem_is_conditional": (
            finite_gns.get("public_koide_promotion_allowed") is False
        ),
        "MCPR_is_historically_target_informed": (
            mcpr.get("provenance_boundary", {})
            .get("historical_selection", {})
            .get("architecture_target_informed")
            is True
        ),
        "circulant_identity_matches_MCPR": abs(q_from_ratio - q_mcpr)
        < Decimal("1e-60"),
        "MCPR_roots_are_in_positive_chamber": min(current_roots) > 0.0,
        "balanced_positive_phase_samples_have_Q_two_thirds": all(
            row["all_roots_nonnegative"]
            and abs(float(row["physical_Q"]) - 2.0 / 3.0) < 1.0e-14
            for row in phase_samples
        ),
        "outside_phase_signed_formula_is_not_physical_Q": (
            min(outside_roots) < 0.0
            and abs(physical_koide(outside_roots) - 2.0 / 3.0) > 1.0e-3
        ),
        "PDG_2026_mass_coordinate_is_exact": (
            measured_masses == EXPECTED_PDG_2026_MASSES_MEV
            and all(
                references["entries"][name]["source"]["edition"] == "2026"
                for name in EXPECTED_PDG_2026_MASSES_MEV
            )
        ),
        "current_PDG_central_Q_recomputed": (
            measured_q == EXPECTED_PDG_2026_Q
        ),
    }

    return {
        "artifact": "oph_koide_circulant_identity",
        "schema": "oph.koide_circulant_identity.v1",
        "status": (
            "CLOSED_CIRCULANT_IDENTITY__FINITE_GNS_BALANCE_CONDITIONAL__"
            "PHYSICAL_ATTACHMENT_OPEN"
        ),
        "claim_class": "conditional_implication",
        "source_only_physical_prediction": False,
        "public_physical_promotion_allowed": False,
        "identity": {
            "carrier": "C = a I + b R + conjugate(b) R^2, R^3 = I",
            "eigenvalues": "lambda_k = a + 2 |b| cos(arg(b) + 2 pi k/3)",
            "domain": (
                "a>0 and all three eigenvalues nonnegative, so physical "
                "sqrt(m_k)=lambda_k"
            ),
            "cosine_sum": "0",
            "cosine_square_sum": "3/2",
            "Q_formula": "Q = 1/3 + (2/3) (|b|/a)^2",
            "equivalence": "Q = 2/3 iff |b|/a = 1/sqrt(2)",
            "phase_content": (
                "Q is phase-independent inside the positive chamber; the "
                "phase remains free and carries the two independent ratios."
            ),
            "outside_chamber_boundary": (
                "When an eigenvalue is negative, physical square roots use "
                "absolute eigenvalues and the signed formula is not physical Q."
            ),
        },
        "conditional_finite_GNS_result": {
            "artifact": "oph_charged_koide_orientation_isometry",
            "status": finite_gns.get("status"),
            "abs_b_over_a": decimal_text(exact_balance),
            "Q": decimal_text(two_thirds),
            "physical_attachment_open": True,
        },
        "current_MCPR_coordinate": {
            "a": decimal_text(a_value),
            "abs_b": decimal_text(abs_b),
            "arg_b": mcpr["regular_C3_shape"]["delta"],
            "abs_b_over_a": decimal_text(modulus_ratio),
            "Q": decimal_text(q_mcpr),
            "abs_b_over_a_minus_exact_balance": decimal_text(
                modulus_ratio_minus_exact_balance
            ),
            "relative_modulus_defect_ppm": decimal_text(
                relative_modulus_defect_ppm
            ),
            "Q_minus_two_thirds": decimal_text(q_mcpr_minus_two_thirds),
            "m_mu_over_m_e": decimal_text(mcpr_mu_over_e),
            "provenance": (
                "retrospective historically target-informed declared-model "
                "coordinate; zero charged targets are read at runtime"
            ),
        },
        "compare_only_PDG_2026_central_coordinate": {
            "masses_MeV": measured_masses,
            "Q": decimal_text(measured_q),
            "Q_minus_two_thirds": decimal_text(measured_q_minus_two_thirds),
            "inferred_abs_b_over_a": decimal_text(measured_ratio),
            "m_mu_over_m_e": decimal_text(measured_mu_over_e),
            "MCPR_Q_minus_measured_central": decimal_text(
                q_mcpr_minus_measured
            ),
            "MCPR_m_mu_over_m_e_minus_measured_central": decimal_text(
                mcpr_mu_over_e_minus_measured
            ),
            "significance_claimed": False,
            "note": (
                "The current PDG central masses give Q near the MCPR "
                "coordinate, not 0.666660511. The MCPR architecture was "
                "historically target-informed, so this near-equality is a "
                "diagnostic and not prospective evidence."
            ),
        },
        "phase_controls": {
            "positive_samples": phase_samples,
            "outside_sample": {
                "phase": format(outside_phase, ".17g"),
                "roots": [format(value, ".17g") for value in outside_roots],
                "physical_Q": format(physical_koide(outside_roots), ".17g"),
            },
        },
        "provenance": {
            "mcpr": {
                "path": MCPR.relative_to(ROOT).as_posix(),
                "sha256": sha256_file(MCPR),
            },
            "finite_GNS": {
                "path": FINITE_GNS.relative_to(ROOT).as_posix(),
                "sha256": sha256_file(FINITE_GNS),
            },
            "face_carrier": {
                "path": FACE_CARRIER.relative_to(ROOT).as_posix(),
                "sha256": sha256_file(FACE_CARRIER),
            },
            "PDG_reference_fixture": {
                "path": PDG_REFERENCES.relative_to(ROOT).as_posix(),
                "sha256": sha256_file(PDG_REFERENCES),
                "role": "compare_only",
            },
        },
        "checks": checks,
        "checks_pass": all(checks.values()),
    }


def canonical_bytes(payload: dict[str, Any]) -> bytes:
    return (json.dumps(payload, indent=2, ensure_ascii=False) + "\n").encode(
        "utf-8"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    artifact = build_artifact()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(canonical_bytes(artifact))
    print(
        json.dumps(
            {
                "status": artifact["status"],
                "checks_pass": artifact["checks_pass"],
                "output": str(args.output),
            },
            indent=2,
        )
    )
    return 0 if artifact["checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
