#!/usr/bin/env python3
"""Fail-closed audit of the clock-to-source-energy conversion lane.

The exact theorem E_star = h*nu_clk/epsilon_clk converts a source-emitted
positive dimensionless clock gap into an absolute energy scale.  This audit
reads the live R_gamma_noG_DAG certificate, recomputes the candidate E_star
and the G identity implied by the stored epsilon_Cs at high precision, and
scans hierarchy/certificates for the five required component producers.

The comparison against the displayed G coordinate is a checksum-identity
taint diagnosis, never a prediction.  The lane verdict is fail-closed:
promotion requires all five component producers and a source-emitted clock
gap, and four component producers are absent from the live repository.

This audit generates no source clock.  It consumes no measured particle
mass.
"""

from __future__ import annotations

import argparse
import json
from decimal import Decimal, localcontext
from pathlib import Path
from typing import Any

import mpmath as mp

WORKING_PREC = 80
WORKING_DPS = 80

ROOT = Path(__file__).resolve().parents[2]
CERTIFICATES_DIR = ROOT / "particles" / "hierarchy" / "certificates"
CLOCK_CERTIFICATE = CERTIFICATES_DIR / "R_gamma_noG_DAG_certificate.json"
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "clock"
    / "clock_source_energy_closure_audit.json"
)

# Exact SI definitions.
H_SI = Decimal("6.62607015e-34")
NU_CS_HZ = Decimal("9192631770")
GEV_J = Decimal("1.602176634e-10")
C_SI = Decimal("299792458")

# CODATA-style display value used only for the checksum-identity diagnosis.
G_DISPLAY_SI = "6.6743e-11"
G_DISPLAY_MATCH_TOL = mp.mpf("1e-6")

REQUIRED_COMPONENTS = (
    "R_U",
    "R_alpha",
    "R_e_abs",
    "R_QCD_nuc_133Cs",
    "R_atom_133Cs",
)
R_U_REQUIRED_FILES = (
    "R_U_interval_certificate.json",
    "R_U_krawczyk_certificate.json",
)
REQUIRED_NEXT_SOURCE_OBJECTS = (
    "R_ALPHA_ATOMIC_SCHEME_PACKET",
    "R_ELECTRON_ABSOLUTE_RATIO_PACKET",
    "R_QCD_NUCLEAR_133CS_PACKET",
    "R_ATOMIC_133CS_HAMILTONIAN_PACKET",
    "R_ATOMIC_133CS_SPECTRAL_ENCLOSURE",
    "R_CLOCK_REFINEMENT_LIMIT",
    "R_NO_G_CLOCK_DAG",
    "R_SOURCE_ENERGY_INTERVAL",
    "PHYSICAL_QF1_TO_QF9_FLAVOR_CARRIER_CERTIFICATE",
    "P_ALPHA_U_STRICT_NO_TARGET_ROOT",
    "COMMON_SCALE_PHYSICAL_YUKAWA_PACKET",
    "FROZEN_RG_THRESHOLD_MATCHING_UNCERTAINTY_PACKET",
)


def candidate_scale_from_gap(epsilon: Decimal) -> tuple[Decimal, Decimal]:
    """Exact SI conversion of a dimensionless clock gap into E_star."""

    if epsilon <= 0:
        raise ValueError("epsilon_Cs must be a positive decimal")
    e_j = H_SI * NU_CS_HZ / epsilon
    e_gev = H_SI * NU_CS_HZ / (epsilon * GEV_J)
    return e_j, e_gev


def inferred_g_from_gap(epsilon: str) -> mp.mpf:
    """Evaluate G_SI = c^5/(4*pi^2*hbar*nu_Cs^2) * epsilon_Cs^2, hbar = h/(2*pi)."""

    with mp.workdps(WORKING_DPS):
        eps = mp.mpf(epsilon)
        h = mp.mpf(str(H_SI))
        hbar = h / (2 * mp.pi)
        c = mp.mpf(str(C_SI))
        nu = mp.mpf(str(NU_CS_HZ))
        return +(c**5 * eps**2 / (4 * mp.pi**2 * hbar * nu**2))


def component_ledger(certificates_dir: Path) -> tuple[dict[str, str], dict[str, list[str]]]:
    """Scan the certificate directory for producer files per required component.

    R_U is supplied when both named R_U certificates exist.  Every other
    component is absent unless a certificate file with a matching name
    exists.
    """

    names = sorted(path.name for path in certificates_dir.glob("*.json"))
    producer_files: dict[str, list[str]] = {}
    ledger: dict[str, str] = {}
    r_u_files = [name for name in R_U_REQUIRED_FILES if (certificates_dir / name).exists()]
    producer_files["R_U"] = r_u_files
    ledger["R_U"] = "supplied" if len(r_u_files) == len(R_U_REQUIRED_FILES) else "absent"
    for component in REQUIRED_COMPONENTS[1:]:
        hits = [name for name in names if component.lower() in name.lower()]
        producer_files[component] = hits
        ledger[component] = "supplied" if hits else "absent"
    return ledger, producer_files


def conversion_theorem_block() -> dict[str, str]:
    """Exact conversion theorem with its interval and quark corollaries."""

    return {
        "theorem": "CLOCK_TO_SOURCE_ENERGY_THEOREM",
        "alias": "CLOCK_YO_SOURCE_ENERGY_THEOREM",
        "hypothesis": (
            "The source emits a positive dimensionless clock gap "
            "epsilon_clk = DeltaE_clk/E_star and the operational clock "
            "convention fixes nu_clk with DeltaE_clk = h*nu_clk."
        ),
        "conclusion": "E_star = h*nu_clk/epsilon_clk",
        "interval_corollary": (
            "For a certified positive interval "
            "0 < eps_lo <= epsilon_clk <= eps_hi, E_star lies in "
            "[h*nu_clk/eps_hi, h*nu_clk/eps_lo]."
        ),
        "quark_corollary": (
            "For a source-emitted ratio r_q = m_q/E_star, "
            "m_q_GeV = r_q*h*nu_clk/(epsilon_clk*J_GeV)."
        ),
        "mathematical_status": (
            "PROVED_CONDITIONALLY_ON_A_SOURCE_EMITTED_POSITIVE_CLOCK_GAP"
        ),
        "application_status": "OPEN_BECAUSE_EPSILON_CS_IS_NOT_SOURCE_EMITTED",
    }


def build_artifact(
    certificate: dict[str, Any] | None = None,
    certificates_dir: Path | None = None,
    certificate_path: str | None = None,
) -> dict[str, Any]:
    """Audit the clock lane against a certificate and the live producer files."""

    with localcontext() as ctx:
        ctx.prec = WORKING_PREC
        return _build_artifact_scoped(certificate, certificates_dir, certificate_path)


def _build_artifact_scoped(
    certificate: dict[str, Any] | None,
    certificates_dir: Path | None,
    certificate_path: str | None,
) -> dict[str, Any]:
    if certificate is None:
        certificate = json.loads(CLOCK_CERTIFICATE.read_text(encoding="utf-8"))
        source_path = str(CLOCK_CERTIFICATE)
    else:
        source_path = certificate_path if certificate_path is not None else "in_memory_certificate"
    scan_dir = CERTIFICATES_DIR if certificates_dir is None else certificates_dir

    epsilon = Decimal(certificate["candidate_values_from_sources"]["epsilon_Cs"])
    e_j, e_gev = candidate_scale_from_gap(epsilon)
    epsilon_round_trip = H_SI * NU_CS_HZ / (e_gev * GEV_J)
    round_trip_residual = abs(epsilon_round_trip - epsilon) / epsilon

    g_inferred = inferred_g_from_gap(str(epsilon))
    g_display = mp.mpf(G_DISPLAY_SI)
    relative_g_difference = abs(g_inferred - g_display) / g_display
    reproduces_displayed_g = bool(relative_g_difference < G_DISPLAY_MATCH_TOL)

    ledger, producer_files = component_ledger(scan_dir)
    components_absent = sorted(
        name for name, state in ledger.items() if state != "supplied"
    )
    final_packet_ready = (
        not components_absent
        and certificate.get("status") == "closed_source_clock_packet"
    )
    promotion_allowed = bool(final_packet_ready)
    verdict = (
        "FINAL_SOURCE_ONLY_PACKET_READY"
        if final_packet_ready
        else "NOT_FINAL_SOURCE_ONLY_PACKET"
    )

    checks = {
        "epsilon_is_positive": bool(epsilon > 0),
        "conversion_round_trip_residual_below_1e_60": bool(
            round_trip_residual < Decimal("1e-60")
        ),
        "component_ledger_covers_required_components": set(ledger) == set(REQUIRED_COMPONENTS),
        "promotion_blocked_when_any_component_absent": (
            not components_absent or promotion_allowed is False
        ),
    }
    checks_pass = all(bool(value) for value in checks.values())

    return {
        "artifact": "oph_clock_source_energy_closure_audit",
        "verdict": verdict,
        "final_packet_ready": final_packet_ready,
        "promotion_allowed": promotion_allowed,
        "application_status": "OPEN_BECAUSE_EPSILON_CS_IS_NOT_SOURCE_EMITTED",
        "exact_conversion_theorem": conversion_theorem_block(),
        "si_constants": {
            "h_J_s": str(H_SI),
            "nu_Cs_Hz": str(NU_CS_HZ),
            "GeV_J": str(GEV_J),
            "c_m_s": str(C_SI),
            "decimal_precision_digits": WORKING_PREC,
            "mpmath_dps": mp.mp.dps,
        },
        "clock_certificate_receipt": {
            "path": source_path,
            "status": certificate.get("status"),
            "components_required": certificate.get("components_required"),
            "component_status": certificate.get("component_status"),
            "formula_if_closed": certificate.get("formula_if_closed"),
            "candidate_epsilon_Cs": str(epsilon),
            "candidate_E_star_J": str(e_j),
            "candidate_E_star_GeV": str(e_gev),
            "conversion_round_trip_residual": str(round_trip_residual),
            "epsilon_role": "calibration_checksum_reproducing_displayed_G",
            "G_algebraically_inferred_from_candidate_epsilon": mp.nstr(g_inferred, 50),
            "displayed_G_SI_for_checksum_diagnosis": G_DISPLAY_SI,
            "relative_difference_from_displayed_G": mp.nstr(relative_g_difference, 20),
            "reproduces_displayed_G_within_1e_6": reproduces_displayed_g,
            "displayed_G_comparison_role": "taint_diagnosis_never_a_prediction",
            "interpretation": (
                "The stored epsilon_Cs algebraically reproduces the "
                "displayed G coordinate through the certificate identity "
                "G_SI = c^5/(4*pi^2*hbar*nu_Cs^2) * epsilon_Cs^2.  "
                "Independent R_alpha, R_e_abs, R_QCD_nuc_133Cs, and "
                "R_atom_133Cs producers are absent, so the stored decimal "
                "functions as a calibration checksum."
            ),
        },
        "component_ledger": ledger,
        "component_producer_files": producer_files,
        "components_absent": components_absent,
        "required_next_source_objects": list(REQUIRED_NEXT_SOURCE_OBJECTS),
        "checks": checks,
        "checks_pass": checks_pass,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, default=CLOCK_CERTIFICATE)
    parser.add_argument("--certificates-dir", type=Path, default=CERTIFICATES_DIR)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    artifact = build_artifact(
        json.loads(args.certificate.read_text(encoding="utf-8")),
        certificates_dir=args.certificates_dir,
        certificate_path=str(args.certificate),
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "verdict": artifact["verdict"],
                "promotion_allowed": artifact["promotion_allowed"],
                "checks_pass": artifact["checks_pass"],
            },
            indent=2,
        )
    )
    return 0 if artifact["checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
