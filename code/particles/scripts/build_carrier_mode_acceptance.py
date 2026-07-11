#!/usr/bin/env python3
"""Emit the fail-closed classical-carrier / quantum-particle gate.

The gate is analytic bookkeeping, not a numerical simulation.  It records
what follows after a Maxwell, perturbative Yang--Mills, or pure-Einstein
quadratic action has been selected, while preventing an abstract symmetry
group or a classical field equation from being promoted to a quantum particle
claim.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[2]
PARTICLES_ROOT = ROOT / "particles"
DEFAULT_JSON_OUT = PARTICLES_ROOT / "runs" / "status" / "carrier_mode_acceptance.json"
DEFAULT_MD_OUT = PARTICLES_ROOT / "CARRIER_MODE_ACCEPTANCE.md"

CLASSICAL_REQUIRED = (
    "action_branch_declared",
    "positive_kinetic_coefficient_declared",
    "phase_and_field_content_declared",
    "quadratic_operator_exhibited",
    "gauge_fixing_exhibited",
    "physical_degrees_of_freedom_counted",
    "massless_characteristic_or_hamiltonian_exhibited",
)

QUANTUM_REQUIRED = (
    "quantization_constructed_from_oph",
    "positive_physical_hilbert_space",
    "physical_two_point_or_hamiltonian_spectrum_constructed",
    "positive_residue_massless_pole",
    "asymptotic_or_deconfined_particle_state",
)


def evaluate_acceptance(evidence: Mapping[str, bool]) -> dict[str, Any]:
    """Evaluate the two gates without inferring any missing receipt."""

    missing_classical = [name for name in CLASSICAL_REQUIRED if evidence.get(name) is not True]
    classical_passed = not missing_classical
    missing_quantum = [name for name in QUANTUM_REQUIRED if evidence.get(name) is not True]
    quantum_passed = classical_passed and not missing_quantum
    return {
        "classical_carrier_gate": {
            "passed": classical_passed,
            "status": (
                "conditional_pass_on_declared_action_phase_branch"
                if classical_passed
                else "not_passed"
            ),
            "missing": missing_classical,
        },
        "quantum_particle_gate": {
            "passed": quantum_passed,
            "status": "passed" if quantum_passed else "not_passed",
            "missing": missing_quantum,
        },
        "particle_promotion_allowed": quantum_passed,
    }


def _base_evidence() -> dict[str, bool]:
    return {
        "abstract_symmetry_group_reconstructed": True,
        "abstract_symmetry_group_alone_sufficient": False,
        "action_branch_declared": True,
        "positive_kinetic_coefficient_declared": True,
        "phase_and_field_content_declared": True,
        "quadratic_operator_exhibited": True,
        "gauge_fixing_exhibited": True,
        "physical_degrees_of_freedom_counted": True,
        "massless_characteristic_or_hamiltonian_exhibited": True,
        # None of the following is constructed by the OPH symmetry/Tannaka
        # step or by the classical field equation alone.
        "quantization_constructed_from_oph": False,
        "positive_physical_hilbert_space": False,
        "physical_two_point_or_hamiltonian_spectrum_constructed": False,
        "positive_residue_massless_pole": False,
        "asymptotic_or_deconfined_particle_state": False,
    }


def _carrier(
    *,
    carrier_id: str,
    label: str,
    branch: str,
    quadratic_action: str,
    gauge_fixing: str,
    quadratic_operator: str,
    physical_degrees_of_freedom: str,
    reduced_classical_hamiltonian: str,
    classical_spectrum: str,
    conditional_propagator: str,
    phase_boundary: str,
    nonclaims: list[str],
) -> dict[str, Any]:
    evidence = _base_evidence()
    gates = evaluate_acceptance(evidence)
    return {
        "carrier_id": carrier_id,
        "label": label,
        "claim_kind": "conditional_classical_carrier_mode",
        "branch": branch,
        "branch_is_additional_input_not_group_output": True,
        "abstract_symmetry_group_alone_sufficient": False,
        "hard_quadratic_mass_parameter_squared": 0,
        "hard_mass_parameter_semantics": (
            "Coefficient in the displayed free quadratic operator on the declared branch; "
            "not a measured rest mass and not by itself a quantum-particle pole claim."
        ),
        "quadratic_action": quadratic_action,
        "gauge_fixing": gauge_fixing,
        "quadratic_operator": quadratic_operator,
        "physical_degrees_of_freedom": physical_degrees_of_freedom,
        "reduced_classical_hamiltonian": reduced_classical_hamiltonian,
        "classical_spectrum": classical_spectrum,
        "conditional_free_two_point_function": conditional_propagator,
        "conditional_two_point_boundary": (
            "Displayed only as the result of a further standard free-field quantization assumption; "
            "the OPH construction has not supplied the required physical Hilbert-space and pole receipts."
        ),
        "phase_boundary": phase_boundary,
        "evidence": evidence,
        **gates,
        "explicit_nonclaims": nonclaims,
    }


def build_payload() -> dict[str, Any]:
    carriers = [
        _carrier(
            carrier_id="photon",
            label="Electromagnetic carrier",
            branch="unbroken_deconfined_Maxwell_branch_with_Z_A_positive",
            quadratic_action=(
                "S2[A]=-(Z_A/4) int F_{mu nu}F^{mu nu} "
                "-(Z_A/(2 xi)) int (partial_mu A^mu)^2"
            ),
            gauge_fixing="Lorenz covariant gauge partial_mu A^mu=0; quotient residual gauge directions",
            quadratic_operator="K_{mu nu}=Z_A[-k^2 eta_{mu nu}+(1-xi^{-1})k_mu k_nu]",
            physical_degrees_of_freedom="two transverse classical polarizations in four dimensions",
            reduced_classical_hamiltonian="H2=(Z_A/2) int (|E_T|^2+|B|^2), positive for Z_A>0",
            classical_spectrum="omega(k)^2=|k|^2; characteristic surface k^2=0",
            conditional_propagator="D^T_{mu nu}(k)=-i P^T_{mu nu}/[Z_A(k^2+i0)]",
            phase_boundary=(
                "Requires an unbroken U(1), a nonzero Maxwell kinetic term, and a deconfined phase. "
                "A Higgs or Stueckelberg realization and medium-induced effective masses are not excluded by symmetry alone."
            ),
            nonclaims=[
                "the compact group alone produces a connection or Maxwell kinetic term",
                "the OPH construction supplies a photon Fock space or LSZ pole",
                "absence of every effective photon mass in every phase or medium",
            ],
        ),
        _carrier(
            carrier_id="gluon",
            label="Color gauge carrier",
            branch="pure_Yang-Mills_quadratic_expansion_before_nonperturbative_confinement",
            quadratic_action=(
                "S2[A]=-(Z_g/4) int (partial_mu A_nu^a-partial_nu A_mu^a)^2 "
                "-(Z_g/(2 xi)) int (partial_mu A^{a mu})^2"
            ),
            gauge_fixing="covariant gauge with the usual constraint/BRST quotient at the perturbative free level",
            quadratic_operator="K^{ab}_{mu nu}=delta^{ab} Z_g[-k^2 eta_{mu nu}+(1-xi^{-1})k_mu k_nu]",
            physical_degrees_of_freedom="two perturbative transverse modes per Lie-algebra generator",
            reduced_classical_hamiltonian="H2=(Z_g/2) sum_a int (|E_T^a|^2+|B^a|^2), positive at the free level for Z_g>0",
            classical_spectrum="omega(k)^2=|k|^2 for the free quadratic modes; characteristic surface k^2=0",
            conditional_propagator="D^{T,ab}_{mu nu}(k)=-i delta^{ab} P^T_{mu nu}/[Z_g(k^2+i0)]",
            phase_boundary=(
                "This is a classical/perturbative pre-confinement statement. On the confining QCD phase, "
                "a colored gauge-potential pole is not an asymptotic gluon particle."
            ),
            nonclaims=[
                "a deconfined asymptotic gluon in four-dimensional confining QCD",
                "a continuum Yang-Mills mass-gap proof",
                "absence of gauge-invariant massive bound states",
            ],
        ),
        _carrier(
            carrier_id="graviton",
            label="Einstein tensor carrier",
            branch="pure_Einstein-Hilbert_linearization_about_a_suitable_Ricci-flat_background",
            quadratic_action=(
                "S2[h]=(Z_h/2) int [-1/2 partial_lambda h_mu_nu partial^lambda h^mu_nu "
                "+ partial_mu h^mu_nu partial^lambda h_lambda_nu - partial_mu h^mu_nu partial_nu h "
                "+ 1/2 partial_lambda h partial^lambda h], the Fierz-Pauli quadratic action from pure Einstein-Hilbert"
            ),
            gauge_fixing="de Donder gauge; on Minkowski reduce to the transverse-traceless sector",
            quadratic_operator="K_TT=Z_h k^2 P^(2)_TT",
            physical_degrees_of_freedom="two transverse-traceless classical tensor polarizations in four dimensions",
            reduced_classical_hamiltonian="H2_TT=(Z_h/2) int (|dot h_TT|^2+|grad h_TT|^2), positive for Z_h>0",
            classical_spectrum="omega(k)^2=|k|^2 on Minkowski; null characteristic k^2=0",
            conditional_propagator="D_TT(k)=i P^(2)_TT/[Z_h(k^2+i0)]",
            phase_boundary=(
                "Requires the pure Einstein kinetic branch and a background admitting the stated mode decomposition. "
                "Other diffeomorphism-invariant field content or higher-derivative terms may add massive modes."
            ),
            nonclaims=[
                "the classical Einstein equation constructs a graviton Hilbert space",
                "an OPH-derived quantum graviton pole or interacting S-matrix",
                "exclusion of every additional massive scalar or spin-two mode in extended theories",
            ],
        ),
    ]
    return {
        "artifact": "oph_massless_carrier_mode_acceptance",
        "schema": "oph_carrier_mode_quantum_particle_gate_v1",
        "github_issue": 536,
        "status": "classical_action_branch_modes_recorded_quantum_particle_gate_open",
        "abstract_symmetry_group_alone_passes_classical_gate": False,
        "abstract_symmetry_group_alone_passes_quantum_gate": False,
        "policy": (
            "A reconstructed compact group is not a particle. Public particle promotion requires every "
            "classical-action receipt and every independent quantum/phase/pole receipt."
        ),
        "classical_required_receipts": list(CLASSICAL_REQUIRED),
        "quantum_required_receipts": list(QUANTUM_REQUIRED),
        "carriers": carriers,
        "common_speed_semantics": {
            "structural_output": "shared invariant null cone / characteristic speed on the declared branch",
            "si_decimal": "299792458 m/s is the exact SI metre/second unit convention, not a predicted decimal magnitude",
        },
    }


def render_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Carrier-Mode and Quantum-Particle Acceptance",
        "",
        "This receipt keeps a classical massless quadratic mode separate from a quantum particle claim.",
        "An abstract compact group or a classical field equation alone passes neither quantum-particle promotion nor a public rest-mass prediction gate.",
        "",
        "| Carrier | Declared branch | Hard quadratic mass parameter | Classical mode gate | Quantum particle gate |",
        "| --- | --- | ---: | --- | --- |",
    ]
    for carrier in payload["carriers"]:
        lines.append(
            f"| {carrier['label']} | `{carrier['branch']}` | "
            f"`{carrier['hard_quadratic_mass_parameter_squared']}` | "
            f"`{carrier['classical_carrier_gate']['status']}` | "
            f"`{carrier['quantum_particle_gate']['status']}` |"
        )
    lines.extend(
        [
            "",
            "The zero in this table is an action parameter, not a `0 GeV` particle prediction. The conditional propagator formulas require a further free-field quantization assumption and do not satisfy the OPH quantum gate in the present corpus.",
            "",
            "For Yang--Mills, the displayed modes are perturbative/pre-confinement modes; confinement blocks an asymptotic colored-gluon interpretation. For gravity, the receipt is restricted to pure Einstein linearization on a suitable background and does not exclude extra modes in extended diffeomorphism-invariant theories.",
            "",
            f"Common-speed semantics: {payload['common_speed_semantics']['structural_output']}. "
            f"The value `{payload['common_speed_semantics']['si_decimal']}`.",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the fail-closed carrier/particle acceptance receipt.")
    parser.add_argument("--json-out", default=str(DEFAULT_JSON_OUT))
    parser.add_argument("--markdown-out", default=str(DEFAULT_MD_OUT))
    args = parser.parse_args()

    payload = build_payload()
    json_out = Path(args.json_out)
    json_out.parent.mkdir(parents=True, exist_ok=True)
    json_out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    markdown_out = Path(args.markdown_out)
    markdown_out.parent.mkdir(parents=True, exist_ok=True)
    markdown_out.write_text(render_markdown(payload) + "\n", encoding="utf-8")
    print(f"saved: {json_out}")
    print(f"saved: {markdown_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
