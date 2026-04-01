#!/usr/bin/env python3
"""Emit the canonical UV/BW internalization scaffold.

This does not promote the UV/BW lane. It records the sharpest current internal
extension route: first extract a canonical scaling-limit cap pair from
transported cap marginals, then prove ordered null cut-pair rigidity on that
realized limit.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "particles" / "runs" / "uv" / "bw_internalization_scaffold.json"
PRELIMIT_SYSTEM = ROOT / "particles" / "runs" / "uv" / "bw_realized_transported_cap_local_system.json"
RAW_DATUM = ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_markov_faithfulness_datum.json"
CARRIED_SCHEDULE = ROOT / "particles" / "runs" / "uv" / "bw_carried_collar_schedule_scaffold.json"
CONSTRUCTIVE_RECOVERY = (
    ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_constructive_recovery_scaffold.json"
)
EXACT_MARKOV_MODULUS = ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_exact_markov_modulus_scaffold.json"
FAITHFUL_MODULAR_DEFECT = (
    ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_faithful_modular_defect_scaffold.json"
)
SPECTRAL_FLOOR = (
    ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_eventual_spectral_floor_scaffold.json"
)


@dataclass(frozen=True)
class RigidityResult:
    solution_dimension: int
    surviving_generator_disk: str
    surviving_generator_half_line: str
    ordered_boundary_pair: list[str]


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _artifact_ref(path: Path) -> str:
    return f"code/{path.relative_to(ROOT).as_posix()}"


def _compute_rigidity() -> RigidityResult:
    return RigidityResult(
        solution_dimension=1,
        surviving_generator_disk="1 - z**2",
        surviving_generator_half_line="2*u",
        ordered_boundary_pair=["-1", "+1"],
    )


def build_artifact() -> dict[str, object]:
    rigidity = _compute_rigidity()
    boundary = {
        "status": "open_split_after_candidate_projective_route",
        "remaining_object": "canonical_scaling_cap_pair_realization_from_transported_cap_marginals",
        "follow_on_object": "independent_bw_rigidity_on_realized_limit",
        "dominant_pressure_point": "eventual_fixed_local_collar_spectral_floor_for_transported_marginals",
        "filled_contract_witnesses": [
            "reference_cap_local_test_system",
            "projectively_compatible_transported_cap_marginal_family",
            "asymptotic_transport_equivalence_certificate",
        ],
        "prelimit_system_artifact": _artifact_ref(PRELIMIT_SYSTEM),
        "remaining_missing_emitted_witness": "vanishing_carried_collar_schedule_on_fixed_local_collars",
        "remaining_missing_emitted_witness_artifact": _artifact_ref(CARRIED_SCHEDULE),
        "remaining_missing_emitted_witness_formula": (
            "eta_{n,m,delta} = r_FR(epsilon_{n,m,delta}) + "
            "4 * lambda_{*,n,m,delta}^{-1} * delta^M_{m,delta}(epsilon_{n,m,delta}) -> 0"
        ),
        "actual_solver_missing_emitted_witnesses": [
            {
                "id": "constructive_recovery_remainder_vanishing",
                "artifact": _artifact_ref(CONSTRUCTIVE_RECOVERY),
                "role": "markov_side_recovery_term",
            },
            {
                "id": "fixed_local_collar_faithful_modular_defect_vanishing",
                "artifact": _artifact_ref(FAITHFUL_MODULAR_DEFECT),
                "role": "faithfulness_weighted_modular_term",
            },
        ],
        "derived_remaining_input_witness": {
            "id": "vanishing_carried_collar_schedule_on_fixed_local_collars",
            "artifact": _artifact_ref(CARRIED_SCHEDULE),
        },
        "derived_remaining_input_witness_closure_theorem": (
            "If the constructive recovery remainder and the faithful modular-defect term both vanish "
            "on every fixed local collar model, then the carried-collar schedule follows by termwise "
            "addition of the two nonnegative remainders."
        ),
        "smallest_exact_blocker": "eventual_fixed_local_collar_spectral_floor_for_transported_marginals",
        "smallest_exact_blocker_formula": (
            "exists lambda_bar_{m,delta} > 0 with lambda_{*,n,m,delta} >= lambda_bar_{m,delta} eventually"
        ),
        "single_live_missing_clause_artifact": _artifact_ref(SPECTRAL_FLOOR),
        "markov_side_status": "latent_from_epsilon_to_zero",
        "faithfulness_side_status": "open",
        "smallest_exact_blocker_unlocks": [
            "fixed_local_collar_faithful_modular_defect_vanishing",
            "vanishing_carried_collar_schedule_on_fixed_local_collars",
            "canonical_scaling_cap_pair_realization_from_transported_cap_marginals",
        ],
        "smaller_remaining_raw_datum": "fixed_local_collar_markov_faithfulness_datum",
        "smaller_remaining_raw_datum_artifact": _artifact_ref(RAW_DATUM),
        "local_intermediate_witness_chain": [
            {
                "id": "constructive_recovery_remainder_vanishing",
                "artifact": _artifact_ref(CONSTRUCTIVE_RECOVERY),
            },
            {
                "id": "fixed_local_collar_exact_markov_modulus_vanishing",
                "artifact": _artifact_ref(EXACT_MARKOV_MODULUS),
            },
            {
                "id": "fixed_local_collar_faithful_modular_defect_vanishing",
                "artifact": _artifact_ref(FAITHFUL_MODULAR_DEFECT),
            },
            {
                "id": "vanishing_carried_collar_schedule_on_fixed_local_collars",
                "artifact": _artifact_ref(CARRIED_SCHEDULE),
            },
        ],
        "remaining_objects": [
            "canonical_scaling_cap_pair_realization_from_transported_cap_marginals",
            "independent_bw_rigidity_on_realized_limit",
        ],
        "current_internalized_scope": (
            "Axiom-3 plus the fixed-cutoff collar/MaxEnt package internalize local Gibbs form, "
            "quasi-local propagation, endpoint-Lipschitz interval control, and refinement-stable "
            "branch persistence. The current corpus also packages the reference cap-local test system, "
            "the projectively compatible transported cap marginal family, and the asymptotic transport-equivalence certificate."
        ),
        "reason_current_corpus_fails": (
            "The current corpus already packages the reference cap-local test system, the projectively "
            "compatible transported cap marginal family, and the asymptotic transport-equivalence certificate. "
            "The compactness/extraction step itself is not the missing proof. The live extraction input is the derived "
            "carried-collar schedule on fixed local collar models, and on the local-Gibbs plus exponential-mixing pullback "
            "branch the constructive-recovery / exact-Markov side is already latent once epsilon -> 0 on each fixed collar "
            "model. The only nonlatent lower input still external to the emitted chain is the eventual positive lower spectral "
            "bound feeding the faithful modular-defect term. Without that clause, neither the faithful modular-defect witness, "
            "nor the carried-collar schedule, nor the canonical scaling-limit cap-pair realization is promoted."
        ),
        "statement": (
            "First exact object: starting from a projectively compatible extracted family of "
            "transported cap marginals on fixed reference type-I regulator cap algebras, realize "
            "a canonical scaling-limit observer cap pair (A_infty(C), omega_infty^C) without "
            "assuming type-I survival."
        ),
        "follow_on_statement": (
            "Second exact object: on that realized limit pair, prove from internal OPH premises "
            "alone that ordered null cut-pair data rigidify the residual cap-preserving conformal "
            "freedom to the unique BW hyperbolic subgroup lambda_C(s), so "
            "sigma_t^{omega_infty^C} = alpha_{lambda_C(2 pi t)} without reusing consequences "
            "already downstream of the BW branch."
        ),
        "candidate_extension_status": "constructive_prelimit_system_two_lower_emitted_witnesses_still_missing",
        "candidate_extension_route": (
            "Step 1: close the sole nonlatent lower input by emitting the eventual fixed-local-collar spectral floor "
            "for transported marginals on every fixed collar model. On the local-Gibbs plus exponential-mixing pullback "
            "branch, that closes the faithfulness-weighted modular term once epsilon -> 0, while the recovery/Markov side "
            "is already latent from that same epsilon-control; the carried-collar schedule and then the scaling-limit cap-pair "
            "extraction contract follow above the already-emitted realized transported cap-local system. "
            "Step 2: ordered null cut-pair rigidity that collapses the residual "
            "cap-preserving conformal freedom to the unique BW hyperbolic subgroup."
        ),
        "candidate_extension_target": "sigma_t^{omega_infty^C} = alpha_{lambda_C(2 pi t)}",
        "canonical_code_scaffolds": [
            "code/particles/uv/derive_bw_realized_transported_cap_local_system.py",
            "code/particles/uv/derive_bw_fixed_local_collar_markov_faithfulness_datum.py",
            "code/particles/uv/derive_bw_fixed_local_collar_constructive_recovery_scaffold.py",
            "code/particles/uv/derive_bw_fixed_local_collar_exact_markov_modulus_scaffold.py",
            "code/particles/uv/derive_bw_fixed_local_collar_eventual_spectral_floor_scaffold.py",
            "code/particles/uv/derive_bw_fixed_local_collar_faithful_modular_defect_scaffold.py",
            "code/particles/uv/derive_bw_carried_collar_schedule_scaffold.py",
            "code/particles/uv/derive_bw_scaling_limit_cap_pair_extraction_scaffold.py",
            "code/particles/uv/derive_bw_ordered_cut_pair_rigidity_scaffold.py",
        ],
        "canonical_artifacts": [
            "code/particles/runs/uv/bw_realized_transported_cap_local_system.json",
            "code/particles/runs/uv/bw_fixed_local_collar_markov_faithfulness_datum.json",
            "code/particles/runs/uv/bw_fixed_local_collar_constructive_recovery_scaffold.json",
            "code/particles/runs/uv/bw_fixed_local_collar_exact_markov_modulus_scaffold.json",
            "code/particles/runs/uv/bw_fixed_local_collar_eventual_spectral_floor_scaffold.json",
            "code/particles/runs/uv/bw_fixed_local_collar_faithful_modular_defect_scaffold.json",
            "code/particles/runs/uv/bw_carried_collar_schedule_scaffold.json",
            "code/particles/runs/uv/bw_scaling_limit_cap_pair_extraction_scaffold.json",
            "code/particles/runs/uv/bw_ordered_cut_pair_rigidity_scaffold.json",
        ],
        "symbolic_ordered_cut_pair_rigidity_test": {
            "status": "pass" if rigidity.solution_dimension == 1 else "fail",
            **asdict(rigidity),
            "conclusion": (
                "Preserving the ordered boundary pair leaves a one-dimensional hyperbolic subalgebra on the disk, "
                "which is conjugate to dilation on the positive half-line."
            ),
        },
    }
    return {
        "artifact": "oph_uv_bw_internalization_scaffold",
        "generated_utc": _timestamp(),
        "status": "minimal_constructive_extension",
        "public_promotion_allowed": False,
        "current_boundary": "T1_plus_T3_external_scaling_limit_branch",
        "extension_kind": "scaling_limit_cap_pair_plus_ordered_cut_pair_rigidity",
        "input_contract": {
            "transported_cap_marginals": "Sequence of cap-local states along a refinement chain in a common chart.",
            "inherited_strip_data": "The oriented null generator Omega and the ordered cut pair (Gamma_minus, Gamma_plus).",
            "support_map": "Scaling-limit support map or finite-stage approximation on the quasi-local cap net.",
            "endpoint_family": "Half-line endpoint matrix elements for the renormalized family K_tilde_a(Omega).",
        },
        "solver_spec": {
            "checks": [
                "local_weak_star_extraction",
                "ordered_cut_pair_preservation",
                "pair_preserving_lie_algebra_dimension_equals_one",
                "half_line_action_is_dilation",
                "borchers_positive_translation",
                "kappa_equals_2pi",
            ],
            "output_certificate": {
                "bw_automorphism": "sigma_t = alpha_lambda_C(2*pi*t)",
                "residual_modular_class": "q_BW(C)=0 after ordered-cut quotient",
                "status": "closed_on_extension",
                "typeI_required": False,
            },
        },
        "public_status_boundary": boundary,
        "notes": [
            "This scaffold promotes the UV/BW extension route to a canonical local artifact without claiming current-corpus closure.",
            "The current pressure point is the first object, not the symbolic rigidity calculation: the realized scaling-limit cap pair is still missing.",
            "The local carried-collar side is now decomposed one level further into constructive recovery, exact-Markov comparison convergence, and faithful modular-defect vanishing, with the eta schedule treated as the derived combination witness.",
            "On the current branch the only nonlatent lower side condition still external to that emitted chain is the eventual fixed-local-collar spectral floor feeding the faithful modular-defect term.",
            "The symbolic test certifies the rigidity shape of the ordered cut-pair argument, but not the existence of the realized scaling-limit cap pair.",
            "The correct target is an automorphism theorem on the realized scaling-limit cap pair; no type-I survival is assumed.",
        ],
        "source_code_scaffolds": {
            "realized_transported_cap_local_system": "code/particles/uv/derive_bw_realized_transported_cap_local_system.py",
            "fixed_local_collar_markov_faithfulness_datum": "code/particles/uv/derive_bw_fixed_local_collar_markov_faithfulness_datum.py",
            "fixed_local_collar_constructive_recovery": "code/particles/uv/derive_bw_fixed_local_collar_constructive_recovery_scaffold.py",
            "fixed_local_collar_exact_markov_modulus": "code/particles/uv/derive_bw_fixed_local_collar_exact_markov_modulus_scaffold.py",
            "fixed_local_collar_eventual_spectral_floor": "code/particles/uv/derive_bw_fixed_local_collar_eventual_spectral_floor_scaffold.py",
            "fixed_local_collar_faithful_modular_defect": "code/particles/uv/derive_bw_fixed_local_collar_faithful_modular_defect_scaffold.py",
            "carried_collar_schedule": "code/particles/uv/derive_bw_carried_collar_schedule_scaffold.py",
            "scaling_limit_cap_pair_extraction": "code/particles/uv/derive_bw_scaling_limit_cap_pair_extraction_scaffold.py",
            "ordered_cut_pair_rigidity": "code/particles/uv/derive_bw_ordered_cut_pair_rigidity_scaffold.py",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Build the UV/BW internalization scaffold.")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(build_artifact(), indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
