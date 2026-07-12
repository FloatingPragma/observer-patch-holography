#!/usr/bin/env python3
"""Derive the dynamics-level obstruction to a nonuniform twelve-port record."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from derive_charged_family_non_singlet_port_attachment import port_moments


ROOT = Path(__file__).resolve().parents[2]
SCREEN = (
    ROOT / "particles" / "hierarchy" / "certificates"
    / "R_screen_sieve_icosahedral_certificate.json"
)
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "leptons"
    / "charged_nonuniform_port_record_dynamics_no_go.json"
)


def convex_cost(weights: list[float]) -> float:
    return sum(value * value for value in weights)


def uniform_minimum(total: float, count: int) -> list[float]:
    return [total / count] * count


def build_artifact(screen: dict[str, Any]) -> dict[str, Any]:
    minimum = screen["convex_defect_minimum"]
    charges = [float(value) for value in minimum["charges"]]
    total = float(minimum["total_charge"])
    uniform = uniform_minimum(total, len(charges))
    first, quadrupole = port_moments(uniform)
    uniform_forced = charges == uniform and len(charges) == 12
    return {
        "artifact": "oph_charged_nonuniform_port_record_dynamics_no_go",
        "status": "CLOSED_NO_GO_CURRENT_OPH_SCREEN_DYNAMICS_FORCE_UNIFORM_RECORD",
        "nonuniform_record_derived": False,
        "public_charged_mass_promotion_allowed": False,
        "declared_branch": {
            "total_curvature_charge": total,
            "port_count": len(charges),
            "convex_cost": "C(q)=sum_i q_i^2",
            "no_marked_point": True,
            "homogeneous_MaxEnt_constraints": True,
            "unique_minimizer": uniform,
            "uniform_minimizer_cost": convex_cost(uniform),
        },
        "theorem": {
            "id": "homogeneous_MaxEnt_icosahedral_dynamics_forbid_deterministic_non_singlet_port_record",
            "statement": (
                "On the declared twelve-port screen branch, fixed total charge 12 and "
                "strictly convex defect cost force q_i=1 at every port. The homogeneous "
                "no-marked-point MaxEnt constraints are A5-invariant; uniqueness from strict "
                "entropy concavity therefore makes the expected port record A5-invariant and "
                "uniform. Its centered first and quadrupole moments vanish. Consequently the "
                "current OPH dynamics cannot derive a deterministic refinement-stable "
                "non-singlet port record."
            ),
            "proof": [
                "Cauchy-Schwarz gives sum_i q_i^2 >= (sum_i q_i)^2/12 = 12, with equality iff every q_i=1.",
                "The screen certificate attains equality with twelve unit defects.",
                "The constraint set and entropy are invariant under the transitive A5 port action.",
                "Strict concavity makes the MaxEnt state unique, so applying any A5 permutation returns the same state.",
                "Every one-port expectation is therefore equal; centering annihilates the record.",
                "Individual fluctuations do not select an ID-independent refinement-stable orientation unless an additional record/conditioning rule is emitted.",
            ],
            "checks": {
                "uniform_minimum_matches_certificate": uniform_forced,
                "first_moment_norm": float(np.linalg.norm(first)),
                "quadrupole_norm": float(np.linalg.norm(quadrupole)),
            },
        },
        "why_fluctuations_do_not_close_the_gap": {
            "ensemble_non_singlet_expectation": 0,
            "sample_orientation": "random_A5_orbit_without_canonical_family_attachment",
            "refinement_problem": "no emitted persistence or spontaneous-symmetry-breaking selection theorem",
        },
        "smallest_consistent_extension": {
            "id": "observer_conditioned_non_equilibrium_port_record_branch",
            "required_objects": [
                "a port-local gauge-invariant observable O_i in the existing screen algebra",
                "source-emitted nonuniform constraint values c_i=<O_i> or local multipliers lambda_i",
                "an observer readback record that selects the A5 orbit without port IDs",
                "a refinement transport theorem preserving the selected non-singlet moment and a simple quadrupole spectrum",
                "a physical argument attaching the resulting spectral lines to the charged sector",
                "a hash-bound no-target-leak packet frozen before charged-mass comparison",
            ],
            "branch_status": (
                "not contained in the homogeneous finite-parameter MaxEnt/refinement branch; "
                "the main paper classifies local multipliers as a non-equilibrium enlarged family"
            ),
            "not_allowed": "choosing c_i or lambda_i from electron, muon, or tau masses",
        },
        "claim_boundary": (
            "A nonuniform twelve-port record is not derivable from the current homogeneous "
            "OPH screen dynamics. It requires a new observer-conditioned non-equilibrium or "
            "spontaneous-symmetry-breaking branch with its own refinement theorem."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--screen", type=Path, default=SCREEN)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    artifact = build_artifact(json.loads(args.screen.read_text(encoding="utf-8")))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(artifact["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
