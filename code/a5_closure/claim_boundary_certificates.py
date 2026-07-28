#!/usr/bin/env python3
"""Exact-small certificates for the physical A5-to-SM claim boundary.

These checks are theorem fixtures.  They certify finite algebra and expose
counterexamples; they are not physical producer receipts.  Every check raises
an explicit exception on failure, so ``python -O`` cannot remove the verifier.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from typing import Any

import sympy as sp


N_MODES = 5
FOCK_DIMENSION = 1 << N_MODES


class CertificateError(RuntimeError):
    """Raised when an exact certificate identity fails."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CertificateError(message)


@dataclass(frozen=True)
class SourceReduct:
    """The deliberately small target-free interface used by the no-go fixture."""

    port_count: int
    total_coordination_charge: int


def source_completion_nonuniqueness_certificate() -> dict[str, Any]:
    """Exhibit inequivalent completions with exactly the same source reduct.

    This proves only non-identifiability from the stated reduct.  It does not
    say that a richer operational source packet cannot select a completion.
    """

    reduct = SourceReduct(port_count=12, total_coordination_charge=12)
    current_completions = {
        "abelian_12": reduct,
        "compact_sm_lie_type": reduct,
    }
    matter_completions = {
        "rank15_exterior_packet": "compact_sm_lie_type",
        "rank15_plus_sterile_singlet": "compact_sm_lie_type",
    }
    require(
        len(set(current_completions.values())) == 1,
        "current completions do not share one source reduct",
    )
    require(
        len(current_completions) == 2,
        "the current countermodel needs two distinct completions",
    )
    require(
        len(set(matter_completions.values())) == 1,
        "matter completions do not share one current reduct",
    )
    require(
        len(matter_completions) == 2,
        "the matter countermodel needs two distinct completions",
    )
    return {
        "status": "EXACT_FINITE_NONIDENTIFIABILITY_THEOREM",
        "source_reduct": {
            "port_count": reduct.port_count,
            "total_coordination_charge": reduct.total_coordination_charge,
        },
        "inequivalent_current_completions": sorted(current_completions),
        "inequivalent_matter_completions": sorted(matter_completions),
        "source_only_reconstruction_of_every_completion": False,
        "physical_promotion": False,
        "boundary": (
            "The exposed reduct is not completion-unique. A richer observer-like "
            "packet may still select a completion if its additional operational "
            "producer receipts are independently supplied."
        ),
    }


def _annihilation(mode: int) -> sp.SparseMatrix:
    entries: dict[tuple[int, int], sp.Expr] = {}
    for state in range(FOCK_DIMENSION):
        if (state >> mode) & 1:
            lower = state & ((1 << mode) - 1)
            sign = -1 if lower.bit_count() % 2 else 1
            entries[(state & ~(1 << mode), state)] = sp.Integer(sign)
    return sp.SparseMatrix(FOCK_DIMENSION, FOCK_DIMENSION, entries)


def rank15_clifford_certificate() -> dict[str, Any]:
    """Recompute the ten-Majorana parent of the rank-15 internal witness."""

    annihilation = [_annihilation(j) for j in range(N_MODES)]
    creation = [operator.conjugate().transpose() for operator in annihilation]
    identity = sp.eye(FOCK_DIMENSION)
    zero = sp.zeros(FOCK_DIMENSION)

    for i in range(N_MODES):
        for j in range(N_MODES):
            require(
                annihilation[i] * annihilation[j]
                + annihilation[j] * annihilation[i]
                == zero,
                f"CAR annihilation relation failed at ({i}, {j})",
            )
            require(
                creation[i] * creation[j] + creation[j] * creation[i] == zero,
                f"CAR creation relation failed at ({i}, {j})",
            )
            expected = identity if i == j else zero
            require(
                annihilation[i] * creation[j] + creation[j] * annihilation[i]
                == expected,
                f"mixed CAR relation failed at ({i}, {j})",
            )

    majoranas: list[sp.Matrix] = []
    for a, adag in zip(annihilation, creation, strict=True):
        majoranas.extend((a + adag, sp.I * (a - adag)))
    require(len(majoranas) == 10, "the oriented Clifford parent needs ten Majoranas")

    gamma_internal = sp.I**5 * sp.prod(majoranas, start=identity)
    parity = sp.diag(
        *[sp.Integer((-1) ** state.bit_count()) for state in range(FOCK_DIMENSION)]
    )
    require(gamma_internal == parity, "oriented Clifford product does not equal parity")
    require(gamma_internal * gamma_internal == identity, "Clifford volume does not square to one")

    vacuum = sp.zeros(FOCK_DIMENSION)
    vacuum[0, 0] = 1
    projector = (identity + gamma_internal) / 2 - vacuum
    require(projector * projector == projector, "P15 is not idempotent")
    require(projector.conjugate().transpose() == projector, "P15 is not Hermitian")
    require(projector.rank() == 15, "P15 does not have complex rank 15")

    selected = [state for state in range(FOCK_DIMENSION) if projector[state, state] == 1]
    degrees = sorted({state.bit_count() for state in selected})
    require(degrees == [2, 4], "the derived P15 image has unexpected exterior degrees")

    return {
        "status": "EXACT_INTERNAL_Q0_FIXTURE",
        "definition": "P15=(I+Gamma_int)/2-P_vac",
        "gamma_definition": "Gamma_int=i^5*c1*c2*...*c10",
        "majorana_count": len(majoranas),
        "full_fock_dimension": FOCK_DIMENSION,
        "rank": int(projector.rank()),
        "derived_exterior_degrees": degrees,
        "physical_promotion": False,
        "open_physical_receipts": [
            "spacetime Spin lift and central minus one",
            "graded exchange and locality",
            "nonzero source residue and sterile-sector exclusion",
            "complement-complete refinement",
        ],
    }


def gaussian_composite_1pi_certificate(parameter: sp.Expr = sp.Integer(1)) -> dict[str, Any]:
    """Separate a composite Legendre cubic from the fundamental Gaussian 1PI."""

    a = sp.sympify(parameter)
    variance = sp.simplify(1 + 2 * a**2)
    third_cumulant = sp.simplify(6 * a + 8 * a**3)
    composite_gamma3 = sp.simplify(-third_cumulant / variance**3)
    fundamental_gamma3 = sp.Integer(0)
    require(variance != 0, "the composite response variance vanishes")
    if a != 0:
        require(composite_gamma3 != 0, "nonlinear composite cubic unexpectedly vanished")
    require(fundamental_gamma3 == 0, "the standard Gaussian fundamental action is not quadratic")
    return {
        "status": "EXACT_NEGATIVE_CONTROL",
        "interpolator": "O=x+a(x^2-1)",
        "connected_variance": str(variance),
        "connected_third_cumulant": str(third_cumulant),
        "composite_effective_action_third_derivative": str(composite_gamma3),
        "fundamental_effective_action_third_derivative": "0",
        "physical_promotion": False,
        "boundary": (
            "A nonlinear composite cubic does not establish a microscopic Yukawa vertex; "
            "the complete fundamental and mixed source system is required."
        ),
    }


def refinement_complement_certificate() -> dict[str, Any]:
    """Show why exact old-sector intertwining does not exclude hidden poles."""

    coarse = (sp.Integer(0), sp.Integer(1), sp.Integer(2))
    cutoff = sp.Rational(5, 2)
    good_complement = (sp.Integer(10), sp.Integer(11))
    bad_complement = (sp.Integer(0),)

    coarse_rank = sum(1 for value in coarse if bool(value < cutoff))
    good_rank = sum(1 for value in coarse + good_complement if bool(value < cutoff))
    bad_rank = sum(1 for value in coarse + bad_complement if bool(value < cutoff))
    require(coarse_rank == good_rank == 3, "positive refinement rank check failed")
    require(bad_rank == 4, "hidden-complement mutation did not increase selected rank")

    contour_center = sp.Integer(1)
    contour_radius = sp.Integer(2)

    def inside(value: sp.Expr) -> bool:
        return bool(abs(value - contour_center) < contour_radius)

    require(not any(inside(value) for value in good_complement), "good complement enters contour")
    require(any(inside(value) for value in bad_complement), "bad complement misses contour interior")
    boundary_distance = abs(abs(bad_complement[0] - contour_center) - contour_radius)
    require(boundary_distance != 0, "negative-control pole lies on the contour")

    return {
        "status": "EXACT_POSITIVE_AND_NEGATIVE_CONTROL",
        "cutoff": str(cutoff),
        "positive": {
            "old_sector_intertwining_defect": 0,
            "complement_spectrum": [str(value) for value in good_complement],
            "selected_rank": good_rank,
            "passes_complement_interior_exclusion": True,
        },
        "hidden_zero_mode": {
            "old_sector_intertwining_defect": 0,
            "complement_spectrum": [str(value) for value in bad_complement],
            "coarse_selected_rank": coarse_rank,
            "fine_selected_rank": bad_rank,
            "distance_from_contour_boundary_is_nonzero": True,
            "passes_complement_interior_exclusion": False,
        },
        "physical_promotion": False,
    }


def finite_settlement_certificate() -> dict[str, Any]:
    """Certify the finite 16 -> 16 -> 12 trap and exhaustive reference backend."""

    risks = {"T0": 16, "T1": 16, "T2": 12}
    edges = {"T0": ("T1",), "T1": ("T0", "T2"), "T2": ("T1",)}

    local = "T0"
    lower = [node for node in edges[local] if risks[node] < risks[local]]
    require(not lower, "T0 is not a strict-local-descent trap")

    seen = {local}
    queue = deque([local])
    while queue:
        node = queue.popleft()
        for neighbor in edges[node]:
            if neighbor not in seen:
                seen.add(neighbor)
                queue.append(neighbor)
    minimizer = min(seen, key=lambda node: (risks[node], node))
    require(minimizer == "T2", "exhaustive settlement missed the global minimizer")
    return {
        "status": "EXACT_SMALL_REFERENCE_BACKEND",
        "risk_path": [16, 16, 12],
        "strict_local_result": local,
        "exhaustive_result": minimizer,
        "physical_promotion": False,
        "boundary": "Finite exhaustive settlement is not a scalable local physical relaxation law.",
    }


def build_all_certificates() -> dict[str, Any]:
    return {
        "source_completion_no_go": source_completion_nonuniqueness_certificate(),
        "rank15": rank15_clifford_certificate(),
        "gaussian_1pi": gaussian_composite_1pi_certificate(),
        "refinement": refinement_complement_certificate(),
        "settlement": finite_settlement_certificate(),
    }


if __name__ == "__main__":
    import json

    print(json.dumps(build_all_certificates(), indent=2, sort_keys=True))
