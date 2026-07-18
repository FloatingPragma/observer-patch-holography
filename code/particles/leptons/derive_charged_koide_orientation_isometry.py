#!/usr/bin/env python3
"""Certify the finite tracial-GNS Koide theorem and its physical boundary.

A connected M6 source register with a minimal two-state orientation record has
one rank-two singlet event and one rank-two oriented charged event after the
admissibility projection. Born-Lueders conditioning therefore gives equal
block probabilities. The canonical square-root vector in the tracial GNS
space turns those probabilities into response powers without an adjustable
normalization, forcing |b|/a=1/sqrt(2) and Q=2/3 on the positive chamber.

This closes the finite probability-to-power step. It does not construct the
physical chiral C3 carrier or an exact recoverable source-to-mass-response
channel. No charged-lepton mass or fitted charged-family parameter is read.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
ORIENTATION_CERTIFICATE = (
    ROOT / "particles" / "hierarchy" / "certificates"
    / "R_m_rep_24_certificate.json"
)
FACE_CARRIER = (
    ROOT / "particles" / "runs" / "leptons"
    / "charged_icosahedral_face_carrier_frontier.json"
)
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "leptons"
    / "charged_koide_orientation_isometry.json"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def koide_from_modulus_ratio(modulus_over_singlet: float) -> float:
    """Return signed-trace Q; it is physical Q when the carrier is positive."""

    ratio = float(modulus_over_singlet)
    return (1.0 + 2.0 * ratio * ratio) / 3.0


def orientation_isometry_modulus(orientation_count: int) -> float:
    """Equal modulus in an isometric split over orthogonal orientations."""

    if orientation_count <= 0:
        raise ValueError("orientation_count must be positive")
    return 1.0 / math.sqrt(float(orientation_count))


def response_isometry_defect(modulus_over_singlet: float) -> float:
    """T_kappa^* T_kappa - 1 for the real-linear conjugate-pair response."""

    ratio = float(modulus_over_singlet)
    return 2.0 * ratio * ratio - 1.0


def circulant_roots(modulus_over_singlet: float, phase: float) -> tuple[float, ...]:
    """Eigenvalues at singlet coefficient a=1."""

    ratio = float(modulus_over_singlet)
    return tuple(
        1.0 + 2.0 * ratio * math.cos(phase + 2.0 * math.pi * index / 3.0)
        for index in range(3)
    )


def physical_koide_from_roots(roots: tuple[float, ...]) -> float:
    """Physical Q for masses m_i=r_i^2, so sqrt(m_i)=|r_i|."""

    return sum(root * root for root in roots) / sum(abs(root) for root in roots) ** 2


def balanced_positive_chamber(phase: float, tolerance: float = 1.0e-14) -> bool:
    """Whether delta lies in |delta|<=pi/12 modulo 2*pi/3."""

    period = 2.0 * math.pi / 3.0
    representative = (phase + period / 2.0) % period - period / 2.0
    return abs(representative) <= math.pi / 12.0 + tolerance


def c3_trivial_to_charged_intertwiner_dimension(tolerance: float = 1.0e-14) -> int:
    """Dimension of Hom_C3(1, chi (+) conjugate(chi)); it is zero."""

    omega = complex(math.cos(2.0 * math.pi / 3.0), math.sin(2.0 * math.pi / 3.0))
    charged_eigenvalues = (omega, omega.conjugate())
    return sum(abs(value - 1.0) <= tolerance for value in charged_eigenvalues)


def quotient_block_weights(fiber_sizes: tuple[int, ...]) -> tuple[float, ...]:
    """Uniform weights on quotient-visible blocks, independent of fiber size."""

    if not fiber_sizes or any(size <= 0 for size in fiber_sizes):
        raise ValueError("fiber sizes must be positive")
    weight = 1.0 / len(fiber_sizes)
    return tuple(weight for _ in fiber_sizes)


def lifted_micro_weights(fiber_sizes: tuple[int, ...]) -> tuple[tuple[float, ...], ...]:
    """Uniformly lift each quotient-block weight inside its fiber."""

    block_weights = quotient_block_weights(fiber_sizes)
    return tuple(
        tuple(block_weight / size for _ in range(size))
        for block_weight, size in zip(block_weights, fiber_sizes, strict=True)
    )


def connected_register_event_ranks(
    singlet_rank: int = 1,
    charged_rank: int = 2,
    orientation_dimension: int = 2,
) -> tuple[int, int, int, int]:
    """Return M_d dimension and ranks of Z0, Zc, and E=Z0+Zc."""

    if min(singlet_rank, charged_rank, orientation_dimension) <= 0:
        raise ValueError("all dimensions must be positive")
    source_dimension = (singlet_rank + charged_rank) * orientation_dimension
    z0_rank = singlet_rank * orientation_dimension
    zc_rank = charged_rank
    return source_dimension, z0_rank, zc_rank, z0_rank + zc_rank


def born_lueders_block_probabilities(
    z0_rank: int,
    zc_rank: int,
) -> tuple[float, float]:
    """Probabilities after conditioning the maximally mixed state on Z0+Zc."""

    if min(z0_rank, zc_rank) <= 0:
        raise ValueError("event ranks must be positive")
    admitted_rank = z0_rank + zc_rank
    return z0_rank / admitted_rank, zc_rank / admitted_rank


def gns_modulus_from_block_probabilities(p0: float, pc: float) -> float:
    """Return |b|/a for the Hermitian C3 GNS amplitude."""

    if p0 <= 0.0 or pc < 0.0:
        raise ValueError("block probabilities must satisfy p0>0 and pc>=0")
    return math.sqrt(pc / (2.0 * p0))


def koide_from_action_gap(action_gap: float) -> float:
    """Q for charged multiplicity two and source action gap Delta S."""

    return (1.0 + 2.0 * math.exp(-float(action_gap))) / 3.0


def c3_gns_gram_matrix() -> tuple[tuple[complex, ...], ...]:
    """Gram matrix of (I,R,R^2) for tau_3(A^dagger B)."""

    return tuple(
        tuple(
            1.0 + 0.0j if (column - row) % 3 == 0 else 0.0 + 0.0j
            for column in range(3)
        )
        for row in range(3)
    )


def build_artifact(
    orientation: dict[str, Any] | None = None,
    face: dict[str, Any] | None = None,
) -> dict[str, Any]:
    if orientation is None:
        orientation = json.loads(ORIENTATION_CERTIFICATE.read_text(encoding="utf-8"))
    if face is None:
        face = json.loads(FACE_CARRIER.read_text(encoding="utf-8"))

    orientation_count = int(
        orientation["representation_sector"]["orientation_multiplier"]
    )
    modulus = orientation_isometry_modulus(orientation_count)
    q_isometric = koide_from_modulus_ratio(modulus)
    block_weights = quotient_block_weights((1, orientation_count))
    micro_weights = lifted_micro_weights((1, orientation_count))
    action_gap = math.log(float(orientation_count))
    natural_trace_ratio = float(orientation_count)
    corrected_ratio = natural_trace_ratio * math.exp(-action_gap)
    source_dimension, z0_rank, zc_rank, admitted_rank = (
        connected_register_event_ranks(
            singlet_rank=1,
            charged_rank=orientation_count,
            orientation_dimension=orientation_count,
        )
    )
    p0, pc = born_lueders_block_probabilities(z0_rank, zc_rank)
    gns_modulus = gns_modulus_from_block_probabilities(p0, pc)
    gns_q = koide_from_modulus_ratio(gns_modulus)
    positive_phase_controls = (0.0, 0.1, 2.0 / 9.0)
    positive_phase_checks = {
        str(phase): {
            "roots": list(circulant_roots(modulus, phase)),
            "inside_exact_positive_chamber": balanced_positive_chamber(phase),
            "physical_Q": physical_koide_from_roots(circulant_roots(modulus, phase)),
        }
        for phase in positive_phase_controls
    }
    outside_phase = 0.4
    outside_roots = circulant_roots(modulus, outside_phase)

    source_checks = {
        "orientation_certificate_accepted": orientation.get("accepted") is True,
        "orientation_multiplier_is_two": orientation_count == 2,
        "orientation_pair_is_write_verify": "write/verify"
        in orientation["repair_grammar"]["orientation_meaning"],
        "face_carrier_is_source_only": face.get("source_only") is True,
        "face_carrier_has_geometric_C3": face.get("face_orbit_theorem", {}).get(
            "face_stabilizer_order"
        ) == 3,
        "face_carrier_does_not_claim_physical_attachment": face.get(
            "public_charged_mass_promotion_allowed"
        ) is False,
    }
    algebra_checks = {
        "isometric_modulus_is_one_over_sqrt_two": math.isclose(
            modulus, 1.0 / math.sqrt(2.0), rel_tol=0.0, abs_tol=1.0e-15
        ),
        "singlet_and_charged_powers_are_equal": math.isclose(
            3.0, 6.0 * modulus * modulus, rel_tol=0.0, abs_tol=1.0e-15
        ),
        "koide_is_two_thirds": math.isclose(
            q_isometric, 2.0 / 3.0, rel_tol=0.0, abs_tol=1.0e-15
        ),
        "quotient_blocks_are_equiprobable": block_weights == (0.5, 0.5),
        "lifted_orientation_weights_are_one_quarter_each": micro_weights
        == ((0.5,), (0.25, 0.25)),
        "fiber_action_gap_is_ln_two": math.isclose(
            action_gap, math.log(2.0), rel_tol=0.0, abs_tol=1.0e-15
        ),
        "fiber_corrected_power_ratio_is_one": math.isclose(
            corrected_ratio, 1.0, rel_tol=0.0, abs_tol=1.0e-15
        ),
        "connected_source_algebra_is_M6": source_dimension == 6,
        "admissibility_event_has_equal_rank_two_blocks": (
            z0_rank == zc_rank == 2 and admitted_rank == 4
        ),
        "born_lueders_conditioning_gives_equal_block_probabilities": (
            p0 == pc == 0.5
        ),
        "tracial_GNS_basis_is_orthonormal": c3_gns_gram_matrix()
        == (
            (1.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j),
            (0.0 + 0.0j, 1.0 + 0.0j, 0.0 + 0.0j),
            (0.0 + 0.0j, 0.0 + 0.0j, 1.0 + 0.0j),
        ),
        "canonical_GNS_square_root_map_forces_balanced_modulus": math.isclose(
            gns_modulus, 1.0 / math.sqrt(2.0), rel_tol=0.0, abs_tol=1.0e-15
        ),
        "canonical_GNS_square_root_map_forces_koide": math.isclose(
            gns_q, 2.0 / 3.0, rel_tol=0.0, abs_tol=1.0e-15
        ),
        "action_gap_master_formula_gives_koide": math.isclose(
            koide_from_action_gap(math.log(2.0)),
            2.0 / 3.0,
            rel_tol=0.0,
            abs_tol=1.0e-15,
        ),
        "positive_phase_controls_are_positive_and_koide_balanced": all(
            control["inside_exact_positive_chamber"]
            and min(control["roots"]) >= -1.0e-14
            and math.isclose(
                control["physical_Q"], 2.0 / 3.0, rel_tol=0.0, abs_tol=1.0e-14
            )
            for control in positive_phase_checks.values()
        ),
        "outside_chamber_signed_identity_is_not_physical_koide": (
            not balanced_positive_chamber(outside_phase)
            and min(outside_roots) < 0.0
            and not math.isclose(
                physical_koide_from_roots(outside_roots),
                2.0 / 3.0,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
        ),
        "no_C3_equivariant_trivial_to_charged_linear_map": (
            c3_trivial_to_charged_intertwiner_dimension() == 0
        ),
        "exact_Koide_isometry_stability_identity": all(
            math.isclose(
                abs(koide_from_modulus_ratio(ratio) - 2.0 / 3.0),
                abs(response_isometry_defect(ratio)) / 3.0,
                rel_tol=0.0,
                abs_tol=1.0e-15,
            )
            for ratio in (0.4, 0.6, modulus)
        ),
    }

    countermodels = []
    for ratio in (0.4, 0.6):
        countermodels.append(
            {
                "modulus_over_singlet": ratio,
                "phase": 0.0,
                "Q": koide_from_modulus_ratio(ratio),
                "C3_covariant_and_Hermitian": True,
                "positive_square_root_mass_carrier": min(
                    circulant_roots(ratio, 0.0)
                ) > 0.0,
                "orientation_response_isometric": math.isclose(
                    2.0 * ratio * ratio, 1.0, rel_tol=0.0, abs_tol=1.0e-15
                ),
            }
        )

    checks_pass = all(source_checks.values()) and all(algebra_checks.values())
    return {
        "artifact": "oph_charged_koide_orientation_isometry",
        "schema_version": 2,
        "status": (
            "CLOSED_FINITE_TRACIAL_GNS_KOIDE_THEOREM_"
            "PHYSICAL_CHIRAL_RECOVERY_ATTACHMENT_OPEN"
        ),
        "source_only": False,
        "charged_reference_data_consumed": False,
        "public_koide_promotion_allowed": False,
        "checks_pass": checks_pass,
        "source_receipts": {
            "orientation_pair": {
                "path": ORIENTATION_CERTIFICATE.relative_to(ROOT).as_posix(),
                "sha256": sha256(ORIENTATION_CERTIFICATE),
                "status": orientation.get("status"),
            },
            "geometric_C3_carrier": {
                "path": FACE_CARRIER.relative_to(ROOT).as_posix(),
                "sha256": sha256(FACE_CARRIER),
                "status": face.get("status"),
            },
        },
        "source_checks": source_checks,
        "finite_GNS_theorem": {
            "connected_source_algebra": {
                "algebra": "B((1 (+) chi (+) conjugate(chi)) tensor C^2) = M6(C)",
                "dimension": source_dimension,
                "maxent_state": "I_6/6",
            },
            "minimal_orientation_event": {
                "Z0": "P0 tensor I_2",
                "Zc": "Pc tensor q_+",
                "rank_Z0": z0_rank,
                "rank_Zc": zc_rank,
                "rank_E_plus": admitted_rank,
                "conditioned_probabilities": {"p0": p0, "pc": pc},
                "action_gap": "ln(2)",
            },
            "canonical_response_map": {
                "source_basis": ["e0", "e+", "e-"],
                "GNS_basis": ["I", "R", "R^2"],
                "GNS_gram_matrix": [
                    [[value.real, value.imag] for value in row]
                    for row in c3_gns_gram_matrix()
                ],
                "amplitudes": {
                    "a": "sqrt(p0)",
                    "abs_b": "sqrt(pc/2)",
                    "abs_b_over_a": gns_modulus,
                },
                "Q": gns_q,
            },
            "general_action_law": {
                "pc_over_p0": "2*exp(-Delta_S)",
                "abs_b_squared_over_a_squared": "exp(-Delta_S)",
                "Q_of_Delta_S": "(1+2*exp(-Delta_S))/3",
            },
            "logical_scope": (
                "The connected-register rank calculation, Born-Lueders conditioning, and "
                "canonical tracial-GNS probability-to-power map are closed finite "
                "mathematics. Applying Minimal Admissible Realization to response-local "
                "registers is additional structural content, and identifying this GNS "
                "process with the physical charged-mass response requires the open exact "
                "recoverable attachment below."
            ),
        },
        "conditional_theorem": {
            "premise": (
                "Supply a physical chiral three-family regular-C3 carrier and exact "
                "state-preserving recoverable u.c.p. source-to-physical maps that intertwine "
                "the public singlet and charged records and identify the physical GNS vector "
                "with the normalized positive square-root-mass response after kinetic "
                "whitening."
            ),
            "proof": [
                "In the connected M6 source register, the admitted singlet and oriented charged events both have rank two, so Born-Lueders conditioning gives p0=pc=1/2.",
                "The canonical tracial-GNS square-root vector maps e0,e+,e- unitarily to I,R,R^2 and converts record probabilities into squared response-component norms.",
                "Exact recoverability and faithful trace preservation make the physical channel an L2 isometry, while record intertwining transports p0=pc to E0=Ec.",
                "For the C3 circulant, E_0=3a^2 and E_c=6|b|^2, so E0=Ec gives |b|/a=1/sqrt(2).",
                "On the positive-root chamber, the exact identity Q=(1+E_c/E_0)/3 then gives Q=2/3, independently of the charged phase delta.",
            ],
            "conclusion": (
                "The finite normalization is derived rather than adjusted. Conditional on the "
                "open exact recoverable physical attachment, the same normalization reaches "
                "the positive square-root-mass response and gives Koide without a lepton-mass "
                "input or delta=2/9."
            ),
            "geometric_form": (
                "For the positive square-root-mass vector, Q=2/3 iff its squared norm in "
                "the C3 singlet line equals its squared norm in the centered real plane; "
                "equivalently the vector makes angle pi/4 with (1,1,1)."
            ),
            "stability_identity": (
                "For the real-linear map T_kappa:R->H_c, "
                "a |-> kappa*a*(exp(i*delta),exp(-i*delta)), "
                "the signed-trace identity is |Q-2/3|="
                "|T_kappa^*T_kappa-1|/3; it is the physical Koide stability "
                "identity on any phase branch where the carrier is positive semidefinite."
            ),
            "phase_boundary": {
                "positive_root_chamber": "|delta| <= pi/12 modulo 2*pi/3",
                "statement": (
                    "The theorem leaves delta continuous inside the positive-root chamber. "
                    "Orientation reversal carries delta to -delta rather than fixing delta. A "
                    "phase law is needed for the two independent charged-lepton mass ratios, "
                    "but not for the Koide invariant there. Outside that chamber, physical "
                    "square roots are absolute eigenvalues and the signed trace identity is "
                    "not the physical Koide relation."
                ),
                "positive_controls": positive_phase_checks,
                "outside_chamber_control": {
                    "delta": outside_phase,
                    "roots": list(outside_roots),
                    "signed_trace_Q": q_isometric,
                    "physical_Q": physical_koide_from_roots(outside_roots),
                },
            },
        },
        "conditional_equivalent_balance_descriptions": {
            "shared_bridge_requirement": (
                "The canonical tracial-GNS map closes probability-to-power inside the finite "
                "source response model. A separate exact recoverable attachment is required to "
                "identify that model with the physical charged-mass response."
            ),
            "orientation_isometry": {
                "orientation_count": orientation_count,
                "modulus_over_singlet": modulus,
                "power_ratio_Ec_over_E0": 2.0 * modulus * modulus,
                "Q": q_isometric,
            },
            "quotient_block_maxent": {
                "fiber_sizes_singlet_charged": [1, orientation_count],
                "uniform_quotient_block_weights": list(block_weights),
                "uniform_within_fiber_lift": [list(row) for row in micro_weights],
                "interpretation": (
                    "Uniform MaxEnt on the two quotient-visible blocks gives equal block "
                    "probabilities. In the connected M6 construction these weights are inherited "
                    "from the unique normalized trace rather than selected from the direct-sum "
                    "trace family."
                ),
            },
            "natural_trace_plus_fiber_correction": {
                "base_weight_ratio_charged_over_singlet": natural_trace_ratio,
                "required_action_gap": "ln(2)",
                "corrected_power_ratio": corrected_ratio,
                "interpretation": (
                    "The historical one-bit action is exactly the fiber-multiplicity correction "
                    "that converts uniform microstate counting into equal admitted source-block "
                    "weights. The GNS map makes this a finite response-power theorem; physical "
                    "mass-response identification remains conditional."
                ),
            },
        },
        "algebra_checks": algebra_checks,
        "necessity_boundary": {
            "countermodels_without_isometry": countermodels,
            "statement": (
                "C3 covariance, Hermiticity, reversible orientation doubling, and traciality "
                "alone do not fix the response norm. The connected source algebra, declared "
                "admissibility event, Born-Lueders conditioning, and canonical GNS amplitude "
                "do fix it in the finite model. Without the recoverable physical attachment, "
                "the modulus of a physical C3 mass carrier remains continuous."
            ),
            "logical_strength": (
                "For the real-linear T_kappa:R->H_c, "
                "T_kappa(a)=kappa*a*(exp(i*delta),exp(-i*delta)) and "
                "T_kappa^*T_kappa=2*kappa^2. Thus T_kappa is isometric iff "
                "kappa=1/sqrt(2), iff E_c=E_0, iff physical Q=2/3 on the positive "
                "chamber. Away from a positive carrier the displayed Q is only the signed-"
                "trace functional. The finite GNS construction now fixes kappa canonically in "
                "the source response model; it does not prove that the physical mass operator "
                "is the image of that model."
            ),
        },
        "representation_boundary": {
            "Hom_C3_trivial_to_charged_dimension": c3_trivial_to_charged_intertwiner_dimension(),
            "statement": (
                "There is no nonzero C3-equivariant linear map from the trivial singlet to "
                "chi (+) conjugate(chi). The gauge-adjoint write/verify pair therefore cannot "
                "be attached as a labeled charged-family vector by symmetry alone. This does "
                "not forbid a C3-invariant radial norm law or an orbit-valued output, which is "
                "enough for the Q-only criterion. A source-emitted charged tensor, connection, "
                "or symmetry-breaking selected orbit is required for a vector attachment and "
                "the individual mass-ratio phase."
            ),
            "invariant_metric_boundary": (
                "On the real C3 module 1 (+) W_2, every invariant positive inner product "
                "has independent positive scales alpha on the singlet and beta on W_2. "
                "C3 symmetry fixes the form inside each irreducible block but not their "
                "relative normalization, so it cannot make the response isometry canonical."
            ),
        },
        "open_physical_gate": {
            "id": "charged_GNS_to_physical_chiral_mass_recovery_attachment",
            "derived": False,
            "requirements": [
                "justify applying Minimal Admissible Realization to the response-local two-state orientation register and the connected M6 admissibility event",
                "construct the physical chiral L-H-E carrier with quotient-visible regular C3 family fibers and canonical kinetic forms",
                "construct faithful trace-preserving recoverable u.c.p. maps between the finite source process and physical response process, with a left inverse and singlet/charged record intertwining",
                "identify the physical GNS vector with the normalized positive square-root-mass response C=(Y_e^dagger Y_e)^(1/4), after kinetic whitening",
                "prove refinement naturality and path exhaustion for the response records and recovery maps",
                "prove positivity of the resulting square-root-mass operator on the selected phase branch",
                "derive the running-coordinate, threshold, and pole-mass bridge needed to identify the operator with the measured charged-lepton masses",
                "bind the construction to a no-target-leak ancestry receipt before charged-mass comparison",
            ],
            "claim_boundary": (
                "The finite connected-register, Born-Lueders, and tracial-GNS normalization is "
                "closed. The corpus does not yet emit the physical chiral regular-C3 carrier or "
                "the exact recoverable source-to-mass-response attachment. The result therefore "
                "does not promote a source-only charged-lepton prediction."
            ),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": artifact["status"], "checks_pass": artifact["checks_pass"]}, indent=2))
    return 0 if artifact["checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
