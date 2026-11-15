"""Bounded exact self-readback experiment for the atomic seam-repair law.

This module tests a finite version of the proposed OPH self-readback closure.
It does not fit a hidden transition matrix.  A candidate dynamics is exercised
through an operational event interface.  The public readback derives only
constraint facts from those events.  A separately fixed reconstruction
algorithm then applies the following finite clauses:

* the primitive mismatch locations are the thirty seams of the oriented
  twelve-port carrier;
* a completed event preserves the endpoint total and lands in the nearest
  integer-agreement shell;
* A3 assigns equal weight to the thirty seam attempts and, conditionally, to
  the unresolved two-outcome parity tie.  The sixty directed labels encode
  those local completion placements; they are not a sixty-letter attempt
  alphabet.

The complete candidate transition kernel is held out until scoring.  A
candidate is a fixed point when its generator is a positive rational multiple
of the reconstructed generator.  The multiplier is the common clock rate.

All promoted decisions use integers and :class:`fractions.Fraction`.  The
state probe covers every nonnegative twelve-port load vector with protected
total zero through four. Constructive audits also cover all 1,352,078
nonnegative total-twelve states and all 531,441 signed states in
``{-1,0,1}^12``. The positive fixed-point receipt concerns the one-step
expectation generator. An IID path law is conditional on proposed
free-composition and temporal-completeness clauses because canonical A3 does
not imply Markovity. The experiment does not certify a complete full-algebra
grammar, refinement compatibility, a physical repair law, protected-sector
selection, or universe closure.
"""

from __future__ import annotations

import argparse
import copy
from collections import deque
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
from functools import lru_cache
from itertools import product
from math import comb
from pathlib import Path
from typing import Any, Mapping, Sequence

REPORT_SCHEMA = "oph.bounded_atomic_self_readback_closure.v1"
VERIFICATION_SCHEMA = "oph.bounded_atomic_self_readback_closure_verification.v1"
REFERENCE_REPORT_PATH = (
    Path(__file__).resolve().parents[2]
    / "data"
    / "repair_closure"
    / "bounded_atomic_self_readback_closure_receipt.json"
)
TOTAL_TWELVE_SOURCE_PROJECTION_PATH = (
    Path(__file__).resolve().parents[2]
    / "data"
    / "repair_closure"
    / "record_counting_source_projection.json"
)
TOTAL_TWELVE_SOURCE_PROJECTION_SHA256 = (
    "sha256:ac3942d586724894daa48e977361e287023d16a3b34080bbce3238ed4c80f628"
)

DIRECTED_SEAM_TORSOR_RECEIPT = "FINITE_DIRECTED_SEAM_TORSOR_RECEIPT"
BOUNDED_SELF_READBACK_RECEIPT = (
    "BOUNDED_EXPECTATION_LEVEL_ATOMIC_SELF_READBACK_FIXED_POINT_RECEIPT"
)
BOUNDED_COUPLED_CLOSURE_RECEIPT = (
    "BOUNDED_COUPLED_STATE_GENERATOR_CLOSURE_RECEIPT"
)
CONDITIONAL_FREE_WORD_LAW_RECEIPT = (
    "CONDITIONAL_FREE_EVENT_WORD_IID_LAW_RECEIPT"
)
GLOBAL_POLICY_RECEIPT = "GLOBAL_A1_A3_POLICY_UNIQUENESS_RECEIPT"
FULL_SELF_READBACK_RECEIPT = "FULL_SELF_READING_UNIVERSE_CLOSURE_RECEIPT"
PHYSICAL_REPAIR_RECEIPT = "PHYSICAL_REPAIR_LAW_RECEIPT"

PORT_COUNT = 12
UNDIRECTED_SEAM_COUNT = 30
# The physical attempt alphabet is the thirty undirected seams.  The sixty
# directed labels encode the two local completion placements and are useful
# for exact combinatorics; they are not sixty independent seam attempts.
SEAM_ATTEMPT_COUNT = UNDIRECTED_SEAM_COUNT
DIRECTED_SEAM_COUNT = 60
MAX_PROBE_TOTAL = 4

# Exact source-side oriented face presentation of the base carrier.  The
# triples are cyclically oriented.  This tuple deliberately replaces the
# floating-coordinate mesh builder in the torsor proof.
ORIENTED_BASE_FACES = (
    (0, 11, 5),
    (0, 5, 1),
    (0, 1, 7),
    (0, 7, 10),
    (0, 10, 11),
    (1, 5, 9),
    (5, 11, 4),
    (11, 10, 2),
    (10, 7, 6),
    (7, 1, 8),
    (3, 9, 4),
    (3, 4, 2),
    (3, 2, 6),
    (3, 6, 8),
    (3, 8, 9),
    (4, 9, 5),
    (2, 4, 11),
    (6, 2, 10),
    (8, 6, 7),
    (9, 8, 1),
)

State = tuple[int, ...]
DirectedEdge = tuple[int, int]
TransitionRow = dict[int, Fraction]
TransitionRows = dict[int, tuple[TransitionRow, ...]]


@dataclass(frozen=True)
class CandidateLaw:
    """One member of the frozen adversarial candidate suite."""

    law_id: str
    local_policy: str
    support_policy: str
    clock: Fraction = Fraction(1)

    def __post_init__(self) -> None:
        if not self.law_id or not isinstance(self.law_id, str):
            raise ValueError("law_id must be a nonempty string")
        if self.local_policy not in {
            "directed_balanced",
            "reverse_directed_balanced",
            "nearest_keep_high_side",
            "one_unit_descent",
            "identity",
            "swap",
            "uniform_all_conserved_splits",
            "lower_index_gets_ceiling",
        }:
            raise ValueError(f"unknown local policy: {self.local_policy}")
        if self.support_policy not in {
            "direct_uniform",
            "direct_edge_biased",
            "direct_direction_biased",
            "distance_two_uniform",
            "antipodal_uniform",
            "direct_plus_distance_two_uniform",
        }:
            raise ValueError(f"unknown support policy: {self.support_policy}")
        if self.clock <= 0 or self.clock > 1:
            raise ValueError("the bounded discrete clock must lie in (0,1]")


def exact_reference_edges() -> tuple[tuple[int, int], ...]:
    """Derive the thirty unoriented seams from the exact face tuple."""

    edges = tuple(
        sorted(
            {
                tuple(sorted((face[index], face[(index + 1) % 3])))
                for face in ORIENTED_BASE_FACES
                for index in range(3)
            }
        )
    )
    if len(edges) != UNDIRECTED_SEAM_COUNT:
        raise ValueError("the exact oriented faces do not define thirty seams")
    degrees = [0] * PORT_COUNT
    for left, right in edges:
        degrees[left] += 1
        degrees[right] += 1
    if degrees != [5] * PORT_COUNT:
        raise ValueError("the exact face presentation is not five-regular")
    return edges


def exact_directed_seam_torsor() -> dict[str, Any]:
    """Construct the proper rotation action using incidence and face cycles.

    The 120 incidence automorphisms are enumerated without coordinates.
    Orientation preservation is decided from the cyclic ordering of the
    twenty source faces.  The surviving sixty permutations act simply
    transitively on the sixty directed seams.
    """

    source_faces = ORIENTED_BASE_FACES
    edges = exact_reference_edges()
    if (
        len(source_faces) != 20
        or len({frozenset(face) for face in source_faces}) != 20
        or len(_oriented_face_set(source_faces)) != 20
    ):
        raise ValueError("the exact oriented face presentation has duplicates")
    directed_face_boundaries: dict[DirectedEdge, int] = {}
    for left, middle, right in source_faces:
        for directed in ((left, middle), (middle, right), (right, left)):
            directed_face_boundaries[directed] = (
                directed_face_boundaries.get(directed, 0) + 1
            )
    if any(
        directed_face_boundaries.get((left, right)) != 1
        or directed_face_boundaries.get((right, left)) != 1
        for left, right in edges
    ):
        raise ValueError(
            "the face cycles do not induce opposite directions on every seam"
        )
    adjacency = [[0] * PORT_COUNT for _ in range(PORT_COUNT)]
    for left, right in edges:
        adjacency[left][right] = 1
        adjacency[right][left] = 1
    automorphisms = _incidence_automorphisms(adjacency)
    if len(automorphisms) != 120:
        raise ValueError("the incidence automorphism group does not have order 120")

    oriented_faces = _oriented_face_set(source_faces)
    rotations = tuple(
        permutation
        for permutation in automorphisms
        if _oriented_face_set(
            tuple(
                (
                    permutation[left],
                    permutation[middle],
                    permutation[right],
                )
                for left, middle, right in source_faces
            )
        )
        == oriented_faces
    )
    if len(rotations) != 60:
        raise ValueError("the oriented incidence action does not have order 60")

    rotation_set = frozenset(rotations)
    if tuple(range(PORT_COUNT)) not in rotation_set:
        raise ValueError("the proper rotation action has no identity")
    if any(
        _compose_permutations(left, right) not in rotation_set
        for left in rotations
        for right in rotations
    ):
        raise ValueError("the proper rotation action is not closed")

    directed = tuple(
        sorted(
            directed_edge
            for left, right in edges
            for directed_edge in ((left, right), (right, left))
        )
    )
    reference = directed[0]
    orbit = {
        (permutation[reference[0]], permutation[reference[1]])
        for permutation in rotations
    }
    stabilizer = [
        permutation
        for permutation in rotations
        if (
            permutation[reference[0]],
            permutation[reference[1]],
        )
        == reference
    ]
    if orbit != set(directed) or len(stabilizer) != 1:
        raise ValueError("the proper action is not simply transitive on directed seams")
    if any(
        (permutation[left], permutation[right]) not in set(directed)
        for permutation in rotations
        for left, right in directed
    ):
        raise ValueError("a proper rotation leaves the directed-seam set")

    atom_hash = _sha256_json([list(edge) for edge in directed])
    action_hash = _sha256_json([list(permutation) for permutation in rotations])
    return {
        "port_count": PORT_COUNT,
        "undirected_seam_count": len(edges),
        "directed_event_atom_count": len(directed),
        "incidence_automorphism_count": len(automorphisms),
        "orientation_preserving_rotation_count": len(rotations),
        "reference_directed_seam": list(reference),
        "directed_seam_orbit_size": len(orbit),
        "reference_stabilizer_order": len(stabilizer),
        "action_closed": True,
        "action_simply_transitive": True,
        "construction": (
            "explicit integer oriented-face tuple; incidence automorphisms "
            "filtered by exact cyclic-face preservation"
        ),
        "floating_coordinate_matching_used": False,
        "floating_geometry_source_used": False,
        "exact_oriented_face_count": len(source_faces),
        "exact_oriented_faces_distinct": True,
        "every_seam_has_two_incident_faces": True,
        "incident_face_directions_are_opposite": True,
        "exact_oriented_face_sha256": _sha256_json(
            [list(face) for face in source_faces]
        ),
        "directed_event_atom_sha256": atom_hash,
        "proper_action_sha256": action_hash,
        "_directed_edges": directed,
        "_rotations": rotations,
    }


def enumerate_protected_sector(total: int) -> tuple[State, ...]:
    """Enumerate all nonnegative twelve-port states with a fixed total."""

    if isinstance(total, bool) or not isinstance(total, int) or total < 0:
        raise ValueError("total must be a nonnegative integer")
    return tuple(_weak_compositions(total, PORT_COUNT))


def frozen_candidate_suite() -> tuple[CandidateLaw, ...]:
    """Return the source-frozen target-free adversarial suite."""

    candidates: list[CandidateLaw] = []
    for policy in ("directed_balanced", "reverse_directed_balanced"):
        for numerator, denominator in ((1, 4), (1, 2), (1, 1)):
            candidates.append(
                CandidateLaw(
                    law_id=f"{policy}_clock_{numerator}_{denominator}",
                    local_policy=policy,
                    support_policy="direct_uniform",
                    clock=Fraction(numerator, denominator),
                )
            )
    candidates.extend(
        (
            CandidateLaw(
                "nearest_keep_high_side_control",
                "nearest_keep_high_side",
                "direct_uniform",
            ),
            CandidateLaw(
                "one_unit_descent_control",
                "one_unit_descent",
                "direct_uniform",
            ),
            CandidateLaw(
                "identity_control",
                "identity",
                "direct_uniform",
            ),
            CandidateLaw(
                "swap_control",
                "swap",
                "direct_uniform",
            ),
            CandidateLaw(
                "uniform_all_splits_control",
                "uniform_all_conserved_splits",
                "direct_uniform",
            ),
            CandidateLaw(
                "port_label_bias_control",
                "lower_index_gets_ceiling",
                "direct_uniform",
            ),
            CandidateLaw(
                "edge_schedule_bias_control",
                "directed_balanced",
                "direct_edge_biased",
            ),
            CandidateLaw(
                "orientation_schedule_bias_control",
                "directed_balanced",
                "direct_direction_biased",
            ),
            CandidateLaw(
                "distance_two_control",
                "directed_balanced",
                "distance_two_uniform",
            ),
            CandidateLaw(
                "antipodal_control",
                "directed_balanced",
                "antipodal_uniform",
            ),
            CandidateLaw(
                "mixed_radius_control",
                "directed_balanced",
                "direct_plus_distance_two_uniform",
            ),
        )
    )
    return tuple(candidates)


def public_constraint_readback(candidate: CandidateLaw) -> dict[str, Any]:
    """Return only the operational constraints visible to reconstruction."""

    context = _experiment_context()
    evaluation = _evaluate_candidate(candidate, context)
    return copy.deepcopy(evaluation["constraint_readback"])


def reconstruct_from_constraint_readback(
    readback: Mapping[str, Any],
) -> CandidateLaw | None:
    """Apply the frozen source reconstruction without fitting observed dynamics.

    The exact field allowlist is part of the anti-circularity boundary.
    Candidate names, update parameters, transition probabilities, and
    held-out transition digests are not accepted by this function. Observed
    support and schedule fields remain in the audit schema, but the
    reconstruction decision does not branch on them. A1 fixes the seam
    support and A3 constructs its counting measure source-side. Wrong-support
    and biased-schedule candidates reach held-out scoring and fail there.
    """

    expected_keys = {
        "schema",
        "grammar_manifest_sha256",
        "protected_totals",
        "max_probe_total",
        "complete_state_event_probe",
        "primitive_probe_count",
        "expected_primitive_probe_count",
        "event_atom_support_is_exact_directed_seam_torsor",
        "event_atom_schedule_is_uniform",
        "event_atom_schedule_has_full_support",
        "all_outcomes_preserve_protected_total",
        "all_outcomes_change_only_event_endpoints",
        "all_outcomes_land_in_nearest_agreement_shell",
        "all_neutral_changes_are_typed_reversible_stutters",
        "local_rule_is_presentation_covariant",
        "endpoint_reversal_covariant",
        "target_fields_read",
        "downstream_measurements_read",
    }
    if not isinstance(readback, Mapping) or set(readback) != expected_keys:
        return None
    if (
        readback.get("schema") != "oph.atomic_repair_constraint_readback.v1"
        or readback.get("protected_totals") != list(range(MAX_PROBE_TOTAL + 1))
        or readback.get("max_probe_total") != MAX_PROBE_TOTAL
        or readback.get("primitive_probe_count")
        != readback.get("expected_primitive_probe_count")
        or readback.get("target_fields_read") != []
        or readback.get("downstream_measurements_read") is not False
    ):
        return None
    required_local_semantics = (
        "complete_state_event_probe",
        "all_outcomes_preserve_protected_total",
        "all_outcomes_change_only_event_endpoints",
        "all_outcomes_land_in_nearest_agreement_shell",
        "all_neutral_changes_are_typed_reversible_stutters",
        "local_rule_is_presentation_covariant",
        "endpoint_reversal_covariant",
    )
    if any(
        readback.get(name) is not True
        for name in required_local_semantics
    ):
        return None
    if readback.get("grammar_manifest_sha256") != _experiment_context()[
        "grammar_manifest_sha256"
    ]:
        return None
    return CandidateLaw(
        law_id="source_reconstructed_directed_seam_law",
        local_policy="directed_balanced",
        support_policy="direct_uniform",
        clock=Fraction(1),
    )


def candidate_transition_rows(candidate: CandidateLaw) -> TransitionRows:
    """Return the exact held-out transition rows for focused audits."""

    context = _experiment_context()
    return copy.deepcopy(_evaluate_candidate(candidate, context)["transition_rows"])


def generator_ray_comparison(
    candidate_rows: TransitionRows,
    reference_rows: TransitionRows,
) -> dict[str, Any]:
    """Compare two finite generators modulo one positive rational clock."""

    if set(candidate_rows) != set(reference_rows):
        return {
            "equivalent": False,
            "clock_ratio": None,
            "reason": "protected_sector_manifest_mismatch",
        }
    ratio: Fraction | None = None
    compared_entries = 0
    for total in sorted(reference_rows):
        left_rows = candidate_rows[total]
        right_rows = reference_rows[total]
        if len(left_rows) != len(right_rows):
            return {
                "equivalent": False,
                "clock_ratio": None,
                "reason": "state_manifest_mismatch",
            }
        for source, (left, right) in enumerate(
            zip(left_rows, right_rows, strict=True)
        ):
            targets = set(left) | set(right) | {source}
            for target in sorted(targets):
                left_generator = left.get(target, Fraction()) - Fraction(
                    int(source == target)
                )
                right_generator = right.get(target, Fraction()) - Fraction(
                    int(source == target)
                )
                compared_entries += 1
                if right_generator == 0:
                    if left_generator != 0:
                        return {
                            "equivalent": False,
                            "clock_ratio": None,
                            "reason": "generator_support_mismatch",
                            "compared_entries": compared_entries,
                        }
                    continue
                candidate_ratio = left_generator / right_generator
                if ratio is None:
                    ratio = candidate_ratio
                elif candidate_ratio != ratio:
                    return {
                        "equivalent": False,
                        "clock_ratio": None,
                        "reason": "nonuniform_generator_ratio",
                        "compared_entries": compared_entries,
                    }
    if ratio is None or ratio <= 0:
        return {
            "equivalent": False,
            "clock_ratio": None,
            "reason": "no_positive_generator_ray",
            "compared_entries": compared_entries,
        }
    return {
        "equivalent": True,
        "clock_ratio": str(ratio),
        "reason": "one_positive_rational_generator_ray",
        "compared_entries": compared_entries,
    }


def self_readback_repair_closure_report() -> dict[str, Any]:
    """Build the deterministic bounded certificate."""

    payload = json.loads(_cached_payload_json())
    report = copy.deepcopy(payload)
    report["certificate_payload_sha256"] = _sha256_json(payload)
    return report


def verify_self_readback_repair_closure_report(
    report: Mapping[str, Any],
) -> dict[str, Any]:
    """Fail closed on mutations, malformed data, or scope promotion."""

    reasons: list[str] = []
    try:
        if not isinstance(report, Mapping):
            raise TypeError("report must be a mapping")
        received = copy.deepcopy(dict(report))
        received_hash = received.pop("certificate_payload_sha256", None)
        if received.get("schema") != REPORT_SCHEMA:
            reasons.append("schema_mismatch")
        if received_hash != _sha256_json(received):
            reasons.append("payload_hash_mismatch")
        expected = json.loads(_cached_payload_json())
        if received != expected:
            reasons.append("producer_replay_mismatch")
        for name in (
            BOUNDED_COUPLED_CLOSURE_RECEIPT,
            GLOBAL_POLICY_RECEIPT,
            FULL_SELF_READBACK_RECEIPT,
            PHYSICAL_REPAIR_RECEIPT,
        ):
            if report.get(name) is not False:
                reasons.append("forbidden_scope_promotion")
    except (
        AttributeError,
        TypeError,
        ValueError,
        OverflowError,
        RecursionError,
    ):
        reasons.append("malformed_or_noncanonical_payload")
    return {
        "schema": VERIFICATION_SCHEMA,
        "receipt": not reasons,
        "status": "PASS" if not reasons else "FAIL",
        "reasons": sorted(set(reasons)),
        "producer_replay": True,
        "independent_implementation": False,
        "claim_boundary": (
            "Verification replays the bounded exact producer. It cannot "
            "promote full-algebra grammar completeness, refinement, physical "
            "repair identification, protected-sector selection, or universe "
            "closure."
        ),
    }


def write_self_readback_repair_closure_report(
    output_path: str | Path,
) -> dict[str, Any]:
    """Write the deterministic certificate using LF line endings."""

    report = self_readback_repair_closure_report()
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return report


@lru_cache(maxsize=1)
def _cached_payload_json() -> str:
    payload = _build_payload()
    _require_no_floats(payload)
    return json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def _build_payload() -> dict[str, Any]:
    context = _experiment_context()
    torsor = context["torsor"]
    reconstructed = CandidateLaw(
        "source_reconstructed_directed_seam_law",
        "directed_balanced",
        "direct_uniform",
    )
    reconstructed_evaluation = _evaluate_candidate(reconstructed, context)
    reference_rows = reconstructed_evaluation["transition_rows"]
    sector_rows = _stationary_closure_rows(
        context["states_by_total"],
        reference_rows,
    )
    mean_bridge = _all_probed_mean_bridge(
        context["states_by_total"],
        reference_rows,
    )
    two_event_path_law = _two_event_path_law_check(
        context["states_by_total"][1],
        reference_rows[1],
        exact_reference_edges(),
    )
    total_twelve = _total_twelve_diagnostics(
        tuple(torsor["_directed_edges"])
    )
    signed_cube = _exhaustive_signed_cube_progress_audit()
    free_word_law = _free_event_word_law_certificate(SEAM_ATTEMPT_COUNT)

    candidate_rows: list[dict[str, Any]] = []
    fixed_ray_digests: set[str] = set()
    for candidate in frozen_candidate_suite():
        evaluated = _evaluate_candidate(candidate, context)
        reconstruction = reconstruct_from_constraint_readback(
            evaluated["constraint_readback"]
        )
        if reconstruction is None:
            comparison = {
                "equivalent": False,
                "clock_ratio": None,
                "reason": "constraint_readback_ineligible",
                "compared_entries": 0,
            }
        else:
            comparison = generator_ray_comparison(
                evaluated["transition_rows"],
                reference_rows,
            )
        fixed = bool(reconstruction is not None and comparison["equivalent"])
        ray_digest = (
            _generator_ray_digest(evaluated["transition_rows"])
            if fixed
            else None
        )
        if ray_digest is not None:
            fixed_ray_digests.add(ray_digest)
        constraint = evaluated["constraint_readback"]
        candidate_rows.append(
            {
                "law_id": candidate.law_id,
                "local_policy": candidate.local_policy,
                "support_policy": candidate.support_policy,
                "declared_clock": str(candidate.clock),
                "constraint_reconstruction_eligible": reconstruction is not None,
                "constraint_checks": {
                    key: constraint[key]
                    for key in (
                        "event_atom_support_is_exact_directed_seam_torsor",
                        "event_atom_schedule_is_uniform",
                        "event_atom_schedule_has_full_support",
                        "all_outcomes_preserve_protected_total",
                        "all_outcomes_change_only_event_endpoints",
                        "all_outcomes_land_in_nearest_agreement_shell",
                        "all_neutral_changes_are_typed_reversible_stutters",
                        "local_rule_is_presentation_covariant",
                        "endpoint_reversal_covariant",
                    )
                },
                "primitive_probe_count": constraint["primitive_probe_count"],
                "held_out_transition_sha256": evaluated[
                    "held_out_transition_sha256"
                ],
                "fixed_modulo_clock": fixed,
                "generator_ray_comparison": comparison,
                "normalized_generator_ray_sha256": ray_digest,
            }
        )

    fixed = [row for row in candidate_rows if row["fixed_modulo_clock"]]
    expected_fixed = {
        f"{policy}_clock_{numerator}_{denominator}"
        for policy in ("directed_balanced", "reverse_directed_balanced")
        for numerator, denominator in ((1, 4), (1, 2), (1, 1))
    }
    fixed_ids = {row["law_id"] for row in fixed}
    bounded_one_step_closure = bool(
        fixed_ids == expected_fixed
        and len(fixed_ray_digests) == 1
        and all(row["unique_stationary_state"] for row in sector_rows)
        and mean_bridge["all_probed_states_exact_identity_verified"]
    )

    by_id = {row["law_id"]: row for row in candidate_rows}
    one_unit_evaluation = _evaluate_candidate(
        CandidateLaw(
            "one_unit_descent_control",
            "one_unit_descent",
            "direct_uniform",
        ),
        context,
    )
    one_unit_stationary = _stationary_closure_rows(
        context["states_by_total"],
        one_unit_evaluation["transition_rows"],
    )

    public_keys = sorted(
        _evaluate_candidate(
            CandidateLaw(
                "public_boundary_probe",
                "directed_balanced",
                "direct_uniform",
            ),
            context,
        )["constraint_readback"]
    )
    payload: dict[str, Any] = {
        "schema": REPORT_SCHEMA,
        "scope": (
            "bounded_nonnegative_atomic_loads_on_one_oriented_twelve_port_carrier"
        ),
        "source_inputs": {
            "carrier": "oriented_twelve_port_icosahedral_boundary_incidence",
            "protected_total_sectors": list(range(MAX_PROBE_TOTAL + 1)),
            "state_count_by_total": {
                str(total): len(states)
                for total, states in context["states_by_total"].items()
            },
            "total_state_count": sum(
                len(states) for states in context["states_by_total"].values()
            ),
            "laboratory_data_used": False,
            "measured_particle_or_cosmology_values_used": False,
            "downstream_target_used": False,
        },
        "directed_seam_torsor": {
            key: value
            for key, value in torsor.items()
            if not key.startswith("_")
        },
        "axiom_clause_specialization": {
            "basis_status": "proposed_a1r_a2r_specialization_not_adopted",
            "canonical_three_axiom_derivation": False,
            "a1_primitive_constraint_locations": (
                "the thirty undirected shared-seam mismatch constraints"
            ),
            "a2_reconciliation_and_completion_atoms": (
                "preserve the endpoint total and minimize integer endpoint "
                "disagreement; the two orientations of each seam label the "
                "two odd-remainder completions"
            ),
            "odd_total_pathwise_a2_equalizer_retraction": False,
            "integer_lift_scope": (
                "exact microhistory and expectation-level lift; the real "
                "pair-average channel is the exact A2 equalizer retraction"
            ),
            "macro_event_to_atomic_record_history_bridge_certified": False,
            "a3_schedule": (
                "maximum randomness gives probability 1/30 per seam and "
                "conditionally 1/2 per unresolved odd tie; the free path "
                "alphabet has thirty seam attempts, while sixty directed "
                "labels encode local completion placements; an IID attempt "
                "law additionally consumes the proposed free-word and "
                "temporal-completeness clauses"
            ),
            "common_clock_selected": False,
            "full_a1_repair_grammar_certified": False,
        },
        "anti_circular_reconstruction": {
            "constraint_readback_schema": (
                "oph.atomic_repair_constraint_readback.v1"
            ),
            "constraint_readback_schema_fields": public_keys,
            "reconstruction_decision_fields": [
                "schema",
                "grammar_manifest_sha256",
                "protected_totals",
                "max_probe_total",
                "complete_state_event_probe",
                "primitive_probe_count",
                "expected_primitive_probe_count",
                "all_outcomes_preserve_protected_total",
                "all_outcomes_change_only_event_endpoints",
                "all_outcomes_land_in_nearest_agreement_shell",
                "all_neutral_changes_are_typed_reversible_stutters",
                "local_rule_is_presentation_covariant",
                "endpoint_reversal_covariant",
                "target_fields_read",
                "downstream_measurements_read",
            ],
            "candidate_law_id_visible_to_reconstructor": False,
            "candidate_update_parameters_visible_to_reconstructor": False,
            "candidate_transition_probabilities_visible_to_reconstructor": False,
            "held_out_transition_digest_visible_to_reconstructor": False,
            "observed_event_support_used_to_construct_law": False,
            "observed_schedule_weights_used_to_construct_law": False,
            "reconstruction_algorithm": (
                "nearest agreement shell plus uniform counting measure on "
                "the source-fixed exact directed-seam torsor"
            ),
            "scoring": (
                "held-out exact continuous-time generator Q=P-I compared "
                "with reconstruction up to one positive rational rate"
            ),
            "clock_boundary": (
                "lazy discrete kernels at different declared rates are not "
                "identified as equal one-tick channels; only their generator "
                "rays, and hence exp(tQ) after time rescaling, are quotiented"
            ),
        },
        "exhaustive_probe": {
            "max_protected_total": MAX_PROBE_TOTAL,
            "state_count": sum(
                len(states) for states in context["states_by_total"].values()
            ),
            "canonical_primitive_state_event_probe_count": (
                sum(
                    len(states)
                    for states in context["states_by_total"].values()
                )
                * DIRECTED_SEAM_COUNT
            ),
            "physical_seam_attempts_per_state": SEAM_ATTEMPT_COUNT,
            "directed_completion_labels_per_state": DIRECTED_SEAM_COUNT,
            "arithmetic": "integers_and_exact_rational_probabilities",
            "random_sampling_used_for_receipt": False,
        },
        "coupled_state_generator_closure": {
            "sector_rows": sector_rows,
            "terminal_shell_rule": (
                "for S=12q+r, each terminal load is q or q+1 and exactly "
                "r ports carry q+1"
            ),
            "general_integer_energy_argument": {
                "energy": "H(x)=sum_i x_i^2",
                "fixed_total_pair_balancing_drop": (
                    "for endpoint gap d>=2, nearest balancing lowers H by "
                    "a positive integer"
                ),
                "neutral_motion": (
                    "for endpoint gap one, the two directed events give a "
                    "wait and a reversible swap"
                ),
                "terminal_shell": (
                    "Euclidean division S=12q+r gives twelve loads in "
                    "{q,q+1}, with r upper loads"
                ),
                "formal_local_integer_theorems_present": True,
                "formal_global_hitting_theorem_present": False,
                "general_signed_state_progress_argument_present": True,
                "maximum_progress_word_length_on_this_carrier": 3,
                "reachable_coordinate_box_is_finite": True,
                "iid_full_support_implies_almost_sure_shell_hitting": True,
                "full_all_total_state_enumeration_performed": False,
            },
            "bounded_shell_dynamics": (
                "lazy_symmetric_exclusion_process_on_the_carrier"
            ),
            "a3_state": "uniform_on_each_fixed_total_terminal_shell",
            "protected_total_sector_selected": False,
        },
        "bounded_signed_state_control": signed_cube,
        "distinguished_total_twelve_sector": total_twelve,
        "exact_conditional_mean_bridge": mean_bridge,
        "two_event_schedule_closure": two_event_path_law,
        "conditional_free_event_word_law": free_word_law,
        "candidate_suite": {
            "frozen_candidate_count": len(candidate_rows),
            "candidate_fixed_point_probe_totals": list(
                range(MAX_PROBE_TOTAL + 1)
            ),
            "total_twelve_candidate_conformance_exhausted": False,
            "complete_grammar_of_every_a1_a3_repair": False,
            "rows": candidate_rows,
            "fixed_candidate_ids": sorted(fixed_ids),
            "fixed_candidate_count": len(fixed),
            "fixed_generator_ray_class_count": len(fixed_ray_digests),
        },
        "adversarial_controls": {
            "nearest_shell_but_nonmaximal_mixing": {
                "law_id": "nearest_keep_high_side_control",
                "constraint_reconstruction_eligible": by_id[
                    "nearest_keep_high_side_control"
                ]["constraint_reconstruction_eligible"],
                "fixed_modulo_clock": by_id[
                    "nearest_keep_high_side_control"
                ]["fixed_modulo_clock"],
                "lesson": (
                    "landing in the nearest shell does not make reconstruction "
                    "an identity; A3 mixing remains discriminating"
                ),
            },
            "one_unit_rule": {
                "law_id": "one_unit_descent_control",
                "constraint_reconstruction_eligible": by_id[
                    "one_unit_descent_control"
                ]["constraint_reconstruction_eligible"],
                "fixed_modulo_clock": by_id[
                    "one_unit_descent_control"
                ]["fixed_modulo_clock"],
                "absorbing_state_count_by_total": {
                    str(row["protected_total"]): row[
                        "absorbing_state_count"
                    ]
                    for row in one_unit_stationary
                },
                "canonical_closed_class_size_by_total": {
                    str(row["protected_total"]): row["closed_class_size"]
                    for row in sector_rows
                },
            },
            "biased_schedule_reaches_held_out_scoring": by_id[
                "edge_schedule_bias_control"
            ]["constraint_reconstruction_eligible"],
            "biased_schedule_rejected_by_held_out_generator": not by_id[
                "edge_schedule_bias_control"
            ]["fixed_modulo_clock"],
            "radius_two_reaches_held_out_scoring": by_id[
                "distance_two_control"
            ]["constraint_reconstruction_eligible"],
            "radius_two_rejected_by_held_out_generator": not by_id[
                "distance_two_control"
            ]["fixed_modulo_clock"],
            "reverse_direction_state_kernel_same_ray_after_event_label_marginalization": by_id[
                "reverse_directed_balanced_clock_1_1"
            ]["fixed_modulo_clock"],
            "clock_quarter_same_ray": by_id[
                "directed_balanced_clock_1_4"
            ]["fixed_modulo_clock"],
        },
        DIRECTED_SEAM_TORSOR_RECEIPT: bool(
            torsor["action_simply_transitive"]
            and torsor["floating_coordinate_matching_used"] is False
        ),
        BOUNDED_SELF_READBACK_RECEIPT: bounded_one_step_closure,
        BOUNDED_COUPLED_CLOSURE_RECEIPT: False,
        CONDITIONAL_FREE_WORD_LAW_RECEIPT: bool(
            free_word_law["a3_optimizer_is_unique_uniform_word_law"]
            and free_word_law["uniform_word_law_factorizes_for_every_length"]
            and free_word_law["finite_word_family_is_prefix_consistent"]
            and two_event_path_law["iid_path_kernel_equals_P_squared"]
            and two_event_path_law[
                "repeat_same_control_differs_from_P_squared"
            ]
        ),
        GLOBAL_POLICY_RECEIPT: False,
        FULL_SELF_READBACK_RECEIPT: False,
        PHYSICAL_REPAIR_RECEIPT: False,
        "status": (
            "BOUNDED_EXPECTATION_LEVEL_REPAIR_FIXED_POINT_UNIQUE_MODULO_CLOCK_IN_"
            "THE_FROZEN_ADVERSARIAL_SUITE"
            if bounded_one_step_closure
            else "BOUNDED_ATOMIC_REPAIR_CLOSURE_NOT_ATTAINED"
        ),
        "claim_boundary": (
            "The exact expectation-level result covers nonnegative atomic "
            "loads with protected "
            "totals zero through four and the complete nonnegative total-12 "
            "sector on one finite oriented carrier. It "
            "shows that a predeclared noncircular reconstruction selects one "
            "one-step mean-generator-ray class inside the frozen adversarial suite "
            "and that "
            "the reconstructed dynamics has one uniform stationary state per "
            "tested protected sector. The total-12 progress audit establishes "
            "the unique all-one closure class on all 1,352,078 nonnegative "
            "compositions. Odd-total paths end in a nearest-balanced shell "
            "rather than the exact agreement equalizer, so the integer "
            "microhistory is not itself an A2-R conditional expectation. It "
            "does not prove that seams are the complete "
            "repair grammar of the full OPH algebra, select a protected "
            "sector, derive temporal independence from the canonical axioms, "
            "establish refinement compatibility, identify a physical clock or "
            "repair law, or select a universe. IID path closure is exact only "
            "conditional on the proposed free-composition and complete "
            "temporal-constraint clauses."
        ),
    }
    return payload


@lru_cache(maxsize=1)
def _experiment_context() -> dict[str, Any]:
    torsor = exact_directed_seam_torsor()
    states_by_total = {
        total: enumerate_protected_sector(total)
        for total in range(MAX_PROBE_TOTAL + 1)
    }
    manifest = {
        "carrier": "oriented_twelve_port_icosahedral_boundary_incidence",
        "port_count": PORT_COUNT,
        "event_atom_count": DIRECTED_SEAM_COUNT,
        "directed_event_atom_sha256": torsor["directed_event_atom_sha256"],
        "proper_action_sha256": torsor["proper_action_sha256"],
        "protected_totals": list(range(MAX_PROBE_TOTAL + 1)),
        "a2_integer_equalizer": "nearest_balanced_endpoint_shell",
    }
    return {
        "torsor": torsor,
        "states_by_total": states_by_total,
        "grammar_manifest_sha256": _sha256_json(manifest),
    }


@lru_cache(maxsize=None)
def _evaluate_candidate_cached(
    candidate: CandidateLaw,
) -> dict[str, Any]:
    return _evaluate_candidate_uncached(candidate, _experiment_context())


def _evaluate_candidate(
    candidate: CandidateLaw,
    context: Mapping[str, Any],
) -> dict[str, Any]:
    if context is _experiment_context():
        return _evaluate_candidate_cached(candidate)
    return _evaluate_candidate_uncached(candidate, context)


def _evaluate_candidate_uncached(
    candidate: CandidateLaw,
    context: Mapping[str, Any],
) -> dict[str, Any]:
    torsor = context["torsor"]
    expected_events = tuple(torsor["_directed_edges"])
    event_weights = _event_weights(candidate.support_policy, expected_events)
    event_atoms = tuple(sorted(event_weights))
    expected_event_set = set(expected_events)

    support_exact = set(event_atoms) == expected_event_set
    schedule_full_support = all(weight > 0 for weight in event_weights.values())
    schedule_uniform = bool(
        event_weights
        and len(set(event_weights.values())) == 1
        and sum(event_weights.values(), Fraction()) == 1
    )

    states_by_total = context["states_by_total"]
    rows_by_total: TransitionRows = {}
    total_preserved = True
    endpoint_local = True
    nearest_shell = True
    neutral_typed = True
    primitive_probe_count = 0
    local_signatures: dict[
        tuple[int, int],
        set[tuple[tuple[int, int, str], ...]],
    ] = {}
    for total, states in states_by_total.items():
        state_index = {state: index for index, state in enumerate(states)}
        transition_rows: list[TransitionRow] = []
        for source_index, state in enumerate(states):
            row: TransitionRow = {
                source_index: Fraction(1) - candidate.clock
            }
            for event, event_weight in event_weights.items():
                primitive_probe_count += 1
                left, right = event
                outcomes = _local_outcomes(candidate.local_policy, state, event)
                if sum(outcomes.values(), Fraction()) != 1:
                    raise ValueError("a local candidate row is not normalized")
                local_signature = tuple(
                    sorted(
                        (
                            target[left],
                            target[right],
                            str(probability),
                        )
                        for target, probability in outcomes.items()
                    )
                )
                local_signatures.setdefault(
                    (state[left], state[right]),
                    set(),
                ).add(local_signature)
                for target, probability in outcomes.items():
                    if sum(target) != total:
                        total_preserved = False
                    if any(
                        target[port] != state[port]
                        for port in range(PORT_COUNT)
                        if port not in event
                    ):
                        endpoint_local = False
                    if abs(target[left] - target[right]) > 1:
                        nearest_shell = False
                    before_energy = state[left] ** 2 + state[right] ** 2
                    after_energy = target[left] ** 2 + target[right] ** 2
                    if after_energy > before_energy:
                        neutral_typed = False
                    elif after_energy == before_energy and target != state:
                        swapped = list(state)
                        swapped[left], swapped[right] = (
                            state[right],
                            state[left],
                        )
                        if (
                            target != tuple(swapped)
                            or abs(state[left] - state[right]) != 1
                        ):
                            neutral_typed = False
                    target_index = state_index.get(target)
                    if target_index is None:
                        raise ValueError("a candidate leaves the protected sector")
                    row[target_index] = row.get(
                        target_index,
                        Fraction(),
                    ) + candidate.clock * event_weight * probability
            row = {
                target: probability
                for target, probability in row.items()
                if probability != 0
            }
            if (
                any(probability < 0 for probability in row.values())
                or sum(row.values(), Fraction()) != 1
            ):
                raise ValueError("a candidate transition row is not stochastic")
            transition_rows.append(row)
        rows_by_total[total] = tuple(transition_rows)

    presentation_covariant = all(
        len(signatures) == 1 for signatures in local_signatures.values()
    )
    reversal_covariant = True
    for left, right in sorted(
        {tuple(sorted(event)) for event in event_atoms}
    ):
        forward_event = (left, right)
        reverse_event = (right, left)
        pair_weight = event_weights.get(
            forward_event,
            Fraction(),
        ) + event_weights.get(reverse_event, Fraction())
        if pair_weight == 0:
            reversal_covariant = False
            continue
        for left_value in range(MAX_PROBE_TOTAL + 1):
            for right_value in range(MAX_PROBE_TOTAL - left_value + 1):
                state = [0] * PORT_COUNT
                state[left] = left_value
                state[right] = right_value
                swapped_state = list(state)
                swapped_state[left], swapped_state[right] = (
                    right_value,
                    left_value,
                )

                forward_distribution: dict[tuple[int, int], Fraction] = {}
                swapped_distribution: dict[tuple[int, int], Fraction] = {}
                for event in (forward_event, reverse_event):
                    event_probability = event_weights.get(event, Fraction())
                    if event_probability == 0:
                        continue
                    for target, probability in _local_outcomes(
                        candidate.local_policy,
                        tuple(state),
                        event,
                    ).items():
                        pair = (target[left], target[right])
                        forward_distribution[pair] = (
                            forward_distribution.get(pair, Fraction())
                            + event_probability * probability / pair_weight
                        )
                    for target, probability in _local_outcomes(
                        candidate.local_policy,
                        tuple(swapped_state),
                        event,
                    ).items():
                        pair = (target[right], target[left])
                        swapped_distribution[pair] = (
                            swapped_distribution.get(pair, Fraction())
                            + event_probability * probability / pair_weight
                        )
                if forward_distribution != swapped_distribution:
                    reversal_covariant = False

    expected_probe_count = sum(
        len(states) for states in states_by_total.values()
    ) * len(event_atoms)
    constraint_readback = {
        "schema": "oph.atomic_repair_constraint_readback.v1",
        "grammar_manifest_sha256": context["grammar_manifest_sha256"],
        "protected_totals": list(range(MAX_PROBE_TOTAL + 1)),
        "max_probe_total": MAX_PROBE_TOTAL,
        "complete_state_event_probe": (
            primitive_probe_count == expected_probe_count
        ),
        "primitive_probe_count": primitive_probe_count,
        "expected_primitive_probe_count": expected_probe_count,
        "event_atom_support_is_exact_directed_seam_torsor": support_exact,
        "event_atom_schedule_is_uniform": schedule_uniform,
        "event_atom_schedule_has_full_support": schedule_full_support,
        "all_outcomes_preserve_protected_total": total_preserved,
        "all_outcomes_change_only_event_endpoints": endpoint_local,
        "all_outcomes_land_in_nearest_agreement_shell": nearest_shell,
        "all_neutral_changes_are_typed_reversible_stutters": neutral_typed,
        "local_rule_is_presentation_covariant": presentation_covariant,
        "endpoint_reversal_covariant": reversal_covariant,
        "target_fields_read": [],
        "downstream_measurements_read": False,
    }
    return {
        "constraint_readback": constraint_readback,
        "transition_rows": rows_by_total,
        "held_out_transition_sha256": _transition_rows_sha256(
            states_by_total,
            rows_by_total,
        ),
    }


def _event_weights(
    support_policy: str,
    direct_events: Sequence[DirectedEdge],
) -> dict[DirectedEdge, Fraction]:
    direct = tuple(direct_events)
    undirected = frozenset(
        tuple(sorted(edge))
        for edge in direct
    )
    distance_two = tuple(
        sorted(
            (left, right)
            for left in range(PORT_COUNT)
            for right in range(PORT_COUNT)
            if left != right
            and _graph_distance(left, right, undirected) == 2
        )
    )
    antipodal = tuple(
        sorted(
            (left, right)
            for left in range(PORT_COUNT)
            for right in range(PORT_COUNT)
            if left != right
            and _graph_distance(left, right, undirected) == 3
        )
    )
    if len(distance_two) != 60 or len(antipodal) != 12:
        raise ValueError("the carrier distance-event census drifted")
    if support_policy == "direct_uniform":
        return _uniform_weights(direct)
    if support_policy == "direct_edge_biased":
        first = tuple(sorted(direct[0]))
        raw = {
            edge: 2 if tuple(sorted(edge)) == first else 1
            for edge in direct
        }
        return _normalized_integer_weights(raw)
    if support_policy == "direct_direction_biased":
        raw = {
            edge: 2 if edge == direct[0] else 1
            for edge in direct
        }
        return _normalized_integer_weights(raw)
    if support_policy == "distance_two_uniform":
        return _uniform_weights(distance_two)
    if support_policy == "antipodal_uniform":
        return _uniform_weights(antipodal)
    if support_policy == "direct_plus_distance_two_uniform":
        return _uniform_weights(tuple(sorted((*direct, *distance_two))))
    raise ValueError(f"unknown support policy: {support_policy}")


def _local_outcomes(
    policy: str,
    state: State,
    event: DirectedEdge,
) -> dict[State, Fraction]:
    left, right = event
    left_value = state[left]
    right_value = state[right]
    total = left_value + right_value
    lower = total // 2
    upper = total - lower

    pairs: dict[tuple[int, int], Fraction]
    if policy == "directed_balanced":
        pairs = {(lower, upper): Fraction(1)}
    elif policy == "reverse_directed_balanced":
        pairs = {(upper, lower): Fraction(1)}
    elif policy == "nearest_keep_high_side":
        if left_value > right_value:
            pairs = {(upper, lower): Fraction(1)}
        elif right_value > left_value:
            pairs = {(lower, upper): Fraction(1)}
        else:
            pairs = {(left_value, right_value): Fraction(1)}
    elif policy == "one_unit_descent":
        if left_value - right_value >= 2:
            pairs = {(left_value - 1, right_value + 1): Fraction(1)}
        elif right_value - left_value >= 2:
            pairs = {(left_value + 1, right_value - 1): Fraction(1)}
        else:
            pairs = {(left_value, right_value): Fraction(1)}
    elif policy == "identity":
        pairs = {(left_value, right_value): Fraction(1)}
    elif policy == "swap":
        pairs = {(right_value, left_value): Fraction(1)}
    elif policy == "uniform_all_conserved_splits":
        pairs = {
            (value, total - value): Fraction(1, total + 1)
            for value in range(total + 1)
        }
    elif policy == "lower_index_gets_ceiling":
        if left < right:
            pairs = {(upper, lower): Fraction(1)}
        else:
            pairs = {(lower, upper): Fraction(1)}
    else:
        raise ValueError(f"unknown local policy: {policy}")

    outcomes: dict[State, Fraction] = {}
    for (new_left, new_right), probability in pairs.items():
        target = list(state)
        target[left] = new_left
        target[right] = new_right
        frozen = tuple(target)
        outcomes[frozen] = outcomes.get(frozen, Fraction()) + probability
    return outcomes


def _stationary_closure_rows(
    states_by_total: Mapping[int, Sequence[State]],
    rows_by_total: TransitionRows,
) -> list[dict[str, Any]]:
    reports: list[dict[str, Any]] = []
    for total in sorted(states_by_total):
        states = tuple(states_by_total[total])
        rows = rows_by_total[total]
        adjacency = [
            tuple(sorted(target for target, probability in row.items() if probability > 0))
            for row in rows
        ]
        components = _strongly_connected_components(adjacency)
        component_of = {
            vertex: component_index
            for component_index, component in enumerate(components)
            for vertex in component
        }
        closed = [
            component
            for component_index, component in enumerate(components)
            if not any(
                component_of[target] != component_index
                for source in component
                for target in adjacency[source]
            )
        ]
        balanced = {
            index
            for index, state in enumerate(states)
            if max(state) - min(state) <= 1
        }
        one_closed_balanced = bool(
            len(closed) == 1 and set(closed[0]) == balanced
        )
        closed_set = set(closed[0]) if len(closed) == 1 else set()
        column_sums = {
            target: sum(
                (
                    rows[source].get(target, Fraction())
                    for source in closed_set
                ),
                Fraction(),
            )
            for target in closed_set
        }
        doubly_stochastic = bool(
            closed_set
            and all(value == 1 for value in column_sums.values())
        )
        reverse = [[] for _ in states]
        for source, targets in enumerate(adjacency):
            for target in targets:
                reverse[target].append(source)
        reachable = _reachable_from_sources(reverse, closed_set)
        energy_nonincrease = all(
            sum(value * value for value in states[target])
            <= sum(value * value for value in states[source])
            for source, row in enumerate(rows)
            for target, probability in row.items()
            if probability > 0
        )
        absorbing = sum(
            set(row) == {source} and row[source] == 1
            for source, row in enumerate(rows)
        )
        reports.append(
            {
                "protected_total": total,
                "state_count": len(states),
                "expected_state_count": comb(total + PORT_COUNT - 1, PORT_COUNT - 1),
                "strong_component_count": len(components),
                "closed_component_count": len(closed),
                "closed_class_size": len(closed_set),
                "expected_balanced_shell_size": comb(PORT_COUNT, total),
                "closed_class_is_exact_balanced_shell": one_closed_balanced,
                "all_states_reach_closed_class": len(reachable) == len(states),
                "energy_nonincreasing_on_every_positive_transition": (
                    energy_nonincrease
                ),
                "closed_class_exactly_doubly_stochastic": doubly_stochastic,
                "unique_stationary_state": bool(
                    one_closed_balanced
                    and len(reachable) == len(states)
                    and doubly_stochastic
                ),
                "stationary_state": (
                    f"uniform_weight_1/{len(closed_set)}_on_balanced_shell"
                    if closed_set
                    else None
                ),
                "absorbing_state_count": absorbing,
            }
        )
    return reports


def _one_atom_mean_bridge(
    states: Sequence[State],
    rows: Sequence[TransitionRow],
) -> dict[str, Any]:
    if len(states) != PORT_COUNT:
        raise ValueError("the one-atom sector does not have twelve states")
    state_to_port = {
        index: state.index(1)
        for index, state in enumerate(states)
    }
    port_to_state = {
        port: state_index
        for state_index, port in state_to_port.items()
    }
    matrix = [
        [Fraction() for _ in range(PORT_COUNT)]
        for _ in range(PORT_COUNT)
    ]
    for source_state, row in enumerate(rows):
        source_port = state_to_port[source_state]
        for target_state, probability in row.items():
            target_port = state_to_port[target_state]
            matrix[source_port][target_port] += probability

    laplacian = _graph_laplacian(PORT_COUNT, exact_reference_edges())
    expected = [
        [
            Fraction(int(left == right))
            - Fraction(laplacian[left][right], 60)
            for right in range(PORT_COUNT)
        ]
        for left in range(PORT_COUNT)
    ]
    exact = matrix == expected and len(port_to_state) == PORT_COUNT
    return {
        "sector_total": 1,
        "state_count": len(states),
        "identity": "P_one_atom = I - L_icosahedron/60",
        "exact_identity_verified": exact,
        "one_atom_generator": "-L_icosahedron/60",
        "physical_time_scale_selected": False,
        "unitary_response_identified": False,
        "claim_boundary": (
            "This is the mean Markov repair kernel. Sharing the carrier "
            "Laplacian with a reversible response does not derive a factor i "
            "or identify unitary physical propagation."
        ),
    }


def _all_probed_mean_bridge(
    states_by_total: Mapping[int, Sequence[State]],
    rows_by_total: TransitionRows,
) -> dict[str, Any]:
    """Verify the exact conditional mean on every exhaustively probed state.

    The full integer kernel is nonlinear. Pairing each directed completion
    with its reverse nevertheless makes its first moment linear. This check
    evaluates that statement from the held-out transition rows rather than
    assuming it from the local formula.
    """

    laplacian = _graph_laplacian(PORT_COUNT, exact_reference_edges())
    states_checked = 0
    coordinate_checks = 0
    for total in sorted(states_by_total):
        states = tuple(states_by_total[total])
        rows = rows_by_total[total]
        if len(states) != len(rows):
            raise ValueError("the mean bridge state and row manifests differ")
        for source_index, (state, row) in enumerate(
            zip(states, rows, strict=True)
        ):
            observed = tuple(
                sum(
                    (
                        probability * states[target_index][coordinate]
                        for target_index, probability in row.items()
                    ),
                    Fraction(),
                )
                for coordinate in range(PORT_COUNT)
            )
            expected = tuple(
                Fraction(state[coordinate])
                - Fraction(
                    sum(
                        laplacian[coordinate][other] * state[other]
                        for other in range(PORT_COUNT)
                    ),
                    60,
                )
                for coordinate in range(PORT_COUNT)
            )
            if observed != expected:
                raise ValueError(
                    "the held-out integer conditional mean is not "
                    "the icosahedral Laplacian step"
                )
            if sum(row.values(), Fraction()) != 1 or source_index not in range(
                len(states)
            ):
                raise ValueError("the held-out mean row is malformed")
            states_checked += 1
            coordinate_checks += PORT_COUNT

    one_atom = _one_atom_mean_bridge(
        states_by_total[1],
        rows_by_total[1],
    )
    return {
        "protected_totals_checked": sorted(states_by_total),
        "states_checked": states_checked,
        "coordinate_expectations_checked": coordinate_checks,
        "identity": "E[X_next | X=x] = (I - L_icosahedron/60) x",
        "all_probed_states_exact_identity_verified": True,
        "full_integer_transition_kernel_is_linear": False,
        "one_atom_restriction": one_atom,
        "physical_time_scale_selected": False,
        "unitary_response_identified": False,
        "claim_boundary": (
            "The first moment is the scalar Laplacian repair step on every "
            "exhaustively probed state. The held-out event kernel remains "
            "nonlinear and is not identified with unitary propagation."
        ),
    }


def _compose_transition_rows(
    left: Sequence[TransitionRow],
    right: Sequence[TransitionRow],
) -> tuple[TransitionRow, ...]:
    if len(left) != len(right):
        raise ValueError("transition kernels have different state counts")
    output: list[TransitionRow] = []
    for source, middle_row in enumerate(left):
        row: TransitionRow = {}
        for middle, first_probability in middle_row.items():
            for target, second_probability in right[middle].items():
                row[target] = (
                    row.get(target, Fraction())
                    + first_probability * second_probability
                )
        row = {
            target: probability
            for target, probability in row.items()
            if probability
        }
        if sum(row.values(), Fraction()) != 1 or source >= len(left):
            raise ValueError("a composed transition row is not stochastic")
        output.append(row)
    return tuple(output)


def _free_event_word_law_certificate(event_count: int) -> dict[str, Any]:
    """Certify the conditional product law for a free finite event grammar.

    This is an algebraic consequence of a *proposed* source clause: words of
    every finite length are freely composable unless the complete temporal
    constraint grammar records a coupling. It is not inferred from the
    canonical A3 statement, which explicitly does not imply Markovity.
    """

    if event_count != SEAM_ATTEMPT_COUNT:
        raise ValueError("the free-word certificate requires thirty seam attempts")
    checked_lengths = tuple(range(5))
    rows: list[dict[str, Any]] = []
    for length in checked_lengths:
        word_count = event_count**length
        probability = Fraction(1, word_count)
        factorized = Fraction(1, event_count) ** length
        if probability != factorized:
            raise ValueError("uniform word counting failed to factorize")
        rows.append(
            {
                "word_length": length,
                "word_count": word_count,
                "uniform_probability_per_word": str(probability),
                "product_probability_per_word": str(factorized),
                "factorization_exact": True,
            }
        )
    prefix_checks = 0
    for longer in checked_lengths:
        for shorter in range(longer + 1):
            continuation_count = event_count ** (longer - shorter)
            marginal = (
                continuation_count * Fraction(1, event_count**longer)
            )
            if marginal != Fraction(1, event_count**shorter):
                raise ValueError("uniform word laws are not prefix-consistent")
            prefix_checks += 1
    return {
        "status": "conditional_exact_theorem",
        "event_alphabet_size": event_count,
        "source_clause": (
            "primitive seam-readback and reconciliation attempts form a free "
            "finite word grammar and remain enabled at every state, including "
            "attempts whose completion is a wait, unless the complete "
            "A2-visible temporal grammar supplies a coupling"
        ),
        "reference": (
            "counting measure on every finite word space, constructed before "
            "candidate dynamics or downstream observations are exposed"
        ),
        "a3_information_projection": (
            "on the full probability simplex over M^k, relative entropy to "
            "the faithful uniform counting reference is nonnegative and "
            "vanishes only at that reference; equivalently KKT stationarity "
            "and strict convexity give the unique optimizer Unif(M^k)"
        ),
        "a3_optimizer_is_unique_uniform_word_law": True,
        "general_identity": (
            "|M^k|=|M|^k and Unif(M^k)(m_1,...,m_k)="
            "product_j Unif(M)(m_j)"
        ),
        "uniform_word_law_factorizes_for_every_length": True,
        "prefix_consistency_identity": (
            "|M|^(l-k) * |M|^(-l) = |M|^(-k)"
        ),
        "finite_prefix_consistency_checks": prefix_checks,
        "finite_word_family_is_prefix_consistent": True,
        "bernoulli_process_follows_from_consistent_finite_marginals": True,
        "conditional_next_attempt_probability": "1/30",
        "finite_length_rows": rows,
        "repeat_same_control_support_at_length_two": event_count,
        "free_word_support_at_length_two": event_count**2,
        "canonical_a3_alone_implies_markovity": False,
        "proposed_a1r_a2r_temporal_clauses_required": True,
        "alphabet_counts_attempts_not_only_successful_repairs": True,
        "conditioning_on_committed_or_strict_repairs_is_uniform": False,
        "claim_boundary": (
            "The factorization theorem is exact after the free-composition "
            "and complete temporal-constraint clauses are supplied. It does "
            "not promote IID scheduling under the canonical three-axiom text."
        ),
    }


def _balanced_seam_attempt_outcomes(
    state: State,
    seam: tuple[int, int],
) -> dict[State, Fraction]:
    """Conditional A3 tie law for one undirected seam attempt.

    The schedule chooses the seam. The local completion then has one outcome
    for an even endpoint total and two endpoint-swapped outcomes with equal
    weight for an odd total.
    """

    left, right = tuple(sorted(seam))
    if (left, right) not in set(exact_reference_edges()):
        raise ValueError("a seam attempt is outside the source carrier")
    outcomes: dict[State, Fraction] = {}
    for directed in ((left, right), (right, left)):
        for target, probability in _local_outcomes(
            "directed_balanced",
            state,
            directed,
        ).items():
            outcomes[target] = (
                outcomes.get(target, Fraction())
                + Fraction(1, 2) * probability
            )
    if sum(outcomes.values(), Fraction()) != 1:
        raise ValueError("the conditional tie kernel is not normalized")
    return outcomes


def _two_event_path_law_check(
    states: Sequence[State],
    one_step_rows: Sequence[TransitionRow],
    seam_attempts: Sequence[tuple[int, int]],
) -> dict[str, Any]:
    """Separate an IID A3 path law from a correlated uniform-marginal control."""

    states = tuple(states)
    seams = tuple(tuple(sorted(seam)) for seam in seam_attempts)
    if (
        len(states) != PORT_COUNT
        or len(seams) != SEAM_ATTEMPT_COUNT
        or len(set(seams)) != SEAM_ATTEMPT_COUNT
    ):
        raise ValueError("the two-event audit requires the S=1 seam alphabet")
    state_index = {state: index for index, state in enumerate(states)}
    iid_rows: list[TransitionRow] = []
    repeated_rows: list[TransitionRow] = []
    for state in states:
        iid_weights: TransitionRow = {}
        repeated_weights: TransitionRow = {}
        for first in seams:
            for middle, first_outcome_probability in (
                _balanced_seam_attempt_outcomes(state, first).items()
            ):
                for repeated, second_outcome_probability in (
                    _balanced_seam_attempt_outcomes(middle, first).items()
                ):
                    repeated_index = state_index[repeated]
                    repeated_weights[repeated_index] = (
                        repeated_weights.get(repeated_index, Fraction())
                        + Fraction(1, SEAM_ATTEMPT_COUNT)
                        * first_outcome_probability
                        * second_outcome_probability
                    )
                for second in seams:
                    for target, second_outcome_probability in (
                        _balanced_seam_attempt_outcomes(middle, second).items()
                    ):
                        target_index = state_index[target]
                        iid_weights[target_index] = (
                            iid_weights.get(target_index, Fraction())
                            + Fraction(1, SEAM_ATTEMPT_COUNT**2)
                            * first_outcome_probability
                            * second_outcome_probability
                        )
        iid_rows.append(
            {
                target: probability
                for target, probability in iid_weights.items()
                if probability
            }
        )
        repeated_rows.append(
            {
                target: probability
                for target, probability in repeated_weights.items()
                if probability
            }
        )

    squared = _compose_transition_rows(one_step_rows, one_step_rows)
    iid = tuple(iid_rows)
    repeated = tuple(repeated_rows)
    if iid != squared:
        raise ValueError("the uniform two-event path law is not P squared")
    differing_entries = sum(
        repeated[source].get(target, Fraction())
        != squared[source].get(target, Fraction())
        for source in range(len(states))
        for target in range(len(states))
    )
    if not differing_entries:
        raise ValueError("the correlated control was not distinguished")
    return {
        "seam_attempt_count": SEAM_ATTEMPT_COUNT,
        "directed_completion_label_count": DIRECTED_SEAM_COUNT,
        "iid_ordered_attempt_pair_count": SEAM_ATTEMPT_COUNT**2,
        "iid_probability_per_ordered_attempt_pair": "1/900",
        "one_attempt_marginal_probability": "1/30",
        "odd_tie_probability_per_completion": "1/2",
        "iid_path_kernel_equals_P_squared": True,
        "repeat_same_control_has_same_one_attempt_marginal": True,
        "repeat_same_control_differs_from_P_squared": True,
        "repeat_same_differing_matrix_entries": differing_entries,
        "iid_two_step_sha256": _transition_rows_sha256(
            {1: states},
            {1: iid},
        ),
        "repeat_same_two_step_sha256": _transition_rows_sha256(
            {1: states},
            {1: repeated},
        ),
        "claim_boundary": (
            "Uniform one-event frequencies do not establish an IID schedule. "
            "The exact two-event table rejects a repeat-the-same-event "
            "control with identical one-event marginals."
        ),
    }


def _construct_integer_progress_witness(
    state: State,
    directed_events: frozenset[DirectedEdge],
    adjacency: Sequence[Sequence[int]],
) -> tuple[DirectedEdge, ...]:
    """Construct a short strict-progress word on an arbitrary integer sector.

    The empty word means that all coordinates differ by at most one, which is
    exactly the balanced shell for the state's conserved total. Otherwise a
    maximum is transported through unit-difference seams until one completed
    balance strictly lowers the quadratic energy.
    """

    if (
        len(state) != PORT_COUNT
        or any(isinstance(value, bool) or not isinstance(value, int) for value in state)
    ):
        raise ValueError("the progress witness requires twelve integer loads")
    if max(state) - min(state) <= 1:
        return ()
    edges = exact_reference_edges()
    for left, right in edges:
        if abs(state[left] - state[right]) >= 2:
            return (
                (left, right)
                if state[left] > state[right]
                else (right, left),
            )

    maximum = max(state)
    sources = [index for index, value in enumerate(state) if value == maximum]
    goals = {
        index
        for index, value in enumerate(state)
        if value <= maximum - 2
    }
    if not goals:
        raise ValueError("a non-shell integer state has no separated loads")
    queue: deque[int] = deque(sources)
    predecessor: dict[int, int | None] = {
        source: None for source in sources
    }
    endpoint: int | None = None
    while queue:
        current = queue.popleft()
        if current in goals:
            endpoint = current
            break
        for neighbor in adjacency[current]:
            if neighbor not in predecessor:
                predecessor[neighbor] = current
                queue.append(neighbor)
    if endpoint is None:
        raise ValueError("the carrier failed the connected progress search")

    path: list[int] = []
    current: int | None = endpoint
    while current is not None:
        path.append(current)
        current = predecessor[current]
    path.reverse()
    events = tuple(zip(path, path[1:]))
    if not events or any(event not in directed_events for event in events):
        raise ValueError("the progress path left the directed seam grammar")
    return events


@lru_cache(maxsize=1)
def _exhaustive_total_twelve_progress_audit() -> dict[str, Any]:
    """Check all 1,352,078 nonnegative total-12 states exactly."""

    edges = exact_reference_edges()
    directed = frozenset(
        event
        for left, right in edges
        for event in ((left, right), (right, left))
    )
    adjacency: list[list[int]] = [[] for _ in range(PORT_COUNT)]
    for left, right in edges:
        adjacency[left].append(right)
        adjacency[right].append(left)
    adjacency = [sorted(row) for row in adjacency]

    state_count = 0
    minimum_count = 0
    old_absorbing_count = 0
    witness_length_counts: dict[int, int] = {}
    maximum_witness_length = 0
    witness_digest = hashlib.sha256()
    for state in _weak_compositions(12, PORT_COUNT):
        state_count += 1
        old_absorbing = all(
            abs(state[left] - state[right]) <= 1
            for left, right in edges
        )
        old_absorbing_count += int(old_absorbing)
        witness = _construct_integer_progress_witness(
            state,
            directed,
            adjacency,
        )
        witness_length_counts[len(witness)] = (
            witness_length_counts.get(len(witness), 0) + 1
        )
        maximum_witness_length = max(maximum_witness_length, len(witness))
        before_energy = sum(value * value for value in state)
        current = state
        for offset, event in enumerate(witness):
            outcomes = _local_outcomes("directed_balanced", current, event)
            if len(outcomes) != 1:
                raise ValueError("a progress event is not deterministic")
            following = next(iter(outcomes))
            following_energy = sum(value * value for value in following)
            if offset + 1 < len(witness) and following_energy != before_energy:
                raise ValueError("a pre-descent transport changed the energy")
            current = following
        after_energy = sum(value * value for value in current)
        if not witness:
            if state != (1,) * PORT_COUNT or before_energy != 12:
                raise ValueError("a nonminimum state lacks a progress witness")
            minimum_count += 1
        elif after_energy >= before_energy:
            raise ValueError("a progress witness did not lower the energy")
        witness_digest.update(
            (
                ",".join(map(str, state))
                + "|"
                + ";".join(f"{left}>{right}" for left, right in witness)
                + "\n"
            ).encode("ascii")
        )

    if state_count != comb(23, 11) or minimum_count != 1:
        raise ValueError("the exhaustive total-12 census drifted")
    return {
        "states_checked": state_count,
        "unique_minimum_state_count": minimum_count,
        "nonminimum_states_with_strict_descent_path": state_count - minimum_count,
        "maximum_events_before_strict_descent": maximum_witness_length,
        "witness_length_counts": {
            str(length): count
            for length, count in sorted(witness_length_counts.items())
        },
        "all_nonminimum_states_have_strict_descent_path": True,
        "old_one_unit_absorbing_state_count": old_absorbing_count,
        "old_one_unit_nonglobal_absorbing_state_count": (
            old_absorbing_count - minimum_count
        ),
        "canonical_unique_closed_class_follows_from_full_support": True,
        "witness_manifest_sha256": "sha256:" + witness_digest.hexdigest(),
    }


@lru_cache(maxsize=1)
def _exhaustive_signed_cube_progress_audit() -> dict[str, Any]:
    """Exercise the general progress construction on all ``{-1,0,1}^12``.

    This is a bounded signed-state control, not the proof of the unbounded
    statement. The proof uses the same explicit construction and finiteness of
    the coordinate box preserved by pair balancing.
    """

    edges = exact_reference_edges()
    directed = frozenset(
        event
        for left, right in edges
        for event in ((left, right), (right, left))
    )
    adjacency: list[list[int]] = [[] for _ in range(PORT_COUNT)]
    for left, right in edges:
        adjacency[left].append(right)
        adjacency[right].append(left)
    adjacency = [sorted(row) for row in adjacency]

    states_checked = 0
    shell_states = 0
    strict_progress_states = 0
    maximum_witness_length = 0
    witness_length_counts: dict[int, int] = {}
    totals: set[int] = set()
    witness_digest = hashlib.sha256()
    for state in product((-1, 0, 1), repeat=PORT_COUNT):
        state = tuple(state)
        states_checked += 1
        totals.add(sum(state))
        witness = _construct_integer_progress_witness(
            state,
            directed,
            adjacency,
        )
        witness_length_counts[len(witness)] = (
            witness_length_counts.get(len(witness), 0) + 1
        )
        maximum_witness_length = max(maximum_witness_length, len(witness))
        before_energy = sum(value * value for value in state)
        current = state
        for offset, event in enumerate(witness):
            following = next(
                iter(_local_outcomes("directed_balanced", current, event))
            )
            following_energy = sum(value * value for value in following)
            if offset + 1 < len(witness) and following_energy != before_energy:
                raise ValueError("signed pre-descent transport changed energy")
            current = following
        after_energy = sum(value * value for value in current)
        if witness:
            if after_energy >= before_energy:
                raise ValueError("signed progress word did not lower energy")
            strict_progress_states += 1
        else:
            if max(state) - min(state) > 1:
                raise ValueError("a signed non-shell state lacks progress")
            shell_states += 1
        witness_digest.update(
            (
                ",".join(map(str, state))
                + "|"
                + ";".join(f"{left}>{right}" for left, right in witness)
                + "\n"
            ).encode("ascii")
        )

    if states_checked != 3**PORT_COUNT or totals != set(range(-12, 13)):
        raise ValueError("the signed-cube state census drifted")
    return {
        "coordinate_alphabet": [-1, 0, 1],
        "states_checked": states_checked,
        "protected_totals_covered": sorted(totals),
        "balanced_shell_states": shell_states,
        "non_shell_states_with_strict_progress_word": strict_progress_states,
        "maximum_events_before_strict_descent": maximum_witness_length,
        "witness_length_counts": {
            str(length): count
            for length, count in sorted(witness_length_counts.items())
        },
        "all_non_shell_states_have_strict_progress_word": True,
        "witness_manifest_sha256": "sha256:" + witness_digest.hexdigest(),
        "claim_boundary": (
            "This exhaustive signed cube is a bounded control. The unbounded "
            "fixed-total result follows from the explicit maximum-transport "
            "construction and the finite coordinate box preserved by every "
            "nearest-balanced event."
        ),
    }


@lru_cache(maxsize=1)
def _total_twelve_source_projection() -> dict[str, Any]:
    try:
        projection = json.loads(
            TOTAL_TWELVE_SOURCE_PROJECTION_PATH.read_text(encoding="utf-8")
        )
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError("the total-twelve source projection is unreadable") from exc
    if (
        not isinstance(projection, dict)
        or projection.get("schema")
        != "oph.record_counting_source_projection.v1"
        or projection.get("source_issue") != 628
        or projection.get("bounded_exit") != "exact_named_realization"
        or projection.get("full_pile", {}).get("total") != 12
        or projection.get("full_pile", {}).get("initial_counts")
        != [12] + [0] * 11
        or projection.get("termination_lyapunov")
        != "load square V(N) = sum_i N_i^2"
        or projection.get("termination_lyapunov_upstream_json_pointer")
        != "/mechanism/termination_lyapunov"
        or _sha256_json(projection)
        != TOTAL_TWELVE_SOURCE_PROJECTION_SHA256
    ):
        raise ValueError("the total-twelve source projection failed closed")
    return projection


def _total_twelve_diagnostics(
    directed_events: Sequence[DirectedEdge],
) -> dict[str, Any]:
    """Replay exact deterministic event streams in the sourced total-12 sector.

    The constructive progress audit exhausts all ``C(23,11)`` states. The
    trajectories separately exercise long event streams from adversarial
    initial states. A small integer linear-congruential stream chooses only
    source-side directed event atoms; it uses no floating randomness.
    """

    if len(directed_events) != DIRECTED_SEAM_COUNT:
        raise ValueError("the total-twelve replay requires sixty event atoms")
    source_projection = _total_twelve_source_projection()
    initial_states = (
        (12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        (6, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0),
        (5, 3, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0),
        (2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0),
        (3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0),
    )
    seeds = (1, 17, 101, 4099, 65537)
    target = (1,) * PORT_COUNT
    rows: list[dict[str, Any]] = []
    for seed, initial in zip(seeds, initial_states, strict=True):
        state = initial
        generator_state = seed
        strict_drops = 0
        neutral_swaps = 0
        waits = 0
        event_digest = hashlib.sha256()
        for step in range(1, 200_001):
            generator_state = (
                1_664_525 * generator_state + 1_013_904_223
            ) % (2**32)
            event = directed_events[generator_state % len(directed_events)]
            outcomes = _local_outcomes(
                "directed_balanced",
                state,
                event,
            )
            if len(outcomes) != 1:
                raise ValueError("a directed atomic event must be deterministic")
            next_state = next(iter(outcomes))
            before_energy = sum(value * value for value in state)
            after_energy = sum(value * value for value in next_state)
            if after_energy < before_energy:
                strict_drops += 1
            elif next_state != state:
                neutral_swaps += 1
            else:
                waits += 1
            event_digest.update(f"{event[0]}>{event[1]}\n".encode("ascii"))
            state = next_state
            if state == target:
                break
        else:
            raise ValueError("a seeded total-twelve replay did not settle")
        rows.append(
            {
                "seed": seed,
                "initial_state": list(initial),
                "initial_energy": sum(value * value for value in initial),
                "steps_to_all_ones": step,
                "strict_energy_drops": strict_drops,
                "neutral_reversible_swaps": neutral_swaps,
                "wait_events": waits,
                "final_state": list(state),
                "final_energy": sum(value * value for value in state),
                "event_prefix_sha256": "sha256:" + event_digest.hexdigest(),
                "settled_to_all_ones": state == target,
            }
        )
    all_settled = all(row["settled_to_all_ones"] for row in rows)
    exhaustive = _exhaustive_total_twelve_progress_audit()
    return {
        "source_context": (
            "the pinned integer port-counting source projection fixes twelve "
            "atomic +1 writes at p00; this campaign tests the nonnegative "
            "total-twelve sector containing that full pile"
        ),
        "source_projection": {
            "path": str(
                TOTAL_TWELVE_SOURCE_PROJECTION_PATH.relative_to(
                    Path(__file__).resolve().parents[2]
                )
            ),
            "canonical_json_sha256": TOTAL_TWELVE_SOURCE_PROJECTION_SHA256,
            "upstream_manifest_sha256": source_projection[
                "upstream_manifest_sha256"
            ],
            "upstream_file_sha256": source_projection[
                "upstream_file_sha256"
            ],
            "source_issue": source_projection["source_issue"],
            "source_generation": source_projection["full_pile"][
                "source_generation"
            ],
            "termination_lyapunov": source_projection[
                "termination_lyapunov"
            ],
            "termination_lyapunov_upstream_json_pointer": source_projection[
                "termination_lyapunov_upstream_json_pointer"
            ],
            "verified": True,
        },
        "source_packet_imported_or_hash_verified_here": True,
        "protected_total": 12,
        "full_nonnegative_composition_count": comb(23, 11),
        "euclidean_division": "12 = 12*1 + 0",
        "minimum_energy": 12,
        "next_energy_floor": 14,
        "minimum_shell_cardinality": 1,
        "unique_minimum_state": [1] * PORT_COUNT,
        "unique_minimum_follows_from_integer_convexity": True,
        "trajectory_arithmetic": "integer_only_lcg_event_schedule",
        "trajectory_rows": rows,
        "all_seeded_trajectories_settled_to_unique_minimum": all_settled,
        "exhaustive_progress_audit": exhaustive,
        "reconstructed_law_progress_only": True,
        "candidate_conformance_exhausted_on_total_twelve": False,
        "all_nonnegative_states_reach_unique_minimum_almost_surely_under_uniform_iid_attempt_law": True,
        "almost_sure_hitting_schedule_premise": (
            "the proposed free-word and temporal-completeness clauses select "
            "the uniform IID thirty-seam attempt law"
        ),
        "exhaustive_total_twelve_transition_graph_built": False,
        "claim_boundary": (
            "The total and quadratic minimizer are supplied by the separate "
            "declared signed-integer counting lane, pinned through a local "
            "source projection. The complete-balancing macro-event is not "
            "identified with that lane's atomic one-unit record move. The "
            "exact progress audit covers all 1,352,078 states in the "
            "nonnegative invariant sector. Extension to every signed state "
            "uses the separate general finite-energy argument."
        ),
    }


def _strongly_connected_components(
    adjacency: Sequence[Sequence[int]],
) -> list[tuple[int, ...]]:
    vertex_count = len(adjacency)
    reverse = [[] for _ in range(vertex_count)]
    for source, targets in enumerate(adjacency):
        for target in targets:
            reverse[target].append(source)

    seen: set[int] = set()
    finish_order: list[int] = []
    for root in range(vertex_count):
        if root in seen:
            continue
        seen.add(root)
        stack: list[tuple[int, int]] = [(root, 0)]
        while stack:
            vertex, offset = stack[-1]
            if offset < len(adjacency[vertex]):
                target = adjacency[vertex][offset]
                stack[-1] = (vertex, offset + 1)
                if target not in seen:
                    seen.add(target)
                    stack.append((target, 0))
            else:
                finish_order.append(vertex)
                stack.pop()

    assigned: set[int] = set()
    components: list[tuple[int, ...]] = []
    for root in reversed(finish_order):
        if root in assigned:
            continue
        assigned.add(root)
        component: list[int] = []
        stack = [root]
        while stack:
            vertex = stack.pop()
            component.append(vertex)
            for target in reverse[vertex]:
                if target not in assigned:
                    assigned.add(target)
                    stack.append(target)
        components.append(tuple(sorted(component)))
    return components


def _reachable_from_sources(
    adjacency: Sequence[Sequence[int]],
    sources: set[int],
) -> set[int]:
    reached = set(sources)
    stack = list(sources)
    while stack:
        vertex = stack.pop()
        for target in adjacency[vertex]:
            if target not in reached:
                reached.add(target)
                stack.append(target)
    return reached


def _generator_ray_digest(rows_by_total: TransitionRows) -> str:
    pivot: Fraction | None = None
    entries: list[tuple[int, int, int, Fraction]] = []
    for total in sorted(rows_by_total):
        for source, row in enumerate(rows_by_total[total]):
            for target in sorted(set(row) | {source}):
                value = row.get(target, Fraction()) - Fraction(
                    int(source == target)
                )
                if value != 0:
                    entries.append((total, source, target, value))
                    if pivot is None and source != target and value > 0:
                        pivot = value
    if pivot is None:
        raise ValueError("a zero generator has no positive ray")
    normalized = [
        [total, source, target, str(value / pivot)]
        for total, source, target, value in entries
    ]
    return _sha256_json(normalized)


def _transition_rows_sha256(
    states_by_total: Mapping[int, Sequence[State]],
    rows_by_total: TransitionRows,
) -> str:
    digest = hashlib.sha256()
    for total in sorted(rows_by_total):
        states = states_by_total[total]
        for source, row in enumerate(rows_by_total[total]):
            for target, probability in sorted(row.items()):
                digest.update(
                    (
                        f"{total}|{','.join(map(str, states[source]))}|"
                        f"{','.join(map(str, states[target]))}|"
                        f"{probability.numerator}/{probability.denominator}\n"
                    ).encode("ascii")
                )
    return "sha256:" + digest.hexdigest()


def _incidence_automorphisms(
    adjacency: Sequence[Sequence[int]],
) -> tuple[tuple[int, ...], ...]:
    vertex_count = len(adjacency)
    assignment: list[int | None] = [None] * vertex_count
    used = [False] * vertex_count
    results: list[tuple[int, ...]] = []

    def consistent(vertex: int, image: int) -> bool:
        return all(
            adjacency[vertex][other]
            == adjacency[image][assignment[other]]  # type: ignore[index]
            for other in range(vertex)
        )

    def search(vertex: int) -> None:
        if vertex == vertex_count:
            results.append(tuple(assignment))  # type: ignore[arg-type]
            return
        for image in range(vertex_count):
            if used[image] or not consistent(vertex, image):
                continue
            assignment[vertex] = image
            used[image] = True
            search(vertex + 1)
            assignment[vertex] = None
            used[image] = False

    search(0)
    return tuple(results)


def _graph_laplacian(
    vertex_count: int,
    edges: Sequence[tuple[int, int]],
) -> tuple[tuple[int, ...], ...]:
    matrix = [[0] * vertex_count for _ in range(vertex_count)]
    for left, right in edges:
        if left == right or not (
            0 <= left < vertex_count and 0 <= right < vertex_count
        ):
            raise ValueError("invalid exact carrier edge")
        matrix[left][left] += 1
        matrix[right][right] += 1
        matrix[left][right] -= 1
        matrix[right][left] -= 1
    return tuple(tuple(row) for row in matrix)


def _oriented_face_set(
    faces: Sequence[Sequence[int]],
) -> frozenset[tuple[int, int, int]]:
    normalized = []
    for face in faces:
        left, middle, right = (int(value) for value in face)
        normalized.append(
            min(
                (left, middle, right),
                (middle, right, left),
                (right, left, middle),
            )
        )
    return frozenset(normalized)


def _compose_permutations(
    left: Sequence[int],
    right: Sequence[int],
) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(len(left)))


def _weak_compositions(total: int, slots: int) -> Sequence[State]:
    if slots == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for remainder in _weak_compositions(total - first, slots - 1):
            yield (first, *remainder)


def _uniform_weights(
    events: Sequence[DirectedEdge],
) -> dict[DirectedEdge, Fraction]:
    if not events or len(set(events)) != len(events):
        raise ValueError("event atoms must be a nonempty unique sequence")
    return {event: Fraction(1, len(events)) for event in events}


def _normalized_integer_weights(
    weights: Mapping[DirectedEdge, int],
) -> dict[DirectedEdge, Fraction]:
    total = sum(weights.values())
    if total <= 0:
        raise ValueError("event weights need positive total")
    return {
        event: Fraction(value, total)
        for event, value in weights.items()
    }


@lru_cache(maxsize=None)
def _graph_distance(
    source: int,
    target: int,
    undirected_edges: frozenset[tuple[int, int]] | set[tuple[int, int]],
) -> int:
    edge_set = set(undirected_edges)
    reached = {source}
    frontier = {source}
    distance = 0
    while frontier:
        if target in frontier:
            return distance
        next_frontier: set[int] = set()
        for vertex in frontier:
            for left, right in edge_set:
                if left == vertex and right not in reached:
                    next_frontier.add(right)
                elif right == vertex and left not in reached:
                    next_frontier.add(left)
        reached.update(next_frontier)
        frontier = next_frontier
        distance += 1
    raise ValueError("carrier graph is disconnected")


def _sha256_json(value: Any) -> str:
    return "sha256:" + hashlib.sha256(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
    ).hexdigest()


def _require_no_floats(value: Any) -> None:
    if isinstance(value, float):
        raise ValueError("receipt payloads may not contain floats")
    if isinstance(value, Mapping):
        for key, item in value.items():
            _require_no_floats(key)
            _require_no_floats(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            _require_no_floats(item)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--verify", type=Path)
    args = parser.parse_args(argv)
    if args.verify is not None:
        try:
            submitted = json.loads(args.verify.read_text(encoding="utf-8"))
            result = verify_self_readback_repair_closure_report(submitted)
        except (OSError, json.JSONDecodeError, TypeError, ValueError):
            result = {
                "schema": VERIFICATION_SCHEMA,
                "receipt": False,
                "status": "FAIL",
                "reasons": ["unreadable_or_malformed_report"],
            }
        if args.output is not None:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(
                json.dumps(result, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
                newline="\n",
            )
        else:
            print(json.dumps(result, indent=2, sort_keys=True))
        return 0 if result["receipt"] else 1

    report = self_readback_repair_closure_report()
    if args.output is not None:
        write_self_readback_repair_closure_report(args.output)
    else:
        print(json.dumps(report, indent=2, sort_keys=True))
    return 0


__all__ = [
    "BOUNDED_COUPLED_CLOSURE_RECEIPT",
    "BOUNDED_SELF_READBACK_RECEIPT",
    "CandidateLaw",
    "CONDITIONAL_FREE_WORD_LAW_RECEIPT",
    "DIRECTED_SEAM_TORSOR_RECEIPT",
    "FULL_SELF_READBACK_RECEIPT",
    "GLOBAL_POLICY_RECEIPT",
    "PHYSICAL_REPAIR_RECEIPT",
    "REFERENCE_REPORT_PATH",
    "candidate_transition_rows",
    "enumerate_protected_sector",
    "exact_reference_edges",
    "exact_directed_seam_torsor",
    "frozen_candidate_suite",
    "generator_ray_comparison",
    "main",
    "public_constraint_readback",
    "reconstruct_from_constraint_readback",
    "self_readback_repair_closure_report",
    "verify_self_readback_repair_closure_report",
    "write_self_readback_repair_closure_report",
]


if __name__ == "__main__":
    raise SystemExit(main())
