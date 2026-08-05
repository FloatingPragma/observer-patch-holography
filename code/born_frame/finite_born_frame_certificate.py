#!/usr/bin/env python3
"""Exact finite Born-frame certificate for issue #687 (B11).

The finite audit has two deliberately separate branches.

* The twelve central port atoms in the source manifest form one classical
  twelve-outcome context.  Their weights determine a state on the
  commutative algebra C^12, but do not determine an ambient M_12 density
  matrix.
* The binary-icosahedral spinor adapter supplies twelve qubit rank-one
  projectors, arranged into six antipodal binary contexts.  This adapter is
  a mathematical construction from source-derived geometry, not a
  source-produced public quantum instrument.

All linear algebra below is exact over Q(sqrt(5)), represented as pairs of
fractions.  No floating-point value enters a verdict.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Sequence


SCHEMA = "oph.finite_born_frame.v1"
ISSUE = 687
VERDICT = (
    "EXACT_FINITE_FRAME_RANK_GAP__"
    "CONTEXT_ADDITIVITY_NOT_DENSITY_FORCING__"
    "PUBLIC_QUANTUM_EFFECT_PRODUCER_MISSING"
)

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "runtime" / "finite_born_frame_certificate.json"
PORT_FRAME_SOURCE = REPOSITORY_ROOT / "Lean" / "Screen" / "PortFrameGram.lean"

# Cross-repository provenance is recorded but is not a runtime dependency.
# The adapter receipt itself says that its Pauli/projector realization is a
# mathematical construction and that physical promotion is forbidden.
EXTERNAL_ADAPTER = {
    "repository": "FloatingPragma/oph-physics-sim",
    "repository_commit_at_projection": "4aa2ce703b5cd172f687c7c9bc3d2d4aa04ed11e",
    "path": "data/quantum/icosahedral_chsh_candidate_receipt.json",
    "file_sha256_at_projection": "9c3ccbe47ec9939647a4f5ba5d18c2f5a37e96b6910702cba174e11ca336c214",
    "internal_receipt_sha256": "f6e8b4e7fbaf5c3539afc17643d125d2691c028dd46cd25564f5af90761bb3b7",
    "schema": "oph.icosahedral_chsh_candidate.v1",
    "verdict": "EXACT_PROJECTIVE_BRANCH_CANDIDATE__TWO_WING_COMPLETED_RECORD_SOURCE_PRODUCER_MISSING",
    "adapter_status": "mathematical construction from the source-derived spin lift",
    "physical_promotion_allowed": False,
}

EXTERNAL_CENTRAL_ATOM_MANIFEST = {
    "repository": "FloatingPragma/oph-physics-sim",
    "repository_commit_at_projection": "4aa2ce703b5cd172f687c7c9bc3d2d4aa04ed11e",
    "path": "tests/fixtures/echosahedral_federation_reference.json",
    "file_sha256_at_projection": "12ec97358ade25f919f9981f1cd7c99c2b27aaaaa9ce14bd2320d73a9c1bfc14",
    "schema": "oph.echosahedral_selector_manifest.v1",
    "declared_port_atom_count": 12,
    "atoms_pairwise_orthogonal": True,
    "atoms_sum_to_one": True,
    "status": "declared source-carrier manifest, not a generated public quantum instrument",
}


class CertificateError(ValueError):
    """Typed fail-closed certificate error."""

    def __init__(self, code: str, message: str) -> None:
        super().__init__(f"{code}: {message}")
        self.code = code


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise CertificateError(code, message)


def _fraction(value: int | Fraction) -> Fraction:
    return value if isinstance(value, Fraction) else Fraction(value)


@dataclass(frozen=True)
class Q5:
    """An exact element ``a + b sqrt(5)`` of Q(sqrt(5))."""

    a: Fraction
    b: Fraction = Fraction(0)

    @staticmethod
    def of(a: int | Fraction = 0, b: int | Fraction = 0) -> "Q5":
        return Q5(_fraction(a), _fraction(b))

    def __add__(self, other: "Q5") -> "Q5":
        return Q5(self.a + other.a, self.b + other.b)

    def __radd__(self, other: int) -> "Q5":
        return self if other == 0 else Q5.of(other) + self

    def __neg__(self) -> "Q5":
        return Q5(-self.a, -self.b)

    def __sub__(self, other: "Q5") -> "Q5":
        return self + (-other)

    def __mul__(self, other: "Q5") -> "Q5":
        return Q5(
            self.a * other.a + 5 * self.b * other.b,
            self.a * other.b + self.b * other.a,
        )

    def __truediv__(self, other: "Q5") -> "Q5":
        return self * other.inverse()

    def inverse(self) -> "Q5":
        norm = self.a * self.a - 5 * self.b * self.b
        require(norm != 0, "Q5_DIVISION", "division by zero")
        return Q5(self.a / norm, -self.b / norm)

    def is_zero(self) -> bool:
        return self.a == 0 and self.b == 0

    def sign(self) -> int:
        """Exact sign, using irrationality of sqrt(5)."""

        if self.b == 0:
            return (self.a > 0) - (self.a < 0)
        if self.a == 0:
            return (self.b > 0) - (self.b < 0)
        if (self.a > 0) == (self.b > 0):
            return 1 if self.a > 0 else -1
        comparison = self.a * self.a - 5 * self.b * self.b
        require(comparison != 0, "Q5_SIGN", "unexpected rational square root of five")
        if self.a > 0:
            return 1 if comparison > 0 else -1
        return 1 if comparison < 0 else -1

    def render(self) -> str:
        def part(value: Fraction) -> str:
            return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"

        if self.b == 0:
            return part(self.a)
        if self.a == 0:
            return f"{part(self.b)}*sqrt(5)"
        return f"{part(self.a)} + {part(self.b)}*sqrt(5)"


ZERO = Q5.of(0)
ONE = Q5.of(1)
HALF = Q5.of(Fraction(1, 2))
PHI = Q5.of(Fraction(1, 2), Fraction(1, 2))
LENGTH_SQUARED = PHI + Q5.of(2)

Vector = tuple[Q5, Q5, Q5]


def vec(x: Q5 | int, y: Q5 | int, z: Q5 | int) -> Vector:
    return tuple(value if isinstance(value, Q5) else Q5.of(value) for value in (x, y, z))  # type: ignore[return-value]


# This is the lexicographically first oriented exact vertex frame after
# relabelling to the PortFrameGram convention, whose antipode is i -> 11-i.
DEFAULT_FRAME: tuple[Vector, ...] = (
    vec(0, 1, PHI),
    vec(0, -1, PHI),
    vec(PHI, 0, 1),
    vec(-PHI, 0, 1),
    vec(1, PHI, 0),
    vec(1, -PHI, 0),
    vec(-1, PHI, 0),
    vec(-1, -PHI, 0),
    vec(PHI, 0, -1),
    vec(-PHI, 0, -1),
    vec(0, 1, -PHI),
    vec(0, -1, -PHI),
)

ANTIPODE = tuple(11 - index for index in range(12))
NEIGHBORS: tuple[tuple[int, ...], ...] = (
    (1, 2, 3, 4, 6),
    (0, 2, 3, 5, 7),
    (0, 1, 4, 5, 8),
    (0, 1, 6, 7, 9),
    (0, 2, 6, 8, 10),
    (1, 2, 7, 8, 11),
    (0, 3, 4, 9, 10),
    (1, 3, 5, 9, 11),
    (2, 4, 5, 10, 11),
    (3, 6, 7, 10, 11),
    (4, 6, 8, 9, 11),
    (5, 7, 8, 9, 10),
)


def dot(left: Sequence[Q5], right: Sequence[Q5]) -> Q5:
    return sum((a * b for a, b in zip(left, right, strict=True)), ZERO)


def matrix_rank(rows: Iterable[Sequence[Q5]]) -> int:
    work = [list(row) for row in rows]
    if not work:
        return 0
    columns = len(work[0])
    require(all(len(row) == columns for row in work), "MATRIX_SHAPE", "ragged matrix")
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, len(work)) if not work[row][column].is_zero()),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column].inverse()
        work[pivot_row] = [entry * scale for entry in work[pivot_row]]
        for row in range(len(work)):
            if row == pivot_row or work[row][column].is_zero():
                continue
            factor = work[row][column]
            work[row] = [
                work[row][index] - factor * work[pivot_row][index]
                for index in range(columns)
            ]
        pivot_row += 1
        if pivot_row == len(work):
            break
    return pivot_row


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_sha256(payload: Any) -> str:
    data = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(data.encode("utf-8")).hexdigest()


def context_matrix() -> list[list[Q5]]:
    rows: list[list[Q5]] = []
    for index in range(6):
        row = [ZERO for _ in range(12)]
        row[index] = ONE
        row[ANTIPODE[index]] = ONE
        rows.append(row)
    return rows


def centered_from_s(frame: Sequence[Vector], s: Vector) -> tuple[Q5, ...]:
    return tuple(dot(axis, s) for axis in frame)


def weights_from_centered(centered: Sequence[Q5]) -> tuple[Q5, ...]:
    return tuple((ONE + value) * HALF for value in centered)


def context_additive(weights: Sequence[Q5]) -> bool:
    return len(weights) == 12 and all(
        (weights[index] + weights[ANTIPODE[index]] - ONE).is_zero()
        for index in range(6)
    )


def box_positive(weights: Sequence[Q5]) -> bool:
    return all(value.sign() >= 0 and (ONE - value).sign() >= 0 for value in weights)


def representative_centered(centered: Sequence[Q5]) -> tuple[Q5, ...]:
    return tuple(centered[index] for index in range(6))


def frame_relations(centered_representatives: Sequence[Q5]) -> tuple[Q5, Q5, Q5]:
    c0, c1, c2, c3, c4, c5 = centered_representatives
    return (
        c0 + c1 - PHI * (c2 + c3),
        c2 - c3 - PHI * (c4 + c5),
        c4 - c5 - PHI * (c0 - c1),
    )


def frame_representable(centered: Sequence[Q5]) -> bool:
    return context_additive(weights_from_centered(centered)) and all(
        relation.is_zero() for relation in frame_relations(representative_centered(centered))
    )


def reconstruct_s(centered: Sequence[Q5]) -> Vector:
    c0, c1, c2, c3, c4, c5 = representative_centered(centered)
    return (
        (c4 + c5) * HALF,
        (c0 - c1) * HALF,
        (c2 + c3) * HALF,
    )


def density_norm_squared(s: Vector) -> Q5:
    return LENGTH_SQUARED * dot(s, s)


def _serialized_vector(values: Sequence[Q5]) -> list[str]:
    return [value.render() for value in values]


def validate_frame(frame: Sequence[Vector], antipode: Sequence[int] = ANTIPODE) -> dict[str, Any]:
    require(len(frame) == 12, "FRAME_CARDINALITY", "expected twelve axes")
    require(len(antipode) == 12, "ANTIPODE_CARDINALITY", "expected twelve antipodes")
    require(
        all(0 <= antipode[index] < 12 and antipode[antipode[index]] == index for index in range(12)),
        "ANTIPODE_INVOLUTION",
        "antipode must be an involution",
    )
    require(
        all(
            all((frame[antipode[index]][coordinate] + frame[index][coordinate]).is_zero() for coordinate in range(3))
            for index in range(12)
        ),
        "ANTIPODE_FRAME",
        "antipodes must be opposite axes",
    )
    require(
        all((dot(axis, axis) - LENGTH_SQUARED).is_zero() for axis in frame),
        "FRAME_NORM",
        "axes must have the common exact norm",
    )
    require(matrix_rank(frame) == 3, "FRAME_RANK", "axis frame must have rank three")

    gram = [[dot(frame[left], frame[right]) for right in range(12)] for left in range(12)]
    require(matrix_rank(gram) == 3, "GRAM_RANK", "Gram matrix must have rank three")
    normalized_adjacent = PHI / LENGTH_SQUARED
    require((normalized_adjacent - Q5.of(0, Fraction(1, 5))).is_zero(), "GRAM_VALUE", "adjacent overlap must be 1/sqrt(5)")
    for left in range(12):
        for right in range(12):
            expected = (
                LENGTH_SQUARED
                if left == right
                else -LENGTH_SQUARED
                if right == antipode[left]
                else PHI
                if right in NEIGHBORS[left]
                else -PHI
            )
            require((gram[left][right] - expected).is_zero(), "GRAM_INCIDENCE", "frame and source incidence disagree")

    transpose_product = [
        [sum((frame[index][row] * frame[index][column] for index in range(12)), ZERO) for column in range(3)]
        for row in range(3)
    ]
    for row in range(3):
        for column in range(3):
            expected = Q5.of(4) * LENGTH_SQUARED if row == column else ZERO
            require((transpose_product[row][column] - expected).is_zero(), "TIGHT_FRAME", "frame operator is not scalar")

    orthogonal_pairs = [
        [left, right]
        for left in range(12)
        for right in range(left + 1, 12)
        if (gram[left][right] + LENGTH_SQUARED).is_zero()
    ]
    require(orthogonal_pairs == [[index, 11 - index] for index in range(6)], "PROJECTIVE_CONTEXTS", "unexpected binary context inventory")
    effect_rows = [[ONE, *axis] for axis in frame]
    require(matrix_rank(effect_rows) == 4, "EFFECT_SPAN", "projectors must span the qubit Hermitian operator system")
    return {
        "axis_rank": 3,
        "gram_rank": 3,
        "effect_operator_system_rank": 4,
        "common_unnormalized_axis_norm_squared": LENGTH_SQUARED.render(),
        "normalized_gram_values": ["1", "1/5*sqrt(5)", "-1/5*sqrt(5)", "-1"],
        "orthogonal_contexts": orthogonal_pairs,
        "context_count": len(orthogonal_pairs),
        "projectors_per_context": 2,
        "contexts_interlock": False,
    }


def _central_atom_branch() -> dict[str, Any]:
    constraint_rank = matrix_rank([[ONE for _ in range(12)]])
    require(constraint_rank == 1, "CENTRAL_CONTEXT_RANK", "central normalization rank drifted")
    rho0_spectrum = [Fraction(1, 12)] * 12
    rho1_spectrum = [Fraction(1, 8), Fraction(1, 24)] + [Fraction(1, 12)] * 10
    require(sum(rho0_spectrum) == sum(rho1_spectrum) == 1, "CENTRAL_DENSITY_TRACE", "density trace drifted")
    require(all(value >= 0 for value in rho0_spectrum + rho1_spectrum), "CENTRAL_DENSITY_POSITIVITY", "density witness is not positive")
    return {
        "object_status": "declared_source_manifest_central_port_atoms",
        "atom_count": 12,
        "context_count": 1,
        "normalization_constraint_rank": constraint_rank,
        "normalized_weight_affine_dimension": 11,
        "nonnegative_weight_region": "the eleven-simplex",
        "state_on_commutative_algebra": {
            "representation": "rho_commutative = diag(w_0,...,w_11)",
            "unique": True,
        },
        "ambient_M12_density_representation": {
            "unique": False,
            "uniform_weight_counterexample": {
                "common_atom_weights": ["1/12"] * 12,
                "rho_0": "I_12/12",
                "rho_1": "I_12/12 + (E_01+E_10)/24",
                "rho_0_spectrum": [str(value) for value in rho0_spectrum],
                "rho_1_spectrum": [str(value) for value in rho1_spectrum],
                "same_declared_atom_weights": True,
                "distinct_density_matrices": True,
            },
        },
        "boundary": "This is a classical central-record context. It is not the noncommuting qubit projector adapter.",
    }


def _projective_branch(frame: Sequence[Vector]) -> dict[str, Any]:
    frame_audit = validate_frame(frame)
    context_rank = matrix_rank(context_matrix())
    require(context_rank == 6, "CONTEXT_RANK", "six binary contexts must have rank six")

    # Exact admissible but non-Born weight: one binary context is deterministic,
    # all other contexts are fair.
    nonborn_centered = (ONE, ZERO, ZERO, ZERO, ZERO, ZERO, -ZERO, -ZERO, -ZERO, -ZERO, -ZERO, -ONE)
    nonborn_weights = weights_from_centered(nonborn_centered)
    nonborn_residuals = frame_relations(representative_centered(nonborn_centered))
    require(context_additive(nonborn_weights), "NONBORN_CONTROL", "control lost context additivity")
    require(box_positive(nonborn_weights), "NONBORN_CONTROL", "control left the probability cube")
    require(any(not value.is_zero() for value in nonborn_residuals), "NONBORN_CONTROL", "control accidentally became Born representable")

    # Exact represented and context-positive, but non-density, witness.  The
    # associated Bloch norm is (phi+2)*9/25 > 1 while every declared port
    # weight remains in [0,1].
    outside_s = vec(Q5.of(Fraction(3, 5)), 0, 0)
    outside_centered = centered_from_s(frame, outside_s)
    outside_weights = weights_from_centered(outside_centered)
    outside_norm = density_norm_squared(outside_s)
    require(context_additive(outside_weights), "NONPOSITIVE_CONTROL", "represented control lost additivity")
    require(box_positive(outside_weights), "NONPOSITIVE_CONTROL", "represented control left the finite probability cube")
    require(frame_representable(outside_centered), "NONPOSITIVE_CONTROL", "represented control left the Born affine slice")
    require((outside_norm - ONE).sign() > 0, "NONPOSITIVE_CONTROL", "represented control is unexpectedly a density matrix")

    inside_s = vec(Q5.of(Fraction(1, 5)), 0, 0)
    inside_centered = centered_from_s(frame, inside_s)
    inside_weights = weights_from_centered(inside_centered)
    inside_norm = density_norm_squared(inside_s)
    require(context_additive(inside_weights) and box_positive(inside_weights), "DENSITY_CONTROL", "valid density control failed finite admissibility")
    require(frame_representable(inside_centered), "DENSITY_CONTROL", "valid density control failed frame relations")
    require(reconstruct_s(inside_centered) == inside_s, "TOMOGRAPHY", "exact reconstruction failed")
    require((ONE - inside_norm).sign() >= 0, "DENSITY_CONTROL", "valid density control is nonpositive")

    relation_rows = [
        [ONE, ONE, -PHI, -PHI, ZERO, ZERO],
        [ZERO, ZERO, ONE, -ONE, -PHI, -PHI],
        [-PHI, PHI, ZERO, ZERO, ONE, -ONE],
    ]
    require(matrix_rank(relation_rows) == 3, "BORN_RELATION_RANK", "Born relation rank drifted")

    return {
        "object_status": "declared_qubit_projector_adapter_from_source_derived_spin_geometry",
        "physical_public_effect_attachment": False,
        "registered_effect_inventory": "zero, identity, and the twelve rank-one projector labels only; no rational effect closure or full positive cone is source-registered",
        "field": "Q(sqrt(5)) with exact rational coefficient pairs",
        "frame": frame_audit,
        "normalization_and_additivity": {
            "constraint_matrix_rank": context_rank,
            "weight_coordinate_count": 12,
            "affine_dimension": 6,
            "nonnegative_region": "the six-cube [0,1]^6 via one probability per antipodal context",
            "parameterization": "w_i=a_i and w_(11-i)=1-a_i for i=0,...,5",
            "noncontextuality_content": "label consistency only; every rank-one projector occurs in exactly one binary context",
        },
        "born_hermitian_slice": {
            "centered_coordinates": "c_i=2w_i-1",
            "representative_ports": [0, 1, 2, 3, 4, 5],
            "exact_relations": [
                "c0+c1=phi*(c2+c3)",
                "c2-c3=phi*(c4+c5)",
                "c4-c5=phi*(c0-c1)",
            ],
            "relation_rank": 3,
            "affine_dimension_inside_normalized_weights": 3,
            "codimension_inside_context_additive_weights": 3,
            "tomographic_coordinates": {
                "x": "(c4+c5)/2",
                "y": "(c0-c1)/2",
                "z": "(c2+c3)/2",
            },
            "weight_formula": "w_i=(1+v_i dot s)/2",
            "unique_trace_one_hermitian_representation_when_it_exists": True,
            "density_positivity_criterion": "(phi+2)*(x^2+y^2+z^2)<=1",
        },
        "decision": {
            "every_context_admissible_weight_has_a_hermitian_trace_representation": False,
            "every_context_admissible_weight_has_a_density_representation": False,
            "density_representation_is_unique_when_it_exists": True,
            "reason": "the context-additive affine space has dimension six, while the informationally complete Born/Hermitian slice has dimension three; finite positivity on the twelve listed projectors does not imply positivity on the full qubit effect cone",
        },
        "exact_controls": {
            "context_additive_nonborn_weight": {
                "weights": _serialized_vector(nonborn_weights),
                "frame_relation_residuals": _serialized_vector(nonborn_residuals),
                "all_weights_in_unit_interval": True,
                "has_trace_one_hermitian_representation": False,
            },
            "born_affine_but_nonpositive_weight": {
                "s": _serialized_vector(outside_s),
                "weights": _serialized_vector(outside_weights),
                "density_norm_squared": outside_norm.render(),
                "all_weights_in_unit_interval": True,
                "has_unique_trace_one_hermitian_representation": True,
                "has_density_representation": False,
            },
            "valid_density_weight": {
                "s": _serialized_vector(inside_s),
                "weights": _serialized_vector(inside_weights),
                "density_norm_squared": inside_norm.render(),
                "exact_tomography_recovers_s": True,
                "has_unique_density_representation": True,
            },
        },
    }


def build_certificate(
    *,
    frame: Sequence[Vector] = DEFAULT_FRAME,
    antipode: Sequence[int] = ANTIPODE,
) -> dict[str, Any]:
    # Validate caller-supplied antipode before using the canonical context
    # equations.  This is an adversarial hook used by the test suite.
    validate_frame(frame, antipode)
    require(tuple(antipode) == ANTIPODE, "ANTIPODE_LINEAGE", "B11 uses the PortFrameGram i -> 11-i labeling")
    require(PORT_FRAME_SOURCE.exists(), "SOURCE_LINEAGE", "PortFrameGram source is missing")
    body: dict[str, Any] = {
        "schema": SCHEMA,
        "issue": ISSUE,
        "verdict": VERDICT,
        "source_lineage": {
            "local_exact_port_frame": {
                "path": "Lean/Screen/PortFrameGram.lean",
                "sha256": _sha256(PORT_FRAME_SOURCE),
                "status": "declared exact icosahedral Gram/incidence packet",
            },
            "external_central_atom_manifest": EXTERNAL_CENTRAL_ATOM_MANIFEST,
            "external_projective_adapter": EXTERNAL_ADAPTER,
            "lineage_decision": "No registered producer emits these qubit projectors as completed public effects. The exact calculation is conditional on the declared spinor adapter; the central port atoms remain a separate classical family.",
        },
        "declared_central_atom_branch": _central_atom_branch(),
        "declared_spinor_projective_branch": _projective_branch(frame),
        "closure_assessment": {
            "finite_enumeration_complete": True,
            "declared_projector_adapter_enumeration_complete": True,
            "literal_source_produced_public_effect_set_available": False,
            "general_gleason_or_busch_theorem_attempted": False,
            "positive_born_derivation_obtained": False,
            "bounded_negative_exit": True,
            "first_missing_positive_antecedent": "a source-produced interlocking public effect family rich enough that context additivity and positivity constrain the full positive cone",
            "issue_687_closure_ready": True,
            "closure_mode": "scoped exact no-go plus typed physical-attachment boundary",
        },
        "claim_boundary": (
            "The declared twelve-axis qubit adapter is informationally complete, so a Born/Hermitian representation is unique whenever it exists. "
            "But its six disjoint binary contexts leave a six-dimensional cube of normalized additive weights, whereas the Born/Hermitian slice is only three-dimensional. "
            "Some additive unit-interval weights have no Hermitian trace representation, and some represented unit-interval weights correspond to a nonpositive matrix. "
            "Therefore this finite family does not derive the Born rule. The source-manifest central atoms instead form one classical context and do not determine an ambient M12 density matrix. "
            "Neither branch supplies a source-produced physical quantum measurement instrument."
        ),
    }
    body["receipt_sha256"] = canonical_sha256(body)
    return body


def validate_certificate(receipt: dict[str, Any]) -> None:
    require(receipt.get("schema") == SCHEMA, "RECEIPT_SCHEMA", "unexpected schema")
    require(receipt.get("issue") == ISSUE, "RECEIPT_ISSUE", "unexpected issue")
    require(receipt.get("verdict") == VERDICT, "RECEIPT_VERDICT", "unexpected verdict")
    body = {key: value for key, value in receipt.items() if key != "receipt_sha256"}
    require(receipt.get("receipt_sha256") == canonical_sha256(body), "RECEIPT_HASH", "receipt hash mismatch")
    require(receipt == build_certificate(), "RECEIPT_REPLAY", "receipt differs from exact replay")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true", help="check the frozen receipt instead of rewriting it")
    args = parser.parse_args()
    receipt = build_certificate()
    if args.check:
        frozen = json.loads(args.out.read_text(encoding="utf-8"))
        validate_certificate(frozen)
        print("FINITE_BORN_FRAME_CERTIFICATE_VALID")
        return 0
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
