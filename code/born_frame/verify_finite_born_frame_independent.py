#!/usr/bin/env python3
"""Independent exact verifier for the issue-687 finite Born-frame receipt.

This module intentionally does not import the producer.  It reconstructs the
frame ranks, context constraints, tomography relations, and both adversarial
weights with a second Q(sqrt(5)) implementation.
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
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RECEIPT = Path(__file__).resolve().parent / "runtime" / "finite_born_frame_certificate.json"


class VerificationError(ValueError):
    pass


def need(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


@dataclass(frozen=True)
class K:
    p: Fraction
    q: Fraction = Fraction(0)

    @staticmethod
    def make(p: int | Fraction = 0, q: int | Fraction = 0) -> "K":
        return K(Fraction(p), Fraction(q))

    def __add__(self, other: "K") -> "K":
        return K(self.p + other.p, self.q + other.q)

    def __radd__(self, other: int) -> "K":
        return self if other == 0 else K.make(other) + self

    def __neg__(self) -> "K":
        return K(-self.p, -self.q)

    def __sub__(self, other: "K") -> "K":
        return self + (-other)

    def __mul__(self, other: "K") -> "K":
        return K(self.p * other.p + 5 * self.q * other.q, self.p * other.q + self.q * other.p)

    def inv(self) -> "K":
        denominator = self.p * self.p - 5 * self.q * self.q
        need(denominator != 0, "independent division by zero")
        return K(self.p / denominator, -self.q / denominator)

    def __truediv__(self, other: "K") -> "K":
        return self * other.inv()

    def zero(self) -> bool:
        return self.p == 0 and self.q == 0

    def sign(self) -> int:
        if self.q == 0:
            return (self.p > 0) - (self.p < 0)
        if self.p == 0:
            return (self.q > 0) - (self.q < 0)
        if (self.p > 0) == (self.q > 0):
            return 1 if self.p > 0 else -1
        comparison = self.p * self.p - 5 * self.q * self.q
        need(comparison != 0, "independent Q(sqrt5) sign degeneracy")
        if self.p > 0:
            return 1 if comparison > 0 else -1
        return 1 if comparison < 0 else -1


Z = K.make()
ONE_K = K.make(1)
H = K.make(Fraction(1, 2))
PHI = K.make(Fraction(1, 2), Fraction(1, 2))
L2 = PHI + K.make(2)


def parse_fraction(text: str) -> Fraction:
    return Fraction(text)


def parse_k(text: str) -> K:
    if " + " in text:
        rational, radical = text.split(" + ", 1)
        need(radical.endswith("*sqrt(5)"), f"malformed radical {text}")
        return K(parse_fraction(rational), parse_fraction(radical[: -len("*sqrt(5)")]))
    if text.endswith("*sqrt(5)"):
        return K.make(0, parse_fraction(text[: -len("*sqrt(5)")]))
    return K.make(parse_fraction(text))


def rank(rows: Iterable[Sequence[K]]) -> int:
    matrix = [list(row) for row in rows]
    if not matrix:
        return 0
    width = len(matrix[0])
    need(all(len(row) == width for row in matrix), "independent ragged matrix")
    active = 0
    for column in range(width):
        pivot = next((row for row in range(active, len(matrix)) if not matrix[row][column].zero()), None)
        if pivot is None:
            continue
        matrix[active], matrix[pivot] = matrix[pivot], matrix[active]
        inverse = matrix[active][column].inv()
        matrix[active] = [entry * inverse for entry in matrix[active]]
        for row in range(len(matrix)):
            if row == active or matrix[row][column].zero():
                continue
            factor = matrix[row][column]
            matrix[row] = [matrix[row][index] - factor * matrix[active][index] for index in range(width)]
        active += 1
    return active


def vector(x: K | int, y: K | int, z: K | int) -> tuple[K, K, K]:
    return tuple(value if isinstance(value, K) else K.make(value) for value in (x, y, z))  # type: ignore[return-value]


FRAME = (
    vector(0, 1, PHI),
    vector(0, -1, PHI),
    vector(PHI, 0, 1),
    vector(-PHI, 0, 1),
    vector(1, PHI, 0),
    vector(1, -PHI, 0),
    vector(-1, PHI, 0),
    vector(-1, -PHI, 0),
    vector(PHI, 0, -1),
    vector(-PHI, 0, -1),
    vector(0, 1, -PHI),
    vector(0, -1, -PHI),
)
ANTI = tuple(11 - index for index in range(12))


def dot(left: Sequence[K], right: Sequence[K]) -> K:
    return sum((x * y for x, y in zip(left, right, strict=True)), Z)


def canonical_hash(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
    return hashlib.sha256(encoded).hexdigest()


def additive(weights: Sequence[K]) -> bool:
    return len(weights) == 12 and all((weights[index] + weights[ANTI[index]] - ONE_K).zero() for index in range(6))


def in_box(weights: Sequence[K]) -> bool:
    return all(value.sign() >= 0 and (ONE_K - value).sign() >= 0 for value in weights)


def relations(centered: Sequence[K]) -> tuple[K, K, K]:
    a, b, c, d, e, f = centered[:6]
    return (
        a + b - PHI * (c + d),
        c - d - PHI * (e + f),
        e - f - PHI * (a - b),
    )


def weights_from_s(s: Sequence[K]) -> list[K]:
    centered = [dot(axis, s) for axis in FRAME]
    return [(ONE_K + value) * H for value in centered]


def verify(receipt: dict[str, Any]) -> dict[str, Any]:
    need(receipt.get("schema") == SCHEMA, "schema mismatch")
    need(receipt.get("issue") == ISSUE, "issue mismatch")
    need(receipt.get("verdict") == VERDICT, "verdict mismatch")
    body = {key: value for key, value in receipt.items() if key != "receipt_sha256"}
    need(receipt.get("receipt_sha256") == canonical_hash(body), "receipt hash mismatch")

    source = receipt["source_lineage"]
    local = source["local_exact_port_frame"]
    local_path = ROOT / local["path"]
    need(local_path.exists(), "local source path missing")
    need(hashlib.sha256(local_path.read_bytes()).hexdigest() == local["sha256"], "local source hash mismatch")
    adapter = source["external_projective_adapter"]
    need(adapter["adapter_status"] == "mathematical construction from the source-derived spin lift", "adapter status promoted")
    need(adapter["physical_promotion_allowed"] is False, "adapter physical promotion")
    central_manifest = source["external_central_atom_manifest"]
    need(central_manifest["declared_port_atom_count"] == 12, "central source atom count mismatch")
    need(central_manifest["atoms_pairwise_orthogonal"] is True, "central source orthogonality missing")
    need(central_manifest["atoms_sum_to_one"] is True, "central source normalization missing")
    need("declared source-carrier manifest" in central_manifest["status"], "central source status promoted")
    need("No registered producer" in source["lineage_decision"], "source boundary missing")

    central = receipt["declared_central_atom_branch"]
    need(central["normalization_constraint_rank"] == 1, "central constraint rank mismatch")
    need(central["normalized_weight_affine_dimension"] == 11, "central affine dimension mismatch")
    need(central["state_on_commutative_algebra"]["unique"] is True, "commutative state uniqueness mismatch")
    ambient = central["ambient_M12_density_representation"]
    need(ambient["unique"] is False, "ambient M12 uniqueness overstated")
    witness = ambient["uniform_weight_counterexample"]
    spectrum0 = [Fraction(value) for value in witness["rho_0_spectrum"]]
    spectrum1 = [Fraction(value) for value in witness["rho_1_spectrum"]]
    need(sum(spectrum0) == sum(spectrum1) == 1, "central density trace failure")
    need(all(value >= 0 for value in spectrum0 + spectrum1), "central density positivity failure")
    need(witness["same_declared_atom_weights"] is True and witness["distinct_density_matrices"] is True, "central ambiguity witness failure")

    gram = [[dot(left, right) for right in FRAME] for left in FRAME]
    need(rank(FRAME) == rank(gram) == 3, "independent frame rank failure")
    need(all((dot(FRAME[index], FRAME[ANTI[index]]) + L2).zero() for index in range(12)), "independent antipode failure")
    need(rank([[ONE_K, *axis] for axis in FRAME]) == 4, "independent effect span failure")
    context_rows = []
    for index in range(6):
        row = [Z] * 12
        row[index] = ONE_K
        row[ANTI[index]] = ONE_K
        context_rows.append(row)
    need(rank(context_rows) == 6, "independent context rank failure")
    relation_rows = (
        (ONE_K, ONE_K, -PHI, -PHI, Z, Z),
        (Z, Z, ONE_K, -ONE_K, -PHI, -PHI),
        (-PHI, PHI, Z, Z, ONE_K, -ONE_K),
    )
    need(rank(relation_rows) == 3, "independent relation rank failure")

    branch = receipt["declared_spinor_projective_branch"]
    need(branch["physical_public_effect_attachment"] is False, "false public effect attachment")
    need("no rational effect closure" in branch["registered_effect_inventory"], "effect inventory overstated")
    frame_report = branch["frame"]
    need(frame_report["axis_rank"] == 3 and frame_report["effect_operator_system_rank"] == 4, "serialized frame rank failure")
    additivity = branch["normalization_and_additivity"]
    need(additivity["constraint_matrix_rank"] == 6 and additivity["affine_dimension"] == 6, "serialized context dimension failure")
    born = branch["born_hermitian_slice"]
    need(born["relation_rank"] == 3 and born["affine_dimension_inside_normalized_weights"] == 3, "serialized Born dimension failure")
    need(born["unique_trace_one_hermitian_representation_when_it_exists"] is True, "tomography uniqueness missing")
    decision = branch["decision"]
    need(decision["every_context_admissible_weight_has_a_hermitian_trace_representation"] is False, "false Hermitian universality")
    need(decision["every_context_admissible_weight_has_a_density_representation"] is False, "false density universality")
    need(decision["density_representation_is_unique_when_it_exists"] is True, "conditional density uniqueness missing")

    controls = branch["exact_controls"]
    nonborn = controls["context_additive_nonborn_weight"]
    nonborn_weights = [parse_k(value) for value in nonborn["weights"]]
    nonborn_centered = [value * K.make(2) - ONE_K for value in nonborn_weights]
    need(additive(nonborn_weights) and in_box(nonborn_weights), "non-Born control not admissible")
    need(any(not value.zero() for value in relations(nonborn_centered)), "non-Born control became representable")
    need(nonborn["has_trace_one_hermitian_representation"] is False, "non-Born control flag promoted")

    outside = controls["born_affine_but_nonpositive_weight"]
    outside_s = tuple(parse_k(value) for value in outside["s"])
    outside_weights = [parse_k(value) for value in outside["weights"]]
    outside_centered = [value * K.make(2) - ONE_K for value in outside_weights]
    outside_norm = L2 * dot(outside_s, outside_s)
    need(outside_weights == weights_from_s(outside_s), "nonpositive control weights do not match frame")
    need(additive(outside_weights) and in_box(outside_weights), "nonpositive control not finite-admissible")
    need(all(value.zero() for value in relations(outside_centered)), "nonpositive control left Born slice")
    need((outside_norm - ONE_K).sign() > 0, "nonpositive control became a density")
    need(parse_k(outside["density_norm_squared"]) == outside_norm, "nonpositive control norm mismatch")
    need(outside["has_density_representation"] is False, "nonpositive control flag promoted")

    inside = controls["valid_density_weight"]
    inside_s = tuple(parse_k(value) for value in inside["s"])
    inside_weights = [parse_k(value) for value in inside["weights"]]
    inside_norm = L2 * dot(inside_s, inside_s)
    need(inside_weights == weights_from_s(inside_s), "density control weights do not match frame")
    need(additive(inside_weights) and in_box(inside_weights), "density control not finite-admissible")
    need((ONE_K - inside_norm).sign() >= 0, "density control became nonpositive")
    need(parse_k(inside["density_norm_squared"]) == inside_norm, "density control norm mismatch")
    need(inside["exact_tomography_recovers_s"] is True and inside["has_unique_density_representation"] is True, "density control uniqueness mismatch")

    closure = receipt["closure_assessment"]
    need(closure["finite_enumeration_complete"] is True, "finite enumeration not closed")
    need(closure["declared_projector_adapter_enumeration_complete"] is True, "adapter enumeration not closed")
    need(closure["literal_source_produced_public_effect_set_available"] is False, "false source effect producer")
    need(closure["positive_born_derivation_obtained"] is False, "false positive Born derivation")
    need(closure["bounded_negative_exit"] is True, "bounded negative exit missing")
    need(closure["issue_687_closure_ready"] is True, "closure assessment drifted")
    need("does not derive the Born rule" in receipt["claim_boundary"], "claim boundary overstates result")
    return {
        "status": "PASS",
        "context_affine_dimension": 6,
        "born_affine_dimension": 3,
        "central_affine_dimension": 11,
        "conditional_density_representation_unique": True,
        "all_context_weights_density_representable": False,
        "physical_public_effect_attachment": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    args = parser.parse_args()
    receipt = json.loads(args.receipt.read_text(encoding="utf-8"))
    print(json.dumps(verify(receipt), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
