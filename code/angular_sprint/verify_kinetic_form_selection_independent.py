#!/usr/bin/env python3
"""Independent exact verifier for the issue-646 kinetic-form receipt.

This verifier deliberately imports neither ``kinetic_form_selection_certificate``
nor any helper under ``code/a5_closure``.  It reads the pinned serialized
carrier and current receipt, rebuilds the four exact adjacency projectors over
Q(sqrt(5)), rebuilds the Hilbert--Schmidt form from its serialized band data,
parses the complete serialized structure-constant table, and checks all 1728
ad-invariance equations.  It then reconstructs the two simple ideals, their
Killing forms, and the Killing-relative coefficients.

The one-loop family/Higgs cancellation is recomputed separately with exact
rational affine arithmetic.  No coupling measurement or comparison payload is
read.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Mapping

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
DEFAULT_RECEIPT = HERE / "runtime" / "kinetic_form_selection_receipt.json"
PORT_RECEIPT_REL = Path(
    "code/a5_closure/receipts/port_current_inner_reference.receipt.json"
)
CARRIER_MANIFEST_REL = Path(
    "code/a5_closure/manifests/echosahedral_federation_reference.json"
)
MATTER_ATTACHMENT_REL = Path(
    "code/a5_closure/manifests/matter_attachment_receipt.json"
)
PORT_RECEIPT = REPO_ROOT / PORT_RECEIPT_REL
CARRIER_MANIFEST = REPO_ROOT / CARRIER_MANIFEST_REL
MATTER_ATTACHMENT = REPO_ROOT / MATTER_ATTACHMENT_REL

SCHEMA = "oph.kinetic_form_selection_independent_verification.v1"
EXPECTED_COMMAND = (
    "python3 code/angular_sprint/"
    "verify_kinetic_form_selection_independent.py --receipt "
    "code/angular_sprint/runtime/kinetic_form_selection_receipt.json"
)
EXPECTED_CLASSIFICATION = (
    "independent exact reconstruction from the pinned serialized carrier "
    "incidence, structure constants, band form, and rank-fifteen matter "
    "table; the verifier imports no producer or a5-closure algebra helper"
)


class IndependentVerificationError(ValueError):
    """The independent verifier rejected malformed exact input."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise IndependentVerificationError(message)


def strict_json_bytes(payload: bytes, *, label: str) -> dict[str, Any]:
    """Decode one JSON object while rejecting duplicate mapping keys."""

    def object_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        out: dict[str, Any] = {}
        for key, value in pairs:
            if key in out:
                raise IndependentVerificationError(
                    f"{label}: duplicate JSON key {key!r}"
                )
            out[key] = value
        return out

    def reject_constant(value: str) -> Any:
        raise IndependentVerificationError(
            f"{label}: non-finite JSON constant {value!r}"
        )

    try:
        value = json.loads(
            payload.decode("ascii"),
            object_pairs_hook=object_pairs,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise IndependentVerificationError(f"{label}: invalid JSON: {exc}") from exc
    require(isinstance(value, dict), f"{label}: top-level value is not an object")
    return value


def strict_load(path: Path) -> dict[str, Any]:
    return strict_json_bytes(path.read_bytes(), label=path.as_posix())


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("ascii")


def tagged_sha256(payload: bytes) -> str:
    return "sha256:" + hashlib.sha256(payload).hexdigest()


def semantic_json_sha256(value: Any) -> str:
    """Match the upstream untagged canonical-JSON binding."""

    payload = json.dumps(
        value,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


@dataclass(frozen=True)
class Q5:
    """Exact element ``a + b sqrt(5)``."""

    a: Fraction
    b: Fraction = Fraction(0)

    @staticmethod
    def of(a: int | str | Fraction, b: int | str | Fraction = 0) -> "Q5":
        return Q5(Fraction(a), Fraction(b))

    def __add__(self, other: "Q5") -> "Q5":
        return Q5(self.a + other.a, self.b + other.b)

    def __sub__(self, other: "Q5") -> "Q5":
        return Q5(self.a - other.a, self.b - other.b)

    def __neg__(self) -> "Q5":
        return Q5(-self.a, -self.b)

    def __mul__(self, other: "Q5") -> "Q5":
        return Q5(
            self.a * other.a + 5 * self.b * other.b,
            self.a * other.b + self.b * other.a,
        )

    def __truediv__(self, other: "Q5") -> "Q5":
        norm = other.a * other.a - 5 * other.b * other.b
        require(norm != 0, "division by zero in Q(sqrt(5))")
        return Q5(
            (self.a * other.a - 5 * self.b * other.b) / norm,
            (self.b * other.a - self.a * other.b) / norm,
        )

    def is_zero(self) -> bool:
        return self.a == 0 and self.b == 0


ZERO = Q5.of(0)
ONE = Q5.of(1)


def parse_q5(text: Any) -> Q5:
    require(isinstance(text, str), "Q(sqrt(5)) value is not a string")
    cleaned = text.replace(" ", "").replace("sqrt5", "sqrt(5)")
    if "*sqrt(5)" not in cleaned:
        return Q5.of(Fraction(cleaned))
    prefix = cleaned.removesuffix("*sqrt(5)").replace("+-", "-")
    split = None
    for index in range(1, len(prefix)):
        if prefix[index] in "+-":
            split = index
    if split is None:
        return Q5.of(0, Fraction(prefix))
    return Q5.of(Fraction(prefix[:split]), Fraction(prefix[split:]))


Matrix = list[list[Q5]]
Vector = list[Q5]


def identity(size: int) -> Matrix:
    return [[ONE if i == j else ZERO for j in range(size)] for i in range(size)]


def matrix_add(left: Matrix, right: Matrix) -> Matrix:
    return [
        [left[i][j] + right[i][j] for j in range(len(left[0]))]
        for i in range(len(left))
    ]


def matrix_scale(factor: Q5, matrix: Matrix) -> Matrix:
    return [[factor * value for value in row] for row in matrix]


def matrix_mul(left: Matrix, right: Matrix) -> Matrix:
    out = [[ZERO for _ in range(len(right[0]))] for _ in range(len(left))]
    for i, row in enumerate(left):
        for k, value in enumerate(row):
            if value.is_zero():
                continue
            for j in range(len(right[0])):
                out[i][j] = out[i][j] + value * right[k][j]
    return out


def matrix_trace(matrix: Matrix) -> Q5:
    total = ZERO
    for index in range(len(matrix)):
        total = total + matrix[index][index]
    return total


def rank(rows: Iterable[Iterable[Q5]]) -> int:
    matrix = [list(row) for row in rows]
    if not matrix:
        return 0
    width = len(matrix[0])
    pivot_row = 0
    for column in range(width):
        pivot = next(
            (row for row in range(pivot_row, len(matrix)) if not matrix[row][column].is_zero()),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        pivot_value = matrix[pivot_row][column]
        matrix[pivot_row] = [value / pivot_value for value in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row == pivot_row or matrix[row][column].is_zero():
                continue
            factor = matrix[row][column]
            matrix[row] = [
                value - factor * base
                for value, base in zip(matrix[row], matrix[pivot_row], strict=True)
            ]
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return pivot_row


def independent_columns(matrix: Matrix, expected_rank: int) -> list[Vector]:
    columns = [
        [matrix[row][column] for row in range(len(matrix))]
        for column in range(len(matrix[0]))
    ]
    basis: list[Vector] = []
    for column in columns:
        if rank(basis + [column]) > len(basis):
            basis.append(column)
    require(len(basis) == expected_rank, f"projector rank != {expected_rank}")
    return basis


def carrier_adjacency(carrier: Mapping[str, Any]) -> Matrix:
    block = carrier.get("carrier")
    require(isinstance(block, Mapping), "carrier block missing")
    ports = block.get("ports")
    edges = block.get("edges")
    require(isinstance(ports, list) and len(ports) == 12, "carrier is not twelve-port")
    require(isinstance(edges, list) and len(edges) == 30, "carrier is not thirty-edge")
    index = {name: i for i, name in enumerate(ports)}
    require(len(index) == 12, "carrier ports are not distinct")
    adjacency = [[ZERO for _ in range(12)] for _ in range(12)]
    for edge in edges:
        require(isinstance(edge, list) and len(edge) == 2, "malformed carrier edge")
        left, right = (index.get(edge[0]), index.get(edge[1]))
        require(left is not None and right is not None and left != right, "invalid carrier edge")
        require(adjacency[left][right].is_zero(), "duplicate carrier edge")
        adjacency[left][right] = adjacency[right][left] = ONE
    require(all(sum(not value.is_zero() for value in row) == 5 for row in adjacency), "carrier is not five-regular")
    return adjacency


def spectral_projectors(adjacency: Matrix) -> dict[str, Matrix]:
    eigenvalues = {
        "unit_band": Q5.of(5),
        "frame_band": Q5.of(0, 1),
        "kernel_band": Q5.of(0, -1),
        "quintet_band": Q5.of(-1),
    }
    unit = identity(12)
    projectors: dict[str, Matrix] = {}
    for name, eigenvalue in eigenvalues.items():
        projector = identity(12)
        denominator = ONE
        for other in eigenvalues.values():
            if other == eigenvalue:
                continue
            shifted = matrix_add(adjacency, matrix_scale(-other, unit))
            projector = matrix_mul(projector, shifted)
            denominator = denominator * (eigenvalue - other)
        projectors[name] = matrix_scale(ONE / denominator, projector)
    total = [[ZERO for _ in range(12)] for _ in range(12)]
    for projector in projectors.values():
        total = matrix_add(total, projector)
    require(total == unit, "spectral projectors are incomplete")
    for left_name, left in projectors.items():
        require(matrix_mul(left, left) == left, f"{left_name} is not idempotent")
        for right_name, right in projectors.items():
            if left_name < right_name:
                require(
                    all(value.is_zero() for row in matrix_mul(left, right) for value in row),
                    f"{left_name}/{right_name} projectors are not orthogonal",
                )
    return projectors


def parse_structure_constants(parent: Mapping[str, Any]) -> list[list[Vector]]:
    closure = parent.get("closure")
    require(isinstance(closure, Mapping), "parent closure block missing")
    serialized = closure.get("structure_constants")
    require(isinstance(serialized, Mapping) and len(serialized) == 66, "structure-constant table is incomplete")
    table = [[[ZERO for _ in range(12)] for _ in range(12)] for _ in range(12)]
    seen: set[tuple[int, int]] = set()
    for key, row in serialized.items():
        require(isinstance(key, str) and key.startswith("[") and key.endswith("]"), "malformed bracket key")
        parts = key[1:-1].split(",")
        require(len(parts) == 2, "malformed bracket index")
        left, right = (int(parts[0]), int(parts[1]))
        require(0 <= left < right < 12 and (left, right) not in seen, "invalid bracket index")
        require(isinstance(row, list) and len(row) == 12, "malformed bracket row")
        parsed = [parse_q5(value) for value in row]
        table[left][right] = parsed
        table[right][left] = [-value for value in parsed]
        seen.add((left, right))
    require(len(seen) == 66, "not every unordered bracket pair is serialized")
    return table


def bracket(table: list[list[Vector]], left: Vector, right: Vector) -> Vector:
    out = [ZERO for _ in range(12)]
    for i, x in enumerate(left):
        if x.is_zero():
            continue
        for j, y in enumerate(right):
            if y.is_zero():
                continue
            factor = x * y
            for k, coefficient in enumerate(table[i][j]):
                out[k] = out[k] + factor * coefficient
    return out


def check_jacobi(table: list[list[Vector]]) -> int:
    standard = [[ONE if i == j else ZERO for i in range(12)] for j in range(12)]
    checks = 0
    for i in range(12):
        for j in range(12):
            for k in range(12):
                total = [ZERO for _ in range(12)]
                for term in (
                    bracket(table, standard[i], bracket(table, standard[j], standard[k])),
                    bracket(table, standard[j], bracket(table, standard[k], standard[i])),
                    bracket(table, standard[k], bracket(table, standard[i], standard[j])),
                ):
                    total = [a + b for a, b in zip(total, term, strict=True)]
                require(all(value.is_zero() for value in total), f"Jacobi failure at {(i, j, k)}")
                checks += 1
    return checks


def form_from_band_data(parent: Mapping[str, Any], projectors: Mapping[str, Matrix]) -> Matrix:
    compactness = parent.get("compactness")
    require(isinstance(compactness, Mapping), "parent compactness block missing")
    coefficients = compactness.get("hilbert_schmidt_pullback_band_coefficients")
    require(isinstance(coefficients, Mapping), "parent band form missing")
    gram = [[ZERO for _ in range(12)] for _ in range(12)]
    for name in ("unit_band", "frame_band", "kernel_band", "quintet_band"):
        require(name in coefficients and name in projectors, f"missing band {name}")
        gram = matrix_add(gram, matrix_scale(parse_q5(coefficients[name]), projectors[name]))
    return gram


def check_ad_invariance(gram: Matrix, table: list[list[Vector]]) -> int:
    checks = 0
    for i in range(12):
        for j in range(12):
            for k in range(12):
                value = ZERO
                for m in range(12):
                    value = value + table[i][j][m] * gram[m][k]
                    value = value + table[i][k][m] * gram[j][m]
                require(value.is_zero(), f"ad-invariance failure at {(i, j, k)}")
                checks += 1
    return checks


def form_on(gram: Matrix, vectors: list[Vector]) -> Matrix:
    out = [[ZERO for _ in vectors] for _ in vectors]
    for a, left in enumerate(vectors):
        for b, right in enumerate(vectors):
            for i, x in enumerate(left):
                for j, y in enumerate(right):
                    out[a][b] = out[a][b] + x * y * gram[i][j]
    return out


def killing_on(table: list[list[Vector]], vectors: list[Vector]) -> Matrix:
    def ad_matrix(vector: Vector) -> Matrix:
        out = [[ZERO for _ in range(12)] for _ in range(12)]
        for i, coefficient in enumerate(vector):
            if coefficient.is_zero():
                continue
            for j in range(12):
                for k in range(12):
                    out[k][j] = out[k][j] + coefficient * table[i][j][k]
        return out

    ads = [ad_matrix(vector) for vector in vectors]
    return [
        [matrix_trace(matrix_mul(left, right)) for right in ads]
        for left in ads
    ]


def relative_coefficient(form: Matrix, killing: Matrix, label: str) -> Q5:
    pivot = next(
        (i for i in range(len(killing)) if not killing[i][i].is_zero()),
        None,
    )
    require(pivot is not None, f"{label}: degenerate Killing form")
    coefficient = form[pivot][pivot] / (-killing[pivot][pivot])
    for i in range(len(form)):
        for j in range(len(form)):
            require(
                form[i][j] == coefficient * (-killing[i][j]),
                f"{label}: form is not Killing-proportional",
            )
    return coefficient


def recompute_exact_core(parent: Mapping[str, Any], carrier: Mapping[str, Any]) -> dict[str, Any]:
    adjacency = carrier_adjacency(carrier)
    projectors = spectral_projectors(adjacency)
    expected_ranks = {"unit_band": 1, "frame_band": 3, "kernel_band": 3, "quintet_band": 5}
    bases = {
        name: independent_columns(projector, expected_ranks[name])
        for name, projector in projectors.items()
    }
    table = parse_structure_constants(parent)
    jacobi_checks = check_jacobi(table)
    gram = form_from_band_data(parent, projectors)
    ad_checks = check_ad_invariance(gram, table)
    su2 = bases["kernel_band"]
    su3 = bases["frame_band"] + bases["quintet_band"]
    require(rank(su2) == 3 and rank(su3) == 8, "simple-ideal dimensions drift")
    c2 = relative_coefficient(form_on(gram, su2), killing_on(table, su2), "su2")
    c3 = relative_coefficient(form_on(gram, su3), killing_on(table, su3), "su3")
    require(c2 == Q5.of(1), "independent su2 coefficient != 1")
    require(c3 == Q5.of(Fraction(1, 6)), "independent su3 coefficient != 1/6")
    return {
        "jacobi_checks": jacobi_checks,
        "ad_invariance_checks": ad_checks,
        "su2": str(c2.a) if c2.b == 0 else None,
        "su3": str(c3.a) if c3.b == 0 else None,
        "ratio": str((c2 / c3).a) if (c2 / c3).b == 0 else None,
    }


def recompute_registered_matter_indices(
    matter: Mapping[str, Any],
) -> tuple[Fraction, Fraction, Fraction]:
    require(
        matter.get("schema") == "oph.local-domain-matter-attachment.v1",
        "matter attachment schema mismatch",
    )
    certificate = matter.get("generation_certificate")
    require(isinstance(certificate, Mapping), "matter generation certificate missing")
    rows = certificate.get("rows")
    require(isinstance(rows, list) and len(rows) == 5, "matter generation table drift")
    u1 = Fraction(0)
    su2 = Fraction(0)
    su3 = Fraction(0)
    state_count = 0
    for row in rows:
        require(isinstance(row, Mapping), "malformed matter row")
        color = int(row.get("color_dimension"))
        weak = int(row.get("weak_dimension"))
        hypercharge = Fraction(row.get("hypercharge"))
        states = int(row.get("weyl_states"))
        require(states == color * weak, "matter state multiplicity mismatch")
        state_count += states
        u1 += Fraction(states) * hypercharge * hypercharge
        if weak == 2:
            su2 += Fraction(color, 2)
        else:
            require(weak == 1, "unsupported weak matter representation")
        if color == 3:
            su3 += Fraction(weak, 2)
        else:
            require(color == 1, "unsupported color matter representation")
    require(
        state_count == certificate.get("weyl_state_count") == 15,
        "matter state count mismatch",
    )
    require(
        (u1, su2, su3) == (Fraction(10, 3), Fraction(2), Fraction(2)),
        "matter representation indices mismatch",
    )
    return u1, su2, su3


def recompute_general_cancellation(
    k: tuple[Fraction, Fraction, Fraction],
) -> dict[str, Any]:
    # Affine triples encode constant + nG coefficient + nH coefficient.
    beta = (
        (Fraction(0), Fraction(20, 9), Fraction(1, 6)),
        (Fraction(-22, 3), Fraction(4, 3), Fraction(1, 6)),
        (Fraction(-11), Fraction(4, 3), Fraction(0)),
    )

    def scaled(factor: Fraction, row: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
        return tuple(factor * value for value in row)

    def minus(left: tuple[Fraction, ...], right: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
        return tuple(a - b for a, b in zip(left, right, strict=True))

    cofactors = (
        minus(scaled(k[1], beta[2]), scaled(k[2], beta[1])),
        minus(scaled(k[2], beta[0]), scaled(k[0], beta[2])),
        minus(scaled(k[0], beta[1]), scaled(k[1], beta[0])),
    )
    expected = (
        (Fraction(-22, 3), Fraction(0), Fraction(-1, 3)),
        (Fraction(110, 3), Fraction(0), Fraction(1, 3)),
        (Fraction(-220, 9), Fraction(0), Fraction(2, 9)),
    )
    require(cofactors == expected, "general determinant cofactor mismatch")
    specialization = tuple(row[0] + 3 * row[1] + row[2] for row in cofactors)
    require(
        specialization == (Fraction(-23, 3), Fraction(37), Fraction(-218, 9)),
        "(nG,nH)=(3,1) specialization mismatch",
    )
    return {
        "cofactors": [[str(value) for value in row] for row in cofactors],
        "specialization": [str(value) for value in specialization],
        "nG_cancels": all(row[1] == 0 for row in cofactors),
    }


def verify_receipt(
    receipt: Mapping[str, Any],
    *,
    repo_root: Path = REPO_ROOT,
) -> dict[str, Any]:
    """Verify one decoded receipt and return a fail-closed result."""

    reasons: list[str] = []
    core: dict[str, Any] | None = None
    cancellation: dict[str, Any] | None = None
    matter_index_strings: list[str] | None = None
    try:
        require(
            receipt.get("schema") == "oph.kinetic_form_selection_receipt.v1",
            "receipt schema mismatch",
        )
        self_payload = dict(receipt)
        stated_hash = self_payload.pop("receipt_sha256", None)
        require(
            stated_hash == tagged_sha256(canonical_json_bytes(self_payload)),
            "receipt self-hash mismatch",
        )

        pins = receipt.get("parent_pins")
        require(isinstance(pins, list) and len(pins) == 3, "parent pin set mismatch")
        expected_pin_paths = {
            "code/a5_closure/receipts/port_current_inner_reference.receipt.json",
            "Lean/Screen/RGRepresentationFrontier.lean",
            "code/a5_closure/manifests/matter_attachment_receipt.json",
        }
        require(
            {pin.get("path") for pin in pins if isinstance(pin, Mapping)}
            == expected_pin_paths,
            "parent pin paths mismatch",
        )
        for pin in pins:
            require(isinstance(pin, Mapping), "malformed parent pin")
            relative = pin.get("path")
            require(isinstance(relative, str) and "\\" not in relative, "invalid pin path")
            parsed = Path(relative)
            require(
                not parsed.is_absolute()
                and "." not in parsed.parts
                and ".." not in parsed.parts,
                "noncanonical parent pin path",
            )
            path = repo_root / relative
            payload = path.read_bytes()
            require(pin.get("bytes") == len(payload), f"parent byte count drift: {relative}")
            require(pin.get("sha256") == tagged_sha256(payload), f"parent hash drift: {relative}")

        parent = strict_load(repo_root / PORT_RECEIPT_REL)
        carrier = strict_load(repo_root / CARRIER_MANIFEST_REL)
        matter = strict_load(repo_root / MATTER_ATTACHMENT_REL)
        require(
            parent.get("carrier_manifest_sha256")
            == semantic_json_sha256(carrier),
            "carrier manifest is not bound by the pinned current receipt",
        )
        core = recompute_exact_core(parent, carrier)
        require(
            receipt.get("ad_invariance", {}).get("verified_basis_triples")
            == core["ad_invariance_checks"]
            == 1728,
            "ad-invariance count mismatch",
        )
        row = receipt.get("killing_relative_coefficients", {})
        require(
            row.get("su2") == core["su2"]
            and row.get("su3") == core["su3"]
            and row.get("ratio_su2_over_su3") == core["ratio"],
            "Killing-relative coefficient mismatch",
        )

        matter_indices = recompute_registered_matter_indices(matter)
        matter_index_strings = [str(value) for value in matter_indices]
        matter_branch = receipt.get("matter_trace_branch", {})
        require(
            matter_branch.get("per_copy_weyl_indices")
            == {
                "u1": str(matter_indices[0]),
                "su2": str(matter_indices[1]),
                "su3": str(matter_indices[2]),
            },
            "matter-index receipt mismatch",
        )
        cancellation = recompute_general_cancellation(matter_indices)
        frozen = matter_branch.get("frozen_rg_statistic", {})
        general = frozen.get("general_family_higgs_cancellation", {})
        require(
            general.get("determinant_cofactors_constant_nG_nH")
            == cancellation["cofactors"],
            "general cofactor certificate mismatch",
        )
        require(
            general.get("nG_coefficients_cancel_exactly") is True
            and cancellation["nG_cancels"],
            "generation cancellation flag mismatch",
        )
        require(
            frozen.get("exact_cofactors") == cancellation["specialization"],
            "declared (nG,nH) specialization cofactor mismatch",
        )

        comparison = receipt.get("comparison_boundary", {})
        require(
            comparison.get("public_measurement_read") is False
            and comparison.get("comparison_permitted") is False,
            "comparison boundary is open",
        )
        physical = receipt.get("physical_selection_boundary", {})
        for key in (
            "current_lift_source_selected",
            "kinetic_form_source_selected",
            "physical_sector_selected",
            "physical_continuum_gauge_action_identified",
        ):
            require(physical.get(key) is False, f"forbidden physical promotion: {key}")
        verification = receipt.get("verification", {})
        require(
            verification.get("command") == EXPECTED_COMMAND
            and verification.get("classification") == EXPECTED_CLASSIFICATION
            and verification.get(
                "reconstructs_all_1728_ad_invariance_equations"
            )
            is True
            and verification.get("reconstructs_killing_relative_coefficients")
            is True
            and verification.get("reconstructs_registered_matter_indices")
            is True
            and verification.get("reconstructs_general_nG_nH_cancellation")
            is True
            and verification.get("comparison_data_read") is False,
            "verification scope or comparison-data boundary drift",
        )
    except (IndependentVerificationError, KeyError, OSError, TypeError, ValueError) as exc:
        reasons.append(str(exc))

    return {
        "schema": SCHEMA,
        "receipt": not reasons,
        "reasons": reasons,
        "exact_core": core,
        "registered_matter_indices_u1_su2_su3": matter_index_strings,
        "general_family_higgs_cancellation": cancellation,
        "verification_scope": (
            "independent exact reconstruction from pinned serialized source "
            "objects; no producer or a5-closure algebra helper imported"
        ),
        "comparison_data_read": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    args = parser.parse_args()
    try:
        receipt = strict_load(args.receipt)
        result = verify_receipt(receipt)
    except (IndependentVerificationError, OSError) as exc:
        result = {
            "schema": SCHEMA,
            "receipt": False,
            "reasons": [str(exc)],
            "comparison_data_read": False,
        }
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["receipt"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
