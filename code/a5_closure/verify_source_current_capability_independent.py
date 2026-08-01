#!/usr/bin/env python3
"""Independent verifier for the bounded issue-566 source capability packet."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


MODULE_DIR = Path(__file__).resolve().parent
DEFAULT_PROJECTION = MODULE_DIR / "manifests" / "source_current_capability_projection.json"
DEFAULT_RECEIPT = MODULE_DIR / "receipts" / "source_current_capability.receipt.json"
PROJECTION_SCHEMA = "oph.source_current_capability_projection.v1"
RECEIPT_SCHEMA = "oph.source_current_capability_receipt.v1"
VERDICT = "BOUNDED_REGISTERED_RESPONSE_ALGEBRA_INSUFFICIENT__ORDER_SENSITIVE_SOURCE_BRIDGE_OPEN"


class VerificationError(ValueError):
    pass


def check(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def load(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    check(isinstance(value, dict), f"{path} is not an object")
    return value


def digest(value: Any) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return "sha256:" + hashlib.sha256(payload).hexdigest()


def file_digest(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


Matrix = list[list[Fraction]]


def eye(n: int) -> Matrix:
    return [[Fraction(i == j) for j in range(n)] for i in range(n)]


def multiply(x: Matrix, y: Matrix) -> Matrix:
    n = len(x)
    return [
        [sum(x[i][k] * y[k][j] for k in range(n)) for j in range(n)]
        for i in range(n)
    ]


def vector(matrix: Matrix) -> list[Fraction]:
    return sum(matrix, [])


def rank(vectors: Sequence[Sequence[Fraction]]) -> int:
    work = [list(map(Fraction, row)) for row in vectors]
    if not work:
        return 0
    pivot_row = 0
    for column in range(len(work[0])):
        pivot = next(
            (row for row in range(pivot_row, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column]
        work[pivot_row] = [entry / scale for entry in work[pivot_row]]
        for row in range(len(work)):
            if row == pivot_row or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [
                work[row][j] - scale * work[pivot_row][j]
                for j in range(len(work[row]))
            ]
        pivot_row += 1
        if pivot_row == len(work):
            break
    return pivot_row


def determinant(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    work = [list(map(Fraction, row)) for row in matrix]
    dimension = len(work)
    check(all(len(row) == dimension for row in work), "square determinant input")
    value = Fraction(1)
    for column in range(dimension):
        pivot = next(
            (row for row in range(column, dimension) if work[row][column]),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            value = -value
        pivot_value = work[column][column]
        value *= pivot_value
        for row in range(column + 1, dimension):
            factor = work[row][column] / pivot_value
            work[row] = [
                work[row][entry] - factor * work[column][entry]
                for entry in range(dimension)
            ]
    return value


def matrix_linear_combination(terms: Sequence[tuple[Fraction, Matrix]]) -> Matrix:
    dimension = len(terms[0][1])
    return [
        [
            sum(coefficient * matrix[row][column] for coefficient, matrix in terms)
            for column in range(dimension)
        ]
        for row in range(dimension)
    ]


def trace_product(left: Matrix, right: Matrix) -> Fraction:
    return sum(
        left[row][column] * right[column][row]
        for row in range(len(left))
        for column in range(len(left))
    )


def render(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def permutation_matrix(permutation: Sequence[int]) -> Matrix:
    n = len(permutation)
    result = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    for source, target in enumerate(permutation):
        result[int(target)][source] = Fraction(1)
    return result


def cyclic(face: Sequence[int]) -> tuple[int, int, int]:
    a, b, c = map(int, face)
    return min((a, b, c), (b, c, a), (c, a, b))


def graph_automorphisms(adjacency: Matrix) -> list[tuple[int, ...]]:
    n = len(adjacency)
    image: list[int | None] = [None] * n
    occupied: set[int] = set()
    answer: list[tuple[int, ...]] = []

    def visit(source: int) -> None:
        if source == n:
            answer.append(tuple(int(value) for value in image))  # type: ignore[arg-type]
            return
        for target in range(n):
            if target in occupied:
                continue
            if any(
                adjacency[source][earlier] != adjacency[target][image[earlier]]  # type: ignore[index]
                for earlier in range(source)
            ):
                continue
            image[source] = target
            occupied.add(target)
            visit(source + 1)
            occupied.remove(target)
            image[source] = None

    visit(0)
    return answer


def verify(projection_path: Path, receipt_path: Path) -> dict[str, Any]:
    projection = load(projection_path)
    receipt = load(receipt_path)
    check(projection.get("schema") == PROJECTION_SCHEMA, "projection schema")
    projection_body = {
        key: value for key, value in projection.items() if key != "projection_sha256"
    }
    check(projection.get("projection_sha256") == digest(projection_body), "projection hash")
    check(receipt.get("schema") == RECEIPT_SCHEMA, "receipt schema")
    receipt_body = {
        key: value for key, value in receipt.items() if key != "receipt_sha256"
    }
    check(receipt.get("receipt_sha256") == digest(receipt_body), "receipt hash")
    check(receipt.get("verdict") == VERDICT, "receipt verdict")
    check(receipt.get("physical_current_source_bridge_attained") is False, "false source gate")
    check(receipt.get("abstract_forced_lie_type_theorem_preserved") is True, "theorem boundary")
    check(
        receipt.get("audit_scope")
        == (
            "the algebra generated by the single registered adjacency recurrence "
            "and its polynomial maximal-distance readback; no bound is asserted for "
            "an unregistered order-sensitive response history"
        ),
        "bounded audit scope",
    )
    check(receipt.get("projection_sha256") == projection["projection_sha256"], "receipt projection pin")

    sources = projection.get("source_files")
    check(isinstance(sources, Mapping), "source pins")
    carrier_path = MODULE_DIR / str(sources.get("carrier"))
    response_path = MODULE_DIR / str(sources.get("response"))
    check(carrier_path.is_file() and response_path.is_file(), "source path")
    check(sources.get("carrier_file_sha256") == file_digest(carrier_path), "carrier file pin")
    check(sources.get("response_file_sha256") == file_digest(response_path), "response file pin")
    carrier_source = load(carrier_path)
    response_source = load(response_path)
    check(sources.get("carrier_payload_sha256") == digest(carrier_source), "carrier payload pin")
    response_body = {
        key: value for key, value in response_source.items() if key != "artifact_sha256"
    }
    check(response_source.get("artifact_sha256") == digest(response_body), "response payload hash")
    check(sources.get("response_payload_sha256") == response_source["artifact_sha256"], "response payload pin")

    carrier = projection.get("carrier")
    check(isinstance(carrier, Mapping), "carrier projection")
    ports = carrier.get("port_order")
    check(isinstance(ports, list) and len(ports) == 12 and len(set(ports)) == 12, "port basis")
    source_ports = carrier_source.get("carrier", {}).get("ports")
    check(ports == source_ports, "projected port order")
    source_index = {port: position for position, port in enumerate(source_ports)}
    source_edges = sorted(
        sorted((source_index[left], source_index[right]))
        for left, right in carrier_source.get("carrier", {}).get("edges", [])
    )
    source_faces = [
        [source_index[port] for port in face]
        for face in carrier_source.get("carrier", {}).get("oriented_faces", [])
    ]
    check(carrier.get("edges") == source_edges, "projected source edges")
    check(carrier.get("oriented_faces") == source_faces, "projected source faces")
    response_binding = response_source.get("carrier_binding", {})
    check(response_binding.get("port_order") == ports, "response port binding")
    check(
        response_binding.get("carrier_manifest_sha256")
        == digest(carrier_source).removeprefix("sha256:"),
        "response carrier binding",
    )
    adjacency = [[Fraction(0) for _ in range(12)] for _ in range(12)]
    for edge in carrier.get("edges", []):
        check(isinstance(edge, list) and len(edge) == 2, "edge row")
        left, right = map(int, edge)
        check(0 <= left < right < 12 and not adjacency[left][right], "edge range")
        adjacency[left][right] = adjacency[right][left] = Fraction(1)
    check(sum(sum(row) for row in adjacency) == 60, "edge count")
    check(all(sum(row) == 5 for row in adjacency), "degree sequence")

    powers = [eye(12)]
    for _ in range(6):
        powers.append(multiply(powers[-1], adjacency))
    traces = projection.get("registered_recurrence", {}).get("steps")
    check(isinstance(traces, list) and len(traces) == 4, "trace count")
    for step, row in enumerate(traces):
        check(isinstance(row, Mapping) and row.get("step") == step, "trace order")
        observed = [[Fraction(entry) for entry in values] for values in row.get("matrix", [])]
        check(observed == powers[step], f"trace {step} replay")
    recurrence = projection["registered_recurrence"]
    check(recurrence.get("raw_runtime_history_serialized") is False, "trace scope")
    source_protocol = response_source.get("source_response", {}).get(
        "impulse_readback_protocol", {}
    )
    check(recurrence.get("filter_coefficients") == source_protocol.get("homogeneous_filter_coefficients"), "filter source")
    check(source_protocol.get("target_labels_used") is False, "target firewall")
    check(source_protocol.get("downstream_labels_used") is False, "downstream firewall")
    check(
        projection.get("registered_refinement_maps")
        == response_source.get("physical_refinement_maps", {}).get(
            "port_persistence_maps"
        ),
        "refinement source",
    )
    runtime = response_source.get("provenance", {}).get("runtime_binding", {})
    unitary_channel = projection.get("registered_unitary_channel", {})
    check(
        unitary_channel
        == {
            "coupling": runtime.get("coupling"),
            "evolution": runtime.get("evolution"),
            "generator": "L = 5*I - A",
            "reversible_response_source": runtime.get(
                "reversible_response_source"
            ),
        },
        "unitary channel source",
    )
    check(
        runtime.get("coupling") == "combinatorial_icosahedral_graph_laplacian"
        and runtime.get("evolution") == "finite_unitary_exp_minus_i_u_gL"
        and runtime.get("reversible_response_source")
        == "finite_unitary_carrier_channel",
        "unitary channel typing",
    )

    response_map = response_source.get("source_response", {}).get("antipode_port_map")
    check(projection.get("readback", {}).get("antipode_port_map") == response_map, "antipode source")
    coefficients = (Fraction(-1), Fraction(1, 2), Fraction(2, 5), Fraction(-1, 10))
    response = [
        [sum(coefficients[k] * powers[k][i][j] for k in range(4)) for j in range(12)]
        for i in range(12)
    ]
    expected = [[Fraction(0) for _ in range(12)] for _ in range(12)]
    for source, target in enumerate(response_map):
        expected[source][int(target)] = Fraction(-1)
    check(response == expected, "response polynomial")
    check(multiply(response, response) == eye(12), "response involution")

    response_rank = rank([vector(matrix) for matrix in powers[:4]])
    closure_rank = rank([vector(matrix) for matrix in powers])
    check(response_rank == closure_rank == 4, "response algebra dimension")
    laplacian = matrix_linear_combination(
        ((Fraction(5), powers[0]), (Fraction(-1), powers[1]))
    )
    check(
        rank([vector(matrix) for matrix in powers[:4]] + [vector(laplacian)])
        == response_rank,
        "unitary generator span",
    )
    polynomial_residual = matrix_linear_combination(
        (
            (Fraction(1), powers[4]),
            (Fraction(-4), powers[3]),
            (Fraction(-10), powers[2]),
            (Fraction(20), powers[1]),
            (Fraction(25), powers[0]),
        )
    )
    check(
        polynomial_residual
        == [[Fraction(0) for _ in range(12)] for _ in range(12)],
        "minimal polynomial",
    )
    check(
        all(multiply(x, y) == multiply(y, x) for x in powers[:4] for y in powers[:4]),
        "response algebra commutator",
    )
    gram = [
        [trace_product(left, right) for right in powers[:4]]
        for left in powers[:4]
    ]
    symmetric_basis = all(
        matrix
        == [list(row) for row in zip(*matrix, strict=True)]
        for matrix in powers[:4]
    )
    check(symmetric_basis, "symmetric response basis")
    principal_minors = [
        determinant([row[:size] for row in gram[:size]])
        for size in range(1, 5)
    ]
    check(
        gram
        == [
            [Fraction(12), Fraction(0), Fraction(60), Fraction(120)],
            [Fraction(0), Fraction(60), Fraction(120), Fraction(780)],
            [Fraction(60), Fraction(120), Fraction(780), Fraction(3120)],
            [Fraction(120), Fraction(780), Fraction(3120), Fraction(16380)],
        ],
        "trace-pairing Gram",
    )
    check(
        principal_minors
        == [Fraction(12), Fraction(720), Fraction(172800), Fraction(207360000)]
        and all(value > 0 for value in principal_minors),
        "positive skew pairing",
    )

    faces = carrier.get("oriented_faces")
    check(isinstance(faces, list) and len(faces) == 20, "face packet")
    oriented = {cyclic(face) for face in faces}
    automorphisms = graph_automorphisms(adjacency)
    rotations = [
        permutation
        for permutation in automorphisms
        if {
            cyclic([permutation[a], permutation[b], permutation[c]])
            for a, b, c in faces
        }
        == oriented
    ]
    check(len(automorphisms) == 120 and len(rotations) == 60, "recharting count")
    basis = [vector(matrix) for matrix in powers[:4]]
    intersection = sum(
        rank(basis + [vector(permutation_matrix(permutation))]) == 4
        for permutation in rotations
    )
    check(intersection == 1, "response/recharting intersection")

    refinement_rows = projection.get("registered_refinement_maps")
    check(isinstance(refinement_rows, list) and refinement_rows, "refinement rows")
    refinement_natural = True
    for row in refinement_rows:
        check(isinstance(row, Mapping), "refinement row")
        row_body = {key: value for key, value in row.items() if key != "map_hash"}
        check(row.get("map_hash") == digest(row_body), "refinement row hash")
        port_map = row.get("port_map")
        check(
            isinstance(port_map, list) and sorted(port_map) == list(range(12)),
            "refinement port map",
        )
        matrix = permutation_matrix(port_map)
        refinement_natural = (
            refinement_natural
            and multiply(matrix, adjacency) == multiply(adjacency, matrix)
        )
    check(refinement_natural, "refinement recurrence naturality")

    audit = receipt.get("audit")
    check(isinstance(audit, Mapping), "audit block")
    algebra = audit.get("response_word_algebra", {})
    check(algebra.get("exact_dimension") == response_rank, "reported algebra dimension")
    check(algebra.get("closure_rank_through_degree_six") == closure_rank, "reported closure rank")
    check(
        algebra.get("minimal_polynomial")
        == "A^4 - 4*A^3 - 10*A^2 + 20*A + 25*I = 0",
        "reported minimal polynomial",
    )
    check(algebra.get("commutator_nonzero_count") == 0, "reported commutator")
    check(algebra.get("commutative") is True, "reported commutativity")
    check(algebra.get("maximum_independent_skew_generator_count") == 4, "reported tangent bound")
    check(
        algebra.get("skew_adjoint_real_form")
        == "i times span_Q{I,A,A^2,A^3}"
        and algebra.get("skew_adjoint_basis_verified") is True,
        "reported skew-adjoint basis",
    )
    check(
        algebra.get("skew_pairing_gram")
        == [[render(value) for value in row] for row in gram],
        "reported trace-pairing Gram",
    )
    check(
        algebra.get("skew_pairing_leading_principal_minors")
        == [render(value) for value in principal_minors]
        and algebra.get("skew_pairing_positive_definite") is True,
        "reported positive skew pairing",
    )
    check(
        audit.get("registered_unitary_channel_audit")
        == {
            "generator": "L = 5*I - A",
            "generator_in_response_word_algebra": True,
            "single_generator_functional_calculus_dimension_upper_bound": 4,
            "order_sensitive_port_indexed_tangent_available": False,
        },
        "reported unitary channel audit",
    )
    recharting = audit.get("recharting_audit", {})
    check(recharting.get("proper_recharting_count") == 60, "reported rotations")
    check(recharting.get("proper_rechartings_in_response_word_algebra") == 1, "reported intersection")
    refinement = audit.get("refinement_audit", {})
    check(
        refinement
        == {
            "registered_map_count": len(refinement_rows),
            "recurrence_natural_on_registered_maps": True,
            "generator_bracket_implementer_intertwining_available": False,
        },
        "reported refinement audit",
    )
    classifications = audit.get("acceptance_classifications", {})
    expected_classifications = {
        "ordered_two_sided_port_response_histories": (
            "PARTIAL_SOURCE_NATIVE_COMMUTATIVE_RECURRENCE_ONLY"
        ),
        "twelve_independent_skew_adjoint_generators": (
            "IMPOSSIBLE_ON_REGISTERED_RESPONSE_WORD_ALGEBRA"
        ),
        "exact_nonabelian_commutator_reconstruction": (
            "IMPOSSIBLE_ON_REGISTERED_RESPONSE_WORD_ALGEBRA"
        ),
        "closed_overlap_words_cover_sixty_rechartings": (
            "MISSING_STATIC_RECHARTINGS_ONLY"
        ),
        "same_word_projective_implementers": "MISSING",
        "response_identity_component_factorizations": (
            "IMPOSSIBLE_ON_REGISTERED_RESPONSE_WORD_ALGEBRA"
        ),
        "refinement_intertwining": "PARTIAL_SOURCE_NATIVE_RECURRENCE_ONLY",
        "adversarial_controls": "PARTIAL_EXACT_PACKET_CONTROLS",
        "target_firewall": "SOURCE_NATIVE_ATTAINED",
    }
    check(
        classifications == expected_classifications,
        "acceptance classifications",
    )
    return {
        "verdict": VERDICT,
        "response_algebra_dimension": response_rank,
        "proper_recharting_count": len(rotations),
        "response_word_recharting_intersection": intersection,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--projection", type=Path, default=DEFAULT_PROJECTION)
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    args = parser.parse_args()
    report = verify(args.projection, args.receipt)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
