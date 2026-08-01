#!/usr/bin/env python3
"""Exact bounded capability audit for the issue-566 source-current bridge.

The registered source response is the four-step adjacency recurrence together
with its maximal-distance readback.  This certificate reconstructs the full
matrix packet and computes the unital response-word algebra over the rationals.
It does not consume the conditional matrix-current fixture.

The result is a bounded obstruction, not a no-go theorem for the A1/A2 current.
The registered recurrence generates a four-dimensional commutative algebra, so
it cannot by itself emit twelve independent generators, a nonzero bracket, or
the sixty proper recharting implementers.  An order-sensitive source producer
can add data outside this algebra and discharge the open bridge.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


MODULE_DIR = Path(__file__).resolve().parent
CARRIER_PATH = MODULE_DIR / "manifests" / "echosahedral_federation_reference.json"
RESPONSE_PATH = MODULE_DIR / "manifests" / "charged_response_semantic_artifact.json"
PROJECTION_PATH = MODULE_DIR / "manifests" / "source_current_capability_projection.json"
RECEIPT_PATH = MODULE_DIR / "receipts" / "source_current_capability.receipt.json"

PROJECTION_SCHEMA = "oph.source_current_capability_projection.v1"
RECEIPT_SCHEMA = "oph.source_current_capability_receipt.v1"
VERDICT = "BOUNDED_REGISTERED_RESPONSE_ALGEBRA_INSUFFICIENT__ORDER_SENSITIVE_SOURCE_BRIDGE_OPEN"
ISSUE = 566
FORBIDDEN_SOURCE_TOKENS = (
    "standard_model",
    "standard model",
    "su(3)",
    "su3",
    "electroweak",
    "hypercharge",
    "particle_mass",
    "measured_coupling",
    "gauge_target",
)


class CapabilityError(ValueError):
    """Typed fail-closed certificate error."""

    def __init__(self, code: str, message: str) -> None:
        super().__init__(f"{code}: {message}")
        self.code = code


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise CapabilityError(code, message)


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), "JSON_OBJECT", f"{path} is not a JSON object")
    return value


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


def canonical_sha256(value: Any) -> str:
    return "sha256:" + hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_sha256(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


Matrix = list[list[Fraction]]


def identity(n: int) -> Matrix:
    return [
        [Fraction(1 if i == j else 0) for j in range(n)]
        for i in range(n)
    ]


def zero(n: int) -> Matrix:
    return [[Fraction(0) for _ in range(n)] for _ in range(n)]


def matmul(left: Matrix, right: Matrix) -> Matrix:
    n = len(left)
    out = zero(n)
    for i in range(n):
        for k in range(n):
            coefficient = left[i][k]
            if coefficient == 0:
                continue
            for j in range(n):
                out[i][j] += coefficient * right[k][j]
    return out


def add(*terms: tuple[Fraction, Matrix]) -> Matrix:
    n = len(terms[0][1])
    out = zero(n)
    for coefficient, matrix in terms:
        for i in range(n):
            for j in range(n):
                out[i][j] += coefficient * matrix[i][j]
    return out


def transpose(matrix: Matrix) -> Matrix:
    return [list(row) for row in zip(*matrix, strict=True)]


def flatten(matrix: Matrix) -> list[Fraction]:
    return [entry for row in matrix for entry in row]


def rational_rank(rows: Sequence[Sequence[Fraction | int]]) -> int:
    matrix = [[Fraction(entry) for entry in row] for row in rows]
    if not matrix:
        return 0
    columns = len(matrix[0])
    active = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(active, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[active], matrix[pivot] = matrix[pivot], matrix[active]
        scale = matrix[active][column]
        matrix[active] = [entry / scale for entry in matrix[active]]
        for row in range(len(matrix)):
            if row == active or matrix[row][column] == 0:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                matrix[row][j] - factor * matrix[active][j]
                for j in range(columns)
            ]
        active += 1
        if active == len(matrix):
            break
    return active


def determinant(matrix: Sequence[Sequence[Fraction | int]]) -> Fraction:
    work = [[Fraction(entry) for entry in row] for row in matrix]
    n = len(work)
    require(all(len(row) == n for row in work), "DETERMINANT", "matrix is not square")
    result = Fraction(1)
    for column in range(n):
        pivot = next((row for row in range(column, n) if work[row][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            result = -result
        pivot_value = work[column][column]
        result *= pivot_value
        for row in range(column + 1, n):
            factor = work[row][column] / pivot_value
            for j in range(column, n):
                work[row][j] -= factor * work[column][j]
    return result


def trace_product(left: Matrix, right: Matrix) -> Fraction:
    return sum(
        left[i][j] * right[j][i]
        for i in range(len(left))
        for j in range(len(left))
    )


def render_fraction(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def render_matrix(matrix: Matrix) -> list[list[int | str]]:
    return [
        [entry.numerator if entry.denominator == 1 else render_fraction(entry) for entry in row]
        for row in matrix
    ]


def parse_matrix(value: Any, n: int) -> Matrix:
    require(isinstance(value, list) and len(value) == n, "TRACE_MATRIX", "wrong row count")
    matrix: Matrix = []
    for row in value:
        require(isinstance(row, list) and len(row) == n, "TRACE_MATRIX", "wrong column count")
        matrix.append([Fraction(entry) for entry in row])
    return matrix


def source_carrier(carrier_manifest: Mapping[str, Any]) -> dict[str, Any]:
    require(
        carrier_manifest.get("schema") == "oph.echosahedral_selector_manifest.v1",
        "CARRIER_SCHEMA",
        "unexpected carrier schema",
    )
    carrier = carrier_manifest.get("carrier")
    require(isinstance(carrier, Mapping), "CARRIER", "carrier block missing")
    ports = list(carrier.get("ports", []))
    require(len(ports) == 12 and len(set(ports)) == 12, "CARRIER_PORTS", "twelve ports required")
    index = {port: position for position, port in enumerate(ports)}
    adjacency = zero(12)
    indexed_edges: list[list[int]] = []
    for edge in carrier.get("edges", []):
        require(
            isinstance(edge, list)
            and len(edge) == 2
            and edge[0] in index
            and edge[1] in index,
            "CARRIER_EDGE",
            "malformed edge",
        )
        left, right = sorted((index[edge[0]], index[edge[1]]))
        require(left != right and adjacency[left][right] == 0, "CARRIER_EDGE", "edge is not simple")
        adjacency[left][right] = adjacency[right][left] = Fraction(1)
        indexed_edges.append([left, right])
    require(len(indexed_edges) == 30, "CARRIER_EDGE", "thirty edges required")
    require(all(sum(row) == 5 for row in adjacency), "CARRIER_DEGREE", "carrier is not five-regular")
    indexed_faces: list[list[int]] = []
    for face in carrier.get("oriented_faces", []):
        require(
            isinstance(face, list) and len(face) == 3 and all(port in index for port in face),
            "CARRIER_FACE",
            "malformed oriented face",
        )
        row = [index[port] for port in face]
        require(
            all(adjacency[row[k]][row[(k + 1) % 3]] == 1 for k in range(3)),
            "CARRIER_FACE",
            "face is not an incidence triangle",
        )
        indexed_faces.append(row)
    require(len(indexed_faces) == 20, "CARRIER_FACE", "twenty oriented faces required")
    return {
        "ports": ports,
        "adjacency": adjacency,
        "edges": sorted(indexed_edges),
        "oriented_faces": indexed_faces,
    }


def cyclic_face_set(faces: Iterable[Sequence[int]]) -> set[tuple[int, int, int]]:
    result: set[tuple[int, int, int]] = set()
    for face in faces:
        a, b, c = (int(value) for value in face)
        result.add(min((a, b, c), (b, c, a), (c, a, b)))
    return result


def incidence_automorphisms(adjacency: Matrix) -> list[tuple[int, ...]]:
    n = len(adjacency)
    assignment: list[int | None] = [None] * n
    used = [False] * n
    result: list[tuple[int, ...]] = []

    def consistent(node: int, image: int) -> bool:
        return all(
            adjacency[node][other] == adjacency[image][assignment[other]]  # type: ignore[index]
            for other in range(node)
        )

    def search(node: int) -> None:
        if node == n:
            result.append(tuple(int(value) for value in assignment))  # type: ignore[arg-type]
            return
        for image in range(n):
            if used[image] or not consistent(node, image):
                continue
            assignment[node] = image
            used[image] = True
            search(node + 1)
            used[image] = False
            assignment[node] = None

    search(0)
    return result


def permutation_matrix(permutation: Sequence[int]) -> Matrix:
    n = len(permutation)
    matrix = zero(n)
    for source, target in enumerate(permutation):
        matrix[int(target)][source] = Fraction(1)
    return matrix


def proper_rechartings(adjacency: Matrix, faces: Sequence[Sequence[int]]) -> list[tuple[int, ...]]:
    oriented = cyclic_face_set(faces)
    result = []
    for permutation in incidence_automorphisms(adjacency):
        image = cyclic_face_set(
            [[permutation[a], permutation[b], permutation[c]] for a, b, c in faces]
        )
        if image == oriented:
            result.append(permutation)
    return result


def walk_strings(value: Any) -> Iterable[str]:
    if isinstance(value, Mapping):
        for key, child in value.items():
            yield str(key)
            yield from walk_strings(child)
    elif isinstance(value, list):
        for child in value:
            yield from walk_strings(child)
    elif isinstance(value, str):
        yield value


def build_projection(
    carrier_manifest: Mapping[str, Any],
    response_artifact: Mapping[str, Any],
) -> dict[str, Any]:
    carrier = source_carrier(carrier_manifest)
    artifact_body = {
        key: value for key, value in response_artifact.items() if key != "artifact_sha256"
    }
    require(
        response_artifact.get("schema") == "oph.charged_response_semantic_artifact.v3"
        and response_artifact.get("artifact_sha256") == canonical_sha256(artifact_body),
        "RESPONSE_HASH",
        "response artifact hash or schema failed",
    )
    binding = response_artifact.get("carrier_binding")
    require(isinstance(binding, Mapping), "RESPONSE_BINDING", "carrier binding missing")
    require(
        binding.get("port_order") == carrier["ports"]
        and binding.get("carrier_manifest_sha256")
        == canonical_sha256(carrier_manifest).removeprefix("sha256:"),
        "RESPONSE_BINDING",
        "response and carrier bindings disagree",
    )
    source_response = response_artifact.get("source_response")
    require(isinstance(source_response, Mapping), "RESPONSE", "source response missing")
    protocol = source_response.get("impulse_readback_protocol")
    require(isinstance(protocol, Mapping), "RESPONSE", "impulse protocol missing")
    require(
        protocol.get("target_labels_used") is False
        and protocol.get("downstream_labels_used") is False
        and protocol.get("unique_solution_rank") == 4,
        "SOURCE_FIREWALL",
        "response protocol does not pass its source firewall",
    )
    require(
        source_response.get("operator") == "negative_graph_antipode_involution"
        and source_response.get("self_adjoint_unitary_involution") is True
        and source_response.get("commutes_with_propagation_generator") is True,
        "RESPONSE",
        "registered response typing failed",
    )
    provenance = response_artifact.get("provenance")
    runtime = provenance.get("runtime_binding") if isinstance(provenance, Mapping) else None
    require(
        isinstance(runtime, Mapping)
        and runtime.get("coupling") == "combinatorial_icosahedral_graph_laplacian"
        and runtime.get("evolution") == "finite_unitary_exp_minus_i_u_gL"
        and runtime.get("reversible_response_source") == "finite_unitary_carrier_channel",
        "RUNTIME_RESPONSE",
        "registered reversible channel binding failed",
    )

    powers = [identity(12)]
    for _ in range(3):
        powers.append(matmul(powers[-1], carrier["adjacency"]))
    projection: dict[str, Any] = {
        "schema": PROJECTION_SCHEMA,
        "issue": ISSUE,
        "source_files": {
            "carrier": CARRIER_PATH.relative_to(MODULE_DIR).as_posix(),
            "response": RESPONSE_PATH.relative_to(MODULE_DIR).as_posix(),
            "carrier_file_sha256": file_sha256(CARRIER_PATH),
            "response_file_sha256": file_sha256(RESPONSE_PATH),
            "carrier_payload_sha256": canonical_sha256(carrier_manifest),
            "response_payload_sha256": response_artifact["artifact_sha256"],
        },
        "carrier": {
            "port_order": carrier["ports"],
            "edges": carrier["edges"],
            "oriented_faces": carrier["oriented_faces"],
        },
        "registered_recurrence": {
            "kind": "exact_reconstruction_from_pinned_incidence_recurrence",
            "raw_runtime_history_serialized": False,
            "source_count": 12,
            "readback_count": 12,
            "steps": [
                {"step": step, "matrix": render_matrix(matrix)}
                for step, matrix in enumerate(powers)
            ],
            "input": protocol.get("input"),
            "recurrence": protocol.get("recurrence"),
            "filter_coefficients": protocol.get("homogeneous_filter_coefficients"),
            "target_labels_used": False,
            "downstream_labels_used": False,
        },
        "readback": {
            "operator": source_response.get("operator"),
            "source": source_response.get("source"),
            "antipode_port_map": source_response.get("antipode_port_map"),
            "polynomial_identity": source_response.get("antipode_polynomial_identity"),
            "self_adjoint_unitary_involution": source_response.get(
                "self_adjoint_unitary_involution"
            ),
            "commutes_with_recurrence_generator": source_response.get(
                "commutes_with_propagation_generator"
            ),
        },
        "registered_unitary_channel": {
            "coupling": runtime.get("coupling"),
            "evolution": runtime.get("evolution"),
            "generator": "L = 5*I - A",
            "reversible_response_source": runtime.get(
                "reversible_response_source"
            ),
        },
        "registered_refinement_maps": response_artifact.get(
            "physical_refinement_maps", {}
        ).get("port_persistence_maps", []),
    }
    hits = sorted(
        {
            token
            for text in walk_strings(projection)
            for token in FORBIDDEN_SOURCE_TOKENS
            if token in text.lower()
        }
    )
    require(not hits, "SOURCE_FIREWALL", f"forbidden target tokens in projection: {hits}")
    projection["projection_sha256"] = canonical_sha256(projection)
    return projection


def audit_projection(projection: Mapping[str, Any]) -> dict[str, Any]:
    require(projection.get("schema") == PROJECTION_SCHEMA, "PROJECTION_SCHEMA", "wrong projection schema")
    body = {key: value for key, value in projection.items() if key != "projection_sha256"}
    require(projection.get("projection_sha256") == canonical_sha256(body), "PROJECTION_HASH", "projection hash failed")
    carrier = projection.get("carrier")
    recurrence = projection.get("registered_recurrence")
    readback = projection.get("readback")
    unitary_channel = projection.get("registered_unitary_channel")
    require(
        isinstance(carrier, Mapping)
        and isinstance(recurrence, Mapping)
        and isinstance(readback, Mapping)
        and isinstance(unitary_channel, Mapping),
        "PROJECTION",
        "projection blocks missing",
    )
    ports = list(carrier.get("port_order", []))
    require(len(ports) == 12 and len(set(ports)) == 12, "PROJECTION", "wrong port basis")
    adjacency = zero(12)
    for edge in carrier.get("edges", []):
        require(isinstance(edge, list) and len(edge) == 2, "PROJECTION", "malformed edge")
        left, right = (int(edge[0]), int(edge[1]))
        require(0 <= left < right < 12 and adjacency[left][right] == 0, "PROJECTION", "bad edge")
        adjacency[left][right] = adjacency[right][left] = Fraction(1)
    require(all(sum(row) == 5 for row in adjacency), "PROJECTION", "bad degree sequence")
    expected_powers = [identity(12)]
    for _ in range(3):
        expected_powers.append(matmul(expected_powers[-1], adjacency))
    traces = recurrence.get("steps")
    require(isinstance(traces, list) and len(traces) == 4, "TRACE_COUNT", "four recurrence traces required")
    observed_powers = []
    for expected_step, row in enumerate(traces):
        require(isinstance(row, Mapping) and row.get("step") == expected_step, "TRACE_ORDER", "trace order changed")
        observed_powers.append(parse_matrix(row.get("matrix"), 12))
    require(observed_powers == expected_powers, "TRACE_REPLAY", "recurrence trace replay failed")
    require(
        recurrence.get("kind") == "exact_reconstruction_from_pinned_incidence_recurrence"
        and recurrence.get("raw_runtime_history_serialized") is False,
        "TRACE_SCOPE",
        "trace scope is overstated",
    )

    power_rank = rational_rank([flatten(matrix) for matrix in observed_powers])
    require(power_rank == 4, "ALGEBRA_RANK", "response power rank is not four")
    a0, a1, a2, a3 = observed_powers
    a4 = matmul(a3, adjacency)
    minimal_residual = add(
        (Fraction(1), a4),
        (Fraction(-4), a3),
        (Fraction(-10), a2),
        (Fraction(20), a1),
        (Fraction(25), a0),
    )
    require(minimal_residual == zero(12), "MINIMAL_POLYNOMIAL", "degree-four relation failed")
    closure_powers = observed_powers[:]
    for _ in range(4, 7):
        closure_powers.append(matmul(closure_powers[-1], adjacency))
    closure_rank = rational_rank([flatten(matrix) for matrix in closure_powers])
    require(closure_rank == 4, "ALGEBRA_CLOSURE", "response words leave the four-dimensional span")
    laplacian = add((Fraction(5), a0), (Fraction(-1), a1))
    require(
        unitary_channel.get("generator") == "L = 5*I - A"
        and unitary_channel.get("evolution") == "finite_unitary_exp_minus_i_u_gL"
        and rational_rank(basis := [flatten(matrix) for matrix in observed_powers])
        == rational_rank(basis + [flatten(laplacian)]),
        "UNITARY_CHANNEL",
        "registered unitary generator leaves the response-word algebra",
    )
    commutator_nonzero = 0
    for left in observed_powers:
        for right in observed_powers:
            if matmul(left, right) != matmul(right, left):
                commutator_nonzero += 1
    require(commutator_nonzero == 0, "ALGEBRA_COMMUTATOR", "registered response algebra is not commutative")

    gram = [
        [trace_product(left, right) for right in observed_powers]
        for left in observed_powers
    ]
    symmetric_basis = all(transpose(matrix) == matrix for matrix in observed_powers)
    require(
        symmetric_basis,
        "PAIRING",
        "multiplication by i does not give a skew-adjoint response basis",
    )
    principal_minors = [
        determinant([row[:size] for row in gram[:size]])
        for size in range(1, 5)
    ]
    require(all(value > 0 for value in principal_minors), "PAIRING", "skew pairing is not positive definite")

    antipode = readback.get("antipode_port_map")
    require(
        isinstance(antipode, list)
        and sorted(int(value) for value in antipode) == list(range(12))
        and all(int(antipode[int(antipode[i])]) == i for i in range(12)),
        "READBACK",
        "antipode is not a permutation involution",
    )
    response = add(
        (Fraction(-1), a0),
        (Fraction(1, 2), a1),
        (Fraction(2, 5), a2),
        (Fraction(-1, 10), a3),
    )
    expected_response = zero(12)
    for source, target in enumerate(antipode):
        expected_response[source][int(target)] = Fraction(-1)
    require(response == expected_response, "READBACK", "negative antipode polynomial failed")
    require(matmul(response, response) == a0 and transpose(response) == response, "READBACK", "response is not a self-adjoint involution")

    faces = carrier.get("oriented_faces")
    require(isinstance(faces, list) and len(faces) == 20, "RECHARTING", "oriented faces missing")
    automorphisms = incidence_automorphisms(adjacency)
    rotations = proper_rechartings(adjacency, faces)
    require(len(automorphisms) == 120 and len(rotations) == 60, "RECHARTING", "carrier automorphism count failed")
    basis_rows = [flatten(matrix) for matrix in observed_powers]
    in_response_algebra = []
    all_commute = True
    for permutation in rotations:
        matrix = permutation_matrix(permutation)
        all_commute = all_commute and matmul(matrix, adjacency) == matmul(adjacency, matrix)
        if rational_rank(basis_rows + [flatten(matrix)]) == power_rank:
            in_response_algebra.append(permutation)
    require(all_commute, "RECHARTING", "a proper recharting does not commute with incidence")
    require(
        in_response_algebra == [tuple(range(12))],
        "RECHARTING",
        "a nonidentity proper recharting lies in the registered response-word algebra",
    )

    refinement_rows = projection.get("registered_refinement_maps")
    require(isinstance(refinement_rows, list) and refinement_rows, "REFINEMENT", "refinement maps missing")
    refinement_natural = True
    for row in refinement_rows:
        require(isinstance(row, Mapping), "REFINEMENT", "malformed refinement row")
        port_map = row.get("port_map")
        require(isinstance(port_map, list) and sorted(port_map) == list(range(12)), "REFINEMENT", "bad port map")
        row_body = {key: value for key, value in row.items() if key != "map_hash"}
        require(
            row.get("map_hash") == canonical_sha256(row_body),
            "REFINEMENT",
            "refinement map hash failed",
        )
        p_matrix = permutation_matrix(port_map)
        refinement_natural = refinement_natural and matmul(p_matrix, adjacency) == matmul(adjacency, p_matrix)
    require(refinement_natural, "REFINEMENT", "registered recurrence is not refinement-natural")

    classifications = {
        "ordered_two_sided_port_response_histories": "PARTIAL_SOURCE_NATIVE_COMMUTATIVE_RECURRENCE_ONLY",
        "twelve_independent_skew_adjoint_generators": "IMPOSSIBLE_ON_REGISTERED_RESPONSE_WORD_ALGEBRA",
        "exact_nonabelian_commutator_reconstruction": "IMPOSSIBLE_ON_REGISTERED_RESPONSE_WORD_ALGEBRA",
        "closed_overlap_words_cover_sixty_rechartings": "MISSING_STATIC_RECHARTINGS_ONLY",
        "same_word_projective_implementers": "MISSING",
        "response_identity_component_factorizations": "IMPOSSIBLE_ON_REGISTERED_RESPONSE_WORD_ALGEBRA",
        "refinement_intertwining": "PARTIAL_SOURCE_NATIVE_RECURRENCE_ONLY",
        "adversarial_controls": "PARTIAL_EXACT_PACKET_CONTROLS",
        "target_firewall": "SOURCE_NATIVE_ATTAINED",
    }
    return {
        "response_word_algebra": {
            "basis": ["I", "A", "A^2", "A^3"],
            "exact_dimension": power_rank,
            "closure_rank_through_degree_six": closure_rank,
            "minimal_polynomial": "A^4 - 4*A^3 - 10*A^2 + 20*A + 25*I = 0",
            "commutator_nonzero_count": commutator_nonzero,
            "commutative": True,
            "skew_adjoint_real_form": "i times span_Q{I,A,A^2,A^3}",
            "skew_adjoint_basis_verified": symmetric_basis,
            "skew_pairing_gram": [[render_fraction(value) for value in row] for row in gram],
            "skew_pairing_leading_principal_minors": [
                render_fraction(value) for value in principal_minors
            ],
            "skew_pairing_positive_definite": True,
            "maximum_independent_skew_generator_count": 4,
        },
        "registered_unitary_channel_audit": {
            "generator": "L = 5*I - A",
            "generator_in_response_word_algebra": True,
            "single_generator_functional_calculus_dimension_upper_bound": 4,
            "order_sensitive_port_indexed_tangent_available": False,
        },
        "recharting_audit": {
            "incidence_automorphism_count": len(automorphisms),
            "proper_recharting_count": len(rotations),
            "proper_rechartings_commute_with_recurrence": all_commute,
            "proper_rechartings_in_response_word_algebra": len(in_response_algebra),
            "nonidentity_proper_rechartings_in_response_word_algebra": 0,
            "static_rechartings_are_ordered_response_words": False,
        },
        "refinement_audit": {
            "registered_map_count": len(refinement_rows),
            "recurrence_natural_on_registered_maps": refinement_natural,
            "generator_bracket_implementer_intertwining_available": False,
        },
        "acceptance_classifications": classifications,
    }


def certificate_payload(projection: Mapping[str, Any]) -> dict[str, Any]:
    audit = audit_projection(projection)
    receipt: dict[str, Any] = {
        "schema": RECEIPT_SCHEMA,
        "issue": ISSUE,
        "projection_sha256": projection["projection_sha256"],
        "source_files": projection["source_files"],
        "verdict": VERDICT,
        "physical_current_source_bridge_attained": False,
        "abstract_forced_lie_type_theorem_preserved": True,
        "audit_scope": (
            "the algebra generated by the single registered adjacency recurrence "
            "and its polynomial maximal-distance readback; no bound is asserted for "
            "an unregistered order-sensitive response history"
        ),
        "audit": audit,
        "claim_boundary": (
            "The pinned adjacency recurrence and negative-antipode readback generate "
            "a four-dimensional commutative response-word algebra. This exact packet "
            "cannot supply twelve independent generators, a nonzero bracket, or the "
            "nonidentity proper rechartings. The result applies only to the registered "
            "recurrence packet. It leaves an order-sensitive A1/A2 response producer "
            "open and does not weaken the abstract conditional Lie-type theorem."
        ),
        "next_source_object": {
            "required": (
                "port-indexed reversible perturbations with both composition orders "
                "recorded on one carrier source"
            ),
            "first_exact_gate": (
                "the twelve first-order port derivatives have real rank twelve; "
                "their exact mixed-order commutators close in that span, their "
                "commutator span has real rank eleven, and the constant linear "
                "combination of the twelve port generators spans the one-dimensional "
                "center"
            ),
            "holonomy_gate": (
                "closed overlap words from those perturbations cover all sixty proper "
                "rechartings and reconstruct their implementers"
            ),
        },
        "verifier_command": (
            "python3 code/a5_closure/verify_source_current_capability_independent.py "
            "--projection code/a5_closure/manifests/source_current_capability_projection.json "
            "--receipt code/a5_closure/receipts/source_current_capability.receipt.json"
        ),
    }
    receipt["receipt_sha256"] = canonical_sha256(receipt)
    return receipt


def verify_receipt(projection: Mapping[str, Any], receipt: Mapping[str, Any]) -> None:
    require(receipt.get("schema") == RECEIPT_SCHEMA, "RECEIPT_SCHEMA", "wrong receipt schema")
    body = {key: value for key, value in receipt.items() if key != "receipt_sha256"}
    require(receipt.get("receipt_sha256") == canonical_sha256(body), "RECEIPT_HASH", "receipt hash failed")
    require(dict(receipt) == certificate_payload(projection), "RECEIPT_REPLAY", "receipt does not replay exactly")


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("all")
    verify = subparsers.add_parser("verify")
    verify.add_argument("--projection", type=Path, default=PROJECTION_PATH)
    verify.add_argument("--receipt", type=Path, default=RECEIPT_PATH)
    args = parser.parse_args()

    if args.command == "all":
        projection = build_projection(load_json(CARRIER_PATH), load_json(RESPONSE_PATH))
        receipt = certificate_payload(projection)
        write_json(PROJECTION_PATH, projection)
        write_json(RECEIPT_PATH, receipt)
        print(VERDICT)
        return 0
    projection = load_json(args.projection)
    receipt = load_json(args.receipt)
    verify_receipt(projection, receipt)
    print("SOURCE_CURRENT_CAPABILITY_VALID")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
