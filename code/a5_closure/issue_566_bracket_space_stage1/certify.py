#!/usr/bin/env python3
"""Exact, target-free first-stage certificate for issue #566.

The sole semantic input is the canonical oriented twelve-port carrier manifest.
From its edges and oriented faces this program reconstructs the full incidence
automorphism group, selects the sixty orientation-preserving permutations,
and certifies the rational invariant space

    Hom_G(exterior_square(Q^12), Q^12).

Here ``target-free'' means target-free after conditioning on the pinned
canonical oriented twelve-port carrier.  No candidate bracket, later response
fixture, measurement, desired gauge algebra, or physical target is accepted as
input.  The packet does not derive the carrier choice.  Its output is only a
basis for the entire equivariant alternating search space; Jacobi and every
later selection gate remain open.
"""

from __future__ import annotations

import argparse
import ast
import copy
import hashlib
import itertools
import json
import sys
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
CARRIER_PATH = REPO / "code/a5_closure/manifests/echosahedral_federation_reference.json"
BASIS_PATH = HERE / "a5_alternating_bracket_reynolds_basis.json"
RECEIPT_PATH = HERE / "a5_alternating_bracket_space_stage1.receipt.json"
VERIFIER_PATH = HERE / "verify.py"
TEST_PATH = HERE / "test_stage1.py"

PORTS = 12
EXPECTED_TOP_KEYS = {
    "architecture",
    "carrier",
    "refinement_tower",
    "schema",
    "source_readback",
}
EXPECTED_CARRIER_KEYS = {
    "atoms_pairwise_orthogonal",
    "atoms_sum_to_one",
    "central_port_atoms",
    "edges",
    "oriented_faces",
    "ports",
}
PRODUCER_ALLOWED_IMPORT_ROOTS = {
    "__future__",
    "argparse",
    "ast",
    "copy",
    "dataclasses",
    "fractions",
    "hashlib",
    "itertools",
    "json",
    "pathlib",
    "sys",
    "typing",
}
VERIFIER_ALLOWED_IMPORT_ROOTS = {
    "__future__",
    "ast",
    "fractions",
    "hashlib",
    "itertools",
    "json",
    "pathlib",
    "sys",
    "typing",
}
TEST_ALLOWED_IMPORT_ROOTS = {
    "__future__",
    "hashlib",
    "json",
    "pathlib",
    "shutil",
    "subprocess",
    "sys",
    "tempfile",
    "unittest",
}

CLAIM_BOUNDARY = (
    "Conditional on the pinned canonical oriented twelve-port carrier, this receipt certifies only the complete "
    "equivariant alternating-bracket search space on its permutation module. Target-free means that no desired "
    "gauge algebra, measurement, or coefficient target enters this stage; it does not mean that this packet derives "
    "the carrier choice. The receipt does not select a bracket, solve Jacobi, establish compactness, identify a "
    "physical current, or close issue #566."
)

Permutation = tuple[int, ...]
Coordinate = tuple[int, int, int]  # output, left input, right input; left < right


class CertificateError(RuntimeError):
    """A fail-closed certificate error with a stable machine code."""

    def __init__(self, code: str, message: str):
        super().__init__(f"{code}: {message}")
        self.code = code


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise CertificateError(code, message)


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")


def sha256_bytes(payload: bytes) -> str:
    return "sha256:" + hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def cyclic_face(face: Sequence[int]) -> tuple[int, int, int]:
    a, b, c = face
    return min((a, b, c), (b, c, a), (c, a, b))


@dataclass(frozen=True)
class Carrier:
    ports: tuple[str, ...]
    edges: frozenset[tuple[int, int]]
    faces: frozenset[tuple[int, int, int]]
    adjacency: tuple[tuple[bool, ...], ...]
    projection: Mapping[str, Any]


def parse_manifest(raw: Mapping[str, Any]) -> Carrier:
    """Validate the fixed schema and retain only incidence/orientation data."""

    require(isinstance(raw, Mapping), "FIREWALL_MANIFEST_TYPE", "carrier manifest must be an object")
    require(set(raw) == EXPECTED_TOP_KEYS, "FIREWALL_TOP_KEYS", "unexpected top-level semantic input")
    require(
        raw.get("schema") == "oph.echosahedral_selector_manifest.v1",
        "FIREWALL_SCHEMA",
        "unexpected carrier schema",
    )
    body = raw.get("carrier")
    require(isinstance(body, Mapping), "FIREWALL_CARRIER_TYPE", "carrier section must be an object")
    require(
        set(body) == EXPECTED_CARRIER_KEYS,
        "FIREWALL_CARRIER_KEYS",
        "unexpected carrier field; target injection is not accepted",
    )

    ports_raw = body.get("ports")
    require(
        isinstance(ports_raw, list) and len(ports_raw) == PORTS,
        "PORT_COUNT",
        "exactly twelve ordered port identifiers are required",
    )
    require(all(type(x) is str for x in ports_raw), "PORT_LABEL", "port identifiers must be strings")
    ports = tuple(ports_raw)
    require(len(set(ports)) == PORTS, "PORT_LABEL", "port identifiers must be distinct")
    index = {name: i for i, name in enumerate(ports)}

    edges_raw = body.get("edges")
    require(isinstance(edges_raw, list) and len(edges_raw) == 30, "EDGE_COUNT", "expected thirty edges")
    edges: set[tuple[int, int]] = set()
    for row in edges_raw:
        require(isinstance(row, list) and len(row) == 2, "EDGE_FORMAT", "each edge must have two ports")
        require(all(type(x) is str and x in index for x in row), "EDGE_PORT", "edge names an unknown port")
        i, j = index[row[0]], index[row[1]]
        require(i != j, "EDGE_LOOP", "carrier edges cannot be loops")
        edge = (min(i, j), max(i, j))
        require(edge not in edges, "EDGE_DUPLICATE", "carrier edge is duplicated")
        edges.add(edge)

    adjacency = [[False for _ in range(PORTS)] for _ in range(PORTS)]
    for i, j in edges:
        adjacency[i][j] = adjacency[j][i] = True
    require(all(sum(row) == 5 for row in adjacency), "EDGE_DEGREE", "carrier graph is not five-regular")
    reached = {0}
    frontier = [0]
    while frontier:
        i = frontier.pop()
        for j in range(PORTS):
            if adjacency[i][j] and j not in reached:
                reached.add(j)
                frontier.append(j)
    require(len(reached) == PORTS, "EDGE_CONNECTED", "carrier graph is disconnected")

    faces_raw = body.get("oriented_faces")
    require(isinstance(faces_raw, list) and len(faces_raw) == 20, "FACE_COUNT", "expected twenty faces")
    faces: set[tuple[int, int, int]] = set()
    directed_edges: dict[tuple[int, int], int] = {}
    for row in faces_raw:
        require(isinstance(row, list) and len(row) == 3, "FACE_FORMAT", "each face must have three ports")
        require(all(type(x) is str and x in index for x in row), "FACE_PORT", "face names an unknown port")
        triple = tuple(index[x] for x in row)
        require(len(set(triple)) == 3, "FACE_REPEAT", "face repeats a port")
        face = cyclic_face(triple)
        require(face not in faces, "FACE_DUPLICATE", "oriented face is duplicated")
        for u, v in ((triple[0], triple[1]), (triple[1], triple[2]), (triple[2], triple[0])):
            require((min(u, v), max(u, v)) in edges, "FACE_NONEDGE", "face contains a non-edge")
            directed_edges[(u, v)] = directed_edges.get((u, v), 0) + 1
        faces.add(face)
    require(
        all(directed_edges.get((i, j), 0) == 1 and directed_edges.get((j, i), 0) == 1 for i, j in edges),
        "FACE_ORIENTATION",
        "oriented faces do not give opposite directions on every shared edge",
    )
    require(PORTS - len(edges) + len(faces) == 2, "FACE_EULER", "carrier triangulation has wrong Euler count")

    projection = {
        "schema": raw["schema"],
        "ports": list(ports),
        "edges": [list(edge) for edge in sorted(edges)],
        "oriented_faces": [list(face) for face in sorted(faces)],
    }
    return Carrier(
        ports=ports,
        edges=frozenset(edges),
        faces=frozenset(faces),
        adjacency=tuple(tuple(row) for row in adjacency),
        projection=projection,
    )


def load_canonical_carrier() -> tuple[dict[str, Any], Carrier]:
    require(CARRIER_PATH.resolve().is_file(), "FIREWALL_PATH", "canonical carrier manifest is missing")
    raw = json.loads(CARRIER_PATH.read_text(encoding="utf-8"))
    require(isinstance(raw, dict), "FIREWALL_MANIFEST_TYPE", "carrier manifest must be an object")
    return raw, parse_manifest(raw)


def compose(p: Permutation, q: Permutation) -> Permutation:
    """Composition p after q."""

    return tuple(p[q[i]] for i in range(len(p)))


def inverse(p: Permutation) -> Permutation:
    out = [0] * len(p)
    for i, image in enumerate(p):
        out[image] = i
    return tuple(out)


def permutation_order(p: Permutation) -> int:
    identity = tuple(range(len(p)))
    value = identity
    for order in range(1, 1 + 2 * len(p) * len(p)):
        value = compose(p, value)
        if value == identity:
            return order
    raise CertificateError("GROUP_ORDER", "permutation order bound exceeded")


def enumerate_incidence_automorphisms(carrier: Carrier) -> list[Permutation]:
    """Backtrack all graph automorphisms without a preloaded permutation list."""

    result: list[Permutation] = []

    def extend(partial: list[int]) -> None:
        source = len(partial)
        if source == PORTS:
            result.append(tuple(partial))
            return
        used = set(partial)
        for target in range(PORTS):
            if target in used:
                continue
            if all(
                carrier.adjacency[source][old_source]
                == carrier.adjacency[target][partial[old_source]]
                for old_source in range(source)
            ):
                extend(partial + [target])

    extend([])
    require(len(result) == 120 and len(set(result)) == 120, "INCIDENCE_AUTOMORPHISMS", "expected 120 graph automorphisms")
    return sorted(result)


def mapped_faces(p: Permutation, faces: Iterable[tuple[int, int, int]]) -> frozenset[tuple[int, int, int]]:
    return frozenset(cyclic_face((p[a], p[b], p[c])) for a, b, c in faces)


def validate_group(group: Sequence[Permutation], carrier: Carrier) -> None:
    require(len(group) == 60 and len(set(group)) == 60, "GROUP_COUNT", "proper action must have sixty distinct rows")
    identity = tuple(range(PORTS))
    require(identity in group, "GROUP_IDENTITY", "identity is absent")
    group_set = set(group)
    for p in group:
        require(len(p) == PORTS and set(p) == set(range(PORTS)), "PERMUTATION_BIJECTION", "group row is not a permutation")
        require(mapped_faces(p, carrier.faces) == carrier.faces, "GROUP_ORIENTATION", "group row reverses carrier orientation")
        for i in range(PORTS):
            for j in range(PORTS):
                require(
                    carrier.adjacency[p[i]][p[j]] == carrier.adjacency[i][j],
                    "GROUP_INCIDENCE",
                    "group row does not preserve incidence",
                )
    for p in group:
        require(inverse(p) in group_set, "GROUP_INVERSE", "inverse is absent")
        for q in group:
            require(compose(p, q) in group_set, "GROUP_CLOSURE", "proper action is not closed")
    require(
        all(any(p[0] == target for p in group) for target in range(PORTS)),
        "GROUP_TRANSITIVITY",
        "proper action is not transitive",
    )


def proper_port_group(carrier: Carrier) -> tuple[list[Permutation], int]:
    full = enumerate_incidence_automorphisms(carrier)
    proper = [p for p in full if mapped_faces(p, carrier.faces) == carrier.faces]
    reversing = [p for p in full if p not in set(proper)]
    require(len(proper) == 60 and len(reversing) == 60, "GROUP_ORIENTATION_SPLIT", "orientation split is not 60 plus 60")
    validate_group(proper, carrier)
    return sorted(proper), len(full)


def generated_group(generators: Sequence[Permutation]) -> set[Permutation]:
    identity = tuple(range(len(generators[0])))
    moves = list(generators) + [inverse(g) for g in generators]
    reached = {identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for move in moves:
            candidate = compose(move, current)
            if candidate not in reached:
                reached.add(candidate)
                frontier.append(candidate)
    return reached


def is_even(p: Permutation) -> bool:
    return sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p))) % 2 == 0


def canonical_generating_pair(group: Sequence[Permutation]) -> tuple[Permutation, Permutation]:
    for a in sorted(group):
        if permutation_order(a) != 2:
            continue
        for b in sorted(group):
            if permutation_order(b) != 3 or permutation_order(compose(a, b)) != 5:
                continue
            if len(generated_group((a, b))) == 60:
                return a, b
    raise CertificateError("A5_GENERATORS", "no canonical (2,3,5) generating pair was found")


def certify_a5_isomorphism(group: Sequence[Permutation]) -> dict[str, Any]:
    """Exhibit an exact isomorphism with the sixty even permutations on five letters."""

    even5 = sorted(tuple(p) for p in itertools.permutations(range(5)) if is_even(tuple(p)))
    require(len(even5) == 60, "A5_MODEL", "even five-letter permutation model has wrong order")
    port_a, port_b = canonical_generating_pair(group)
    model_a, model_b = canonical_generating_pair(even5)
    moves = [
        (port_a, model_a),
        (port_b, model_b),
        (inverse(port_a), inverse(model_a)),
        (inverse(port_b), inverse(model_b)),
    ]
    port_identity = tuple(range(PORTS))
    model_identity = tuple(range(5))
    mapping: dict[Permutation, Permutation] = {port_identity: model_identity}
    reverse: dict[Permutation, Permutation] = {model_identity: port_identity}
    frontier = [port_identity]
    while frontier:
        current = frontier.pop()
        current_model = mapping[current]
        for port_move, model_move in moves:
            next_port = compose(port_move, current)
            next_model = compose(model_move, current_model)
            if next_port in mapping:
                require(mapping[next_port] == next_model, "A5_RELATION", "generator word has inconsistent image")
            else:
                require(next_model not in reverse, "A5_INJECTIVITY", "two port rows map to one model row")
                mapping[next_port] = next_model
                reverse[next_model] = next_port
                frontier.append(next_port)
    require(set(mapping) == set(group) and set(mapping.values()) == set(even5), "A5_SURJECTIVITY", "isomorphism does not cover both groups")
    for p in group:
        for q in group:
            require(
                mapping[compose(p, q)] == compose(mapping[p], mapping[q]),
                "A5_HOMOMORPHISM",
                "displayed bijection does not preserve multiplication",
            )
    serialized_map = [[list(p), list(mapping[p])] for p in sorted(mapping)]
    return {
        "model": "even permutations on five unlabeled symbols",
        "model_order": len(even5),
        "bijection_size": len(mapping),
        "port_generators": [list(port_a), list(port_b)],
        "model_generators": [list(model_a), list(model_b)],
        "generator_orders": [2, 3, 5],
        "homomorphism_checked_pairs": len(group) ** 2,
        "mapping_sha256": sha256_bytes(canonical_bytes(serialized_map)),
        "isomorphic": True,
    }


def all_coordinates() -> list[Coordinate]:
    return [(k, i, j) for k in range(PORTS) for i in range(PORTS) for j in range(i + 1, PORTS)]


def act_coordinate(p: Permutation, coordinate: Coordinate) -> tuple[Coordinate, int]:
    k, i, j = coordinate
    left, right = p[i], p[j]
    sign = 1 if left < right else -1
    return (p[k], min(left, right), max(left, right)), sign


def character_certificate(group: Sequence[Permutation]) -> dict[str, Any]:
    rows: dict[tuple[int, int, int], int] = {}
    numerator = 0
    for p in group:
        chi = sum(p[i] == i for i in range(PORTS))
        square = compose(p, p)
        chi_square = sum(square[i] == i for i in range(PORTS))
        wedge_character = (chi * chi - chi_square) // 2
        key = (permutation_order(p), chi, chi_square)
        rows[key] = rows.get(key, 0) + 1
        numerator += chi * wedge_character
    require(numerator % len(group) == 0, "CHARACTER_INTEGRAL", "projector trace is not integral")
    dimension = numerator // len(group)
    table = []
    for (order, chi, chi_square), count in sorted(rows.items()):
        wedge_character = (chi * chi - chi_square) // 2
        contribution = chi * wedge_character
        table.append(
            {
                "element_order": order,
                "count": count,
                "chi_V": chi,
                "chi_V_of_square": chi_square,
                "chi_exterior_square_V": wedge_character,
                "projector_trace_contribution_each": contribution,
                "projector_trace_contribution_total": count * contribution,
            }
        )
    return {
        "formula": "sum_g chi_V(g)*(chi_V(g)^2-chi_V(g^2))/2 divided by group_order",
        "class_rows": table,
        "projector_trace_numerator": numerator,
        "projector_trace_denominator": len(group),
        "dimension": dimension,
    }


def build_reynolds_basis(group: Sequence[Permutation]) -> list[dict[Coordinate, Fraction]]:
    coordinates = all_coordinates()
    unassigned = set(coordinates)
    vectors: list[dict[Coordinate, Fraction]] = []
    while unassigned:
        seed = min(unassigned)
        sums: dict[Coordinate, Fraction] = {}
        signs: dict[Coordinate, int] = {}
        conflict = False
        for p in group:
            image, sign = act_coordinate(p, seed)
            sums[image] = sums.get(image, Fraction(0)) + Fraction(sign, len(group))
            if image in signs and signs[image] != sign:
                conflict = True
            signs[image] = sign
        orbit = set(signs)
        unassigned -= orbit
        if conflict:
            require(all(value == 0 for value in sums.values()), "REYNOLDS_CONFLICT", "signed stabilizer cancellation is incomplete")
            continue
        vector = {coordinate: value for coordinate, value in sums.items() if value}
        require(set(vector) == orbit, "REYNOLDS_SUPPORT", "Reynolds vector lost an orientable orbit coordinate")
        orbit_size = len(orbit)
        require(
            all(abs(value) == Fraction(1, orbit_size) for value in vector.values()),
            "REYNOLDS_NORMALIZATION",
            "unexpected Reynolds coefficient",
        )
        vectors.append(dict(sorted(vector.items())))
    return vectors


def check_covariance(vectors: Sequence[Mapping[Coordinate, Fraction]], group: Sequence[Permutation]) -> int:
    checks = 0
    coordinates = all_coordinates()
    for vector in vectors:
        for p in group:
            for coordinate in coordinates:
                image, sign = act_coordinate(p, coordinate)
                require(
                    vector.get(image, Fraction(0)) == sign * vector.get(coordinate, Fraction(0)),
                    "COVARIANCE",
                    "basis vector is not equivariant",
                )
                checks += 1
    return checks


def check_alternation(vectors: Sequence[Mapping[Coordinate, Fraction]]) -> int:
    """Expand the wedge convention and check B(e_i,e_j)=-B(e_j,e_i)."""

    checks = 0
    for vector in vectors:
        for k in range(PORTS):
            for i in range(PORTS):
                for j in range(PORTS):
                    def coefficient(left: int, right: int) -> Fraction:
                        if left == right:
                            return Fraction(0)
                        if left < right:
                            return vector.get((k, left, right), Fraction(0))
                        return -vector.get((k, right, left), Fraction(0))

                    require(coefficient(i, j) == -coefficient(j, i), "ALTERNATION", "wedge expansion is not alternating")
                    checks += 1
    return checks


def rational_rank(matrix: Sequence[Sequence[Fraction]]) -> int:
    if not matrix:
        return 0
    work = [list(row) for row in matrix]
    row = 0
    columns = len(work[0])
    for column in range(columns):
        pivot = next((r for r in range(row, len(work)) if work[r][column]), None)
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        scale = work[row][column]
        work[row] = [value / scale for value in work[row]]
        for r in range(len(work)):
            if r != row and work[r][column]:
                factor = work[r][column]
                work[r] = [x - factor * y for x, y in zip(work[r], work[row])]
        row += 1
        if row == len(work):
            break
    return row


def determinant(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    require(all(len(row) == len(matrix) for row in matrix), "DETERMINANT_SHAPE", "minor is not square")
    work = [list(row) for row in matrix]
    result = Fraction(1)
    for column in range(len(work)):
        pivot = next((r for r in range(column, len(work)) if work[r][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            result = -result
        value = work[column][column]
        result *= value
        for r in range(column + 1, len(work)):
            if work[r][column]:
                factor = work[r][column] / value
                for c in range(column, len(work)):
                    work[r][c] -= factor * work[column][c]
    return result


def fraction_json(value: Fraction) -> dict[str, int]:
    return {"numerator": value.numerator, "denominator": value.denominator}


def rank_certificate(vectors: Sequence[Mapping[Coordinate, Fraction]]) -> dict[str, Any]:
    coordinates = all_coordinates()
    matrix = [[vector.get(coordinate, Fraction(0)) for vector in vectors] for coordinate in coordinates]
    rank = rational_rank(matrix)
    require(rank == 14, "BASIS_RANK", f"expected rank fourteen, got {rank}")
    seeds = [min(vector) for vector in vectors]
    minor = [[vectors[column].get(seed, Fraction(0)) for column in range(len(vectors))] for seed in seeds]
    minor_det = determinant(minor)
    require(minor_det != 0, "BASIS_MINOR", "canonical seed minor is singular")
    gram = [
        [
            sum((left.get(c, Fraction(0)) * right.get(c, Fraction(0)) for c in coordinates), Fraction(0))
            for right in vectors
        ]
        for left in vectors
    ]
    require(all(gram[i][j] == 0 for i in range(len(gram)) for j in range(len(gram)) if i != j), "BASIS_GRAM", "orbit basis is not orthogonal")
    gram_diagonal = [gram[i][i] for i in range(len(gram))]
    gram_det = Fraction(1)
    for value in gram_diagonal:
        gram_det *= value
    require(gram_det != 0, "BASIS_GRAM", "Gram determinant vanishes")
    return {
        "ambient_coordinate_count": len(coordinates),
        "matrix_shape": [len(coordinates), len(vectors)],
        "rank_over_Q": rank,
        "pivot_coordinates": [list(seed) for seed in seeds],
        "pivot_minor_determinant": fraction_json(minor_det),
        "gram_diagonal": [fraction_json(value) for value in gram_diagonal],
        "gram_determinant": fraction_json(gram_det),
    }


def encode_basis(vectors: Sequence[Mapping[Coordinate, Fraction]], group_sha256: str) -> dict[str, Any]:
    rows = []
    for number, vector in enumerate(vectors):
        seed = min(vector)
        orbit_size = len(vector)
        entries = [
            [coordinate[0], coordinate[1], coordinate[2], value.numerator, value.denominator]
            for coordinate, value in sorted(vector.items())
        ]
        rows.append(
            {
                "basis_id": f"R{number:02d}",
                "seed_coordinate": list(seed),
                "orbit_size": orbit_size,
                "stabilizer_size": 60 // orbit_size,
                "integral_orbit_scale": orbit_size,
                "entries": entries,
            }
        )
    return {
        "schema": "oph.a5_alternating_bracket_reynolds_basis.v1",
        "issue": 566,
        "field": "Q",
        "port_dimension": PORTS,
        "domain_dimension": PORTS * (PORTS - 1) // 2,
        "ambient_coordinate_count": len(all_coordinates()),
        "coordinate_convention": "[output,left,right] with left<right; swapping inputs negates the coefficient",
        "group_action_convention": "(k,i,j) maps to (g(k),sort(g(i),g(j))) with the sorting sign",
        "normalization": "each vector is the exact Reynolds average (1/60)*sum_g g(seed)",
        "proper_port_group_sha256": group_sha256,
        "basis": rows,
    }


def audit_imports(paths: Sequence[Path]) -> dict[str, list[str]]:
    allowlists = {
        Path(__file__).resolve(): PRODUCER_ALLOWED_IMPORT_ROOTS,
        VERIFIER_PATH.resolve(): VERIFIER_ALLOWED_IMPORT_ROOTS,
        TEST_PATH.resolve(): TEST_ALLOWED_IMPORT_ROOTS,
    }
    require(set(map(Path.resolve, paths)) == set(allowlists), "FIREWALL_IMPORT_PATH", "unexpected audited source set")
    result: dict[str, list[str]] = {}
    for path in paths:
        imports: set[str] = set()
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imports.add(node.module.split(".")[0])
        allowed = allowlists[path.resolve()]
        require(
            imports <= allowed,
            "FIREWALL_IMPORT",
            f"non-allowlisted imports in {path.name}: {sorted(imports - allowed)}",
        )
        result[path.name] = sorted(imports)
    return dict(sorted(result.items()))


def expect_error(name: str, expected: str, action: Any) -> dict[str, Any]:
    actual = "NO_ERROR"
    try:
        action()
    except CertificateError as exc:
        actual = exc.code
    return {"name": name, "expected_error": expected, "actual_error": actual, "passed": actual == expected}


def run_mutation_tests(
    raw: Mapping[str, Any], carrier: Carrier, group: Sequence[Permutation], vectors: Sequence[Mapping[Coordinate, Fraction]]
) -> list[dict[str, Any]]:
    tests: list[dict[str, Any]] = []

    missing_edge = copy.deepcopy(raw)
    missing_edge["carrier"]["edges"].pop()
    tests.append(expect_error("remove_carrier_edge", "EDGE_COUNT", lambda: parse_manifest(missing_edge)))

    reversed_face = copy.deepcopy(raw)
    reversed_face["carrier"]["oriented_faces"][0][1], reversed_face["carrier"]["oriented_faces"][0][2] = (
        reversed_face["carrier"]["oriented_faces"][0][2],
        reversed_face["carrier"]["oriented_faces"][0][1],
    )
    tests.append(expect_error("reverse_one_oriented_face", "FACE_ORIENTATION", lambda: parse_manifest(reversed_face)))

    injected = copy.deepcopy(raw)
    injected["carrier"]["target_payload"] = {"accepted": True}
    tests.append(expect_error("inject_unregistered_semantic_target", "FIREWALL_CARRIER_KEYS", lambda: parse_manifest(injected)))

    bad_group = list(group)
    row = list(bad_group[1])
    row[-1] = row[-2]
    bad_group[1] = tuple(row)
    tests.append(expect_error("break_permutation_bijection", "PERMUTATION_BIJECTION", lambda: validate_group(bad_group, carrier)))

    sign_flip = [dict(vector) for vector in vectors]
    coordinate = min(sign_flip[0])
    sign_flip[0][coordinate] = -sign_flip[0][coordinate]
    tests.append(expect_error("flip_one_basis_coefficient", "COVARIANCE", lambda: check_covariance(sign_flip, group)))

    dropped = [dict(vector) for vector in vectors[:-1]]
    tests.append(expect_error("drop_one_basis_vector", "BASIS_RANK", lambda: rank_certificate(dropped)))

    changed_scale = [dict(vector) for vector in vectors]
    coordinate = min(changed_scale[0])
    changed_scale[0][coordinate] += Fraction(1, 60)
    tests.append(expect_error("alter_reynolds_normalization", "COVARIANCE", lambda: check_covariance(changed_scale, group)))

    require(all(row["passed"] for row in tests), "MUTATION_TEST", "one or more mutation controls failed")
    return tests


def build_artifacts() -> tuple[dict[str, Any], dict[str, Any]]:
    raw, carrier = load_canonical_carrier()
    group, full_order = proper_port_group(carrier)
    group_rows = [list(row) for row in group]
    group_sha256 = sha256_bytes(canonical_bytes(group_rows))
    a5_witness = certify_a5_isomorphism(group)
    character = character_certificate(group)
    require(character["dimension"] == 14, "CHARACTER_DIMENSION", "equivariant dimension is not fourteen")

    vectors = build_reynolds_basis(group)
    require(len(vectors) == 14, "BASIS_COUNT", "Reynolds construction did not yield fourteen vectors")
    require(sum(len(vector) for vector in vectors) == len(all_coordinates()), "BASIS_PARTITION", "basis supports do not partition all coordinates")
    alternation_checks = check_alternation(vectors)
    covariance_checks = check_covariance(vectors, group)
    rank = rank_certificate(vectors)
    mutation_tests = run_mutation_tests(raw, carrier, group, vectors)
    imported_by_file = audit_imports((Path(__file__).resolve(), VERIFIER_PATH.resolve(), TEST_PATH.resolve()))
    imported = sorted({root for roots in imported_by_file.values() for root in roots})

    basis = encode_basis(vectors, group_sha256)
    basis_sha256 = sha256_bytes(canonical_bytes(basis))
    orbit_sizes = sorted(len(vector) for vector in vectors)
    receipt: dict[str, Any] = {
        "schema": "oph.a5_alternating_bracket_space_stage1.receipt.v1",
        "issue": 566,
        "status": "EXACT_TARGET_FREE_SEARCH_SPACE_STAGE1_ONLY",
        "claim_boundary": CLAIM_BOUNDARY,
        "semantic_input": {
            "path": "code/a5_closure/manifests/echosahedral_federation_reference.json",
            "file_sha256": sha256_file(CARRIER_PATH),
            "carrier_projection_sha256": sha256_bytes(canonical_bytes(carrier.projection)),
            "accessed_value_paths": [
                "/schema",
                "/carrier/ports",
                "/carrier/edges",
                "/carrier/oriented_faces",
            ],
            "ignored_value_sections": [
                "/architecture",
                "/carrier/atoms_pairwise_orthogonal",
                "/carrier/atoms_sum_to_one",
                "/carrier/central_port_atoms",
                "/refinement_tower",
                "/source_readback",
            ],
        },
        "target_firewall": {
            "enabled": True,
            "semantic_input_count": 1,
            "arbitrary_input_path_allowed": False,
            "command_line_target_parameters_allowed": False,
            "desired_coefficient_input_allowed": False,
            "measurement_input_allowed": False,
            "environment_target_read": False,
            "network_target_read": False,
            "nonstdlib_imports": [],
            "audited_import_roots": imported,
            "audited_import_roots_by_file": imported_by_file,
            "exact_manifest_keysets_enforced": True,
            "basis_seed_rule": "lexicographically first unassigned tensor coordinate",
            "target_free_scope": "conditional_on_pinned_canonical_oriented_twelve_port_carrier",
        },
        "proper_port_action": {
            "port_count": PORTS,
            "incidence_automorphism_count": full_order,
            "orientation_preserving_count": len(group),
            "orientation_reversing_count": full_order - len(group),
            "permutation_rows_sha256": group_sha256,
            "permutation_rows": group_rows,
            "identity_present": True,
            "all_rows_bijective": True,
            "incidence_preserved": True,
            "oriented_faces_preserved": True,
            "closure_checked_pairs": len(group) ** 2,
            "inverses_present": True,
            "transitive_on_ports": True,
            "a5_isomorphism": a5_witness,
        },
        "representation": {
            "field": "Q",
            "V_dimension": PORTS,
            "exterior_square_V_dimension": PORTS * (PORTS - 1) // 2,
            "hom_ambient_dimension": len(all_coordinates()),
            "action": "simultaneous output and signed unordered-input permutation",
        },
        "dimension_certificate": character,
        "reynolds_basis": {
            "path": "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_reynolds_basis.json",
            "canonical_json_sha256": basis_sha256,
            "basis_count": len(vectors),
            "orbit_sizes": orbit_sizes,
            "support_union_size": sum(len(vector) for vector in vectors),
            "supports_pairwise_disjoint": True,
            "all_signed_orbits_orientable": True,
            "normalization": "R(seed)=(1/60)*sum_g g(seed); nonzero entries are signed reciprocals of orbit size",
            "alternation_checks": alternation_checks,
            "exact_alternation_passed": True,
            "covariance_checks": covariance_checks,
            "exact_covariance_passed": True,
            "rank_certificate": rank,
        },
        "mutation_tests": mutation_tests,
        "later_gates": {
            "carrier_choice_derived_by_this_packet": False,
            "bracket_selected": False,
            "source_selection": False,
            "jacobi_solution_found": False,
            "jacobi_solution_classification_complete": False,
            "compactness_established": False,
            "physical_current_identified": False,
            "issue_566_closed": False,
        },
        "implementation_pins": {
            "producer_sha256": sha256_file(Path(__file__).resolve()),
            "independent_verifier_sha256": sha256_file(VERIFIER_PATH.resolve()),
            "test_sha256": sha256_file(TEST_PATH.resolve()),
        },
    }
    receipt["receipt_sha256"] = sha256_bytes(canonical_bytes(receipt))
    return basis, receipt


def serialized(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True) + "\n"


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="verify committed artifacts are byte-for-byte current")
    args = parser.parse_args(argv)
    basis, receipt = build_artifacts()
    basis_text = serialized(basis)
    receipt_text = serialized(receipt)
    if args.check:
        require(BASIS_PATH.read_text(encoding="utf-8") == basis_text, "ARTIFACT_STALE", "basis artifact is stale")
        require(RECEIPT_PATH.read_text(encoding="utf-8") == receipt_text, "ARTIFACT_STALE", "receipt artifact is stale")
    else:
        BASIS_PATH.write_text(basis_text, encoding="utf-8")
        RECEIPT_PATH.write_text(receipt_text, encoding="utf-8")
    print(
        json.dumps(
            {
                "status": receipt["status"],
                "proper_action_order": receipt["proper_port_action"]["orientation_preserving_count"],
                "dimension": receipt["dimension_certificate"]["dimension"],
                "basis_rank": receipt["reynolds_basis"]["rank_certificate"]["rank_over_Q"],
                "later_gates_all_false": not any(receipt["later_gates"].values()),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except CertificateError as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(1)
