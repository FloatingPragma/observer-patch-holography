#!/usr/bin/env python3
"""Independent verifier for the issue #566 alternating-space stage.

This file deliberately imports no producer code.  It reconstructs the proper
port action by oriented-triangle propagation (a different algorithm from the
producer's graph-automorphism backtracking), recomputes the character trace and
all Reynolds vectors, and performs its own exact rational row reduction.

The target-free claim is conditional on the pinned canonical oriented
twelve-port carrier.  This verifier also reconstructs every claim-bearing
receipt and basis metadata field, including the statement that the packet does
not derive the carrier, select a bracket, or close the issue.
"""

from __future__ import annotations

import ast
import hashlib
import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Sequence


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
MANIFEST_PATH = REPO / "code/a5_closure/manifests/echosahedral_federation_reference.json"
BASIS_PATH = HERE / "a5_alternating_bracket_reynolds_basis.json"
RECEIPT_PATH = HERE / "a5_alternating_bracket_space_stage1.receipt.json"
PRODUCER_PATH = HERE / "certify.py"
TEST_PATH = HERE / "test_stage1.py"

TOP_KEYS = {"architecture", "carrier", "refinement_tower", "schema", "source_readback"}
CARRIER_KEYS = {
    "atoms_pairwise_orthogonal",
    "atoms_sum_to_one",
    "central_port_atoms",
    "edges",
    "oriented_faces",
    "ports",
}
PRODUCER_ALLOWED_IMPORTS = {
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
VERIFIER_ALLOWED_IMPORTS = {
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
TEST_ALLOWED_IMPORTS = {
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

STATUS = "EXACT_TARGET_FREE_SEARCH_SPACE_STAGE1_ONLY"
CLAIM_BOUNDARY = (
    "Conditional on the pinned canonical oriented twelve-port carrier, this receipt certifies only the complete "
    "equivariant alternating-bracket search space on its permutation module. Target-free means that no desired "
    "gauge algebra, measurement, or coefficient target enters this stage; it does not mean that this packet derives "
    "the carrier choice. The receipt does not select a bracket, solve Jacobi, establish compactness, identify a "
    "physical current, or close issue #566."
)
SEMANTIC_PATH = "code/a5_closure/manifests/echosahedral_federation_reference.json"
BASIS_ARTIFACT_PATH = "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_reynolds_basis.json"
ACCESSED_PATHS = ["/schema", "/carrier/ports", "/carrier/edges", "/carrier/oriented_faces"]
IGNORED_PATHS = [
    "/architecture",
    "/carrier/atoms_pairwise_orthogonal",
    "/carrier/atoms_sum_to_one",
    "/carrier/central_port_atoms",
    "/refinement_tower",
    "/source_readback",
]
LATER_GATES = {
    "bracket_selected": False,
    "carrier_choice_derived_by_this_packet": False,
    "compactness_established": False,
    "issue_566_closed": False,
    "jacobi_solution_classification_complete": False,
    "jacobi_solution_found": False,
    "physical_current_identified": False,
    "source_selection": False,
}
PRODUCER_MUTATIONS = [
    {"name": "remove_carrier_edge", "expected_error": "EDGE_COUNT", "actual_error": "EDGE_COUNT", "passed": True},
    {"name": "reverse_one_oriented_face", "expected_error": "FACE_ORIENTATION", "actual_error": "FACE_ORIENTATION", "passed": True},
    {
        "name": "inject_unregistered_semantic_target",
        "expected_error": "FIREWALL_CARRIER_KEYS",
        "actual_error": "FIREWALL_CARRIER_KEYS",
        "passed": True,
    },
    {
        "name": "break_permutation_bijection",
        "expected_error": "PERMUTATION_BIJECTION",
        "actual_error": "PERMUTATION_BIJECTION",
        "passed": True,
    },
    {"name": "flip_one_basis_coefficient", "expected_error": "COVARIANCE", "actual_error": "COVARIANCE", "passed": True},
    {"name": "drop_one_basis_vector", "expected_error": "BASIS_RANK", "actual_error": "BASIS_RANK", "passed": True},
    {"name": "alter_reynolds_normalization", "expected_error": "COVARIANCE", "actual_error": "COVARIANCE", "passed": True},
]

RECEIPT_KEYS = {
    "claim_boundary",
    "dimension_certificate",
    "implementation_pins",
    "issue",
    "later_gates",
    "mutation_tests",
    "proper_port_action",
    "receipt_sha256",
    "representation",
    "reynolds_basis",
    "schema",
    "semantic_input",
    "status",
    "target_firewall",
}
SEMANTIC_KEYS = {"path", "file_sha256", "carrier_projection_sha256", "accessed_value_paths", "ignored_value_sections"}
FIREWALL_KEYS = {
    "enabled",
    "semantic_input_count",
    "arbitrary_input_path_allowed",
    "command_line_target_parameters_allowed",
    "desired_coefficient_input_allowed",
    "measurement_input_allowed",
    "environment_target_read",
    "network_target_read",
    "nonstdlib_imports",
    "audited_import_roots",
    "audited_import_roots_by_file",
    "exact_manifest_keysets_enforced",
    "basis_seed_rule",
    "target_free_scope",
}
ACTION_KEYS = {
    "port_count",
    "incidence_automorphism_count",
    "orientation_preserving_count",
    "orientation_reversing_count",
    "permutation_rows_sha256",
    "permutation_rows",
    "identity_present",
    "all_rows_bijective",
    "incidence_preserved",
    "oriented_faces_preserved",
    "closure_checked_pairs",
    "inverses_present",
    "transitive_on_ports",
    "a5_isomorphism",
}
A5_KEYS = {
    "model",
    "model_order",
    "bijection_size",
    "port_generators",
    "model_generators",
    "generator_orders",
    "homomorphism_checked_pairs",
    "mapping_sha256",
    "isomorphic",
}
REPRESENTATION_KEYS = {"field", "V_dimension", "exterior_square_V_dimension", "hom_ambient_dimension", "action"}
DIMENSION_KEYS = {
    "formula",
    "class_rows",
    "projector_trace_numerator",
    "projector_trace_denominator",
    "dimension",
}
CLASS_ROW_KEYS = {
    "element_order",
    "count",
    "chi_V",
    "chi_V_of_square",
    "chi_exterior_square_V",
    "projector_trace_contribution_each",
    "projector_trace_contribution_total",
}
REYNOLDS_RECEIPT_KEYS = {
    "path",
    "canonical_json_sha256",
    "basis_count",
    "orbit_sizes",
    "support_union_size",
    "supports_pairwise_disjoint",
    "all_signed_orbits_orientable",
    "normalization",
    "alternation_checks",
    "exact_alternation_passed",
    "covariance_checks",
    "exact_covariance_passed",
    "rank_certificate",
}
RANK_KEYS = {
    "ambient_coordinate_count",
    "matrix_shape",
    "rank_over_Q",
    "pivot_coordinates",
    "pivot_minor_determinant",
    "gram_diagonal",
    "gram_determinant",
}
FRACTION_KEYS = {"numerator", "denominator"}
PIN_KEYS = {"producer_sha256", "independent_verifier_sha256", "test_sha256"}
BASIS_KEYS = {
    "schema",
    "issue",
    "field",
    "port_dimension",
    "domain_dimension",
    "ambient_coordinate_count",
    "coordinate_convention",
    "group_action_convention",
    "normalization",
    "proper_port_group_sha256",
    "basis",
}
BASIS_ROW_KEYS = {"basis_id", "seed_coordinate", "orbit_size", "stabilizer_size", "integral_orbit_scale", "entries"}
MUTATION_ROW_KEYS = {"name", "expected_error", "actual_error", "passed"}

Permutation = tuple[int, ...]
Coordinate = tuple[int, int, int]


class VerificationError(RuntimeError):
    def __init__(self, code: str, message: str):
        super().__init__(f"{code}: {message}")
        self.code = code


def check(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise VerificationError(code, message)


def check_keys(raw: Any, expected: set[str], code: str) -> None:
    check(isinstance(raw, Mapping), code, "value is not an object")
    check(set(raw) == expected, code, f"keyset differs: {sorted(set(raw) ^ expected)}")


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")


def digest(value: Any) -> str:
    return "sha256:" + hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_digest(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def cycle(face: Sequence[int]) -> tuple[int, int, int]:
    a, b, c = face
    return min((a, b, c), (b, c, a), (c, a, b))


def parse_carrier(raw: Mapping[str, Any]) -> dict[str, Any]:
    check(isinstance(raw, Mapping), "FIREWALL_MANIFEST_TYPE", "manifest is not an object")
    check(set(raw) == TOP_KEYS, "FIREWALL_TOP_KEYS", "unexpected top-level field")
    check(raw.get("schema") == "oph.echosahedral_selector_manifest.v1", "FIREWALL_SCHEMA", "wrong schema")
    carrier = raw.get("carrier")
    check(isinstance(carrier, Mapping), "FIREWALL_CARRIER_TYPE", "carrier is not an object")
    check(set(carrier) == CARRIER_KEYS, "FIREWALL_CARRIER_KEYS", "unexpected carrier field")
    ports_raw = carrier.get("ports")
    check(isinstance(ports_raw, list) and len(ports_raw) == 12, "PORT_COUNT", "wrong port count")
    check(all(type(x) is str for x in ports_raw) and len(set(ports_raw)) == 12, "PORT_LABEL", "bad port labels")
    ports = tuple(ports_raw)
    index = {name: i for i, name in enumerate(ports)}

    edge_rows = carrier.get("edges")
    check(isinstance(edge_rows, list) and len(edge_rows) == 30, "EDGE_COUNT", "wrong edge count")
    edges: set[tuple[int, int]] = set()
    for row in edge_rows:
        check(isinstance(row, list) and len(row) == 2, "EDGE_FORMAT", "bad edge row")
        check(all(type(x) is str and x in index for x in row), "EDGE_PORT", "edge has unknown port")
        i, j = index[row[0]], index[row[1]]
        check(i != j, "EDGE_LOOP", "edge loop")
        edges.add((min(i, j), max(i, j)))
    check(len(edges) == 30, "EDGE_DUPLICATE", "duplicate edge")

    face_rows = carrier.get("oriented_faces")
    check(isinstance(face_rows, list) and len(face_rows) == 20, "FACE_COUNT", "wrong face count")
    faces: set[tuple[int, int, int]] = set()
    directed_third: dict[tuple[int, int], int] = {}
    for row in face_rows:
        check(isinstance(row, list) and len(row) == 3, "FACE_FORMAT", "bad face row")
        check(all(type(x) is str and x in index for x in row), "FACE_PORT", "face has unknown port")
        a, b, c = (index[x] for x in row)
        check(len({a, b, c}) == 3, "FACE_REPEAT", "face repeats a port")
        canonical = cycle((a, b, c))
        check(canonical not in faces, "FACE_DUPLICATE", "duplicate face")
        faces.add(canonical)
        for u, v, w in ((a, b, c), (b, c, a), (c, a, b)):
            check((min(u, v), max(u, v)) in edges, "FACE_NONEDGE", "face contains non-edge")
            check((u, v) not in directed_third, "FACE_ORIENTATION", "directed edge occurs twice")
            directed_third[(u, v)] = w
    check(
        all((i, j) in directed_third and (j, i) in directed_third for i, j in edges),
        "FACE_ORIENTATION",
        "orientation does not cover both directions",
    )
    adjacency = [[False] * 12 for _ in range(12)]
    for i, j in edges:
        adjacency[i][j] = adjacency[j][i] = True
    check(all(sum(row) == 5 for row in adjacency), "EDGE_DEGREE", "graph is not five-regular")
    projection = {
        "schema": raw["schema"],
        "ports": list(ports),
        "edges": [list(edge) for edge in sorted(edges)],
        "oriented_faces": [list(face) for face in sorted(faces)],
    }
    return {
        "ports": ports,
        "edges": edges,
        "faces": faces,
        "third": directed_third,
        "adjacency": adjacency,
        "projection": projection,
    }


def propagate_action(carrier: Mapping[str, Any], target_root: tuple[int, int, int], orientation: int) -> Permutation:
    source_root = min(carrier["faces"])
    mapping = dict(zip(source_root, target_root))
    changed = True
    while changed:
        changed = False
        for (u, v), w in carrier["third"].items():
            if u not in mapping or v not in mapping:
                continue
            key = (mapping[u], mapping[v]) if orientation == 1 else (mapping[v], mapping[u])
            check(key in carrier["third"], "PROPAGATION_EDGE", "mapped directed edge is absent")
            image = carrier["third"][key]
            if w in mapping:
                check(mapping[w] == image, "PROPAGATION_CONFLICT", "triangle propagation conflicts")
            else:
                mapping[w] = image
                changed = True
    check(len(mapping) == 12, "PROPAGATION_INCOMPLETE", "triangle propagation did not cover all ports")
    row = tuple(mapping[i] for i in range(12))
    check(set(row) == set(range(12)), "PERMUTATION_BIJECTION", "propagated row is not bijective")
    return row


def reconstruct_actions(carrier: Mapping[str, Any]) -> tuple[list[Permutation], list[Permutation]]:
    proper: list[Permutation] = []
    reversing: list[Permutation] = []
    for face in sorted(carrier["faces"]):
        a, b, c = face
        positive = (a, b, c)
        negative = (a, c, b)
        for root in (positive, positive[1:] + positive[:1], positive[2:] + positive[:2]):
            proper.append(propagate_action(carrier, root, 1))
        for root in (negative, negative[1:] + negative[:1], negative[2:] + negative[:2]):
            reversing.append(propagate_action(carrier, root, -1))
    proper = sorted(set(proper))
    reversing = sorted(set(reversing))
    check(len(proper) == 60, "GROUP_COUNT", "proper action count is not sixty")
    check(len(reversing) == 60, "FULL_GROUP_COUNT", "reversing action count is not sixty")
    check(set(proper).isdisjoint(reversing), "ORIENTATION_SPLIT", "orientation classes overlap")
    return proper, reversing


def compose(p: Permutation, q: Permutation) -> Permutation:
    return tuple(p[q[i]] for i in range(len(p)))


def invert(p: Permutation) -> Permutation:
    row = [0] * len(p)
    for i, image in enumerate(p):
        row[image] = i
    return tuple(row)


def order(p: Permutation) -> int:
    identity = tuple(range(len(p)))
    value = identity
    for n in range(1, 100):
        value = compose(p, value)
        if value == identity:
            return n
    raise VerificationError("GROUP_ORDER", "order bound exceeded")


def closure(generators: Sequence[Permutation]) -> set[Permutation]:
    identity = tuple(range(len(generators[0])))
    moves = list(generators) + [invert(p) for p in generators]
    result = {identity}
    pending = [identity]
    while pending:
        current = pending.pop()
        for move in moves:
            candidate = compose(move, current)
            if candidate not in result:
                result.add(candidate)
                pending.append(candidate)
    return result


def verify_group(group: Sequence[Permutation], carrier: Mapping[str, Any]) -> None:
    check(len(group) == 60 and len(set(group)) == 60, "GROUP_COUNT", "wrong group size")
    group_set = set(group)
    for p in group:
        check(len(p) == 12 and set(p) == set(range(12)), "PERMUTATION_BIJECTION", "bad group row")
        for i in range(12):
            for j in range(12):
                check(
                    carrier["adjacency"][i][j] == carrier["adjacency"][p[i]][p[j]],
                    "GROUP_INCIDENCE",
                    "incidence is not preserved",
                )
        check(invert(p) in group_set, "GROUP_INVERSE", "inverse is absent")
        for q in group:
            check(compose(p, q) in group_set, "GROUP_CLOSURE", "closure fails")


def is_even(p: Permutation) -> bool:
    return sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p))) % 2 == 0


def verify_a5_witness(group: Sequence[Permutation], witness: Mapping[str, Any]) -> None:
    check_keys(witness, A5_KEYS, "A5_WITNESS_KEYS")
    even5 = sorted(tuple(p) for p in itertools.permutations(range(5)) if is_even(tuple(p)))
    check(len(even5) == 60, "A5_MODEL", "even model has wrong size")
    port_generators = tuple(tuple(row) for row in witness.get("port_generators", []))
    model_generators = tuple(tuple(row) for row in witness.get("model_generators", []))
    check(len(port_generators) == 2 and all(p in group for p in port_generators), "A5_GENERATORS", "bad port generators")
    check(len(model_generators) == 2 and all(p in even5 for p in model_generators), "A5_GENERATORS", "bad model generators")
    check(
        [order(port_generators[0]), order(port_generators[1]), order(compose(*port_generators))] == [2, 3, 5],
        "A5_GENERATOR_ORDERS",
        "port generator orders are wrong",
    )
    check(
        [order(model_generators[0]), order(model_generators[1]), order(compose(*model_generators))] == [2, 3, 5],
        "A5_GENERATOR_ORDERS",
        "model generator orders are wrong",
    )
    check(closure(port_generators) == set(group), "A5_GENERATION", "port generators do not generate")
    check(closure(model_generators) == set(even5), "A5_GENERATION", "model generators do not generate")
    moves = [
        (port_generators[0], model_generators[0]),
        (port_generators[1], model_generators[1]),
        (invert(port_generators[0]), invert(model_generators[0])),
        (invert(port_generators[1]), invert(model_generators[1])),
    ]
    p_identity, m_identity = tuple(range(12)), tuple(range(5))
    mapping = {p_identity: m_identity}
    pending = [p_identity]
    while pending:
        current = pending.pop()
        for p_move, m_move in moves:
            p_next = compose(p_move, current)
            m_next = compose(m_move, mapping[current])
            if p_next in mapping:
                check(mapping[p_next] == m_next, "A5_RELATION", "word image conflicts")
            else:
                mapping[p_next] = m_next
                pending.append(p_next)
    check(set(mapping) == set(group) and set(mapping.values()) == set(even5), "A5_BIJECTION", "map is not bijective")
    for p in group:
        for q in group:
            check(mapping[compose(p, q)] == compose(mapping[p], mapping[q]), "A5_HOMOMORPHISM", "map is not a homomorphism")
    serialized = [[list(p), list(mapping[p])] for p in sorted(mapping)]
    check(witness["model"] == "even permutations on five unlabeled symbols", "A5_METADATA", "model label differs")
    check(witness["model_order"] == 60, "A5_METADATA", "model order differs")
    check(witness["bijection_size"] == len(mapping) == 60, "A5_METADATA", "bijection size differs")
    check(witness["generator_orders"] == [2, 3, 5], "A5_METADATA", "generator-order declaration differs")
    check(witness["homomorphism_checked_pairs"] == len(group) ** 2, "A5_METADATA", "pair count differs")
    check(witness["isomorphic"] is True, "A5_METADATA", "isomorphism verdict is not true")
    check(digest(serialized) == witness.get("mapping_sha256"), "A5_MAPPING_HASH", "A5 map digest differs")


def coordinates() -> list[Coordinate]:
    return [(k, i, j) for k in range(12) for i in range(12) for j in range(i + 1, 12)]


def act(p: Permutation, coordinate: Coordinate) -> tuple[Coordinate, int]:
    k, i, j = coordinate
    a, b = p[i], p[j]
    return (p[k], min(a, b), max(a, b)), (1 if a < b else -1)


def character_data(group: Sequence[Permutation]) -> dict[str, Any]:
    bins: dict[tuple[int, int, int], int] = {}
    total = 0
    for p in group:
        fixed = sum(p[i] == i for i in range(12))
        fixed_square = sum(compose(p, p)[i] == i for i in range(12))
        wedge = (fixed * fixed - fixed_square) // 2
        bins[(order(p), fixed, fixed_square)] = bins.get((order(p), fixed, fixed_square), 0) + 1
        total += fixed * wedge
    table = []
    for (element_order, fixed, fixed_square), count in sorted(bins.items()):
        wedge = (fixed * fixed - fixed_square) // 2
        each = fixed * wedge
        table.append(
            {
                "element_order": element_order,
                "count": count,
                "chi_V": fixed,
                "chi_V_of_square": fixed_square,
                "chi_exterior_square_V": wedge,
                "projector_trace_contribution_each": each,
                "projector_trace_contribution_total": count * each,
            }
        )
    check(total == 840 and total // 60 == 14, "CHARACTER_DIMENSION", "character trace is not fourteen")
    return {"rows": table, "numerator": total, "dimension": total // 60}


def reynolds(seed: Coordinate, group: Sequence[Permutation]) -> dict[Coordinate, Fraction]:
    vector: dict[Coordinate, Fraction] = {}
    for p in group:
        image, sign = act(p, seed)
        vector[image] = vector.get(image, Fraction(0)) + Fraction(sign, len(group))
    return {coordinate: value for coordinate, value in vector.items() if value}


def recompute_basis(group: Sequence[Permutation]) -> list[dict[Coordinate, Fraction]]:
    unseen = set(coordinates())
    result: list[dict[Coordinate, Fraction]] = []
    while unseen:
        seed = min(unseen)
        orbit = {act(p, seed)[0] for p in group}
        vector = reynolds(seed, group)
        unseen -= orbit
        if vector:
            check(set(vector) == orbit, "REYNOLDS_SUPPORT", "partial signed orbit survived")
            result.append(dict(sorted(vector.items())))
    check(len(result) == 14, "BASIS_COUNT", "wrong Reynolds basis count")
    return result


def expected_basis_artifact(
    vectors: Sequence[Mapping[Coordinate, Fraction]], group_sha256: str
) -> dict[str, Any]:
    rows = []
    for number, vector in enumerate(vectors):
        seed = min(vector)
        orbit_size = len(vector)
        rows.append(
            {
                "basis_id": f"R{number:02d}",
                "seed_coordinate": list(seed),
                "orbit_size": orbit_size,
                "stabilizer_size": 60 // orbit_size,
                "integral_orbit_scale": orbit_size,
                "entries": [
                    [coordinate[0], coordinate[1], coordinate[2], value.numerator, value.denominator]
                    for coordinate, value in sorted(vector.items())
                ],
            }
        )
    return {
        "schema": "oph.a5_alternating_bracket_reynolds_basis.v1",
        "issue": 566,
        "field": "Q",
        "port_dimension": 12,
        "domain_dimension": 66,
        "ambient_coordinate_count": 792,
        "coordinate_convention": "[output,left,right] with left<right; swapping inputs negates the coefficient",
        "group_action_convention": "(k,i,j) maps to (g(k),sort(g(i),g(j))) with the sorting sign",
        "normalization": "each vector is the exact Reynolds average (1/60)*sum_g g(seed)",
        "proper_port_group_sha256": group_sha256,
        "basis": rows,
    }


def decode_basis(raw: Mapping[str, Any]) -> list[dict[Coordinate, Fraction]]:
    check_keys(raw, BASIS_KEYS, "BASIS_KEYS")
    check(raw.get("schema") == "oph.a5_alternating_bracket_reynolds_basis.v1", "BASIS_SCHEMA", "wrong basis schema")
    rows = raw.get("basis")
    check(isinstance(rows, list), "BASIS_FORMAT", "basis is not a list")
    result = []
    for number, row in enumerate(rows):
        check_keys(row, BASIS_ROW_KEYS, "BASIS_ROW_KEYS")
        check(isinstance(row, Mapping) and row.get("basis_id") == f"R{number:02d}", "BASIS_FORMAT", "bad basis row")
        entries = row.get("entries")
        check(isinstance(entries, list), "BASIS_FORMAT", "entries are not a list")
        vector: dict[Coordinate, Fraction] = {}
        for entry in entries:
            check(
                isinstance(entry, list) and len(entry) == 5 and all(type(value) is int for value in entry),
                "BASIS_ENTRY",
                "basis entry must contain five integers",
            )
            k, i, j, numerator, denominator = entry
            coordinate = (k, i, j)
            check(coordinate in set(coordinates()) and denominator > 0, "BASIS_ENTRY", "bad coordinate or denominator")
            check(coordinate not in vector, "BASIS_ENTRY", "duplicate coordinate")
            vector[coordinate] = Fraction(numerator, denominator)
        check(list(row.get("seed_coordinate", [])) == list(min(vector)), "BASIS_SEED", "seed differs")
        check(row.get("orbit_size") == len(vector), "BASIS_ORBIT", "orbit size differs")
        result.append(vector)
    return result


def verify_covariance(vectors: Sequence[Mapping[Coordinate, Fraction]], group: Sequence[Permutation]) -> int:
    count = 0
    for vector in vectors:
        for p in group:
            for coordinate in coordinates():
                image, sign = act(p, coordinate)
                check(vector.get(image, Fraction(0)) == sign * vector.get(coordinate, Fraction(0)), "COVARIANCE", "coefficient covariance fails")
                count += 1
    return count


def verify_alternation(vectors: Sequence[Mapping[Coordinate, Fraction]]) -> int:
    count = 0
    for vector in vectors:
        for k in range(12):
            for i in range(12):
                for j in range(12):
                    def coefficient(left: int, right: int) -> Fraction:
                        if left == right:
                            return Fraction(0)
                        if left < right:
                            return vector.get((k, left, right), Fraction(0))
                        return -vector.get((k, right, left), Fraction(0))

                    check(coefficient(i, j) == -coefficient(j, i), "ALTERNATION", "antisymmetry fails")
                    count += 1
    return count


def exact_rank(vectors: Sequence[Mapping[Coordinate, Fraction]]) -> int:
    matrix = [[vector.get(coordinate, Fraction(0)) for vector in vectors] for coordinate in coordinates()]
    row = 0
    for column in range(len(vectors)):
        pivot = next((r for r in range(row, len(matrix)) if matrix[r][column]), None)
        if pivot is None:
            continue
        matrix[row], matrix[pivot] = matrix[pivot], matrix[row]
        pivot_value = matrix[row][column]
        matrix[row] = [value / pivot_value for value in matrix[row]]
        for r in range(row + 1, len(matrix)):
            if matrix[r][column]:
                factor = matrix[r][column]
                matrix[r] = [x - factor * y for x, y in zip(matrix[r], matrix[row])]
        row += 1
    return row


def fraction_from_json(raw: Mapping[str, Any]) -> Fraction:
    check_keys(raw, FRACTION_KEYS, "FRACTION_KEYS")
    check(type(raw["numerator"]) is int and type(raw["denominator"]) is int, "FRACTION_FORMAT", "fraction is not integral")
    check(raw["denominator"] > 0, "FRACTION_FORMAT", "fraction denominator is not positive")
    return Fraction(raw["numerator"], raw["denominator"])


def verify_rank_receipt(vectors: Sequence[Mapping[Coordinate, Fraction]], raw: Mapping[str, Any]) -> None:
    check_keys(raw, RANK_KEYS, "RANK_KEYS")
    rank = exact_rank(vectors)
    check(rank == 14, "BASIS_RANK", "exact row reduction did not give fourteen")
    check(raw["ambient_coordinate_count"] == 792, "RANK_RECEIPT", "ambient coordinate count differs")
    check(raw["matrix_shape"] == [792, 14], "RANK_RECEIPT", "matrix shape differs")
    check(raw["rank_over_Q"] == rank, "RANK_RECEIPT", "declared rank differs")
    seeds = [min(vector) for vector in vectors]
    check(raw.get("pivot_coordinates") == [list(seed) for seed in seeds], "RANK_RECEIPT", "pivot rows differ")
    pivot_det = Fraction(1)
    for number, seed in enumerate(seeds):
        check(all(vectors[j].get(seed, 0) == 0 for j in range(len(vectors)) if j != number), "RANK_RECEIPT", "seed minor is not diagonal")
        pivot_det *= vectors[number][seed]
    check(pivot_det == fraction_from_json(raw["pivot_minor_determinant"]), "RANK_RECEIPT", "minor determinant differs")
    diagonal = []
    for vector in vectors:
        diagonal.append(sum((value * value for value in vector.values()), Fraction(0)))
    check(diagonal == [fraction_from_json(row) for row in raw["gram_diagonal"]], "RANK_RECEIPT", "Gram diagonal differs")
    gram_det = Fraction(1)
    for value in diagonal:
        gram_det *= value
    check(gram_det == fraction_from_json(raw["gram_determinant"]), "RANK_RECEIPT", "Gram determinant differs")


def audit_imports() -> dict[str, list[str]]:
    paths = (PRODUCER_PATH.resolve(), Path(__file__).resolve(), TEST_PATH.resolve())
    allowlists = {
        PRODUCER_PATH.resolve(): PRODUCER_ALLOWED_IMPORTS,
        Path(__file__).resolve(): VERIFIER_ALLOWED_IMPORTS,
        TEST_PATH.resolve(): TEST_ALLOWED_IMPORTS,
    }
    result: dict[str, list[str]] = {}
    for path in paths:
        imports: set[str] = set()
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imports.add(node.module.split(".")[0])
        allowed = allowlists[path]
        check(imports <= allowed, "FIREWALL_IMPORT", f"non-allowlisted import in {path.name}: {sorted(imports - allowed)}")
        result[path.name] = sorted(imports)
    return dict(sorted(result.items()))


def expect_error(expected: str, action: Any) -> bool:
    try:
        action()
    except VerificationError as exc:
        return exc.code == expected
    return False


def independent_mutations(
    manifest: Mapping[str, Any], carrier: Mapping[str, Any], group: Sequence[Permutation], vectors: Sequence[Mapping[Coordinate, Fraction]]
) -> dict[str, bool]:
    injected = json.loads(json.dumps(manifest))
    injected["carrier"]["target_payload"] = True
    target_caught = expect_error("FIREWALL_CARRIER_KEYS", lambda: parse_carrier(injected))

    bad_group = list(group)
    bad_group[0] = tuple([bad_group[0][0]] * 12)
    group_caught = expect_error("PERMUTATION_BIJECTION", lambda: verify_group(bad_group, carrier))

    bad_vectors = [dict(vector) for vector in vectors]
    seed = min(bad_vectors[0])
    bad_vectors[0][seed] = -bad_vectors[0][seed]
    coefficient_caught = expect_error("COVARIANCE", lambda: verify_covariance(bad_vectors, group))

    rank_caught = exact_rank(vectors[:-1]) == 13
    result = {
        "target_injection_caught": target_caught,
        "group_row_tamper_caught": group_caught,
        "basis_coefficient_tamper_caught": coefficient_caught,
        "basis_vector_drop_caught": rank_caught,
    }
    check(all(result.values()), "MUTATION_TEST", "independent mutation control failed")
    return result


def verify() -> dict[str, Any]:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    basis_raw = json.loads(BASIS_PATH.read_text(encoding="utf-8"))
    receipt = json.loads(RECEIPT_PATH.read_text(encoding="utf-8"))
    check(isinstance(manifest, dict) and isinstance(basis_raw, dict) and isinstance(receipt, dict), "JSON_TYPE", "artifact root must be an object")
    check_keys(receipt, RECEIPT_KEYS, "RECEIPT_KEYS")

    stored_receipt_hash = receipt.get("receipt_sha256")
    unhashed = dict(receipt)
    unhashed.pop("receipt_sha256", None)
    check(stored_receipt_hash == digest(unhashed), "RECEIPT_HASH", "receipt self-hash differs")
    check(receipt.get("schema") == "oph.a5_alternating_bracket_space_stage1.receipt.v1", "RECEIPT_SCHEMA", "wrong receipt schema")
    check(receipt.get("issue") == 566, "RECEIPT_ISSUE", "wrong issue")
    check(receipt.get("status") == STATUS, "RECEIPT_STATUS", "stage-one status differs")
    check(receipt.get("claim_boundary") == CLAIM_BOUNDARY, "CLAIM_BOUNDARY", "claim boundary differs")

    carrier = parse_carrier(manifest)
    proper, reversing = reconstruct_actions(carrier)
    verify_group(proper, carrier)
    group_hash = digest([list(row) for row in proper])
    action_receipt = receipt["proper_port_action"]
    check_keys(action_receipt, ACTION_KEYS, "ACTION_KEYS")
    check(len(proper) + len(reversing) == action_receipt["incidence_automorphism_count"] == 120, "GROUP_RECEIPT", "full action count differs")
    check(action_receipt["port_count"] == 12, "GROUP_RECEIPT", "port count differs")
    check(action_receipt["orientation_preserving_count"] == len(proper) == 60, "GROUP_RECEIPT", "proper count differs")
    check(action_receipt["orientation_reversing_count"] == len(reversing) == 60, "GROUP_RECEIPT", "reversing count differs")
    check(group_hash == action_receipt["permutation_rows_sha256"], "GROUP_HASH", "proper action hash differs")
    check(action_receipt["permutation_rows"] == [list(row) for row in proper], "GROUP_ROWS", "serialized proper rows differ")
    expected_action_flags = {
        "identity_present": True,
        "all_rows_bijective": True,
        "incidence_preserved": True,
        "oriented_faces_preserved": True,
        "closure_checked_pairs": len(proper) ** 2,
        "inverses_present": True,
        "transitive_on_ports": True,
    }
    for key, expected_value in expected_action_flags.items():
        check(action_receipt[key] == expected_value, "GROUP_RECEIPT", f"{key} differs")
    verify_a5_witness(proper, action_receipt["a5_isomorphism"])

    semantic = receipt["semantic_input"]
    check_keys(semantic, SEMANTIC_KEYS, "SEMANTIC_KEYS")
    check(semantic["path"] == SEMANTIC_PATH, "INPUT_PATH", "semantic input path differs")
    check(file_digest(MANIFEST_PATH) == semantic["file_sha256"], "INPUT_HASH", "manifest hash differs")
    check(digest(carrier["projection"]) == semantic["carrier_projection_sha256"], "INPUT_HASH", "carrier projection hash differs")
    check(semantic["accessed_value_paths"] == ACCESSED_PATHS, "INPUT_SCOPE", "accessed value paths differ")
    check(semantic["ignored_value_sections"] == IGNORED_PATHS, "INPUT_SCOPE", "ignored value sections differ")
    imports_by_file = audit_imports()
    imports = sorted({root for roots in imports_by_file.values() for root in roots})
    firewall = receipt["target_firewall"]
    check_keys(firewall, FIREWALL_KEYS, "FIREWALL_KEYS")
    expected_firewall = {
        "enabled": True,
        "semantic_input_count": 1,
        "arbitrary_input_path_allowed": False,
        "command_line_target_parameters_allowed": False,
        "desired_coefficient_input_allowed": False,
        "measurement_input_allowed": False,
        "environment_target_read": False,
        "network_target_read": False,
        "nonstdlib_imports": [],
        "audited_import_roots": imports,
        "audited_import_roots_by_file": imports_by_file,
        "exact_manifest_keysets_enforced": True,
        "basis_seed_rule": "lexicographically first unassigned tensor coordinate",
        "target_free_scope": "conditional_on_pinned_canonical_oriented_twelve_port_carrier",
    }
    check(firewall == expected_firewall, "FIREWALL_RECEIPT", "firewall metadata differs")

    representation = receipt["representation"]
    check_keys(representation, REPRESENTATION_KEYS, "REPRESENTATION_KEYS")
    check(
        representation
        == {
            "field": "Q",
            "V_dimension": 12,
            "exterior_square_V_dimension": 66,
            "hom_ambient_dimension": 792,
            "action": "simultaneous output and signed unordered-input permutation",
        },
        "REPRESENTATION_RECEIPT",
        "representation metadata differs",
    )

    character = character_data(proper)
    dimension = receipt["dimension_certificate"]
    check_keys(dimension, DIMENSION_KEYS, "DIMENSION_KEYS")
    check(
        dimension["formula"] == "sum_g chi_V(g)*(chi_V(g)^2-chi_V(g^2))/2 divided by group_order",
        "CHARACTER_RECEIPT",
        "character formula differs",
    )
    check(all(isinstance(row, Mapping) and set(row) == CLASS_ROW_KEYS for row in dimension["class_rows"]), "CHARACTER_KEYS", "character row keyset differs")
    check(dimension["class_rows"] == character["rows"], "CHARACTER_RECEIPT", "character table differs")
    check(
        dimension["projector_trace_numerator"] == character["numerator"]
        and dimension["projector_trace_denominator"] == 60
        and dimension["dimension"] == character["dimension"] == 14,
        "CHARACTER_RECEIPT",
        "dimension receipt differs",
    )

    expected = recompute_basis(proper)
    check(basis_raw == expected_basis_artifact(expected, group_hash), "BASIS_ARTIFACT", "basis metadata or coefficients differ")
    check(digest(basis_raw) == receipt["reynolds_basis"]["canonical_json_sha256"], "BASIS_HASH", "basis hash differs")
    vectors = decode_basis(basis_raw)
    check(vectors == expected, "REYNOLDS_BASIS", "serialized basis differs from independent Reynolds averages")
    reynolds_receipt = receipt["reynolds_basis"]
    check_keys(reynolds_receipt, REYNOLDS_RECEIPT_KEYS, "REYNOLDS_RECEIPT_KEYS")
    supports = [set(vector) for vector in vectors]
    support_union = set().union(*supports)
    supports_pairwise_disjoint = sum(map(len, supports)) == len(support_union)
    check(reynolds_receipt["path"] == BASIS_ARTIFACT_PATH, "REYNOLDS_RECEIPT", "basis path differs")
    check(reynolds_receipt["basis_count"] == len(vectors) == 14, "REYNOLDS_RECEIPT", "basis count differs")
    check(reynolds_receipt["orbit_sizes"] == sorted(map(len, vectors)), "REYNOLDS_RECEIPT", "orbit sizes differ")
    check(reynolds_receipt["support_union_size"] == len(support_union) == 792, "REYNOLDS_RECEIPT", "support union differs")
    check(reynolds_receipt["supports_pairwise_disjoint"] is supports_pairwise_disjoint is True, "REYNOLDS_RECEIPT", "support partition differs")
    check(reynolds_receipt["all_signed_orbits_orientable"] is True, "REYNOLDS_RECEIPT", "signed-orbit verdict differs")
    check(
        reynolds_receipt["normalization"]
        == "R(seed)=(1/60)*sum_g g(seed); nonzero entries are signed reciprocals of orbit size",
        "REYNOLDS_RECEIPT",
        "normalization declaration differs",
    )
    alternation_checks = verify_alternation(vectors)
    check(
        reynolds_receipt["exact_alternation_passed"] is True
        and alternation_checks == reynolds_receipt["alternation_checks"],
        "ALTERNATION_RECEIPT",
        "alternation receipt differs",
    )
    covariance_checks = verify_covariance(vectors, proper)
    check(
        reynolds_receipt["exact_covariance_passed"] is True
        and covariance_checks == reynolds_receipt["covariance_checks"],
        "COVARIANCE_RECEIPT",
        "covariance receipt differs",
    )
    verify_rank_receipt(vectors, reynolds_receipt["rank_certificate"])

    check(
        isinstance(receipt["mutation_tests"], list)
        and all(isinstance(row, Mapping) and set(row) == MUTATION_ROW_KEYS for row in receipt["mutation_tests"]),
        "MUTATION_KEYS",
        "producer mutation keyset differs",
    )
    check(receipt["mutation_tests"] == PRODUCER_MUTATIONS, "MUTATION_RECEIPT", "producer mutation receipt differs")
    mutations = independent_mutations(manifest, carrier, proper, vectors)
    later = receipt.get("later_gates")
    check_keys(later, set(LATER_GATES), "LATER_GATE_KEYS")
    check(later == LATER_GATES, "LATER_GATE", "a later gate or boundary differs")
    pins = receipt["implementation_pins"]
    check_keys(pins, PIN_KEYS, "IMPLEMENTATION_PIN_KEYS")
    check(pins["producer_sha256"] == file_digest(PRODUCER_PATH), "IMPLEMENTATION_PIN", "producer pin differs")
    check(pins["independent_verifier_sha256"] == file_digest(Path(__file__).resolve()), "IMPLEMENTATION_PIN", "verifier pin differs")
    check(pins["test_sha256"] == file_digest(TEST_PATH), "IMPLEMENTATION_PIN", "test pin differs")
    return {
        "verified": True,
        "proper_action_order": len(proper),
        "a5_isomorphism_verified": True,
        "dimension": character["dimension"],
        "basis_rank": exact_rank(vectors),
        "exact_alternation_checks": alternation_checks,
        "exact_covariance_checks": covariance_checks,
        "independent_mutations": mutations,
        "later_gates_all_false": not any(later.values()),
    }


if __name__ == "__main__":
    try:
        print(json.dumps(verify(), sort_keys=True))
    except VerificationError as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(1)
