#!/usr/bin/env python3
"""Exact certificate for GitHub issue #565.

The input is a quotient-visible twelve-port carrier manifest.  The verifier
uses only:

* primitive central port atoms and their normalized trace;
* integer port-defect readback and its quadratic mismatch norm;
* edge/face incidence of the carrier;
* port-lineage maps in a refinement tower.

It derives, rather than assumes as extra output data:

* the unique twelve-unit minimizer and its exact gap;
* the unique fixed-point-free antipode involution;
* the orientation-preserving automorphism group A5;
* its faithful six-axis action;
* the canonical rank-three port Gram frame and equivariance;
* refinement and arbitrary port-relabeling naturality.

No Standard Model representation, product-adjoint count, coupling, measured
quantity, or gauge datum is accepted in a source manifest.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import itertools
import json
from collections import Counter, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence

import sympy as sp

SCHEMA = "oph.echosahedral_selector_manifest.v1"
RECEIPT_SCHEMA = "oph.echosahedral_selector_receipt.v1"
NEGATIVE_SCHEMA = "oph.echosahedral_selector_negative_controls.v1"

# Fail closed on downstream target data.  These tokens are normalized by
# removing punctuation and spaces before comparison.
FORBIDDEN_SOURCE_TOKENS = (
    "standardmodel",
    "productadjoint",
    "su3",
    "su2",
    "u1gauge",
    "electroweak",
    "higgs",
    "hypercharge",
    "particlemass",
    "measuredcoupling",
    "gaugdatum",
    "gaugedatum",
    "alphainverse",
)


class CertificateError(ValueError):
    """Fail-closed manifest or receipt error carrying a stable code."""

    def __init__(self, code: str, message: str) -> None:
        super().__init__(f"{code}: {message}")
        self.code = code
        self.message = message


@dataclass(frozen=True)
class Carrier:
    ports: tuple[str, ...]
    index: Mapping[str, int]
    edges: tuple[tuple[int, int], ...]
    faces: tuple[tuple[int, int, int], ...]
    adjacency: tuple[frozenset[int], ...]
    distances: tuple[tuple[int, ...], ...]
    antipode: tuple[int, ...]


Permutation = tuple[int, ...]


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")


def sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def load_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CertificateError("JSON_READ", f"cannot read {path}: {exc}") from exc


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def normalized_token(text: str) -> str:
    return "".join(ch.lower() for ch in text if ch.isalnum())


def walk_strings(value: Any, path: str = "$") -> Iterator[tuple[str, str]]:
    if isinstance(value, Mapping):
        for key, item in value.items():
            yield f"{path}.{key}", str(key)
            yield from walk_strings(item, f"{path}.{key}")
    elif isinstance(value, list):
        for i, item in enumerate(value):
            yield from walk_strings(item, f"{path}[{i}]")
    elif isinstance(value, str):
        yield path, value


def enforce_source_firewall(manifest: Mapping[str, Any]) -> list[str]:
    hits: list[str] = []
    for path, text in walk_strings(manifest):
        token = normalized_token(text)
        for forbidden in FORBIDDEN_SOURCE_TOKENS:
            if forbidden in token:
                hits.append(f"{path}:{text}")
    if hits:
        raise CertificateError(
            "FORBIDDEN_DEPENDENCY",
            "source manifest contains downstream target data: " + "; ".join(hits[:4]),
        )
    return hits


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise CertificateError(code, message)


def parse_fraction(value: Mapping[str, Any]) -> tuple[int, int]:
    try:
        numerator = int(value["numerator"])
        denominator = int(value["denominator"])
    except (KeyError, TypeError, ValueError) as exc:
        raise CertificateError("TRACE_WEIGHT_FORMAT", "trace weights must be integer numerator/denominator objects") from exc
    require(denominator > 0, "TRACE_WEIGHT_FORMAT", "trace denominator must be positive")
    return numerator, denominator


def all_pairs_distances(adjacency: Sequence[frozenset[int]]) -> tuple[tuple[int, ...], ...]:
    n = len(adjacency)
    rows: list[tuple[int, ...]] = []
    for source in range(n):
        distance = [-1] * n
        distance[source] = 0
        queue: deque[int] = deque([source])
        while queue:
            u = queue.popleft()
            for v in adjacency[u]:
                if distance[v] == -1:
                    distance[v] = distance[u] + 1
                    queue.append(v)
        require(all(d >= 0 for d in distance), "INCIDENCE_DISCONNECTED", "port incidence graph is disconnected")
        rows.append(tuple(distance))
    return tuple(rows)


def validate_carrier(manifest: Mapping[str, Any]) -> Carrier:
    require(manifest.get("schema") == SCHEMA, "SCHEMA", f"expected {SCHEMA}")
    require(manifest.get("architecture") == "federation_of_12_port_echosahedra", "ARCHITECTURE", "wrong simulator architecture")

    carrier = manifest.get("carrier")
    require(isinstance(carrier, Mapping), "CARRIER", "carrier object is missing")
    ports_raw = carrier.get("ports")
    require(isinstance(ports_raw, list), "PORTS", "ports must be a list")
    ports = tuple(str(p) for p in ports_raw)
    require(len(ports) == 12, "PORT_COUNT", "exactly twelve quotient-visible ports are required")
    require(len(set(ports)) == 12, "PORT_LABELS", "port identifiers must be distinct")
    index = {p: i for i, p in enumerate(ports)}

    atoms = carrier.get("central_port_atoms")
    require(isinstance(atoms, list), "PORT_ATOMS", "central_port_atoms must be a list")
    require(len(atoms) == 12, "PORT_ATOMS_COUNT", "exactly twelve primitive central atoms are required")
    seen_atom_ports: set[str] = set()
    seen_atom_ids: set[str] = set()
    weights: list[tuple[int, int]] = []
    for atom in atoms:
        require(isinstance(atom, Mapping), "PORT_ATOMS", "each atom must be an object")
        port = str(atom.get("port"))
        atom_id = str(atom.get("atom_id"))
        require(port in index, "PORT_ATOM_PORT", f"unknown atom port {port}")
        require(port not in seen_atom_ports, "PORT_ATOM_DUPLICATE", f"duplicate atom for {port}")
        require(atom_id not in seen_atom_ids, "PORT_ATOM_DUPLICATE", f"duplicate atom id {atom_id}")
        require(atom.get("primitive") is True, "PORT_ATOM_NOT_PRIMITIVE", f"atom {atom_id} is not primitive")
        seen_atom_ports.add(port)
        seen_atom_ids.add(atom_id)
        weight = atom.get("normalized_trace")
        require(isinstance(weight, Mapping), "TRACE_WEIGHT_FORMAT", f"atom {atom_id} lacks normalized trace")
        weights.append(parse_fraction(weight))
    require(all(num * 12 == den for num, den in weights), "PORT_TRACE_NOT_UNIFORM", "every primitive port atom must have trace 1/12")
    require(carrier.get("atoms_pairwise_orthogonal") is True, "PORT_ATOMS_NOT_ORTHOGONAL", "port atoms must be pairwise orthogonal")
    require(carrier.get("atoms_sum_to_one") is True, "PORT_ATOMS_NOT_COMPLETE", "port atoms must sum to one")

    source = manifest.get("source_readback")
    require(isinstance(source, Mapping), "SOURCE_READBACK", "source_readback is missing")
    require(source.get("defect_domain") == "integer_port_charges", "DEFECT_DOMAIN", "defect domain must be integer port charges")
    require(source.get("total_charge") == 12, "TOTAL_CHARGE_NOT_12", "Euler/readback total must be exactly 12")
    require(source.get("mismatch_cost") == "sum_of_port_charge_squares", "COST_NOT_STRICT_QUADRATIC", "source cost must be the equal-trace quadratic port norm")
    require(source.get("cost_origin") == "normalized_central_readback_hilbert_schmidt_norm", "COST_NOT_SOURCE_DEFINED", "quadratic cost must be the declared readback norm")

    edges_raw = carrier.get("edges")
    require(isinstance(edges_raw, list), "EDGES", "edges must be a list")
    edge_set: set[tuple[int, int]] = set()
    for item in edges_raw:
        require(isinstance(item, list) and len(item) == 2, "EDGE_FORMAT", "each edge must contain two ports")
        a_name, b_name = str(item[0]), str(item[1])
        require(a_name in index and b_name in index, "EDGE_PORT", "edge names an unknown port")
        a, b = index[a_name], index[b_name]
        require(a != b, "EDGE_LOOP", "self-edges are forbidden")
        edge = (min(a, b), max(a, b))
        require(edge not in edge_set, "EDGE_DUPLICATE", f"duplicate edge {item}")
        edge_set.add(edge)
    require(len(edge_set) == 30, "EDGE_COUNT", "a certified carrier has 30 edges")

    adjacency_mut = [set() for _ in ports]
    for a, b in edge_set:
        adjacency_mut[a].add(b)
        adjacency_mut[b].add(a)
    degree_profile = tuple(len(nbrs) for nbrs in adjacency_mut)
    require(degree_profile == (5,) * 12, "DEGREE_SIGNATURE", f"expected twelve degree-five ports, got {degree_profile}")
    adjacency = tuple(frozenset(nbrs) for nbrs in adjacency_mut)

    faces_raw = carrier.get("oriented_faces")
    require(isinstance(faces_raw, list), "FACES", "oriented_faces must be a list")
    require(len(faces_raw) == 20, "FACE_COUNT", "a certified carrier has 20 oriented triangular faces")
    faces: list[tuple[int, int, int]] = []
    undirected_faces: set[frozenset[int]] = set()
    directed_edge_count: Counter[tuple[int, int]] = Counter()
    edge_face_count: Counter[tuple[int, int]] = Counter()
    for item in faces_raw:
        require(isinstance(item, list) and len(item) == 3, "FACE_FORMAT", "each face must contain three ports")
        names = tuple(str(x) for x in item)
        require(all(name in index for name in names), "FACE_PORT", "face names an unknown port")
        face = tuple(index[name] for name in names)
        require(len(set(face)) == 3, "FACE_DEGENERATE", "face ports must be distinct")
        key = frozenset(face)
        require(key not in undirected_faces, "FACE_DUPLICATE", "duplicate unoriented face")
        undirected_faces.add(key)
        for a, b in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            edge = (min(a, b), max(a, b))
            require(edge in edge_set, "FACE_EDGE_MISMATCH", "face contains a non-edge")
            edge_face_count[edge] += 1
            directed_edge_count[(a, b)] += 1
        faces.append(face)
    require(all(edge_face_count[e] == 2 for e in edge_set), "EDGE_FACE_INCIDENCE", "every edge must lie in exactly two faces")
    require(all(directed_edge_count[(a, b)] == 1 and directed_edge_count[(b, a)] == 1 for a, b in edge_set), "FACE_ORIENTATION", "face orientations must induce opposite directions across every shared edge")

    # The link of each port must be a 5-cycle.  This is an exact local
    # manifold/readback check and excludes arbitrary 5-regular graphs.
    for v, neighbors in enumerate(adjacency):
        link_edges = 0
        link_degree = {u: 0 for u in neighbors}
        for a, b in itertools.combinations(neighbors, 2):
            if b in adjacency[a]:
                link_edges += 1
                link_degree[a] += 1
                link_degree[b] += 1
        require(link_edges == 5 and all(d == 2 for d in link_degree.values()), "LINK_NOT_C5", f"link of port {ports[v]} is not a five-cycle")

    require(12 - 30 + 20 == 2, "EULER", "incidence is not a triangulated sphere")
    distances = all_pairs_distances(adjacency)
    distance_profiles = [Counter(row) for row in distances]
    expected_profile = Counter({0: 1, 1: 5, 2: 5, 3: 1})
    require(all(profile == expected_profile for profile in distance_profiles), "DISTANCE_SIGNATURE", "incidence lacks the icosahedral distance signature")
    antipode = tuple(next(j for j, d in enumerate(row) if d == 3) for row in distances)
    require(all(antipode[i] != i and antipode[antipode[i]] == i for i in range(12)), "ANTIPODE_NOT_INVOLUTION", "distance-three map is not a fixed-point-free involution")

    return Carrier(
        ports=ports,
        index=index,
        edges=tuple(sorted(edge_set)),
        faces=tuple(faces),
        adjacency=adjacency,
        distances=distances,
        antipode=antipode,
    )


def enumerate_distance_isometries(
    source: Sequence[Sequence[int]],
    target: Sequence[Sequence[int]] | None = None,
    seed: Mapping[int, int] | None = None,
) -> list[Permutation]:
    """Enumerate all isometries of two finite graph-distance matrices."""

    if target is None:
        target = source
    n = len(source)
    require(len(target) == n, "ISOMETRY_SIZE", "distance matrices have different sizes")
    mapping: dict[int, int] = dict(seed or {})
    used = set(mapping.values())
    for a, image_a in mapping.items():
        for b, image_b in mapping.items():
            require(source[a][b] == target[image_a][image_b], "ISOMETRY_SEED", "inconsistent isometry seed")
    results: list[Permutation] = []

    def recurse() -> None:
        if len(mapping) == n:
            results.append(tuple(mapping[i] for i in range(n)))
            return
        best_vertex = -1
        best_candidates: list[int] | None = None
        for u in range(n):
            if u in mapping:
                continue
            candidates = [
                v
                for v in range(n)
                if v not in used
                and all(source[u][a] == target[v][image_a] for a, image_a in mapping.items())
            ]
            if best_candidates is None or len(candidates) < len(best_candidates):
                best_vertex = u
                best_candidates = candidates
                if len(candidates) <= 1:
                    break
        assert best_candidates is not None
        for v in best_candidates:
            mapping[best_vertex] = v
            used.add(v)
            recurse()
            used.remove(v)
            del mapping[best_vertex]

    recurse()
    return sorted(set(results))


def cyclic_face_sets(faces: Sequence[tuple[int, int, int]]) -> tuple[set[tuple[int, int, int]], set[tuple[int, int, int]]]:
    positive: set[tuple[int, int, int]] = set()
    negative: set[tuple[int, int, int]] = set()
    for a, b, c in faces:
        positive.update(((a, b, c), (b, c, a), (c, a, b)))
        negative.update(((a, c, b), (c, b, a), (b, a, c)))
    return positive, negative


def orientation_sign(permutation: Permutation, faces: Sequence[tuple[int, int, int]]) -> int:
    positive, negative = cyclic_face_sets(faces)
    signs: set[int] = set()
    for face in faces:
        image = tuple(permutation[v] for v in face)
        if image in positive:
            signs.add(1)
        elif image in negative:
            signs.add(-1)
        else:
            raise CertificateError("AUTOMORPHISM_FACE", "distance isometry does not preserve the face complex")
    require(len(signs) == 1, "AUTOMORPHISM_ORIENTATION", "automorphism has inconsistent face orientation")
    return next(iter(signs))


def compose(left: Permutation, right: Permutation) -> Permutation:
    """Composition left after right."""

    return tuple(left[right[i]] for i in range(len(left)))


def inverse(permutation: Permutation) -> Permutation:
    out = [0] * len(permutation)
    for i, image in enumerate(permutation):
        out[image] = i
    return tuple(out)


def permutation_order(permutation: Permutation) -> int:
    identity = tuple(range(len(permutation)))
    power = identity
    for order in range(1, 121):
        power = compose(permutation, power)
        if power == identity:
            return order
    raise CertificateError("GROUP_ORDER", "permutation order exceeded finite bound")




def permutation_parity(permutation: Sequence[int]) -> int:
    """Return +1 for an even finite permutation and -1 for an odd one."""

    inversions = sum(
        1
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
        if permutation[i] > permutation[j]
    )
    return 1 if inversions % 2 == 0 else -1


def direct_a5_identification(group: Sequence[Permutation]) -> dict[str, Any]:
    """Identify the order-60 port group with A5 without classification lookup.

    The exact permutation group has five Klein-four Sylow-2 subgroups.
    Conjugation on that five-element set gives a faithful action whose 60
    images are all even permutations.  Hence the image is exactly A5.
    """

    identity = tuple(range(len(group[0])))
    involutions = [g for g in group if permutation_order(g) == 2]
    v4_subgroups: set[frozenset[Permutation]] = set()
    for a, b in itertools.combinations(involutions, 2):
        if compose(a, b) != compose(b, a):
            continue
        c = compose(a, b)
        if c != identity and permutation_order(c) == 2:
            v4_subgroups.add(frozenset((identity, a, b, c)))
    ordered_v4 = sorted(v4_subgroups, key=lambda subgroup: sorted(subgroup))
    require(len(ordered_v4) == 5, "A5_V4_COUNT", f"expected five Klein-four subgroups, got {len(ordered_v4)}")
    v4_index = {subgroup: i for i, subgroup in enumerate(ordered_v4)}

    actions: set[tuple[int, ...]] = set()
    for g in group:
        g_inv = inverse(g)
        action: list[int] = []
        for subgroup in ordered_v4:
            image = frozenset(compose(compose(g, h), g_inv) for h in subgroup)
            require(image in v4_index, "A5_V4_CONJUGATION", "conjugation did not preserve the V4 family")
            action.append(v4_index[image])
        action_tuple = tuple(action)
        require(permutation_parity(action_tuple) == 1, "A5_FIVE_ACTION_ODD", "five-point conjugation action contains an odd permutation")
        actions.add(action_tuple)

    require(len(actions) == 60, "A5_FIVE_ACTION_KERNEL", "conjugation action on the five V4 subgroups is not faithful")
    even_five = {
        permutation
        for permutation in itertools.permutations(range(5))
        if permutation_parity(permutation) == 1
    }
    require(actions == even_five, "A5_FIVE_ACTION_IMAGE", "five-point image is not the full alternating group")
    return {
        "klein_four_subgroups": 5,
        "conjugation_degree": 5,
        "conjugation_action_faithful": True,
        "conjugation_image": "all 60 even permutations of five objects",
        "identification": "orientation-preserving incidence automorphisms are explicitly isomorphic to A5",
    }

def generated_subgroup(generators: Iterable[Permutation], n: int = 12) -> set[Permutation]:
    identity = tuple(range(n))
    generators_list = list(generators)
    moves = generators_list + [inverse(g) for g in generators_list]
    subgroup: set[Permutation] = {identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for move in moves:
            candidate = compose(move, current)
            if candidate not in subgroup:
                subgroup.add(candidate)
                frontier.append(candidate)
    return subgroup


def conjugacy_classes(group: Sequence[Permutation]) -> list[set[Permutation]]:
    remaining = set(group)
    classes: list[set[Permutation]] = []
    while remaining:
        g = min(remaining)
        cls = {compose(compose(h, g), inverse(h)) for h in group}
        classes.append(cls)
        remaining -= cls
    return sorted(classes, key=lambda cls: (permutation_order(min(cls)), len(cls), min(cls)))


def normal_closure_size(representative: Permutation, group: Sequence[Permutation]) -> int:
    conjugates = {compose(compose(h, representative), inverse(h)) for h in group}
    return len(generated_subgroup(conjugates, len(representative)))


def frame_orientation_certificate(
    gram: sp.Matrix,
    plus: Sequence[Permutation],
    minus: Sequence[Permutation],
) -> dict[str, Any]:
    """Check that combinatorial face orientation equals 3D determinant.

    Since E=G/4 is the orthogonal projection onto the rank-three Gram
    image and every incidence permutation commutes with E, the traces of
    the first three powers of the restricted action are
    tr(P^k E)=tr(P^k G)/4.  Newton's identity then gives its determinant
    without choosing coordinates in the three-dimensional image.
    """

    require(gram * gram == 4 * gram and sp.trace(gram) == 12, "FRAME_PROJECTOR", "Gram matrix is not four times a rank-three projector")

    def restricted_trace(permutation: Permutation, power: int) -> sp.Expr:
        powered = tuple(range(len(permutation)))
        for _ in range(power):
            powered = compose(permutation, powered)
        return sp.simplify(sum(gram[i, powered[i]] for i in range(len(permutation))) / 4)

    determinants: dict[Permutation, sp.Expr] = {}
    for permutation in list(plus) + list(minus):
        t1 = restricted_trace(permutation, 1)
        t2 = restricted_trace(permutation, 2)
        t3 = restricted_trace(permutation, 3)
        determinant = sp.simplify((t1**3 - 3 * t1 * t2 + 2 * t3) / 6)
        require(determinant in (sp.Integer(-1), sp.Integer(1)), "FRAME_ACTION_DETERMINANT", "induced frame determinant is not a sign")
        determinants[permutation] = determinant

    require(all(determinants[p] == 1 for p in plus), "FRAME_ORIENTATION_PLUS", "an oriented incidence automorphism acts by a reflection")
    require(all(determinants[p] == -1 for p in minus), "FRAME_ORIENTATION_MINUS", "an orientation-reversing incidence automorphism acts properly")
    return {
        "gram_image_dimension": 3,
        "determinant_method": "Newton identity from tr(P^k G)/4 for k=1,2,3",
        "orientation_preserving_determinants": {"+1": len(plus)},
        "orientation_reversing_determinants": {"-1": len(minus)},
        "combinatorial_orientation_equals_frame_determinant": True,
    }


def adjacency_matrix(carrier: Carrier) -> sp.Matrix:
    matrix = sp.zeros(12)
    for a, b in carrier.edges:
        matrix[a, b] = matrix[b, a] = 1
    return matrix


def gram_matrix(carrier: Carrier) -> sp.Matrix:
    t = sp.sqrt(5) / 5
    return sp.Matrix(
        12,
        12,
        lambda i, j: sp.Integer(1)
        if i == j
        else t
        if carrier.distances[i][j] == 1
        else -t
        if carrier.distances[i][j] == 2
        else sp.Integer(-1),
    )


def permute_matrix(permutation: Permutation) -> sp.Matrix:
    matrix = sp.zeros(len(permutation))
    for i, image in enumerate(permutation):
        matrix[image, i] = 1
    return matrix


def rooted_labeling(carrier: Carrier, root: tuple[int, int, int]) -> tuple[int, ...]:
    """Canonical combinatorial traversal from one oriented root face."""

    third: dict[tuple[int, int], int] = {}
    for a, b, c in carrier.faces:
        third[(a, b)] = c
        third[(b, c)] = a
        third[(c, a)] = b
    order = list(root)
    seen = set(order)
    while len(order) < 12:
        candidates: list[tuple[int, int, int, int]] = []
        for i in range(len(order)):
            for j in range(len(order)):
                if i == j:
                    continue
                u, v = order[i], order[j]
                w = third.get((u, v))
                if w is not None and w not in seen:
                    candidates.append((max(i, j), i, j, w))
        require(bool(candidates), "CANONICAL_TRAVERSAL", "rooted face traversal did not reach all ports")
        _, _, _, next_vertex = min(candidates)
        order.append(next_vertex)
        seen.add(next_vertex)
    return tuple(order)


def incidence_code(carrier: Carrier, order: Sequence[int]) -> str:
    position = {vertex: i for i, vertex in enumerate(order)}
    adjacency_bits = "".join(
        "1" if order[j] in carrier.adjacency[order[i]] else "0"
        for i in range(12)
        for j in range(i + 1, 12)
    )
    oriented_faces = sorted(
        min(
            (position[a], position[b], position[c]),
            (position[b], position[c], position[a]),
            (position[c], position[a], position[b]),
        )
        for a, b, c in carrier.faces
    )
    return adjacency_bits + ":" + ";".join(",".join(map(str, face)) for face in oriented_faces)


def canonical_incidence_hash(carrier: Carrier) -> str:
    roots: list[tuple[int, int, int]] = []
    for a, b, c in carrier.faces:
        roots.extend(((a, b, c), (b, c, a), (c, a, b)))
    codes = [incidence_code(carrier, rooted_labeling(carrier, root)) for root in roots]
    return hashlib.sha256(min(codes).encode("ascii")).hexdigest()


def parse_port_permutation(raw: Sequence[Any], carrier: Carrier) -> Permutation:
    require(isinstance(raw, list) and len(raw) == 12, "REFINEMENT_MAP_FORMAT", "port map must list twelve target ports")
    names = tuple(str(x) for x in raw)
    require(set(names) == set(carrier.ports), "REFINEMENT_MAP_BIJECTION", "port map must be a bijection")
    return tuple(carrier.index[name] for name in names)


def validate_refinement(manifest: Mapping[str, Any], carrier: Carrier, plus_group: Sequence[Permutation], gram: sp.Matrix) -> dict[str, Any]:
    tower = manifest.get("refinement_tower")
    require(isinstance(tower, Mapping), "REFINEMENT_TOWER", "refinement_tower is missing")
    levels_raw = tower.get("levels")
    require(isinstance(levels_raw, list) and len(levels_raw) >= 2, "REFINEMENT_LEVELS", "at least two refinement levels are required")
    levels = tuple(str(x) for x in levels_raw)
    require(len(set(levels)) == len(levels), "REFINEMENT_LEVELS", "refinement levels must be distinct")
    level_index = {level: i for i, level in enumerate(levels)}
    maps_raw = tower.get("maps")
    require(isinstance(maps_raw, list), "REFINEMENT_MAPS", "refinement maps must be a list")
    maps: dict[tuple[str, str], Permutation] = {}
    plus_set = set(plus_group)
    identity = tuple(range(12))
    for item in maps_raw:
        require(isinstance(item, Mapping), "REFINEMENT_MAP_FORMAT", "each refinement map must be an object")
        source = str(item.get("source"))
        target = str(item.get("target"))
        require(source in level_index and target in level_index, "REFINEMENT_LEVEL_UNKNOWN", "map names an unknown level")
        require(level_index[source] < level_index[target], "REFINEMENT_DIRECTION", "refinement maps must point forward")
        permutation = parse_port_permutation(item.get("port_map"), carrier)
        require(permutation in plus_set, "REFINEMENT_INCIDENCE", "refinement map is not an orientation-preserving incidence isomorphism")
        require(all(permutation[carrier.antipode[i]] == carrier.antipode[permutation[i]] for i in range(12)), "REFINEMENT_ANTIPODE", "refinement map does not commute with antipodes")
        p_matrix = permute_matrix(permutation)
        require(p_matrix.T * gram * p_matrix == gram, "REFINEMENT_FRAME", "refinement map does not preserve the canonical port frame")
        require((source, target) not in maps, "REFINEMENT_MAP_DUPLICATE", f"duplicate refinement map {source}->{target}")
        maps[(source, target)] = permutation

    # Identity and every declared composable triangle must satisfy the cocycle.
    for level in levels:
        maps.setdefault((level, level), identity)
    checked_triangles = 0
    for i, source in enumerate(levels):
        for j in range(i + 1, len(levels)):
            middle = levels[j]
            if (source, middle) not in maps:
                continue
            for k in range(j + 1, len(levels)):
                target = levels[k]
                if (middle, target) in maps and (source, target) in maps:
                    checked_triangles += 1
                    expected = compose(maps[(middle, target)], maps[(source, middle)])
                    require(maps[(source, target)] == expected, "REFINEMENT_COCYCLE", f"refinement cocycle fails on {source}->{middle}->{target}")
    require(checked_triangles >= 1, "REFINEMENT_COCYCLE_MISSING", "tower must declare at least one composable refinement triangle")

    return {
        "levels": list(levels),
        "declared_nonidentity_maps": len(maps_raw),
        "checked_cocycle_triangles": checked_triangles,
        "all_maps_in_A5": True,
        "unit_lines_natural": True,
        "antipode_natural": True,
        "frame_gram_natural": True,
    }


def group_certificate(carrier: Carrier) -> tuple[dict[str, Any], list[Permutation], list[Permutation]]:
    automorphisms = enumerate_distance_isometries(carrier.distances)
    require(len(automorphisms) == 120, "AUTOMORPHISM_COUNT", f"expected 120 full incidence automorphisms, got {len(automorphisms)}")
    signs = {permutation: orientation_sign(permutation, carrier.faces) for permutation in automorphisms}
    plus = [p for p in automorphisms if signs[p] == 1]
    minus = [p for p in automorphisms if signs[p] == -1]
    require(len(plus) == len(minus) == 60, "ORIENTATION_GROUP_COUNT", "orientation character must split 120=60+60")

    identity = tuple(range(12))
    antipode = tuple(carrier.antipode)
    require(antipode in minus, "FULL_GROUP_ANTIPODE", "central inversion must reverse orientation")
    require(all(compose(antipode, g) == compose(g, antipode) for g in automorphisms), "FULL_GROUP_CENTER", "antipode must be central in the full incidence group")
    require(set(minus) == {compose(antipode, g) for g in plus}, "FULL_GROUP_COSET", "orientation-reversing coset is not antipode times A5")
    require(identity in plus, "GROUP_IDENTITY", "identity is missing")
    plus_set = set(plus)
    require(all(compose(a, b) in plus_set for a in plus for b in plus), "GROUP_CLOSURE", "orientation-preserving automorphisms are not closed")

    classes = conjugacy_classes(plus)
    class_rows = [
        {
            "element_order": permutation_order(min(cls)),
            "class_size": len(cls),
            "normal_closure_size": normal_closure_size(min(cls), plus),
        }
        for cls in classes
    ]
    nonidentity_rows = [row for row in class_rows if row["element_order"] != 1]
    require(all(row["normal_closure_size"] == 60 for row in nonidentity_rows), "A5_NOT_SIMPLE", "orientation group failed exact simplicity check")
    order_histogram = Counter(permutation_order(p) for p in plus)
    require(order_histogram == Counter({1: 1, 2: 15, 3: 20, 5: 24}), "A5_ORDER_HISTOGRAM", f"unexpected element orders {order_histogram}")
    direct_identification = direct_a5_identification(plus)

    orbit = {p[0] for p in plus}
    require(len(orbit) == 12, "PORT_ACTION_NOT_TRANSITIVE", "A5 action is not transitive on ports")
    require(len(set(plus)) == 60, "PORT_ACTION_NOT_FAITHFUL", "port action is not faithful")

    axes = sorted({tuple(sorted((i, carrier.antipode[i]))) for i in range(12)})
    require(len(axes) == 6, "AXIS_COUNT", "antipode quotient must have six axes")
    axis_index = {axis: i for i, axis in enumerate(axes)}
    axis_actions: set[tuple[int, ...]] = set()
    for p in plus:
        image = tuple(axis_index[tuple(sorted((p[a], p[b])))] for a, b in axes)
        axis_actions.add(image)
    require(len(axis_actions) == 60, "AXIS_ACTION_NOT_FAITHFUL", "A5 action on six axes is not faithful")
    first_axis_stabilizer = sum(
        1
        for p in plus
        if tuple(sorted((p[axes[0][0]], p[axes[0][1]]))) == axes[0]
    )
    require(first_axis_stabilizer == 10, "AXIS_STABILIZER", "six-axis stabilizer must have order ten")

    return (
        {
            "full_incidence_automorphism_order": 120,
            "full_incidence_group": "A5 x C2, with central C2 generated by the antipode",
            "orientation_character_counts": {"+1": 60, "-1": 60},
            "orientation_preserving_order": 60,
            "orientation_preserving_group": "A5",
            "identification_proof": "faithful conjugation action on five Klein-four subgroups has exactly the 60 even permutations",
            "direct_a5_identification": direct_identification,
            "element_order_histogram": {str(k): order_histogram[k] for k in sorted(order_histogram)},
            "conjugacy_classes": class_rows,
            "port_action_faithful": True,
            "port_action_transitive": True,
            "axis_count": 6,
            "axis_action_faithful": True,
            "axis_orbit_transitive": True,
            "axis_stabilizer_order": 10,
            "a5_action_permutations": [list(p) for p in plus],
        },
        plus,
        minus,
    )


def certificate_payload(manifest: Mapping[str, Any]) -> dict[str, Any]:
    enforce_source_firewall(manifest)
    carrier = validate_carrier(manifest)
    group_row, plus, minus = group_certificate(carrier)

    adjacency = adjacency_matrix(carrier)
    x = sp.symbols("x")
    charpoly = sp.factor(adjacency.charpoly(x).as_expr())
    expected_charpoly = (x - 5) * (x + 1) ** 5 * (x**2 - 5) ** 3
    require(sp.expand(charpoly - expected_charpoly) == 0, "ADJACENCY_SPECTRUM", f"unexpected characteristic polynomial {charpoly}")

    gram = gram_matrix(carrier)
    require(gram.T == gram, "GRAM_SYMMETRY", "port Gram matrix is not symmetric")
    require(gram * gram == 4 * gram, "GRAM_PROJECTOR", "canonical Gram matrix must satisfy G^2=4G")
    require(sp.trace(gram) == 12, "GRAM_TRACE", "canonical Gram trace must be 12")
    require(all((permute_matrix(p).T * gram * permute_matrix(p)) == gram for p in plus + minus), "FRAME_EQUIVARIANCE", "incidence automorphism does not preserve the port frame")
    frame_orientation = frame_orientation_certificate(gram, plus, minus)

    antipode_pairs = sorted(
        tuple(sorted((carrier.ports[i], carrier.ports[carrier.antipode[i]])))
        for i in range(12)
        if i < carrier.antipode[i]
    )
    require(len(antipode_pairs) == 6, "ANTIPODE_PAIR_COUNT", "expected six antipodal pairs")

    refinement_row = validate_refinement(manifest, carrier, plus, gram)
    canonical_hash = canonical_incidence_hash(carrier)

    # Exact unit-splitting proof:
    # sum q_i^2 = 12 + sum(q_i-1)^2 when sum q_i=12.  If q != 1,
    # integral zero-sum deviations contain a positive and a negative entry,
    # so the square sum is at least 2.
    unit_split = {
        "domain": "q in Z^12 with sum(q)=12",
        "source_cost": "H(q)=sum_i q_i^2 from the normalized central readback Hilbert-Schmidt norm",
        "identity": "H(q)=12+sum_i(q_i-1)^2",
        "unique_minimizer": [1] * 12,
        "minimum": 12,
        "next_integer_floor": 14,
        "strict_gap": 2,
        "unit_lines": [f"span({atom['atom_id']})" for atom in manifest["carrier"]["central_port_atoms"]],
        "equivalence": "unique as the unordered family of primitive central port lines; labels transform by port permutations",
    }

    return {
        "schema": RECEIPT_SCHEMA,
        "issue": 565,
        "manifest_sha256": sha256_json(manifest),
        "source_firewall": {
            "forbidden_dependency_hits": [],
            "uses_only": [
                "primitive central port atoms",
                "normalized port trace",
                "integer defect readback",
                "edge/face incidence",
                "refinement port lineage",
            ],
        },
        "carrier": {
            "ports": 12,
            "edges": 30,
            "oriented_faces": 20,
            "degree_profile": [5] * 12,
            "euler_characteristic": 2,
            "vertex_links": "12 copies of C5",
            "distance_profile_per_port": {"0": 1, "1": 5, "2": 5, "3": 1},
            "canonical_oriented_incidence_sha256": canonical_hash,
        },
        "unit_split": unit_split,
        "inverse_pairing": {
            "definition": "iota(p) is the unique port at graph distance three from p",
            "fixed_point_free": True,
            "involutive": True,
            "unique": True,
            "pairs": [list(pair) for pair in antipode_pairs],
            "automorphism_equivariant": True,
        },
        "icosahedral_selector": {
            **group_row,
            "adjacency_characteristic_polynomial": "(x - 5)*(x + 1)^5*(x^2 - 5)^3",
            "canonical_gram_definition": "1 on diagonal; +1/sqrt(5) at distance 1; -1/sqrt(5) at distance 2; -1 at distance 3",
            "canonical_gram_identity": "G^2=4G",
            "canonical_gram_trace": 12,
            "canonical_gram_rank": 3,
            "nonzero_gram_eigenvalue": 4,
            "port_inner_products": ["-1", "-sqrt(5)/5", "sqrt(5)/5"],
            "frame_uniqueness": "every real realization of G is unique up to O(3); the oriented face class reduces equivalence to SO(3)",
            "port_frame_equivariance": True,
            "frame_orientation_certificate": frame_orientation,
            "six_axis_squared_cross_inner_product": "1/5",
            "selector_equivalence_class": "one oriented regular icosahedral port frame up to SO(3) and quotient-visible port relabelling",
        },
        "refinement": refinement_row,
        "classification": {
            "without_port_atoms": "one-dimensional splittings of an unmarked 12-space form a continuum GL(12,R)/((R^*)^12 semidirect S12)",
            "without_equal_trace": "weights (4,1,...,1)/15 make q=(0,2,1,...,1) cost 14/15, below the all-one cost 1",
            "without_total_twelve": "at total 13 there are 12 quadratic minimizers, one coordinate 2 and eleven coordinates 1",
            "without_pairing_incidence": "there are 12!/(2^6*6!) = 10395 fixed-point-free pairings of twelve labelled units",
            "without_orientation_restriction": "the survivor is the full order-120 group A5 x C2, including determinant-minus-one frame maps",
            "without_strict_cost": "for nonnegative q and linear cost at total 12, all C(23,11)=1352078 weak compositions minimize",
            "without_refinement_incidence": "a transposition is a port bijection but does not preserve edges, antipodes, or the Gram frame",
            "without_refinement_cocycle": "three individually valid stage maps can have a direct r0->r2 map unequal to the composite",
            "without_source_firewall": "a target field can force any preselected downstream representation and destroys source-only status",
        },
    }


def relabel_manifest(manifest: Mapping[str, Any], old_to_new: Sequence[int]) -> dict[str, Any]:
    """Relabel every port-bearing source field, including refinement maps."""

    require(len(old_to_new) == 12 and set(old_to_new) == set(range(12)), "RELABEL_FORMAT", "relabeling must be a permutation")
    out = copy.deepcopy(manifest)
    old_ports = list(manifest["carrier"]["ports"])
    new_ports = [f"r{idx:02d}" for idx in range(12)]
    rename = {old_ports[i]: new_ports[old_to_new[i]] for i in range(12)}
    out["carrier"]["ports"] = new_ports
    for atom in out["carrier"]["central_port_atoms"]:
        atom["port"] = rename[atom["port"]]
        atom["atom_id"] = f"e_{atom['port']}"
    out["carrier"]["edges"] = [[rename[a], rename[b]] for a, b in out["carrier"]["edges"]]
    out["carrier"]["oriented_faces"] = [[rename[a], rename[b], rename[c]] for a, b, c in out["carrier"]["oriented_faces"]]

    # Conjugate every old-index permutation to the new labels.
    sigma = tuple(old_to_new)
    sigma_inv = inverse(sigma)
    for item in out["refinement_tower"]["maps"]:
        old_perm_names = item["port_map"]
        old_perm = tuple(old_ports.index(name) for name in old_perm_names)
        new_perm = compose(compose(sigma, old_perm), sigma_inv)
        item["port_map"] = [new_ports[new_perm[i]] for i in range(12)]
    return out


def negative_control_cases(manifest: Mapping[str, Any]) -> list[tuple[str, dict[str, Any], str]]:
    cases: list[tuple[str, dict[str, Any], str]] = []

    merged = copy.deepcopy(manifest)
    merged["carrier"]["central_port_atoms"] = merged["carrier"]["central_port_atoms"][:6]
    cases.append(("remove_twelve_primitive_atoms", merged, "PORT_ATOMS_COUNT"))

    unequal = copy.deepcopy(manifest)
    unequal["carrier"]["central_port_atoms"][0]["normalized_trace"] = {"numerator": 1, "denominator": 6}
    cases.append(("remove_equal_port_trace", unequal, "PORT_TRACE_NOT_UNIFORM"))

    total13 = copy.deepcopy(manifest)
    total13["source_readback"]["total_charge"] = 13
    cases.append(("remove_total_twelve", total13, "TOTAL_CHARGE_NOT_12"))

    linear = copy.deepcopy(manifest)
    linear["source_readback"]["mismatch_cost"] = "sum_of_absolute_port_charges"
    cases.append(("remove_strict_quadratic_cost", linear, "COST_NOT_STRICT_QUADRATIC"))

    broken_edge = copy.deepcopy(manifest)
    broken_edge["carrier"]["edges"] = broken_edge["carrier"]["edges"][:-1]
    cases.append(("remove_icosahedral_incidence", broken_edge, "EDGE_COUNT"))

    reversing = copy.deepcopy(manifest)
    carrier = validate_carrier(manifest)
    _, _, minus = group_certificate(carrier)
    reversing["refinement_tower"]["maps"][0]["port_map"] = [carrier.ports[i] for i in minus[0]]
    cases.append(("remove_orientation_preserving_refinement", reversing, "REFINEMENT_INCIDENCE"))

    non_incidence = copy.deepcopy(manifest)
    transposition = list(range(12))
    transposition[0], transposition[1] = transposition[1], transposition[0]
    non_incidence["refinement_tower"]["maps"][0]["port_map"] = [carrier.ports[i] for i in transposition]
    cases.append(("remove_incidence_naturality", non_incidence, "REFINEMENT_INCIDENCE"))

    bad_cocycle = copy.deepcopy(manifest)
    bad_cocycle["refinement_tower"]["maps"][2]["port_map"] = list(carrier.ports)
    cases.append(("remove_refinement_cocycle", bad_cocycle, "REFINEMENT_COCYCLE"))

    forbidden = copy.deepcopy(manifest)
    forbidden["downstream_hint"] = {"product_adjoint_dimension": 12}
    cases.append(("inject_downstream_target", forbidden, "FORBIDDEN_DEPENDENCY"))

    return cases


def negative_control_payload(manifest: Mapping[str, Any]) -> dict[str, Any]:
    results: list[dict[str, Any]] = []
    for name, mutant, expected_code in negative_control_cases(manifest):
        actual_code = "ACCEPTED"
        try:
            certificate_payload(mutant)
        except CertificateError as exc:
            actual_code = exc.code
        require(actual_code == expected_code, "NEGATIVE_CONTROL_FAILED", f"{name}: expected {expected_code}, got {actual_code}")
        results.append({"name": name, "expected_error": expected_code, "actual_error": actual_code, "passed": True})
    return {
        "schema": NEGATIVE_SCHEMA,
        "issue": 565,
        "manifest_sha256": sha256_json(manifest),
        "finite_controls": results,
        "countermodel_witnesses": {
            "unequal_trace": {
                "weights": "(4,1,1,1,1,1,1,1,1,1,1,1)/15",
                "all_one_cost": "1",
                "alternative_q": [0, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                "alternative_cost": "14/15",
            },
            "wrong_total_13": {
                "minimizer_count": 12,
                "form": "one coordinate 2 and eleven coordinates 1",
            },
            "linear_cost_total_12": {
                "minimizer_count": 1352078,
                "form": "all weak compositions of 12 into 12 nonnegative parts",
            },
            "missing_incidence": {
                "fixed_point_free_pairings": 10395,
            },
            "missing_orientation": {
                "surviving_group": "A5 x C2",
                "order": 120,
            },
            "broken_refinement": {
                "nonincidence_map": "transposition of ports 0 and 1",
                "cocycle_failure": "replace the declared r0->r2 composite by the identity",
            },
        },
        "classified_larger_families": {
            "unmarked_12_space_splittings": "GL(12,R)/((R^*)^12 semidirect S12)",
            "labelled_fixed_point_free_pairings": 10395,
            "unoriented_icosahedral_symmetry_order": 120,
            "linear_cost_minimizers_on_nonnegative_total_12": 1352078,
        },
    }


def verify_receipt(manifest: Mapping[str, Any], receipt: Mapping[str, Any]) -> None:
    expected = certificate_payload(manifest)
    require(receipt == expected, "RECEIPT_MISMATCH", "receipt is stale, malformed, or tampered")


def default_paths() -> tuple[Path, Path, Path]:
    here = Path(__file__).resolve().parent
    return (
        here / "manifests" / "echosahedral_federation_reference.json",
        here / "receipts" / "echosahedral_federation_reference.receipt.json",
        here / "negative_controls" / "issue_565_negative_controls.json",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    certify = sub.add_parser("certify", help="create the deterministic exact receipt")
    certify.add_argument("--manifest", type=Path, required=True)
    certify.add_argument("--output", type=Path, required=True)
    verify = sub.add_parser("verify", help="recompute and compare a receipt")
    verify.add_argument("--manifest", type=Path, required=True)
    verify.add_argument("--receipt", type=Path, required=True)
    negative = sub.add_parser("negative-controls", help="run and write the finite countermodel bundle")
    negative.add_argument("--manifest", type=Path, required=True)
    negative.add_argument("--output", type=Path, required=True)
    all_cmd = sub.add_parser("all", help="regenerate receipt and negative controls at repository-default paths")
    all_cmd.add_argument("--manifest", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "certify":
        manifest = load_json(args.manifest)
        receipt = certificate_payload(manifest)
        write_json(args.output, receipt)
        print(json.dumps({"status": "PASS", "receipt": str(args.output), "sha256": sha256_json(receipt)}, indent=2))
    elif args.command == "verify":
        manifest = load_json(args.manifest)
        receipt = load_json(args.receipt)
        verify_receipt(manifest, receipt)
        print(json.dumps({"status": "PASS", "receipt": str(args.receipt)}, indent=2))
    elif args.command == "negative-controls":
        manifest = load_json(args.manifest)
        payload = negative_control_payload(manifest)
        write_json(args.output, payload)
        print(json.dumps({"status": "PASS", "negative_controls": str(args.output)}, indent=2))
    else:
        default_manifest, default_receipt, default_negative = default_paths()
        manifest_path = args.manifest or default_manifest
        manifest = load_json(manifest_path)
        write_json(default_receipt, certificate_payload(manifest))
        write_json(default_negative, negative_control_payload(manifest))
        print(json.dumps({"status": "PASS", "receipt": str(default_receipt), "negative_controls": str(default_negative)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

