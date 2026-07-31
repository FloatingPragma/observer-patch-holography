#!/usr/bin/env python3
"""Exact directed-seam repair bridge on the twelve-port carrier.

This certificate closes one finite algebraic bridge and nothing larger.
Starting from the unoriented incidence complex, it:

* constructs a coherent face orientation without coordinates;
* enumerates the 120 incidence automorphisms and proves that the 60
  orientation-preserving automorphisms act simply transitively on the 60
  directed seams;
* identifies the terminal of the #628 integer unit-transfer rule as one
  member of the two nearest-balanced directed placements;
* proves that the equal mixture of opposite directed placements is the
  rational pair-average update; and
* proves on the total-load-one sector that the uniform directed-seam orbit
  gives the exact channel ``I - Delta/60``.

For an odd endpoint total, neither integer placement lies in the exact
agreement equalizer.  Their mixture agrees only in expectation.  The
certificate therefore concerns a finite expected working-readback channel.
It does not identify a pathwise record repair, a physical clock, a
refinement-compatible semigroup, or a universe-selecting law.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import deque
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any, Mapping, Sequence

MODULE_DIR = Path(__file__).resolve().parent
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

import a3_scheduler_kernel_certificate as scheduler  # noqa: E402
import echosahedral_selector_certificate as carrier_cert  # noqa: E402
import record_counting_mechanism_certificate as counting  # noqa: E402

CertificateError = carrier_cert.CertificateError
require = carrier_cert.require
load_json = carrier_cert.load_json
sha256_json = carrier_cert.sha256_json
write_json = carrier_cert.write_json

SCHEMA = "oph.directed_seam_repair_certificate.v1"
PORT_COUNT = 12
EDGE_COUNT = 30
DIRECTED_SEAM_COUNT = 60
AUDIT_BOUND = 12
MANIFEST_PATH = MODULE_DIR / "manifests" / "directed_seam_repair_reference.json"

Edge = tuple[int, int]
DirectedSeam = tuple[int, int]
Permutation = tuple[int, ...]
Matrix = tuple[tuple[Fraction, ...], ...]


def port_label(index: int) -> str:
    return f"p{index:02d}"


def canonical_cycle(face: Sequence[int]) -> tuple[int, int, int]:
    """Choose a canonical cyclic representative without reversing a face."""

    require(len(face) == 3, "FACE_ARITY", "a face must contain three vertices")
    a, b, c = face
    require(len({a, b, c}) == 3, "FACE_DEGENERATE", "face vertices must differ")
    return min((a, b, c), (b, c, a), (c, a, b))


def unoriented_triangles(
    adjacency: Sequence[Sequence[int]],
) -> tuple[tuple[int, int, int], ...]:
    """Enumerate the triangular 2-cells from graph incidence alone."""

    n = len(adjacency)
    neighbors = [frozenset(row) for row in adjacency]
    faces = tuple(
        (a, b, c)
        for a, b, c in combinations(range(n), 3)
        if b in neighbors[a] and c in neighbors[a] and c in neighbors[b]
    )
    require(
        len(faces) == 20,
        "FACE_COUNT",
        f"the unoriented carrier must have twenty triangular faces, got {len(faces)}",
    )
    return faces


def coherent_orientation(
    adjacency: Sequence[Sequence[int]],
) -> tuple[tuple[int, int, int], ...]:
    """Orient the closed triangulation by propagation across shared edges.

    The lexicographically first unoriented face supplies an arbitrary global
    sign.  Every neighboring face is oriented so that a shared edge receives
    the opposite boundary direction.  Reversing the seed reverses every face
    and leaves the orientation-preserving automorphism subgroup unchanged.
    """

    faces = unoriented_triangles(adjacency)
    face_keys = tuple(frozenset(face) for face in faces)
    edge_to_faces: dict[Edge, list[frozenset[int]]] = {}
    for key in face_keys:
        for u, v in combinations(sorted(key), 2):
            edge_to_faces.setdefault((u, v), []).append(key)
    require(
        len(edge_to_faces) == EDGE_COUNT,
        "EDGE_COUNT",
        f"the face complex must have thirty seams, got {len(edge_to_faces)}",
    )
    require(
        all(len(incident) == 2 for incident in edge_to_faces.values()),
        "EDGE_FACE_INCIDENCE",
        "every seam must belong to exactly two triangular faces",
    )

    seed = face_keys[0]
    oriented: dict[frozenset[int], tuple[int, int, int]] = {
        seed: canonical_cycle(tuple(sorted(seed)))
    }
    queue: deque[frozenset[int]] = deque([seed])
    while queue:
        key = queue.popleft()
        a, b, c = oriented[key]
        for u, v in ((a, b), (b, c), (c, a)):
            edge = (min(u, v), max(u, v))
            neighbor = next(item for item in edge_to_faces[edge] if item != key)
            third = next(vertex for vertex in neighbor if vertex not in edge)
            proposed = canonical_cycle((v, u, third))
            if neighbor in oriented:
                require(
                    oriented[neighbor] == proposed,
                    "NONORIENTABLE",
                    "orientation propagation is inconsistent",
                )
            else:
                oriented[neighbor] = proposed
                queue.append(neighbor)

    require(
        len(oriented) == len(faces),
        "FACE_CONNECTIVITY",
        "the face-adjacency graph must be connected",
    )
    directed_counts: dict[DirectedSeam, int] = {}
    for a, b, c in oriented.values():
        for directed in ((a, b), (b, c), (c, a)):
            directed_counts[directed] = directed_counts.get(directed, 0) + 1
    for u, v in edge_to_faces:
        require(
            directed_counts.get((u, v)) == 1
            and directed_counts.get((v, u)) == 1,
            "ORIENTATION_BOUNDARY",
            "the two incident faces must induce opposite seam directions",
        )
    return tuple(sorted(oriented.values()))


def orientation_sign(
    permutation: Permutation,
    oriented_faces: Sequence[tuple[int, int, int]],
) -> int:
    """Return the combinatorial surface-orientation character."""

    positive: set[tuple[int, int, int]] = set()
    negative: set[tuple[int, int, int]] = set()
    for a, b, c in oriented_faces:
        positive.update(((a, b, c), (b, c, a), (c, a, b)))
        negative.update(((a, c, b), (c, b, a), (b, a, c)))
    signs: set[int] = set()
    for face in oriented_faces:
        image = tuple(permutation[vertex] for vertex in face)
        if image in positive:
            signs.add(1)
        elif image in negative:
            signs.add(-1)
        else:
            raise CertificateError(
                "AUTOMORPHISM_FACE",
                "an incidence automorphism failed to preserve the face complex",
            )
    require(
        len(signs) == 1,
        "AUTOMORPHISM_ORIENTATION",
        "an automorphism has no consistent orientation character",
    )
    return next(iter(signs))


def directed_seams(edges: Sequence[Edge]) -> tuple[DirectedSeam, ...]:
    require(
        len(edges) == EDGE_COUNT,
        "EDGE_COUNT",
        "the directed primitive orbit requires thirty unoriented seams",
    )
    result = tuple(sorted((u, v) for edge in edges for u, v in (edge, edge[::-1])))
    require(
        len(result) == DIRECTED_SEAM_COUNT and len(set(result)) == len(result),
        "DIRECTED_SEAM_COUNT",
        "the carrier must have sixty distinct directed seams",
    )
    return result


def incidence_and_orbit() -> dict[str, Any]:
    """Build the combinatorial orientation and regular directed-seam action."""

    source_manifest = load_json(
        MODULE_DIR / "manifests" / counting.CARRIER_MANIFEST_NAME
    )
    carrier = carrier_cert.validate_carrier(source_manifest)
    faces = coherent_orientation(carrier.adjacency)
    automorphisms = carrier_cert.enumerate_distance_isometries(carrier.distances)
    require(
        len(automorphisms) == 120,
        "AUTOMORPHISM_COUNT",
        "the incidence automorphism group must have order 120",
    )
    signs = {
        permutation: orientation_sign(permutation, faces)
        for permutation in automorphisms
    }
    proper = tuple(sorted(p for p in automorphisms if signs[p] == 1))
    reversing = tuple(sorted(p for p in automorphisms if signs[p] == -1))
    require(
        len(proper) == len(reversing) == 60,
        "PROPER_GROUP_COUNT",
        "the orientation character must split the automorphisms as 60 plus 60",
    )

    seams = directed_seams(carrier.edges)
    seam_set = set(seams)
    identity = tuple(range(PORT_COUNT))
    require(identity in proper, "PROPER_IDENTITY", "the proper group lacks identity")
    proper_set = set(proper)
    require(
        all(
            carrier_cert.compose(left, right) in proper_set
            for left in proper
            for right in proper
        ),
        "PROPER_CLOSURE",
        "the proper incidence automorphisms are not closed",
    )

    representative = seams[0]
    representative_orbit = tuple(
        (permutation[representative[0]], permutation[representative[1]])
        for permutation in proper
    )
    require(
        set(representative_orbit) == seam_set
        and len(set(representative_orbit)) == DIRECTED_SEAM_COUNT,
        "DIRECTED_ORBIT",
        "the proper group is not transitive and free on directed seams",
    )
    stabilizer = tuple(
        permutation
        for permutation in proper
        if (
            permutation[representative[0]],
            permutation[representative[1]],
        )
        == representative
    )
    require(
        stabilizer == (identity,),
        "DIRECTED_STABILIZER",
        "a directed seam must have trivial proper stabilizer",
    )

    transporter_checks = 0
    for source in seams:
        images = [
            (permutation[source[0]], permutation[source[1]])
            for permutation in proper
        ]
        require(
            len(set(images)) == DIRECTED_SEAM_COUNT and set(images) == seam_set,
            "SIMPLE_TRANSITIVITY",
            "some directed seam lacks a unique proper transporter to every target",
        )
        transporter_checks += len(images)

    return {
        "carrier": carrier,
        "source_manifest": source_manifest,
        "faces": faces,
        "automorphisms": tuple(automorphisms),
        "proper": proper,
        "reversing": reversing,
        "seams": seams,
        "representative": representative,
        "transporter_checks": transporter_checks,
    }


def nearest_balanced_pair(left: int, right: int) -> tuple[int, int]:
    """Floor/ceiling shell for the conserved endpoint total."""

    total = left + right
    lower = total // 2
    upper = total - lower
    require(
        upper - lower in (0, 1),
        "BALANCED_SHELL",
        "integer floor and ceiling must differ by at most one",
    )
    return lower, upper


def directed_balanced_pair(left: int, right: int) -> tuple[int, int]:
    """Put the floor at the directed source and the ceiling at its target."""

    return nearest_balanced_pair(left, right)


def opposite_directed_balanced_pair(left: int, right: int) -> tuple[int, int]:
    lower, upper = nearest_balanced_pair(left, right)
    return upper, lower


def directed_pair_expectation(left: int, right: int) -> tuple[Fraction, Fraction]:
    """Equal mixture of the two opposite directed balanced placements."""

    forward = directed_balanced_pair(left, right)
    reverse = opposite_directed_balanced_pair(left, right)
    return (
        Fraction(forward[0] + reverse[0], 2),
        Fraction(forward[1] + reverse[1], 2),
    )


def pair_average(left: int | Fraction, right: int | Fraction) -> tuple[Fraction, Fraction]:
    mean = (Fraction(left) + Fraction(right)) / 2
    return mean, mean


def settle_pair_via_628(
    left: int,
    right: int,
) -> tuple[tuple[int, int], tuple[tuple[int, int], ...]]:
    """Replay #628 on one selected seam until its balanced shell is reached."""

    state = (left, right)
    history: list[tuple[int, int]] = []
    edge = [(0, 1)]
    while True:
        moves = counting.admissible_moves(state, edge)
        if not moves:
            break
        require(
            len(moves) == 1,
            "MICRO_MOVE_UNIQUENESS",
            "one selected seam must have at most one descending direction",
        )
        move = moves[0]
        before = counting.descent_potential(state)
        difference = state[move[0]] - state[move[1]]
        state = counting.apply_move(state, move)
        require(
            counting.descent_potential(state) - before == -2 * (difference - 1),
            "MICRO_DESCENT",
            "the #628 microstep must obey its exact descent identity",
        )
        history.append(move)
    return state, tuple(history)


def integer_macro_bridge() -> dict[str, Any]:
    """Exact bounded replay plus the parity-complete algebraic identities."""

    pair_checks = 0
    microstep_checks = 0
    odd_checks = 0
    even_checks = 0
    for left in range(-AUDIT_BOUND, AUDIT_BOUND + 1):
        for right in range(-AUDIT_BOUND, AUDIT_BOUND + 1):
            lower, upper = nearest_balanced_pair(left, right)
            terminal, history = settle_pair_via_628(left, right)
            require(
                terminal in ((lower, upper), (upper, lower)),
                "MICRO_TERMINAL",
                "the #628 microdynamics did not reach the nearest-balanced shell",
            )
            require(
                len(history) == abs(left - right) // 2,
                "MICRO_LENGTH",
                "the #628 microhistory length differs from floor(|a-b|/2)",
            )
            require(
                sum(terminal) == left + right,
                "MICRO_CONSERVATION",
                "the #628 microhistory failed to conserve endpoint total",
            )
            require(
                directed_pair_expectation(left, right)
                == pair_average(left, right),
                "PAIR_EXPECTATION",
                "opposite balanced placements did not average to pair agreement",
            )
            pair_checks += 1
            microstep_checks += len(history)
            if (left + right) % 2:
                require(
                    lower != upper
                    and abs(terminal[0] - terminal[1]) == 1,
                    "ODD_SHELL",
                    "an odd endpoint total must retain pathwise mismatch one",
                )
                odd_checks += 1
            else:
                require(
                    lower == upper and terminal[0] == terminal[1],
                    "EVEN_EQUALIZER",
                    "an even endpoint total must reach exact integer agreement",
                )
                even_checks += 1

    return {
        "audit_integer_range": [-AUDIT_BOUND, AUDIT_BOUND],
        "audited_input_pairs": pair_checks,
        "audited_628_microsteps": microstep_checks,
        "even_total_cases": even_checks,
        "odd_total_cases": odd_checks,
        "micro_rule": "transfer one unit from the higher endpoint when the difference is at least two",
        "micro_history_length": "floor(abs(a-b)/2)",
        "terminal_shell": "{(floor((a+b)/2),ceil((a+b)/2)), (ceil((a+b)/2),floor((a+b)/2))}",
        "directed_placement": "(floor((a+b)/2),ceil((a+b)/2))",
        "opposite_placement": "(ceil((a+b)/2),floor((a+b)/2))",
        "equal_direction_expectation": "((a+b)/2,(a+b)/2)",
        "expectation_equals_pair_average": True,
        "pathwise_exact_agreement": "even endpoint totals only",
        "odd_total_statement": "the two directed placements label the nearest-balanced shell; a fixed #628 descent history reaches the member selected by its initial load ordering",
    }


def identity_matrix(size: int) -> Matrix:
    return tuple(
        tuple(Fraction(int(i == j)) for j in range(size)) for i in range(size)
    )


def pair_average_matrix(size: int, seam: DirectedSeam) -> Matrix:
    u, v = seam
    require(
        0 <= u < size and 0 <= v < size and u != v,
        "SEAM_RANGE",
        "a seam must join distinct in-range ports",
    )
    rows = [list(row) for row in identity_matrix(size)]
    rows[u][u] = Fraction(1, 2)
    rows[u][v] = Fraction(1, 2)
    rows[v][u] = Fraction(1, 2)
    rows[v][v] = Fraction(1, 2)
    return tuple(tuple(row) for row in rows)


def graph_laplacian(size: int, edges: Sequence[Edge]) -> Matrix:
    rows = [[Fraction(0) for _ in range(size)] for _ in range(size)]
    for u, v in edges:
        rows[u][u] += 1
        rows[v][v] += 1
        rows[u][v] -= 1
        rows[v][u] -= 1
    return tuple(tuple(row) for row in rows)


def matrix_add(left: Matrix, right: Matrix) -> Matrix:
    require(
        len(left) == len(right)
        and all(len(a) == len(b) for a, b in zip(left, right)),
        "MATRIX_SHAPE",
        "matrix shapes differ",
    )
    return tuple(
        tuple(a + b for a, b in zip(left_row, right_row))
        for left_row, right_row in zip(left, right)
    )


def matrix_scale(scalar: Fraction, matrix: Matrix) -> Matrix:
    return tuple(tuple(scalar * entry for entry in row) for row in matrix)


def matrix_vector(matrix: Matrix, vector: Sequence[Fraction]) -> tuple[Fraction, ...]:
    require(
        all(len(row) == len(vector) for row in matrix),
        "MATRIX_VECTOR_SHAPE",
        "matrix and vector shapes differ",
    )
    return tuple(
        sum((entry * value for entry, value in zip(row, vector)), Fraction(0))
        for row in matrix
    )


def expected_directed_update(
    counts: Sequence[int],
    seam: DirectedSeam,
) -> tuple[Fraction, ...]:
    """Expected full state under opposite balanced placements of one seam."""

    u, v = seam
    require(
        len(counts) == PORT_COUNT and u != v,
        "STATE_SHAPE",
        "the directed update needs twelve loads and distinct endpoints",
    )
    out = [Fraction(value) for value in counts]
    out[u], out[v] = directed_pair_expectation(counts[u], counts[v])
    return tuple(out)


def directed_balanced_update(
    counts: Sequence[int],
    seam: DirectedSeam,
) -> tuple[int, ...]:
    """Apply one deterministic directed nearest-balanced placement."""

    u, v = seam
    require(
        len(counts) == PORT_COUNT and u != v,
        "STATE_SHAPE",
        "the directed update needs twelve loads and distinct endpoints",
    )
    out = list(counts)
    out[u], out[v] = directed_balanced_pair(counts[u], counts[v])
    return tuple(out)


def directed_s1_matrix(seam: DirectedSeam) -> Matrix:
    """Linear matrix induced by a directed placement on the S=1 basis."""

    columns = []
    for source in range(PORT_COUNT):
        basis = tuple(1 if index == source else 0 for index in range(PORT_COUNT))
        columns.append(
            tuple(Fraction(value) for value in directed_balanced_update(basis, seam))
        )
    return tuple(
        tuple(columns[column][row] for column in range(PORT_COUNT))
        for row in range(PORT_COUNT)
    )


def s1_channel_certificate(
    edges: Sequence[Edge],
    seams: Sequence[DirectedSeam],
) -> dict[str, Any]:
    """Prove the directed-orbit channel identity on total load one."""

    require(
        tuple(sorted(seams)) == directed_seams(edges),
        "DIRECTED_ORBIT_SUPPORT",
        "the channel support is not the complete directed-seam orbit",
    )
    matrices: list[Matrix] = []
    basis_checks = 0
    opposite_pair_checks = 0
    for seam in seams:
        matrix = directed_s1_matrix(seam)
        reverse_matrix = directed_s1_matrix((seam[1], seam[0]))
        require(
            matrix_scale(
                Fraction(1, 2),
                matrix_add(matrix, reverse_matrix),
            )
            == pair_average_matrix(PORT_COUNT, seam),
            "OPPOSITE_DIRECTION_MEAN",
            "opposite directed S=1 placements did not average to pairAverage",
        )
        opposite_pair_checks += 1
        for source in range(PORT_COUNT):
            basis = tuple(1 if index == source else 0 for index in range(PORT_COUNT))
            vector = tuple(Fraction(value) for value in basis)
            require(
                tuple(
                    Fraction(value)
                    for value in directed_balanced_update(basis, seam)
                )
                == matrix_vector(matrix, vector),
                "S1_PRIMITIVE",
                "a directed balanced primitive disagrees with its S=1 matrix",
            )
            basis_checks += 1
        matrices.append(matrix)

    channel = matrix_scale(
        Fraction(1, DIRECTED_SEAM_COUNT),
        sum_matrices(matrices),
    )
    laplacian = graph_laplacian(PORT_COUNT, edges)
    expected = matrix_add(
        identity_matrix(PORT_COUNT),
        matrix_scale(Fraction(-1, 60), laplacian),
    )
    require(
        channel == expected,
        "UNIFORM_CHANNEL",
        "the uniform directed orbit did not give I - Delta/60",
    )
    require(
        all(
            sum((channel[row][column] for row in range(PORT_COUNT)), Fraction(0))
            == 1
            for column in range(PORT_COUNT)
        )
        and all(entry >= 0 for row in channel for entry in row),
        "S1_STOCHASTIC",
        "the total-load-one channel must be exactly stochastic",
    )
    require(
        all(channel[i][i] == Fraction(11, 12) for i in range(PORT_COUNT)),
        "S1_DIAGONAL",
        "the five-regular channel must have diagonal 11/12",
    )
    for u, v in edges:
        require(
            channel[u][v] == channel[v][u] == Fraction(1, 60),
            "S1_NEIGHBOR",
            "each adjacent transition must have weight 1/60",
        )

    projected_weights = {
        edge: sum(
            Fraction(1, DIRECTED_SEAM_COUNT)
            for seam in seams
            if frozenset(seam) == frozenset(edge)
        )
        for edge in edges
    }
    require(
        all(weight == Fraction(1, EDGE_COUNT) for weight in projected_weights.values()),
        "UNDIRECTED_PROJECTION",
        "uniform directed scheduling must project to the #614 1/30 seam kernel",
    )

    return {
        "sector": "S=1 total integer working load",
        "primitive_basis_checks": basis_checks,
        "opposite_direction_pair_average_checks": opposite_pair_checks,
        "directed_probability": "1/60",
        "projected_undirected_probability": "1/30",
        "opposite_direction_expectation": "pairAverage",
        "channel_identity": "R = I - Delta_ico/60",
        "positive_generator_identity": "I - R = Delta_ico/60",
        "diagonal_entry": "11/12",
        "adjacent_entry": "1/60",
        "nonadjacent_entry": "0",
        "column_stochastic": True,
        "exact_matrix": fraction_matrix_json(channel),
        "exact_laplacian": fraction_matrix_json(laplacian),
    }


def sum_matrices(matrices: Sequence[Matrix]) -> Matrix:
    require(bool(matrices), "MATRIX_EMPTY", "at least one matrix is required")
    total = tuple(
        tuple(Fraction(0) for _ in range(len(matrices[0][0])))
        for _ in range(len(matrices[0]))
    )
    for matrix in matrices:
        total = matrix_add(total, matrix)
    return total


def fraction_matrix_json(matrix: Matrix) -> list[list[str]]:
    return [[str(entry) for entry in row] for row in matrix]


def validate_upstream_manifest(
    path: Path,
    schema: str,
    issue: int | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    value = load_json(path)
    require(
        isinstance(value, dict) and value.get("schema") == schema,
        "UPSTREAM_SCHEMA",
        f"{path.name} does not carry {schema}",
    )
    if issue is not None:
        require(
            value.get("issue") == issue,
            "UPSTREAM_ISSUE",
            f"{path.name} is not the expected issue-{issue} packet",
        )
    if "manifest_sha256" in value:
        body = {key: item for key, item in value.items() if key != "manifest_sha256"}
        require(
            value["manifest_sha256"] == "sha256:" + sha256_json(body),
            "UPSTREAM_HASH",
            f"{path.name} has an invalid manifest hash",
        )
    pin = {
        "path": f"manifests/{path.name}",
        "schema": schema,
        "canonical_json_sha256": "sha256:" + sha256_json(value),
    }
    return value, pin


def build_payload() -> dict[str, Any]:
    orbit = incidence_and_orbit()
    carrier = orbit["carrier"]

    _, carrier_pin = validate_upstream_manifest(
        MODULE_DIR / "manifests" / counting.CARRIER_MANIFEST_NAME,
        carrier_cert.SCHEMA,
    )
    scheduler_manifest, scheduler_pin = validate_upstream_manifest(
        MODULE_DIR / "manifests" / "a3_scheduler_kernel_reference.json",
        scheduler.SCHEMA,
        scheduler.ISSUE,
    )
    counting_manifest, counting_pin = validate_upstream_manifest(
        counting.MANIFEST_PATH,
        counting.SCHEMA,
        counting.ISSUE,
    )
    require(
        scheduler_manifest["a3_selection"]["selected_kernel_probability_per_seam"]
        == "1/30",
        "SCHEDULER_PIN",
        "the pinned scheduler must select the exact undirected 1/30 kernel",
    )
    require(
        counting_manifest["bounded_exit"] == "exact_named_realization",
        "COUNTING_PIN",
        "the pinned #628 mechanism must retain its bounded exact scope",
    )

    faces = orbit["faces"]
    representative = orbit["representative"]
    macro = integer_macro_bridge()
    s1 = s1_channel_certificate(carrier.edges, orbit["seams"])

    return {
        "schema": SCHEMA,
        "scope": {
            "object": "finite expected scalar working-readback repair on one twelve-port carrier",
            "arithmetic": "integers and exact rational Fractions only",
            "coordinate_geometry_used": False,
            "source_issues": [614, 628],
            "producer_replay": True,
            "independent_implementation": False,
            "claim_boundary": (
                "The certificate proves a one-step expected scalar channel "
                "and a local integer macro bridge. It does not prove an IID "
                "path law, exact odd-total agreement, a record instrument, "
                "refinement, a physical clock, or universe selection."
            ),
        },
        "upstream_pins": {
            "carrier": carrier_pin,
            "undirected_scheduler": scheduler_pin,
            "integer_record_counting_mechanism": counting_pin,
            "lean_pair_average": {
                "path": "Lean/ObserverPatchHolography/ScalarSeamRepair.lean",
                "theorems": [
                    "eq_pairAverage_of_supported_agreeing_sum_preserving",
                    "pairAverage_eq_id_sub_half_edgeLaplacian",
                    "uniform_thirty_seam_repair",
                ],
                "role": "formal counterpart; this certificate recomputes its scalar identities independently in exact rational arithmetic",
            },
        },
        "combinatorial_orientation": {
            "source": "unoriented adjacency and triangular incidence",
            "seed_rule": "lexicographically first triangular face with increasing vertex labels",
            "propagation_rule": "shared seams receive opposite boundary directions",
            "oriented_faces": [
                [port_label(a), port_label(b), port_label(c)] for a, b, c in faces
            ],
            "faces": len(faces),
            "seams": len(carrier.edges),
            "coherent_closed_orientation": True,
            "global_sign_is_conventional": True,
        },
        "directed_seam_primitive_orbit": {
            "full_incidence_automorphisms": len(orbit["automorphisms"]),
            "orientation_preserving_automorphisms": len(orbit["proper"]),
            "orientation_reversing_automorphisms": len(orbit["reversing"]),
            "directed_seams": len(orbit["seams"]),
            "representative": [
                port_label(representative[0]),
                port_label(representative[1]),
            ],
            "representative_orbit_size": len(orbit["seams"]),
            "representative_stabilizer_order": 1,
            "unique_transporter_checks": orbit["transporter_checks"],
            "simply_transitive": True,
            "unique_invariant_probability_on_the_orbit": "1/60",
        },
        "integer_macro_bridge": macro,
        "uniform_s1_channel": s1,
        "verdict": {
            "combinatorial_directed_orbit": "exact",
            "integer_micro_to_balanced_shell": "exact",
            "opposite_direction_expectation_to_pair_average": "exact",
            "uniform_s1_channel": "exact_named_realization",
            "pathwise_odd_total_agreement": "not attained",
            "opposite_shell_placement_as_628_descent": "not generally attained",
            "directed_orbit_schedule_as_physical_selection": "open",
            "protected_record_and_checkpoint_instrument": "open",
            "noncommutative_common_repair_workspace": "open",
            "refinement_semigroup_naturality": "open",
            "physical_clock_and_laboratory_attachment": "open",
            "full_self_readback_and_universe_selection": "open",
        },
    }


def signed_manifest() -> dict[str, Any]:
    payload = build_payload()
    return {**payload, "manifest_sha256": "sha256:" + sha256_json(payload)}


def verify_manifest(manifest: Mapping[str, Any]) -> dict[str, Any]:
    require(
        isinstance(manifest, Mapping),
        "MANIFEST",
        "the directed-seam manifest must be an object",
    )
    body = {key: value for key, value in manifest.items() if key != "manifest_sha256"}
    require(
        manifest.get("manifest_sha256") == "sha256:" + sha256_json(body),
        "MANIFEST_HASH",
        "the directed-seam manifest hash is invalid",
    )
    require(
        body == build_payload(),
        "MANIFEST_REPLAY",
        "the directed-seam manifest does not match exact recomputation",
    )
    return {
        "status": "PASS",
        "schema": SCHEMA,
        "manifest_sha256": manifest["manifest_sha256"],
        "producer_replay": True,
        "independent_implementation": False,
        "claim_boundary": body["scope"]["claim_boundary"],
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    certify = subparsers.add_parser("certify")
    certify.add_argument("--output", type=Path, default=MANIFEST_PATH)
    verify = subparsers.add_parser("verify")
    verify.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "certify":
            manifest = signed_manifest()
            write_json(args.output, manifest)
            print(
                json.dumps(
                    {
                        "status": "PASS",
                        "schema": SCHEMA,
                        "output": str(args.output),
                        "manifest_sha256": manifest["manifest_sha256"],
                    },
                    sort_keys=True,
                )
            )
        else:
            print(json.dumps(verify_manifest(load_json(args.manifest)), sort_keys=True))
    except (CertificateError, OSError, KeyError, TypeError, ValueError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
