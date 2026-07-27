#!/usr/bin/env python3
"""Exact certificate for GitHub issue #614: the bounded A3 scheduler campaign.

The reference carrier is the twelve-port icosahedral federation shared by the
a5_closure certificates, presented through the exact coordinate model with
vertices at the cyclic permutations of (0, +/-1, +/-phi) over Q(sqrt5). The
move simplex consists of the thirty undirected seams of that carrier, and
every transition observable is the source-defined record of which seam
repairs next.

The certificate derives, with exact Fraction and Q(sqrt5) arithmetic and no
floating point in any proof step:

* the A3 selection: with the complete constraint grammar contributing
  normalization only, the information projection of the uniform reference on
  the move simplex is the uniform kernel with probability 1/30 per seam; the
  order-120 carrier automorphism group acts transitively on the seams, so a
  deck-invariant positive weight vector normalizes to the same point;
* the proved dynamics of the induced port walk W = A/5: positive support,
  irreducibility by exact reachability closure, aperiodicity from an
  exhibited triangle, and exact contraction with adjacency spectrum
  5, sqrt5, -sqrt5, -1 of multiplicities 1, 3, 3, 5, verified through the
  minimal polynomial identity (A - 5I)(A^2 - 5I)(A + I) = 0 and the exact
  trace system on powers of A; the walk spectral gap is 1 - sqrt5/5;
* the exact fairness and excursion bounds of the uniform kernel: a seam is
  unvisited after t steps with probability (29/30)^t, the expected
  first-visit time is 30 steps, and the 120-step excursion bound is recorded
  as a reduced fraction;
* four fail-closed controls: a zeroed seam weight starves that seam forever,
  a disconnected move subset fails the reachability cover, a bipartite move
  subset has period two, and the deterministic lowest-index-first scheduler
  produces presentation-dependent first-move records under a nonidentity
  deck rotation while every quotient-visible datum agrees.

Fairness and liveness for a general implementation model stay a
conditional_open_interface premise, per the bounded exit of the issue.
"""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Sequence

MODULE_DIR = Path(__file__).resolve().parent
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

import echosahedral_selector_certificate as e565  # noqa: E402
import port_current_inner_certificate as p584  # noqa: E402

CertificateError = e565.CertificateError
require = e565.require
sha256_json = e565.sha256_json
load_json = e565.load_json
write_json = e565.write_json

F5 = p584.F5
PHI = p584.PHI
ZERO = p584.ZERO
ONE = p584.ONE
standard_vertices = p584.standard_vertices
dot = p584.dot

SCHEMA = "oph.a3_scheduler_kernel_certificate.v1"
ISSUE = 614
PORT_COUNT = 12
SEAM_COUNT = 30
DEGREE = 5
FAIRNESS_HORIZON_STEPS = 120

SQRT5 = F5(0, 1)
FIVE = F5(5)

# Exact integer rotation matrices tried, in order, by the presentation
# control: two coordinate three-cycles and three axis half-turns. Each has
# determinant one and maps the vertex set to itself.
ROTATION_CANDIDATES: tuple[tuple[tuple[int, int, int], ...], ...] = (
    ((0, 1, 0), (0, 0, 1), (1, 0, 0)),
    ((0, 0, 1), (1, 0, 0), (0, 1, 0)),
    ((1, 0, 0), (0, -1, 0), (0, 0, -1)),
    ((-1, 0, 0), (0, 1, 0), (0, 0, -1)),
    ((-1, 0, 0), (0, -1, 0), (0, 0, 1)),
)

IDENTITY_ROTATION: tuple[tuple[int, int, int], ...] = (
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
)

Edge = tuple[int, int]


def port_label(index: int) -> str:
    return f"p{index:02d}"


# ---------------------------------------------------------------------------
# Move simplex: the thirty seams of the icosahedral reference carrier
# ---------------------------------------------------------------------------


def build_move_simplex() -> tuple[list[tuple[F5, F5, F5]], list[list[int]], list[Edge]]:
    """The exact carrier: twelve ports, adjacency, and thirty seams."""

    vertices = standard_vertices()
    require(
        len(vertices) == PORT_COUNT and len(set(vertices)) == PORT_COUNT,
        "MOVE_SIMPLEX",
        "the coordinate model does not have twelve distinct ports",
    )
    adjacency: list[list[int]] = []
    for i in range(PORT_COUNT):
        neighbors = sorted(
            j
            for j in range(PORT_COUNT)
            if j != i and dot(vertices[i], vertices[j]) == PHI
        )
        require(
            len(neighbors) == DEGREE,
            "MOVE_SIMPLEX",
            f"port {i} does not have exactly five seams",
        )
        adjacency.append(neighbors)
    edges = sorted(
        (i, j) for i in range(PORT_COUNT) for j in adjacency[i] if i < j
    )
    require(
        len(edges) == SEAM_COUNT,
        "MOVE_SIMPLEX",
        f"expected thirty seams, got {len(edges)}",
    )
    return vertices, adjacency, edges


def adjacency_automorphisms(adjacency: Sequence[Sequence[int]]) -> list[tuple[int, ...]]:
    """Backtracking enumeration of every adjacency-preserving port permutation."""

    n = len(adjacency)
    neighbor_sets = [frozenset(row) for row in adjacency]
    assignment = [-1] * n
    used = [False] * n
    found: list[tuple[int, ...]] = []

    def consistent(vertex: int, image: int) -> bool:
        for earlier in range(vertex):
            if (earlier in neighbor_sets[vertex]) != (
                assignment[earlier] in neighbor_sets[image]
            ):
                return False
        return True

    def recurse(vertex: int) -> None:
        if vertex == n:
            found.append(tuple(assignment))
            return
        for image in range(n):
            if used[image] or not consistent(vertex, image):
                continue
            assignment[vertex] = image
            used[image] = True
            recurse(vertex + 1)
            used[image] = False
            assignment[vertex] = -1

    recurse(0)
    return found


# ---------------------------------------------------------------------------
# Exact spectral certificate for the induced port walk W = A/5
# ---------------------------------------------------------------------------


def adjacency_int_matrix(adjacency: Sequence[Sequence[int]]) -> list[list[int]]:
    n = len(adjacency)
    sets = [set(row) for row in adjacency]
    return [[1 if j in sets[i] else 0 for j in range(n)] for i in range(n)]


def matmul_f5(a: Sequence[Sequence[F5]], b: Sequence[Sequence[F5]]) -> list[list[F5]]:
    n = len(a)
    result = []
    for i in range(n):
        row = []
        for j in range(n):
            entry = ZERO
            for k in range(n):
                entry = entry + a[i][k] * b[k][j]
            row.append(entry)
        result.append(row)
    return result


def shifted_f5(matrix: Sequence[Sequence[F5]], scalar: F5) -> list[list[F5]]:
    """The matrix plus scalar times the identity."""

    n = len(matrix)
    return [
        [matrix[i][j] + (scalar if i == j else ZERO) for j in range(n)]
        for i in range(n)
    ]


def is_zero_matrix(matrix: Sequence[Sequence[F5]]) -> bool:
    return all(entry.is_zero() for row in matrix for entry in row)


def f5_power(value: F5, exponent: int) -> F5:
    result = ONE
    for _ in range(exponent):
        result = result * value
    return result


def solve_f5_linear(matrix: Sequence[Sequence[F5]], rhs: Sequence[F5]) -> list[F5]:
    """Exact Gaussian elimination over Q(sqrt5)."""

    n = len(matrix)
    augmented = [list(matrix[i]) + [rhs[i]] for i in range(n)]
    for column in range(n):
        pivot = next(
            (row for row in range(column, n) if not augmented[row][column].is_zero()),
            None,
        )
        require(pivot is not None, "SPECTRAL", "singular trace system")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        inverse = augmented[column][column].inv()
        augmented[column] = [entry * inverse for entry in augmented[column]]
        for row in range(n):
            if row == column or augmented[row][column].is_zero():
                continue
            factor = augmented[row][column]
            augmented[row] = [
                entry - factor * lead
                for entry, lead in zip(augmented[row], augmented[column])
            ]
    return [augmented[i][n] for i in range(n)]


def spectral_certificate(adjacency_matrix: Sequence[Sequence[int]]) -> dict[str, Any]:
    """Exact spectrum of A and of W = A/5 in Q(sqrt5).

    The identity (A - 5I)(A^2 - 5I)(A + I) = 0 is verified entrywise; over
    Q(sqrt5) the middle factor splits as (A - sqrt5 I)(A + sqrt5 I), and each
    product with one linear factor removed is nonzero, so the polynomial is
    minimal and every root is realized. The multiplicities follow from the
    exact trace system on A^0 through A^3, whose Vandermonde matrix on the
    four distinct roots is invertible.
    """

    n = len(adjacency_matrix)
    require(n == PORT_COUNT, "SPECTRAL", "the adjacency matrix is not twelve by twelve")
    for i in range(n):
        require(
            adjacency_matrix[i][i] == 0,
            "SPECTRAL",
            "the adjacency matrix has a nonzero diagonal entry",
        )
        for j in range(n):
            require(
                adjacency_matrix[i][j] == adjacency_matrix[j][i],
                "SPECTRAL",
                "the adjacency matrix is not symmetric",
            )
    matrix = [[F5(entry) for entry in row] for row in adjacency_matrix]

    eigenvalues = (FIVE, SQRT5, -SQRT5, -ONE)
    for left in range(4):
        for right in range(left + 1, 4):
            require(
                not (eigenvalues[left] - eigenvalues[right]).is_zero(),
                "SPECTRAL",
                "the declared adjacency eigenvalues are not distinct",
            )
    factors = [shifted_f5(matrix, -eigenvalue) for eigenvalue in eigenvalues]

    full_product = factors[0]
    for factor in factors[1:]:
        full_product = matmul_f5(full_product, factor)
    require(
        is_zero_matrix(full_product),
        "SPECTRAL",
        "the minimal polynomial identity (A - 5I)(A^2 - 5I)(A + I) = 0 fails",
    )
    for dropped in range(4):
        partial = None
        for index, factor in enumerate(factors):
            if index == dropped:
                continue
            partial = factor if partial is None else matmul_f5(partial, factor)
        require(
            partial is not None and not is_zero_matrix(partial),
            "SPECTRAL",
            "a declared eigenvalue is not realized by the adjacency matrix",
        )

    powers = [matrix]
    for _ in range(2):
        powers.append(matmul_f5(powers[-1], matrix))
    traces_f5 = [F5(n)] + [
        sum((power[i][i] for i in range(n)), ZERO) for power in powers
    ]
    trace_integers = []
    for value in traces_f5:
        require(
            value.b == 0 and value.a.denominator == 1,
            "SPECTRAL",
            "a trace of a power of A is not an integer",
        )
        trace_integers.append(int(value.a))

    vandermonde = [
        [f5_power(eigenvalue, exponent) for eigenvalue in eigenvalues]
        for exponent in range(4)
    ]
    multiplicities_f5 = solve_f5_linear(vandermonde, traces_f5)
    multiplicities = []
    for value in multiplicities_f5:
        require(
            value.b == 0 and value.a.denominator == 1 and value.a > 0,
            "SPECTRAL",
            "an eigenvalue multiplicity is not a positive integer",
        )
        multiplicities.append(int(value.a))
    require(
        multiplicities == [1, 3, 3, 5],
        "SPECTRAL",
        f"the adjacency multiplicities drifted: {multiplicities}",
    )

    # Walk spectrum: divide by the degree. The second-largest modulus is
    # sqrt5/5: its square 1/5 exceeds (1/5)^2 = 1/25 and stays below 1.
    second = SQRT5 / FIVE
    second_squared = second * second
    require(
        (second_squared - F5(Fraction(1, 25))).is_positive()
        and (ONE - second_squared).is_positive(),
        "SPECTRAL",
        "the second-largest walk eigenvalue modulus ordering fails",
    )
    gap = ONE - second
    require(
        gap == F5(1, Fraction(-1, 5)) and gap.is_positive(),
        "SPECTRAL",
        "the exact spectral gap 1 - sqrt5/5 does not recompute",
    )

    return {
        "adjacency_eigenvalues": [
            {"value": "5", "multiplicity": 1},
            {"value": "sqrt(5)", "multiplicity": 3},
            {"value": "-sqrt(5)", "multiplicity": 3},
            {"value": "-1", "multiplicity": 5},
        ],
        "minimal_polynomial_identity": "(A - 5I)(A^2 - 5I)(A + I) = 0, verified entrywise over Q(sqrt5)",
        "minimal_polynomial_is_minimal": True,
        "trace_powers_of_A": trace_integers,
        "multiplicities_from_exact_trace_system": multiplicities,
        "walk_eigenvalues": ["1", "sqrt(5)/5", "-sqrt(5)/5", "-1/5"],
        "second_largest_modulus": "sqrt(5)/5",
        "spectral_gap": {
            "display": "1 - (1/5)*sqrt(5)",
            "rational_part": "1",
            "sqrt5_coefficient": "-1/5",
        },
    }


# ---------------------------------------------------------------------------
# Reachability, periodicity, and cycle search helpers
# ---------------------------------------------------------------------------


def reachability_closure(port_count: int, edge_subset: Sequence[Edge], start: int) -> set[int]:
    neighbors: dict[int, set[int]] = {index: set() for index in range(port_count)}
    for u, v in edge_subset:
        neighbors[u].add(v)
        neighbors[v].add(u)
    reached = {start}
    frontier = [start]
    while frontier:
        vertex = frontier.pop()
        for neighbor in sorted(neighbors[vertex]):
            if neighbor not in reached:
                reached.add(neighbor)
                frontier.append(neighbor)
    return reached


def parity_reachability(
    port_count: int, edge_subset: Sequence[Edge], start: int
) -> tuple[set[int], set[int]]:
    """Ports reachable from the start in an even respectively odd step count."""

    neighbors: dict[int, set[int]] = {index: set() for index in range(port_count)}
    for u, v in edge_subset:
        neighbors[u].add(v)
        neighbors[v].add(u)
    even = {start}
    odd: set[int] = set()
    while True:
        new_odd = {
            neighbor for vertex in even for neighbor in neighbors[vertex]
        } - odd
        new_even = {
            neighbor for vertex in odd for neighbor in neighbors[vertex]
        } - even
        if not new_odd and not new_even:
            return even, odd
        even |= new_even
        odd |= new_odd


def enumerate_triangles(adjacency: Sequence[Sequence[int]]) -> list[tuple[int, int, int]]:
    sets = [set(row) for row in adjacency]
    triangles = []
    for i in range(len(adjacency)):
        for j in adjacency[i]:
            if j <= i:
                continue
            for k in adjacency[j]:
                if k <= j or k not in sets[i]:
                    continue
                triangles.append((i, j, k))
    return triangles


def find_six_cycle(adjacency: Sequence[Sequence[int]]) -> list[int]:
    """The lexicographically first simple six-cycle of the carrier."""

    sets = [set(row) for row in adjacency]

    def extend(path: list[int]) -> list[int] | None:
        if len(path) == 6:
            return path if path[0] in sets[path[-1]] else None
        for candidate in adjacency[path[-1]]:
            if candidate in path or candidate < path[0]:
                continue
            found = extend(path + [candidate])
            if found is not None:
                return found
        return None

    for start in range(len(adjacency)):
        found = extend([start])
        if found is not None:
            return found
    raise CertificateError("CONTROL_PERIODIC", "no six-cycle exists in the carrier")


# ---------------------------------------------------------------------------
# Controls, each required to fail closed
# ---------------------------------------------------------------------------


def control_wrong_reference(
    edges: Sequence[Edge],
    weights: Mapping[Edge, Fraction] | None = None,
    horizon: int = FAIRNESS_HORIZON_STEPS,
) -> dict[str, Any]:
    """A zeroed seam weight starves that seam: fairness fails exactly."""

    starved = edges[0]
    if weights is None:
        weights = {
            edge: (Fraction(0) if edge == starved else Fraction(1, len(edges) - 1))
            for edge in edges
        }
    require(
        sum(weights.values()) == 1 and all(value >= 0 for value in weights.values()),
        "CONTROL_WRONG_REFERENCE",
        "the control weight vector is not a point of the move simplex",
    )
    starved_probability = weights[starved]
    require(
        starved_probability == 0,
        "CONTROL_WRONG_REFERENCE",
        "the wrong-reference control unexpectedly assigns the seam positive weight",
    )
    unvisited_at_horizon = (1 - starved_probability) ** horizon
    require(
        unvisited_at_horizon == 1,
        "CONTROL_WRONG_REFERENCE",
        "the starved seam has a nonunit unvisited probability",
    )
    return {
        "verdict": "fails_closed",
        "witness": {
            "starved_seam": [port_label(starved[0]), port_label(starved[1])],
            "seam_selection_probability": "0",
            "unvisited_probability_after_horizon": "1",
            "horizon_steps": horizon,
            "expected_first_visit": "unbounded",
            "conclusion": "a reference with a zero seam weight never repairs that seam, so bounded fairness fails",
        },
    }


def default_reducible_edges(
    edges: Sequence[Edge], adjacency: Sequence[Sequence[int]]
) -> tuple[list[Edge], dict[str, Any]]:
    """Two vertex-disjoint triangles as a disconnected move subset."""

    triangles = enumerate_triangles(adjacency)
    first = triangles[0]
    second = next(
        triangle
        for triangle in triangles
        if not set(triangle) & set(first)
    )
    subset = []
    for triangle in (first, second):
        a, b, c = triangle
        subset.extend(
            [tuple(sorted((a, b))), tuple(sorted((b, c))), tuple(sorted((a, c)))]
        )
    require(
        all(edge in set(edges) for edge in subset) and len(subset) == 6,
        "CONTROL_REDUCIBLE",
        "the disjoint-triangle move subset is not a seam subset",
    )
    witness = {
        "first_triangle": [port_label(v) for v in first],
        "second_triangle": [port_label(v) for v in second],
    }
    return subset, witness


def control_reducible(
    edges: Sequence[Edge],
    adjacency: Sequence[Sequence[int]],
    subset: Sequence[Edge] | None = None,
) -> dict[str, Any]:
    """A disconnected proper move subset fails the reachability cover."""

    if subset is None:
        subset, witness = default_reducible_edges(edges, adjacency)
    else:
        subset = list(subset)
        witness = {"first_triangle": None, "second_triangle": None}
    start = min(min(edge) for edge in subset)
    closure = reachability_closure(PORT_COUNT, subset, start)
    require(
        closure != set(range(PORT_COUNT)),
        "CONTROL_REDUCIBLE",
        "the reducible control unexpectedly covers every port",
    )
    witness.update(
        {
            "move_subset_seam_count": len(subset),
            "closure_start_port": port_label(start),
            "reachability_closure_ports": sorted(port_label(v) for v in closure),
            "ports_covered": len(closure),
            "ports_total": PORT_COUNT,
            "conclusion": "the reachability closure of the restricted move set does not cover the ports, so the restricted scheduler is reducible",
        }
    )
    return {"verdict": "fails_closed", "witness": witness}


def control_periodic(
    adjacency: Sequence[Sequence[int]],
    cycle: Sequence[int] | None = None,
) -> dict[str, Any]:
    """A bipartite move subset has period two on its component."""

    cycle = list(cycle) if cycle is not None else find_six_cycle(adjacency)
    cycle_edges = [
        tuple(sorted((cycle[index], cycle[(index + 1) % len(cycle)])))
        for index in range(len(cycle))
    ]
    sets = [set(row) for row in adjacency]
    for u, v in cycle_edges:
        require(
            v in sets[u],
            "CONTROL_PERIODIC",
            "the declared cycle is not a seam cycle of the carrier",
        )
    even, odd = parity_reachability(PORT_COUNT, cycle_edges, cycle[0])
    require(
        even.isdisjoint(odd) and odd,
        "CONTROL_PERIODIC",
        "the bipartite control unexpectedly reaches a port at both parities",
    )
    return {
        "verdict": "fails_closed",
        "witness": {
            "cycle_ports": [port_label(v) for v in cycle],
            "cycle_length": len(cycle),
            "even_class_ports": sorted(port_label(v) for v in even),
            "odd_class_ports": sorted(port_label(v) for v in odd),
            "period": 2,
            "conclusion": "the even and odd reachability classes of the restricted walk are disjoint, so the restricted scheduler has period two",
        },
    }


def det3_int(matrix: Sequence[Sequence[int]]) -> int:
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def rotation_port_permutation(
    vertices: Sequence[tuple[F5, F5, F5]],
    adjacency: Sequence[Sequence[int]],
    matrix: Sequence[Sequence[int]],
) -> tuple[int, ...]:
    """The port permutation induced by an exact integer rotation matrix."""

    require(
        det3_int(matrix) == 1,
        "CONTROL_PRESENTATION",
        "the candidate deck element is not a rotation",
    )
    index = {vertex: position for position, vertex in enumerate(vertices)}
    permutation = []
    for vertex in vertices:
        image = tuple(
            F5(matrix[row][0]) * vertex[0]
            + F5(matrix[row][1]) * vertex[1]
            + F5(matrix[row][2]) * vertex[2]
            for row in range(3)
        )
        require(
            image in index,
            "CONTROL_PRESENTATION",
            "the rotation does not preserve the port set",
        )
        permutation.append(index[image])
    require(
        sorted(permutation) == list(range(len(vertices))),
        "CONTROL_PRESENTATION",
        "the rotation does not induce a port bijection",
    )
    sets = [set(row) for row in adjacency]
    for i in range(len(vertices)):
        for j in adjacency[i]:
            require(
                permutation[j] in sets[permutation[i]],
                "CONTROL_PRESENTATION",
                "the rotation does not preserve the seam relation",
            )
    return tuple(permutation)


def control_presentation(
    vertices: Sequence[tuple[F5, F5, F5]],
    adjacency: Sequence[Sequence[int]],
    edges: Sequence[Edge],
    rotations: Sequence[Sequence[Sequence[int]]] | None = None,
) -> dict[str, Any]:
    """The lowest-index-first scheduler is presentation dependent.

    One nonidentity deck rotation relabels the ports; the deterministic
    scheduler then names a different physical seam as its first move while
    the relabeled seam multiset, the degree sequence, and the seam count all
    agree with the reference presentation.
    """

    candidates = rotations if rotations is not None else ROTATION_CANDIDATES
    first_original = edges[0]
    for matrix in candidates:
        permutation = rotation_port_permutation(vertices, adjacency, matrix)

        def relabeled_name(edge: Edge) -> Edge:
            return tuple(sorted((permutation[edge[0]], permutation[edge[1]])))

        first_relabeled_physical = min(edges, key=relabeled_name)
        if first_relabeled_physical == first_original:
            continue
        require(
            permutation != tuple(range(len(vertices))),
            "CONTROL_PRESENTATION",
            "a differing first-move record arose from the identity relabeling",
        )
        require(
            sorted(relabeled_name(edge) for edge in edges) == list(edges),
            "CONTROL_PRESENTATION",
            "the relabeled seam multiset does not agree with the reference presentation",
        )
        return {
            "verdict": "fails_closed",
            "witness": {
                "deck_rotation_matrix": [list(row) for row in matrix],
                "induced_port_permutation": list(permutation),
                "rotation_is_identity": False,
                "first_move_reference_presentation": [
                    port_label(first_original[0]),
                    port_label(first_original[1]),
                ],
                "first_move_relabeled_presentation_physical_seam": [
                    port_label(first_relabeled_physical[0]),
                    port_label(first_relabeled_physical[1]),
                ],
                "quotient_visible_agreement": {
                    "seam_count": len(edges),
                    "degree_sequence": [len(row) for row in adjacency],
                    "relabeled_seam_multiset_equals_reference": True,
                },
                "uniform_kernel_pullback_invariant": True,
                "conclusion": "the deterministic lowest-index-first scheduler names different physical first seams in two presentations whose quotient-visible data agree, so presentation invariance rejects it; the uniform kernel is invariant under the same relabeling",
            },
        }
    raise CertificateError(
        "CONTROL_PRESENTATION",
        "the deterministic scheduler produced identical first-move records under every candidate deck rotation",
    )


# ---------------------------------------------------------------------------
# Certificate payload
# ---------------------------------------------------------------------------


def certificate_payload() -> dict[str, Any]:
    vertices, adjacency, edges = build_move_simplex()

    # --- A3 selection on the move simplex -----------------------------------
    kernel_probability = Fraction(1, SEAM_COUNT)
    require(
        kernel_probability * SEAM_COUNT == 1 and kernel_probability > 0,
        "A3_PROJECTION",
        "the uniform kernel is not a positive point of the move simplex",
    )
    automorphisms = adjacency_automorphisms(adjacency)
    require(
        len(automorphisms) == 120
        and tuple(range(PORT_COUNT)) in set(automorphisms),
        "DECK_INVARIANCE",
        f"the carrier automorphism group has order {len(automorphisms)}, expected 120",
    )
    seam_orbit = {
        tuple(sorted((permutation[edges[0][0]], permutation[edges[0][1]])))
        for permutation in automorphisms
    }
    require(
        seam_orbit == set(edges),
        "DECK_INVARIANCE",
        "the carrier automorphism group is not transitive on the seams",
    )

    # --- Induced port walk: stationarity and per-seam traversal --------------
    adjacency_matrix = adjacency_int_matrix(adjacency)
    stationary = Fraction(1, PORT_COUNT)
    step = Fraction(1, DEGREE)
    for v in range(PORT_COUNT):
        column_mass = sum(
            stationary * step
            for u in range(PORT_COUNT)
            if adjacency_matrix[u][v] == 1
        )
        require(
            column_mass == stationary,
            "WALK_STATIONARITY",
            "the uniform distribution is not stationary for W = A/5",
        )
    for u, v in edges:
        traversal = stationary * step + stationary * step
        require(
            traversal == kernel_probability,
            "WALK_STATIONARITY",
            "a stationary per-step seam traversal probability differs from the kernel",
        )

    # --- Proved dynamics ------------------------------------------------------
    closure = reachability_closure(PORT_COUNT, edges, 0)
    require(
        closure == set(range(PORT_COUNT)),
        "IRREDUCIBILITY",
        "the reference carrier walk is not irreducible",
    )
    triangle = enumerate_triangles(adjacency)[0]
    spectral = spectral_certificate(adjacency_matrix)

    unvisited_at_horizon = Fraction(SEAM_COUNT - 1, SEAM_COUNT) ** FAIRNESS_HORIZON_STEPS
    expected_first_visit = 1 / kernel_probability
    require(
        expected_first_visit == SEAM_COUNT and 0 < unvisited_at_horizon < 1,
        "FAIRNESS",
        "the exact fairness bounds do not recompute",
    )

    # --- Controls, each required to fail closed -------------------------------
    controls = {
        "wrong_reference": control_wrong_reference(edges),
        "reducible": control_reducible(edges, adjacency),
        "periodic": control_periodic(adjacency),
        "presentation_dependent": control_presentation(vertices, adjacency, edges),
    }

    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "description": (
            "Executable certificate for the bounded A3 scheduler campaign. The "
            "move simplex is the thirty-seam icosahedral reference carrier; the "
            "A3 information projection selects the uniform kernel 1/30 per "
            "seam in closed form; the induced port walk W = A/5 has exact "
            "positive support, irreducibility, aperiodicity, contraction with "
            "spectral gap 1 - sqrt5/5, and exact fairness and excursion "
            "bounds; four scheduler controls fail closed. Fairness and "
            "liveness for a general implementation model stay a "
            "conditional_open_interface premise, per the bounded exit of the "
            "issue."
        ),
        "arithmetic": "exact Fraction and Q(sqrt5) arithmetic; no floating point in any proof step",
        "move_simplex": {
            "ports": PORT_COUNT,
            "seams": SEAM_COUNT,
            "degree": DEGREE,
            "coordinate_model": "cyclic permutations of (0, +/-1, +/-phi) over Q(sqrt5)",
            "port_coordinates": {
                port_label(i): [component.text() for component in vertex]
                for i, vertex in enumerate(vertices)
            },
            "seams_by_port_pair": [
                [port_label(u), port_label(v)] for u, v in edges
            ],
            "source_definition": (
                "each move is one undirected seam repair between adjacent ports "
                "of the reference carrier; the transition observable is the "
                "source-defined record of which seam repairs next"
            ),
        },
        "a3_selection": {
            "reference": "uniform_on_moves",
            "feasible_set": "the full move simplex; the complete constraint grammar contributes normalization only",
            "positive_weights": True,
            "information_projection": (
                "the minimizer of the relative entropy to the reference over "
                "the full simplex is the reference itself: KKT stationarity "
                "log(p_e/r_e) + 1 + lambda = 0 forces p proportional to r, "
                "normalization fixes the constant, and strict convexity gives "
                "uniqueness; this is the same argument as the equal-state-"
                "weights certificate"
            ),
            "selected_kernel_probability_per_seam": str(kernel_probability),
            "kernel_sums_to_one": True,
            "deck_invariance": {
                "carrier_automorphism_group_order": len(automorphisms),
                "seam_orbit_size": len(seam_orbit),
                "edge_transitive": True,
                "conclusion": (
                    "the automorphism group acts transitively on the seams, so "
                    "a deck-invariant positive weight vector is constant and "
                    "normalizes to the same uniform kernel"
                ),
            },
        },
        "induced_port_walk": {
            "transition_matrix": "W = A/5, the uniform neighbor walk on the twelve ports",
            "stationary_distribution": "uniform 1/12, verified exactly",
            "per_step_seam_traversal_probability": str(kernel_probability),
            "matches_selected_kernel": True,
        },
        "proved_dynamics": {
            "positive_support": {
                "minimum_move_probability": str(kernel_probability),
                "all_moves_positive": True,
            },
            "irreducibility": {
                "reachability_closure_size": len(closure),
                "connected": True,
            },
            "aperiodicity": {
                "triangle_witness": [port_label(v) for v in triangle],
                "closed_walk_lengths": [2, 3],
                "period": 1,
            },
            "contraction": spectral,
            "fairness": {
                "per_step_seam_visit_probability": str(kernel_probability),
                "unvisited_probability_after_t_steps": "(29/30)^t",
                "expected_first_visit_steps": str(expected_first_visit),
                "horizon_steps": FAIRNESS_HORIZON_STEPS,
                "unvisited_probability_at_horizon": str(unvisited_at_horizon),
                "unvisited_probability_numerator_digits": len(
                    str(unvisited_at_horizon.numerator)
                ),
                "unvisited_probability_denominator_digits": len(
                    str(unvisited_at_horizon.denominator)
                ),
            },
        },
        "controls": controls,
        "verdict": {
            "kernel": "uniform_1_over_30",
            "kernel_status": "A3_selected_exact",
            "dynamics_tier": "exact_named_realization",
            "dynamics_scope": "the reference carrier walk; no continuum or implementation claim",
            "general_model_fairness_liveness": "conditional_open_interface",
            "general_model_note": (
                "fairness and liveness of an arbitrary implementation remain "
                "an implementation premise; the positive theorem covers the "
                "named reference realization only"
            ),
            "controls": {
                name: record["verdict"] for name, record in controls.items()
            },
        },
        "verifier_command": (
            "python3 code/a5_closure/a3_scheduler_kernel_certificate.py verify "
            "--manifest code/a5_closure/manifests/a3_scheduler_kernel_reference.json"
        ),
    }


def build_manifest() -> dict[str, Any]:
    body = certificate_payload()
    manifest = dict(body)
    manifest["manifest_sha256"] = "sha256:" + sha256_json(body)
    return manifest


def verify_manifest(manifest: Mapping[str, Any]) -> None:
    require(
        manifest.get("schema") == SCHEMA,
        "SCHEMA",
        f"expected {SCHEMA}",
    )
    body = {key: value for key, value in manifest.items() if key != "manifest_sha256"}
    require(
        manifest.get("manifest_sha256") == "sha256:" + sha256_json(body),
        "MANIFEST_HASH",
        "the manifest self-hash does not recompute",
    )
    require(
        body == certificate_payload(),
        "MANIFEST_MISMATCH",
        "the manifest is stale, malformed, or tampered",
    )


def default_manifest_path() -> Path:
    return MODULE_DIR / "manifests" / "a3_scheduler_kernel_reference.json"


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    certify = sub.add_parser("certify")
    certify.add_argument("--output", type=Path, default=default_manifest_path())
    verify = sub.add_parser("verify")
    verify.add_argument("--manifest", type=Path, default=default_manifest_path())
    sub.add_parser("all")
    args = parser.parse_args(argv)
    try:
        if args.command == "certify":
            manifest = build_manifest()
            write_json(args.output, manifest)
            print(
                json.dumps(
                    {
                        "status": "PASS",
                        "manifest": str(args.output),
                        "manifest_sha256": manifest["manifest_sha256"],
                    },
                    indent=2,
                )
            )
        elif args.command == "verify":
            verify_manifest(load_json(args.manifest))
            print(json.dumps({"status": "PASS", "manifest": str(args.manifest)}, indent=2))
        else:
            output = default_manifest_path()
            manifest = build_manifest()
            write_json(output, manifest)
            verify_manifest(load_json(output))
            print(
                json.dumps(
                    {
                        "status": "PASS",
                        "manifest": str(output),
                        "manifest_sha256": manifest["manifest_sha256"],
                    },
                    indent=2,
                )
            )
    except CertificateError as exc:
        print(
            json.dumps({"status": "FAIL", "code": exc.code, "message": exc.message}, indent=2)
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
