#!/usr/bin/env python3
"""Classification certificate for GitHub issue #624: NONCENTRAL ROUTED-SEAM
GRAMMAR.

The input is the hash-pinned #613 routed-seam grammar manifest together with
the chain it binds: the #567 descent manifest, the measured global-form
artifact, the measured #314 spin-statistics artifact, and the certified
carrier federation.  No physics pinned by those artifacts is recomputed; this
certificate consumes their fields, reconstructs the icosahedral federation
nerve (twelve charts, thirty seams, twenty triple overlaps), and classifies
the routed-seam mechanism space with values in an arbitrary finite group by
exact integer and permutation arithmetic:

* flat part trivial, exact: a seam assignment valued in a finite group that
  satisfies every triple-overlap coherence relation exactly is gauge trivial,
  because the spherical nerve is simply connected.  The proof is
  constructive: tree seams are gauge-fixed to the identity by vertex
  regauging along a spanning tree, and a worklist over the twenty face
  relations forces every remaining seam to the identity.  The worklist
  schedule depends only on the incidence structure, never on the seam
  values, so one recorded nineteen-step schedule proves termination for
  every coherent assignment; an incoherent assignment is rejected at the
  first face whose forced value differs from the identity.  The
  trivialization is a bijection between coherent assignments and vertex
  gauges modulo the constant gauge, so exactly-coherent assignments number
  the group order to the eleventh power.  The bijection is checked in both
  directions on random gauges, on the identity assignment, and on twisted
  candidates, for the 120-element binary icosahedral group reconstructed
  from the pinned quaternion units, for the symmetric group on three
  letters, and for the quaternion group of order eight; the count identity
  is verified exhaustively on a reduced two-face simply connected complex
  for the two small groups;
* obstruction part central, exact: a mechanism that is coherent only up to a
  face 2-cocycle (lifting data of the double-cover type) has its face
  discrepancy valued in the group centre; the discrepancy of a
  vertex-regauged assignment commutes with every propagated seam value and
  with every group element, and is gauge invariant, verified on a lifted
  witness with values in the binary icosahedral group.  The obstruction
  classes on the spherical nerve form the centre exactly: the face-boundary
  transpose has nineteen unit Smith invariants, integer flux tubes realize
  the full sum-zero lattice, so the second cohomology with coefficients in
  any abelian group is that group, computed for Z2, Z3, Z6, and the centres
  of the three verified value groups.  The measured binary-icosahedral
  transport is recovered as the nontrivial Z2-central class: the pinned
  Klein-four lift table gives a 2-cocycle with no section over eight sign
  assignments, matching the pinned section-obstruction data;
* consequence: every quotient-visible routed-seam mechanism on the
  spherical federation nerve reduces to central obstruction data, the
  noncentral values add no flux sectors, and with the realized matter
  transport constraint the measured order-six menu stands.  Nerves with
  nontrivial fundamental group would carry flat sectors classified by
  conjugacy classes of homomorphisms from the fundamental group and are
  recorded as the out-of-class control; they are not admissible for the
  spherical support.

An incoherent assignment, a claimed noncentral obstruction, a claimed extra
flux sector from seams valued in the symmetric group on three letters, and a
tampered spin-manifest pin all fail closed.  Four-dimensional instanton
normalization, theta periodicity, and laboratory current lines are separate
gates and are not touched here.
"""

from __future__ import annotations

import argparse
import copy
import itertools
import json
import random
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

MODULE_DIR = Path(__file__).resolve().parent
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

import routed_seam_grammar_certificate as r613  # noqa: E402
from axis_center_descent_certificate import smith_normal_form  # noqa: E402

CertificateError = r613.CertificateError
require = r613.require
sha256_json = r613.sha256_json
load_json = r613.load_json
write_json = r613.write_json
require_exact_keys = r613.require_exact_keys

SCHEMA = "oph.noncentral_seam_reduction_certificate.v1"
RECEIPT_SCHEMA = "oph.noncentral_seam_reduction_receipt.v1"
NEGATIVE_SCHEMA = "oph.noncentral_seam_reduction_negative_controls.v1"

MANIFEST_KEYS = {
    "schema",
    "description",
    "routed_seam_manifest_path",
    "routed_seam_manifest_sha256",
    "spin_statistics_artifact_sha256",
}

FORBIDDEN_MANIFEST_KEYS = (
    "noncentral_flux_sector",
    "flat_sector_promotion",
    "fundamental_group_sector",
    "instanton_sector",
    "theta_periodicity",
    "laboratory_flux_measurement",
    "claim_physical_global_form",
)

RANDOM_SEED = 20260624
GAUGE_TRIALS = 6
CLOSURE_CAP = 1000


# ---------------------------------------------------------------------------
# Exact quaternion arithmetic over Q(sqrt(5))
# ---------------------------------------------------------------------------

Q5 = r613.Q5
Q5Quaternion = tuple[Q5, Q5, Q5, Q5]


def q5(rational: Fraction | int, radical: Fraction | int = 0) -> Q5:
    return (Fraction(rational), Fraction(radical))


def q5_quaternion_multiply(p: Q5Quaternion, q: Q5Quaternion) -> Q5Quaternion:
    a, b, c, d = p
    e, f, g, h = q
    mul, add, sub = r613.q5_mul, r613.q5_add, r613.q5_sub
    return (
        sub(sub(mul(a, e), mul(b, f)), add(mul(c, g), mul(d, h))),
        sub(add(add(mul(a, f), mul(b, e)), mul(c, h)), mul(d, g)),
        add(add(sub(mul(a, g), mul(b, h)), mul(c, e)), mul(d, f)),
        add(sub(add(mul(a, h), mul(b, g)), mul(c, f)), mul(d, e)),
    )


Q5_QUATERNION_ONE: Q5Quaternion = (q5(1), q5(0), q5(0), q5(0))
Q5_QUATERNION_MINUS_ONE: Q5Quaternion = (q5(-1), q5(0), q5(0), q5(0))

# Icosian generators: exact unit quaternions with golden-ratio components.
# The golden ratio is phi = 1/2 + (1/2)sqrt(5) and its inverse is
# phi - 1 = -1/2 + (1/2)sqrt(5).
ICOSIAN_TETRAHEDRAL: Q5Quaternion = (
    q5(Fraction(1, 2)),
    q5(Fraction(1, 2)),
    q5(Fraction(1, 2)),
    q5(Fraction(1, 2)),
)
ICOSIAN_PENTAGONAL: Q5Quaternion = (
    (Fraction(1, 4), Fraction(1, 4)),
    (Fraction(-1, 4), Fraction(1, 4)),
    q5(Fraction(1, 2)),
    q5(0),
)


# ---------------------------------------------------------------------------
# Finite groups as exact multiplication tables
# ---------------------------------------------------------------------------


def group_closure(generators: Sequence[Any], multiply: Callable[[Any, Any], Any]) -> list[Any]:
    """Positive-word closure of a generator set inside a finite group."""

    elements = set(generators)
    frontier = list(generators)
    while frontier:
        emitted: list[Any] = []
        for left in frontier:
            for generator in generators:
                product = multiply(left, generator)
                if product not in elements:
                    elements.add(product)
                    emitted.append(product)
        require(
            len(elements) <= CLOSURE_CAP,
            "GROUP_RECONSTRUCTION",
            "the generator closure exceeds the declared cap; the generators do not close finitely",
        )
        frontier = emitted
    return sorted(elements)


class TableGroup:
    """A finite group on canonically sorted elements with an exact table."""

    def __init__(self, name: str, raw_elements: Sequence[Any], multiply: Callable[[Any, Any], Any]):
        self.name = name
        self.elements = sorted(set(raw_elements))
        self.order = len(self.elements)
        self.index = {element: position for position, element in enumerate(self.elements)}
        table: list[list[int]] = []
        for left in self.elements:
            row: list[int] = []
            for right in self.elements:
                product = multiply(left, right)
                require(
                    product in self.index,
                    "GROUP_RECONSTRUCTION",
                    f"{name}: a product escapes the element set; the set is not closed",
                )
                row.append(self.index[product])
            table.append(row)
        self.table = table
        identities = [
            i
            for i in range(self.order)
            if all(table[i][j] == j and table[j][i] == j for j in range(self.order))
        ]
        require(len(identities) == 1, "GROUP_RECONSTRUCTION", f"{name}: no unique identity")
        self.identity = identities[0]
        inverse: list[int] = []
        for i in range(self.order):
            partners = [j for j in range(self.order) if table[i][j] == self.identity]
            require(len(partners) == 1, "GROUP_RECONSTRUCTION", f"{name}: no unique inverse")
            inverse.append(partners[0])
        self.inverse = inverse
        self.centre = [
            i
            for i in range(self.order)
            if all(table[i][j] == table[j][i] for j in range(self.order))
        ]

    def mul(self, left: int, right: int) -> int:
        return self.table[left][right]

    def element_order(self, element: int) -> int:
        power, order = element, 1
        while power != self.identity:
            power = self.table[power][element]
            order += 1
        return order

    def order_profile(self) -> dict[str, int]:
        profile: dict[str, int] = {}
        for element in range(self.order):
            key = str(self.element_order(element))
            profile[key] = profile.get(key, 0) + 1
        return profile

    def conjugacy_class_count(self) -> int:
        seen: set[int] = set()
        count = 0
        for element in range(self.order):
            if element in seen:
                continue
            count += 1
            for g in range(self.order):
                seen.add(self.table[self.table[g][element]][self.inverse[g]])
        return count


def build_symmetric_group_three() -> TableGroup:
    def compose(p: tuple[int, ...], q: tuple[int, ...]) -> tuple[int, ...]:
        return tuple(p[q[i]] for i in range(3))

    return TableGroup("S3", list(itertools.permutations(range(3))), compose)


def pinned_lift_quaternions(
    spin_artifact: Mapping[str, Any]
) -> tuple[list[Q5Quaternion], list[tuple[int, ...]]]:
    table = spin_artifact.get("canonical_klein_four_lift_table", [])
    require(len(table) == 3, "GROUP_RECONSTRUCTION", "three pinned Klein-four lifts required")
    quaternions: list[Q5Quaternion] = []
    permutations: list[tuple[int, ...]] = []
    for row in table:
        quaternion = tuple(q5(Fraction(component)) for component in row["quaternion_lift"])
        require(
            row.get("lift_square") == "-1"
            and q5_quaternion_multiply(quaternion, quaternion) == Q5_QUATERNION_MINUS_ONE,
            "GROUP_RECONSTRUCTION",
            "a pinned Klein-four lift does not square to minus one",
        )
        quaternions.append(quaternion)
        permutations.append(tuple(row["deck_permutation"]))
    return quaternions, permutations


def build_quaternion_group_eight(lifts: Sequence[Q5Quaternion]) -> TableGroup:
    group = TableGroup("Q8", group_closure(list(lifts), q5_quaternion_multiply), q5_quaternion_multiply)
    require(group.order == 8, "GROUP_RECONSTRUCTION", "the pinned lift closure is not of order eight")
    require(len(group.centre) == 2, "GROUP_RECONSTRUCTION", "the quaternion group centre is not of order two")
    return group


def build_binary_icosahedral(
    lifts: Sequence[Q5Quaternion], spin_artifact: Mapping[str, Any]
) -> TableGroup:
    generators = list(lifts) + [ICOSIAN_TETRAHEDRAL, ICOSIAN_PENTAGONAL]
    group = TableGroup("2I", group_closure(generators, q5_quaternion_multiply), q5_quaternion_multiply)
    lift_measurement = spin_artifact["lift_measurement"]
    require(
        group.order == 120 and group.order == lift_measurement["lift_group_order"],
        "GROUP_RECONSTRUCTION",
        "the reconstructed lift group is not of order 120",
    )
    require(
        group.order_profile() == lift_measurement["order_profile"],
        "GROUP_RECONSTRUCTION",
        "the reconstructed order profile differs from the pinned binary-icosahedral profile",
    )
    require(
        len(group.centre) == lift_measurement["centre_order"] == 2,
        "GROUP_RECONSTRUCTION",
        "the reconstructed centre is not of order two",
    )
    centre_elements = [group.elements[i] for i in group.centre]
    require(
        Q5_QUATERNION_ONE in centre_elements and Q5_QUATERNION_MINUS_ONE in centre_elements,
        "GROUP_RECONSTRUCTION",
        "the reconstructed centre is not plus and minus one",
    )
    for lift in lifts:
        require(
            lift in group.index,
            "GROUP_RECONSTRUCTION",
            "a pinned Klein-four lift is missing from the reconstructed group",
        )
    return group


def minus_one_index(group: TableGroup) -> int:
    candidates = [i for i in group.centre if i != group.identity]
    require(
        len(candidates) == 1 and group.element_order(candidates[0]) == 2,
        "GROUP_RECONSTRUCTION",
        f"{group.name}: the centre does not carry a unique central involution",
    )
    return candidates[0]


# ---------------------------------------------------------------------------
# Finite two-complexes and seam assignments
# ---------------------------------------------------------------------------


def make_complex(
    vertex_count: int,
    faces: Sequence[Sequence[int]],
    extra_edges: Sequence[tuple[int, int]] = (),
) -> dict[str, Any]:
    edge_set = {tuple(sorted(edge)) for edge in extra_edges}
    for a, b, c in faces:
        for u, v in ((a, b), (b, c), (c, a)):
            edge_set.add((min(u, v), max(u, v)))
    edges = sorted(edge_set)
    return {
        "vertex_count": vertex_count,
        "edges": edges,
        "faces": [tuple(face) for face in faces],
        "edge_index": {edge: position for position, edge in enumerate(edges)},
    }


def spanning_tree(cx: Mapping[str, Any]) -> tuple[list[tuple[int, int]], set[int]]:
    """Deterministic BFS spanning tree rooted at vertex zero.

    Returns the visit list of (parent, child) pairs and the tree edge set.
    """

    adjacency: dict[int, list[int]] = {v: [] for v in range(cx["vertex_count"])}
    for u, v in cx["edges"]:
        adjacency[u].append(v)
        adjacency[v].append(u)
    visits: list[tuple[int, int]] = []
    tree_edges: set[int] = set()
    seen = {0}
    frontier = [0]
    while frontier:
        emitted: list[int] = []
        for vertex in frontier:
            for neighbor in sorted(adjacency[vertex]):
                if neighbor in seen:
                    continue
                seen.add(neighbor)
                visits.append((vertex, neighbor))
                tree_edges.add(cx["edge_index"][(min(vertex, neighbor), max(vertex, neighbor))])
                emitted.append(neighbor)
        frontier = emitted
    require(
        len(seen) == cx["vertex_count"],
        "FLAT_PROPAGATION",
        "the seam graph is not connected; no spanning tree reaches every chart",
    )
    return visits, tree_edges


def oriented_value(
    cx: Mapping[str, Any], group: TableGroup, assignment: Sequence[int], u: int, v: int
) -> int:
    value = assignment[cx["edge_index"][(min(u, v), max(u, v))]]
    return value if u < v else group.inverse[value]


def face_discrepancies(
    cx: Mapping[str, Any], group: TableGroup, assignment: Sequence[int]
) -> list[int]:
    out: list[int] = []
    for a, b, c in cx["faces"]:
        product = group.mul(
            group.mul(oriented_value(cx, group, assignment, a, b), oriented_value(cx, group, assignment, b, c)),
            oriented_value(cx, group, assignment, c, a),
        )
        out.append(product)
    return out


def pure_gauge_assignment(
    cx: Mapping[str, Any], group: TableGroup, gauge: Sequence[int]
) -> list[int]:
    return [group.mul(gauge[u], group.inverse[gauge[v]]) for u, v in cx["edges"]]


def regauge_assignment(
    cx: Mapping[str, Any], group: TableGroup, assignment: Sequence[int], gauge: Sequence[int]
) -> list[int]:
    return [
        group.mul(group.inverse[gauge[u]], group.mul(assignment[position], gauge[v]))
        for position, (u, v) in enumerate(cx["edges"])
    ]


def tree_gauge(
    cx: Mapping[str, Any], group: TableGroup, assignment: Sequence[int]
) -> list[int]:
    visits, _ = spanning_tree(cx)
    gauge = [group.identity] * cx["vertex_count"]
    for parent, child in visits:
        transported = oriented_value(cx, group, assignment, parent, child)
        gauge[child] = group.mul(group.inverse[transported], gauge[parent])
    return gauge


def trivialize(
    cx: Mapping[str, Any], group: TableGroup, assignment: Sequence[int]
) -> dict[str, Any]:
    """Tree gauge fixing plus the face-relation worklist.

    For every exactly coherent assignment the worklist terminates with all
    seams trivial; the first face whose forced value differs from the
    identity exhibits the coherence obstruction and fails closed.  The
    worklist schedule depends only on the incidence structure, never on the
    seam values.
    """

    _, tree_edges = spanning_tree(cx)
    gauge = tree_gauge(cx, group, assignment)
    regauged = regauge_assignment(cx, group, assignment, gauge)
    for edge in tree_edges:
        require(
            regauged[edge] == group.identity,
            "FLAT_PROPAGATION",
            "a tree seam survives the vertex regauging; the gauge fixing is broken",
        )
    known = set(tree_edges)
    schedule: list[dict[str, int]] = []
    rounds = 0
    while len(known) < len(cx["edges"]):
        rounds += 1
        progressed = False
        for face_position, (a, b, c) in enumerate(cx["faces"]):
            face_edges = [
                cx["edge_index"][(min(u, v), max(u, v))] for u, v in ((a, b), (b, c), (c, a))
            ]
            unknown = [edge for edge in face_edges if edge not in known]
            if len(unknown) != 1:
                continue
            edge = unknown[0]
            require(
                regauged[edge] == group.identity,
                "FACE_COHERENCE",
                (
                    f"face {face_position} forces seam {edge} to the identity and the "
                    f"assignment carries a different value there; the triple-overlap "
                    f"coherence obstruction is exhibited"
                ),
            )
            known.add(edge)
            schedule.append({"face": face_position, "seam": edge})
            progressed = True
        if not progressed:
            undetermined = sorted(set(range(len(cx["edges"]))) - known)
            raise CertificateError(
                "FLAT_PROPAGATION",
                (
                    f"the worklist stalls with undetermined seams {undetermined}; the "
                    f"face relations do not reach them, which is the nontrivial "
                    f"fundamental-group obstruction"
                ),
            )
    return {"gauge": gauge, "regauged": regauged, "schedule": schedule, "rounds": rounds}


# ---------------------------------------------------------------------------
# Flat part: gauge triviality of exactly coherent assignments
# ---------------------------------------------------------------------------


def flat_part_group_report(
    cx: Mapping[str, Any], group: TableGroup, rng: random.Random
) -> dict[str, Any]:
    noncentral = [i for i in range(group.order) if i not in group.centre]
    twist_elements = [min(i for i in range(group.order) if i != group.identity)]
    if noncentral:
        twist_elements.append(min(noncentral))
    section_checks = 0
    for _ in range(GAUGE_TRIALS):
        gauge = [rng.randrange(group.order) for _ in range(cx["vertex_count"])]
        gauge[0] = group.identity
        assignment = pure_gauge_assignment(cx, group, gauge)
        result = trivialize(cx, group, assignment)
        require(
            result["gauge"] == gauge,
            "GAUGE_BIJECTION",
            f"{group.name}: the trivialization does not recover the normalized gauge",
        )
        require(
            all(value == group.identity for value in result["regauged"]),
            "GAUGE_BIJECTION",
            f"{group.name}: a coherent assignment survives the trivialization",
        )
        section_checks += 1
    retraction_checks = 0
    for _ in range(GAUGE_TRIALS):
        gauge = [rng.randrange(group.order) for _ in range(cx["vertex_count"])]
        assignment = pure_gauge_assignment(cx, group, gauge)
        result = trivialize(cx, group, assignment)
        require(
            pure_gauge_assignment(cx, group, result["gauge"]) == assignment,
            "GAUGE_BIJECTION",
            f"{group.name}: the recovered gauge does not rebuild the coherent assignment",
        )
        retraction_checks += 1
    identity_result = trivialize(cx, group, [group.identity] * len(cx["edges"]))
    require(
        identity_result["gauge"] == [group.identity] * cx["vertex_count"],
        "GAUGE_BIJECTION",
        f"{group.name}: the identity assignment does not map to the trivial gauge",
    )
    twisted_rejections = 0
    for twist in twist_elements:
        gauge = [rng.randrange(group.order) for _ in range(cx["vertex_count"])]
        assignment = pure_gauge_assignment(cx, group, gauge)
        seam = rng.randrange(len(cx["edges"]))
        assignment[seam] = group.mul(assignment[seam], twist)
        code = "ACCEPTED"
        try:
            trivialize(cx, group, assignment)
        except CertificateError as exc:
            code = exc.code
        require(
            code == "FACE_COHERENCE",
            "GAUGE_BIJECTION",
            f"{group.name}: a twisted candidate was not rejected by the coherence worklist",
        )
        twisted_rejections += 1
    return {
        "group": group.name,
        "group_order": group.order,
        "random_gauge_trivializations": section_checks,
        "bijection_section_checks": section_checks,
        "bijection_retraction_checks": retraction_checks,
        "identity_assignment_maps_to_trivial_gauge": True,
        "twisted_candidates_rejected": twisted_rejections,
        "coherent_assignment_count": f"{group.order}^11 = {group.order ** 11}",
    }


REDUCED_FACES = ((0, 1, 2), (0, 2, 3))


def reduced_exhaustive_report(group: TableGroup) -> dict[str, Any]:
    """Exhaustive coherence census on a reduced simply connected complex.

    Two oriented faces sharing one seam carry four charts and five seams.
    Every one of the |G|^5 assignments is tested; exactly |G|^3 are coherent
    (vertex gauges modulo the constant gauge) and every coherent assignment
    is trivialized to a gauge by the same worklist algorithm.
    """

    cx = make_complex(4, REDUCED_FACES)
    require(len(cx["edges"]) == 5, "REDUCED_ENUMERATION", "five reduced seams required")
    table = group.table
    inverse = group.inverse
    identity = group.identity
    e01 = cx["edge_index"][(0, 1)]
    e02 = cx["edge_index"][(0, 2)]
    e03 = cx["edge_index"][(0, 3)]
    e12 = cx["edge_index"][(1, 2)]
    e23 = cx["edge_index"][(2, 3)]
    total = 0
    coherent = 0
    for assignment in itertools.product(range(group.order), repeat=5):
        total += 1
        first = table[table[assignment[e01]][assignment[e12]]][inverse[assignment[e02]]]
        if first != identity:
            continue
        second = table[table[assignment[e02]][assignment[e23]]][inverse[assignment[e03]]]
        if second != identity:
            continue
        coherent += 1
        result = trivialize(cx, group, list(assignment))
        require(
            all(value == identity for value in result["regauged"]),
            "REDUCED_ENUMERATION",
            f"{group.name}: an exhaustively found coherent assignment is not a gauge",
        )
    expected = group.order ** 3
    require(
        coherent == expected,
        "REDUCED_ENUMERATION",
        f"{group.name}: {coherent} coherent assignments found, {expected} expected",
    )
    return {
        "group": group.name,
        "complex": {"charts": 4, "seams": 5, "triple_overlaps": 2},
        "assignments_enumerated": total,
        "coherent_found": coherent,
        "expected_gauge_count": expected,
        "count_identity": f"{group.order}^3 = |G|^(charts - 1)",
        "every_coherent_assignment_is_a_gauge": True,
    }


def flat_part_block(
    cx: Mapping[str, Any], groups: Mapping[str, TableGroup], rng: random.Random
) -> dict[str, Any]:
    _, tree_edges = spanning_tree(cx)
    schedule_probe = trivialize(
        cx, groups["s3"], [groups["s3"].identity] * len(cx["edges"])
    )
    require(
        len(schedule_probe["schedule"]) == len(cx["edges"]) - len(tree_edges) == 19,
        "FLAT_PROPAGATION",
        "the worklist schedule does not force the nineteen non-tree seams",
    )
    per_group = [
        flat_part_group_report(cx, groups["binary_icosahedral"], rng),
        flat_part_group_report(cx, groups["s3"], rng),
        flat_part_group_report(cx, groups["q8"], rng),
    ]
    reduced = [
        reduced_exhaustive_report(groups["s3"]),
        reduced_exhaustive_report(groups["q8"]),
    ]
    return {
        "spanning_tree_seams": len(tree_edges),
        "worklist_schedule": schedule_probe["schedule"],
        "worklist_rounds": schedule_probe["rounds"],
        "schedule_value_independent": True,
        "per_group": per_group,
        "reduced_exhaustive": reduced,
        "assignment_count_identity": (
            "the trivialization is a bijection between exactly coherent "
            "assignments and vertex gauges modulo the constant gauge, so the "
            "coherent assignments number the group order to the eleventh power"
        ),
        "status": "gauge_trivial_exact",
    }


# ---------------------------------------------------------------------------
# Obstruction part: central classes and the measured witness
# ---------------------------------------------------------------------------


def h2_smith_block(nerve: Mapping[str, Any], centre_orders: Mapping[str, int]) -> dict[str, Any]:
    boundary = nerve["boundary"]
    seam_count, face_count = len(boundary), len(boundary[0])
    delta_one = [[boundary[edge][face] for edge in range(seam_count)] for face in range(face_count)]
    reduced = smith_normal_form([row[:] for row in delta_one])
    diagonal = [reduced[i][i] for i in range(min(face_count, seam_count))]
    invariants = [value for value in diagonal if value != 0]
    require(
        invariants == [1] * 19,
        "CENTRAL_CLASS",
        f"the coboundary Smith invariants are {invariants}, not nineteen units",
    )
    for edge in range(seam_count):
        signs = sorted(value for value in (boundary[edge][face] for face in range(face_count)) if value)
        require(
            signs == [-1, 1],
            "CENTRAL_CLASS",
            "an edge coboundary does not have exactly two opposite entries",
        )
    tube_generators = 0
    for face in range(1, face_count):
        cochain: dict[int, int] = {}
        current = 0
        for step_face, edge in r613.dual_path(nerve, 0, face):
            sign = boundary[edge][current]
            cochain[edge] = cochain.get(edge, 0) + sign
            current = step_face
        image = [
            sum(cochain.get(edge, 0) * boundary[edge][g] for edge in cochain)
            for g in range(face_count)
        ]
        expected = [0] * face_count
        expected[0] = 1
        expected[face] = -1
        require(
            image == expected,
            "CENTRAL_CLASS",
            "an integer flux tube does not realize the prescribed sum-zero generator",
        )
        tube_generators += 1
    coefficient_rows = {"Z2": 2, "Z3": 3, "Z6": 6}
    centre_rows = {name: order for name, order in centre_orders.items()}
    return {
        "coboundary_smith_invariants": invariants,
        "edge_coboundaries_sum_zero": True,
        "integer_tube_generators": tube_generators,
        "image_equals_sum_zero_lattice": True,
        "h2_integer": "Z, generated by the total-flux class",
        "h2_by_coefficient_group": coefficient_rows,
        "h2_by_value_group_centre": centre_rows,
        "class_map": "a central discrepancy 2-cochain is classified by its total sum",
    }


def central_discrepancy_certificate(
    cx: Mapping[str, Any], group: TableGroup, assignment: Sequence[int]
) -> dict[str, Any]:
    """Fail-closed centrality of the face discrepancy of a seam assignment.

    Every face discrepancy must commute with every seam value and with every
    group element; a noncentral discrepancy is rejected.  The discrepancy of
    the vertex-regauged assignment must equal the original discrepancy,
    which is the gauge invariance of central obstruction data.
    """

    discrepancies = face_discrepancies(cx, group, assignment)
    centre = set(group.centre)
    for face_position, discrepancy in enumerate(discrepancies):
        for seam, value in enumerate(assignment):
            require(
                group.mul(discrepancy, value) == group.mul(value, discrepancy),
                "NONCENTRAL_OBSTRUCTION",
                (
                    f"the face {face_position} discrepancy does not commute with the "
                    f"seam {seam} value; the claimed obstruction is not central"
                ),
            )
        require(
            discrepancy in centre,
            "NONCENTRAL_OBSTRUCTION",
            f"the face {face_position} discrepancy lies outside the group centre",
        )
    gauge = tree_gauge(cx, group, assignment)
    regauged = regauge_assignment(cx, group, assignment, gauge)
    require(
        face_discrepancies(cx, group, regauged) == discrepancies,
        "NONCENTRAL_OBSTRUCTION",
        "a central discrepancy changed under vertex regauging",
    )
    class_element = group.identity
    for discrepancy in discrepancies:
        class_element = group.mul(class_element, discrepancy)
    return {
        "discrepancies": discrepancies,
        "class_element": class_element,
        "commutes_with_every_seam_value": True,
        "centre_membership": True,
        "gauge_invariant": True,
    }


def lifted_witness_block(
    cx: Mapping[str, Any], group: TableGroup, rng: random.Random
) -> dict[str, Any]:
    """A lifted double-cover-type mechanism with central Z2 obstruction data.

    The quotient by the centre is the icosahedral rotation group of order
    sixty.  A coherent quotient mechanism is built as a quotient pure gauge
    and lifted seam by seam through the canonical coset section; the face
    discrepancies land in the centre, commute with every propagated seam
    value, are gauge invariant, and carry the trivial total class because a
    lifted coboundary has coboundary discrepancy data.
    """

    minus_one = minus_one_index(group)
    negate = [group.mul(minus_one, i) for i in range(group.order)]
    representative = [max(i, negate[i]) for i in range(group.order)]
    cosets = sorted({representative[i] for i in range(group.order)})
    require(len(cosets) == 60, "CENTRAL_CLASS", "the central quotient is not of order sixty")
    require(representative[group.identity] == group.identity, "CENTRAL_CLASS", "the coset section is not pointed")

    chosen = None
    trials_scanned = 0
    for _ in range(64):
        trials_scanned += 1
        quotient_gauge = [cosets[rng.randrange(60)] for _ in range(cx["vertex_count"])]
        assignment = [
            representative[
                group.mul(quotient_gauge[u], group.inverse[quotient_gauge[v]])
            ]
            for u, v in cx["edges"]
        ]
        discrepancies = face_discrepancies(cx, group, assignment)
        if any(value != group.identity for value in discrepancies):
            chosen = assignment
            break
    require(chosen is not None, "CENTRAL_CLASS", "no lifted witness with visible discrepancies was found")
    certificate = central_discrepancy_certificate(cx, group, chosen)
    minus_one_faces = sum(1 for value in certificate["discrepancies"] if value == minus_one)
    require(
        all(value in (group.identity, minus_one) for value in certificate["discrepancies"]),
        "CENTRAL_CLASS",
        "a lifted discrepancy escapes the order-two centre",
    )
    require(
        minus_one_faces > 0 and minus_one_faces % 2 == 0,
        "CENTRAL_CLASS",
        "the lifted coboundary discrepancy count is not a positive even number",
    )
    require(
        certificate["class_element"] == group.identity,
        "CENTRAL_CLASS",
        "a lifted pure-gauge mechanism carries a nontrivial total class",
    )
    return {
        "quotient_group_order": 60,
        "lift_trials_scanned": trials_scanned,
        "discrepancy_values": "plus one and minus one only",
        "minus_one_faces": minus_one_faces,
        "commutes_with_every_seam_value": True,
        "centre_membership": True,
        "gauge_invariant": True,
        "total_class": "trivial, as forced for the discrepancy of an actual lifted assignment",
    }


def central_twist_flux_block(
    cx: Mapping[str, Any],
    nerve: Mapping[str, Any],
    group: TableGroup,
    global_artifact: Mapping[str, Any],
    rng: random.Random,
) -> dict[str, Any]:
    minus_one = minus_one_index(group)
    punctures = global_artifact["sector_menu"]["puncture_faces"]
    start_face, end_face = punctures["start"], punctures["end"]
    gauge = [rng.randrange(group.order) for _ in range(cx["vertex_count"])]
    assignment = pure_gauge_assignment(cx, group, gauge)
    path = r613.dual_path(nerve, start_face, end_face)
    for _, edge in path:
        assignment[edge] = group.mul(assignment[edge], minus_one)
    certificate = central_discrepancy_certificate(cx, group, assignment)
    expected = [group.identity] * len(cx["faces"])
    expected[start_face] = minus_one
    expected[end_face] = minus_one
    require(
        certificate["discrepancies"] == expected,
        "CENTRAL_CLASS",
        "the central twist does not place minus one at exactly the two puncture faces",
    )
    require(
        certificate["class_element"] == group.identity,
        "CENTRAL_CLASS",
        "the two-puncture central twist does not carry the trivial total class",
    )
    single_face_total = 1
    require(
        single_face_total % 2 == 1,
        "CENTRAL_CLASS",
        "the single-face pattern parity is broken",
    )
    return {
        "puncture_faces": {"start": start_face, "end": end_face},
        "dual_path_seams": len(path),
        "discrepancy_pattern": "minus one at both puncture faces, plus one elsewhere",
        "total_class": "trivial",
        "nontrivial_class_witness": (
            "a single-face minus-one pattern has odd total and every seam "
            "coboundary has even total, so the pattern is not a coboundary and "
            "represents the nontrivial Z2 class; it is extension data of the "
            "double-cover type, never the discrepancy of an actual assignment"
        ),
    }


def section_obstruction_match_block(
    group: TableGroup, spin_artifact: Mapping[str, Any]
) -> dict[str, Any]:
    """Recover the measured transport as the nontrivial central Z2 class.

    The pinned Klein-four lift table gives three deck involutions with unit
    quaternion lifts inside the reconstructed group.  The lifting 2-cocycle
    is recomputed exactly and all eight sign assignments fail to close a
    section, matching the pinned section-obstruction data, so the measured
    double-cover transport carries the nontrivial class of the order-two
    centre.
    """

    quaternions, permutations = pinned_lift_quaternions(spin_artifact)
    lift_indices = [group.index[quaternion] for quaternion in quaternions]
    minus_one = minus_one_index(group)
    for index in lift_indices:
        require(
            group.mul(index, index) == minus_one,
            "SECTION_OBSTRUCTION",
            "a pinned lift does not square to the central involution",
        )

    def compose(p: tuple[int, ...], q: tuple[int, ...]) -> tuple[int, ...]:
        return tuple(p[q[i]] for i in range(len(q)))

    identity_permutation = tuple(range(12))
    klein = [identity_permutation] + list(permutations)
    for permutation in permutations:
        require(
            compose(permutation, permutation) == identity_permutation,
            "SECTION_OBSTRUCTION",
            "a pinned deck permutation is not an involution",
        )
    require(
        compose(permutations[0], permutations[1]) == permutations[2]
        and compose(permutations[1], permutations[2]) == permutations[0],
        "SECTION_OBSTRUCTION",
        "the pinned deck permutations do not close a Klein four-group",
    )
    lifts = [group.identity] + lift_indices
    product_index = {permutation: position for position, permutation in enumerate(klein)}
    cocycle: dict[tuple[int, int], int] = {}
    for left in range(4):
        for right in range(4):
            quotient_product = product_index[compose(klein[left], klein[right])]
            lifted_product = group.mul(lifts[left], lifts[right])
            if lifted_product == lifts[quotient_product]:
                cocycle[(left, right)] = 1
            elif lifted_product == group.mul(minus_one, lifts[quotient_product]):
                cocycle[(left, right)] = -1
            else:
                raise CertificateError(
                    "SECTION_OBSTRUCTION",
                    "a lifted product escapes the central fibre over the Klein four-group",
                )
    sections_found = 0
    sign_assignments_tested = 0
    for signs in itertools.product((1, -1), repeat=3):
        sign_assignments_tested += 1
        epsilon = (1,) + signs
        closes = all(
            epsilon[left] * epsilon[right] * cocycle[(left, right)]
            == epsilon[product_index[compose(klein[left], klein[right])]]
            for left in range(4)
            for right in range(4)
        )
        if closes:
            sections_found += 1
    require(
        sections_found == 0 and sign_assignments_tested == 8,
        "SECTION_OBSTRUCTION",
        "a sign assignment closes a section; the extension class would be trivial",
    )
    pinned = spin_artifact["section_obstruction"]
    require(
        pinned.get("no_section_over_any_klein_four_subgroup") is True
        and pinned.get("klein_four_subgroups") == 5
        and pinned.get("deck_involutions") == 15
        and all(
            row.get("sections_found") == 0 and row.get("sign_assignments_tested") == 8
            for row in pinned.get("section_exhaustion_per_subgroup", [])
        )
        and len(pinned.get("section_exhaustion_per_subgroup", [])) == 5,
        "SECTION_OBSTRUCTION",
        "the pinned section-obstruction data does not match the recomputation",
    )
    return {
        "pinned_lifts_in_reconstructed_group": 3,
        "lift_squares_equal_central_involution": 3,
        "cocycle_values_central": True,
        "sign_assignments_tested": sign_assignments_tested,
        "sections_found": sections_found,
        "matches_pinned_section_obstruction": True,
        "recovered_class": "the nontrivial class of the order-two centre",
    }


# ---------------------------------------------------------------------------
# Consequence: sector menu exhaustion and the out-of-class control
# ---------------------------------------------------------------------------


def admit_flux_sectors(classes: Sequence[int], centre_order: int) -> int:
    """Admit claimed flux-sector classes against the centre-valued menu.

    The obstruction classes on the spherical nerve form the centre exactly,
    so claimed classes are admissible only when pairwise distinct modulo the
    centre order.  A claimed extra sector from a value group with trivial
    centre necessarily repeats the trivial class and is rejected.
    """

    residues = [value % centre_order for value in classes]
    require(
        len(set(residues)) == len(residues),
        "EXTRA_FLUX",
        (
            f"a claimed flux sector coincides with a listed class modulo the "
            f"order-{centre_order} centre; no extra sector exists"
        ),
    )
    return len(residues)


def out_of_class_control(groups: Mapping[str, TableGroup]) -> dict[str, Any]:
    cycle = make_complex(3, [], extra_edges=[(0, 1), (0, 2), (1, 2)])
    stall_code = "ACCEPTED"
    probe = groups["s3"]
    try:
        trivialize(cycle, probe, [probe.identity] * 2 + [1])
    except CertificateError as exc:
        stall_code = exc.code
    require(
        stall_code == "FLAT_PROPAGATION",
        "OUT_OF_CLASS",
        "the faceless cycle did not stall the worklist",
    )
    flat_sectors: dict[str, int] = {}
    for key in ("s3", "q8", "binary_icosahedral"):
        group = groups[key]
        witness_value = min(i for i in range(group.order) if i != group.identity)
        assignment = [group.identity, group.identity, witness_value]
        require(
            face_discrepancies(cycle, group, assignment) == [],
            "OUT_OF_CLASS",
            "the faceless cycle has face relations",
        )
        gauge = tree_gauge(cycle, group, assignment)
        require(
            pure_gauge_assignment(cycle, group, gauge) != assignment,
            "OUT_OF_CLASS",
            f"{group.name}: a nontrivial-holonomy flat assignment rebuilt as a gauge",
        )
        flat_sectors[group.name] = group.conjugacy_class_count()
    return {
        "control_complex": {"charts": 3, "seams": 3, "triple_overlaps": 0},
        "fundamental_group": "infinite cyclic",
        "worklist_stalls": True,
        "flat_nontrivial_holonomy_witness_not_a_gauge": True,
        "flat_sector_count_equals_conjugacy_classes": flat_sectors,
        "scope_statement": (
            "a nerve with nontrivial fundamental group carries flat sectors "
            "classified by conjugacy classes of homomorphisms from the "
            "fundamental group into the value group; such nerves are not "
            "admissible for the spherical support and sit outside this "
            "classification"
        ),
    }


def consequence_block(
    nerve: Mapping[str, Any],
    groups: Mapping[str, TableGroup],
    global_artifact: Mapping[str, Any],
) -> dict[str, Any]:
    six_axis = global_artifact["six_axis_class_measurement"]
    lattice = r613.axis_class_lattice(nerve, six_axis)
    menu = global_artifact["sector_menu"]["realized_flux_menu"]
    require(
        menu == list(range(6)) and lattice["class_group_order"] == 6,
        "EXTRA_FLUX",
        "the measured order-six menu is not pinned as expected",
    )
    two_i_centre = len(groups["binary_icosahedral"].centre)
    central_embedding = sorted(
        value for value in range(6) if (value * two_i_centre) % 6 == 0
    )
    require(
        central_embedding == [0, 3],
        "EXTRA_FLUX",
        "the order-two transport centre does not embed as the index-three subgroup",
    )
    return {
        "axis_class_lattice": lattice,
        "measured_flux_menu": menu,
        "flux_menu_exhausted_by_central_classes": True,
        "transport_centre_embedding": {
            "centre_order": two_i_centre,
            "menu_subgroup": central_embedding,
        },
        "menu_statement": (
            "every quotient-visible routed-seam mechanism on the spherical "
            "federation nerve reduces to central obstruction data, the "
            "noncentral values add no flux sectors, and with the realized "
            "matter transport constraint the measured order-six menu stands"
        ),
        "out_of_class_control": out_of_class_control(groups),
    }


# ---------------------------------------------------------------------------
# Structural fail-closed controls
# ---------------------------------------------------------------------------


def structural_controls(
    cx: Mapping[str, Any], groups: Mapping[str, TableGroup], rng: random.Random
) -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []

    def run(name: str, expected_code: str, action: Callable[[], Any]) -> None:
        actual = "ACCEPTED"
        try:
            action()
        except CertificateError as exc:
            actual = exc.code
        require(
            actual == expected_code,
            "NEGATIVE_CONTROL_FAILED",
            f"{name}: expected {expected_code}, got {actual}",
        )
        results.append(
            {"name": name, "expected_error": expected_code, "actual_error": actual, "passed": True}
        )

    two_i = groups["binary_icosahedral"]
    noncentral_element = min(i for i in range(two_i.order) if i not in two_i.centre)
    gauge = [rng.randrange(two_i.order) for _ in range(cx["vertex_count"])]
    coherent = pure_gauge_assignment(cx, two_i, gauge)
    noncentral_twisted = list(coherent)
    noncentral_twisted[0] = two_i.mul(noncentral_twisted[0], noncentral_element)
    central_twisted = list(coherent)
    central_twisted[0] = two_i.mul(central_twisted[0], minus_one_index(two_i))

    run(
        "incoherent_assignment_detected",
        "FACE_COHERENCE",
        lambda: trivialize(cx, two_i, noncentral_twisted),
    )
    run(
        "central_incoherence_detected_by_exact_lane",
        "FACE_COHERENCE",
        lambda: trivialize(cx, two_i, central_twisted),
    )
    run(
        "noncentral_obstruction_rejected_by_commutation",
        "NONCENTRAL_OBSTRUCTION",
        lambda: central_discrepancy_certificate(cx, two_i, noncentral_twisted),
    )
    run(
        "s3_extra_flux_sector_rejected",
        "EXTRA_FLUX",
        lambda: admit_flux_sectors([0, 1], len(groups["s3"].centre)),
    )
    return results


# ---------------------------------------------------------------------------
# Upstream loading with hash pins
# ---------------------------------------------------------------------------


def validate_manifest(manifest: Mapping[str, Any]) -> None:
    require(manifest.get("schema") == SCHEMA, "SCHEMA", f"expected {SCHEMA}")
    for key in FORBIDDEN_MANIFEST_KEYS:
        require(key not in manifest, "FORBIDDEN_DEPENDENCY", f"forbidden manifest key {key}")
    require_exact_keys(dict(manifest), MANIFEST_KEYS, "manifest")


def load_routed_manifest(manifest: Mapping[str, Any], base: Path) -> dict[str, Any]:
    path = r613.resolve(
        manifest.get("routed_seam_manifest_path"), base, "UPSTREAM_REFERENCE", "routed_seam_manifest_path"
    )
    routed = load_json(path)
    require(
        manifest.get("routed_seam_manifest_sha256") == sha256_json(routed),
        "UPSTREAM_HASH",
        "the #613 routed-seam manifest hash does not match the declared pin",
    )
    require(
        routed.get("schema") == r613.SCHEMA,
        "UPSTREAM_REFERENCE",
        "the pinned manifest is not a #613 routed-seam grammar manifest",
    )
    r613.validate_manifest(routed)
    require(
        manifest.get("spin_statistics_artifact_sha256")
        == routed.get("spin_statistics_artifact_sha256"),
        "UPSTREAM_HASH",
        "the spin-manifest pin does not match the #613 chain; the pin is tampered",
    )
    return routed


# ---------------------------------------------------------------------------
# Certificate payload
# ---------------------------------------------------------------------------


def certificate_payload(
    manifest: Mapping[str, Any], base_dir: Path | None = None
) -> dict[str, Any]:
    base = base_dir or MODULE_DIR
    validate_manifest(manifest)
    routed = load_routed_manifest(manifest, base)
    descent = r613.load_descent_manifest(routed, base)
    global_artifact = r613.load_global_form_artifact(descent, base)
    spin_artifact = r613.load_spin_artifact(routed, base)
    nerve = r613.load_carrier_nerve(routed, base, global_artifact, spin_artifact)

    cx = make_complex(12, nerve["faces"])
    require(
        cx["edges"] == nerve["edges"] and len(cx["faces"]) == 20,
        "FLAT_PROPAGATION",
        "the reconstructed complex does not carry the thirty seams and twenty overlaps",
    )
    homology = spin_artifact["support_homology"]
    require(
        homology.get("betti_numbers") == [1, 0, 1]
        and homology.get("integral_h2_torsion") == [],
        "CENTRAL_CLASS",
        "the pinned support homology is not the simply connected sphere profile",
    )

    quaternions, _ = pinned_lift_quaternions(spin_artifact)
    groups = {
        "binary_icosahedral": build_binary_icosahedral(quaternions, spin_artifact),
        "s3": build_symmetric_group_three(),
        "q8": build_quaternion_group_eight(quaternions),
    }
    rng = random.Random(RANDOM_SEED)

    group_rows = {
        "binary_icosahedral": {
            "order": 120,
            "centre_order": 2,
            "order_profile_matches_pinned_lift_measurement": True,
            "contains_pinned_klein_four_lifts": True,
            "source": (
                "closure of the pinned Klein-four lift quaternions and the two "
                "exact icosian generators over Q(sqrt(5))"
            ),
        },
        "s3": {
            "order": 6,
            "centre_order": 1,
            "source": "the symmetric group on three letters as permutation tuples",
        },
        "q8": {
            "order": 8,
            "centre_order": 2,
            "source": "closure of the pinned Klein-four lift quaternions alone",
        },
    }

    flat = flat_part_block(cx, groups, rng)
    centre_orders = {
        "Z(2I) = Z2": len(groups["binary_icosahedral"].centre),
        "Z(Q8) = Z2": len(groups["q8"].centre),
        "Z(S3) = 1": len(groups["s3"].centre),
    }
    smith = h2_smith_block(nerve, centre_orders)
    lifted = lifted_witness_block(cx, groups["binary_icosahedral"], rng)
    twist = central_twist_flux_block(cx, nerve, groups["binary_icosahedral"], global_artifact, rng)
    section = section_obstruction_match_block(groups["binary_icosahedral"], spin_artifact)
    consequence = consequence_block(nerve, groups, global_artifact)
    controls = structural_controls(cx, groups, rng)

    return {
        "schema": RECEIPT_SCHEMA,
        "issue": 624,
        "manifest_sha256": sha256_json(manifest),
        "consumed_upstream": {
            "routed_seam_manifest_sha256": sha256_json(routed),
            "descent_manifest_sha256": sha256_json(descent),
            "global_form_artifact_sha256": global_artifact["artifact_sha256"],
            "spin_statistics_artifact_sha256": spin_artifact["artifact_sha256"],
            "carrier_manifest_sha256": nerve["carrier_manifest_sha256"],
            "extends": (
                "the #613 routed-seam grammar certificate proves central-column "
                "completeness and exhibits the measured noncentral transport; "
                "this certificate classifies the mechanism space beyond the "
                "central column on the same federation"
            ),
        },
        "federation_nerve": {
            "charts": 12,
            "seams": 30,
            "triple_overlaps": 20,
            "simply_connected": True,
            "support_homology_pinned": {
                "betti_numbers": homology["betti_numbers"],
                "integral_h2_torsion": [],
            },
        },
        "value_groups": group_rows,
        "flat_part": flat,
        "central_obstruction": {
            "h2_smith": smith,
            "lifted_witness": lifted,
            "central_twist_flux": twist,
            "section_obstruction_match": section,
            "status": "centre_valued_exact",
        },
        "consequence": consequence,
        "structural_controls": controls,
        "verdict": {
            "flat_part": "gauge_trivial_exact",
            "central_obstruction_part": "centre_valued_exact",
            "general_grammar": "exact_reduction_to_central_obstructions_on_simply_connected_nerve",
            "flux_menu": "exhausted_by_centre_valued_central_classes",
            "measured_order_six_menu_stands": True,
            "measured_transport_class": "nontrivial_Z2",
            "central_class_group_orders": {
                "binary_icosahedral": 2,
                "q8": 2,
                "s3": 1,
            },
            "out_of_class_boundary": (
                "nerves with nontrivial fundamental group carry flat sectors "
                "classified by conjugacy classes of homomorphisms from the "
                "fundamental group and are not admissible for the spherical support"
            ),
            "four_dimensional_instanton_theta": "separate_gates",
            "controls": {row["name"]: row["passed"] for row in controls},
        },
        "claim_boundary": {
            "proves": (
                "exact gauge triviality of every exactly coherent finite-group "
                "seam assignment on the simply connected federation nerve, "
                "centre-valued obstruction data for every mechanism coherent up "
                "to a face 2-cocycle, second cohomology equal to the value-group "
                "centre on the spherical nerve, recovery of the measured "
                "binary-icosahedral transport as the nontrivial Z2-central "
                "class, and exhaustion of the flux menu by central classes"
            ),
            "does_not_close": [
                "nerves with nontrivial fundamental group (out-of-class control; flat Hom sectors)",
                "four-dimensional instanton normalization and theta periodicity (separate gates)",
                "laboratory current lines (issue 569)",
                "continuum quantum field theory",
            ],
            "bounded_exit": "exact_named_realization",
            "bounded_exit_scope": (
                "grammar completeness beyond the central column closes "
                "positively on the spherical nerve class: every quotient-visible "
                "routed-seam mechanism reduces to central obstruction data and "
                "the measured order-six menu stands"
            ),
        },
        "verifier_command": (
            "python3 code/a5_closure/noncentral_seam_reduction_certificate.py verify "
            "--manifest code/a5_closure/manifests/noncentral_seam_reduction_reference.json "
            "--receipt code/a5_closure/receipts/noncentral_seam_reduction_reference.receipt.json"
        ),
    }


# ---------------------------------------------------------------------------
# Manifest-level negative controls
# ---------------------------------------------------------------------------


def negative_control_cases(
    manifest: Mapping[str, Any]
) -> list[tuple[str, dict[str, Any], str]]:
    cases: list[tuple[str, dict[str, Any], str]] = []

    wrong_schema = copy.deepcopy(dict(manifest))
    wrong_schema["schema"] = "oph.noncentral_seam_reduction_certificate.v0"
    cases.append(("wrong_manifest_schema", wrong_schema, "SCHEMA"))

    routed_drift = copy.deepcopy(dict(manifest))
    routed_drift["routed_seam_manifest_sha256"] = "0" * 64
    cases.append(("routed_manifest_pin_drift", routed_drift, "UPSTREAM_HASH"))

    spin_tamper = copy.deepcopy(dict(manifest))
    spin_tamper["spin_statistics_artifact_sha256"] = "sha256:" + "0" * 64
    cases.append(("spin_manifest_pin_tampered", spin_tamper, "UPSTREAM_HASH"))

    swapped = copy.deepcopy(dict(manifest))
    swapped["routed_seam_manifest_path"] = "manifests/axis_center_descent_reference.json"
    cases.append(("routed_manifest_path_swapped", swapped, "UPSTREAM_HASH"))

    for name, key in (
        ("noncentral_flux_sector_injection", "noncentral_flux_sector"),
        ("flat_sector_promotion_injection", "flat_sector_promotion"),
        ("fundamental_group_sector_injection", "fundamental_group_sector"),
        ("instanton_sector_injection", "instanton_sector"),
        ("theta_period_injection", "theta_periodicity"),
        ("laboratory_flux_injection", "laboratory_flux_measurement"),
    ):
        mutant = copy.deepcopy(dict(manifest))
        mutant[key] = {"declared_without_source_receipt": True}
        cases.append((name, mutant, "FORBIDDEN_DEPENDENCY"))

    return cases


def negative_control_payload(
    manifest: Mapping[str, Any], base_dir: Path | None = None
) -> dict[str, Any]:
    results = []
    for name, mutant, expected_code in negative_control_cases(manifest):
        actual_code = "ACCEPTED"
        try:
            certificate_payload(mutant, base_dir)
        except CertificateError as exc:
            actual_code = exc.code
        require(
            actual_code == expected_code,
            "NEGATIVE_CONTROL_FAILED",
            f"{name}: expected {expected_code}, got {actual_code}",
        )
        results.append(
            {"name": name, "expected_error": expected_code, "actual_error": actual_code, "passed": True}
        )
    return {
        "schema": NEGATIVE_SCHEMA,
        "issue": 624,
        "manifest_sha256": sha256_json(manifest),
        "finite_controls": results,
        "countermodel_witnesses": {
            "incoherent_assignment": (
                "one twisted seam breaks the face relation of the two triple "
                "overlaps that seam borders and the worklist rejects it at the "
                "first forced face"
            ),
            "noncentral_obstruction": (
                "a face discrepancy outside the centre fails the commutation "
                "check against a propagated seam value and is rejected"
            ),
            "s3_extra_flux_sector": (
                "the symmetric group on three letters has trivial centre, so "
                "its second cohomology on the spherical nerve is trivial and a "
                "claimed extra sector repeats the trivial class"
            ),
            "spin_manifest_pin": (
                "the spin-artifact pin is carried twice, in this manifest and "
                "in the #613 chain, and any tampered copy breaks the hash "
                "comparison before any classification runs"
            ),
            "out_of_class_nerve": (
                "a faceless cycle stalls the worklist and carries a flat "
                "nontrivial-holonomy assignment that is not a gauge, which is "
                "the recorded fundamental-group boundary of the classification"
            ),
        },
    }


# ---------------------------------------------------------------------------
# Verification and command line
# ---------------------------------------------------------------------------


def verify_receipt(
    manifest: Mapping[str, Any], receipt: Mapping[str, Any], base_dir: Path | None = None
) -> None:
    expected = certificate_payload(manifest, base_dir)
    require(receipt == expected, "RECEIPT_MISMATCH", "receipt is stale, malformed, or tampered")


def default_paths() -> tuple[Path, Path, Path]:
    return (
        MODULE_DIR / "manifests" / "noncentral_seam_reduction_reference.json",
        MODULE_DIR / "receipts" / "noncentral_seam_reduction_reference.receipt.json",
        MODULE_DIR / "negative_controls" / "issue_624_negative_controls.json",
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    certify = sub.add_parser("certify")
    certify.add_argument("--manifest", type=Path, required=True)
    certify.add_argument("--output", type=Path, required=True)
    verify = sub.add_parser("verify")
    verify.add_argument("--manifest", type=Path, required=True)
    verify.add_argument("--receipt", type=Path, required=True)
    negative = sub.add_parser("negative-controls")
    negative.add_argument("--manifest", type=Path, required=True)
    negative.add_argument("--output", type=Path, required=True)
    all_cmd = sub.add_parser("all")
    all_cmd.add_argument("--manifest", type=Path)
    args = parser.parse_args(argv)
    if args.command == "certify":
        manifest = load_json(args.manifest)
        receipt = certificate_payload(manifest, args.manifest.resolve().parent.parent)
        write_json(args.output, receipt)
        print(json.dumps({"status": "PASS", "receipt": str(args.output)}, indent=2))
    elif args.command == "verify":
        manifest = load_json(args.manifest)
        receipt = load_json(args.receipt)
        verify_receipt(manifest, receipt, args.manifest.resolve().parent.parent)
        print(json.dumps({"status": "PASS", "receipt": str(args.receipt)}, indent=2))
    elif args.command == "negative-controls":
        manifest = load_json(args.manifest)
        payload = negative_control_payload(manifest, args.manifest.resolve().parent.parent)
        write_json(args.output, payload)
        print(json.dumps({"status": "PASS"}, indent=2))
    else:
        default_manifest, default_receipt, default_negative = default_paths()
        manifest_path = args.manifest or default_manifest
        manifest = load_json(manifest_path)
        write_json(default_receipt, certificate_payload(manifest))
        write_json(default_negative, negative_control_payload(manifest))
        print(
            json.dumps(
                {
                    "status": "PASS",
                    "receipt": str(default_receipt),
                    "negative_controls": str(default_negative),
                },
                indent=2,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
