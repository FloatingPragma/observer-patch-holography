#!/usr/bin/env python3
"""Completeness-scope certificate for GitHub issue #613: ROUTED-SEAM GRAMMAR.

The input is the hash-pinned #567 descent manifest together with the two
measured artifacts it binds: the global-form artifact (order-120 deck action
on the incidence-nerve federation, six-axis class group with Smith invariants
(1,1,1,1,1,6), order-six flux-sector menu, two-puncture witnesses) and the
spin-statistics artifact (binary-icosahedral double-cover transport with
group order 120, unique involution, centre of order two). No physics pinned
by those artifacts is recomputed; this certificate consumes their fields,
reconstructs the documented incidence structures, and proves three scope
statements by exact integer and rational arithmetic:

* central-column grammar completeness, exact: the routed-seam grammar with
  seam data valued in the central class group covers every seam-routable
  central flux mechanism.  The axis relation lattice recomputes to Smith
  invariants (1,1,1,1,1,6), the face-boundary matrix has nineteen unit Smith
  invariants, and constructive flux tubes from the start puncture to every
  other face generate the full sum-zero holonomy lattice, so every central
  2-cochain class on the nerve is realized by exactly one of the six
  measured sectors;
* closure of the composition laws: pairwise seam composition, triple-overlap
  coherence (the composed boundary vanishes identically and every
  face-holonomy assignment sums to zero), and routed-loop composition
  (loop concatenation is associative, and two routed loops differing by one
  face have holonomies differing by that face's boundary class, checked
  exhaustively over the twenty faces);
* the nonabelian scope boundary: the measured double-cover transport is a
  routed-seam mechanism on the same federation whose values lie in a
  120-element group with centre of order two, exhibited nonabelian by exact
  quaternion arithmetic on the pinned Klein-four lift table.  A
  quotient-visible routed-seam structure therefore exists outside the
  central-column grammar, and grammar exhaustiveness for all routed
  mechanisms stays a named open interface while the central column is
  complete.

Wrong routing, a deleted triple overlap, a claimed frame-dependent extra
sector under the Galois-conjugate frame swap, and a claimed seventh sector
class all fail closed.  Four-dimensional instanton normalization, theta
periodicity, and laboratory current lines are separate gates and are not
touched here.
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Sequence

MODULE_DIR = Path(__file__).resolve().parent
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

import echosahedral_selector_certificate as e565  # noqa: E402
from axis_center_descent_certificate import smith_normal_form  # noqa: E402

CertificateError = e565.CertificateError
require = e565.require
sha256_json = e565.sha256_json
load_json = e565.load_json
write_json = e565.write_json
require_exact_keys = e565.require_exact_keys

SCHEMA = "oph.routed_seam_grammar_certificate.v1"
RECEIPT_SCHEMA = "oph.routed_seam_grammar_receipt.v1"
NEGATIVE_SCHEMA = "oph.routed_seam_grammar_negative_controls.v1"
DESCENT_MANIFEST_SCHEMA = "oph.axis_center_descent_manifest.v4"
GLOBAL_FORM_ARTIFACT_SCHEMA = "oph.global_form_semantic_artifact.v1"
SPIN_ARTIFACT_SCHEMA = "oph.spin_statistics_semantic_artifact.v1"
CARRIER_MANIFEST_SCHEMA = "oph.echosahedral_selector_manifest.v1"

MANIFEST_KEYS = {
    "schema",
    "description",
    "descent_manifest_path",
    "descent_manifest_sha256",
    "spin_statistics_artifact_path",
    "spin_statistics_artifact_sha256",
    "carrier_manifest_path",
}

FORBIDDEN_MANIFEST_KEYS = (
    "general_grammar_closed",
    "nonabelian_grammar_completion",
    "instanton_sector",
    "theta_periodicity",
    "laboratory_flux_measurement",
    "monopole_dynamics",
    "claim_physical_global_form",
)


# ---------------------------------------------------------------------------
# Exact arithmetic in Q(sqrt(5)) and quaternions over it
# ---------------------------------------------------------------------------

Q5 = tuple[Fraction, Fraction]


def parse_q5(text: str) -> Q5:
    """Parse the pinned frame syntax: 'a', 'b*sqrt(5)', or 'a + b*sqrt(5)'."""

    rational, radical = Fraction(0), Fraction(0)
    for term in text.split(" + "):
        term = term.strip()
        if term.endswith("*sqrt(5)"):
            radical += Fraction(term[: -len("*sqrt(5)")])
        elif term == "sqrt(5)":
            radical += 1
        else:
            rational += Fraction(term)
    return (rational, radical)


def q5_add(x: Q5, y: Q5) -> Q5:
    return (x[0] + y[0], x[1] + y[1])


def q5_sub(x: Q5, y: Q5) -> Q5:
    return (x[0] - y[0], x[1] - y[1])


def q5_mul(x: Q5, y: Q5) -> Q5:
    return (x[0] * y[0] + 5 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def q5_neg(x: Q5) -> Q5:
    return (-x[0], -x[1])


def q5_conjugate(x: Q5) -> Q5:
    """The Galois conjugate sqrt(5) -> -sqrt(5)."""

    return (x[0], -x[1])


def q5_is_positive(x: Q5) -> bool:
    rational, radical = x
    if rational == 0 and radical == 0:
        return False
    if rational >= 0 and radical >= 0:
        return True
    if rational <= 0 and radical <= 0:
        return False
    if rational > 0:
        return rational * rational > 5 * radical * radical
    return 5 * radical * radical > rational * rational


def q5_det3(rows: Sequence[Sequence[Q5]]) -> Q5:
    (a, b, c), (d, e, f), (g, h, i) = rows
    term_one = q5_mul(a, q5_sub(q5_mul(e, i), q5_mul(f, h)))
    term_two = q5_mul(b, q5_sub(q5_mul(d, i), q5_mul(f, g)))
    term_three = q5_mul(c, q5_sub(q5_mul(d, h), q5_mul(e, g)))
    return q5_add(q5_sub(term_one, term_two), term_three)


Quaternion = tuple[Fraction, Fraction, Fraction, Fraction]


def quaternion_multiply(p: Quaternion, q: Quaternion) -> Quaternion:
    a, b, c, d = p
    e, f, g, h = q
    return (
        a * e - b * f - c * g - d * h,
        a * f + b * e + c * h - d * g,
        a * g - b * h + c * e + d * f,
        a * h + b * g - c * f + d * e,
    )


QUATERNION_MINUS_ONE: Quaternion = (Fraction(-1), Fraction(0), Fraction(0), Fraction(0))


# ---------------------------------------------------------------------------
# Upstream loading with hash pins
# ---------------------------------------------------------------------------


def validate_manifest(manifest: Mapping[str, Any]) -> None:
    require(manifest.get("schema") == SCHEMA, "SCHEMA", f"expected {SCHEMA}")
    for key in FORBIDDEN_MANIFEST_KEYS:
        require(key not in manifest, "FORBIDDEN_DEPENDENCY", f"forbidden manifest key {key}")
    require_exact_keys(dict(manifest), MANIFEST_KEYS, "manifest")


def resolve(path_raw: Any, base: Path, code: str, name: str) -> Path:
    require(isinstance(path_raw, str), code, f"{name} is missing")
    path = Path(path_raw)
    if not path.is_absolute():
        path = base / path
    return path


def load_descent_manifest(manifest: Mapping[str, Any], base: Path) -> dict[str, Any]:
    path = resolve(manifest.get("descent_manifest_path"), base, "UPSTREAM_REFERENCE", "descent_manifest_path")
    descent = load_json(path)
    require(
        manifest.get("descent_manifest_sha256") == sha256_json(descent),
        "UPSTREAM_HASH",
        "the #567 descent manifest hash does not match the declared pin",
    )
    require(
        descent.get("schema") == DESCENT_MANIFEST_SCHEMA,
        "UPSTREAM_REFERENCE",
        "the pinned manifest is not a #567 descent manifest",
    )
    return descent


def load_global_form_artifact(descent: Mapping[str, Any], base: Path) -> dict[str, Any]:
    """Load the measured global-form artifact through the descent manifest's pin."""

    path = resolve(
        descent.get("global_form_artifact_path"), base, "GLOBAL_FORM_ARTIFACT", "global_form_artifact_path"
    )
    artifact = load_json(path)
    declared = descent.get("global_form_artifact_sha256")
    require(
        isinstance(declared, str) and declared == artifact.get("artifact_sha256"),
        "UPSTREAM_HASH",
        "the global form artifact hash does not match the descent manifest's pin",
    )
    body = {key: value for key, value in artifact.items() if key != "artifact_sha256"}
    require(
        artifact.get("artifact_sha256") == "sha256:" + sha256_json(body),
        "GLOBAL_FORM_ARTIFACT",
        "the global form artifact self-hash does not recompute",
    )
    require(
        artifact.get("schema") == GLOBAL_FORM_ARTIFACT_SCHEMA and artifact.get("issue") == 567,
        "GLOBAL_FORM_ARTIFACT",
        "the pinned artifact is not a #567 global form artifact",
    )
    axis = artifact.get("six_axis_class_measurement", {})
    require(
        axis.get("axis_count") == 6
        and axis.get("class_group_order") == 6
        and axis.get("smith_invariants") == [1, 1, 1, 1, 1, 6]
        and axis.get("rotation_action_transitive") is True
        and axis.get("antipode_reverses_every_oriented_axis") is True,
        "GLOBAL_FORM_ARTIFACT",
        "the measured six-axis class group block is incomplete",
    )
    deck = artifact.get("federation_deck_action", {})
    require(
        deck.get("deck_group_order") == 120
        and deck.get("charts") == 12
        and deck.get("seams") == 30
        and deck.get("triple_overlaps") == 20,
        "GLOBAL_FORM_ARTIFACT",
        "the measured federation deck action block is incomplete",
    )
    menu = artifact.get("sector_menu", {})
    require(
        menu.get("realized_flux_menu") == [0, 1, 2, 3, 4, 5]
        and menu.get("boundary_smith_invariants_all_unit") is True
        and menu.get("complex") == {"vertices": 12, "seams": 30, "faces": 20}
        and menu.get("puncture_faces", {}).get("antipodal") is True
        and menu.get("single_puncture_impossibility", {}).get(
            "single_puncture_nonzero_flux_impossible"
        )
        is True
        and len(menu.get("flux_tube_witnesses", [])) == 6,
        "GLOBAL_FORM_ARTIFACT",
        "the measured sector menu block is incomplete",
    )
    require(
        artifact.get("federation_sector_class", {}).get("measured_sector_class") == 0
        and artifact.get("refined_sector_menu", {}).get("refinement_natural_sector_menu") is True,
        "GLOBAL_FORM_ARTIFACT",
        "the sector class or refinement naturality block is incomplete",
    )
    return artifact


def load_spin_artifact(manifest: Mapping[str, Any], base: Path) -> dict[str, Any]:
    path = resolve(
        manifest.get("spin_statistics_artifact_path"), base, "SPIN_ARTIFACT", "spin_statistics_artifact_path"
    )
    artifact = load_json(path)
    declared = manifest.get("spin_statistics_artifact_sha256")
    require(
        isinstance(declared, str) and declared == artifact.get("artifact_sha256"),
        "UPSTREAM_HASH",
        "the spin statistics artifact hash does not match the declared pin",
    )
    body = {key: value for key, value in artifact.items() if key != "artifact_sha256"}
    require(
        artifact.get("artifact_sha256") == "sha256:" + sha256_json(body),
        "SPIN_ARTIFACT",
        "the spin statistics artifact self-hash does not recompute",
    )
    require(
        artifact.get("schema") == SPIN_ARTIFACT_SCHEMA and artifact.get("issue") == 314,
        "SPIN_ARTIFACT",
        "the pinned artifact is not a #314 spin statistics artifact",
    )
    return artifact


# ---------------------------------------------------------------------------
# Nerve reconstruction from the pinned carrier manifest
# ---------------------------------------------------------------------------


def edge_list(face_rows: Sequence[Sequence[int]]) -> list[tuple[int, int]]:
    edges = set()
    for a, b, c in face_rows:
        for u, v in ((a, b), (b, c), (c, a)):
            edges.add((min(u, v), max(u, v)))
    return sorted(edges)


def boundary_two(
    face_rows: Sequence[Sequence[int]], edges: Sequence[tuple[int, int]]
) -> list[list[int]]:
    edge_index = {edge: position for position, edge in enumerate(edges)}
    matrix = [[0] * len(face_rows) for _ in range(len(edges))]
    for face_position, (a, b, c) in enumerate(face_rows):
        for u, v in ((a, b), (b, c), (c, a)):
            key = (min(u, v), max(u, v))
            sign = 1 if (u, v) == key else -1
            matrix[edge_index[key]][face_position] += sign
    return matrix


def build_nerve(face_rows: Sequence[Sequence[int]]) -> dict[str, Any]:
    """The seam/face incidence structure of an oriented triangulated nerve.

    Fails closed when any seam does not border exactly two faces with
    opposite orientation signs, which is exactly the missing-triple-overlap
    detection: deleting one face from the twenty leaves its three seams with
    a single bordering face each.
    """

    edges = edge_list(face_rows)
    boundary = boundary_two(face_rows, edges)
    for position, row in enumerate(boundary):
        signs = sorted(value for value in row if value != 0)
        require(
            signs == [-1, 1],
            "NERVE_INCIDENCE",
            f"seam {position} does not border exactly two faces with opposite signs",
        )
    return {"faces": [tuple(face) for face in face_rows], "edges": edges, "boundary": boundary}


def load_carrier_nerve(
    manifest: Mapping[str, Any],
    base: Path,
    global_artifact: Mapping[str, Any],
    spin_artifact: Mapping[str, Any],
) -> dict[str, Any]:
    path = resolve(manifest.get("carrier_manifest_path"), base, "CARRIER_BINDING", "carrier_manifest_path")
    carrier_manifest = load_json(path)
    carrier_sha = sha256_json(carrier_manifest)
    for name, artifact in (("global form", global_artifact), ("spin statistics", spin_artifact)):
        require(
            artifact.get("carrier_binding", {}).get("carrier_manifest_sha256") == carrier_sha,
            "CARRIER_BINDING",
            f"the {name} artifact does not bind the loaded carrier manifest; "
            "the two consumed artifacts must share one A1 federation",
        )
    require(
        carrier_manifest.get("schema") == CARRIER_MANIFEST_SCHEMA,
        "CARRIER_BINDING",
        "the pinned carrier manifest is not a certified selector manifest",
    )
    carrier = carrier_manifest.get("carrier", {})
    ports = list(carrier.get("ports", []))
    require(len(ports) == 12 and len(set(ports)) == 12, "CARRIER_BINDING", "twelve distinct ports required")
    for artifact in (global_artifact, spin_artifact):
        require(
            artifact.get("carrier_binding", {}).get("port_order") == ports,
            "CARRIER_BINDING",
            "an artifact port order differs from the carrier manifest port order",
        )
    index = {port: position for position, port in enumerate(ports)}
    face_rows = [[index[p] for p in face] for face in carrier.get("oriented_faces", [])]
    require(len(face_rows) == 20, "CARRIER_BINDING", "twenty oriented faces required")
    nerve = build_nerve(face_rows)
    declared_edges = sorted(
        (min(index[a], index[b]), max(index[a], index[b])) for a, b in carrier.get("edges", [])
    )
    require(
        declared_edges == nerve["edges"] and len(nerve["edges"]) == 30,
        "CARRIER_BINDING",
        "the declared carrier edges are not the thirty seams of the oriented faces",
    )
    antipode_ports = spin_artifact.get("carrier_binding", {}).get("antipode", {})
    antipode = [index[antipode_ports[port]] for port in ports]
    require(
        all(antipode[antipode[i]] == i and antipode[i] != i for i in range(12)),
        "CARRIER_BINDING",
        "the pinned antipode is not a fixed-point-free involution",
    )
    nerve.update(
        {
            "ports": ports,
            "index": index,
            "antipode": antipode,
            "carrier_manifest_sha256": carrier_sha,
        }
    )
    return nerve


# ---------------------------------------------------------------------------
# Holonomy, routing, and flux-tube checking
# ---------------------------------------------------------------------------


def face_holonomies(
    seam_values: Mapping[int, int],
    boundary: Sequence[Sequence[int]],
    face_count: int,
    modulus: int,
) -> list[int]:
    holonomies = []
    for face in range(face_count):
        total = 0
        for edge, value in seam_values.items():
            total += boundary[edge][face] * value
        holonomies.append(total % modulus)
    return holonomies


def validate_flux_tube(
    nerve: Mapping[str, Any],
    seam_values: Mapping[int, int],
    flux: int,
    start_face: int,
    end_face: int,
    modulus: int,
) -> list[int]:
    """Fail-closed flux-tube check: routing, coherence, and antipodal support.

    A seam assignment passes exactly when its face holonomies are +flux at the
    start puncture, -flux at the end puncture, and zero on every interior
    face, and its support seams form a connected dual path between the two
    punctures.  Tampering any one seam value breaks the holonomy of the
    triple overlaps that seam borders, so a wrongly routed assignment is
    rejected here.
    """

    boundary = nerve["boundary"]
    face_count = len(nerve["faces"])
    for edge, value in seam_values.items():
        require(
            0 <= edge < len(nerve["edges"]) and 0 <= value < modulus,
            "SEAM_ROUTING",
            "a seam value lies outside the seam list or the class group",
        )
    holonomies = face_holonomies(seam_values, boundary, face_count, modulus)
    expected = [0] * face_count
    expected[start_face] = flux % modulus
    expected[end_face] = (-flux) % modulus
    require(
        holonomies == expected,
        "SEAM_ROUTING",
        "the seam assignment violates triple-overlap coherence: its face "
        "holonomies are not +flux at the start puncture, -flux at the end "
        "puncture, and zero on every interior face",
    )
    require(
        sum(holonomies) % modulus == 0,
        "SEAM_ROUTING",
        "face holonomies of seam data must sum to zero on the closed support",
    )
    if flux % modulus != 0:
        degrees: dict[int, int] = {}
        for edge in seam_values:
            for face in range(face_count):
                if boundary[edge][face] != 0:
                    degrees[face] = degrees.get(face, 0) + 1
        endpoints = sorted(face for face, degree in degrees.items() if degree == 1)
        require(
            endpoints == sorted((start_face, end_face))
            and all(degree in (1, 2) for degree in degrees.values()),
            "SEAM_ROUTING",
            "the support seams do not form a dual path between the punctures",
        )
    return holonomies


def dual_path(
    nerve: Mapping[str, Any], start: int, goal: int
) -> list[tuple[int, int]]:
    """A deterministic shortest dual path as (face, shared_edge) steps."""

    boundary = nerve["boundary"]
    neighbors: dict[int, list[tuple[int, int]]] = {}
    for edge in range(len(nerve["edges"])):
        bordering = [face for face in range(len(nerve["faces"])) if boundary[edge][face] != 0]
        left, right = bordering
        neighbors.setdefault(left, []).append((right, edge))
        neighbors.setdefault(right, []).append((left, edge))
    previous: dict[int, tuple[int, int]] = {}
    frontier = [start]
    seen = {start}
    while frontier and goal not in seen:
        nxt: list[int] = []
        for face in frontier:
            for other, edge in sorted(neighbors[face]):
                if other in seen:
                    continue
                seen.add(other)
                previous[other] = (face, edge)
                nxt.append(other)
        frontier = nxt
    require(goal in seen, "SEAM_ROUTING", "no dual path joins the requested faces")
    path: list[tuple[int, int]] = []
    cursor = goal
    while cursor != start:
        parent, edge = previous[cursor]
        path.append((cursor, edge))
        cursor = parent
    path.reverse()
    return path


def flux_tube(
    nerve: Mapping[str, Any], start: int, goal: int, flux: int, modulus: int
) -> dict[int, int]:
    """An exact seam assignment carrying +flux at start and -flux at goal."""

    boundary = nerve["boundary"]
    seam_values: dict[int, int] = {}
    carried = flux % modulus
    current = start
    for step_face, edge in dual_path(nerve, start, goal):
        sign = boundary[edge][current]
        require(sign in (1, -1), "SEAM_ROUTING", "path edge sign must be +-1")
        seam_values[edge] = (sign * carried) % modulus
        current = step_face
    return seam_values


# ---------------------------------------------------------------------------
# Central-column grammar completeness
# ---------------------------------------------------------------------------


def axis_class_lattice(
    nerve: Mapping[str, Any], six_axis: Mapping[str, Any]
) -> dict[str, Any]:
    """Recompute the axis relation lattice the descent chain pins.

    The lattice is free of rank six on the antipodal axes modulo the diagonal
    vector and the zero-sum sublattice; the relation matrix columns are the
    all-ones vector and the difference basis of the zero-sum lattice.  Its
    Smith invariants must recompute to the pinned (1,1,1,1,1,6).
    """

    antipode = nerve["antipode"]
    axes = sorted({tuple(sorted((i, antipode[i]))) for i in range(12)})
    require(len(axes) == 6, "AXIS_CLASS", "six antipodal axes required")
    pinned_axes = [
        sorted(pair) for pair in six_axis.get("axes", [])
    ]
    named_axes = [sorted((nerve["ports"][a], nerve["ports"][b])) for a, b in axes]
    require(
        sorted(named_axes) == sorted(pinned_axes),
        "AXIS_CLASS",
        "the reconstructed antipodal axes differ from the pinned axis list",
    )
    relations = [[1] * 6]
    for i in range(5):
        column = [0] * 6
        column[i] = 1
        column[i + 1] = -1
        relations.append(column)
    relation_matrix = [[relations[j][i] for j in range(6)] for i in range(6)]
    reduced = smith_normal_form(relation_matrix)
    invariants = [reduced[i][i] for i in range(6)]
    require(
        invariants == six_axis.get("smith_invariants"),
        "AXIS_CLASS",
        f"the recomputed axis Smith invariants {invariants} differ from the pinned values",
    )
    require(
        invariants[:-1] == [1] * 5 and invariants[-1] == six_axis.get("class_group_order"),
        "AXIS_CLASS",
        "the Smith normal form does not have unit invariants with a final six",
    )
    return {
        "axes": named_axes,
        "relation_matrix_columns": "all-ones vector and the difference basis of the zero-sum lattice",
        "smith_invariants": invariants,
        "class_group_order": invariants[-1],
        "matches_pinned_measurement": True,
    }


def admit_sector_classes(classes: Sequence[int], class_group_order: int) -> int:
    """Admit a claimed sector-class list; a duplicate class fails closed.

    The Smith computation fixes the class group order, so a list of claimed
    classes is admissible exactly when its residues are pairwise distinct.
    A claimed seventh sector necessarily repeats one of the six residues and
    is rejected here.
    """

    residues = [value % class_group_order for value in classes]
    require(
        len(set(residues)) == len(residues),
        "EXTRA_FLUX",
        "a claimed extra sector class coincides with a listed class modulo "
        f"the order-{class_group_order} class group fixed by the Smith computation",
    )
    return len(residues)


def central_column_completeness(
    nerve: Mapping[str, Any],
    global_artifact: Mapping[str, Any],
) -> dict[str, Any]:
    six_axis = global_artifact["six_axis_class_measurement"]
    sector_menu = global_artifact["sector_menu"]
    class_order = six_axis["class_group_order"]
    lattice = axis_class_lattice(nerve, six_axis)

    boundary = nerve["boundary"]
    reduced = smith_normal_form([row[:] for row in boundary])
    diagonal = [reduced[i][i] for i in range(min(len(reduced), len(reduced[0])))]
    boundary_invariants = [value for value in diagonal if value != 0]
    require(
        boundary_invariants == [1] * 19,
        "SECTOR_SMITH",
        f"the face-boundary Smith invariants are {boundary_invariants}, not nineteen units",
    )
    require(
        sector_menu.get("boundary_smith_invariants_all_unit") is True,
        "SECTOR_SMITH",
        "the pinned menu does not record all-unit boundary Smith invariants",
    )
    column_sums = [
        sum(boundary[edge][face] for face in range(20)) for edge in range(30)
    ]
    require(
        all(value == 0 for value in column_sums),
        "SECTOR_SMITH",
        "each seam must appear with opposite signs in its two faces",
    )

    # Pinned two-puncture witnesses: every class of the order-six group is
    # realized by exactly one witness, and each witness recomputes exactly.
    punctures = sector_menu["puncture_faces"]
    start_face, end_face = punctures["start"], punctures["end"]
    antipode = nerve["antipode"]
    require(
        frozenset(antipode[v] for v in nerve["faces"][start_face])
        == frozenset(nerve["faces"][end_face]),
        "SEAM_ROUTING",
        "the pinned puncture faces are not antipodal on the reconstructed nerve",
    )
    realized_classes = []
    witness_rows = []
    for witness in sector_menu["flux_tube_witnesses"]:
        seam_values = {int(edge): value for edge, value in witness["seam_values"].items()}
        validate_flux_tube(
            nerve, seam_values, witness["flux"], witness["start_face"], witness["end_face"], class_order
        )
        require(
            witness["start_face"] == start_face
            and witness["end_face"] == end_face
            and witness.get("interior_faces_flat") is True
            and witness.get("dual_path_length") == len(seam_values),
            "SEAM_ROUTING",
            "a pinned witness is inconsistent with the pinned puncture data",
        )
        realized_classes.append(witness["flux"] % class_order)
        witness_rows.append({"flux": witness["flux"], "seams_used": len(seam_values)})
    require(
        sorted(realized_classes) == list(range(class_order)),
        "SEAM_ROUTING",
        "the pinned witnesses do not realize every class exactly once",
    )
    admit_sector_classes(realized_classes, class_order)

    # Constructive completeness: a flux tube from the start puncture to every
    # other face realizes the difference generators of the sum-zero holonomy
    # lattice; with the linearity of the holonomy map this realizes every
    # sum-zero central face-holonomy assignment over the class group.
    generator_tubes: dict[int, dict[int, int]] = {}
    for face in range(20):
        if face == start_face:
            continue
        tube = flux_tube(nerve, start_face, face, 1, class_order)
        validate_flux_tube(nerve, tube, 1, start_face, face, class_order)
        generator_tubes[face] = tube
    pair_checks = 0
    for face_a, tube_a in generator_tubes.items():
        for face_b, tube_b in generator_tubes.items():
            if face_a == face_b:
                continue
            for scale in range(1, class_order):
                combined: dict[int, int] = {}
                for edge, value in tube_b.items():
                    combined[edge] = (combined.get(edge, 0) + scale * value) % class_order
                for edge, value in tube_a.items():
                    combined[edge] = (combined.get(edge, 0) - scale * value) % class_order
                holonomies = face_holonomies(combined, boundary, 20, class_order)
                expected = [0] * 20
                expected[face_a] = scale % class_order
                expected[face_b] = (-scale) % class_order
                require(
                    holonomies == expected,
                    "SEAM_ROUTING",
                    "a composed generator tube does not realize the prescribed pair holonomy",
                )
                pair_checks += 1

    # Subgroup obstruction menu: recomputed from the divisors of the class
    # group order and required to equal the pinned menu.
    divisors = [d for d in range(1, class_order + 1) if class_order % d == 0]
    recomputed_menu = {
        f"order_{d}": {
            "liftable_fluxes": [c for c in range(class_order) if (c * d) % class_order == 0],
            "obstructed_fluxes": [c for c in range(class_order) if (c * d) % class_order != 0],
        }
        for d in divisors
    }
    require(
        recomputed_menu == sector_menu["subgroup_obstruction_menu"],
        "SECTOR_SMITH",
        "the recomputed subgroup obstruction menu differs from the pinned menu",
    )

    return {
        "axis_class_lattice": lattice,
        "boundary_smith_invariants": boundary_invariants,
        "face_holonomy_sum_identically_zero": True,
        "puncture_faces_antipodal": True,
        "pinned_witnesses_recomputed": witness_rows,
        "realized_class_menu": sorted(realized_classes),
        "generator_tubes_constructed": len(generator_tubes),
        "pair_realization_checks": pair_checks,
        "subgroup_obstruction_menu_recomputed": True,
        "completeness_statement": (
            "the Smith normal form of the axis relation lattice has unit "
            "invariants except the final six, so the central class group is "
            "cyclic of order six; the face-boundary matrix has nineteen unit "
            "Smith invariants and zero column sums, so the routed-seam "
            "holonomy map surjects onto the sum-zero central 2-cochain "
            "lattice; every central 2-cochain class on the nerve is realized "
            "by exactly one of the six measured sectors, and every "
            "seam-routable central flux mechanism lies in the measured "
            "order-six menu"
        ),
        "status": "exact",
    }


# ---------------------------------------------------------------------------
# Composition laws on the finite nerve
# ---------------------------------------------------------------------------


def composition_laws(
    nerve: Mapping[str, Any],
    witnesses: Sequence[Mapping[str, Any]],
    class_order: int,
) -> dict[str, Any]:
    boundary = nerve["boundary"]
    edges = nerve["edges"]
    faces = nerve["faces"]
    cochains = [
        {int(edge): value for edge, value in witness["seam_values"].items()}
        for witness in witnesses
    ]

    def holonomy(cochain: Mapping[int, int]) -> list[int]:
        return face_holonomies(cochain, boundary, len(faces), class_order)

    # Pairwise seam composition: the holonomy map is additive on seam data.
    pairwise_checks = 0
    for left in cochains:
        for right in cochains:
            summed: dict[int, int] = dict(left)
            for edge, value in right.items():
                summed[edge] = (summed.get(edge, 0) + value) % class_order
            left_h, right_h, sum_h = holonomy(left), holonomy(right), holonomy(summed)
            require(
                sum_h == [(left_h[f] + right_h[f]) % class_order for f in range(len(faces))],
                "COMPOSITION",
                "pairwise seam composition does not commute with face holonomy",
            )
            pairwise_checks += 1
    scalar_checks = 0
    for edge in range(len(edges)):
        base = holonomy({edge: 1})
        for scale in range(class_order):
            scaled = holonomy({edge: scale % class_order})
            require(
                scaled == [(scale * base[f]) % class_order for f in range(len(faces))],
                "COMPOSITION",
                "seam scaling does not commute with face holonomy",
            )
            scalar_checks += 1

    # Triple-overlap coherence: the composed boundary vanishes identically,
    # so the coboundary of every seam assignment is a closed face-holonomy
    # assignment, and every face-holonomy assignment sums to zero.
    vertex_boundary = [[0] * len(edges) for _ in range(12)]
    for position, (u, v) in enumerate(edges):
        vertex_boundary[u][position] -= 1
        vertex_boundary[v][position] += 1
    composed_entries = 0
    for vertex in range(12):
        for face in range(len(faces)):
            total = sum(
                vertex_boundary[vertex][edge] * boundary[edge][face]
                for edge in range(len(edges))
            )
            require(
                total == 0,
                "COMPOSITION",
                "the composed boundary does not vanish: a face boundary is not closed",
            )
            composed_entries += 1

    # Routed loops: the boundary loop of each face has holonomy equal to that
    # face's holonomy, loop concatenation is associative, and two routed
    # loops differing by one face have holonomies differing by that face's
    # boundary class, exhaustively over the twenty faces.
    face_chains = []
    for face in range(len(faces)):
        chain = tuple(boundary[edge][face] for edge in range(len(edges)))
        face_chains.append(chain)

    def loop_holonomy(chain: Sequence[int], cochain: Mapping[int, int]) -> int:
        return sum(chain[edge] * value for edge, value in cochain.items()) % class_order

    loop_checks = 0
    for face, chain in enumerate(face_chains):
        for cochain in cochains:
            require(
                loop_holonomy(chain, cochain) == holonomy(cochain)[face],
                "ROUTED_LOOP",
                "a face boundary loop holonomy differs from the face holonomy",
            )
            loop_checks += 1
    base_chain = face_chains[0]
    homotopy_checks = 0
    for face, chain in enumerate(face_chains):
        moved = tuple(base_chain[e] + chain[e] for e in range(len(edges)))
        for cochain in cochains:
            require(
                loop_holonomy(moved, cochain)
                == (loop_holonomy(base_chain, cochain) + holonomy(cochain)[face]) % class_order,
                "ROUTED_LOOP",
                "a homotopy move across a face does not shift the loop holonomy "
                "by that face's boundary class",
            )
            homotopy_checks += 1
    associativity_checks = 0
    for chain_a in face_chains:
        for chain_b in face_chains:
            left_pair = tuple(chain_a[e] + chain_b[e] for e in range(len(edges)))
            for chain_c in face_chains:
                left = tuple(left_pair[e] + chain_c[e] for e in range(len(edges)))
                right_pair = tuple(chain_b[e] + chain_c[e] for e in range(len(edges)))
                right = tuple(chain_a[e] + right_pair[e] for e in range(len(edges)))
                require(
                    left == right,
                    "ROUTED_LOOP",
                    "loop concatenation is not associative on the finite nerve",
                )
                associativity_checks += 1

    return {
        "pairwise_composition_checks": pairwise_checks,
        "seam_scaling_checks": scalar_checks,
        "triple_overlap_coherence_entries": composed_entries,
        "coboundary_lands_in_class_group": True,
        "face_loop_holonomy_checks": loop_checks,
        "homotopy_move_checks": homotopy_checks,
        "loop_associativity_checks": associativity_checks,
        "status": "closed",
    }


# ---------------------------------------------------------------------------
# Frame reconstruction and Galois-conjugate control
# ---------------------------------------------------------------------------


def frame_and_conjugate(
    nerve: Mapping[str, Any], spin_artifact: Mapping[str, Any]
) -> dict[str, Any]:
    """Exact frame checks and the conjugate-frame invariance of the menu.

    The pinned port-vertex frame is parsed over Q(sqrt(5)).  Adjacency at
    squared distance four reproduces the thirty seams, the antipode is exact
    negation, and all twenty oriented face determinants are positive.  The
    Galois conjugation sqrt(5) -> -sqrt(5) leaves the seam set, the face set,
    and the orientation signs invariant, so the conjugate frame carries the
    identical boundary matrix and the identical class menu.
    """

    frame_raw = spin_artifact["port_vertex_frame"]
    coordinates = [
        tuple(parse_q5(component) for component in frame_raw[port]) for port in nerve["ports"]
    ]
    antipode = nerve["antipode"]
    for i in range(12):
        require(
            all(coordinates[antipode[i]][k] == q5_neg(coordinates[i][k]) for k in range(3)),
            "FRAME",
            "the pinned frame antipode is not exact negation",
        )

    def squared_distance(u: Sequence[Q5], v: Sequence[Q5]) -> Q5:
        total: Q5 = (Fraction(0), Fraction(0))
        for k in range(3):
            difference = q5_sub(u[k], v[k])
            total = q5_add(total, q5_mul(difference, difference))
        return total

    four: Q5 = (Fraction(4), Fraction(0))

    def frame_edges(points: Sequence[Sequence[Q5]]) -> list[tuple[int, int]]:
        return sorted(
            (i, j)
            for i in range(12)
            for j in range(i + 1, 12)
            if squared_distance(points[i], points[j]) == four
        )

    require(
        frame_edges(coordinates) == nerve["edges"],
        "FRAME",
        "adjacency at squared distance four does not reproduce the thirty seams",
    )
    determinants = [q5_det3([coordinates[v] for v in face]) for face in nerve["faces"]]
    require(
        all(q5_is_positive(det) for det in determinants),
        "FRAME",
        "an oriented face determinant is not positive in the pinned frame",
    )

    conjugate = [tuple(q5_conjugate(c) for c in point) for point in coordinates]
    require(
        frame_edges(conjugate) == nerve["edges"],
        "CONJUGATE_FRAME",
        "the Galois-conjugate frame does not reproduce the same thirty seams",
    )
    conjugate_determinants = [q5_det3([conjugate[v] for v in face]) for face in nerve["faces"]]
    require(
        all(q5_is_positive(det) for det in conjugate_determinants),
        "CONJUGATE_FRAME",
        "the Galois-conjugate frame does not preserve the face orientations",
    )
    return {
        "frame_arithmetic": "exact_q_sqrt5",
        "antipode_is_exact_negation": True,
        "seams_from_squared_distance_four": 30,
        "oriented_face_determinants_positive": 20,
        "conjugate_frame_same_seams": True,
        "conjugate_frame_same_oriented_faces": True,
        "conjugate_frame_same_boundary_matrix": True,
        "class_menu_invariant_under_conjugation": True,
    }


def reject_frame_dependent_sector(
    claimed_class: int, base_menu: Sequence[int], class_group_order: int
) -> None:
    """A claimed conjugate-frame-only sector fails closed.

    The conjugate frame carries the identical boundary matrix, so its class
    menu equals the base menu; any claimed frame-dependent extra sector has a
    class residue that the base menu carries and the claim is rejected.
    """

    residue = claimed_class % class_group_order
    require(
        residue not in list(base_menu),
        "CONJUGATE_FRAME",
        f"the claimed frame-dependent sector reduces to class {residue}, which "
        "the frame-independent menu carries; no extra sector exists",
    )


# ---------------------------------------------------------------------------
# Nonabelian scope boundary
# ---------------------------------------------------------------------------


def nonabelian_boundary(spin_artifact: Mapping[str, Any]) -> dict[str, Any]:
    """The measured double-cover transport as the mandatory scope countermodel.

    The pinned lift measurement carries a 120-element unit-quaternion group
    with centre of order two.  The pinned Klein-four lift table supplies three
    unit quaternions; exact quaternion arithmetic shows each squares to minus
    one and no two commute, so the transport values do not lie in any abelian
    group and in particular not in the central order-six column.
    """

    lift = spin_artifact.get("lift_measurement", {})
    require(
        lift.get("lift_group_order") == 120
        and lift.get("centre_order") == 2
        and lift.get("centre_elements") == ["+1", "-1"]
        and lift.get("unique_nontrivial_involution") == "-1",
        "NONABELIAN_WITNESS",
        "the pinned lift measurement does not carry the binary-icosahedral profile",
    )
    order_profile = lift.get("order_profile", {})
    require(
        sum(order_profile.values()) == 120 and order_profile.get("2") == 1,
        "NONABELIAN_WITNESS",
        "the pinned order profile does not sum to 120 with a unique involution",
    )
    table = spin_artifact.get("canonical_klein_four_lift_table", [])
    require(len(table) == 3, "NONABELIAN_WITNESS", "three Klein-four lifts required")
    lifts = []
    for row in table:
        quaternion = tuple(Fraction(component) for component in row["quaternion_lift"])
        require(
            row.get("lift_square") == "-1"
            and quaternion_multiply(quaternion, quaternion) == QUATERNION_MINUS_ONE,
            "NONABELIAN_WITNESS",
            "a pinned Klein-four lift does not square to minus one",
        )
        lifts.append(quaternion)
    noncommuting_pairs = 0
    for i in range(3):
        for j in range(i + 1, 3):
            require(
                quaternion_multiply(lifts[i], lifts[j])
                != quaternion_multiply(lifts[j], lifts[i]),
                "NONABELIAN_WITNESS",
                "two pinned Klein-four lifts commute; the countermodel fails",
            )
            noncommuting_pairs += 1
    require(
        spin_artifact.get("section_obstruction", {}).get(
            "no_section_over_any_klein_four_subgroup"
        )
        is True,
        "NONABELIAN_WITNESS",
        "the pinned section obstruction is missing",
    )
    return {
        "artifact_sha256": spin_artifact["artifact_sha256"],
        "transport_group_order": 120,
        "centre_order": 2,
        "unique_nontrivial_involution": "-1",
        "klein_four_lift_squares_minus_one": 3,
        "noncommuting_lift_pairs": noncommuting_pairs,
        "double_cover_non_split": True,
        "scope_statement": (
            "the measured double-cover transport is a routed-seam mechanism "
            "on the same federation with values in a nonabelian 120-element "
            "group whose centre has order two; those values do not lie in "
            "the central order-six column, so a quotient-visible routed-seam "
            "structure exists outside the central-column grammar and grammar "
            "exhaustiveness for all routed mechanisms is a named open "
            "interface"
        ),
    }


# ---------------------------------------------------------------------------
# Structural fail-closed controls
# ---------------------------------------------------------------------------


def structural_controls(
    nerve: Mapping[str, Any], global_artifact: Mapping[str, Any]
) -> list[dict[str, Any]]:
    sector_menu = global_artifact["sector_menu"]
    class_order = global_artifact["six_axis_class_measurement"]["class_group_order"]
    menu = sector_menu["realized_flux_menu"]
    results = []

    def run(name: str, expected_code: str, action) -> None:
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

    witness = sector_menu["flux_tube_witnesses"][1]
    tampered = {int(edge): value for edge, value in witness["seam_values"].items()}
    first_edge = sorted(tampered)[0]
    tampered[first_edge] = (tampered[first_edge] + 1) % class_order

    def wrong_routing() -> None:
        validate_flux_tube(
            nerve, tampered, witness["flux"], witness["start_face"], witness["end_face"], class_order
        )

    run("wrong_routing_one_seam_tampered", "SEAM_ROUTING", wrong_routing)

    def missing_triple_overlap() -> None:
        build_nerve([list(face) for face in nerve["faces"][:-1]])

    run("missing_triple_overlap_deleted_face", "NERVE_INCIDENCE", missing_triple_overlap)

    def frame_dependent_sector() -> None:
        reject_frame_dependent_sector(3, menu, class_order)

    run("conjugate_frame_extra_sector_claim", "CONJUGATE_FRAME", frame_dependent_sector)

    def seventh_sector() -> None:
        admit_sector_classes(list(menu) + [6], class_order)

    run("extra_flux_seventh_sector_claim", "EXTRA_FLUX", seventh_sector)

    return results


# ---------------------------------------------------------------------------
# Certificate payload
# ---------------------------------------------------------------------------


def certificate_payload(
    manifest: Mapping[str, Any], base_dir: Path | None = None
) -> dict[str, Any]:
    base = base_dir or MODULE_DIR
    validate_manifest(manifest)
    descent = load_descent_manifest(manifest, base)
    global_artifact = load_global_form_artifact(descent, base)
    spin_artifact = load_spin_artifact(manifest, base)
    nerve = load_carrier_nerve(manifest, base, global_artifact, spin_artifact)

    frame = frame_and_conjugate(nerve, spin_artifact)
    class_order = global_artifact["six_axis_class_measurement"]["class_group_order"]
    central = central_column_completeness(nerve, global_artifact)
    laws = composition_laws(
        nerve, global_artifact["sector_menu"]["flux_tube_witnesses"], class_order
    )
    countermodel = nonabelian_boundary(spin_artifact)
    controls = structural_controls(nerve, global_artifact)

    return {
        "schema": RECEIPT_SCHEMA,
        "issue": 613,
        "manifest_sha256": sha256_json(manifest),
        "consumed_upstream": {
            "descent_manifest_sha256": sha256_json(descent),
            "global_form_artifact_sha256": global_artifact["artifact_sha256"],
            "spin_statistics_artifact_sha256": spin_artifact["artifact_sha256"],
            "carrier_manifest_sha256": nerve["carrier_manifest_sha256"],
            "shared_federation": (
                "both consumed artifacts bind the same certified carrier, so "
                "the central column and the nonabelian countermodel live on "
                "one A1 federation of twelve charts, thirty seams, and "
                "twenty triple overlaps"
            ),
            "global_form_fields": [
                "six_axis_class_measurement",
                "federation_deck_action",
                "federation_sector_class",
                "sector_menu",
                "refined_sector_menu.refinement_natural_sector_menu",
                "carrier_binding",
            ],
            "spin_statistics_fields": [
                "carrier_binding",
                "port_vertex_frame",
                "lift_measurement",
                "canonical_klein_four_lift_table",
                "section_obstruction",
            ],
        },
        "federation_nerve": {
            "charts": 12,
            "seams": 30,
            "triple_overlaps": 20,
            "deck_group_order": global_artifact["federation_deck_action"]["deck_group_order"],
            "reference_sector_class": global_artifact["federation_sector_class"][
                "measured_sector_class"
            ],
            "refinement_natural_sector_menu": True,
        },
        "frame_reconstruction": frame,
        "central_column_grammar": central,
        "composition_laws": laws,
        "nonabelian_scope_boundary": countermodel,
        "structural_controls": controls,
        "verdict": {
            "central_column_completeness": "exact",
            "central_class_group_smith_invariants": central["axis_class_lattice"][
                "smith_invariants"
            ],
            "composition_laws": "closed",
            "general_grammar": "conditional_open_interface",
            "general_grammar_witness": {
                "artifact_sha256": countermodel["artifact_sha256"],
                "transport_group_order": countermodel["transport_group_order"],
                "centre_order": countermodel["centre_order"],
                "noncommuting_lift_pairs": countermodel["noncommuting_lift_pairs"],
            },
            "four_dimensional_instanton_theta": "separate_gates",
            "controls": {
                row["name"]: row["passed"] for row in controls
            },
        },
        "claim_boundary": {
            "proves": (
                "exact completeness of the routed-seam grammar with seam data "
                "valued in the central order-six class group, closure of the "
                "pairwise, triple-overlap, and routed-loop composition laws "
                "on the finite nerve, and the existence of a measured "
                "routed-seam mechanism with noncentral nonabelian values on "
                "the same federation"
            ),
            "does_not_close": [
                "grammar exhaustiveness for routed mechanisms with values outside the central column",
                "four-dimensional instanton normalization and theta periodicity (separate gates)",
                "laboratory current lines (issue 569)",
                "continuum quantum field theory",
            ],
            "bounded_exit": "conditional_open_interface",
            "bounded_exit_scope": (
                "the central column exits exact_named_realization inside the "
                "measured order-six menu; the general grammar exits "
                "conditional_open_interface on the nonabelian witness"
            ),
        },
        "verifier_command": (
            "python3 code/a5_closure/routed_seam_grammar_certificate.py verify "
            "--manifest code/a5_closure/manifests/routed_seam_grammar_reference.json "
            "--receipt code/a5_closure/receipts/routed_seam_grammar_reference.receipt.json"
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
    wrong_schema["schema"] = "oph.routed_seam_grammar_certificate.v0"
    cases.append(("wrong_manifest_schema", wrong_schema, "SCHEMA"))

    descent_drift = copy.deepcopy(dict(manifest))
    descent_drift["descent_manifest_sha256"] = "0" * 64
    cases.append(("descent_manifest_pin_drift", descent_drift, "UPSTREAM_HASH"))

    spin_drift = copy.deepcopy(dict(manifest))
    spin_drift["spin_statistics_artifact_sha256"] = "sha256:" + "0" * 64
    cases.append(("spin_artifact_pin_drift", spin_drift, "UPSTREAM_HASH"))

    swapped = copy.deepcopy(dict(manifest))
    swapped["spin_statistics_artifact_path"] = str(manifest.get("descent_manifest_path"))
    cases.append(("spin_artifact_path_swapped", swapped, "UPSTREAM_HASH"))

    wrong_carrier = copy.deepcopy(dict(manifest))
    wrong_carrier["carrier_manifest_path"] = str(manifest.get("descent_manifest_path"))
    cases.append(("carrier_manifest_swapped", wrong_carrier, "CARRIER_BINDING"))

    for name, key in (
        ("general_grammar_promotion", "general_grammar_closed"),
        ("nonabelian_completion_injection", "nonabelian_grammar_completion"),
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
        "issue": 613,
        "manifest_sha256": sha256_json(manifest),
        "finite_controls": results,
        "countermodel_witnesses": {
            "wrong_routing": (
                "one tampered seam value breaks the holonomy of the triple "
                "overlaps that seam borders and the composition checker rejects it"
            ),
            "missing_triple_overlap": (
                "deleting one face leaves three seams with a single bordering "
                "face each and the incidence check rejects the nerve"
            ),
            "conjugate_frame": (
                "the Galois-conjugate frame carries the identical seam set, "
                "oriented faces, and boundary matrix, so a claimed "
                "frame-dependent extra sector reduces to a class the "
                "frame-independent menu carries and is rejected"
            ),
            "extra_flux": (
                "a claimed seventh sector class repeats one of the six "
                "residues fixed by the Smith computation and is rejected"
            ),
            "nonabelian_witness": (
                "the pinned binary-icosahedral transport keeps general "
                "grammar exhaustiveness a named open interface; a manifest "
                "cannot promote the central-column result to a general "
                "completeness claim"
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
        MODULE_DIR / "manifests" / "routed_seam_grammar_reference.json",
        MODULE_DIR / "receipts" / "routed_seam_grammar_reference.receipt.json",
        MODULE_DIR / "negative_controls" / "issue_613_negative_controls.json",
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
