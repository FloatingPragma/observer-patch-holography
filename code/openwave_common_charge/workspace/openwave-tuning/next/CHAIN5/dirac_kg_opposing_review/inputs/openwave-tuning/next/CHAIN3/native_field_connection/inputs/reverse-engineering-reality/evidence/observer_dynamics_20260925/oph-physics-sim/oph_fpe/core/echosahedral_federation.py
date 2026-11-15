"""Two-level echosahedral carrier/federation conformance contracts.

The microscopic carrier and the inter-carrier federation are deliberately
different objects.  Each carrier has a hidden regular-icosahedron
presentation with twelve ports; seams connect typed, connected subsets of
those local ports.  This module validates that finite implementation surface.

It does *not* identify one carrier with an observer, screen cap, S2 point, H3
point, event, or BW/KMS source.  Those remain downstream, independently gated
outputs.  Existing screen-port reports can be checked through a fail-closed
reference bridge, and a compact JSON instrument-bundle API is provided, but
neither path promotes the current engine to a physical federation source.
"""

from __future__ import annotations

import argparse
import copy
from dataclasses import dataclass, replace
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Iterable, Literal, Mapping, Sequence

import networkx as nx
import numpy as np

from oph_fpe.core.icosahedral import (
    build_geodesic_icosahedral_tower,
    icosahedral_a5_port_permutations,
)
from oph_fpe.core.screen_ports import (
    echosahedral_patch_architecture_report,
    echosahedral_port_names,
)


CollarKind = Literal[
    "single_port",
    "antipodal_pair",
    "edge_bundle",
    "face_collar",
    "connected_bundle",
]
BoundaryCondition = Literal[
    "open_external",
    "fixed_external",
    "measured_external",
]

_SHA256_RE = re.compile(r"^sha256:[0-9a-f]{64}$")
_PORT_COUNT = 12
_REFERENCE_CARRIER_TEMPLATE_ID = "__shared_regular_icosahedron_template__"
_REFERENCE_STRUCTURE_TOKEN = "shared_regular_icosahedron_12_30_20_a5_v1"
_REPORT_DETAIL_LIMIT = 64
_CARDINALITY_SEMANTICS = (
    "exact_declared_finite_source_carrier_cardinality_separate_from_support_regulator"
)
_INTERFACE_HASH_CHECK_SCOPE = (
    "content_addressed_canonical_checked_matrix_schema_no_higher_overlap_proof"
)
_EMBEDDED_LOCAL_TEMPLATE_FIELDS = frozenset(
    {
        "port_coordinates",
        "hidden_port_coordinates",
        "port_names",
        "local_port_names",
        "edges",
        "faces",
        "antipode",
        "a5_actions",
        "local_a5_frame",
    }
)

HIDDEN_PRESENTATION_FIELDS = frozenset(
    {
        "port_coordinates",
        "hidden_port_coordinates",
        "port_names",
        "local_port_names",
        "local_a5_frame",
        "a5_actions",
        "worker_id",
        "worker_assignment",
        "shard_id",
        "shard_assignment",
        "queue_order",
        "queue_position",
        "memory_layout",
        "repair_iteration",
        "repair_depth",
        "retry_count",
        "candidate_kappa",
        "candidate_clock_scale",
        "target_h3_label",
        "local_outward_normal",
    }
)

CARRIER_FORBIDDEN_PROMOTION_FIELDS = HIDDEN_PRESENTATION_FIELDS | frozenset(
    {
        "s2_point",
        "s2_cap",
        "h3_point",
        "event_position",
        "event_id_from_carrier_id",
        "bw_source_from_carrier_coordinates",
        "cap_normal_from_local_carrier",
    }
)

QUOTIENT_VISIBLE_CARRIER_FIELDS = frozenset(
    {
        "carrier_id",
        "port_response",
        "seam_interface_packet",
        "semantic_record",
        "checkpoint_continuation",
        "repair_normal_form",
        "interface_algebra_sha256",
    }
)


@dataclass(frozen=True)
class EchosahedralCarrier:
    """One finite twelve-port hidden carrier presentation."""

    carrier_id: str
    port_names: tuple[str, ...]
    port_coordinates: tuple[tuple[float, float, float], ...]
    edges: tuple[tuple[int, int], ...]
    faces: tuple[tuple[int, int, int], ...]
    antipode: tuple[int, ...]
    a5_actions: tuple[tuple[int, ...], ...]


@dataclass(frozen=True)
class InterfaceAlgebraBinding:
    """Hashes binding the same overlap-visible algebra at both endpoints."""

    interface_algebra_id: str
    interface_algebra_sha256: str
    left_interface_algebra_sha256: str
    right_interface_algebra_sha256: str


@dataclass(frozen=True)
class SeamBundle:
    """A typed connected collar with an explicit orientation-reversing bijection."""

    seam_id: str
    left_carrier_id: str
    right_carrier_id: str
    left_ports: tuple[int, ...]
    right_ports: tuple[int, ...]
    left_to_right_ports: tuple[int, ...]
    right_to_left_ports: tuple[int, ...]
    left_to_right_orientation: tuple[int, ...]
    right_to_left_orientation: tuple[int, ...]
    collar_kind: CollarKind
    interface_algebra: InterfaceAlgebraBinding


@dataclass(frozen=True)
class ExternalBoundaryBundle:
    """An explicit declaration for ports not sewn to another carrier."""

    boundary_id: str
    carrier_id: str
    ports: tuple[int, ...]
    boundary_condition: BoundaryCondition
    boundary_algebra_sha256: str


@dataclass(frozen=True)
class ObserverSupport:
    """One observer token supported on a connected carrier subfederation."""

    observer_token: str
    carrier_ids: frozenset[str]
    visible_seam_ids: frozenset[str]
    record_algebra_sha256: str
    checkpoint_cut_sha256: str


@dataclass(frozen=True)
class TripleOverlapBundle:
    """One explicitly sourced nonempty three-chart overlap.

    ``oriented_carrier_ids`` lists the charts in positive cyclic order and
    ``oriented_seam_ids`` lists the corresponding pairwise overlaps
    ``(01, 12, 20)``.  Each pairwise collar restricts to the same canonical
    matrix algebra on the triple overlap.  In the incidence-nerve producer
    below that algebra is ``M_1(C)``: the three restriction maps are therefore
    the unique unital star homomorphisms and their composite is exactly the
    identity.  The explicit bundle matters because a triangle in the carrier
    adjacency graph alone does not say that the three pairwise overlaps have
    a nonempty common restriction.
    """

    overlap_id: str
    oriented_carrier_ids: tuple[str, str, str]
    oriented_seam_ids: tuple[str, str, str]
    restriction_algebra_id: str
    restriction_algebra_sha256: str


@dataclass(frozen=True)
class SupportTowerBinding:
    """Source-visible simplicial map from federation nerve to support tower."""

    geometry_family: str
    carrier_to_base_vertex: tuple[tuple[str, int], ...]
    orientation: str
    refinement_rule: str
    source_incidence_sha256: str


@dataclass(frozen=True)
class EchosahedralFederation:
    """Carrier set, typed seams, external boundary, and observer supports."""

    federation_id: str
    carriers: tuple[EchosahedralCarrier, ...]
    seams: tuple[SeamBundle, ...]
    external_boundaries: tuple[ExternalBoundaryBundle, ...]
    observer_supports: tuple[ObserverSupport, ...] = ()
    triple_overlaps: tuple[TripleOverlapBundle, ...] = ()
    support_tower_binding: SupportTowerBinding | None = None


def interface_algebra_sha256(schema: Mapping[str, Any] | Sequence[Any] | str) -> str:
    """Hash a JSON-serializable interface-algebra schema."""

    encoded = json.dumps(
        schema, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return "sha256:" + hashlib.sha256(encoded).hexdigest()


def finite_matrix_interface_algebra_schema(
    interface_algebra_id: str, matrix_size: int
) -> dict[str, Any]:
    """Canonical schema for the matrix algebra actually checked on a seam.

    The seam homomorphism verifier below proves conjugation on ``M_n(C)``,
    where ``n`` is the collar-port count.  Binding only three equal opaque
    hashes did not establish that the hashed object was this algebra (or even
    that the hash had a preimage available to the verifier).  The canonical
    schema makes the checked algebra and its dimension part of the hash
    preimage.
    """

    size = int(matrix_size)
    if size <= 0:
        raise ValueError("matrix_size must be positive")
    return {
        "schema": "oph.finite-interface-matrix-star-algebra.v1",
        "interface_algebra_id": str(interface_algebra_id),
        "scalar_field": "complex",
        "matrix_size": size,
    }


@lru_cache(maxsize=1)
def _reference_echosahedral_carrier_template() -> EchosahedralCarrier:
    """Build the immutable 12-port template once for all reference carriers."""

    base = build_geodesic_icosahedral_tower(0).levels[0]
    coordinates = np.asarray(base.vertices, dtype=float)
    antipode = tuple(
        int(np.argmin(np.linalg.norm(coordinates + coordinate, axis=1)))
        for coordinate in coordinates
    )
    return EchosahedralCarrier(
        carrier_id=_REFERENCE_CARRIER_TEMPLATE_ID,
        port_names=tuple(echosahedral_port_names(12)),
        port_coordinates=tuple(
            tuple(float(value) for value in row) for row in coordinates
        ),
        edges=tuple(tuple(int(value) for value in row) for row in base.edges),
        faces=tuple(tuple(int(value) for value in row) for row in base.faces),
        antipode=antipode,
        a5_actions=icosahedral_a5_port_permutations(),
    )


def reference_echosahedral_carrier(carrier_id: str) -> EchosahedralCarrier:
    """Instantiate one ID over the shared immutable 12-port template.

    Every tuple-valued presentation field is shared with the single cached
    reference template.  Only ``carrier_id`` is replaced, so a federation of
    canonical carriers does not rebuild or duplicate the 12/30/20/A5 data.
    """

    if not isinstance(carrier_id, str) or not carrier_id:
        raise ValueError("carrier_id must be a nonempty string")
    return replace(_reference_echosahedral_carrier_template(), carrier_id=carrier_id)


def _reference_incidence_sha256() -> str:
    """Bind the oriented incidence data used by the federation-nerve producer."""

    template = _reference_echosahedral_carrier_template()
    return interface_algebra_sha256(
        {
            "schema": "oph.echosahedral_oriented_incidence.v1",
            "port_count": len(template.port_names),
            "edges": [list(edge) for edge in sorted(template.edges)],
            "oriented_faces": [list(face) for face in template.faces],
            "antipode": list(template.antipode),
        }
    )


def reference_incidence_nerve_federation() -> EchosahedralFederation:
    """Derive a nonvacuous federation cover from the carrier incidence itself.

    There is one chart carrier for each of the twelve source ports, one
    pairwise seam for each of the thirty source edges, and one declared triple
    restriction for each of the twenty oriented source faces.  For an edge
    ``{u,v}``, chart ``u`` spends local port ``v`` and chart ``v`` spends local
    port ``u``.  Every seam endpoint is consequently owned exactly once.

    This is the nerve of the certified oriented icosahedral boundary, not an
    SM-shaped or target-labelled cover.  The construction is deterministic
    from the already-audited ``(12,30,20)`` incidence packet.
    """

    template = _reference_echosahedral_carrier_template()
    carrier_ids = tuple(f"incidence-chart-{index:02d}" for index in range(_PORT_COUNT))
    carriers = tuple(reference_echosahedral_carrier(item) for item in carrier_ids)
    interface_id = "incidence-nerve-scalar-restriction-v1"
    interface_hash = interface_algebra_sha256(
        finite_matrix_interface_algebra_schema(interface_id, 1)
    )
    binding = InterfaceAlgebraBinding(
        interface_algebra_id=interface_id,
        interface_algebra_sha256=interface_hash,
        left_interface_algebra_sha256=interface_hash,
        right_interface_algebra_sha256=interface_hash,
    )

    def seam_id(left: int, right: int) -> str:
        low, high = sorted((int(left), int(right)))
        return f"incidence-seam-{low:02d}-{high:02d}"

    seams = tuple(
        SeamBundle(
            seam_id=seam_id(left, right),
            left_carrier_id=carrier_ids[left],
            right_carrier_id=carrier_ids[right],
            left_ports=(right,),
            right_ports=(left,),
            left_to_right_ports=(left,),
            right_to_left_ports=(right,),
            left_to_right_orientation=(-1,),
            right_to_left_orientation=(-1,),
            collar_kind="single_port",
            interface_algebra=binding,
        )
        for left, right in sorted(template.edges)
    )
    triple_overlaps = tuple(
        TripleOverlapBundle(
            overlap_id=f"incidence-face-{face_index:02d}",
            oriented_carrier_ids=tuple(carrier_ids[index] for index in face),
            oriented_seam_ids=(
                seam_id(face[0], face[1]),
                seam_id(face[1], face[2]),
                seam_id(face[2], face[0]),
            ),
            restriction_algebra_id=interface_id,
            restriction_algebra_sha256=interface_hash,
        )
        for face_index, face in enumerate(template.faces)
    )

    neighbors: dict[int, set[int]] = {index: set() for index in range(_PORT_COUNT)}
    for left, right in template.edges:
        neighbors[left].add(right)
        neighbors[right].add(left)
    boundary_hash = interface_algebra_sha256(
        {
            "schema": "oph.incidence_nerve.external_boundary.v1",
            "source_incidence_sha256": _reference_incidence_sha256(),
        }
    )
    external_boundaries: list[ExternalBoundaryBundle] = []
    for carrier_index, carrier in enumerate(carriers):
        remaining = set(range(_PORT_COUNT)) - neighbors[carrier_index]
        graph = nx.Graph()
        graph.add_nodes_from(remaining)
        graph.add_edges_from(
            (left, right)
            for left, right in carrier.edges
            if left in remaining and right in remaining
        )
        for component_index, component in enumerate(
            sorted(
                (tuple(sorted(item)) for item in nx.connected_components(graph)),
                key=lambda item: (item[0], len(item), item),
            )
        ):
            external_boundaries.append(
                ExternalBoundaryBundle(
                    boundary_id=(
                        f"incidence-boundary-{carrier_index:02d}-"
                        f"{component_index:02d}"
                    ),
                    carrier_id=carrier.carrier_id,
                    ports=component,
                    boundary_condition="open_external",
                    boundary_algebra_sha256=boundary_hash,
                )
            )

    record_hash = interface_algebra_sha256(
        {
            "schema": "oph.incidence_nerve.observer_record.v1",
            "source_incidence_sha256": _reference_incidence_sha256(),
        }
    )
    checkpoint_hash = interface_algebra_sha256(
        {
            "schema": "oph.incidence_nerve.observer_checkpoint.v1",
            "source_incidence_sha256": _reference_incidence_sha256(),
        }
    )
    observer = ObserverSupport(
        observer_token="incidence-nerve-observer",
        carrier_ids=frozenset(carrier_ids),
        visible_seam_ids=frozenset(seam.seam_id for seam in seams),
        record_algebra_sha256=record_hash,
        checkpoint_cut_sha256=checkpoint_hash,
    )
    support_binding = SupportTowerBinding(
        geometry_family="nested_geodesic_icosahedral",
        carrier_to_base_vertex=tuple(
            (carrier_id, index) for index, carrier_id in enumerate(carrier_ids)
        ),
        orientation="source_oriented_faces_outward",
        refinement_rule="normalized_edge_midpoint_four_child",
        source_incidence_sha256=_reference_incidence_sha256(),
    )
    return EchosahedralFederation(
        federation_id="source-derived-incidence-nerve-v1",
        carriers=carriers,
        seams=seams,
        external_boundaries=tuple(external_boundaries),
        observer_supports=(observer,),
        triple_overlaps=triple_overlaps,
        support_tower_binding=support_binding,
    )


def echosahedral_carrier_conformance_report(
    carrier: EchosahedralCarrier,
    *,
    tolerance: float = 5.0e-11,
) -> dict[str, Any]:
    """Return an isolated copy of an ID-independent cached local audit.

    Carrier identity is not part of local icosahedral conformance.  Normalizing
    the ID before cache lookup makes 256k canonical carriers share one exact
    local audit instead of populating a cache with ID-distinguished copies.
    """

    try:
        if _uses_shared_reference_template(carrier):
            report = copy.deepcopy(_shared_reference_conformance_cached(tolerance))
        else:
            normalized = replace(carrier, carrier_id=_REFERENCE_CARRIER_TEMPLATE_ID)
            report = copy.deepcopy(
                _echosahedral_carrier_conformance_cached(normalized, tolerance)
            )
    except (IndexError, OverflowError, TypeError, ValueError) as exc:
        report = _malformed_carrier_conformance_report(carrier, exc)
    identifier_valid = bool(isinstance(carrier.carrier_id, str) and carrier.carrier_id)
    report["carrier_id"] = carrier.carrier_id
    report["carrier_identifier_valid"] = identifier_valid
    if not identifier_valid:
        report["blockers"] = sorted(
            set(report["blockers"]) | {"carrier_id_must_be_nonempty_string"}
        )
        report["ECHOSAHEDRAL_CARRIER_CONFORMANCE"] = False
        report["ECHOSAHEDRAL_CARRIER_CONFORMANCE_RECEIPT"] = False
    return report


def _malformed_carrier_conformance_report(
    carrier: EchosahedralCarrier, exc: Exception
) -> dict[str, Any]:
    """Emit a stable fail-closed row when a typed carrier is malformed at runtime."""

    blocker = f"malformed_carrier_presentation:{type(exc).__name__}"
    return {
        "schema": "oph.echosahedral_carrier.conformance.v1",
        "carrier_id": carrier.carrier_id,
        "carrier_identifier_valid": bool(
            isinstance(carrier.carrier_id, str) and carrier.carrier_id
        ),
        "port_count": None,
        "edge_count": None,
        "face_count": None,
        "hidden_presentation_coordinates": True,
        "hidden_coordinates_eligible_for_promoted_geometry": False,
        "structural_class_sha256": interface_algebra_sha256(
            {
                "schema": "oph.echosahedral_carrier.malformed.v1",
                "exception_type": type(exc).__name__,
            }
        ),
        "blockers": [blocker],
        "ECHOSAHEDRAL_CARRIER_CONFORMANCE": False,
        "ECHOSAHEDRAL_CARRIER_CONFORMANCE_RECEIPT": False,
        "LOCAL_RESPONSE_1_3_3PRIME_5_DECOMPOSITION_RECEIPT": False,
        "claim_boundary": (
            "Malformed runtime presentation data cannot earn local conformance."
        ),
    }


@lru_cache(maxsize=512)
def _echosahedral_carrier_conformance_cached(
    carrier: EchosahedralCarrier,
    tolerance: float = 5.0e-11,
) -> dict[str, Any]:
    """Recompute complete local 12/30/20, antipode, and A5 conformance."""

    blockers: list[str] = []
    coordinates = np.asarray(carrier.port_coordinates, dtype=float)
    edges = tuple(tuple(int(value) for value in edge) for edge in carrier.edges)
    faces = tuple(tuple(int(value) for value in face) for face in carrier.faces)
    antipode = tuple(int(value) for value in carrier.antipode)
    actions = tuple(
        tuple(int(value) for value in action) for action in carrier.a5_actions
    )
    if (
        len(carrier.port_names) != _PORT_COUNT
        or len(set(carrier.port_names)) != _PORT_COUNT
    ):
        blockers.append("port_names_are_not_twelve_unique_hidden_labels")
    if coordinates.shape != (_PORT_COUNT, 3) or not np.all(np.isfinite(coordinates)):
        blockers.append("hidden_coordinates_are_not_finite_12_by_3")
        coordinates = np.zeros((_PORT_COUNT, 3), dtype=float)
    coordinate_norm_residual = float(
        np.max(np.abs(np.linalg.norm(coordinates, axis=1) - 1.0))
    )
    if coordinate_norm_residual > tolerance:
        blockers.append("port_coordinates_are_not_on_unit_reference_sphere")

    edge_set: set[tuple[int, int]] = set()
    edge_shape_valid = len(edges) == 30
    for edge in edges:
        if (
            len(edge) != 2
            or edge[0] == edge[1]
            or min(edge) < 0
            or max(edge) >= _PORT_COUNT
        ):
            edge_shape_valid = False
            continue
        edge_set.add(tuple(sorted(edge)))
    edge_shape_valid = edge_shape_valid and len(edge_set) == 30
    if not edge_shape_valid:
        blockers.append("edge_incidence_is_not_30_distinct_simple_edges")
    degrees = np.zeros(_PORT_COUNT, dtype=int)
    for left, right in edge_set:
        degrees[left] += 1
        degrees[right] += 1
    degree_five = bool(np.all(degrees == 5))
    if not degree_five:
        blockers.append("not_every_port_has_local_degree_five")

    face_shape_valid = len(faces) == 20
    oriented_face_set: set[tuple[int, int, int]] = set()
    unoriented_face_set: set[tuple[int, int, int]] = set()
    outward_residual = math.inf
    outward_values: list[float] = []
    for face in faces:
        if (
            len(face) != 3
            or len(set(face)) != 3
            or min(face) < 0
            or max(face) >= _PORT_COUNT
        ):
            face_shape_valid = False
            continue
        cyclic = _cyclic_face_key(face)
        unoriented = tuple(sorted(face))
        oriented_face_set.add(cyclic)
        unoriented_face_set.add(unoriented)
        a, b, c = face
        outward_values.append(
            float(
                np.dot(
                    np.cross(
                        coordinates[b] - coordinates[a], coordinates[c] - coordinates[a]
                    ),
                    coordinates[a] + coordinates[b] + coordinates[c],
                )
            )
        )
    face_shape_valid = bool(
        face_shape_valid
        and len(oriented_face_set) == 20
        and len(unoriented_face_set) == 20
    )
    if outward_values:
        outward_residual = float(min(outward_values))
    outward_oriented = bool(face_shape_valid and outward_residual > tolerance)
    if not face_shape_valid:
        blockers.append("face_incidence_is_not_20_distinct_triangles")
    if not outward_oriented:
        blockers.append("faces_are_not_consistently_outward_oriented")

    face_edges_valid = True
    edge_face_counts: dict[tuple[int, int], int] = {edge: 0 for edge in edge_set}
    for face in faces:
        if len(face) != 3:
            face_edges_valid = False
            continue
        for pair in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            key = tuple(sorted(pair))
            if key not in edge_set:
                face_edges_valid = False
            else:
                edge_face_counts[key] += 1
    closed_surface = bool(
        face_edges_valid
        and len(edge_face_counts) == 30
        and all(count == 2 for count in edge_face_counts.values())
        and _PORT_COUNT - len(edge_set) + len(unoriented_face_set) == 2
    )
    if not closed_surface:
        blockers.append("edge_face_incidence_is_not_a_closed_euler_two_surface")

    antipode_valid = bool(
        len(antipode) == _PORT_COUNT
        and sorted(antipode) == list(range(_PORT_COUNT))
        and all(antipode[antipode[index]] == index for index in range(_PORT_COUNT))
        and all(antipode[index] != index for index in range(_PORT_COUNT))
    )
    antipodal_residual = math.inf
    if antipode_valid:
        antipodal_residual = float(
            np.max(
                np.linalg.norm(
                    coordinates + coordinates[np.asarray(antipode, dtype=int)], axis=1
                )
            )
        )
        antipode_valid = antipodal_residual <= tolerance
    if not antipode_valid:
        blockers.append("antipode_is_not_a_fixed_point_free_geometric_involution")

    a5 = _a5_local_action_audit(
        actions,
        edge_set=edge_set,
        oriented_face_set=oriented_face_set,
        antipode=antipode,
        coordinates=coordinates,
        tolerance=tolerance,
    )
    if not a5["receipt"]:
        blockers.extend(a5["blockers"])

    reference_matches = _reference_icosahedron_isomorphism_receipt(
        coordinates,
        edge_set=edge_set,
        oriented_face_set=oriented_face_set,
        antipode=antipode,
        tolerance=tolerance,
    )
    if not reference_matches["receipt"]:
        blockers.append("carrier_is_not_an_exact_relabeling_of_reference_icosahedron")

    adjacency = np.zeros((_PORT_COUNT, _PORT_COUNT), dtype=float)
    for left, right in edge_set:
        adjacency[left, right] = 1.0
        adjacency[right, left] = 1.0
    spectrum = np.linalg.eigvalsh(adjacency)
    multiplicities = _eigenvalue_multiplicities(spectrum, tolerance=1.0e-8)
    irrep_multiplicities = sorted(multiplicity for _, multiplicity in multiplicities)
    response_decomposition = irrep_multiplicities == [1, 3, 3, 5]
    if not response_decomposition:
        blockers.append("adjacency_response_does_not_have_1_3_3prime_5_sectors")

    reference = echosahedral_patch_architecture_report(12)
    passed = not blockers
    structural_payload = {
        "schema": "oph.echosahedral_carrier.structural_class.v1",
        "port_count": 12,
        "edge_count": 30,
        "face_count": 20,
        "degree_profile": [5] * 12,
        "antipodal_pair_count": 6,
        "a5_order": 60,
        "a5_order_profile": {"1": 1, "2": 15, "3": 20, "5": 24},
        "response_multiplicities": [1, 3, 3, 5],
    }
    return {
        "schema": "oph.echosahedral_carrier.conformance.v1",
        "carrier_id": carrier.carrier_id,
        "port_count": len(carrier.port_names),
        "edge_count": len(edge_set),
        "face_count": len(oriented_face_set),
        "vertex_degree_profile": {
            str(value): int(np.sum(degrees == value)) for value in sorted(set(degrees))
        },
        "maximum_unit_norm_residual": coordinate_norm_residual,
        "minimum_outward_face_orientation_witness": outward_residual,
        "maximum_antipodal_residual": antipodal_residual,
        "antipode_fixed_point_count": sum(
            index == partner for index, partner in enumerate(antipode)
        )
        if len(antipode) == _PORT_COUNT
        else None,
        "closed_euler_two_surface": closed_surface,
        "a5_action": a5,
        "reference_isomorphism": reference_matches,
        "adjacency_spectrum": [float(value) for value in spectrum],
        "adjacency_eigenvalue_multiplicities": [
            {"eigenvalue": value, "multiplicity": multiplicity}
            for value, multiplicity in multiplicities
        ],
        "declared_response_sector_dimensions": [1, 3, 3, 5],
        "hidden_presentation_coordinates": True,
        "hidden_coordinates_eligible_for_promoted_geometry": False,
        "shared_reference_template_hash": reference.get("template_hash"),
        "structural_class_sha256": interface_algebra_sha256(structural_payload),
        "blockers": sorted(set(blockers)),
        "ECHOSAHEDRAL_CARRIER_CONFORMANCE": passed,
        "ECHOSAHEDRAL_CARRIER_CONFORMANCE_RECEIPT": passed,
        "LOCAL_RESPONSE_1_3_3PRIME_5_DECOMPOSITION_RECEIPT": response_decomposition,
        "claim_boundary": (
            "This certifies a hidden finite carrier presentation only. The local XYZ "
            "coordinates and port frame are ineligible as S2, H3, event, cap-normal, "
            "clock, or BW source coordinates."
        ),
    }


@lru_cache(maxsize=16)
def _shared_reference_conformance_cached(tolerance: float) -> dict[str, Any]:
    """Cache the canonical audit by tolerance without rehashing its large tuples."""

    return _echosahedral_carrier_conformance_cached(
        _reference_echosahedral_carrier_template(), tolerance
    )


class _BoundedBlockerLedger:
    """Count every blocker while retaining only bounded distinct examples."""

    def __init__(self, limit: int = _REPORT_DETAIL_LIMIT) -> None:
        self.limit = limit
        self.occurrence_count = 0
        self._examples: list[str] = []
        self._seen_examples: set[str] = set()

    def add(self, blocker: str) -> None:
        self.occurrence_count += 1
        if blocker not in self._seen_examples and len(self._examples) < self.limit:
            self._seen_examples.add(blocker)
            self._examples.append(blocker)

    def extend(self, blockers: Sequence[str]) -> None:
        for blocker in blockers:
            self.add(str(blocker))

    @property
    def examples(self) -> list[str]:
        return sorted(self._examples)

    @property
    def truncated(self) -> bool:
        return self.occurrence_count > len(self._examples)


class _CanonicalRowDigest:
    """Incrementally bind all verified rows without retaining them in a report."""

    def __init__(self) -> None:
        self._hasher = hashlib.sha256()
        self.row_count = 0

    def update(self, row: Mapping[str, Any]) -> None:
        encoded = json.dumps(
            row, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode("utf-8")
        self._hasher.update(len(encoded).to_bytes(8, "big"))
        self._hasher.update(encoded)
        self.row_count += 1

    def hexdigest(self) -> str:
        return "sha256:" + self._hasher.hexdigest()


class _DisjointSet:
    """Small streaming connectivity audit for the carrier seam graph."""

    def __init__(self, nodes: Iterable[str]) -> None:
        self.parent = {node: node for node in nodes}
        self.rank = {node: 0 for node in nodes}
        self.component_count = len(self.parent)

    def find(self, node: str) -> str:
        parent = self.parent[node]
        while parent != self.parent[parent]:
            self.parent[parent] = self.parent[self.parent[parent]]
            parent = self.parent[parent]
        while node != parent:
            next_node = self.parent[node]
            self.parent[node] = parent
            node = next_node
        return parent

    def union(self, left: str, right: str) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        if self.rank[left_root] < self.rank[right_root]:
            left_root, right_root = right_root, left_root
        self.parent[right_root] = left_root
        if self.rank[left_root] == self.rank[right_root]:
            self.rank[left_root] += 1
        self.component_count -= 1


def _carrier_structural_key(carrier: EchosahedralCarrier) -> EchosahedralCarrier:
    """Return an exact hashable presentation key with identity removed."""

    return replace(carrier, carrier_id=_REFERENCE_CARRIER_TEMPLATE_ID)


def _uses_shared_reference_template(carrier: EchosahedralCarrier) -> bool:
    """Recognize the O(1) canonical shared-template representation."""

    template = _reference_echosahedral_carrier_template()
    return bool(
        carrier.port_names is template.port_names
        and carrier.port_coordinates is template.port_coordinates
        and carrier.edges is template.edges
        and carrier.faces is template.faces
        and carrier.antipode is template.antipode
        and carrier.a5_actions is template.a5_actions
    )


def federation_sewing_report(
    federation: EchosahedralFederation,
) -> dict[str, Any]:
    """Validate a finite federation with bounded, content-bound diagnostics.

    Validation still visits every carrier, seam, boundary and observer support,
    but repeated carrier presentations are audited once and output rows are
    bounded.  Twelve-bit occupancy masks replace the former 12N endpoint set.
    Consequently a valid 256k-carrier input does not create a 256k-entry local
    report or a three-million-tuple ``all_ports`` allocation.
    """

    blockers = _BoundedBlockerLedger()
    carrier_by_id: dict[str, EchosahedralCarrier] = {}
    carrier_reports: dict[str, dict[str, Any]] = {}
    structure_audits: dict[Any, dict[str, Any]] = {}
    structure_counts: dict[Any, int] = {}
    structure_first_ids: dict[Any, str] = {}
    reference_structure = _REFERENCE_STRUCTURE_TOKEN
    shared_reference_template_carrier_count = 0
    every_carrier_conforms = True

    for carrier_index, carrier in enumerate(federation.carriers):
        identifier_valid = bool(
            isinstance(carrier.carrier_id, str) and carrier.carrier_id
        )
        display_id = str(carrier.carrier_id)
        if not identifier_valid:
            blockers.add(f"invalid_carrier_id:{display_id}")
        elif carrier.carrier_id in carrier_by_id:
            blockers.add(f"duplicate_carrier_id:{carrier.carrier_id}")
        else:
            carrier_by_id[carrier.carrier_id] = carrier

        structure: Any
        if _uses_shared_reference_template(carrier):
            structure = _REFERENCE_STRUCTURE_TOKEN
            audit_carrier = _reference_echosahedral_carrier_template()
        else:
            audit_carrier = _carrier_structural_key(carrier)
            try:
                hash(audit_carrier)
                structure = audit_carrier
            except TypeError:
                structure = ("malformed_unhashable_presentation", carrier_index)
        if structure not in structure_audits:
            structure_audits[structure] = echosahedral_carrier_conformance_report(
                audit_carrier
            )
            structure_counts[structure] = 0
            structure_first_ids[structure] = display_id
        structure_counts[structure] += 1
        if structure == reference_structure:
            shared_reference_template_carrier_count += 1
        structural_report = structure_audits[structure]
        carrier_conforms = bool(
            identifier_valid and structural_report["ECHOSAHEDRAL_CARRIER_CONFORMANCE"]
        )
        every_carrier_conforms = every_carrier_conforms and carrier_conforms
        if not carrier_conforms:
            for item in structural_report["blockers"] or ["carrier_identifier_invalid"]:
                blockers.add(f"carrier:{display_id}:{item}")

        if len(carrier_reports) < _REPORT_DETAIL_LIMIT:
            example = copy.deepcopy(structural_report)
            example["carrier_id"] = carrier.carrier_id
            example["carrier_identifier_valid"] = identifier_valid
            if not identifier_valid:
                example["blockers"] = sorted(
                    set(example["blockers"]) | {"carrier_id_must_be_nonempty_string"}
                )
                example["ECHOSAHEDRAL_CARRIER_CONFORMANCE"] = False
                example["ECHOSAHEDRAL_CARRIER_CONFORMANCE_RECEIPT"] = False
            report_key = (
                carrier.carrier_id
                if identifier_valid and carrier.carrier_id not in carrier_reports
                else f"entry:{len(carrier_reports)}:{display_id}"
            )
            carrier_reports[report_key] = example

    if not federation.carriers:
        blockers.add("federation_has_no_carriers")
        every_carrier_conforms = False
    exact_carrier_cardinality = bool(
        federation.carriers and len(carrier_by_id) == len(federation.carriers)
    )
    carrier_conformance = bool(exact_carrier_cardinality and every_carrier_conforms)
    structure_examples: list[dict[str, Any]] = []
    for structure, report in structure_audits.items():
        if len(structure_examples) == _REPORT_DETAIL_LIMIT:
            break
        structure_examples.append(
            {
                "first_carrier_id": structure_first_ids[structure],
                "carrier_count": structure_counts[structure],
                "structural_class_sha256": report["structural_class_sha256"],
                "conforming": report["ECHOSAHEDRAL_CARRIER_CONFORMANCE"],
                "blockers": report["blockers"],
            }
        )
    reference_audit = structure_audits.get(reference_structure)
    carrier_summary = {
        "schema": "oph.echosahedral_federation.carrier_conformance_summary.v1",
        "verification_mode": "exact_structural_deduplication_with_bounded_examples",
        "carrier_entry_count": len(federation.carriers),
        "unique_carrier_id_count": len(carrier_by_id),
        "unique_structural_presentation_count": len(structure_audits),
        "local_conformance_audit_count": len(structure_audits),
        "shared_reference_template_carrier_count": (
            shared_reference_template_carrier_count
        ),
        "shared_reference_template_conformance_verified_once": bool(
            shared_reference_template_carrier_count
            and reference_audit is not None
            and reference_audit["ECHOSAHEDRAL_CARRIER_CONFORMANCE"]
        ),
        "carrier_report_example_count": len(carrier_reports),
        "carrier_reports_truncated": (len(federation.carriers) > len(carrier_reports)),
        "structural_class_examples": structure_examples,
        "structural_class_examples_truncated": (
            len(structure_audits) > len(structure_examples)
        ),
        "all_carriers_conform": carrier_conformance,
    }

    seam_rows: list[dict[str, Any]] = []
    seam_failure_examples: list[dict[str, Any]] = []
    seam_digest = _CanonicalRowDigest()
    all_seams_pass = True
    all_interface_schema_hashes_agree = bool(federation.seams)
    all_interface_algebra_maps_pass = bool(federation.seams)
    occupied_masks: dict[str, int] = {}
    unknown_port_examples: list[dict[str, Any]] = []
    unknown_port_reference_count = 0
    connectivity = _DisjointSet(carrier_by_id)
    seam_by_id: dict[str, SeamBundle] = {}
    for seam in federation.seams:
        row = _seam_bundle_report(seam, carrier_by_id)
        all_interface_schema_hashes_agree = bool(
            all_interface_schema_hashes_agree
            and row["endpoint_interface_algebra_hashes_agree"]
        )
        all_interface_algebra_maps_pass = bool(
            all_interface_algebra_maps_pass
            and row["INTERFACE_ALGEBRA_MAP_HOMOMORPHISM_RECEIPT"]
        )
        seam_digest.update(row)
        if len(seam_rows) < _REPORT_DETAIL_LIMIT:
            seam_rows.append(row)
        if not row["SEAM_BUNDLE_RECEIPT"]:
            all_seams_pass = False
            if len(seam_failure_examples) < _REPORT_DETAIL_LIMIT:
                seam_failure_examples.append(row)
        if not seam.seam_id or seam.seam_id in seam_by_id:
            blockers.add(f"duplicate_or_empty_seam_id:{seam.seam_id}")
        seam_by_id.setdefault(seam.seam_id, seam)
        if not row["SEAM_BUNDLE_RECEIPT"]:
            blockers.extend([f"seam:{seam.seam_id}:{item}" for item in row["blockers"]])
        if (
            seam.left_carrier_id in carrier_by_id
            and seam.right_carrier_id in carrier_by_id
        ):
            connectivity.union(seam.left_carrier_id, seam.right_carrier_id)
        for carrier_id, ports in (
            (seam.left_carrier_id, row["left_ports"]),
            (seam.right_carrier_id, row["right_ports"]),
        ):
            for raw_port in ports:
                if type(raw_port) is not int:
                    unknown_port_reference_count += 1
                    if len(unknown_port_examples) < _REPORT_DETAIL_LIMIT:
                        unknown_port_examples.append(
                            {
                                "carrier_id": str(carrier_id),
                                "port": repr(raw_port),
                                "source": "seam_noninteger",
                            }
                        )
                    continue
                port = raw_port
                if carrier_id not in carrier_by_id or not 0 <= port < _PORT_COUNT:
                    unknown_port_reference_count += 1
                    if len(unknown_port_examples) < _REPORT_DETAIL_LIMIT:
                        unknown_port_examples.append(
                            {"carrier_id": carrier_id, "port": port, "source": "seam"}
                        )
                    continue
                bit = 1 << port
                prior = occupied_masks.get(carrier_id, 0)
                if prior & bit:
                    blockers.add(
                        f"local_port_used_by_more_than_one_seam:{carrier_id}:P{port}"
                    )
                occupied_masks[carrier_id] = prior | bit

    triple_overlap_ids: set[str] = set()
    triple_overlap_rows: list[dict[str, Any]] = []
    triple_overlap_all_rows: list[dict[str, Any]] = []
    triple_overlap_failure_examples: list[dict[str, Any]] = []
    triple_overlap_digest = _CanonicalRowDigest()
    all_triple_overlaps_pass = True
    for overlap in federation.triple_overlaps:
        row = _triple_overlap_report(overlap, carrier_by_id, seam_by_id)
        triple_overlap_all_rows.append(row)
        triple_overlap_digest.update(row)
        if len(triple_overlap_rows) < _REPORT_DETAIL_LIMIT:
            triple_overlap_rows.append(row)
        if not overlap.overlap_id or overlap.overlap_id in triple_overlap_ids:
            blockers.add(f"duplicate_or_empty_triple_overlap_id:{overlap.overlap_id}")
            all_triple_overlaps_pass = False
        triple_overlap_ids.add(overlap.overlap_id)
        if not row["TRIPLE_OVERLAP_RESTRICTION_RECEIPT"]:
            all_triple_overlaps_pass = False
            if len(triple_overlap_failure_examples) < _REPORT_DETAIL_LIMIT:
                triple_overlap_failure_examples.append(row)
            blockers.extend(
                [
                    f"triple_overlap:{overlap.overlap_id}:{item}"
                    for item in row["blockers"]
                ]
            )

    boundary_ids: set[str] = set()
    boundary_rows: list[dict[str, Any]] = []
    boundary_failure_examples: list[dict[str, Any]] = []
    boundary_digest = _CanonicalRowDigest()
    all_boundaries_pass = True
    declared_external_masks: dict[str, int] = {}
    for boundary in federation.external_boundaries:
        row = _external_boundary_report(boundary, carrier_by_id)
        boundary_digest.update(row)
        if len(boundary_rows) < _REPORT_DETAIL_LIMIT:
            boundary_rows.append(row)
        if not row["EXPLICIT_EXTERNAL_BOUNDARY_RECEIPT"]:
            all_boundaries_pass = False
            if len(boundary_failure_examples) < _REPORT_DETAIL_LIMIT:
                boundary_failure_examples.append(row)
        if not boundary.boundary_id or boundary.boundary_id in boundary_ids:
            blockers.add(
                f"duplicate_or_empty_external_boundary_id:{boundary.boundary_id}"
            )
        boundary_ids.add(boundary.boundary_id)
        if not row["EXPLICIT_EXTERNAL_BOUNDARY_RECEIPT"]:
            blockers.extend(
                [f"boundary:{boundary.boundary_id}:{item}" for item in row["blockers"]]
            )
        for raw_port in row["ports"]:
            if type(raw_port) is not int:
                unknown_port_reference_count += 1
                if len(unknown_port_examples) < _REPORT_DETAIL_LIMIT:
                    unknown_port_examples.append(
                        {
                            "carrier_id": str(boundary.carrier_id),
                            "port": repr(raw_port),
                            "source": "boundary_noninteger",
                        }
                    )
                continue
            port = raw_port
            carrier_id = boundary.carrier_id
            if carrier_id not in carrier_by_id or not 0 <= port < _PORT_COUNT:
                unknown_port_reference_count += 1
                if len(unknown_port_examples) < _REPORT_DETAIL_LIMIT:
                    unknown_port_examples.append(
                        {"carrier_id": carrier_id, "port": port, "source": "boundary"}
                    )
                continue
            bit = 1 << port
            prior = declared_external_masks.get(carrier_id, 0)
            if prior & bit:
                blockers.add(
                    f"local_port_declared_external_more_than_once:{carrier_id}:P{port}"
                )
            if occupied_masks.get(carrier_id, 0) & bit:
                blockers.add(f"sewn_port_also_declared_external:{carrier_id}:P{port}")
            declared_external_masks[carrier_id] = prior | bit

    full_mask = (1 << _PORT_COUNT) - 1
    dangling_port_count = 0
    dangling_port_examples: list[dict[str, Any]] = []
    for carrier_id in carrier_by_id:
        declared = occupied_masks.get(carrier_id, 0) | declared_external_masks.get(
            carrier_id, 0
        )
        missing = full_mask & ~declared
        dangling_port_count += missing.bit_count()
        if missing and len(dangling_port_examples) < _REPORT_DETAIL_LIMIT:
            for port in range(_PORT_COUNT):
                if missing & (1 << port):
                    dangling_port_examples.append(
                        {"carrier_id": carrier_id, "port": port}
                    )
                    if len(dangling_port_examples) == _REPORT_DETAIL_LIMIT:
                        break
    if dangling_port_count:
        blockers.add("dangling_ports_lack_explicit_external_boundary_declaration")
    if unknown_port_reference_count:
        blockers.add("seam_or_boundary_references_unknown_local_port")

    federation_connected = bool(
        len(carrier_by_id) == 1
        or (len(carrier_by_id) > 1 and connectivity.component_count == 1)
    )
    if not federation_connected:
        blockers.add("intercarrier_seam_graph_is_disconnected")

    observer_rows: list[dict[str, Any]] = []
    observer_failure_examples: list[dict[str, Any]] = []
    observer_digest = _CanonicalRowDigest()
    all_observers_pass = True
    observer_tokens: set[str] = set()
    for support in federation.observer_supports:
        row = observer_support_report(support, carrier_by_id, seam_by_id)
        observer_digest.update(row)
        if len(observer_rows) < _REPORT_DETAIL_LIMIT:
            observer_rows.append(row)
        if support.observer_token in observer_tokens:
            blockers.add("duplicate_observer_token")
        observer_tokens.add(support.observer_token)
        if not row["CONNECTED_OBSERVER_SUPPORT_RECEIPT"]:
            all_observers_pass = False
            if len(observer_failure_examples) < _REPORT_DETAIL_LIMIT:
                observer_failure_examples.append(row)
            blockers.extend(
                [f"observer:{row['observer_token']}:{item}" for item in row["blockers"]]
            )

    passed = bool(
        carrier_conformance
        and blockers.occurrence_count == 0
        and all_seams_pass
        and all_triple_overlaps_pass
        and all_boundaries_pass
        and all_observers_pass
    )
    cocycle = _higher_overlap_cocycle_report(
        federation.seams,
        carrier_by_id,
        explicit_triple_rows=triple_overlap_all_rows,
    )
    cocycle_receipt = bool(cocycle["HIGHER_OVERLAP_COCYCLE_RECEIPT"])
    nonvacuous_cocycle_witness = bool(
        cocycle["NONVACUOUS_HIGHER_OVERLAP_COCYCLE_WITNESS"]
    )
    full_interface_sewing = bool(
        all_interface_schema_hashes_agree
        and all_interface_algebra_maps_pass
        and cocycle_receipt
    )
    physical_realization = bool(
        passed
        and full_interface_sewing
        and nonvacuous_cocycle_witness
        and federation.observer_supports
    )
    sewn_port_count = sum(mask.bit_count() for mask in occupied_masks.values())
    external_port_count = sum(
        mask.bit_count() for mask in declared_external_masks.values()
    )
    return {
        "schema": "oph.echosahedral_federation.sewing.v1",
        "federation_id": federation.federation_id,
        "cardinality_semantics": _CARDINALITY_SEMANTICS,
        "carrier_count": len(federation.carriers),
        "exact_source_carrier_count": (
            len(federation.carriers) if exact_carrier_cardinality else None
        ),
        "carrier_count_is_exact_declared_federation_cardinality": (
            exact_carrier_cardinality
        ),
        "support_regulator_count": None,
        "carrier_count_is_support_regulator_count": False,
        "carrier_count_is_support_chart_cell_count": False,
        "carrier_count_is_screen_entropy_capacity_N_star": False,
        "carrier_count_is_primitive_observer_count": False,
        "seam_count": len(federation.seams),
        "declared_triple_overlap_count": len(federation.triple_overlaps),
        "external_boundary_bundle_count": len(federation.external_boundaries),
        "observer_support_count": len(federation.observer_supports),
        "sewn_local_port_count": sewn_port_count,
        "declared_external_local_port_count": external_port_count,
        "undeclared_dangling_port_count": dangling_port_count,
        "undeclared_dangling_ports": dangling_port_examples,
        "undeclared_dangling_port_examples_truncated": (
            dangling_port_count > len(dangling_port_examples)
        ),
        "unknown_local_port_reference_count": unknown_port_reference_count,
        "unknown_local_port_reference_examples": unknown_port_examples,
        "carrier_conformance": carrier_reports,
        "carrier_conformance_summary": carrier_summary,
        "seams": seam_rows,
        "seam_rows_sha256": seam_digest.hexdigest(),
        "seam_rows_reported_count": len(seam_rows),
        "seam_rows_truncated": len(federation.seams) > len(seam_rows),
        "seam_failure_examples": seam_failure_examples,
        "triple_overlaps": triple_overlap_rows,
        "triple_overlap_rows_sha256": triple_overlap_digest.hexdigest(),
        "triple_overlap_rows_reported_count": len(triple_overlap_rows),
        "triple_overlap_rows_truncated": (
            len(federation.triple_overlaps) > len(triple_overlap_rows)
        ),
        "triple_overlap_failure_examples": triple_overlap_failure_examples,
        "external_boundaries": boundary_rows,
        "external_boundary_rows_sha256": boundary_digest.hexdigest(),
        "external_boundary_rows_reported_count": len(boundary_rows),
        "external_boundary_rows_truncated": (
            len(federation.external_boundaries) > len(boundary_rows)
        ),
        "external_boundary_failure_examples": boundary_failure_examples,
        "observer_supports": observer_rows,
        "observer_support_rows_sha256": observer_digest.hexdigest(),
        "observer_support_rows_reported_count": len(observer_rows),
        "observer_support_rows_truncated": (
            len(federation.observer_supports) > len(observer_rows)
        ),
        "observer_support_failure_examples": observer_failure_examples,
        "report_detail_limit": _REPORT_DETAIL_LIMIT,
        "federation_carrier_graph_component_count": connectivity.component_count,
        "federation_carrier_graph_connected": federation_connected,
        "blocker_occurrence_count": blockers.occurrence_count,
        "blockers": blockers.examples,
        "blocker_examples_truncated": blockers.truncated,
        "higher_overlap_cocycle": cocycle,
        "ECHOSAHEDRAL_CARRIER_CONFORMANCE": carrier_conformance,
        "STRUCTURAL_FEDERATION_SEWING_RECEIPT": passed,
        "FEDERATION_SEWING_RECEIPT": passed,
        "INTERFACE_SCHEMA_HASH_BINDING_RECEIPT": (all_interface_schema_hashes_agree),
        "INTERFACE_ALGEBRA_MAP_HOMOMORPHISM_RECEIPT": all_interface_algebra_maps_pass,
        "HIGHER_OVERLAP_COCYCLE_RECEIPT": cocycle_receipt,
        "HIGHER_OVERLAP_COCYCLE_CONDITION_RECEIPT": cocycle_receipt,
        "NONVACUOUS_HIGHER_OVERLAP_COCYCLE_WITNESS": (
            nonvacuous_cocycle_witness
        ),
        "CONNECTED_OBSERVER_SUPPORT_WITNESS": bool(
            federation.observer_supports and all_observers_pass
        ),
        "FULL_INTERFACE_ALGEBRA_SEWING_RECEIPT": full_interface_sewing,
        "PHYSICAL_ECHOSAHEDRAL_FEDERATION_REALIZATION_RECEIPT": physical_realization,
        "CARRIER_TO_SUPPORT_CHART_REALIZATION_RECEIPT": False,
        "claim_boundary": (
            "carrier_count is the exact cardinality of the declared finite source "
            "federation when the ID receipt passes; it is separate from every support "
            "regulator, S2 chart-cell count, N_star, observer count, H3-point count, "
            "and event count. Sewing checks port bijections, orientations, boundary "
            "coverage, and content-addressed schema identity. The interface algebra "
            "map and higher-overlap cocycle condition are computed from exact finite "
            "checks on the typed seams. A separate nonvacuous witness records whether "
            "the declared cover contains a composable seam triangle; absence of one "
            "does not falsify the logical Cech condition, but it cannot promote "
            "the physical federation realization receipt: that receipt also "
            "requires a nonvacuous seam-triangle witness. None of these "
            "receipts implies emergent geometry, which requires independent "
            "downstream receipts."
        ),
    }


def observer_support_report(
    support: ObserverSupport,
    carrier_by_id: Mapping[str, EchosahedralCarrier],
    seam_by_id: Mapping[str, SeamBundle],
) -> dict[str, Any]:
    """Check that an observer token has connected visible carrier support."""

    blockers: list[str] = []
    if not support.observer_token:
        blockers.append("empty_observer_token")
    if not support.carrier_ids:
        blockers.append("observer_support_has_no_carriers")
    unknown_carriers = sorted(set(support.carrier_ids) - set(carrier_by_id))
    if unknown_carriers:
        blockers.append("observer_support_references_unknown_carrier")
    unknown_seams = sorted(set(support.visible_seam_ids) - set(seam_by_id))
    if unknown_seams:
        blockers.append("observer_support_references_unknown_visible_seam")
    if not _is_sha256(support.record_algebra_sha256):
        blockers.append("invalid_record_algebra_sha256")
    if not _is_sha256(support.checkpoint_cut_sha256):
        blockers.append("invalid_checkpoint_cut_sha256")

    support_carrier_ids = set(support.carrier_ids)
    connectivity = _DisjointSet(support_carrier_ids)
    for seam_id in support.visible_seam_ids:
        seam = seam_by_id.get(seam_id)
        if seam is None:
            continue
        endpoints = {seam.left_carrier_id, seam.right_carrier_id}
        if not endpoints <= support_carrier_ids:
            blockers.append("visible_seam_leaves_observer_carrier_support")
            continue
        connectivity.union(seam.left_carrier_id, seam.right_carrier_id)
    connected = bool(
        len(support_carrier_ids) == 1
        or (len(support_carrier_ids) > 1 and connectivity.component_count == 1)
    )
    if not connected:
        blockers.append("observer_carrier_support_is_disconnected")
    passed = not blockers
    sorted_carrier_ids = sorted(support.carrier_ids)
    sorted_seam_ids = sorted(support.visible_seam_ids)
    carrier_digest = _CanonicalRowDigest()
    for carrier_id in sorted_carrier_ids:
        carrier_digest.update({"carrier_id": carrier_id})
    seam_digest = _CanonicalRowDigest()
    for seam_id in sorted_seam_ids:
        seam_digest.update({"seam_id": seam_id})
    return {
        "schema": "oph.echosahedral_federation.observer_support.v1",
        "observer_token": support.observer_token,
        "carrier_ids": sorted_carrier_ids[:_REPORT_DETAIL_LIMIT],
        "carrier_ids_sha256": carrier_digest.hexdigest(),
        "carrier_ids_truncated": len(sorted_carrier_ids) > _REPORT_DETAIL_LIMIT,
        "visible_seam_ids": sorted_seam_ids[:_REPORT_DETAIL_LIMIT],
        "visible_seam_ids_sha256": seam_digest.hexdigest(),
        "visible_seam_ids_truncated": len(sorted_seam_ids) > _REPORT_DETAIL_LIMIT,
        "carrier_count": len(support.carrier_ids),
        "visible_seam_count": len(support.visible_seam_ids),
        "one_carrier_support_allowed": True,
        "connected": connected,
        "blockers": sorted(set(blockers)),
        "CONNECTED_OBSERVER_SUPPORT_RECEIPT": passed,
        "OBSERVER_EQUALS_ONE_CARRIER_ASSUMPTION": False,
    }


def _runtime_exact_int_tuple(value: Any) -> tuple[tuple[int, ...], bool]:
    """Reject Python bool/string coercions at the typed-runtime boundary."""

    if not isinstance(value, tuple) or not all(type(item) is int for item in value):
        return (), False
    return value, True


def _seam_bundle_report(
    seam: SeamBundle,
    carrier_by_id: Mapping[str, EchosahedralCarrier],
) -> dict[str, Any]:
    blockers: list[str] = []
    if not seam.seam_id:
        blockers.append("empty_seam_id")
    if seam.left_carrier_id == seam.right_carrier_id:
        blockers.append("seam_endpoints_must_be_distinct_carriers")
    left = carrier_by_id.get(seam.left_carrier_id)
    right = carrier_by_id.get(seam.right_carrier_id)
    if left is None:
        blockers.append("unknown_left_carrier")
    if right is None:
        blockers.append("unknown_right_carrier")
    left_ports, left_ports_exact = _runtime_exact_int_tuple(seam.left_ports)
    right_ports, right_ports_exact = _runtime_exact_int_tuple(seam.right_ports)
    forward, forward_exact = _runtime_exact_int_tuple(seam.left_to_right_ports)
    backward, backward_exact = _runtime_exact_int_tuple(seam.right_to_left_ports)
    forward_orientation, forward_orientation_exact = _runtime_exact_int_tuple(
        seam.left_to_right_orientation
    )
    backward_orientation, backward_orientation_exact = _runtime_exact_int_tuple(
        seam.right_to_left_orientation
    )
    exact_integer_arrays = bool(
        left_ports_exact
        and right_ports_exact
        and forward_exact
        and backward_exact
        and forward_orientation_exact
        and backward_orientation_exact
    )
    if not exact_integer_arrays:
        blockers.append("seam_port_and_orientation_arrays_must_be_exact_integer_tuples")
    bundle_sizes_valid = bool(
        left_ports
        and len(left_ports) == len(set(left_ports))
        and len(right_ports) == len(set(right_ports)) == len(left_ports)
        and len(forward) == len(left_ports)
        and len(backward) == len(right_ports)
    )
    if not bundle_sizes_valid:
        blockers.append("port_bundles_and_gluing_maps_do_not_have_equal_nonzero_size")
    port_ranges_valid = bool(
        all(
            0 <= port < _PORT_COUNT
            for port in left_ports + right_ports + forward + backward
        )
    )
    if not port_ranges_valid:
        blockers.append("seam_port_out_of_range")
    bijection_valid = bool(
        bundle_sizes_valid
        and set(forward) == set(right_ports)
        and set(backward) == set(left_ports)
    )
    inverse_composition = False
    if bijection_valid:
        forward_map = dict(zip(left_ports, forward, strict=True))
        backward_map = dict(zip(right_ports, backward, strict=True))
        inverse_composition = bool(
            all(backward_map[forward_map[port]] == port for port in left_ports)
            and all(forward_map[backward_map[port]] == port for port in right_ports)
        )
    if not bijection_valid:
        blockers.append("gluing_map_is_not_a_bijection_between_endpoint_bundles")
    if not inverse_composition:
        blockers.append("forward_reverse_gluing_composition_is_not_identity")

    orientation_valid = bool(
        exact_integer_arrays
        and len(forward_orientation) == len(left_ports)
        and len(backward_orientation) == len(right_ports)
        and all(sign == -1 for sign in forward_orientation)
        and all(sign == -1 for sign in backward_orientation)
    )
    orientation_inverse_composition = False
    if orientation_valid and bijection_valid:
        forward_sign = dict(zip(left_ports, forward_orientation, strict=True))
        backward_sign = dict(zip(right_ports, backward_orientation, strict=True))
        forward_map = dict(zip(left_ports, forward, strict=True))
        backward_map = dict(zip(right_ports, backward, strict=True))
        orientation_inverse_composition = bool(
            all(
                forward_sign[port] * backward_sign[forward_map[port]] == 1
                and backward_map[forward_map[port]] == port
                for port in left_ports
            )
        )
    if not orientation_valid:
        blockers.append("seam_orientation_maps_are_not_explicit_reversals")
    if not orientation_inverse_composition:
        blockers.append("forward_reverse_orientation_composition_is_not_identity")

    left_connected = bool(
        left is not None
        and (
            _port_subset_connected(left, left_ports)
            or (
                seam.collar_kind == "antipodal_pair"
                and _collar_kind_matches(left, left_ports, seam.collar_kind)
            )
        )
    )
    right_connected = bool(
        right is not None
        and (
            _port_subset_connected(right, right_ports)
            or (
                seam.collar_kind == "antipodal_pair"
                and _collar_kind_matches(right, right_ports, seam.collar_kind)
            )
        )
    )
    if not left_connected:
        blockers.append("left_port_bundle_is_not_connected_in_local_incidence")
    if not right_connected:
        blockers.append("right_port_bundle_is_not_connected_in_local_incidence")
    collar_type_valid = bool(
        left is not None
        and right is not None
        and _collar_kind_matches(left, left_ports, seam.collar_kind)
        and _collar_kind_matches(right, right_ports, seam.collar_kind)
    )
    if not collar_type_valid:
        blockers.append("declared_collar_kind_does_not_match_local_port_incidence")

    binding = seam.interface_algebra
    hashes_valid = all(
        _is_sha256(value)
        for value in (
            binding.interface_algebra_sha256,
            binding.left_interface_algebra_sha256,
            binding.right_interface_algebra_sha256,
        )
    )
    schema_hashes_agree = bool(
        hashes_valid
        and binding.interface_algebra_sha256
        == binding.left_interface_algebra_sha256
        == binding.right_interface_algebra_sha256
    )
    checked_schema: dict[str, Any] | None = None
    checked_schema_sha256: str | None = None
    schema_hash_binds_checked_algebra = False
    if binding.interface_algebra_id and left_ports:
        checked_schema = finite_matrix_interface_algebra_schema(
            binding.interface_algebra_id, len(left_ports)
        )
        checked_schema_sha256 = interface_algebra_sha256(checked_schema)
        schema_hash_binds_checked_algebra = bool(
            schema_hashes_agree
            and binding.interface_algebra_sha256 == checked_schema_sha256
        )
    if not binding.interface_algebra_id:
        blockers.append("empty_interface_algebra_id")
    if not hashes_valid:
        blockers.append("invalid_interface_algebra_sha256")
    if not schema_hashes_agree:
        blockers.append("endpoint_interface_algebra_hashes_do_not_agree")
    if not schema_hash_binds_checked_algebra:
        blockers.append("interface_algebra_hash_does_not_bind_checked_matrix_algebra")

    algebra_map = _interface_algebra_homomorphism_check(
        left_ports=left_ports,
        right_ports=right_ports,
        forward=forward,
    )
    homomorphism_receipt = bool(
        bijection_valid
        and inverse_composition
        and orientation_valid
        and orientation_inverse_composition
        and schema_hash_binds_checked_algebra
        and algebra_map["receipt"]
    )

    passed = not blockers
    return {
        "schema": "oph.echosahedral_federation.seam_bundle.v1",
        "seam_id": seam.seam_id,
        "left_carrier_id": seam.left_carrier_id,
        "right_carrier_id": seam.right_carrier_id,
        "left_ports": list(left_ports),
        "right_ports": list(right_ports),
        "bundle_size": len(left_ports),
        "collar_kind": seam.collar_kind,
        "left_bundle_connected": left_connected,
        "right_bundle_connected": right_connected,
        "gluing_map_bijective": bijection_valid,
        "forward_reverse_composition_identity": inverse_composition,
        "orientation_reversing": orientation_valid,
        "forward_reverse_orientation_composition_identity": (
            orientation_inverse_composition
        ),
        "bundle_connectivity_mode": (
            "declared_antipode_involution"
            if seam.collar_kind == "antipodal_pair"
            else "induced_local_edge_incidence"
        ),
        "interface_algebra_id": binding.interface_algebra_id,
        "interface_algebra_sha256": binding.interface_algebra_sha256,
        "checked_interface_algebra_schema": checked_schema,
        "checked_interface_algebra_schema_sha256": checked_schema_sha256,
        "interface_hash_check_scope": _INTERFACE_HASH_CHECK_SCOPE,
        "interface_schema_hashes_agree": schema_hashes_agree,
        "interface_schema_hash_binds_checked_algebra": (
            schema_hash_binds_checked_algebra
        ),
        # Compatibility alias.  This means hash equality only; it is not a
        # mathematical assertion that an algebra map has been constructed.
        "endpoint_interface_algebra_hashes_agree": schema_hashes_agree,
        "interface_algebra_map": algebra_map,
        "blockers": sorted(set(blockers)),
        "SEAM_BUNDLE_RECEIPT": passed,
        "INTERFACE_SCHEMA_HASH_BINDING_RECEIPT": schema_hashes_agree,
        "INTERFACE_ALGEBRA_MAP_HOMOMORPHISM_RECEIPT": homomorphism_receipt,
        "FULL_INTERFACE_ALGEBRA_SEAM_RECEIPT": homomorphism_receipt,
        "claim_boundary": (
            "The endpoint hashes bind the canonical schema of the matrix "
            "algebra actually checked on this collar. The "
            "interface algebra map receipt verifies the permutation-induced "
            "finite matrix-algebra isomorphism exactly on the collar bundle. "
            "Triple-overlap cocycle coherence is certified at the federation "
            "level, and neither receipt implies emergent geometry."
        ),
    }


def _interface_algebra_homomorphism_check(
    *,
    left_ports: tuple[int, ...],
    right_ports: tuple[int, ...],
    forward: tuple[int, ...],
) -> dict[str, Any]:
    """Verify the induced finite interface-algebra map exactly.

    The orientation-reversing port bijection induces the conjugation map
    ``A -> P A P^T`` between the endpoint collar matrix algebras.  Every check
    below runs in exact integer arithmetic: bijectivity of the induced
    position permutation, the basis product law on elementary matrices, and
    preservation of sums, products, unit, and star on dense integer matrices.
    """

    blockers: list[str] = []
    size = len(left_ports)
    sigma: list[int] | None = None
    if (
        size == 0
        or len(right_ports) != size
        or len(forward) != size
        or len(set(left_ports)) != size
        or len(set(right_ports)) != size
    ):
        blockers.append("collar_bundles_do_not_define_a_finite_algebra_pair")
    else:
        right_position = {port: index for index, port in enumerate(right_ports)}
        positions = [right_position.get(port) for port in forward]
        if any(position is None for position in positions):
            blockers.append("gluing_image_is_not_contained_in_right_collar")
        else:
            sigma = [int(position) for position in positions]  # type: ignore[arg-type]
            if len(set(sigma)) != size:
                blockers.append("induced_position_map_is_not_injective")
                sigma = None
    basis_product_law = False
    linear_structure_preserved = False
    unit_preserved = False
    star_preserved = False
    products_preserved = False
    if sigma is not None:
        basis_product_law = all(
            (j == k) == (sigma[j] == sigma[k])
            for j in range(size)
            for k in range(size)
        )
        permutation = np.zeros((size, size), dtype=np.int64)
        for index, image in enumerate(sigma):
            permutation[image, index] = 1

        def conjugate(matrix: np.ndarray) -> np.ndarray:
            return permutation @ matrix @ permutation.T

        first = np.arange(1, size * size + 1, dtype=np.int64).reshape(size, size)
        second = first.T + 7
        linear_structure_preserved = bool(
            np.array_equal(conjugate(first + second), conjugate(first) + conjugate(second))
        )
        products_preserved = bool(
            np.array_equal(conjugate(first @ second), conjugate(first) @ conjugate(second))
        )
        unit_preserved = bool(
            np.array_equal(conjugate(np.eye(size, dtype=np.int64)), np.eye(size, dtype=np.int64))
        )
        star_preserved = bool(
            np.array_equal(conjugate(first.T), conjugate(first).T)
        )
        basis_images_exact = all(
            np.array_equal(
                conjugate(_elementary_matrix(size, i, j)),
                _elementary_matrix(size, sigma[i], sigma[j]),
            )
            for i in range(size)
            for j in range(size)
        )
        if not basis_images_exact:
            blockers.append("elementary_matrix_images_do_not_match_induced_map")
        if not basis_product_law:
            blockers.append("basis_product_law_fails_for_induced_map")
        if not (
            linear_structure_preserved
            and products_preserved
            and unit_preserved
            and star_preserved
        ):
            blockers.append("dense_algebra_operations_are_not_preserved_exactly")
    receipt = bool(sigma is not None and not blockers)
    return {
        "collar_dimension": size,
        "induced_position_permutation": sigma,
        "map_model": "conjugation_by_collar_permutation_matrix",
        "bijective_algebra_map": bool(sigma is not None),
        "basis_product_law_exact": basis_product_law,
        "sums_preserved_exact": linear_structure_preserved,
        "products_preserved_exact": products_preserved,
        "unit_preserved_exact": unit_preserved,
        "star_preserved_exact": star_preserved,
        "blockers": blockers,
        "receipt": receipt,
    }


def _elementary_matrix(size: int, row: int, column: int) -> np.ndarray:
    matrix = np.zeros((size, size), dtype=np.int64)
    matrix[row, column] = 1
    return matrix


def _seam_partial_port_maps(
    seam: SeamBundle,
) -> tuple[dict[int, int], dict[int, int]] | None:
    """Return exact forward/backward partial port maps or None if malformed."""

    left_ports, left_exact = _runtime_exact_int_tuple(seam.left_ports)
    right_ports, right_exact = _runtime_exact_int_tuple(seam.right_ports)
    forward, forward_exact = _runtime_exact_int_tuple(seam.left_to_right_ports)
    backward, backward_exact = _runtime_exact_int_tuple(seam.right_to_left_ports)
    if not (
        left_exact
        and right_exact
        and forward_exact
        and backward_exact
        and left_ports
        and len(left_ports) == len(set(left_ports)) == len(forward)
        and len(right_ports) == len(set(right_ports)) == len(backward)
        and len(left_ports) == len(right_ports)
    ):
        return None
    return (
        dict(zip(left_ports, forward, strict=True)),
        dict(zip(right_ports, backward, strict=True)),
    )


def _triple_overlap_report(
    overlap: TripleOverlapBundle,
    carrier_by_id: Mapping[str, EchosahedralCarrier],
    seam_by_id: Mapping[str, SeamBundle],
) -> dict[str, Any]:
    """Verify an explicit common restriction of three pairwise seam algebras."""

    blockers: list[str] = []
    carriers = overlap.oriented_carrier_ids
    seam_ids = overlap.oriented_seam_ids
    exact_shapes = bool(
        isinstance(carriers, tuple)
        and len(carriers) == 3
        and all(isinstance(item, str) and item for item in carriers)
        and len(set(carriers)) == 3
        and isinstance(seam_ids, tuple)
        and len(seam_ids) == 3
        and all(isinstance(item, str) and item for item in seam_ids)
        and len(set(seam_ids)) == 3
    )
    if not overlap.overlap_id:
        blockers.append("empty_triple_overlap_id")
    if not exact_shapes:
        blockers.append("triple_overlap_requires_three_distinct_ordered_carriers_and_seams")
    if any(carrier_id not in carrier_by_id for carrier_id in carriers):
        blockers.append("triple_overlap_references_unknown_carrier")
    if any(seam_id not in seam_by_id for seam_id in seam_ids):
        blockers.append("triple_overlap_references_unknown_seam")

    canonical_schema: dict[str, Any] | None = None
    canonical_hash: str | None = None
    restriction_hash_binds_checked_algebra = False
    if overlap.restriction_algebra_id:
        canonical_schema = finite_matrix_interface_algebra_schema(
            overlap.restriction_algebra_id, 1
        )
        canonical_hash = interface_algebra_sha256(canonical_schema)
        restriction_hash_binds_checked_algebra = bool(
            overlap.restriction_algebra_sha256 == canonical_hash
        )
    else:
        blockers.append("empty_triple_restriction_algebra_id")
    if not _is_sha256(overlap.restriction_algebra_sha256):
        blockers.append("invalid_triple_restriction_algebra_sha256")
    if not restriction_hash_binds_checked_algebra:
        blockers.append("triple_restriction_hash_does_not_bind_canonical_M1_algebra")

    directed_legs: list[dict[str, Any]] = []
    pairwise_restrictions_valid = exact_shapes
    if exact_shapes:
        for index, seam_id in enumerate(seam_ids):
            source = carriers[index]
            target = carriers[(index + 1) % 3]
            seam = seam_by_id.get(seam_id)
            source_port: int | None = None
            target_port: int | None = None
            leg_valid = False
            if seam is not None:
                maps = _seam_partial_port_maps(seam)
                if (
                    maps is not None
                    and seam.interface_algebra.interface_algebra_id
                    == overlap.restriction_algebra_id
                    and seam.interface_algebra.interface_algebra_sha256
                    == overlap.restriction_algebra_sha256
                    and len(seam.left_ports) == len(seam.right_ports) == 1
                ):
                    forward, backward = maps
                    if (
                        seam.left_carrier_id == source
                        and seam.right_carrier_id == target
                    ):
                        source_port = next(iter(forward))
                        target_port = forward[source_port]
                        leg_valid = True
                    elif (
                        seam.right_carrier_id == source
                        and seam.left_carrier_id == target
                    ):
                        source_port = next(iter(backward))
                        target_port = backward[source_port]
                        leg_valid = True
            if not leg_valid:
                blockers.append(
                    f"triple_restriction_leg_{index}_does_not_match_named_seam"
                )
                pairwise_restrictions_valid = False
            directed_legs.append(
                {
                    "seam_id": seam_id,
                    "source_carrier_id": source,
                    "target_carrier_id": target,
                    "source_port": source_port,
                    "target_port": target_port,
                    "restriction_to_common_M1_exact": leg_valid,
                }
            )

    local_restriction_pairs: list[dict[str, Any]] = []
    if len(directed_legs) == 3:
        for index, carrier_id in enumerate(carriers):
            incoming = directed_legs[(index - 1) % 3]["target_port"]
            outgoing = directed_legs[index]["source_port"]
            local_restriction_pairs.append(
                {
                    "carrier_id": carrier_id,
                    "incoming_seam_port": incoming,
                    "outgoing_seam_port": outgoing,
                    "common_scalar_restriction": bool(
                        incoming is not None and outgoing is not None
                    ),
                }
            )
    local_common_restrictions_valid = bool(
        len(local_restriction_pairs) == 3
        and all(row["common_scalar_restriction"] for row in local_restriction_pairs)
    )
    if exact_shapes and not local_common_restrictions_valid:
        blockers.append("local_pairwise_collars_lack_common_triple_restrictions")

    # M_1(C) has one matrix unit.  Every unital complex star-homomorphism
    # fixes it, so six alternating seam/local restriction maps send E_00
    # exactly back to E_00.  This is an algebra cocycle, not a comparison of
    # unrelated local integer port labels.
    identity_composition = bool(
        exact_shapes
        and pairwise_restrictions_valid
        and local_common_restrictions_valid
        and restriction_hash_binds_checked_algebra
    )
    if exact_shapes and not identity_composition:
        blockers.append("triple_overlap_algebra_composition_is_not_identity")
    receipt = not blockers
    return {
        "schema": "oph.echosahedral_federation.triple_overlap.v1",
        "overlap_id": overlap.overlap_id,
        "oriented_carrier_ids": list(carriers),
        "oriented_seam_ids": list(seam_ids),
        "restriction_algebra_id": overlap.restriction_algebra_id,
        "restriction_algebra_sha256": overlap.restriction_algebra_sha256,
        "checked_restriction_algebra_schema": canonical_schema,
        "checked_restriction_algebra_schema_sha256": canonical_hash,
        "restriction_hash_binds_checked_algebra": (
            restriction_hash_binds_checked_algebra
        ),
        "directed_seam_restrictions": directed_legs,
        "local_common_restriction_pairs": local_restriction_pairs,
        "matrix_unit_count": 1,
        "matrix_unit_identity_composition_exact": identity_composition,
        "blockers": sorted(set(blockers)),
        "TRIPLE_OVERLAP_RESTRICTION_RECEIPT": receipt,
        "EXACT_IDENTITY_COCYCLE_ON_TRIPLE_RESTRICTION": identity_composition,
        "claim_boundary": (
            "This verifies a declared nonempty common M_1(C) restriction of "
            "three typed pairwise collars. It does not infer a triple overlap "
            "from carrier-graph adjacency alone."
        ),
    }


def _higher_overlap_cocycle_report(
    seams: Sequence[SeamBundle],
    carrier_ids: Iterable[str],
    *,
    explicit_triple_rows: Sequence[Mapping[str, Any]] = (),
) -> dict[str, Any]:
    """Check every composable seam triangle for exact identity holonomy.

    A triangle is a triple of federated carriers with a seam on each pair.
    Around each triangle every combination of the three partial port
    bijections is composed in cyclic order; wherever the composition is
    defined it must return the starting port.  Empty composition domains are
    counted separately and impose no condition, mirroring the Cech cocycle
    convention on empty triple overlaps.
    """

    known = set(carrier_ids)
    neighbor_map: dict[str, set[str]] = {}
    pair_maps: dict[tuple[str, str], list[tuple[dict[int, int], dict[int, int], str]]] = {}
    malformed_seam_count = 0
    for seam in seams:
        if (
            seam.left_carrier_id not in known
            or seam.right_carrier_id not in known
            or seam.left_carrier_id == seam.right_carrier_id
        ):
            malformed_seam_count += 1
            continue
        maps = _seam_partial_port_maps(seam)
        if maps is None:
            malformed_seam_count += 1
            continue
        forward_map, backward_map = maps
        left_id = seam.left_carrier_id
        right_id = seam.right_carrier_id
        neighbor_map.setdefault(left_id, set()).add(right_id)
        neighbor_map.setdefault(right_id, set()).add(left_id)
        pair_maps.setdefault((left_id, right_id), []).append(
            (forward_map, backward_map, seam.seam_id)
        )
    triangles: set[tuple[str, str, str]] = set()
    for first_id, second_id in pair_maps:
        common = neighbor_map.get(first_id, set()) & neighbor_map.get(second_id, set())
        for third_id in common:
            triangles.add(tuple(sorted((first_id, second_id, third_id))))

    def _leg_maps(source: str, target: str) -> list[tuple[dict[int, int], str]]:
        rows: list[tuple[dict[int, int], str]] = []
        for forward_map, _, seam_id in pair_maps.get((source, target), []):
            rows.append((forward_map, seam_id))
        for _, backward_map, seam_id in pair_maps.get((target, source), []):
            rows.append((backward_map, seam_id))
        return rows

    composable_loop_count = 0
    violations: list[dict[str, Any]] = []
    for first_id, second_id, third_id in sorted(triangles):
        for first_map, first_seam in _leg_maps(first_id, second_id):
            for second_map, second_seam in _leg_maps(second_id, third_id):
                for third_map, third_seam in _leg_maps(third_id, first_id):
                    for port, mid in first_map.items():
                        if mid not in second_map:
                            continue
                        far = second_map[mid]
                        if far not in third_map:
                            continue
                        composable_loop_count += 1
                        if third_map[far] != port:
                            if len(violations) < _REPORT_DETAIL_LIMIT:
                                violations.append(
                                    {
                                        "triangle": [first_id, second_id, third_id],
                                        "seam_ids": [
                                            first_seam,
                                            second_seam,
                                            third_seam,
                                        ],
                                        "start_port": port,
                                        "returned_port": third_map[far],
                                    }
                                )
    inferred_cocycle_condition = bool(
        pair_maps
        and malformed_seam_count == 0
        and not violations
    )
    explicit_triple_count = len(explicit_triple_rows)
    valid_explicit_triple_count = sum(
        row.get("TRIPLE_OVERLAP_RESTRICTION_RECEIPT") is True
        and row.get("EXACT_IDENTITY_COCYCLE_ON_TRIPLE_RESTRICTION") is True
        for row in explicit_triple_rows
    )
    explicit_triples_valid = bool(
        valid_explicit_triple_count == explicit_triple_count
    )
    cocycle_condition = bool(
        inferred_cocycle_condition and explicit_triples_valid
    )
    positive_triangle_witness = bool(
        cocycle_condition
        and explicit_triple_count > 0
        and valid_explicit_triple_count == explicit_triple_count
    )
    return {
        "schema": "oph.echosahedral_federation.higher_overlap_cocycle.v1",
        "triangle_count": len(triangles),
        "composable_port_loop_count": composable_loop_count,
        "composable_triple_restriction_algebra_loop_count": (
            valid_explicit_triple_count
        ),
        "declared_nonempty_triple_overlap_count": explicit_triple_count,
        "valid_declared_triple_overlap_count": valid_explicit_triple_count,
        "inferred_port_label_cocycle_condition": inferred_cocycle_condition,
        "positive_triangle_witness": positive_triangle_witness,
        "vacuous_triangle_convention": (
            "empty_triple_overlap_domains_impose_no_condition"
        ),
        "declared_deck_elements": "identity_only",
        "malformed_seam_count": malformed_seam_count,
        "identity_violation_count": len(violations),
        "identity_violations": violations,
        "HIGHER_OVERLAP_COCYCLE_CONDITION_RECEIPT": cocycle_condition,
        "NONVACUOUS_HIGHER_OVERLAP_COCYCLE_WITNESS": positive_triangle_witness,
        # Compatibility alias: this is the logical Cech condition, which is
        # vacuously true when the declared cover has no triple overlap.
        "HIGHER_OVERLAP_COCYCLE_RECEIPT": cocycle_condition,
    }


def _external_boundary_report(
    boundary: ExternalBoundaryBundle,
    carrier_by_id: Mapping[str, EchosahedralCarrier],
) -> dict[str, Any]:
    blockers: list[str] = []
    carrier = carrier_by_id.get(boundary.carrier_id)
    ports, ports_exact = _runtime_exact_int_tuple(boundary.ports)
    if not boundary.boundary_id:
        blockers.append("empty_boundary_id")
    if carrier is None:
        blockers.append("unknown_boundary_carrier")
    if not ports_exact:
        blockers.append("external_boundary_ports_must_be_an_exact_integer_tuple")
    if not ports or len(ports) != len(set(ports)):
        blockers.append("external_boundary_ports_must_be_nonempty_and_unique")
    if any(port < 0 or port >= _PORT_COUNT for port in ports):
        blockers.append("external_boundary_port_out_of_range")
    connected = bool(carrier is not None and _port_subset_connected(carrier, ports))
    if not connected:
        blockers.append("external_boundary_port_bundle_is_not_connected")
    if boundary.boundary_condition not in {
        "open_external",
        "fixed_external",
        "measured_external",
    }:
        blockers.append("unknown_external_boundary_condition")
    if not _is_sha256(boundary.boundary_algebra_sha256):
        blockers.append("invalid_boundary_algebra_sha256")
    return {
        "schema": "oph.echosahedral_federation.external_boundary.v1",
        "boundary_id": boundary.boundary_id,
        "carrier_id": boundary.carrier_id,
        "ports": list(ports),
        "port_count": len(ports),
        "connected": connected,
        "boundary_condition": boundary.boundary_condition,
        "boundary_algebra_sha256": boundary.boundary_algebra_sha256,
        "blockers": sorted(set(blockers)),
        "EXPLICIT_EXTERNAL_BOUNDARY_RECEIPT": not blockers,
    }


def presentation_firewall_report(
    federation: EchosahedralFederation,
    promoted_payload: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Reject hidden carrier presentation data from promoted payloads."""

    payload = {} if promoted_payload is None else promoted_payload
    encountered_keys = _recursive_mapping_keys(payload)
    forbidden_keys = sorted(encountered_keys & CARRIER_FORBIDDEN_PROMOTION_FIELDS)
    coordinate_leaks = _find_hidden_coordinate_values(payload, federation.carriers)
    port_label_leaks = _find_hidden_port_label_values(payload, federation.carriers)
    passed = not forbidden_keys and not coordinate_leaks and not port_label_leaks
    return {
        "schema": "oph.echosahedral_federation.presentation_firewall.v1",
        "hidden_presentation_fields": sorted(HIDDEN_PRESENTATION_FIELDS),
        "quotient_visible_carrier_fields": sorted(QUOTIENT_VISIBLE_CARRIER_FIELDS),
        "forbidden_promoted_fields": sorted(CARRIER_FORBIDDEN_PROMOTION_FIELDS),
        "encountered_payload_keys": sorted(encountered_keys),
        "forbidden_keys_present": forbidden_keys,
        "hidden_coordinate_value_paths": coordinate_leaks,
        "hidden_port_label_value_paths": port_label_leaks,
        "CARRIER_PRESENTATION_FIREWALL_RECEIPT": passed,
        "HIDDEN_COORDINATES_EXCLUDED_FROM_PROMOTED_PAYLOAD_RECEIPT": bool(
            not coordinate_leaks
            and not (
                {"port_coordinates", "hidden_port_coordinates"} & set(forbidden_keys)
            )
        ),
        "LOCAL_PORT_NAMES_EXCLUDED_FROM_PROMOTED_PAYLOAD_RECEIPT": bool(
            not port_label_leaks
            and not ({"port_names", "local_port_names"} & set(forbidden_keys))
        ),
        "claim_boundary": (
            "The schema firewall excludes known presentation fields and exact hidden "
            "coordinate/name values. It does not by itself prove quotient invariance "
            "of a future dynamical emergence evaluator."
        ),
    }


def relabel_federation_ports(
    federation: EchosahedralFederation,
    port_permutations: Mapping[str, Sequence[int]],
) -> EchosahedralFederation:
    """Co-transform every carrier, seam, and external-boundary port label."""

    carrier_ids = {carrier.carrier_id for carrier in federation.carriers}
    if set(port_permutations) != carrier_ids:
        raise ValueError("one explicit port permutation is required per carrier")
    normalized = {
        carrier_id: _validated_port_permutation(permutation)
        for carrier_id, permutation in port_permutations.items()
    }
    carriers = tuple(
        _relabel_carrier(carrier, normalized[carrier.carrier_id])
        for carrier in federation.carriers
    )
    seams = []
    for seam in federation.seams:
        left_permutation = normalized[seam.left_carrier_id]
        right_permutation = normalized[seam.right_carrier_id]
        seams.append(
            replace(
                seam,
                left_ports=tuple(left_permutation[port] for port in seam.left_ports),
                right_ports=tuple(right_permutation[port] for port in seam.right_ports),
                left_to_right_ports=tuple(
                    right_permutation[port] for port in seam.left_to_right_ports
                ),
                right_to_left_ports=tuple(
                    left_permutation[port] for port in seam.right_to_left_ports
                ),
            )
        )
    boundaries = tuple(
        replace(
            boundary,
            ports=tuple(
                normalized[boundary.carrier_id][port] for port in boundary.ports
            ),
        )
        for boundary in federation.external_boundaries
    )
    return replace(
        federation,
        carriers=carriers,
        seams=tuple(seams),
        external_boundaries=boundaries,
    )


def carrier_quotient_invariance_report(
    source: EchosahedralFederation,
    transformed: EchosahedralFederation,
    port_permutations: Mapping[str, Sequence[int]],
) -> dict[str, Any]:
    """Verify one exact presentation relabeling and its quotient-visible export."""

    blockers: list[str] = []
    try:
        expected = relabel_federation_ports(source, port_permutations)
    except (KeyError, TypeError, ValueError) as exc:
        expected = None
        blockers.append(f"invalid_presentation_transform:{type(exc).__name__}")
    exact_cotransformation = bool(expected is not None and expected == transformed)
    if not exact_cotransformation:
        blockers.append("transformed_federation_is_not_exact_co_transformation")
    source_sewing = federation_sewing_report(source)
    transformed_sewing = federation_sewing_report(transformed)
    if not source_sewing["FEDERATION_SEWING_RECEIPT"]:
        blockers.append("source_federation_sewing_invalid")
    if not transformed_sewing["FEDERATION_SEWING_RECEIPT"]:
        blockers.append("transformed_federation_sewing_invalid")
    source_payload = _quotient_visible_contract_payload(source)
    transformed_payload = _quotient_visible_contract_payload(transformed)
    source_hash = interface_algebra_sha256(source_payload)
    transformed_hash = interface_algebra_sha256(transformed_payload)
    quotient_exports_equal = source_hash == transformed_hash
    if not quotient_exports_equal:
        blockers.append("quotient_visible_contract_export_changed_under_relabeling")
    source_firewall = presentation_firewall_report(source)
    transformed_firewall = presentation_firewall_report(transformed)
    passed = bool(
        not blockers
        and exact_cotransformation
        and quotient_exports_equal
        and source_firewall["CARRIER_PRESENTATION_FIREWALL_RECEIPT"]
        and transformed_firewall["CARRIER_PRESENTATION_FIREWALL_RECEIPT"]
    )
    return {
        "schema": "oph.echosahedral_federation.quotient_invariance.v1",
        "source_federation_id": source.federation_id,
        "transformed_federation_id": transformed.federation_id,
        "transformation_class": "explicit_local_port_relabeling_with_seams_and_boundaries_cotransformed",
        "carrier_permutation_count": len(port_permutations),
        "exact_cotransformation_verified": exact_cotransformation,
        "source_quotient_visible_contract_sha256": source_hash,
        "transformed_quotient_visible_contract_sha256": transformed_hash,
        "quotient_visible_contract_exports_equal": quotient_exports_equal,
        "blockers": sorted(set(blockers)),
        "CARRIER_QUOTIENT_INVARIANCE_RECEIPT": passed,
        "FULL_DYNAMICAL_QUOTIENT_INVARIANCE_RECEIPT": False,
        "claim_boundary": (
            "This receipt proves invariance of the finite carrier/sewing contract under "
            "the supplied exact relabeling. Engine dynamics, repair schedules, and "
            "downstream geometry evaluators are not yet integrated or certified."
        ),
    }


_REFINEMENT_DEFECT_MATCH_TOLERANCE = 1.0e-9
_REFINEMENT_IDENTITY_TOLERANCE = 5.0e-14
_REFINEMENT_NATURALITY_TOLERANCE = 5.0e-12
_REFINEMENT_MAX_LEVEL = 2


def carrier_refinement_naturality_report(
    *,
    max_level: int = _REFINEMENT_MAX_LEVEL,
    embedding_override: Sequence[int] | None = None,
    port_permutation_override: Sequence[Sequence[int]] | None = None,
) -> dict[str, Any]:
    """Recompute the carrier-to-support refinement naturality contract.

    The default path recomputes every map from the shared reference template
    and the exact geodesic tower.  The override arguments exist for adversarial
    verification: a supplied embedding or action family is validated against
    the recomputed reference data and any mismatch fails closed with named
    blockers.
    """

    if embedding_override is None and port_permutation_override is None:
        return copy.deepcopy(_carrier_refinement_naturality_cached(int(max_level)))
    return _compute_carrier_refinement_naturality(
        max_level=int(max_level),
        embedding_override=embedding_override,
        port_permutation_override=port_permutation_override,
    )


def controlled_oriented_s2_limit_report(
    federation: EchosahedralFederation,
    *,
    max_level: int = 3,
) -> dict[str, Any]:
    """Verify that one source-incidence nerve has controlled oriented-S2 limit.

    The source-side assertion is a simplicial one: carriers, seams, and
    declared nonempty triple restrictions must map bijectively to the
    vertices, edges, and outward faces of the certified icosahedral boundary.
    The support-side assertion uses normalized edge-midpoint subdivision.

    If ``delta_n`` is the largest geodesic edge at level ``n``, a child edge
    joining an old endpoint to a midpoint has length at most ``delta_n/2``.
    For an edge joining two normalized midpoints, the normalization inequality

        ||u/||u|| - v/||v|||| <= 2 ||u-v|| / (||u|| + ||v||)

    gives ``2 asin(tan(delta_n/2)/2)``.  On
    ``[0, delta_0]``, with ``delta_0 = acos(1/sqrt(5))``, elementary
    monotonicity/concavity bounds give this quantity at most
    ``(3/4) delta_n``.  Hence ``delta_n <= delta_0 (3/4)^n -> 0``.  Every
    level is an outward-oriented Euler-two spherical triangulation, so the
    controlled limit is the oriented unit two-sphere.
    """

    blockers: list[str] = []
    binding = federation.support_tower_binding
    binding_present = binding is not None
    if binding is None:
        blockers.append("source_support_tower_binding_missing")
        binding_rows: tuple[tuple[str, int], ...] = ()
    else:
        binding_rows = binding.carrier_to_base_vertex
        if binding.geometry_family != "nested_geodesic_icosahedral":
            blockers.append("support_geometry_family_mismatch")
        if binding.orientation != "source_oriented_faces_outward":
            blockers.append("support_orientation_contract_mismatch")
        if binding.refinement_rule != "normalized_edge_midpoint_four_child":
            blockers.append("support_refinement_rule_mismatch")
        if binding.source_incidence_sha256 != _reference_incidence_sha256():
            blockers.append("support_source_incidence_hash_mismatch")

    map_exact = bool(
        isinstance(binding_rows, tuple)
        and len(binding_rows) == _PORT_COUNT
        and all(
            isinstance(row, tuple)
            and len(row) == 2
            and isinstance(row[0], str)
            and row[0]
            and type(row[1]) is int
            for row in binding_rows
        )
    )
    carrier_to_vertex = dict(binding_rows) if map_exact else {}
    carrier_ids = {carrier.carrier_id for carrier in federation.carriers}
    map_bijective = bool(
        map_exact
        and len(carrier_to_vertex) == _PORT_COUNT
        and set(carrier_to_vertex) == carrier_ids
        and set(carrier_to_vertex.values()) == set(range(_PORT_COUNT))
    )
    if not map_bijective:
        blockers.append("carrier_to_base_vertex_map_is_not_a_bijection")

    template = _reference_echosahedral_carrier_template()
    expected_edges = {tuple(sorted(edge)) for edge in template.edges}
    observed_edges: set[tuple[int, int]] = set()
    seam_source_law = map_bijective
    if map_bijective:
        for seam in federation.seams:
            if (
                seam.left_carrier_id not in carrier_to_vertex
                or seam.right_carrier_id not in carrier_to_vertex
            ):
                seam_source_law = False
                continue
            left_vertex = carrier_to_vertex[seam.left_carrier_id]
            right_vertex = carrier_to_vertex[seam.right_carrier_id]
            observed_edges.add(tuple(sorted((left_vertex, right_vertex))))
            if not (
                seam.left_ports == (right_vertex,)
                and seam.right_ports == (left_vertex,)
                and seam.left_to_right_ports == (left_vertex,)
                and seam.right_to_left_ports == (right_vertex,)
            ):
                seam_source_law = False
    edge_nerve_exact = bool(
        seam_source_law
        and len(federation.seams) == len(expected_edges) == 30
        and observed_edges == expected_edges
    )
    if not edge_nerve_exact:
        blockers.append("federation_seams_are_not_the_source_incidence_edges")

    expected_faces = {_cyclic_face_key(face) for face in template.faces}
    observed_faces: set[tuple[int, int, int]] = set()
    if map_bijective:
        for overlap in federation.triple_overlaps:
            try:
                face = tuple(
                    carrier_to_vertex[carrier_id]
                    for carrier_id in overlap.oriented_carrier_ids
                )
            except KeyError:
                continue
            if len(face) == 3:
                observed_faces.add(_cyclic_face_key(face))
    face_nerve_exact = bool(
        len(federation.triple_overlaps) == len(expected_faces) == 20
        and observed_faces == expected_faces
    )
    if not face_nerve_exact:
        blockers.append(
            "declared_triple_overlaps_are_not_the_source_oriented_faces"
        )

    level_cap = max(2, int(max_level))
    tower = build_geodesic_icosahedral_tower(level_cap)
    tower_receipt = tower.receipt()
    level_rows: list[dict[str, Any]] = []
    base_delta = math.acos(1.0 / math.sqrt(5.0))
    contraction_factor = 0.75
    prefix_contraction = True
    oriented_sphere_prefix = True
    for mesh in tower.levels:
        endpoints = mesh.vertices[mesh.edges]
        edge_dots = np.sum(endpoints[:, 0, :] * endpoints[:, 1, :], axis=1)
        edge_angles = np.arccos(np.clip(edge_dots, -1.0, 1.0))
        maximum_angle = float(np.max(edge_angles))
        bound = base_delta * contraction_factor**mesh.level
        within_bound = bool(maximum_angle <= bound + 5.0e-14)
        receipt = mesh.receipt()
        oriented = bool(
            receipt["GEODESIC_ICOSAHEDRAL_GEOMETRY_RECEIPT"]
            and receipt["outward_oriented"]
            and receipt["euler_characteristic"] == 2
        )
        prefix_contraction = prefix_contraction and within_bound
        oriented_sphere_prefix = oriented_sphere_prefix and oriented
        level_rows.append(
            {
                "level": mesh.level,
                "vertex_count": mesh.vertex_count,
                "edge_count": mesh.edge_count,
                "face_count": mesh.face_count,
                "maximum_geodesic_edge_angle": maximum_angle,
                "analytic_mesh_bound": bound,
                "mesh_bound_verified": within_bound,
                "outward_oriented_euler_two_sphere": oriented,
                "geometry_hash": mesh.geometry_hash,
            }
        )

    half_base = base_delta / 2.0
    tangent_bound_coefficient = 1.0 / (2.0 * math.cos(half_base))
    sine_concavity_coefficient = math.sin(3.0 * half_base / 4.0) / half_base
    analytic_endpoint_inequality = bool(
        0.0 < base_delta < math.pi
        and tangent_bound_coefficient <= sine_concavity_coefficient
        and contraction_factor < 1.0
    )
    analytic_all_levels_contraction = analytic_endpoint_inequality
    if not analytic_all_levels_contraction:
        blockers.append("analytic_mesh_contraction_inequality_failed")
    if not prefix_contraction:
        blockers.append("computed_refinement_prefix_exceeds_analytic_mesh_bound")
    if not oriented_sphere_prefix:
        blockers.append("refinement_prefix_is_not_outward_oriented_sphere")

    receipt = bool(
        binding_present
        and map_bijective
        and edge_nerve_exact
        and face_nerve_exact
        and tower_receipt["GEODESIC_ICOSAHEDRAL_TOWER_RECEIPT"]
        and prefix_contraction
        and oriented_sphere_prefix
        and analytic_all_levels_contraction
        and not blockers
    )
    return {
        "schema": "oph.echosahedral_federation.controlled_oriented_s2_limit.v1",
        "source_incidence_sha256": _reference_incidence_sha256(),
        "support_tower_binding_present": binding_present,
        "carrier_to_base_vertex_bijection": map_bijective,
        "source_edge_nerve_exact": edge_nerve_exact,
        "source_oriented_face_nerve_exact": face_nerve_exact,
        "nerve_simplex_counts": {
            "vertices": len(carrier_to_vertex),
            "edges": len(observed_edges),
            "oriented_faces": len(observed_faces),
        },
        "topological_identification": (
            "simplicial_isomorphism_to_oriented_icosahedral_boundary"
        ),
        "refinement_family": "normalized_edge_midpoint_four_child",
        "base_maximum_geodesic_edge_angle": base_delta,
        "analytic_contraction_factor": contraction_factor,
        "analytic_mesh_bound": "delta_n <= delta_0 * (3/4)^n",
        "normalization_inequality": (
            "norm(normalize(u)-normalize(v)) <= "
            "2*norm(u-v)/(norm(u)+norm(v))"
        ),
        "analytic_endpoint_coefficients": {
            "tangent_upper_coefficient": tangent_bound_coefficient,
            "sine_concavity_lower_coefficient": sine_concavity_coefficient,
        },
        "analytic_all_levels_contraction_proved": (
            analytic_all_levels_contraction
        ),
        "mesh_bound_tends_to_zero": analytic_all_levels_contraction,
        "computed_prefix_level_count": len(tower.levels),
        "computed_prefix_mesh_bounds_verified": prefix_contraction,
        "computed_prefix_outward_oriented": oriented_sphere_prefix,
        "levels": level_rows,
        "blockers": sorted(set(blockers)),
        "INCIDENCE_NERVE_TO_SUPPORT_SIMPLICIAL_ISOMORPHISM_RECEIPT": bool(
            map_bijective and edge_nerve_exact and face_nerve_exact
        ),
        "CONTROLLED_ORIENTED_S2_LIMIT_RECEIPT": receipt,
        "claim_boundary": (
            "This is a source-incidence theorem for the commutative oriented "
            "support: the explicit federation nerve is the icosahedral "
            "boundary and its normalized-midpoint mesh size tends to zero. "
            "It supplies no length scale, H3 frame, event manifold, BW/KMS "
            "clock, continuum QFT, or laboratory attachment."
        ),
    }


@lru_cache(maxsize=4)
def _carrier_refinement_naturality_cached(max_level: int) -> dict[str, Any]:
    return _compute_carrier_refinement_naturality(
        max_level=max_level,
        embedding_override=None,
        port_permutation_override=None,
    )


def _nearest_vertex_matching(
    points: np.ndarray, vertices: np.ndarray
) -> tuple[list[int], float]:
    """Match each point to its nearest vertex; return indices and max distance."""

    squared = np.sum(
        (points[:, None, :] - vertices[None, :, :]) ** 2,
        axis=-1,
    )
    nearest = np.argmin(squared, axis=1)
    residual = float(np.sqrt(np.max(np.min(squared, axis=1)))) if points.size else 0.0
    return [int(value) for value in nearest], residual


def _compute_carrier_refinement_naturality(
    *,
    max_level: int,
    embedding_override: Sequence[int] | None,
    port_permutation_override: Sequence[Sequence[int]] | None,
) -> dict[str, Any]:
    blockers: list[str] = []
    template = _reference_echosahedral_carrier_template()
    coordinates = np.asarray(template.port_coordinates, dtype=float)
    twelve_port_type = bool(
        len(template.port_names) == _PORT_COUNT
        and coordinates.shape == (_PORT_COUNT, 3)
    )
    level_cap = max(2, int(max_level))
    tower = build_geodesic_icosahedral_tower(level_cap)
    base = tower.levels[0]

    # (a) Carrier embedding: ports match the degree-five defect vertices of
    # tower level zero by nearest-vertex bijection.
    degrees = np.bincount(
        np.asarray(base.edges, dtype=np.int64).reshape(-1),
        minlength=base.vertex_count,
    )
    defect_ids = np.flatnonzero(degrees == 5).astype(np.int64)
    embedding: list[int] = []
    embedding_residual = math.inf
    if defect_ids.size != _PORT_COUNT:
        blockers.append("level_zero_defect_vertex_count_is_not_twelve")
    else:
        defect_coordinates = base.vertices[defect_ids]
        matched, embedding_residual = _nearest_vertex_matching(
            coordinates, defect_coordinates
        )
        embedding = [int(defect_ids[index]) for index in matched]
    embedding_bijective = bool(
        len(embedding) == _PORT_COUNT and len(set(embedding)) == _PORT_COUNT
    )
    if not embedding_bijective:
        blockers.append("carrier_port_to_defect_vertex_matching_is_not_a_bijection")
    if embedding_residual > _REFINEMENT_DEFECT_MATCH_TOLERANCE:
        blockers.append("carrier_port_to_defect_vertex_residual_exceeds_tolerance")
    declared_embedding = embedding
    if embedding_override is not None:
        declared_embedding = [int(value) for value in embedding_override]
        if declared_embedding != embedding:
            blockers.append(
                "declared_embedding_differs_from_recomputed_nearest_vertex_matching"
            )
    embedding_supplied = bool(
        embedding_bijective
        and embedding_residual <= _REFINEMENT_DEFECT_MATCH_TOLERANCE
        and declared_embedding == embedding
    )

    # (b) Coarse-graining maps: adjacent-level conditional expectations with
    # exact normalization and state preservation receipts.
    refinements = tower.cell_refinements
    refinement_rows: list[dict[str, Any]] = []
    refinement_maps_valid = len(refinements) >= 2
    if len(refinements) < 2:
        blockers.append("refinement_tower_depth_below_two_levels")
    for mapping in refinements:
        row_valid = bool(
            mapping.normalization_residual <= _REFINEMENT_IDENTITY_TOLERANCE
            and mapping.state_preservation_residual <= _REFINEMENT_IDENTITY_TOLERANCE
        )
        refinement_maps_valid = refinement_maps_valid and row_valid
        if not row_valid:
            blockers.append(
                "cell_refinement_receipt_failed_for_levels_"
                f"{mapping.coarse_level}_{mapping.fine_level}"
            )
        refinement_rows.append(
            {
                "coarse_level": mapping.coarse_level,
                "fine_level": mapping.fine_level,
                "normalization_residual": mapping.normalization_residual,
                "state_preservation_residual": mapping.state_preservation_residual,
                "map_hash": mapping.map_hash,
            }
        )

    # Port-supported realization: each port observable spreads with equal
    # weight over the level-zero cells incident to its embedded defect vertex.
    identity_residual = math.inf
    realization: np.ndarray | None = None
    if embedding_supplied:
        vertex_to_port = {vertex: port for port, vertex in enumerate(embedding)}
        realization = np.zeros((base.face_count, _PORT_COUNT), dtype=float)
        for face_index, face in enumerate(base.faces):
            for vertex in face:
                realization[face_index, vertex_to_port[int(vertex)]] += 1.0 / 3.0
        if refinement_maps_valid:
            round_trip = tower.conditional_expectation_cells(
                tower.embed_cells(realization, coarse_level=0, fine_level=2),
                fine_level=2,
                coarse_level=0,
            )
            adjacent = refinements[0].conditional_expectation(
                refinements[0].embed(realization)
            )
            identity_residual = float(
                max(
                    np.max(np.abs(round_trip - realization)),
                    np.max(np.abs(adjacent - realization)),
                )
            )
    identity_verified = bool(
        realization is not None
        and refinement_maps_valid
        and identity_residual <= _REFINEMENT_IDENTITY_TOLERANCE
    )
    if not identity_verified:
        blockers.append("coarse_after_embedding_identity_residual_exceeds_tolerance")

    # (c) Seam-law naturality: every registered A5 port action lifts to face
    # permutations at each level that intertwine lineage, expectation weights,
    # and the port-supported realization.
    if port_permutation_override is None:
        actions = tuple(
            tuple(int(value) for value in row) for row in template.a5_actions
        )
    else:
        actions = tuple(
            tuple(int(value) for value in row) for row in port_permutation_override
        )
    family_valid = bool(
        len(actions) == 60
        and len(set(actions)) == 60
        and all(sorted(row) == list(range(_PORT_COUNT)) for row in actions)
    )
    if not family_valid:
        blockers.append("a5_action_family_is_not_60_unique_port_bijections")
    face_lookups = []
    for mesh in tower.levels:
        face_lookups.append(
            {
                tuple(sorted(int(value) for value in face)): face_index
                for face_index, face in enumerate(mesh.faces)
            }
        )
    action_violation_examples: list[dict[str, Any]] = []
    seam_law_violation_count = 0
    realization_equivariance_residual = 0.0
    chart_equivariance_residuals = {1: 0.0, 2: 0.0}
    chart_equivariant = bool(family_valid and realization is not None)
    embedded_realizations: dict[int, np.ndarray] = {}
    if realization is not None:
        for fine in (1, 2):
            embedded_realizations[fine] = tower.embed_cells(
                realization, coarse_level=0, fine_level=fine
            )
    if family_valid:
        for action_index, row in enumerate(actions):
            row_array = np.asarray(row, dtype=np.int64)
            failure: str | None = None
            rotation, fit_residual = _proper_rotation_fit(
                coordinates, coordinates[row_array]
            )
            if (
                fit_residual > _REFINEMENT_NATURALITY_TOLERANCE
                or float(np.linalg.det(rotation)) <= 0.0
            ):
                failure = "port_permutation_has_no_proper_rotation_realization"
            face_perms: list[np.ndarray] = []
            if failure is None:
                for level_index, mesh in enumerate(tower.levels):
                    mapped_vertices, vertex_residual = _nearest_vertex_matching(
                        np.asarray(mesh.vertices, dtype=float) @ rotation,
                        np.asarray(mesh.vertices, dtype=float),
                    )
                    if (
                        vertex_residual > _REFINEMENT_NATURALITY_TOLERANCE
                        or len(set(mapped_vertices)) != mesh.vertex_count
                    ):
                        failure = (
                            "rotation_does_not_permute_level_"
                            f"{level_index}_vertices"
                        )
                        break
                    lookup = face_lookups[level_index]
                    mapped_faces = np.empty(mesh.face_count, dtype=np.int64)
                    for face_index, face in enumerate(mesh.faces):
                        key = tuple(
                            sorted(mapped_vertices[int(vertex)] for vertex in face)
                        )
                        target = lookup.get(key)
                        if target is None:
                            failure = (
                                "rotation_does_not_permute_level_"
                                f"{level_index}_faces"
                            )
                            break
                        mapped_faces[face_index] = target
                    if failure is not None:
                        break
                    face_perms.append(mapped_faces)
            if failure is None:
                for mapping in refinements:
                    coarse_perm = face_perms[mapping.coarse_level]
                    fine_perm = face_perms[mapping.fine_level]
                    lineage_preserved = bool(
                        np.array_equal(
                            mapping.child_to_parent[fine_perm],
                            coarse_perm[mapping.child_to_parent],
                        )
                    )
                    weight_residual = float(
                        np.max(
                            np.abs(
                                mapping.conditional_expectation_weights[fine_perm]
                                - mapping.conditional_expectation_weights
                            )
                        )
                    )
                    if (
                        not lineage_preserved
                        or weight_residual > _REFINEMENT_NATURALITY_TOLERANCE
                    ):
                        failure = (
                            "refinement_lineage_not_equivariant_for_levels_"
                            f"{mapping.coarse_level}_{mapping.fine_level}"
                        )
                        break
            if failure is None and realization is not None:
                inverse_base_faces = _inverse_permutation(
                    tuple(int(value) for value in face_perms[0])
                )
                base_residual = float(
                    np.max(
                        np.abs(
                            realization[:, row_array]
                            - realization[np.asarray(inverse_base_faces), :]
                        )
                    )
                )
                realization_equivariance_residual = max(
                    realization_equivariance_residual, base_residual
                )
                if base_residual > _REFINEMENT_NATURALITY_TOLERANCE:
                    failure = "port_realization_a5_equivariance_failed_at_level_0"
                else:
                    for fine in (1, 2):
                        inverse_fine = _inverse_permutation(
                            tuple(int(value) for value in face_perms[fine])
                        )
                        fine_residual = float(
                            np.max(
                                np.abs(
                                    embedded_realizations[fine][:, row_array]
                                    - embedded_realizations[fine][
                                        np.asarray(inverse_fine), :
                                    ]
                                )
                            )
                        )
                        chart_equivariance_residuals[fine] = max(
                            chart_equivariance_residuals[fine], fine_residual
                        )
                        if fine_residual > _REFINEMENT_NATURALITY_TOLERANCE:
                            failure = (
                                "chart_realization_a5_equivariance_failed_at_level_"
                                f"{fine}"
                            )
                            chart_equivariant = False
                            break
            if failure is not None:
                seam_law_violation_count += 1
                if len(action_violation_examples) < 8:
                    action_violation_examples.append(
                        {"action_index": action_index, "failure": failure}
                    )
    seam_law_verified = bool(
        family_valid and embedding_supplied and seam_law_violation_count == 0
    )
    if family_valid and seam_law_violation_count:
        blockers.append(
            f"seam_law_naturality_violated_for_{seam_law_violation_count}_of_60_actions"
        )

    # (d) Quotient commutation: presentation relabelings co-transform the
    # embedding and the port-supported realization, and the typed quotient
    # invariance machinery accepts the same relabelings.
    quotient_checks: list[dict[str, Any]] = []
    quotient_verified = embedding_supplied and realization is not None
    if embedding_supplied and realization is not None:
        relabelings = (
            tuple(reversed(range(_PORT_COUNT))),
            tuple(int(value) for value in template.a5_actions[7]),
        )
        for relabeling in relabelings:
            relabeled = _relabel_carrier(template, relabeling)
            relabeled_coordinates = np.asarray(
                relabeled.port_coordinates, dtype=float
            )
            matched, relabeled_residual = _nearest_vertex_matching(
                relabeled_coordinates, base.vertices[defect_ids]
            )
            relabeled_embedding = [int(defect_ids[index]) for index in matched]
            embedding_commutes = bool(
                relabeled_residual <= _REFINEMENT_DEFECT_MATCH_TOLERANCE
                and all(
                    relabeled_embedding[relabeling[port]] == embedding[port]
                    for port in range(_PORT_COUNT)
                )
            )
            realization_commutes = False
            if embedding_commutes:
                relabeled_vertex_to_port = {
                    vertex: port for port, vertex in enumerate(relabeled_embedding)
                }
                relabeled_realization = np.zeros_like(realization)
                for face_index, face in enumerate(base.faces):
                    for vertex in face:
                        relabeled_realization[
                            face_index, relabeled_vertex_to_port[int(vertex)]
                        ] += 1.0 / 3.0
                realization_commutes = bool(
                    np.max(
                        np.abs(
                            relabeled_realization[:, np.asarray(relabeling)]
                            - realization
                        )
                    )
                    <= _REFINEMENT_NATURALITY_TOLERANCE
                )
            witness_boundary = ExternalBoundaryBundle(
                boundary_id="refinement-quotient-boundary",
                carrier_id="refinement-quotient-carrier",
                ports=tuple(range(_PORT_COUNT)),
                boundary_condition="open_external",
                boundary_algebra_sha256=interface_algebra_sha256(
                    {"schema": "oph.refinement_naturality.boundary.v1"}
                ),
            )
            witness_federation = EchosahedralFederation(
                federation_id="refinement-quotient-witness",
                carriers=(
                    reference_echosahedral_carrier("refinement-quotient-carrier"),
                ),
                seams=(),
                external_boundaries=(witness_boundary,),
            )
            permutation_map = {"refinement-quotient-carrier": relabeling}
            transformed = relabel_federation_ports(
                witness_federation, permutation_map
            )
            invariance = carrier_quotient_invariance_report(
                witness_federation, transformed, permutation_map
            )
            typed_invariance = bool(
                invariance["CARRIER_QUOTIENT_INVARIANCE_RECEIPT"]
            )
            quotient_checks.append(
                {
                    "relabeling": list(relabeling),
                    "embedding_commutes": embedding_commutes,
                    "realization_commutes": realization_commutes,
                    "typed_quotient_invariance": typed_invariance,
                }
            )
            quotient_verified = bool(
                quotient_verified
                and embedding_commutes
                and realization_commutes
                and typed_invariance
            )
    if not quotient_verified:
        blockers.append("quotient_relabeling_refinement_commutation_failed")

    # (e) Chart realization: the composite port-to-cell map is well defined,
    # A5-equivariant, and rides on receipts for levels 0->1->2.
    realization_verified = bool(
        embedding_supplied
        and refinement_maps_valid
        and identity_verified
        and chart_equivariant
        and seam_law_violation_count == 0
        and family_valid
    )
    if not realization_verified:
        blockers.append("carrier_to_support_chart_realization_checks_failed")

    naturality = bool(
        twelve_port_type
        and embedding_supplied
        and refinement_maps_valid
        and identity_verified
        and seam_law_verified
        and quotient_verified
    )
    chart = bool(naturality and realization_verified)
    return {
        "schema": "oph.echosahedral_federation.refinement_naturality.v2",
        "local_twelve_port_type_declared": twelve_port_type,
        "carrier_embedding_map_supplied": embedding_supplied,
        "carrier_coarse_graining_map_supplied": refinement_maps_valid,
        "coarse_after_embedding_identity_verified": identity_verified,
        "seam_law_naturality_verified": seam_law_verified,
        "quotient_map_commutation_verified": quotient_verified,
        "carrier_to_support_realization_verified": realization_verified,
        "embedding_port_to_defect_vertex": list(declared_embedding),
        "defect_vertex_ids": [int(value) for value in defect_ids],
        "embedding_match_residual": embedding_residual,
        "embedding_match_tolerance": _REFINEMENT_DEFECT_MATCH_TOLERANCE,
        "coarse_after_embedding_identity_residual": identity_residual,
        "identity_tolerance": _REFINEMENT_IDENTITY_TOLERANCE,
        "naturality_tolerance": _REFINEMENT_NATURALITY_TOLERANCE,
        "checked_action_count": len(actions),
        "seam_law_violation_count": seam_law_violation_count,
        "seam_law_violation_examples": action_violation_examples,
        "port_realization_equivariance_residual": (
            realization_equivariance_residual
        ),
        "chart_realization_equivariance_residuals": {
            str(level): value
            for level, value in sorted(chart_equivariance_residuals.items())
        },
        "cell_refinement_receipts": refinement_rows,
        "quotient_commutation_checks": quotient_checks,
        "tower_levels_checked": level_cap + 1,
        "blockers": sorted(set(blockers)),
        "CARRIER_REFINEMENT_NATURALITY_RECEIPT": naturality,
        "CARRIER_TO_SUPPORT_CHART_REALIZATION_RECEIPT": chart,
        "claim_boundary": (
            "This receipt certifies the finite natural transformation between "
            "the twelve-port reference carrier and the exact geodesic tower "
            "through level two: defect-vertex embedding, state-preserving "
            "coarse graining, A5 seam-law naturality, and quotient "
            "commutation. It selects no physical support chart and carries no "
            "H3, event, BW, or KMS claim."
        ),
    }


def echosahedral_federation_receipt(
    federation: EchosahedralFederation,
    *,
    promoted_payload: Mapping[str, Any] | None = None,
    equivalent_presentation: EchosahedralFederation | None = None,
    presentation_port_permutations: Mapping[str, Sequence[int]] | None = None,
) -> dict[str, Any]:
    """Emit the four non-promoting carrier-level parent receipts."""

    sewing = federation_sewing_report(federation)
    firewall = presentation_firewall_report(federation, promoted_payload)
    if (
        equivalent_presentation is not None
        and presentation_port_permutations is not None
    ):
        quotient = carrier_quotient_invariance_report(
            federation,
            equivalent_presentation,
            presentation_port_permutations,
        )
    else:
        quotient = {
            "schema": "oph.echosahedral_federation.quotient_invariance.v1",
            "blockers": ["independent_equivalent_presentation_witness_missing"],
            "CARRIER_QUOTIENT_INVARIANCE_RECEIPT": False,
            "FULL_DYNAMICAL_QUOTIENT_INVARIANCE_RECEIPT": False,
        }
    refinement = carrier_refinement_naturality_report()
    support_limit = controlled_oriented_s2_limit_report(federation)
    carrier_receipt = bool(sewing["ECHOSAHEDRAL_CARRIER_CONFORMANCE"])
    sewing_receipt = bool(sewing["FEDERATION_SEWING_RECEIPT"])
    quotient_receipt = bool(
        quotient["CARRIER_QUOTIENT_INVARIANCE_RECEIPT"]
        and firewall["CARRIER_PRESENTATION_FIREWALL_RECEIPT"]
    )
    refinement_receipt = bool(refinement["CARRIER_REFINEMENT_NATURALITY_RECEIPT"])
    chart_realization_receipt = bool(
        refinement["CARRIER_TO_SUPPORT_CHART_REALIZATION_RECEIPT"]
    )
    homomorphism_receipt = bool(
        sewing["INTERFACE_ALGEBRA_MAP_HOMOMORPHISM_RECEIPT"]
    )
    cocycle_receipt = bool(sewing["HIGHER_OVERLAP_COCYCLE_RECEIPT"])
    nonvacuous_cocycle_witness = bool(
        sewing["NONVACUOUS_HIGHER_OVERLAP_COCYCLE_WITNESS"]
    )
    full_sewing_receipt = bool(sewing["FULL_INTERFACE_ALGEBRA_SEWING_RECEIPT"])
    physical_realization_receipt = bool(
        sewing["PHYSICAL_ECHOSAHEDRAL_FEDERATION_REALIZATION_RECEIPT"]
    )
    source_instrument_valid = bool(
        carrier_receipt
        and sewing_receipt
        and quotient_receipt
        and refinement_receipt
        and full_sewing_receipt
        and physical_realization_receipt
        and support_limit["CONTROLLED_ORIENTED_S2_LIMIT_RECEIPT"]
    )
    s2_support_chart_receipt = bool(
        chart_realization_receipt
        and carrier_receipt
        and sewing_receipt
        and physical_realization_receipt
        and support_limit["CONTROLLED_ORIENTED_S2_LIMIT_RECEIPT"]
    )
    return {
        "schema": "oph.echosahedral_federation.parent_receipts.v1",
        "instrument_scope": "finite_two_level_echosahedral_carrier_federation_contract",
        "federation_id": federation.federation_id,
        "cardinality_semantics": _CARDINALITY_SEMANTICS,
        "carrier_count": len(federation.carriers),
        "exact_source_carrier_count": sewing["exact_source_carrier_count"],
        "support_regulator_count": None,
        "carrier_count_is_support_regulator_count": False,
        "support_chart_cell_count": None,
        "carrier_to_support_chart_realization": (
            "verified_finite_reference_tower_levels_0_2"
            if chart_realization_receipt
            else "unproved"
        ),
        "sewing": sewing,
        "presentation_firewall": firewall,
        "quotient_invariance": quotient,
        "refinement_naturality": refinement,
        "controlled_oriented_s2_limit": support_limit,
        "ECHOSAHEDRAL_CARRIER_CONFORMANCE": carrier_receipt,
        "FEDERATION_SEWING_RECEIPT": sewing_receipt,
        "CARRIER_QUOTIENT_INVARIANCE_RECEIPT": quotient_receipt,
        "CARRIER_REFINEMENT_NATURALITY_RECEIPT": refinement_receipt,
        "INTERFACE_ALGEBRA_MAP_HOMOMORPHISM_RECEIPT": homomorphism_receipt,
        "HIGHER_OVERLAP_COCYCLE_RECEIPT": cocycle_receipt,
        "HIGHER_OVERLAP_COCYCLE_CONDITION_RECEIPT": cocycle_receipt,
        "NONVACUOUS_HIGHER_OVERLAP_COCYCLE_WITNESS": (
            nonvacuous_cocycle_witness
        ),
        "FULL_INTERFACE_ALGEBRA_SEWING_RECEIPT": full_sewing_receipt,
        "PHYSICAL_ECHOSAHEDRAL_FEDERATION_REALIZATION_RECEIPT": (
            physical_realization_receipt
        ),
        "CARRIER_TO_SUPPORT_CHART_REALIZATION_RECEIPT": chart_realization_receipt,
        "INCIDENCE_NERVE_TO_SUPPORT_SIMPLICIAL_ISOMORPHISM_RECEIPT": bool(
            support_limit[
                "INCIDENCE_NERVE_TO_SUPPORT_SIMPLICIAL_ISOMORPHISM_RECEIPT"
            ]
        ),
        "CONTROLLED_ORIENTED_S2_LIMIT_RECEIPT": bool(
            support_limit["CONTROLLED_ORIENTED_S2_LIMIT_RECEIPT"]
        ),
        "ECHOSAHEDRAL_FEDERATION_SOURCE_INSTRUMENT_VALID": source_instrument_valid,
        "S2_SUPPORT_CHART_EMERGENCE_RECEIPT": s2_support_chart_receipt,
        "H3_FRAME_EMERGENCE_RECEIPT": False,
        "EVENT_MANIFOLD_RECEIPT": False,
        "BW_KMS_CLOCK_RECEIPT": False,
        "PHYSICAL_H3_KMS_EMERGENCE_RECEIPT": False,
        "claim_boundary": (
            "These parent receipts stop at the declared finite carrier contract "
            "and its exact refinement bridge. carrier_count is the exact "
            "source-federation cardinality when unique IDs verify, and is "
            "separate from the support regulator, S2 chart-cell count, N_star, "
            "observer count, H3-point count, and event count. The interface "
            "algebra map, logical higher-overlap cocycle condition, refinement "
            "naturality, and "
            "support-chart realization receipts are computed finite checks on "
            "the reference structures. The oriented-S2 promotion additionally "
            "requires the full source-incidence nerve and the analytic "
            "mesh-contraction theorem, not merely a finite chart embedding. "
            "A separate flag records whether a nonvacuous explicitly declared "
            "triple restriction was present. H3 frames, events, BW/KMS clocks, and "
            "the 2pi normalization require independent downstream receipts."
        ),
    }


def reference_federation_instrument_bundle(
    federation: EchosahedralFederation,
    *,
    promoted_payload: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Serialize a compact shared-template bundle without duplicating XYZ data.

    This exporter intentionally accepts only carriers in the canonical shared
    reference presentation.  Presentation-relabeling witnesses use the typed
    API above; a compact bundle must first be canonicalized explicitly.
    """

    if any(
        carrier != reference_echosahedral_carrier(carrier.carrier_id)
        for carrier in federation.carriers
    ):
        raise ValueError(
            "reference instrument bundles require canonical shared-template carriers"
        )
    return {
        "schema": "oph.echosahedral_federation.instrument_bundle.v1",
        "federation_id": federation.federation_id,
        "cardinality_semantics": _CARDINALITY_SEMANTICS,
        "exact_source_carrier_count": len(federation.carriers),
        "support_regulator_count": None,
        "local_carrier_template": "regular_icosahedron_12_30_20_antipode_a5_v1",
        "local_carrier_template_structural_class_sha256": (
            echosahedral_carrier_conformance_report(
                _reference_echosahedral_carrier_template()
            )["structural_class_sha256"]
        ),
        "interface_hash_check_scope": _INTERFACE_HASH_CHECK_SCOPE,
        "carrier_ids": [carrier.carrier_id for carrier in federation.carriers],
        "seams": [
            {
                "seam_id": seam.seam_id,
                "left_carrier_id": seam.left_carrier_id,
                "right_carrier_id": seam.right_carrier_id,
                "left_ports": list(seam.left_ports),
                "right_ports": list(seam.right_ports),
                "left_to_right_ports": list(seam.left_to_right_ports),
                "right_to_left_ports": list(seam.right_to_left_ports),
                "left_to_right_orientation": list(seam.left_to_right_orientation),
                "right_to_left_orientation": list(seam.right_to_left_orientation),
                "collar_kind": seam.collar_kind,
                "interface_algebra": {
                    "interface_algebra_id": (
                        seam.interface_algebra.interface_algebra_id
                    ),
                    "interface_algebra_sha256": (
                        seam.interface_algebra.interface_algebra_sha256
                    ),
                    "left_interface_algebra_sha256": (
                        seam.interface_algebra.left_interface_algebra_sha256
                    ),
                    "right_interface_algebra_sha256": (
                        seam.interface_algebra.right_interface_algebra_sha256
                    ),
                },
            }
            for seam in federation.seams
        ],
        "triple_overlaps": [
            {
                "overlap_id": overlap.overlap_id,
                "oriented_carrier_ids": list(overlap.oriented_carrier_ids),
                "oriented_seam_ids": list(overlap.oriented_seam_ids),
                "restriction_algebra_id": overlap.restriction_algebra_id,
                "restriction_algebra_sha256": (
                    overlap.restriction_algebra_sha256
                ),
            }
            for overlap in federation.triple_overlaps
        ],
        "support_tower_binding": (
            None
            if federation.support_tower_binding is None
            else {
                "geometry_family": (
                    federation.support_tower_binding.geometry_family
                ),
                "carrier_to_base_vertex": [
                    [carrier_id, vertex]
                    for carrier_id, vertex in (
                        federation.support_tower_binding.carrier_to_base_vertex
                    )
                ],
                "orientation": federation.support_tower_binding.orientation,
                "refinement_rule": (
                    federation.support_tower_binding.refinement_rule
                ),
                "source_incidence_sha256": (
                    federation.support_tower_binding.source_incidence_sha256
                ),
            }
        ),
        "presentation_relabeling_witness": (
            "reverse_every_local_port_order_v1"
            if federation.support_tower_binding is not None
            else None
        ),
        "external_boundaries": [
            {
                "boundary_id": boundary.boundary_id,
                "carrier_id": boundary.carrier_id,
                "ports": list(boundary.ports),
                "boundary_condition": boundary.boundary_condition,
                "boundary_algebra_sha256": boundary.boundary_algebra_sha256,
            }
            for boundary in federation.external_boundaries
        ],
        "observer_supports": [
            {
                "observer_token": support.observer_token,
                "carrier_ids": sorted(support.carrier_ids),
                "visible_seam_ids": sorted(support.visible_seam_ids),
                "record_algebra_sha256": support.record_algebra_sha256,
                "checkpoint_cut_sha256": support.checkpoint_cut_sha256,
            }
            for support in federation.observer_supports
        ],
        "promoted_payload": {} if promoted_payload is None else promoted_payload,
    }


def verify_reference_federation_instrument_bundle(
    bundle: Mapping[str, Any],
) -> dict[str, Any]:
    """Parse and verify the compact JSON-facing carrier instrument bundle."""

    try:
        if not isinstance(bundle, Mapping):
            raise TypeError("instrument bundle must be an object")
        if bundle.get("schema") != "oph.echosahedral_federation.instrument_bundle.v1":
            raise ValueError("instrument bundle schema mismatch")
        if (
            bundle.get("local_carrier_template")
            != "regular_icosahedron_12_30_20_antipode_a5_v1"
        ):
            raise ValueError("unknown local carrier template")
        expected_template_hash = echosahedral_carrier_conformance_report(
            _reference_echosahedral_carrier_template()
        )["structural_class_sha256"]
        if (
            bundle.get(
                "local_carrier_template_structural_class_sha256",
                expected_template_hash,
            )
            != expected_template_hash
        ):
            raise ValueError("local carrier template structural hash mismatch")
        if bundle.get("interface_hash_check_scope", _INTERFACE_HASH_CHECK_SCOPE) != (
            _INTERFACE_HASH_CHECK_SCOPE
        ):
            raise ValueError("interface_hash_check_scope mismatch")
        federation_id = bundle["federation_id"]
        if not isinstance(federation_id, str) or not federation_id:
            raise ValueError("federation_id must be a nonempty string")
        raw_carrier_ids = bundle["carrier_ids"]
        if isinstance(raw_carrier_ids, (str, bytes)) or not isinstance(
            raw_carrier_ids, Sequence
        ):
            raise TypeError("carrier_ids must be an array of strings")
        carrier_ids = tuple(raw_carrier_ids)
        if not all(isinstance(value, str) and value for value in carrier_ids):
            raise TypeError("carrier_ids must contain only nonempty strings")
        if not carrier_ids or len(set(carrier_ids)) != len(carrier_ids):
            raise ValueError("carrier_ids must be nonempty and unique")
        declared_count = bundle.get("exact_source_carrier_count", len(carrier_ids))
        if type(declared_count) is not int or declared_count != len(carrier_ids):
            raise ValueError("exact_source_carrier_count does not match carrier_ids")
        declared_semantics = bundle.get("cardinality_semantics", _CARDINALITY_SEMANTICS)
        if declared_semantics != _CARDINALITY_SEMANTICS:
            raise ValueError("cardinality_semantics mismatch")
        support_regulator_count = bundle.get("support_regulator_count")
        if support_regulator_count is not None:
            raise ValueError(
                "support_regulator_count must remain separate and null in this bundle"
            )
        instrument_surface = {
            key: value for key, value in bundle.items() if key != "promoted_payload"
        }
        embedded_template_fields = sorted(
            _recursive_mapping_keys(instrument_surface)
            & _EMBEDDED_LOCAL_TEMPLATE_FIELDS
        )
        if embedded_template_fields:
            raise ValueError(
                "per-carrier local template fields are forbidden:"
                + ",".join(embedded_template_fields)
            )
        carriers = tuple(reference_echosahedral_carrier(value) for value in carrier_ids)
        seams = tuple(_seam_from_bundle_row(row) for row in bundle.get("seams", ()))
        triple_overlaps = tuple(
            _triple_overlap_from_bundle_row(row)
            for row in bundle.get("triple_overlaps", ())
        )
        support_binding = _support_binding_from_bundle_value(
            bundle.get("support_tower_binding")
        )
        boundaries = tuple(
            _boundary_from_bundle_row(row)
            for row in bundle.get("external_boundaries", ())
        )
        supports = tuple(
            _observer_from_bundle_row(row)
            for row in bundle.get("observer_supports", ())
        )
        federation = EchosahedralFederation(
            federation_id=federation_id,
            carriers=carriers,
            seams=seams,
            external_boundaries=boundaries,
            observer_supports=supports,
            triple_overlaps=triple_overlaps,
            support_tower_binding=support_binding,
        )
        witness_rule = bundle.get("presentation_relabeling_witness")
        equivalent_presentation = None
        presentation_port_permutations = None
        if witness_rule is not None:
            if witness_rule != "reverse_every_local_port_order_v1":
                raise ValueError("unknown presentation_relabeling_witness")
            presentation_port_permutations = {
                carrier_id: tuple(reversed(range(_PORT_COUNT)))
                for carrier_id in carrier_ids
            }
            equivalent_presentation = relabel_federation_ports(
                federation, presentation_port_permutations
            )
        report = echosahedral_federation_receipt(
            federation,
            promoted_payload=bundle.get("promoted_payload", {}),
            equivalent_presentation=equivalent_presentation,
            presentation_port_permutations=presentation_port_permutations,
        )
    except (KeyError, TypeError, ValueError) as exc:
        return {
            "schema": "oph.echosahedral_federation.instrument_bundle_verification.v1",
            "parse_error": f"{type(exc).__name__}:{exc}",
            "INSTRUMENT_BUNDLE_SCHEMA_RECEIPT": False,
            "ECHOSAHEDRAL_CARRIER_CONFORMANCE": False,
            "FEDERATION_SEWING_RECEIPT": False,
            "CARRIER_QUOTIENT_INVARIANCE_RECEIPT": False,
            "CARRIER_REFINEMENT_NATURALITY_RECEIPT": False,
            "INTERFACE_ALGEBRA_MAP_HOMOMORPHISM_RECEIPT": False,
            "HIGHER_OVERLAP_COCYCLE_RECEIPT": False,
            "FULL_INTERFACE_ALGEBRA_SEWING_RECEIPT": False,
            "PHYSICAL_ECHOSAHEDRAL_FEDERATION_REALIZATION_RECEIPT": False,
            "CARRIER_TO_SUPPORT_CHART_REALIZATION_RECEIPT": False,
            "ECHOSAHEDRAL_FEDERATION_SOURCE_INSTRUMENT_VALID": False,
        }
    report = dict(report)
    report["bundle_schema"] = bundle["schema"]
    report["cardinality_semantics"] = _CARDINALITY_SEMANTICS
    report["exact_source_carrier_count"] = len(carrier_ids)
    report["support_regulator_count"] = None
    report["shared_template_encoded_once"] = True
    report["per_carrier_coordinate_tables_embedded"] = False
    report["shared_template_conformance_verified_once"] = bool(
        report["sewing"]["carrier_conformance_summary"][
            "shared_reference_template_conformance_verified_once"
        ]
    )
    report["INSTRUMENT_BUNDLE_SCHEMA_RECEIPT"] = True
    return report


def screen_port_map_carrier_bridge_report(
    screen_port_report: Mapping[str, Any],
) -> dict[str, Any]:
    """Audit current ScreenPortMap output without promoting it to full sewing."""

    local = bool(
        screen_port_report.get("ECHOSAHEDRAL_CARRIER_CONFORMANCE_RECEIPT") is True
        and screen_port_report.get("ports_per_patch") == 12
    )
    singleton_reference = bool(
        screen_port_report.get("REFERENCE_SINGLETON_FEDERATION_SEWING_RECEIPT") is True
    )
    federation_section = screen_port_report.get("federation_sewing", {})
    general_bundle_schema = bool(
        isinstance(federation_section, Mapping)
        and federation_section.get("general_bundle_schema_implemented") is True
    )
    interface_hashes = bool(
        isinstance(federation_section, Mapping)
        and federation_section.get("interface_algebra_hashes_bound") is True
    )
    refinement = carrier_refinement_naturality_report()
    refinement_naturality = bool(
        refinement["CARRIER_REFINEMENT_NATURALITY_RECEIPT"]
    )
    chart_realization = bool(
        refinement["CARRIER_TO_SUPPORT_CHART_REALIZATION_RECEIPT"]
    )
    blockers = [
        item
        for condition, item in (
            (not local, "local_twelve_port_carrier_conformance_missing"),
            (not singleton_reference, "singleton_reference_sewing_missing"),
            (not general_bundle_schema, "typed_general_seam_bundle_schema_missing"),
            (not interface_hashes, "interface_algebra_hash_binding_missing"),
            (True, "explicit_external_boundary_bundle_ledger_missing"),
            (True, "connected_observer_support_ledger_missing"),
            (not chart_realization, "carrier_to_support_chart_realization_missing"),
            (not refinement_naturality, "carrier_refinement_naturality_missing"),
        )
        if condition
    ]
    return {
        "schema": "oph.echosahedral_federation.screen_port_map_bridge.v1",
        "bridge_scope": "existing_screen_port_map_reference_audit_only",
        "local_carrier_conformance_imported": local,
        "singleton_reference_sewing_imported": singleton_reference,
        "typed_general_bundle_schema_present": general_bundle_schema,
        "interface_algebra_hashes_bound": interface_hashes,
        "blockers": blockers,
        "REFERENCE_SCREEN_PORT_MAP_CARRIER_BRIDGE_RECEIPT": bool(
            local and singleton_reference
        ),
        "ECHOSAHEDRAL_CARRIER_CONFORMANCE": local,
        "FEDERATION_SEWING_RECEIPT": False,
        "CARRIER_QUOTIENT_INVARIANCE_RECEIPT": False,
        "CARRIER_REFINEMENT_NATURALITY_RECEIPT": refinement_naturality,
        "CARRIER_TO_SUPPORT_CHART_REALIZATION_RECEIPT": chart_realization,
        "ECHOSAHEDRAL_FEDERATION_SOURCE_INSTRUMENT_VALID": False,
        "S2_SUPPORT_CHART_EMERGENCE_RECEIPT": False,
        "H3_FRAME_EMERGENCE_RECEIPT": False,
        "EVENT_MANIFOLD_RECEIPT": False,
        "BW_KMS_CLOCK_RECEIPT": False,
        "claim_boundary": (
            "The current engine report may certify the shared local carrier "
            "template and singleton routing reference, and it imports the "
            "recomputed reference-template refinement naturality. Typed collar "
            "bundles, interface hashes, explicit external boundaries, observer "
            "supports, and quotient dynamics are open items for this engine."
        ),
    }


def main(argv: Sequence[str] | None = None) -> int:
    """Verify a compact bundle with ``python -m ...echosahedral_federation``."""

    parser = argparse.ArgumentParser(
        description="Verify an OPH echosahedral federation instrument bundle"
    )
    parser.add_argument("bundle", type=Path, help="JSON instrument bundle path")
    parser.add_argument(
        "--require",
        choices=("schema", "carrier", "sewing", "quotient", "refinement", "source"),
        default="sewing",
        help="receipt tier controlling the process exit status",
    )
    args = parser.parse_args(argv)
    try:
        payload = json.loads(args.bundle.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        report = {
            "schema": "oph.echosahedral_federation.instrument_bundle_verification.v1",
            "parse_error": f"{type(exc).__name__}:{exc}",
            "INSTRUMENT_BUNDLE_SCHEMA_RECEIPT": False,
        }
    else:
        report = verify_reference_federation_instrument_bundle(payload)
    print(json.dumps(report, sort_keys=True, indent=2))
    receipt_by_tier = {
        "schema": "INSTRUMENT_BUNDLE_SCHEMA_RECEIPT",
        "carrier": "ECHOSAHEDRAL_CARRIER_CONFORMANCE",
        "sewing": "FEDERATION_SEWING_RECEIPT",
        "quotient": "CARRIER_QUOTIENT_INVARIANCE_RECEIPT",
        "refinement": "CARRIER_REFINEMENT_NATURALITY_RECEIPT",
        "source": "ECHOSAHEDRAL_FEDERATION_SOURCE_INSTRUMENT_VALID",
    }
    return 0 if report.get(receipt_by_tier[args.require]) is True else 1


def _a5_local_action_audit(
    actions: Sequence[Sequence[int]],
    *,
    edge_set: set[tuple[int, int]],
    oriented_face_set: set[tuple[int, int, int]],
    antipode: Sequence[int],
    coordinates: np.ndarray,
    tolerance: float,
) -> dict[str, Any]:
    blockers: list[str] = []
    permutations = tuple(tuple(int(value) for value in row) for row in actions)
    shape_valid = bool(
        len(permutations) == 60
        and all(
            len(row) == _PORT_COUNT and sorted(row) == list(range(_PORT_COUNT))
            for row in permutations
        )
        and len(set(permutations)) == 60
    )
    if not shape_valid:
        blockers.append("a5_action_is_not_60_unique_twelve_port_permutations")
    identity = tuple(range(_PORT_COUNT))
    group_closed = False
    inverses_present = False
    order_profile: dict[int, int] = {}
    if shape_valid:
        action_set = set(permutations)
        group_closed = all(
            _compose_permutations(left, right) in action_set
            for left in permutations
            for right in permutations
        )
        inverses_present = all(
            _inverse_permutation(row) in action_set for row in permutations
        )
        for row in permutations:
            try:
                order = _permutation_order(row)
            except ValueError:
                blockers.append("registered_permutation_order_exceeds_a5_bound")
                continue
            order_profile[order] = order_profile.get(order, 0) + 1
    expected_order_profile = {1: 1, 2: 15, 3: 20, 5: 24}
    if not group_closed:
        blockers.append("a5_permutation_family_not_closed")
    if not inverses_present:
        blockers.append("a5_permutation_inverses_missing")
    if identity not in set(permutations) or order_profile != expected_order_profile:
        blockers.append("a5_element_order_profile_invalid")

    edge_preserved = bool(shape_valid)
    face_preserved = bool(shape_valid)
    antipode_commutes = bool(
        shape_valid
        and len(antipode) == _PORT_COUNT
        and sorted(int(value) for value in antipode) == list(range(_PORT_COUNT))
    )
    maximum_gram_residual = math.inf
    maximum_rotation_residual = math.inf
    minimum_rotation_determinant = -math.inf
    if shape_valid:
        gram = coordinates @ coordinates.T
        gram_residuals: list[float] = []
        rotation_residuals: list[float] = []
        determinants: list[float] = []
        for row in permutations:
            edge_preserved = edge_preserved and all(
                tuple(sorted((row[left], row[right]))) in edge_set
                for left, right in edge_set
            )
            face_preserved = face_preserved and all(
                _cyclic_face_key((row[a], row[b], row[c])) in oriented_face_set
                for a, b, c in oriented_face_set
            )
            if antipode_commutes:
                antipode_commutes = all(
                    row[int(antipode[index])] == int(antipode[row[index]])
                    for index in range(_PORT_COUNT)
                )
            indexed = np.asarray(row, dtype=int)
            gram_residuals.append(
                float(np.max(np.abs(gram - gram[np.ix_(indexed, indexed)])))
            )
            rotation, residual = _proper_rotation_fit(coordinates, coordinates[indexed])
            determinants.append(float(np.linalg.det(rotation)))
            rotation_residuals.append(residual)
        maximum_gram_residual = max(gram_residuals, default=math.inf)
        maximum_rotation_residual = max(rotation_residuals, default=math.inf)
        minimum_rotation_determinant = min(determinants, default=-math.inf)
    if not edge_preserved:
        blockers.append("a5_action_does_not_preserve_all_edges")
    if not face_preserved:
        blockers.append("a5_action_does_not_preserve_all_oriented_faces")
    if not antipode_commutes:
        blockers.append("a5_action_does_not_commute_with_antipode")
    if maximum_gram_residual > tolerance:
        blockers.append("a5_action_does_not_preserve_coordinate_gram_matrix")
    if maximum_rotation_residual > tolerance or minimum_rotation_determinant <= 0.0:
        blockers.append("a5_action_is_not_realized_by_proper_coordinate_rotations")
    return {
        "registered_action_count": len(permutations),
        "unique_action_count": len(set(permutations)),
        "group_closed": group_closed,
        "inverses_present": inverses_present,
        "element_order_profile": {
            str(order): count for order, count in sorted(order_profile.items())
        },
        "all_actions_preserve_edges": edge_preserved,
        "all_actions_preserve_oriented_faces": face_preserved,
        "all_actions_commute_with_antipode": antipode_commutes,
        "maximum_coordinate_gram_residual": maximum_gram_residual,
        "maximum_proper_rotation_fit_residual": maximum_rotation_residual,
        "minimum_rotation_determinant": minimum_rotation_determinant,
        "blockers": sorted(set(blockers)),
        "receipt": not blockers,
        "A5_ORDER_60_LOCAL_ACTION_RECEIPT": not blockers,
    }


def _reference_icosahedron_isomorphism_receipt(
    coordinates: np.ndarray,
    *,
    edge_set: set[tuple[int, int]],
    oriented_face_set: set[tuple[int, int, int]],
    antipode: Sequence[int],
    tolerance: float,
) -> dict[str, Any]:
    reference = reference_echosahedral_carrier("reference")
    reference_coordinates = np.asarray(reference.port_coordinates, dtype=float)
    reference_edges = {tuple(sorted(edge)) for edge in reference.edges}
    reference_faces = {_cyclic_face_key(face) for face in reference.faces}
    candidate_graph = nx.Graph()
    candidate_graph.add_nodes_from(range(_PORT_COUNT))
    candidate_graph.add_edges_from(edge_set)
    reference_graph = nx.Graph()
    reference_graph.add_nodes_from(range(_PORT_COUNT))
    reference_graph.add_edges_from(reference_edges)
    valid_maps = 0
    minimum_coordinate_residual = math.inf
    if len(antipode) == _PORT_COUNT and sorted(
        int(value) for value in antipode
    ) == list(range(_PORT_COUNT)):
        matcher = nx.algorithms.isomorphism.GraphMatcher(
            candidate_graph, reference_graph
        )
        for mapping in matcher.isomorphisms_iter():
            if not all(
                _cyclic_face_key((mapping[face[0]], mapping[face[1]], mapping[face[2]]))
                in reference_faces
                for face in oriented_face_set
            ):
                continue
            if not all(
                mapping[int(antipode[index])] == reference.antipode[mapping[index]]
                for index in range(_PORT_COUNT)
            ):
                continue
            target = np.zeros_like(reference_coordinates)
            for candidate_port, reference_port in mapping.items():
                target[candidate_port] = reference_coordinates[reference_port]
            rotation, residual = _proper_rotation_fit(coordinates, target)
            if np.linalg.det(rotation) > 0.0 and residual <= tolerance:
                valid_maps += 1
                minimum_coordinate_residual = min(minimum_coordinate_residual, residual)
    receipt = valid_maps > 0
    return {
        "proper_oriented_reference_isomorphism_count": valid_maps,
        "minimum_coordinate_rotation_residual": minimum_coordinate_residual,
        "receipt": receipt,
        "EXACT_REFERENCE_ICOSAHEDRON_ISOMORPHISM_RECEIPT": receipt,
    }


def _proper_rotation_fit(
    source: np.ndarray, target: np.ndarray
) -> tuple[np.ndarray, float]:
    covariance = np.asarray(source, dtype=float).T @ np.asarray(target, dtype=float)
    left, _, right_t = np.linalg.svd(covariance)
    rotation = left @ right_t
    if float(np.linalg.det(rotation)) < 0.0:
        left[:, -1] *= -1.0
        rotation = left @ right_t
    residual = float(
        np.max(
            np.linalg.norm(
                np.asarray(source, dtype=float) @ rotation
                - np.asarray(target, dtype=float),
                axis=1,
            )
        )
    )
    return rotation, residual


def _eigenvalue_multiplicities(
    eigenvalues: np.ndarray, *, tolerance: float
) -> list[tuple[float, int]]:
    result: list[tuple[float, int]] = []
    for value in np.sort(np.asarray(eigenvalues, dtype=float)):
        if result and abs(value - result[-1][0]) <= tolerance:
            mean, count = result[-1]
            result[-1] = ((mean * count + float(value)) / (count + 1), count + 1)
        else:
            result.append((float(value), 1))
    return result


def _cyclic_face_key(face: Sequence[int]) -> tuple[int, int, int]:
    a, b, c = (int(value) for value in face)
    return min((a, b, c), (b, c, a), (c, a, b))


def _compose_permutations(left: Sequence[int], right: Sequence[int]) -> tuple[int, ...]:
    return tuple(int(left[int(right[index])]) for index in range(len(left)))


def _inverse_permutation(permutation: Sequence[int]) -> tuple[int, ...]:
    inverse = [0] * len(permutation)
    for source, target in enumerate(permutation):
        inverse[int(target)] = source
    return tuple(inverse)


def _permutation_order(permutation: Sequence[int]) -> int:
    identity = tuple(range(len(permutation)))
    current = identity
    for order in range(1, 61):
        current = _compose_permutations(permutation, current)
        if current == identity:
            return order
    raise ValueError("local permutation order exceeds A5 bound")


def _port_subset_connected(carrier: EchosahedralCarrier, ports: Sequence[int]) -> bool:
    port_set = {int(port) for port in ports}
    if not port_set or not port_set <= set(range(_PORT_COUNT)):
        return False
    pending = [next(iter(port_set))]
    visited: set[int] = set()
    while pending:
        port = pending.pop()
        if port in visited:
            continue
        visited.add(port)
        for left, right in carrier.edges:
            if left == port and right in port_set and right not in visited:
                pending.append(right)
            elif right == port and left in port_set and left not in visited:
                pending.append(left)
    return visited == port_set


@lru_cache(maxsize=1)
def _reference_incidence_sets() -> tuple[
    frozenset[tuple[int, int]], frozenset[tuple[int, int, int]]
]:
    template = _reference_echosahedral_carrier_template()
    return (
        frozenset(tuple(sorted(edge)) for edge in template.edges),
        frozenset(tuple(sorted(face)) for face in template.faces),
    )


def _collar_kind_matches(
    carrier: EchosahedralCarrier,
    ports: Sequence[int],
    collar_kind: str,
) -> bool:
    port_tuple = tuple(int(port) for port in ports)
    port_set = set(port_tuple)
    if not port_tuple or not port_set <= set(range(_PORT_COUNT)):
        return False
    if _uses_shared_reference_template(carrier):
        edge_set, face_set = _reference_incidence_sets()
    else:
        edge_set = frozenset(tuple(sorted(edge)) for edge in carrier.edges)
        face_set = frozenset(tuple(sorted(face)) for face in carrier.faces)
    if collar_kind == "single_port":
        return len(port_tuple) == 1
    if collar_kind == "antipodal_pair":
        return bool(
            len(port_tuple) == 2
            and len(carrier.antipode) == _PORT_COUNT
            and carrier.antipode[port_tuple[0]] == port_tuple[1]
            and carrier.antipode[port_tuple[1]] == port_tuple[0]
        )
    if collar_kind == "edge_bundle":
        return len(port_tuple) == 2 and tuple(sorted(port_tuple)) in edge_set
    if collar_kind == "face_collar":
        return len(port_tuple) == 3 and tuple(sorted(port_tuple)) in face_set
    if collar_kind == "connected_bundle":
        return bool(port_set and _port_subset_connected(carrier, port_tuple))
    return False


def _validated_port_permutation(permutation: Sequence[int]) -> tuple[int, ...]:
    normalized = tuple(int(value) for value in permutation)
    if len(normalized) != _PORT_COUNT or sorted(normalized) != list(range(_PORT_COUNT)):
        raise ValueError("port permutation must be a bijection of 0..11")
    return normalized


def _relabel_carrier(
    carrier: EchosahedralCarrier,
    permutation: Sequence[int],
) -> EchosahedralCarrier:
    mapping = _validated_port_permutation(permutation)
    inverse = _inverse_permutation(mapping)
    names: list[str] = [""] * _PORT_COUNT
    coordinates: list[tuple[float, float, float]] = [(0.0, 0.0, 0.0)] * _PORT_COUNT
    antipode = [0] * _PORT_COUNT
    for old, new in enumerate(mapping):
        names[new] = carrier.port_names[old]
        coordinates[new] = carrier.port_coordinates[old]
        antipode[new] = mapping[carrier.antipode[old]]
    actions = []
    for action in carrier.a5_actions:
        conjugated = tuple(
            mapping[action[inverse[new_port]]] for new_port in range(_PORT_COUNT)
        )
        actions.append(conjugated)
    return replace(
        carrier,
        port_names=tuple(names),
        port_coordinates=tuple(coordinates),
        edges=tuple((mapping[left], mapping[right]) for left, right in carrier.edges),
        faces=tuple(
            (mapping[first], mapping[second], mapping[third])
            for first, second, third in carrier.faces
        ),
        antipode=tuple(antipode),
        a5_actions=tuple(actions),
    )


def _quotient_visible_contract_payload(
    federation: EchosahedralFederation,
) -> dict[str, Any]:
    """Return compact digests of the quotient-visible contract rows."""

    carrier_digest = _CanonicalRowDigest()
    for carrier in sorted(federation.carriers, key=lambda item: item.carrier_id):
        report = echosahedral_carrier_conformance_report(carrier)
        carrier_digest.update(
            {
                "carrier_id": carrier.carrier_id,
                "structural_class_sha256": report["structural_class_sha256"],
                "conforming": report["ECHOSAHEDRAL_CARRIER_CONFORMANCE"],
            }
        )
    seam_digest = _CanonicalRowDigest()
    for seam in sorted(federation.seams, key=lambda item: item.seam_id):
        seam_digest.update(
            {
                "seam_id": seam.seam_id,
                "endpoint_carrier_ids": sorted(
                    (seam.left_carrier_id, seam.right_carrier_id)
                ),
                "collar_kind": seam.collar_kind,
                "bundle_size": len(seam.left_ports),
                "orientation_reversing": bool(
                    all(sign == -1 for sign in seam.left_to_right_orientation)
                    and all(sign == -1 for sign in seam.right_to_left_orientation)
                ),
                "interface_algebra_id": seam.interface_algebra.interface_algebra_id,
                "interface_algebra_sha256": seam.interface_algebra.interface_algebra_sha256,
            }
        )
    boundary_digest = _CanonicalRowDigest()
    for boundary in sorted(
        federation.external_boundaries, key=lambda item: item.boundary_id
    ):
        boundary_digest.update(
            {
                "boundary_id": boundary.boundary_id,
                "carrier_id": boundary.carrier_id,
                "port_count": len(boundary.ports),
                "boundary_condition": boundary.boundary_condition,
                "boundary_algebra_sha256": boundary.boundary_algebra_sha256,
            }
        )
    observer_digest = _CanonicalRowDigest()
    for support in sorted(
        federation.observer_supports, key=lambda item: item.observer_token
    ):
        observer_digest.update(
            {
                "observer_token": support.observer_token,
                "carrier_ids": sorted(support.carrier_ids),
                "visible_seam_ids": sorted(support.visible_seam_ids),
                "record_algebra_sha256": support.record_algebra_sha256,
                "checkpoint_cut_sha256": support.checkpoint_cut_sha256,
            }
        )
    triple_digest = _CanonicalRowDigest()
    for overlap in sorted(
        federation.triple_overlaps, key=lambda item: item.overlap_id
    ):
        triple_digest.update(
            {
                "overlap_id": overlap.overlap_id,
                "oriented_carrier_ids": list(overlap.oriented_carrier_ids),
                "oriented_seam_ids": list(overlap.oriented_seam_ids),
                "restriction_algebra_id": overlap.restriction_algebra_id,
                "restriction_algebra_sha256": (
                    overlap.restriction_algebra_sha256
                ),
            }
        )
    support_binding = federation.support_tower_binding
    return {
        "schema": "oph.echosahedral_federation.quotient_visible_contract.v1",
        "federation_id": federation.federation_id,
        "cardinality_semantics": _CARDINALITY_SEMANTICS,
        "exact_source_carrier_count": len(federation.carriers),
        "support_regulator_count": None,
        "carrier_row_count": carrier_digest.row_count,
        "carrier_rows_sha256": carrier_digest.hexdigest(),
        "seam_row_count": seam_digest.row_count,
        "seam_rows_sha256": seam_digest.hexdigest(),
        "external_boundary_row_count": boundary_digest.row_count,
        "external_boundary_rows_sha256": boundary_digest.hexdigest(),
        "observer_support_row_count": observer_digest.row_count,
        "observer_support_rows_sha256": observer_digest.hexdigest(),
        "triple_overlap_row_count": triple_digest.row_count,
        "triple_overlap_rows_sha256": triple_digest.hexdigest(),
        "support_tower_binding": (
            None
            if support_binding is None
            else {
                "geometry_family": support_binding.geometry_family,
                "carrier_to_base_vertex": [
                    [carrier_id, vertex]
                    for carrier_id, vertex in support_binding.carrier_to_base_vertex
                ],
                "orientation": support_binding.orientation,
                "refinement_rule": support_binding.refinement_rule,
                "source_incidence_sha256": support_binding.source_incidence_sha256,
            }
        ),
    }


def _recursive_mapping_keys(value: Any) -> set[str]:
    keys: set[str] = set()
    if isinstance(value, Mapping):
        for key, child in value.items():
            keys.add(str(key))
            keys.update(_recursive_mapping_keys(child))
    elif isinstance(value, (list, tuple)):
        for child in value:
            keys.update(_recursive_mapping_keys(child))
    return keys


def _walk_payload_values(value: Any, path: str = "$") -> list[tuple[str, Any]]:
    rows = [(path, value)]
    if isinstance(value, Mapping):
        for key, child in value.items():
            rows.extend(_walk_payload_values(child, f"{path}.{key}"))
    elif isinstance(value, (list, tuple)):
        for index, child in enumerate(value):
            rows.extend(_walk_payload_values(child, f"{path}[{index}]"))
    return rows


def _find_hidden_coordinate_values(
    payload: Mapping[str, Any], carriers: Sequence[EchosahedralCarrier]
) -> list[str]:
    # Canonical 256k carriers share one tuple by identity. Deduplicate on that
    # O(1) identity before computing the comparatively expensive value signature.
    coordinate_tables_by_identity: dict[
        int, tuple[tuple[float, float, float], ...]
    ] = {}
    for carrier in carriers:
        coordinate_tables_by_identity.setdefault(
            id(carrier.port_coordinates), carrier.port_coordinates
        )
    reference_signatures = {
        _coordinate_set_signature(np.asarray(coordinates, dtype=float))
        for coordinates in coordinate_tables_by_identity.values()
    }
    leaks: list[str] = []
    for path, value in _walk_payload_values(payload):
        if isinstance(value, (str, bytes, Mapping)):
            continue
        try:
            array = np.asarray(value, dtype=float)
        except (TypeError, ValueError):
            continue
        if (
            array.shape == (_PORT_COUNT, 3)
            and _coordinate_set_signature(array) in reference_signatures
        ):
            leaks.append(path)
    return sorted(set(leaks))


def _coordinate_set_signature(coordinates: np.ndarray) -> tuple[tuple[float, ...], ...]:
    rounded = np.round(np.asarray(coordinates, dtype=float), decimals=12)
    return tuple(sorted(tuple(float(value) for value in row) for row in rounded))


def _find_hidden_port_label_values(
    payload: Mapping[str, Any], carriers: Sequence[EchosahedralCarrier]
) -> list[str]:
    label_tables_by_identity: dict[int, tuple[str, ...]] = {}
    for carrier in carriers:
        label_tables_by_identity.setdefault(id(carrier.port_names), carrier.port_names)
    label_sets = {
        frozenset(port_names) for port_names in label_tables_by_identity.values()
    }
    leaks: list[str] = []
    for path, value in _walk_payload_values(payload):
        if not isinstance(value, (list, tuple)):
            continue
        if value and all(isinstance(item, str) for item in value):
            if frozenset(value) in label_sets:
                leaks.append(path)
    return sorted(set(leaks))


def _seam_from_bundle_row(row: Mapping[str, Any]) -> SeamBundle:
    if not isinstance(row, Mapping):
        raise TypeError("seam row must be an object")
    binding = row["interface_algebra"]
    if not isinstance(binding, Mapping):
        raise TypeError("interface_algebra must be an object")
    return SeamBundle(
        seam_id=_bundle_string(row, "seam_id"),
        left_carrier_id=_bundle_string(row, "left_carrier_id"),
        right_carrier_id=_bundle_string(row, "right_carrier_id"),
        left_ports=_bundle_int_tuple(row, "left_ports"),
        right_ports=_bundle_int_tuple(row, "right_ports"),
        left_to_right_ports=_bundle_int_tuple(row, "left_to_right_ports"),
        right_to_left_ports=_bundle_int_tuple(row, "right_to_left_ports"),
        left_to_right_orientation=_bundle_int_tuple(row, "left_to_right_orientation"),
        right_to_left_orientation=_bundle_int_tuple(row, "right_to_left_orientation"),
        collar_kind=_bundle_string(row, "collar_kind"),
        interface_algebra=InterfaceAlgebraBinding(
            interface_algebra_id=_bundle_string(binding, "interface_algebra_id"),
            interface_algebra_sha256=_bundle_string(
                binding, "interface_algebra_sha256"
            ),
            left_interface_algebra_sha256=_bundle_string(
                binding, "left_interface_algebra_sha256"
            ),
            right_interface_algebra_sha256=_bundle_string(
                binding, "right_interface_algebra_sha256"
            ),
        ),
    )


def _boundary_from_bundle_row(row: Mapping[str, Any]) -> ExternalBoundaryBundle:
    if not isinstance(row, Mapping):
        raise TypeError("external boundary row must be an object")
    return ExternalBoundaryBundle(
        boundary_id=_bundle_string(row, "boundary_id"),
        carrier_id=_bundle_string(row, "carrier_id"),
        ports=_bundle_int_tuple(row, "ports"),
        boundary_condition=_bundle_string(row, "boundary_condition"),
        boundary_algebra_sha256=_bundle_string(row, "boundary_algebra_sha256"),
    )


def _triple_overlap_from_bundle_row(
    row: Mapping[str, Any],
) -> TripleOverlapBundle:
    if not isinstance(row, Mapping):
        raise TypeError("triple overlap row must be an object")
    carriers = _bundle_string_tuple(row, "oriented_carrier_ids")
    seams = _bundle_string_tuple(row, "oriented_seam_ids")
    if len(carriers) != 3 or len(seams) != 3:
        raise ValueError("triple overlap rows require exactly three carriers and seams")
    return TripleOverlapBundle(
        overlap_id=_bundle_string(row, "overlap_id"),
        oriented_carrier_ids=(carriers[0], carriers[1], carriers[2]),
        oriented_seam_ids=(seams[0], seams[1], seams[2]),
        restriction_algebra_id=_bundle_string(row, "restriction_algebra_id"),
        restriction_algebra_sha256=_bundle_string(
            row, "restriction_algebra_sha256"
        ),
    )


def _support_binding_from_bundle_value(
    value: Any,
) -> SupportTowerBinding | None:
    if value is None:
        return None
    if not isinstance(value, Mapping):
        raise TypeError("support_tower_binding must be an object or null")
    raw_rows = value["carrier_to_base_vertex"]
    if isinstance(raw_rows, (str, bytes)) or not isinstance(raw_rows, Sequence):
        raise TypeError("carrier_to_base_vertex must be an array")
    rows: list[tuple[str, int]] = []
    for row in raw_rows:
        if (
            isinstance(row, (str, bytes))
            or not isinstance(row, Sequence)
            or len(row) != 2
            or not isinstance(row[0], str)
            or type(row[1]) is not int
        ):
            raise TypeError(
                "carrier_to_base_vertex rows must be [string, exact integer]"
            )
        rows.append((row[0], row[1]))
    return SupportTowerBinding(
        geometry_family=_bundle_string(value, "geometry_family"),
        carrier_to_base_vertex=tuple(rows),
        orientation=_bundle_string(value, "orientation"),
        refinement_rule=_bundle_string(value, "refinement_rule"),
        source_incidence_sha256=_bundle_string(value, "source_incidence_sha256"),
    )


def _observer_from_bundle_row(row: Mapping[str, Any]) -> ObserverSupport:
    if not isinstance(row, Mapping):
        raise TypeError("observer support row must be an object")
    return ObserverSupport(
        observer_token=_bundle_string(row, "observer_token"),
        carrier_ids=frozenset(_bundle_string_tuple(row, "carrier_ids")),
        visible_seam_ids=frozenset(_bundle_string_tuple(row, "visible_seam_ids")),
        record_algebra_sha256=_bundle_string(row, "record_algebra_sha256"),
        checkpoint_cut_sha256=_bundle_string(row, "checkpoint_cut_sha256"),
    )


def _bundle_string(row: Mapping[str, Any], key: str) -> str:
    value = row[key]
    if not isinstance(value, str):
        raise TypeError(f"{key} must be a string")
    return value


def _bundle_int_tuple(row: Mapping[str, Any], key: str) -> tuple[int, ...]:
    value = row[key]
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise TypeError(f"{key} must be an integer array")
    if not all(type(item) is int for item in value):
        raise TypeError(f"{key} must contain exact integers")
    return tuple(value)


def _bundle_string_tuple(row: Mapping[str, Any], key: str) -> tuple[str, ...]:
    value = row[key]
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise TypeError(f"{key} must be a string array")
    if not all(isinstance(item, str) for item in value):
        raise TypeError(f"{key} must contain only strings")
    return tuple(value)


def _is_sha256(value: str) -> bool:
    return bool(isinstance(value, str) and _SHA256_RE.fullmatch(value))


if __name__ == "__main__":
    raise SystemExit(main())
