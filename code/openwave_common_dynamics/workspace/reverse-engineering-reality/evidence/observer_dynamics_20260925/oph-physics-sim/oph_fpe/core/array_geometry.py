"""Screen geometry routing for the large-array engines.

Three geometry families are routed here.

``icosahedral_tower`` is the production carrier family.  Every patch is one
outward-oriented spherical triangular cell of the committed geodesic
icosahedral refinement tower (:mod:`oph_fpe.core.icosahedral`), every patch
adjacency is the dual seam across one shared mesh edge, and the twelve
persistent ``q(v) = 6 - degree(v) = 1`` primal vertices are the global
defect-port sites.  At refinement level zero the complex is exactly the
regular icosahedron: 12 vertices, 30 edges, 20 faces.  Exact cell counts are
``20 * 4**level``; arbitrary counts are rejected rather than truncated or
padded.

Port-slot reconciliation for the production family: the twelve-slot local
port template of :mod:`oph_fpe.core.screen_ports` (``ports_per_patch: 12``)
is the hidden per-patch echosahedral presentation interface, not the seam
adjacency.  The carrier adjacency graph is exactly 3-regular at every level,
so exactly three of the twelve local slots route seams and the remaining nine
stay exposed/reserved local state.  The port map records the split per patch
(``node_degree``, ``unused_port_slots``,
``unrouted_exposed_or_reserved_port_slot_count``) and the geometry report
records the exact degree histogram; no slot is padded with a fabricated
neighbor and any degree above twelve fails closed through the port-map
overflow counter.

The ``nested_geodesic_icosahedral`` chart families supply the same tower as a
strict global support chart, typed separately from any carrier federation.
``fibonacci_sphere`` (and any other unlisted family) is the legacy KNN
control cellulation.  Identifying a legacy chart cell with a carrier requires
a realization theorem that this adapter does not provide; the legacy reports
expose that compatibility shortcut and never promote it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np

from oph_fpe.core.graph import fibonacci_sphere_points
from oph_fpe.core.icosahedral import (
    GeodesicIcosahedralTower,
    build_geodesic_icosahedral_tower,
    geodesic_icosahedral_patch_arrays,
    icosahedral_a5_equivariance_report,
    icosahedral_count_bracket,
    icosahedral_defect_port_report,
    resolve_icosahedral_level,
    supported_icosahedral_count,
)


ICOSAHEDRAL_ARRAY_FAMILIES = frozenset(
    {
        "subdivided_icosahedral_screen",
        "nested_geodesic_icosahedral",
        "geodesic_icosahedral_refinement",
        "nested_icosahedral_cap_net",
    }
)

ICOSAHEDRAL_TOWER_FAMILY = "icosahedral_tower"
ICOSAHEDRAL_TOWER_PATCH_BASIS = "cells"
ICOSAHEDRAL_TOWER_PORTS_PER_PATCH_TEMPLATE = 12
ICOSAHEDRAL_TOWER_SEAM_SLOTS_PER_PATCH = 3


@dataclass(frozen=True)
class ArrayScreenGeometry:
    points: np.ndarray
    edge_left: np.ndarray
    edge_right: np.ndarray
    report: dict[str, Any]

    @property
    def patch_count(self) -> int:
        return int(self.points.shape[0])

    @property
    def edge_count(self) -> int:
        return int(self.edge_left.size)


def icosahedral_patch_degree_histogram(
    level: int,
    *,
    patch_basis: str = "cells",
) -> dict[int, int]:
    """Exact adjacency-degree histogram of one tower level's patch graph.

    Cell-basis patches are dual-adjacent across shared mesh edges, so every
    outward-oriented spherical triangular cell has exactly three neighbors at
    every level: the histogram is ``{3: 20 * 4**level}``.  Vertex-basis
    patches follow the primal triangulation: the twelve base vertices keep
    degree five at every level and every refinement vertex has degree six, so
    the histogram is ``{5: 12}`` at level zero and
    ``{5: 12, 6: 10 * 4**level - 10}`` above it.  The histogram is computed
    from the built arrays, not asserted from the formulas.
    """

    points, left, right = geodesic_icosahedral_patch_arrays(
        level,
        patch_basis=patch_basis,
    )
    degree = np.bincount(
        np.concatenate(
            (np.asarray(left, dtype=np.int64), np.asarray(right, dtype=np.int64))
        ),
        minlength=int(points.shape[0]),
    )
    values, counts = np.unique(degree, return_counts=True)
    return {int(value): int(count) for value, count in zip(values, counts, strict=True)}


def resolve_icosahedral_tower_rung(graph_config: Mapping[str, Any] | None) -> dict[str, Any]:
    """Validate an ``icosahedral_tower`` graph section without building meshes.

    The production family accepts only exact tower rungs: an explicit
    ``refinement_level``, an exact cell ``patch_count``, or both in agreement.
    ``nominal_patch_count`` stays a campaign label.  All returned counts are
    closed-form (cells ``20 * 4**L``, seams ``30 * 4**L``, port-site vertices
    ``10 * 4**L + 2``), so configuration validation at any rung costs no mesh
    construction.
    """

    graph = dict(graph_config or {})
    family = str(graph.get("family", ""))
    if family != ICOSAHEDRAL_TOWER_FAMILY:
        raise ValueError(
            f"resolve_icosahedral_tower_rung requires family: {ICOSAHEDRAL_TOWER_FAMILY}, "
            f"got {family!r}"
        )
    patch_basis = str(graph.get("patch_basis", ICOSAHEDRAL_TOWER_PATCH_BASIS))
    if patch_basis != ICOSAHEDRAL_TOWER_PATCH_BASIS:
        raise ValueError(
            "the icosahedral_tower carrier family places patch state on the "
            "spherical triangular cells of the committed complex; vertex "
            "support charts stay available under the nested_geodesic_icosahedral "
            "chart families"
        )
    if "patch_count_policy" in graph and str(graph["patch_count_policy"]) != "exact":
        raise ValueError(
            "the icosahedral_tower carrier family supports only exact rungs; "
            "nearest/floor/ceil count policies belong to the chart families"
        )
    level: int
    if "refinement_level" in graph:
        level = int(graph["refinement_level"])
        expected = supported_icosahedral_count(level, ICOSAHEDRAL_TOWER_PATCH_BASIS)
        if "patch_count" in graph and int(graph["patch_count"]) != expected:
            raise ValueError(
                f"refinement_level={level} has exactly {expected} cells, not "
                f"patch_count={int(graph['patch_count'])}; use nominal_patch_count "
                "for a campaign rung label"
            )
    elif "patch_count" in graph:
        level, _ = resolve_icosahedral_level(
            int(graph["patch_count"]),
            patch_basis=ICOSAHEDRAL_TOWER_PATCH_BASIS,
            policy="exact",
        )
    else:
        raise ValueError(
            "icosahedral_tower requires refinement_level or an exact cell patch_count"
        )
    cell_count = supported_icosahedral_count(level, "cells")
    vertex_count = supported_icosahedral_count(level, "vertices")
    seam_count = 30 * 4**level
    nominal_raw = graph.get("nominal_patch_count")
    nominal_count = int(nominal_raw) if nominal_raw is not None else cell_count
    return {
        "schema": "oph.icosahedral_tower_rung.v1",
        "family": ICOSAHEDRAL_TOWER_FAMILY,
        "patch_basis": ICOSAHEDRAL_TOWER_PATCH_BASIS,
        "refinement_level": int(level),
        "patch_count": cell_count,
        "cell_count": cell_count,
        "seam_count": seam_count,
        "port_site_vertex_count": vertex_count,
        "euler_characteristic": vertex_count - seam_count + cell_count,
        "base_carrier_counts": {"vertices": 12, "edges": 30, "faces": 20},
        "expected_degree_histogram": {"3": cell_count},
        "ports_per_patch_template": ICOSAHEDRAL_TOWER_PORTS_PER_PATCH_TEMPLATE,
        "routed_seam_slots_per_patch": ICOSAHEDRAL_TOWER_SEAM_SLOTS_PER_PATCH,
        "exposed_reserved_slots_per_patch": (
            ICOSAHEDRAL_TOWER_PORTS_PER_PATCH_TEMPLATE
            - ICOSAHEDRAL_TOWER_SEAM_SLOTS_PER_PATCH
        ),
        "nominal_patch_count": nominal_count,
        "nominal_count_is_exact": nominal_count == cell_count,
        "ICOSAHEDRAL_TOWER_RUNG_RESOLUTION_RECEIPT": True,
    }


def icosahedral_tower_carrier_verification(
    points: np.ndarray,
    edge_left: np.ndarray,
    edge_right: np.ndarray,
    *,
    refinement_level: int,
) -> dict[str, Any]:
    """Fail-closed comparison of a supplied carrier against the committed tower.

    The reference is rebuilt from :func:`geodesic_icosahedral_patch_arrays` in
    the canonical cell basis and the supplied arrays must match it exactly:
    identical cell-center coordinates, identical seam endpoint arrays in
    canonical order, the exact 3-regular degree histogram, and the level-zero
    base complex equal to the regular icosahedron (12/30/20, Euler
    characteristic two).  A shuffled, rewired, dropped, or padded adjacency
    fails the receipt.
    """

    level = int(refinement_level)
    reference_points, reference_left, reference_right = geodesic_icosahedral_patch_arrays(
        level,
        patch_basis=ICOSAHEDRAL_TOWER_PATCH_BASIS,
    )
    supplied_points = np.asarray(points, dtype=float)
    supplied_left = np.asarray(edge_left, dtype=np.int64)
    supplied_right = np.asarray(edge_right, dtype=np.int64)
    expected_cells = supported_icosahedral_count(level, "cells")
    expected_seams = 30 * 4**level

    patch_count_exact = bool(supplied_points.shape == reference_points.shape)
    seam_count_exact = bool(
        supplied_left.shape == reference_left.shape
        and supplied_right.shape == reference_right.shape
    )
    endpoint_indices_in_range = bool(
        seam_count_exact
        and (
            supplied_left.size == 0
            or (
                int(min(supplied_left.min(), supplied_right.min())) >= 0
                and int(max(supplied_left.max(), supplied_right.max())) < expected_cells
            )
        )
    )
    points_match_committed_builder = bool(
        patch_count_exact and np.array_equal(supplied_points, reference_points)
    )
    seams_match_committed_builder = bool(
        seam_count_exact
        and endpoint_indices_in_range
        and np.array_equal(supplied_left, reference_left)
        and np.array_equal(supplied_right, reference_right)
    )
    if patch_count_exact and endpoint_indices_in_range:
        degree = np.bincount(
            np.concatenate((supplied_left, supplied_right)),
            minlength=expected_cells,
        )
        values, counts = np.unique(degree, return_counts=True)
        degree_histogram = {
            str(int(value)): int(count)
            for value, count in zip(values, counts, strict=True)
        }
    else:
        degree_histogram = {}
    degree_histogram_exact = degree_histogram == {"3": expected_cells}

    base = build_geodesic_icosahedral_tower(0).levels[0]
    base_carrier_is_regular_icosahedron = bool(
        base.vertex_count == 12
        and base.edge_count == 30
        and base.face_count == 20
        and base.euler_characteristic == 2
    )
    tower_receipt = build_geodesic_icosahedral_tower(level).receipt()[
        "GEODESIC_ICOSAHEDRAL_TOWER_RECEIPT"
    ]
    receipt = bool(
        patch_count_exact
        and seam_count_exact
        and endpoint_indices_in_range
        and points_match_committed_builder
        and seams_match_committed_builder
        and degree_histogram_exact
        and base_carrier_is_regular_icosahedron
        and tower_receipt
    )
    return {
        "schema": "oph.icosahedral_tower_carrier_verification.v1",
        "refinement_level": level,
        "expected_cell_count": expected_cells,
        "expected_seam_count": expected_seams,
        "patch_count_exact": patch_count_exact,
        "seam_count_exact": seam_count_exact,
        "endpoint_indices_in_range": endpoint_indices_in_range,
        "points_match_committed_builder": points_match_committed_builder,
        "seams_match_committed_builder_in_canonical_order": seams_match_committed_builder,
        "degree_histogram": degree_histogram,
        "degree_histogram_exact": degree_histogram_exact,
        "base_carrier_is_regular_icosahedron": base_carrier_is_regular_icosahedron,
        "tower_geometry_receipt": bool(tower_receipt),
        "TRUE_TOWER_CARRIER_VERIFICATION_RECEIPT": receipt,
    }


def array_screen_geometry_from_config(
    graph_config: Mapping[str, Any] | None,
    *,
    knn_builder: Any,
) -> ArrayScreenGeometry:
    """Instantiate the production tower carrier, a strict chart, or a control.

    ``family: icosahedral_tower`` builds the production carrier federation on
    the committed refinement tower.  The chart families build the same tower
    as a separately typed support chart.  Every other family is the legacy
    Fibonacci/KNN control; legacy engines still place one local state row on
    each returned chart point, and the legacy reports expose that
    compatibility shortcut without promoting it.
    """

    graph = dict(graph_config or {})
    family = str(graph.get("family", "fibonacci_sphere"))
    if family == ICOSAHEDRAL_TOWER_FAMILY:
        return _icosahedral_tower_screen_geometry(graph)
    requested = int(graph.get("patch_count", graph.get("nodes", 65_536)))
    if family not in ICOSAHEDRAL_ARRAY_FAMILIES:
        neighbors = int(graph.get("neighbors", graph.get("degree", 8)))
        points = fibonacci_sphere_points(requested)
        left, right = knn_builder(points, neighbors)
        report = {
            "schema": "oph.array_screen_geometry.v1",
            "mesh_family": family,
            "geometry_family": "legacy_fibonacci_knn_control",
            "nominal_patch_count": requested,
            "actual_patch_count": int(points.shape[0]),
            "edge_count": int(left.size),
            "patch_basis": "points",
            "TRUE_ICOSAHEDRAL_REFINEMENT_TOWER_RECEIPT": False,
            "NESTED_ICOSAHEDRAL_LINEAGE_RECEIPT": False,
            "A5_EQUIVARIANT_REFINEMENT_RECEIPT": False,
            "PAPER_MULTIRESOLUTION_REGULATOR_CERTIFICATE": False,
            "CARRIER_TO_SUPPORT_REALIZATION_RECEIPT": False,
            "CARRIER_SUPPORT_CONFLATION_PRESENT_IN_LEGACY_ENGINE": True,
            "claim_boundary": (
                "Legacy Fibonacci/KNN global cellulation retained only as a control. "
                "The legacy engine may attach local state rows to these points, but "
                "that compatibility shortcut is not a carrier federation or an "
                "earned support chart."
            ),
        }
        return ArrayScreenGeometry(points, left, right, report)

    patch_basis = str(graph.get("patch_basis", graph.get("node_basis", "cells")))
    if patch_basis not in {"cells", "vertices"}:
        raise ValueError("icosahedral array patch_basis must be 'cells' or 'vertices'")
    nominal_raw = graph.get("nominal_patch_count")
    nominal_count = int(nominal_raw) if nominal_raw is not None else None
    if "refinement_level" in graph:
        level = int(graph["refinement_level"])
        actual_count = supported_icosahedral_count(level, patch_basis)
        if "patch_count" in graph and int(graph["patch_count"]) != actual_count:
            raise ValueError(
                f"refinement_level={level} on {patch_basis} has {actual_count} "
                f"patches, not patch_count={int(graph['patch_count'])}; use "
                "nominal_patch_count for the campaign rung label"
            )
    else:
        policy = str(graph.get("patch_count_policy", "exact"))
        level, bracket = resolve_icosahedral_level(
            requested,
            patch_basis=patch_basis,
            policy=policy,
        )
        actual_count = supported_icosahedral_count(level, patch_basis)
        if not bracket["exact_supported_count"] and nominal_count is None:
            nominal_count = requested
    if nominal_count is None and "patch_count" not in graph:
        nominal_count = actual_count
    points, left, right = geodesic_icosahedral_patch_arrays(
        level,
        patch_basis=patch_basis,
    )
    tower = build_geodesic_icosahedral_tower(level)
    tower_receipt = tower.receipt()
    a5_receipt = icosahedral_a5_equivariance_report(level)
    defect_ports = icosahedral_defect_port_report(level)
    levels = _tower_level_rows(tower, patch_basis)
    expectations = _tower_expectation_rows(tower)
    report = {
        "schema": "oph.array_screen_geometry.v1",
        "mesh_family": "nested_geodesic_icosahedral",
        "geometry_family": "nested_geodesic_icosahedral",
        "patch_basis": patch_basis,
        "refinement_level": level,
        "nominal_patch_count": nominal_count,
        "actual_patch_count": int(points.shape[0]),
        "edge_count": int(left.size),
        "nominal_count_mapping": (
            icosahedral_count_bracket(nominal_count, patch_basis)
            if nominal_count is not None
            else None
        ),
        "levels": levels,
        "conditional_expectations": expectations,
        "nested_lineage_receipt": tower_receipt[
            "GEODESIC_ICOSAHEDRAL_TOWER_RECEIPT"
        ],
        "conditional_expectations_receipt": bool(
            expectations and all(row["state_preserving"] for row in expectations)
        ),
        "TRUE_ICOSAHEDRAL_REFINEMENT_TOWER_RECEIPT": tower_receipt[
            "GEODESIC_ICOSAHEDRAL_TOWER_RECEIPT"
        ],
        "NESTED_ICOSAHEDRAL_LINEAGE_RECEIPT": tower_receipt[
            "GEODESIC_ICOSAHEDRAL_TOWER_RECEIPT"
        ],
        "A5_EQUIVARIANT_REFINEMENT_RECEIPT": a5_receipt[
            "A5_EQUIVARIANT_REFINEMENT_RECEIPT"
        ],
        "REFERENCE_ICOSAHEDRAL_A5_GEOMETRY_RECEIPT": a5_receipt[
            "REFERENCE_ICOSAHEDRAL_A5_GEOMETRY_RECEIPT"
        ],
        "TWELVE_PERSISTENT_COMBINATORIAL_DEFECT_PORTS_RECEIPT": defect_ports[
            "TWELVE_PERSISTENT_COMBINATORIAL_DEFECT_PORTS_RECEIPT"
        ],
        "PAPER_MULTIRESOLUTION_REGULATOR_CERTIFICATE": False,
        "CARRIER_TO_SUPPORT_REALIZATION_RECEIPT": False,
        "CARRIER_SUPPORT_CONFLATION_PRESENT_IN_LEGACY_ENGINE": True,
        "tower": tower_receipt,
        "a5_action": a5_receipt,
        "global_screen_sieve_defect_ports": defect_ports,
        "full_multiresolution_blockers": [
            "finite_detail_algebras_not_instantiated",
            "faithful_prime_cap_states_not_instantiated",
            "local_presentation_circuits_not_instantiated",
            "noncommutative_state_preserving_expectations_not_instantiated",
            "mixed_gns_reference_tower_not_instantiated",
        ],
        "claim_boundary": (
            "Reference global nested icosahedral support chart plus exact "
            "integer-permutation A5 equivariance and the commutative cell-area "
            "expectation scaffold. It is separate from the exact-count microscopic "
            "carrier federation: no global cell is identified with one carrier by "
            "this report. The carrier-to-support realization and full "
            "noncommutative multiresolution/BW certificates remain false."
        ),
    }
    return ArrayScreenGeometry(points, left, right, report)


def _icosahedral_tower_screen_geometry(graph: dict[str, Any]) -> ArrayScreenGeometry:
    """Build the production carrier federation on the committed tower."""

    rung = resolve_icosahedral_tower_rung(graph)
    level = int(rung["refinement_level"])
    points, left, right = geodesic_icosahedral_patch_arrays(
        level,
        patch_basis=ICOSAHEDRAL_TOWER_PATCH_BASIS,
    )
    tower = build_geodesic_icosahedral_tower(level)
    tower_receipt = tower.receipt()
    a5_receipt = icosahedral_a5_equivariance_report(level)
    defect_ports = icosahedral_defect_port_report(level)
    verification = icosahedral_tower_carrier_verification(
        points,
        left,
        right,
        refinement_level=level,
    )
    levels = _tower_level_rows(tower, ICOSAHEDRAL_TOWER_PATCH_BASIS)
    expectations = _tower_expectation_rows(tower)
    degree_histogram = verification["degree_histogram"]
    max_degree = max((int(key) for key in degree_histogram), default=0)
    port_slot_receipt = bool(
        verification["degree_histogram_exact"]
        and max_degree <= ICOSAHEDRAL_TOWER_PORTS_PER_PATCH_TEMPLATE
        and 2 * int(left.size) == 3 * int(points.shape[0])
    )
    true_tower_receipt = bool(
        tower_receipt["GEODESIC_ICOSAHEDRAL_TOWER_RECEIPT"]
        and verification["TRUE_TOWER_CARRIER_VERIFICATION_RECEIPT"]
    )
    report = {
        "schema": "oph.array_screen_geometry.v1",
        "mesh_family": ICOSAHEDRAL_TOWER_FAMILY,
        "geometry_family": ICOSAHEDRAL_TOWER_FAMILY,
        "patch_basis": ICOSAHEDRAL_TOWER_PATCH_BASIS,
        "refinement_level": level,
        "nominal_patch_count": rung["nominal_patch_count"],
        "actual_patch_count": int(points.shape[0]),
        "edge_count": int(left.size),
        "nominal_count_mapping": icosahedral_count_bracket(
            int(rung["nominal_patch_count"]),
            ICOSAHEDRAL_TOWER_PATCH_BASIS,
        ),
        "rung_resolution": rung,
        "carrier_semantics": {
            "patch": "outward_oriented_spherical_triangular_cell_of_the_committed_tower",
            "patch_adjacency": "dual_seam_across_one_shared_mesh_edge",
            "seam_count_equals_primal_mesh_edge_count": int(left.size) == 30 * 4**level,
            "global_defect_port_sites": "twelve_persistent_q_equals_1_primal_vertices",
            "base_complex_at_level_0": {"vertices": 12, "edges": 30, "faces": 20},
            "a5_action": "exact_integer_cell_permutations_certified_per_level",
        },
        "carrier_degree_histogram": degree_histogram,
        "carrier_verification": verification,
        "port_slot_reconciliation": {
            "local_port_template": "regular_icosahedron_twelve_vertex_directions",
            "ports_per_patch_template": ICOSAHEDRAL_TOWER_PORTS_PER_PATCH_TEMPLATE,
            "template_role": (
                "hidden per-patch echosahedral presentation interface; the seam "
                "adjacency is the 3-regular cell-dual graph, not a degree-12 KNN graph"
            ),
            "carrier_adjacency_degree_histogram": degree_histogram,
            "routed_seam_slots_per_patch": ICOSAHEDRAL_TOWER_SEAM_SLOTS_PER_PATCH,
            "exposed_reserved_slots_per_patch": (
                ICOSAHEDRAL_TOWER_PORTS_PER_PATCH_TEMPLATE
                - ICOSAHEDRAL_TOWER_SEAM_SLOTS_PER_PATCH
            ),
            "unused_slot_marking": (
                "the screen-port map records node_degree, unused_port_slots, and "
                "unrouted_exposed_or_reserved_port_slot_count per federation"
            ),
            "adjacency_padding_free": verification["degree_histogram_exact"],
            "overflow_fail_closed": (
                "any patch degree above the twelve-slot template raises the "
                "screen-port overflow counter and fails the port-state receipts"
            ),
        },
        "levels": levels,
        "conditional_expectations": expectations,
        "nested_lineage_receipt": tower_receipt["GEODESIC_ICOSAHEDRAL_TOWER_RECEIPT"],
        "conditional_expectations_receipt": bool(
            expectations and all(row["state_preserving"] for row in expectations)
        ),
        "TRUE_ICOSAHEDRAL_REFINEMENT_TOWER_RECEIPT": true_tower_receipt,
        "NESTED_ICOSAHEDRAL_LINEAGE_RECEIPT": tower_receipt[
            "GEODESIC_ICOSAHEDRAL_TOWER_RECEIPT"
        ],
        "A5_EQUIVARIANT_REFINEMENT_RECEIPT": a5_receipt[
            "A5_EQUIVARIANT_REFINEMENT_RECEIPT"
        ],
        "REFERENCE_ICOSAHEDRAL_A5_GEOMETRY_RECEIPT": a5_receipt[
            "REFERENCE_ICOSAHEDRAL_A5_GEOMETRY_RECEIPT"
        ],
        "TWELVE_PERSISTENT_COMBINATORIAL_DEFECT_PORTS_RECEIPT": defect_ports[
            "TWELVE_PERSISTENT_COMBINATORIAL_DEFECT_PORTS_RECEIPT"
        ],
        "BASE_CARRIER_IS_REGULAR_ICOSAHEDRON_RECEIPT": verification[
            "base_carrier_is_regular_icosahedron"
        ],
        "CARRIER_ADJACENCY_DEGREE_HISTOGRAM_RECEIPT": verification[
            "degree_histogram_exact"
        ],
        "PORT_SLOT_DEGREE_RECONCILIATION_RECEIPT": port_slot_receipt,
        "TRUE_TOWER_CARRIER_VERIFICATION_RECEIPT": verification[
            "TRUE_TOWER_CARRIER_VERIFICATION_RECEIPT"
        ],
        "PAPER_MULTIRESOLUTION_REGULATOR_CERTIFICATE": False,
        "CARRIER_TO_SUPPORT_REALIZATION_RECEIPT": False,
        "CARRIER_SUPPORT_CONFLATION_PRESENT_IN_LEGACY_ENGINE": False,
        "PHYSICAL_A5_PORT_EMERGENCE_RECEIPT": False,
        "tower": tower_receipt,
        "a5_action": a5_receipt,
        "global_screen_sieve_defect_ports": defect_ports,
        "full_multiresolution_blockers": [
            "finite_detail_algebras_not_instantiated",
            "faithful_prime_cap_states_not_instantiated",
            "local_presentation_circuits_not_instantiated",
            "noncommutative_state_preserving_expectations_not_instantiated",
            "mixed_gns_reference_tower_not_instantiated",
        ],
        "claim_boundary": (
            "Production carrier family: every patch is one outward-oriented "
            "spherical triangular cell of the committed geodesic icosahedral "
            "refinement tower, every patch adjacency is the dual seam across one "
            "shared mesh edge, and the twelve persistent q=1 primal vertices are "
            "the global defect-port sites. At level 0 the complex is exactly the "
            "regular icosahedron (12 vertices, 30 edges, 20 faces). The "
            "twelve-slot local port template is per-patch presentation data: "
            "exactly three slots route seams and nine stay exposed/reserved, "
            "with the split recorded in the port map. The identification of "
            "engine state rows with tower cells is this family's declared "
            "architecture, typed as one object rather than a chart silently "
            "standing in for a carrier federation; it is not a theorem that "
            "physical source dynamics selected or realized this federation, and "
            "the physical emergence, carrier-to-support realization, and "
            "noncommutative multiresolution/BW certificates stay false."
        ),
    }
    return ArrayScreenGeometry(points, left, right, report)


def _tower_level_rows(
    tower: GeodesicIcosahedralTower,
    patch_basis: str,
) -> list[dict[str, Any]]:
    """Emit per-level lineage rows shared by the tower-backed families."""

    rows: list[dict[str, Any]] = []
    for item in tower.levels:
        rows.append(
            {
                "level_id": f"ico-L{item.level}",
                "refinement_level": item.level,
                "patch_basis": patch_basis,
                "patch_count": (
                    item.face_count if patch_basis == "cells" else item.vertex_count
                ),
                "vertex_count": item.vertex_count,
                "cell_count": item.face_count,
                "geometry_hash": item.geometry_hash,
                **(
                    {}
                    if item.level == 0
                    else {
                        "parent_level_id": f"ico-L{item.level - 1}",
                        "lineage_hash": tower.cell_refinements[
                            item.level - 1
                        ].map_hash,
                    }
                ),
            }
        )
    return rows


def _tower_expectation_rows(tower: GeodesicIcosahedralTower) -> list[dict[str, Any]]:
    """Emit conditional-expectation rows shared by the tower-backed families."""

    rows: list[dict[str, Any]] = []
    for mapping in tower.cell_refinements:
        rows.append(
            {
                "fine_level_id": f"ico-L{mapping.fine_level}",
                "coarse_level_id": f"ico-L{mapping.coarse_level}",
                "operator_hash": mapping.map_hash,
                "scope": "commutative_cell_geometry_scaffold",
                "unital": True,
                "positive": True,
                "state_preserving": mapping.state_preservation_residual <= 5.0e-14,
                "cap_isotony_compatible": True,
                "noncommutative_prime_cap_expectation": False,
            }
        )
    return rows
