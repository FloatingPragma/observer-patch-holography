"""Nested geodesic icosahedral screen geometry.

This module implements the *geometric* part of the multiresolution regulator
used in the OPH papers.  It deliberately does not claim to construct the
finite detail algebras, presentation circuits, or faithful noncommutative
states required by the full multiresolution regulator certificate.

The canonical patch basis is the set of outward-oriented spherical triangular
cells.  A refinement replaces every cell by four geodesic children.  On the
commutative cell-value scaffold, the refinement embedding is pullback along
``child -> parent`` and the conditional expectation is spherical-area-weighted
averaging over the four children.  These maps are unital, positive, left
inverse, and preserve the normalized spherical-area state.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from functools import lru_cache
from typing import Iterable, Literal

import networkx as nx
import numpy as np
from scipy.spatial import cKDTree


PatchBasis = Literal["cells", "vertices"]


def _readonly(array: np.ndarray) -> np.ndarray:
    result = np.ascontiguousarray(array)
    result.setflags(write=False)
    return result


@dataclass(frozen=True)
class CellRefinementMap:
    """Adjacent-level lineage and maps for the commutative cell scaffold."""

    coarse_level: int
    fine_level: int
    child_to_parent: np.ndarray
    children_by_parent: tuple[tuple[int, int, int, int], ...]
    conditional_expectation_weights: np.ndarray
    coarse_reference_weights: np.ndarray
    fine_reference_weights: np.ndarray
    normalization_residual: float
    state_preservation_residual: float
    map_hash: str

    def embed(self, coarse_values: np.ndarray) -> np.ndarray:
        """Embed a coarse cell observable as a child-wise constant observable."""

        values = np.asarray(coarse_values)
        if values.ndim == 0 or values.shape[0] != len(self.children_by_parent):
            raise ValueError(
                "coarse_values first dimension must equal the coarse cell count "
                f"({len(self.children_by_parent)})"
            )
        return values[self.child_to_parent]

    def conditional_expectation(self, fine_values: np.ndarray) -> np.ndarray:
        """Area-average a fine cell observable over each coarse parent."""

        values = np.asarray(fine_values)
        if values.ndim == 0 or values.shape[0] != self.child_to_parent.size:
            raise ValueError(
                "fine_values first dimension must equal the fine cell count "
                f"({self.child_to_parent.size})"
            )
        output = np.zeros((len(self.children_by_parent),) + values.shape[1:], dtype=np.result_type(values, float))
        for parent, children in enumerate(self.children_by_parent):
            child_ids = np.asarray(children, dtype=np.int64)
            weights = self.conditional_expectation_weights[child_ids]
            reshape = (weights.size,) + (1,) * (values.ndim - 1)
            output[parent] = np.sum(values[child_ids] * weights.reshape(reshape), axis=0)
        return output

    def receipt(self) -> dict:
        return {
            "schema": "oph.icosahedral_cell_refinement.v1",
            "scope": "commutative_cell_geometry_scaffold",
            "coarse_level": self.coarse_level,
            "fine_level": self.fine_level,
            "coarse_cell_count": len(self.children_by_parent),
            "fine_cell_count": int(self.child_to_parent.size),
            "children_per_parent": 4,
            "embedding": "childwise_constant_pullback",
            "conditional_expectation": "spherical_area_weighted_child_average",
            "reference_state": "normalized_spherical_area",
            "normalization_residual": self.normalization_residual,
            "state_preservation_residual": self.state_preservation_residual,
            "map_hash": self.map_hash,
            "GEOMETRIC_CELL_REFINEMENT_RECEIPT": bool(
                self.normalization_residual <= 5.0e-14
                and self.state_preservation_residual <= 5.0e-14
            ),
            "PAPER_MULTIRESOLUTION_REGULATOR_CERTIFICATE": False,
            "certificate_blockers": [
                "finite_detail_algebras_not_instantiated",
                "faithful_detail_states_not_instantiated",
                "local_presentation_circuits_not_instantiated",
                "noncommutative_state_preserving_expectations_not_instantiated",
            ],
        }


@dataclass(frozen=True)
class GeodesicIcosahedralLevel:
    """One deterministic level of a dyadically refined icosahedral sphere."""

    level: int
    frequency: int
    vertices: np.ndarray
    faces: np.ndarray
    edges: np.ndarray
    spherical_face_areas: np.ndarray
    vertex_birth_levels: np.ndarray
    vertex_parent_support: tuple[tuple[tuple[int, float], ...], ...]
    parent_face: np.ndarray | None
    children_by_parent_face: tuple[tuple[int, int, int, int], ...] | None
    geometry_hash: str

    @property
    def vertex_count(self) -> int:
        return int(self.vertices.shape[0])

    @property
    def edge_count(self) -> int:
        return int(self.edges.shape[0])

    @property
    def face_count(self) -> int:
        return int(self.faces.shape[0])

    @property
    def euler_characteristic(self) -> int:
        return self.vertex_count - self.edge_count + self.face_count

    def receipt(self) -> dict:
        radii = np.linalg.norm(self.vertices, axis=1)
        total_area = float(np.sum(self.spherical_face_areas))
        return {
            "schema": "oph.geodesic_icosahedral_level.v1",
            "geometry_family": "nested_geodesic_icosahedral",
            "level": self.level,
            "frequency": self.frequency,
            "vertex_count": self.vertex_count,
            "edge_count": self.edge_count,
            "face_count": self.face_count,
            "expected_counts": {
                "vertices": supported_icosahedral_count(self.level, "vertices"),
                "edges": 30 * self.frequency * self.frequency,
                "cells": supported_icosahedral_count(self.level, "cells"),
            },
            "euler_characteristic": self.euler_characteristic,
            "maximum_unit_sphere_residual": float(np.max(np.abs(radii - 1.0))),
            "total_spherical_area": total_area,
            "four_pi_area_residual": abs(total_area - 4.0 * math.pi),
            "outward_oriented": bool(_faces_are_outward(self.vertices, self.faces)),
            "geometry_hash": self.geometry_hash,
            "GEODESIC_ICOSAHEDRAL_GEOMETRY_RECEIPT": bool(
                self.euler_characteristic == 2
                and np.max(np.abs(radii - 1.0)) <= 5.0e-14
                and abs(total_area - 4.0 * math.pi) <= 5.0e-12
                and _faces_are_outward(self.vertices, self.faces)
            ),
        }


@dataclass(frozen=True)
class GeodesicIcosahedralTower:
    """A finite prefix of the nested icosahedral refinement tower."""

    levels: tuple[GeodesicIcosahedralLevel, ...]
    cell_refinements: tuple[CellRefinementMap, ...]

    @property
    def max_level(self) -> int:
        return len(self.levels) - 1

    def embed_cells(
        self,
        coarse_values: np.ndarray,
        *,
        coarse_level: int,
        fine_level: int,
    ) -> np.ndarray:
        """Apply the refinement embeddings between any two stored levels."""

        self._validate_level_interval(coarse_level, fine_level)
        values = np.asarray(coarse_values)
        for level in range(coarse_level, fine_level):
            values = self.cell_refinements[level].embed(values)
        return values

    def conditional_expectation_cells(
        self,
        fine_values: np.ndarray,
        *,
        fine_level: int,
        coarse_level: int,
    ) -> np.ndarray:
        """Coarse-grain cell observables between any two stored levels."""

        self._validate_level_interval(coarse_level, fine_level)
        values = np.asarray(fine_values)
        for level in range(fine_level - 1, coarse_level - 1, -1):
            values = self.cell_refinements[level].conditional_expectation(values)
        return values

    def cell_ancestor_ids(self, *, fine_level: int, coarse_level: int) -> np.ndarray:
        """Return each fine cell's ancestor identifier at ``coarse_level``."""

        self._validate_level_interval(coarse_level, fine_level)
        ancestors = np.arange(self.levels[fine_level].face_count, dtype=np.int64)
        for level in range(fine_level - 1, coarse_level - 1, -1):
            ancestors = self.cell_refinements[level].child_to_parent[ancestors]
        return ancestors

    def _validate_level_interval(self, coarse_level: int, fine_level: int) -> None:
        if not (0 <= coarse_level <= fine_level <= self.max_level):
            raise ValueError(
                f"require 0 <= coarse_level <= fine_level <= {self.max_level}"
            )

    def receipt(self) -> dict:
        geometry_receipts = [level.receipt() for level in self.levels]
        refinement_receipts = [mapping.receipt() for mapping in self.cell_refinements]
        return {
            "schema": "oph.geodesic_icosahedral_tower.v1",
            "geometry_family": "nested_geodesic_icosahedral",
            "levels": geometry_receipts,
            "cell_refinements": refinement_receipts,
            "embedding_composition": "adjacent_childwise_pullbacks_composed_by_tower",
            "conditional_expectation_composition": (
                "adjacent_spherical_area_expectations_composed_by_tower"
            ),
            "GEODESIC_ICOSAHEDRAL_TOWER_RECEIPT": bool(
                all(item["GEODESIC_ICOSAHEDRAL_GEOMETRY_RECEIPT"] for item in geometry_receipts)
                and all(item["GEOMETRIC_CELL_REFINEMENT_RECEIPT"] for item in refinement_receipts)
            ),
            "PAPER_MULTIRESOLUTION_REGULATOR_CERTIFICATE": False,
            "claim_boundary": (
                "This receipt certifies the nested geodesic mesh, cell lineage, and the "
                "state-preserving commutative spherical-area scaffold only. It does not "
                "instantiate the papers' detail algebras, faithful detail states, "
                "presentation circuits, or noncommutative conditional expectations."
            ),
        }


class UnsupportedIcosahedralPatchCount(ValueError):
    """Raised when an arbitrary count is presented as an exact geodesic mesh."""


def supported_icosahedral_count(level: int, basis: PatchBasis = "cells") -> int:
    """Return the exact number of cells or vertices at dyadic level ``level``."""

    if int(level) != level or level < 0:
        raise ValueError("refinement level must be a nonnegative integer")
    frequency_squared = 4 ** int(level)
    if basis == "cells":
        return 20 * frequency_squared
    if basis == "vertices":
        return 10 * frequency_squared + 2
    raise ValueError(f"unknown icosahedral patch basis: {basis!r}")


def icosahedral_count_bracket(requested_count: int, basis: PatchBasis = "cells") -> dict:
    """Describe the supported counts surrounding an arbitrary requested count."""

    requested = int(requested_count)
    if requested < 1:
        raise ValueError("requested_count must be positive")
    supported: list[tuple[int, int]] = []
    level = 0
    while True:
        count = supported_icosahedral_count(level, basis)
        supported.append((level, count))
        if count >= requested:
            break
        level += 1
    lower_items = [item for item in supported if item[1] <= requested]
    lower = lower_items[-1] if lower_items else None
    upper_items = [item for item in supported if item[1] >= requested]
    upper = upper_items[0]
    candidates = [item for item in (lower, upper) if item is not None]
    nearest = min(candidates, key=lambda item: (abs(item[1] - requested), item[1]))
    exact = nearest[1] == requested
    return {
        "basis": basis,
        "requested_count": requested,
        "exact_supported_count": exact,
        "lower": None if lower is None else {"level": lower[0], "count": lower[1]},
        "upper": {"level": upper[0], "count": upper[1]},
        "nearest": {
            "level": nearest[0],
            "count": nearest[1],
            "absolute_difference": abs(nearest[1] - requested),
            "relative_difference": abs(nearest[1] - requested) / requested,
        },
        "nominal_label_only_if_not_exact": not exact,
    }


def nominal_campaign_rung_mapping(
    requested_counts: Iterable[int] = (4096, 16384, 65536, 262144),
) -> dict:
    """Report honest cell/vertex mesh brackets for nominal campaign rungs."""

    rungs = []
    for requested in requested_counts:
        cell_bracket = icosahedral_count_bracket(int(requested), "cells")
        vertex_bracket = icosahedral_count_bracket(int(requested), "vertices")
        rungs.append(
            {
                "nominal_patch_count": int(requested),
                "cell_basis": cell_bracket,
                "vertex_basis": vertex_bracket,
                "recommended_geometry": {
                    "patch_basis": "cells",
                    "refinement_level": cell_bracket["nearest"]["level"],
                    "actual_patch_count": cell_bracket["nearest"]["count"],
                    "reason": (
                        "The papers count canonical geometric screen cells; the nominal "
                        "campaign rung remains a budget label, not an exact mesh count."
                    ),
                },
            }
        )
    return {
        "schema": "oph.icosahedral_nominal_rung_mapping.v1",
        "geometry_family": "nested_geodesic_icosahedral",
        "refinement_rule": "each_outward_spherical_triangle_has_four_children",
        "rungs": rungs,
        "NO_ARBITRARY_COUNT_TRUNCATION_RECEIPT": True,
    }


def icosahedral_a5_port_permutations() -> tuple[tuple[int, ...], ...]:
    """Return the 60 proper-rotation permutations of the base twelve ports.

    Each row maps an old port index to its rotated port index.  The immutable
    tuple is the canonical reference action used by the local echosahedral
    patch template; it is not a claim that physical source dynamics selected
    or coherently trivialized this action across a patch federation.
    """

    return _base_rotation_group()[1]


def icosahedral_a5_equivariance_report(max_level: int) -> dict:
    """Certify the exact discrete rotational action and its refinement lift.

    The action is first recovered as the 60 orientation-preserving rotations
    of the regular base icosahedron.  All group operations below are then
    checked on integer permutations, not on approximate character fits.  Two
    generators (of orders three and five) are lifted to every stored mesh
    level and checked against vertex incidence, face incidence, parent-face
    lineage, spherical-area weights, embeddings, and conditional
    expectations.

    This is an exact *regulator/reference-geometry* receipt.  It does not say
    that target-free source dynamics selected the icosahedron; that stronger
    physical producer claim has its own deliberately false flag below.
    """

    tower = build_geodesic_icosahedral_tower(max_level)
    rotations, base_permutations = _base_rotation_group()
    identity = tuple(range(tower.levels[0].vertex_count))
    permutation_set = set(base_permutations)
    group_closed = all(
        _compose_permutations(left, right) in permutation_set
        for left in base_permutations
        for right in base_permutations
    )
    inverse_closed = all(
        _inverse_permutation(permutation) in permutation_set
        for permutation in base_permutations
    )
    generator_indices = _find_a5_generator_pair(base_permutations)
    generator_rows: list[dict] = []
    level_rows: list[dict] = []
    prior_face_permutations: dict[str, np.ndarray] = {}
    all_level_checks = True

    for generator_name, generator_index in zip(
        ("order_3_generator", "order_5_generator"),
        generator_indices,
        strict=True,
    ):
        permutation = base_permutations[generator_index]
        generator_rows.append(
            {
                "generator_id": generator_name,
                "base_vertex_permutation": list(permutation),
                "order": _permutation_order(permutation),
                "permutation_hash": _hash_arrays(
                    np.asarray(permutation, dtype=np.int64)
                ),
            }
        )

    for level_index, mesh in enumerate(tower.levels):
        face_lookup = {
            tuple(sorted(int(value) for value in face)): face_id
            for face_id, face in enumerate(mesh.faces)
        }
        generator_checks: list[dict] = []
        current_face_permutations: dict[str, np.ndarray] = {}
        for generator_name, generator_index in zip(
            ("order_3_generator", "order_5_generator"),
            generator_indices,
            strict=True,
        ):
            rotation = rotations[generator_index]
            vertex_permutation, coordinate_residual = _coordinate_permutation(
                mesh.vertices, rotation
            )
            mapped_faces = np.asarray(
                [
                    face_lookup[
                        tuple(
                            sorted(
                                int(vertex_permutation[int(vertex)])
                                for vertex in face
                            )
                        )
                    ]
                    for face in mesh.faces
                ],
                dtype=np.int64,
            )
            current_face_permutations[generator_name] = mapped_faces
            vertex_edges = {
                (min(int(left), int(right)), max(int(left), int(right)))
                for left, right in mesh.edges
            }
            incidence_preserved = all(
                (
                    min(
                        int(vertex_permutation[int(left)]),
                        int(vertex_permutation[int(right)]),
                    ),
                    max(
                        int(vertex_permutation[int(left)]),
                        int(vertex_permutation[int(right)]),
                    ),
                )
                in vertex_edges
                for left, right in mesh.edges
            )
            face_area_residual = float(
                np.max(
                    np.abs(
                        mesh.spherical_face_areas[mapped_faces]
                        - mesh.spherical_face_areas
                    )
                )
            )
            parent_lineage_preserved = True
            expectation_weight_residual = 0.0
            if level_index > 0:
                parent_face = mesh.parent_face
                assert parent_face is not None
                coarse_faces = prior_face_permutations[generator_name]
                parent_lineage_preserved = bool(
                    np.array_equal(
                        parent_face[mapped_faces],
                        coarse_faces[parent_face],
                    )
                )
                refinement = tower.cell_refinements[level_index - 1]
                expectation_weight_residual = float(
                    np.max(
                        np.abs(
                            refinement.conditional_expectation_weights[mapped_faces]
                            - refinement.conditional_expectation_weights
                        )
                    )
                )
            passed = bool(
                coordinate_residual <= 5.0e-12
                and incidence_preserved
                and face_area_residual <= 5.0e-12
                and parent_lineage_preserved
                and expectation_weight_residual <= 5.0e-12
            )
            all_level_checks = all_level_checks and passed
            generator_checks.append(
                {
                    "generator_id": generator_name,
                    "vertex_permutation_hash": _hash_arrays(vertex_permutation),
                    "face_permutation_hash": _hash_arrays(mapped_faces),
                    "maximum_coordinate_residual": coordinate_residual,
                    "vertex_edge_incidence_preserved": incidence_preserved,
                    "maximum_face_area_residual": face_area_residual,
                    "parent_face_lineage_preserved": parent_lineage_preserved,
                    "conditional_expectation_weight_residual": expectation_weight_residual,
                    "passed": passed,
                }
            )
        prior_face_permutations = current_face_permutations
        level_rows.append(
            {
                "level": level_index,
                "vertex_count": mesh.vertex_count,
                "cell_count": mesh.face_count,
                "generator_checks": generator_checks,
                "passed": all(row["passed"] for row in generator_checks),
            }
        )

    base_group_passed = bool(
        len(base_permutations) == 60
        and identity in permutation_set
        and group_closed
        and inverse_closed
        and {_permutation_order(base_permutations[index]) for index in generator_indices}
        == {3, 5}
        and len(
            _generated_permutation_subgroup(
                [base_permutations[index] for index in generator_indices]
            )
        )
        == 60
    )
    receipt = bool(base_group_passed and all_level_checks)
    return {
        "schema": "oph.icosahedral_a5_equivariance.v1",
        "geometry_family": "nested_geodesic_icosahedral",
        "group_description": (
            "orientation-preserving regular-icosahedron rotation group, "
            "classified analytically as A5"
        ),
        "base_rotation_count": len(base_permutations),
        "integer_permutation_group_closed": group_closed,
        "integer_permutation_inverses_present": inverse_closed,
        "faithful_base_vertex_action": len(set(base_permutations)) == 60,
        "generator_rows": generator_rows,
        "level_rows": level_rows,
        "A5_ROTATION_GROUP_ORDER_60_RECEIPT": base_group_passed,
        "A5_EQUIVARIANT_REFINEMENT_RECEIPT": receipt,
        "REFERENCE_ICOSAHEDRAL_A5_GEOMETRY_RECEIPT": receipt,
        "PHYSICAL_A5_PORT_EMERGENCE_RECEIPT": False,
        "physical_emergence_blockers": [
            "calibrated_curvature_KL_risk_not_emitted_by_geometry_constructor",
            "ground_state_complete_settlement_not_emitted_by_geometry_constructor",
            "atomic_defect_projection_not_emitted_by_geometry_constructor",
            "pairwise_Fisher_position_risk_not_emitted_by_geometry_constructor",
        ],
        "classification_dependency": (
            "The classical theorem identifying the order-60 proper rotational "
            "symmetry group of the regular icosahedron with A5 is used only after "
            "the discrete faithful order-60 action has been certified."
        ),
        "claim_boundary": (
            "Exact integer-permutation equivariance of the configured reference "
            "mesh. This is not evidence that repair/readback dynamics selected the "
            "icosahedron and cannot by itself promote an A5-to-SM result."
        ),
    }


def icosahedral_defect_port_report(max_level: int) -> dict:
    """Extract the twelve persistent combinatorial defect-port channels.

    The charge ``q(v)=6-degree(v)`` belongs to the *primal vertex
    triangulation*.  It must not be evaluated on the dual cell graph used for
    the large observer-patch arrays.  The base twelve vertices have ``q=1``;
    regular edgewise subdivision preserves them and introduces only degree-six
    (``q=0``) vertices.
    """

    tower = build_geodesic_icosahedral_tower(max_level)
    rows: list[dict] = []
    persistent_ids: list[int] | None = None
    passed = True
    for mesh in tower.levels:
        degree = np.bincount(mesh.edges.reshape(-1), minlength=mesh.vertex_count)
        charge = 6 - degree
        unit_ids = np.flatnonzero(charge == 1).astype(np.int64)
        nonzero_ids = np.flatnonzero(charge != 0).astype(np.int64)
        if persistent_ids is None:
            persistent_ids = [int(value) for value in unit_ids]
        level_passed = bool(
            int(np.sum(charge)) == 12
            and len(unit_ids) == 12
            and np.array_equal(unit_ids, nonzero_ids)
            and [int(value) for value in unit_ids] == persistent_ids
            and np.all(charge[unit_ids] == 1)
        )
        passed = passed and level_passed
        rows.append(
            {
                "level": mesh.level,
                "primal_vertex_count": mesh.vertex_count,
                "dual_cell_patch_count": mesh.face_count,
                "total_combinatorial_charge": int(np.sum(charge)),
                "unit_defect_vertex_ids": [int(value) for value in unit_ids],
                "nonzero_defect_vertex_ids": [int(value) for value in nonzero_ids],
                "new_vertices_are_degree_six": bool(
                    mesh.level == 0 or np.all(degree[12:] == 6)
                ),
                "charge_vector_hash": _hash_arrays(charge.astype(np.int64)),
                "passed": level_passed,
            }
        )
    assert persistent_ids is not None
    return {
        "schema": "oph.icosahedral_defect_ports.v1",
        "charge_definition": "q(v)=6-degree(v)",
        "charge_domain": "primal_vertex_triangulation",
        "large_array_patch_domain": "dual_spherical_triangular_cells",
        "domains_are_distinct": True,
        "persistent_port_vertex_ids": persistent_ids,
        "persistent_port_count": len(persistent_ids),
        "levels": rows,
        "TWELVE_PERSISTENT_COMBINATORIAL_DEFECT_PORTS_RECEIPT": passed,
        "COMBINATORIAL_DEFECT_CHARGE_IS_PHYSICAL_CURVATURE_RECEIPT": False,
        "PHYSICAL_ATOMIC_DEFECT_PROJECTION_RECEIPT": False,
        "physical_bridge_blockers": [
            "combinatorial_Regge_charge_to_measured_curvature_bridge_missing",
            "atomic_observer_visible_projection_missing",
            "source_selection_and_complete_settlement_missing",
        ],
        "claim_boundary": (
            "The twelve persistent q=1 channels are an exact combinatorial "
            "property of the configured regular refinement. Projected geodesic "
            "curvature and physical observer-visible atomic ports require separate "
            "bridge receipts."
        ),
    }


@lru_cache(maxsize=2)
def build_geodesic_icosahedral_tower(max_level: int) -> GeodesicIcosahedralTower:
    """Build and cache levels ``0..max_level`` of the deterministic tower."""

    if int(max_level) != max_level or max_level < 0:
        raise ValueError("max_level must be a nonnegative integer")
    vertices, faces = _base_icosahedron()
    birth_levels = np.zeros(vertices.shape[0], dtype=np.int16)
    parent_support: tuple[tuple[tuple[int, float], ...], ...] = tuple(
        ((index, 1.0),) for index in range(vertices.shape[0])
    )
    levels: list[GeodesicIcosahedralLevel] = [
        _make_level(
            level=0,
            vertices=vertices,
            faces=faces,
            birth_levels=birth_levels,
            parent_support=parent_support,
            parent_face=None,
            children_by_parent=None,
        )
    ]
    refinements: list[CellRefinementMap] = []
    for level_index in range(1, int(max_level) + 1):
        fine_data = _refine_level(levels[-1], level_index)
        fine_level = _make_level(level=level_index, **fine_data)
        refinements.append(_make_cell_refinement(levels[-1], fine_level))
        levels.append(fine_level)
    return GeodesicIcosahedralTower(tuple(levels), tuple(refinements))


def geodesic_icosahedral_graph(
    level: int,
    *,
    patch_basis: PatchBasis = "cells",
    nominal_patch_count: int | None = None,
) -> nx.Graph:
    """Build the cell-dual or vertex-edge graph for an exact tower level."""

    tower = build_geodesic_icosahedral_tower(level)
    mesh = tower.levels[level]
    graph = nx.Graph()
    if patch_basis == "vertices":
        for node, point in enumerate(mesh.vertices):
            graph.add_node(
                node,
                screen_xyz=tuple(float(value) for value in point),
                geometry_entity="vertex",
                refinement_level=mesh.level,
                birth_level=int(mesh.vertex_birth_levels[node]),
                parent_support=tuple(mesh.vertex_parent_support[node]),
            )
        for left, right in mesh.edges:
            dot = float(np.clip(np.dot(mesh.vertices[left], mesh.vertices[right]), -1.0, 1.0))
            graph.add_edge(int(left), int(right), angular_length=float(math.acos(dot)))
    elif patch_basis == "cells":
        for node, face in enumerate(mesh.faces):
            center = np.sum(mesh.vertices[face], axis=0)
            center /= np.linalg.norm(center)
            attributes = {
                "screen_xyz": tuple(float(value) for value in center),
                "geometry_entity": "spherical_triangular_cell",
                "refinement_level": mesh.level,
                "mesh_vertices": tuple(int(value) for value in face),
                "spherical_area": float(mesh.spherical_face_areas[node]),
            }
            if mesh.parent_face is not None:
                attributes["parent_cell"] = int(mesh.parent_face[node])
                attributes["conditional_expectation_weight"] = float(
                    tower.cell_refinements[-1].conditional_expectation_weights[node]
                )
            graph.add_node(node, **attributes)
        edge_to_faces: dict[tuple[int, int], list[int]] = {}
        for face_index, (a, b, c) in enumerate(mesh.faces):
            for left, right in ((a, b), (b, c), (c, a)):
                key = (int(min(left, right)), int(max(left, right)))
                edge_to_faces.setdefault(key, []).append(face_index)
        for mesh_edge, incident_faces in sorted(edge_to_faces.items()):
            if len(incident_faces) != 2:
                raise AssertionError(
                    f"closed icosahedral sphere edge {mesh_edge} has {len(incident_faces)} incident faces"
                )
            graph.add_edge(
                int(incident_faces[0]),
                int(incident_faces[1]),
                shared_mesh_edge=mesh_edge,
            )
    else:
        raise ValueError(f"unknown icosahedral patch basis: {patch_basis!r}")

    actual_count = graph.number_of_nodes()
    graph.graph.update(
        {
            "geometry_family": "nested_geodesic_icosahedral",
            "patch_basis": patch_basis,
            "refinement_level": mesh.level,
            "frequency": mesh.frequency,
            "vertex_count": mesh.vertex_count,
            "mesh_edge_count": mesh.edge_count,
            "cell_count": mesh.face_count,
            "actual_patch_count": actual_count,
            "nominal_patch_count": nominal_patch_count,
            "nominal_count_is_exact": nominal_patch_count in (None, actual_count),
            "geometry_hash": mesh.geometry_hash,
            "geometry_receipt": mesh.receipt(),
            "paper_multiresolution_regulator_certificate": False,
            "claim_boundary": (
                "Exact nested geodesic geometry and commutative cell-lineage scaffold; "
                "finite detail algebras and presentation circuits remain uninstantiated."
            ),
        }
    )
    if level > 0:
        graph.graph["cell_refinement_receipt"] = tower.cell_refinements[-1].receipt()
    if nominal_patch_count is not None:
        graph.graph["nominal_count_mapping"] = icosahedral_count_bracket(
            nominal_patch_count, patch_basis
        )
    return graph


def geodesic_icosahedral_patch_arrays(
    level: int,
    *,
    patch_basis: PatchBasis = "cells",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(screen_xyz, edge_left, edge_right)`` without NetworkX overhead.

    This is the intended adapter for the large array simulator.  Cell edges
    are dual adjacencies across shared mesh edges; vertex edges are the primal
    triangular-mesh edges.
    """

    mesh = build_geodesic_icosahedral_tower(level).levels[level]
    if patch_basis == "vertices":
        return mesh.vertices, mesh.edges[:, 0], mesh.edges[:, 1]
    if patch_basis != "cells":
        raise ValueError(f"unknown icosahedral patch basis: {patch_basis!r}")

    points = np.sum(mesh.vertices[mesh.faces], axis=1)
    points /= np.linalg.norm(points, axis=1, keepdims=True)
    faces = mesh.faces
    mesh_edges = np.concatenate(
        (faces[:, (0, 1)], faces[:, (1, 2)], faces[:, (2, 0)]),
        axis=0,
    )
    low = np.minimum(mesh_edges[:, 0], mesh_edges[:, 1])
    high = np.maximum(mesh_edges[:, 0], mesh_edges[:, 1])
    edge_keys = low * mesh.vertex_count + high
    face_ids = np.tile(np.arange(mesh.face_count, dtype=np.int64), 3)
    order = np.argsort(edge_keys, kind="stable")
    sorted_keys = edge_keys[order]
    sorted_faces = face_ids[order]
    if sorted_keys.size % 2 or not np.array_equal(sorted_keys[0::2], sorted_keys[1::2]):
        raise AssertionError("every closed spherical mesh edge must have exactly two incident cells")
    left = np.minimum(sorted_faces[0::2], sorted_faces[1::2])
    right = np.maximum(sorted_faces[0::2], sorted_faces[1::2])
    return _readonly(points), _readonly(left), _readonly(right)


def resolve_icosahedral_level(
    requested_count: int,
    *,
    patch_basis: PatchBasis,
    policy: Literal["exact", "nearest", "floor", "ceil"] = "exact",
) -> tuple[int, dict]:
    """Resolve a count under an explicit non-fabricating level policy."""

    bracket = icosahedral_count_bracket(requested_count, patch_basis)
    if policy == "exact":
        if not bracket["exact_supported_count"]:
            raise UnsupportedIcosahedralPatchCount(
                f"{requested_count} is not an exact {patch_basis} count for a dyadic "
                "geodesic icosahedral mesh; supported bracket is "
                f"lower={bracket['lower']} upper={bracket['upper']}. Set an explicit "
                "refinement_level and use nominal_patch_count as a label, or choose "
                "patch_count_policy: nearest/floor/ceil."
            )
        return int(bracket["nearest"]["level"]), bracket
    if policy == "nearest":
        return int(bracket["nearest"]["level"]), bracket
    if policy == "floor":
        if bracket["lower"] is None:
            raise UnsupportedIcosahedralPatchCount(
                f"no supported {patch_basis} count is at or below {requested_count}"
            )
        return int(bracket["lower"]["level"]), bracket
    if policy == "ceil":
        return int(bracket["upper"]["level"]), bracket
    raise ValueError(f"unknown icosahedral patch count policy: {policy!r}")


def _base_icosahedron() -> tuple[np.ndarray, np.ndarray]:
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    vertices = np.asarray(
        [
            (-1.0, phi, 0.0),
            (1.0, phi, 0.0),
            (-1.0, -phi, 0.0),
            (1.0, -phi, 0.0),
            (0.0, -1.0, phi),
            (0.0, 1.0, phi),
            (0.0, -1.0, -phi),
            (0.0, 1.0, -phi),
            (phi, 0.0, -1.0),
            (phi, 0.0, 1.0),
            (-phi, 0.0, -1.0),
            (-phi, 0.0, 1.0),
        ],
        dtype=float,
    )
    vertices /= np.linalg.norm(vertices, axis=1, keepdims=True)
    faces = np.asarray(
        [
            (0, 11, 5),
            (0, 5, 1),
            (0, 1, 7),
            (0, 7, 10),
            (0, 10, 11),
            (1, 5, 9),
            (5, 11, 4),
            (11, 10, 2),
            (10, 7, 6),
            (7, 1, 8),
            (3, 9, 4),
            (3, 4, 2),
            (3, 2, 6),
            (3, 6, 8),
            (3, 8, 9),
            (4, 9, 5),
            (2, 4, 11),
            (6, 2, 10),
            (8, 6, 7),
            (9, 8, 1),
        ],
        dtype=np.int64,
    )
    return _clean_coordinates(vertices), _orient_faces_outward(vertices, faces)


def _refine_level(coarse: GeodesicIcosahedralLevel, level: int) -> dict:
    vertices = [point.copy() for point in coarse.vertices]
    birth_levels = [int(value) for value in coarse.vertex_birth_levels]
    # Inherited vertices embed identically at every adjacent refinement.  New
    # vertices carry their two coarse edge endpoints as immediate parents.
    parent_support: list[tuple[tuple[int, float], ...]] = [
        ((index, 1.0),) for index in range(coarse.vertex_count)
    ]
    midpoint_ids: dict[tuple[int, int], int] = {}

    def midpoint(left: int, right: int) -> int:
        key = (min(left, right), max(left, right))
        if key not in midpoint_ids:
            point = coarse.vertices[key[0]] + coarse.vertices[key[1]]
            point /= np.linalg.norm(point)
            midpoint_ids[key] = len(vertices)
            vertices.append(point)
            birth_levels.append(level)
            parent_support.append(((key[0], 0.5), (key[1], 0.5)))
        return midpoint_ids[key]

    fine_faces: list[tuple[int, int, int]] = []
    parent_face: list[int] = []
    children_by_parent: list[tuple[int, int, int, int]] = []
    for parent, (a_value, b_value, c_value) in enumerate(coarse.faces):
        a, b, c = int(a_value), int(b_value), int(c_value)
        ab = midpoint(a, b)
        bc = midpoint(b, c)
        ca = midpoint(c, a)
        child_start = len(fine_faces)
        fine_faces.extend(((a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)))
        parent_face.extend((parent, parent, parent, parent))
        children_by_parent.append(
            (child_start, child_start + 1, child_start + 2, child_start + 3)
        )

    return {
        "vertices": _clean_coordinates(np.asarray(vertices, dtype=float)),
        "faces": np.asarray(fine_faces, dtype=np.int64),
        "birth_levels": np.asarray(birth_levels, dtype=np.int16),
        "parent_support": tuple(parent_support),
        "parent_face": np.asarray(parent_face, dtype=np.int64),
        "children_by_parent": tuple(children_by_parent),
    }


def _make_level(
    *,
    level: int,
    vertices: np.ndarray,
    faces: np.ndarray,
    birth_levels: np.ndarray,
    parent_support: tuple[tuple[tuple[int, float], ...], ...],
    parent_face: np.ndarray | None,
    children_by_parent: tuple[tuple[int, int, int, int], ...] | None,
) -> GeodesicIcosahedralLevel:
    vertices = _clean_coordinates(vertices)
    faces = _orient_faces_outward(vertices, faces)
    edges = _mesh_edges(faces)
    areas = _spherical_face_areas(vertices, faces)
    geometry_hash = _geometry_hash(level, vertices, faces, edges, parent_face)
    return GeodesicIcosahedralLevel(
        level=level,
        frequency=2**level,
        vertices=_readonly(vertices),
        faces=_readonly(faces.astype(np.int64)),
        edges=_readonly(edges.astype(np.int64)),
        spherical_face_areas=_readonly(areas),
        vertex_birth_levels=_readonly(birth_levels.astype(np.int16)),
        vertex_parent_support=parent_support,
        parent_face=None if parent_face is None else _readonly(parent_face.astype(np.int64)),
        children_by_parent_face=children_by_parent,
        geometry_hash=geometry_hash,
    )


def _make_cell_refinement(
    coarse: GeodesicIcosahedralLevel,
    fine: GeodesicIcosahedralLevel,
) -> CellRefinementMap:
    if fine.parent_face is None or fine.children_by_parent_face is None:
        raise ValueError("fine level has no parent-face lineage")
    weights = np.zeros(fine.face_count, dtype=float)
    for parent, children in enumerate(fine.children_by_parent_face):
        child_ids = np.asarray(children, dtype=np.int64)
        child_areas = fine.spherical_face_areas[child_ids]
        weights[child_ids] = child_areas / np.sum(child_areas)
    normalization_residual = float(
        max(
            abs(float(np.sum(weights[np.asarray(children, dtype=np.int64)])) - 1.0)
            for children in fine.children_by_parent_face
        )
    )
    coarse_reference = coarse.spherical_face_areas / (4.0 * math.pi)
    fine_reference = fine.spherical_face_areas / (4.0 * math.pi)
    aggregated = np.bincount(
        fine.parent_face,
        weights=fine_reference,
        minlength=coarse.face_count,
    )
    state_residual = float(np.max(np.abs(aggregated - coarse_reference)))
    payload = {
        "coarse_level": coarse.level,
        "fine_level": fine.level,
        "child_to_parent": fine.parent_face.tolist(),
        "weights": np.round(weights, 15).tolist(),
    }
    map_hash = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return CellRefinementMap(
        coarse_level=coarse.level,
        fine_level=fine.level,
        child_to_parent=fine.parent_face,
        children_by_parent=fine.children_by_parent_face,
        conditional_expectation_weights=_readonly(weights),
        coarse_reference_weights=_readonly(coarse_reference),
        fine_reference_weights=_readonly(fine_reference),
        normalization_residual=normalization_residual,
        state_preservation_residual=state_residual,
        map_hash=map_hash,
    )


def _clean_coordinates(vertices: np.ndarray) -> np.ndarray:
    result = np.asarray(vertices, dtype=float).copy()
    result /= np.linalg.norm(result, axis=1, keepdims=True)
    result[np.abs(result) < 5.0e-16] = 0.0
    return result


def _orient_faces_outward(vertices: np.ndarray, faces: np.ndarray) -> np.ndarray:
    result = np.asarray(faces, dtype=np.int64).copy()
    for index, (a, b, c) in enumerate(result):
        normal = np.cross(vertices[b] - vertices[a], vertices[c] - vertices[a])
        center = vertices[a] + vertices[b] + vertices[c]
        if float(np.dot(normal, center)) < 0.0:
            result[index, 1], result[index, 2] = result[index, 2], result[index, 1]
    return result


def _faces_are_outward(vertices: np.ndarray, faces: np.ndarray) -> bool:
    for a, b, c in faces:
        normal = np.cross(vertices[b] - vertices[a], vertices[c] - vertices[a])
        center = vertices[a] + vertices[b] + vertices[c]
        if float(np.dot(normal, center)) <= 0.0:
            return False
    return True


def _mesh_edges(faces: np.ndarray) -> np.ndarray:
    edge_set: set[tuple[int, int]] = set()
    for a_value, b_value, c_value in faces:
        a, b, c = int(a_value), int(b_value), int(c_value)
        edge_set.update(
            {
                (min(a, b), max(a, b)),
                (min(b, c), max(b, c)),
                (min(c, a), max(c, a)),
            }
        )
    return np.asarray(sorted(edge_set), dtype=np.int64)


def _spherical_face_areas(vertices: np.ndarray, faces: np.ndarray) -> np.ndarray:
    areas = np.zeros(faces.shape[0], dtype=float)
    for index, (a, b, c) in enumerate(faces):
        u, v, w = vertices[a], vertices[b], vertices[c]
        numerator = abs(float(np.dot(u, np.cross(v, w))))
        denominator = 1.0 + float(np.dot(u, v) + np.dot(v, w) + np.dot(w, u))
        areas[index] = 2.0 * math.atan2(numerator, denominator)
    return areas


def _geometry_hash(
    level: int,
    vertices: np.ndarray,
    faces: np.ndarray,
    edges: np.ndarray,
    parent_face: np.ndarray | None,
) -> str:
    digest = hashlib.sha256()
    digest.update(b"oph.nested_geodesic_icosahedral.v1\0")
    digest.update(int(level).to_bytes(4, "little", signed=False))
    rounded = np.round(vertices, decimals=15)
    rounded[np.abs(rounded) < 0.5e-15] = 0.0
    digest.update(np.asarray(rounded, dtype="<f8").tobytes())
    digest.update(np.asarray(faces, dtype="<i8").tobytes())
    digest.update(np.asarray(edges, dtype="<i8").tobytes())
    if parent_face is not None:
        digest.update(np.asarray(parent_face, dtype="<i8").tobytes())
    return digest.hexdigest()


def _hash_arrays(*arrays: np.ndarray) -> str:
    digest = hashlib.sha256()
    digest.update(b"oph.icosahedral.discrete-action.v1\0")
    for array in arrays:
        values = np.ascontiguousarray(array)
        digest.update(str(values.dtype).encode("ascii"))
        digest.update(b"\0")
        digest.update(json.dumps(list(values.shape)).encode("ascii"))
        digest.update(b"\0")
        digest.update(values.tobytes())
    return digest.hexdigest()


@lru_cache(maxsize=1)
def _base_rotation_group() -> tuple[tuple[np.ndarray, ...], tuple[tuple[int, ...], ...]]:
    """Recover all proper base-icosahedron rotations deterministically."""

    base = build_geodesic_icosahedral_tower(0).levels[0]
    reference = base.vertices[base.faces[0]].T
    inverse_reference = np.linalg.inv(reference)
    rotations: list[np.ndarray] = []
    permutations: list[tuple[int, ...]] = []
    for face in base.faces:
        for shift in range(3):
            target = base.vertices[np.roll(face, -shift)].T
            rotation = target @ inverse_reference
            if float(np.linalg.det(rotation)) <= 0.0:
                continue
            permutation, residual = _coordinate_permutation(base.vertices, rotation)
            candidate = tuple(int(value) for value in permutation)
            if residual > 5.0e-12:
                raise AssertionError("candidate icosahedral rotation misses the base vertices")
            if np.max(np.abs(rotation.T @ rotation - np.eye(3))) > 5.0e-12:
                raise AssertionError("candidate icosahedral rotation is not orthogonal")
            if candidate not in permutations:
                rotations.append(_readonly(rotation))
                permutations.append(candidate)
    if len(permutations) != 60:
        raise AssertionError(
            f"regular icosahedron must have 60 proper rotations, found {len(permutations)}"
        )
    return tuple(rotations), tuple(permutations)


def _coordinate_permutation(
    points: np.ndarray,
    rotation: np.ndarray,
) -> tuple[np.ndarray, float]:
    mapped = np.asarray(points, dtype=float) @ np.asarray(rotation, dtype=float).T
    distances, indices = cKDTree(np.asarray(points, dtype=float)).query(mapped, k=1)
    permutation = np.asarray(indices, dtype=np.int64)
    if len(np.unique(permutation)) != points.shape[0]:
        raise AssertionError("rotation did not induce a vertex permutation")
    return permutation, float(np.max(distances)) if distances.size else 0.0


def _compose_permutations(
    left: tuple[int, ...],
    right: tuple[int, ...],
) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(len(left)))


def _inverse_permutation(permutation: tuple[int, ...]) -> tuple[int, ...]:
    inverse = [0] * len(permutation)
    for source, target in enumerate(permutation):
        inverse[target] = source
    return tuple(inverse)


def _permutation_order(permutation: tuple[int, ...]) -> int:
    identity = tuple(range(len(permutation)))
    value = identity
    for order in range(1, 121):
        value = _compose_permutations(permutation, value)
        if value == identity:
            return order
    raise AssertionError("icosahedral permutation order exceeded the group bound")


def _generated_permutation_subgroup(
    generators: Iterable[tuple[int, ...]],
) -> set[tuple[int, ...]]:
    generators = tuple(generators)
    if not generators:
        return set()
    identity = tuple(range(len(generators[0])))
    seen = {identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            for candidate in (
                _compose_permutations(generator, current),
                _compose_permutations(current, generator),
            ):
                if candidate not in seen:
                    seen.add(candidate)
                    frontier.append(candidate)
    return seen


@lru_cache(maxsize=1)
def _find_a5_generator_pair(
    permutations: tuple[tuple[int, ...], ...],
) -> tuple[int, int]:
    for left_index, left in enumerate(permutations):
        if _permutation_order(left) != 3:
            continue
        for right_index, right in enumerate(permutations):
            if _permutation_order(right) != 5:
                continue
            if len(_generated_permutation_subgroup((left, right))) == 60:
                return left_index, right_index
    raise AssertionError("failed to find order-3/order-5 generators of the A5 action")
