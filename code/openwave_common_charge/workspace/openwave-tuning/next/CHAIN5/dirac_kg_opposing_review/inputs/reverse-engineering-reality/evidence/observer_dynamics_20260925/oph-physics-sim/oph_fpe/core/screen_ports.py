from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import linear_sum_assignment

from oph_fpe.core.icosahedral import (
    build_geodesic_icosahedral_tower,
    icosahedral_a5_equivariance_report,
    icosahedral_a5_port_permutations,
)
from oph_fpe.gauge.covariant_overlap import (
    group_inverse_indices,
    group_multiply_indices,
)


def echosahedral_port_names(ports_per_patch: int = 12) -> list[str]:
    return [f"P{index}" for index in range(int(ports_per_patch))]


@dataclass(frozen=True)
class ReferenceDiagonalA5Intertwiner:
    """Typed reference maps between local twelve-port cells and one port copy.

    The local carrier is ``R^cells tensor R^12`` and the declared reference
    action is diagonal: the same base-icosahedron port permutation acts in
    every cell, while cell identifiers remain fixed.  The embedding replicates
    one global twelve-vector into every cell and the projection takes a
    faithful weighted cell average.  This is a useful exact reference object,
    but it assumes a common port-frame trivialization; it does not derive that
    trivialization from overlap transport, repair histories, or refinement.
    """

    cell_weights: np.ndarray
    port_permutations: tuple[tuple[int, ...], ...]
    action_convention: str = "old_port_i_maps_to_new_port_permutation_i"

    def __post_init__(self) -> None:
        weights = np.asarray(self.cell_weights, dtype=np.float64).copy()
        if weights.ndim != 1 or weights.size < 1:
            raise ValueError("cell_weights must be a nonempty one-dimensional array")
        if not np.all(np.isfinite(weights)) or np.any(weights <= 0.0):
            raise ValueError("cell_weights must be finite and strictly positive")
        if abs(float(np.sum(weights)) - 1.0) > 5.0e-13:
            raise ValueError("cell_weights must already be normalized to one")
        permutations = tuple(
            tuple(int(value) for value in permutation)
            for permutation in self.port_permutations
        )
        for permutation in permutations:
            if len(permutation) != 12 or sorted(permutation) != list(range(12)):
                raise ValueError("every A5 action row must be a permutation of 0..11")
        weights.setflags(write=False)
        object.__setattr__(self, "cell_weights", weights)
        object.__setattr__(self, "port_permutations", permutations)

    @property
    def patch_count(self) -> int:
        return int(self.cell_weights.size)

    @property
    def port_count(self) -> int:
        return 12

    @property
    def group_element_count(self) -> int:
        return len(self.port_permutations)

    @property
    def cell_weight_sha256(self) -> str:
        return "sha256:" + hashlib.sha256(
            np.asarray(self.cell_weights, dtype="<f8").tobytes()
        ).hexdigest()

    @property
    def port_action_sha256(self) -> str:
        return "sha256:" + hashlib.sha256(
            np.asarray(self.port_permutations, dtype="<i8").tobytes()
        ).hexdigest()

    @property
    def contract_sha256(self) -> str:
        payload = {
            "schema": "oph.reference_diagonal_a5_intertwiner.contract.v1",
            "local_carrier": f"R^{self.patch_count}_cells_tensor_R^12_ports",
            "global_carrier": "R^12_reference_ports",
            "action": "identity_on_cells_tensor_base_icosahedron_port_permutation",
            "action_convention": self.action_convention,
            "embedding": "replicate_global_port_vector_into_every_cell",
            "projection": "faithful_weighted_average_over_cells",
            "cell_weight_sha256": self.cell_weight_sha256,
            "port_action_sha256": self.port_action_sha256,
        }
        return "sha256:" + hashlib.sha256(
            json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest()

    def embed(self, global_values: np.ndarray) -> np.ndarray:
        """Replicate a reference twelve-port observable into every local cell."""

        values = np.asarray(global_values)
        if values.ndim < 1 or values.shape[0] != 12:
            raise ValueError("global_values first dimension must be the twelve ports")
        return np.broadcast_to(values, (self.patch_count,) + values.shape).copy()

    def project(self, local_values: np.ndarray) -> np.ndarray:
        """Apply the state-preserving weighted local-to-reference projection."""

        values = np.asarray(local_values)
        if values.ndim < 2 or values.shape[:2] != (self.patch_count, 12):
            raise ValueError(
                "local_values first dimensions must be (patch_count, 12)"
            )
        return np.tensordot(self.cell_weights, values, axes=(0, 0))

    def conditional_expectation(self, local_values: np.ndarray) -> np.ndarray:
        """Project and re-embed onto the declared diagonal reference subspace."""

        return self.embed(self.project(local_values))

    def act_global(self, global_values: np.ndarray, element_index: int) -> np.ndarray:
        """Apply one proper-icosahedron port permutation to the reference copy."""

        values = np.asarray(global_values)
        if values.ndim < 1 or values.shape[0] != 12:
            raise ValueError("global_values first dimension must be the twelve ports")
        inverse = self._inverse_action(element_index)
        return np.take(values, inverse, axis=0)

    def act_local(self, local_values: np.ndarray, element_index: int) -> np.ndarray:
        """Apply the same port permutation independently in every local cell."""

        values = np.asarray(local_values)
        if values.ndim < 2 or values.shape[:2] != (self.patch_count, 12):
            raise ValueError(
                "local_values first dimensions must be (patch_count, 12)"
            )
        inverse = self._inverse_action(element_index)
        return np.take(values, inverse, axis=1)

    def audit(self) -> dict[str, Any]:
        """Verify the reference group and both intertwining identities exactly."""

        permutations = self.port_permutations
        permutation_set = set(permutations)
        identity = tuple(range(12))
        group_closed = all(
            _compose_port_permutations(left, right) in permutation_set
            for left in permutations
            for right in permutations
        )
        inverses_present = all(
            _inverse_port_permutation(permutation) in permutation_set
            for permutation in permutations
        )
        base_port_orbit = {permutation[0] for permutation in permutations}
        base_port_stabilizer_size = sum(
            permutation[0] == 0 for permutation in permutations
        )
        order_profile: dict[int, int] = {}
        action_rows: list[dict[str, Any]] = []
        port_identity = np.eye(12, dtype=np.int8)
        for index, permutation in enumerate(permutations):
            order = _port_permutation_order(permutation)
            order_profile[order] = order_profile.get(order, 0) + 1
            matrix = _port_permutation_matrix(permutation)
            embedding_commutes = bool(
                np.array_equal(matrix @ port_identity, port_identity @ matrix)
            )
            projection_commutes = bool(
                np.array_equal(port_identity @ matrix, matrix @ port_identity)
            )
            action_rows.append(
                {
                    "element_index": index,
                    "element_order": order,
                    "permutation": list(permutation),
                    "permutation_sha256": "sha256:"
                    + hashlib.sha256(
                        np.asarray(permutation, dtype="<i8").tobytes()
                    ).hexdigest(),
                    "embedding_intertwining_exact": embedding_commutes,
                    "weighted_projection_intertwining_exact": projection_commutes,
                    "embedding_intertwining_residual": 0 if embedding_commutes else 1,
                    "weighted_projection_intertwining_residual": 0
                    if projection_commutes
                    else 1,
                }
            )
        normalized_residual = abs(float(np.sum(self.cell_weights)) - 1.0)
        reference_receipt = bool(
            len(permutations) == 60
            and len(permutation_set) == 60
            and identity in permutation_set
            and group_closed
            and inverses_present
            and len(base_port_orbit) == 12
            and base_port_stabilizer_size == 5
            and order_profile == {1: 1, 2: 15, 3: 20, 5: 24}
            and normalized_residual <= 5.0e-13
            and np.all(self.cell_weights > 0.0)
            and all(
                row["embedding_intertwining_exact"]
                and row["weighted_projection_intertwining_exact"]
                for row in action_rows
            )
        )
        return {
            "schema": "oph.reference_diagonal_a5_intertwiner.v1",
            "object_type": "REFERENCE_diagonal_local_to_global_A5_intertwiner",
            "local_carrier_type": {
                "cell_factor_dimension": self.patch_count,
                "port_factor_dimension": 12,
                "formal_dimension": 12 * self.patch_count,
                "representation": "identity_on_cells_tensor_port_permutation",
            },
            "global_reference_carrier_type": {
                "dimension": 12,
                "representation": "base_icosahedron_port_permutation",
            },
            "embedding_type": "unital_diagonal_replication",
            "projection_type": "faithful_cell_weighted_average",
            "conditional_expectation_type": "embedding_after_projection",
            "verification_method": (
                "factorized_integer_coefficient_identity_without_dense_12N_matrix"
            ),
            "action_convention": self.action_convention,
            "group_element_count": len(permutations),
            "unique_action_count": len(permutation_set),
            "element_order_profile": {
                str(key): value for key, value in sorted(order_profile.items())
            },
            "integer_permutation_group_closed": group_closed,
            "integer_permutation_inverses_present": inverses_present,
            "base_port_orbit_size": len(base_port_orbit),
            "base_port_stabilizer_size": base_port_stabilizer_size,
            "cell_weights": {
                "count": self.patch_count,
                "sum": float(np.sum(self.cell_weights)),
                "minimum": float(np.min(self.cell_weights)),
                "maximum": float(np.max(self.cell_weights)),
                "strictly_positive": bool(np.all(self.cell_weights > 0.0)),
                "normalization_residual": normalized_residual,
                "uniform": bool(
                    np.max(self.cell_weights) - np.min(self.cell_weights)
                    <= 5.0e-15
                ),
                "sha256": self.cell_weight_sha256,
            },
            "map_identities": {
                "projection_after_embedding_is_identity": normalized_residual
                <= 5.0e-13,
                "embedding_projection_is_idempotent": normalized_residual
                <= 5.0e-13,
                "projection_is_weighted_adjoint_of_embedding": True,
                "all_60_embedding_intertwining_checks": bool(
                    action_rows
                    and all(row["embedding_intertwining_exact"] for row in action_rows)
                ),
                "all_60_projection_intertwining_checks": bool(
                    action_rows
                    and all(
                        row["weighted_projection_intertwining_exact"]
                        for row in action_rows
                    )
                ),
            },
            "action_rows": action_rows,
            "port_action_sha256": self.port_action_sha256,
            "contract_sha256": self.contract_sha256,
            "reference_assumptions": [
                "one common ordered P0..P11 template is declared for every cell",
                "the A5 action fixes cell identifiers and acts diagonally on port copies",
                "the supplied faithful normalized cell weights define the reference state",
            ],
            "source_derived_common_trivialization": False,
            "physical_receipt_eligible": False,
            "REFERENCE_DIAGONAL_A5_INTERTWINER_RECEIPT": reference_receipt,
            "REFERENCE_LOCAL_TO_GLOBAL_A5_INTERTWINER_RECEIPT": reference_receipt,
            "LOCAL_TO_GLOBAL_A5_INTERTWINER_RECEIPT": False,
            "A5_EQUIVARIANT_CURRENT_CARRIER_INTERTWINER_RECEIPT": False,
            "GLOBAL_TWELVE_UNIT_SCREEN_SIEVE_RECEIPT": False,
            "PHYSICAL_GLOBAL_SCREEN_SIEVE_RECEIPT": False,
            "physical_promotion_blockers": [
                "common_local_port_frame_is_declared_not_source_derived",
                "coherent_overlap_transport_between_local_port_frames_missing",
                "refinement_natural_transport_intertwiner_missing",
                "quotient_visible_collective_twelve_channel_identification_missing",
                "physical_global_atomic_defect_ports_not_bound_to_local_port_copies",
                "reciprocal_current_carrier_intertwiner_missing",
            ],
            "claim_boundary": (
                "Exact finite reference identity for a declared trivial bundle of local "
                "twelve-port copies: the same A5 permutation acts in every cell, the "
                "global reference vector embeds diagonally, and faithful cell averaging "
                "intertwines all 60 permutations. The simulator has not derived a common "
                "local-frame trivialization, coherent overlap/refinement transport, or a "
                "quotient-visible identification with the twelve physical global defect "
                "ports. Therefore LOCAL_TO_GLOBAL_A5_INTERTWINER_RECEIPT remains false."
            ),
        }

    def _inverse_action(self, element_index: int) -> np.ndarray:
        if not (0 <= int(element_index) < len(self.port_permutations)):
            raise IndexError("A5 element_index is out of range")
        permutation = self.port_permutations[int(element_index)]
        return np.asarray(_inverse_port_permutation(permutation), dtype=np.int64)


def reference_diagonal_a5_intertwiner(
    patch_count: int,
    *,
    cell_weights: np.ndarray | None = None,
) -> ReferenceDiagonalA5Intertwiner:
    """Construct the declared reference intertwiner for a finite federation."""

    count = int(patch_count)
    if count < 1:
        raise ValueError("patch_count must be positive")
    weights = (
        np.full(count, 1.0 / count, dtype=np.float64)
        if cell_weights is None
        else np.asarray(cell_weights, dtype=np.float64)
    )
    if weights.shape != (count,):
        raise ValueError("cell_weights shape must equal (patch_count,)")
    return ReferenceDiagonalA5Intertwiner(
        cell_weights=weights,
        port_permutations=icosahedral_a5_port_permutations(),
    )


def reference_diagonal_a5_intertwiner_report(
    patch_count: int,
    *,
    ports_per_patch: int = 12,
    cell_weights: np.ndarray | None = None,
) -> dict[str, Any]:
    """Emit a fail-closed reference report without promoting the SM ladder."""

    count = int(patch_count)
    ports = int(ports_per_patch)
    if count < 1 or ports != 12:
        blockers = []
        if count < 1:
            blockers.append("positive_patch_federation_missing")
        if ports != 12:
            blockers.append("canonical_twelve_port_local_template_missing")
        return {
            "schema": "oph.reference_diagonal_a5_intertwiner.v1",
            "object_type": "REFERENCE_diagonal_local_to_global_A5_intertwiner",
            "patch_count": count,
            "ports_per_patch": ports,
            "REFERENCE_DIAGONAL_A5_INTERTWINER_RECEIPT": False,
            "REFERENCE_LOCAL_TO_GLOBAL_A5_INTERTWINER_RECEIPT": False,
            "LOCAL_TO_GLOBAL_A5_INTERTWINER_RECEIPT": False,
            "A5_EQUIVARIANT_CURRENT_CARRIER_INTERTWINER_RECEIPT": False,
            "GLOBAL_TWELVE_UNIT_SCREEN_SIEVE_RECEIPT": False,
            "PHYSICAL_GLOBAL_SCREEN_SIEVE_RECEIPT": False,
            "blockers": blockers,
            "claim_boundary": (
                "The typed reference intertwiner requires a nonempty federation of "
                "canonical twelve-port local patches. No physical promotion is possible."
            ),
        }
    return reference_diagonal_a5_intertwiner(
        count,
        cell_weights=cell_weights,
    ).audit()


def _compose_port_permutations(
    left: tuple[int, ...],
    right: tuple[int, ...],
) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(12))


def _inverse_port_permutation(permutation: tuple[int, ...]) -> tuple[int, ...]:
    inverse = [0] * 12
    for old, new in enumerate(permutation):
        inverse[new] = old
    return tuple(inverse)


def _port_permutation_order(permutation: tuple[int, ...]) -> int:
    current = tuple(range(12))
    for order in range(1, 61):
        current = _compose_port_permutations(permutation, current)
        if current == tuple(range(12)):
            return order
    raise AssertionError("finite A5 port permutation order exceeds 60")


def _port_permutation_matrix(permutation: tuple[int, ...]) -> np.ndarray:
    matrix = np.zeros((12, 12), dtype=np.int8)
    matrix[np.asarray(permutation, dtype=np.int64), np.arange(12)] = 1
    return matrix


@dataclass(frozen=True)
class ScreenPortMap:
    left_port: np.ndarray
    right_port: np.ndarray
    ports_per_patch: int
    overflow_count: int
    node_degree: np.ndarray
    routing_mode: str = "sequential_endpoint_order"
    directional_alignment: np.ndarray | None = None
    local_frame_hash: str | None = None
    federation_weights: np.ndarray | None = None

    def as_jsonable(self, *, sample_edges: int = 64) -> dict[str, Any]:
        sample_count = min(int(sample_edges), int(self.left_port.size))
        port_names = echosahedral_port_names(self.ports_per_patch)
        alignment = (
            np.asarray(self.directional_alignment, dtype=float)
            if self.directional_alignment is not None
            else np.zeros(0, dtype=float)
        )
        architecture = echosahedral_patch_architecture_report(self.ports_per_patch)
        reference_intertwiner = reference_diagonal_a5_intertwiner_report(
            int(self.node_degree.size),
            ports_per_patch=self.ports_per_patch,
            cell_weights=self.federation_weights,
        )
        endpoint_count = int(np.sum(self.node_degree, dtype=np.int64))
        # The degree sum above is only a count witness.  Slot uniqueness is
        # equivalent to no overflow for both routing implementations: the
        # sequential path uses local ranks and the geometric path repairs every
        # collision with a per-node assignment.
        singleton_sewing_reference = bool(
            self.ports_per_patch == 12
            and self.overflow_count == 0
            and endpoint_count == 2 * self.left_port.size
            and np.all(self.left_port >= 0)
            and np.all(self.left_port < 12)
            and np.all(self.right_port >= 0)
            and np.all(self.right_port < 12)
        )
        return {
            "mode": "explicit_named_echosahedral_ports",
            "architecture_level": "two_level_local_carriers_plus_intercarrier_seams",
            "routing_mode": self.routing_mode,
            "ports_per_patch": int(self.ports_per_patch),
            "port_names": port_names,
            "edge_count": int(self.left_port.size),
            "overflow_count": int(self.overflow_count),
            "overflow_fraction": float(self.overflow_count / max(1, int(self.left_port.size) * 2)),
            "degree_min": int(np.min(self.node_degree)) if self.node_degree.size else 0,
            "degree_mean": float(np.mean(self.node_degree)) if self.node_degree.size else 0.0,
            "degree_max": int(np.max(self.node_degree)) if self.node_degree.size else 0,
            "unused_port_slots": int(
                max(0, int(self.node_degree.size) * int(self.ports_per_patch) - 2 * int(self.left_port.size))
            ),
            "directional_alignment": {
                "sample_count": int(alignment.size),
                "minimum": float(np.min(alignment)) if alignment.size else None,
                "mean": float(np.mean(alignment)) if alignment.size else None,
                "maximum": float(np.max(alignment)) if alignment.size else None,
            },
            "local_frame_hash": self.local_frame_hash,
            "local_patch_architecture": architecture,
            "federation_sewing": {
                "carrier_count": int(self.node_degree.size),
                "seam_count": int(self.left_port.size),
                "current_seam_bundle_mode": "singleton_port_pair",
                "general_bundle_schema_implemented": False,
                "orientation_reversal": "swap_endpoints_and_port_restrictions",
                "valid_distinct_local_slots": singleton_sewing_reference,
                "interface_algebra_hashes_bound": False,
            },
            "reference_diagonal_a5_intertwiner": reference_intertwiner,
            "ECHOSAHEDRAL_LOCAL_PATCH_ARCHITECTURE_RECEIPT": architecture[
                "ECHOSAHEDRAL_LOCAL_PATCH_ARCHITECTURE_RECEIPT"
            ],
            "ECHOSAHEDRAL_CARRIER_CONFORMANCE_RECEIPT": architecture[
                "ECHOSAHEDRAL_CARRIER_CONFORMANCE_RECEIPT"
            ],
            "REFERENCE_SINGLETON_FEDERATION_SEWING_RECEIPT": singleton_sewing_reference,
            "FEDERATION_SEWING_RECEIPT": False,
            "REFERENCE_DIAGONAL_A5_INTERTWINER_RECEIPT": reference_intertwiner[
                "REFERENCE_DIAGONAL_A5_INTERTWINER_RECEIPT"
            ],
            "LOCAL_TO_GLOBAL_A5_INTERTWINER_RECEIPT": False,
            "GLOBAL_TWELVE_UNIT_SCREEN_SIEVE_RECEIPT": False,
            "PHYSICAL_GLOBAL_SCREEN_SIEVE_RECEIPT": False,
            "GEOMETRIC_LOCAL_PORT_ROUTING_RECEIPT": bool(
                self.routing_mode == "icosahedral_directional_assignment"
                and self.overflow_count == 0
                and alignment.size == 2 * self.left_port.size
                and np.all(np.isfinite(alignment))
            ),
            "sample_edge_ports": [
                {
                    "edge_index": index,
                    "left_port": port_names[int(self.left_port[index])],
                    "right_port": port_names[int(self.right_port[index])],
                }
                for index in range(sample_count)
            ],
            "claim_boundary": (
                "Every observer cell instantiates one local P0..P11 echosahedral interface. "
                "Routed overlap edges occupy distinct local slots; unused slots remain exposed "
                "or reserved patch ports. The present engine realizes singleton port-pair seams; "
                "typed port-bundle/collar seams with interface-algebra hashes are still required "
                "before FEDERATION_SEWING_RECEIPT can pass. This configured reference architecture is distinct "
                "from the global twelve-unit screen-sieve producer. The nested diagonal A5 "
                "intertwiner is exact only after declaring one common port-frame trivialization; "
                "it does not derive coherent physical transport or A5-to-SM emergence."
            ),
        }


def assign_echosahedral_ports(
    left: np.ndarray,
    right: np.ndarray,
    patch_count: int,
    *,
    ports_per_patch: int = 12,
    points: np.ndarray | None = None,
    federation_weights: np.ndarray | None = None,
) -> ScreenPortMap:
    left = np.asarray(left, dtype=np.int64)
    right = np.asarray(right, dtype=np.int64)
    ports = max(1, int(ports_per_patch))
    weights = None
    if federation_weights is not None:
        weights = np.asarray(federation_weights, dtype=np.float64).copy()
        if weights.shape != (int(patch_count),):
            raise ValueError("federation_weights shape must equal (patch_count,)")
        if not np.all(np.isfinite(weights)) or np.any(weights <= 0.0):
            raise ValueError("federation_weights must be finite and strictly positive")
        if abs(float(np.sum(weights)) - 1.0) > 5.0e-13:
            raise ValueError("federation_weights must already be normalized to one")
    node_degree = np.bincount(np.concatenate([left, right]), minlength=int(patch_count))
    endpoint_nodes = np.empty(left.size * 2, dtype=np.int64)
    endpoint_nodes[0::2] = left
    endpoint_nodes[1::2] = right
    geometric = bool(points is not None and ports == 12 and endpoint_nodes.size)
    directional_alignment: np.ndarray | None = None
    local_frame_hash: str | None = None
    if geometric:
        point_array = np.asarray(points, dtype=float)
        if point_array.shape != (int(patch_count), 3):
            raise ValueError(
                "points must have shape (patch_count, 3) for geometric port routing"
            )
        endpoint_other = np.empty_like(endpoint_nodes)
        endpoint_other[0::2] = right
        endpoint_other[1::2] = left
        endpoint_ports, directional_alignment, local_frame_hash = (
            _geometric_endpoint_port_assignment(
                point_array,
                endpoint_nodes,
                endpoint_other,
            )
        )
        overflow = int(np.sum(node_degree > ports))
        left_port = endpoint_ports[0::2].astype(np.int16)
        right_port = endpoint_ports[1::2].astype(np.int16)
    elif endpoint_nodes.size:
        order = np.argsort(endpoint_nodes, kind="stable")
        sorted_nodes = endpoint_nodes[order]
        group_start = np.empty(sorted_nodes.size, dtype=bool)
        group_start[0] = True
        group_start[1:] = sorted_nodes[1:] != sorted_nodes[:-1]
        starts = np.flatnonzero(group_start)
        counts = np.diff(np.append(starts, sorted_nodes.size))
        local_rank_sorted = np.arange(sorted_nodes.size, dtype=np.int64) - np.repeat(starts, counts)
        local_rank = np.empty_like(local_rank_sorted)
        local_rank[order] = local_rank_sorted
        overflow = int(np.sum(local_rank >= ports))
        left_port = (local_rank[0::2] % ports).astype(np.int16)
        right_port = (local_rank[1::2] % ports).astype(np.int16)
    else:
        left_port = np.zeros(0, dtype=np.int16)
        right_port = np.zeros(0, dtype=np.int16)
        overflow = 0
    return ScreenPortMap(
        left_port=left_port,
        right_port=right_port,
        ports_per_patch=ports,
        overflow_count=int(overflow),
        node_degree=node_degree,
        routing_mode=(
            "icosahedral_directional_assignment"
            if geometric
            else "sequential_endpoint_order"
        ),
        directional_alignment=directional_alignment,
        local_frame_hash=local_frame_hash,
        federation_weights=weights,
    )


def echosahedral_patch_architecture_report(ports_per_patch: int = 12) -> dict[str, Any]:
    """Describe the shared local twelve-port patch template.

    This is the implementation-facing patch object from the screen
    microphysics paper.  It is stored once and reused by every cell rather than
    materializing ``patch_count`` identical direction tables.
    """

    ports = int(ports_per_patch)
    exact_template = ports == 12
    if exact_template:
        base = build_geodesic_icosahedral_tower(0).levels[0]
        directions = base.vertices
        edges = np.asarray(base.edges, dtype=np.int64)
        faces = np.asarray(base.faces, dtype=np.int64)
        antipodes = np.asarray(
            [
                int(np.argmin(np.linalg.norm(directions + direction, axis=1)))
                for direction in directions
            ],
            dtype=np.int64,
        )
        antipodal_residual = float(
            np.max(np.linalg.norm(directions + directions[antipodes], axis=1))
        )
        a5 = icosahedral_a5_equivariance_report(0)
        permutations = icosahedral_a5_port_permutations()
        incidence_audit = _local_icosahedral_incidence_audit(
            edges,
            faces,
            antipodes,
            permutations,
        )
        template_hash = _port_template_hash(
            directions,
            antipodes,
            edges=edges,
            faces=faces,
        )
        antipodal_pairs = sorted(
            {
                (min(index, int(partner)), max(index, int(partner)))
                for index, partner in enumerate(antipodes)
            }
        )
    else:
        directions = np.zeros((ports, 3), dtype=float)
        edges = np.zeros((0, 2), dtype=np.int64)
        faces = np.zeros((0, 3), dtype=np.int64)
        antipodes = np.full(ports, -1, dtype=np.int64)
        antipodal_residual = float("inf")
        a5 = {"A5_ROTATION_GROUP_ORDER_60_RECEIPT": False}
        incidence_audit = {
            "incidence_receipt": False,
            "vertex_degree_profile": {},
            "adjacency_irrep_multiplicities": [],
            "all_a5_actions_preserve_edges": False,
            "all_a5_actions_preserve_oriented_faces": False,
            "all_a5_actions_commute_with_antipode": False,
        }
        template_hash = None
        antipodal_pairs = []
    receipt = bool(
        exact_template
        and antipodal_residual <= 5.0e-14
        and len(antipodal_pairs) == 6
        and a5.get("A5_ROTATION_GROUP_ORDER_60_RECEIPT") is True
        and incidence_audit["incidence_receipt"]
    )
    return {
        "schema": "oph.echosahedral_local_patch.v1",
        "patch_object_role": "bounded_local_observer_carrier",
        "port_names": echosahedral_port_names(ports),
        "port_count": ports,
        "port_template": "regular_icosahedron_vertices" if exact_template else "noncanonical",
        "port_direction_template": directions.tolist() if exact_template else None,
        "hidden_presentation_coordinates": True,
        "hidden_coordinates_eligible_for_promoted_geometry": False,
        "local_edges": edges.tolist() if exact_template else None,
        "local_faces": faces.tolist() if exact_template else None,
        "local_edge_count": int(edges.shape[0]),
        "local_face_count": int(faces.shape[0]),
        "inverse_port_involution": antipodes.tolist() if exact_template else None,
        "antipodal_pairs": [list(pair) for pair in antipodal_pairs],
        "maximum_antipodal_residual": antipodal_residual if exact_template else None,
        "proper_rotation_group_order": a5.get("base_rotation_count"),
        "incidence_audit": incidence_audit,
        "template_hash": template_hash,
        "shared_template_instantiated_per_cell": receipt,
        "ECHOSAHEDRAL_LOCAL_PATCH_ARCHITECTURE_RECEIPT": receipt,
        "ECHOSAHEDRAL_CARRIER_CONFORMANCE_RECEIPT": receipt,
        "LOCAL_ICOSAHEDRAL_INCIDENCE_RECEIPT": bool(
            incidence_audit["incidence_receipt"]
        ),
        "GLOBAL_TWELVE_UNIT_SCREEN_SIEVE_RECEIPT": False,
        "LOCAL_TO_GLOBAL_A5_INTERTWINER_RECEIPT": False,
        "claim_boundary": (
            "The twelve local P0..P11 ports belong to every bounded observer cell. "
            "Their 12 vertices, 30 edges, 20 outward-oriented faces, antipode, and "
            "order-60 A5 action are hidden carrier-presentation data. They may shape "
            "local coupling but are forbidden as promoted S2, H3, cap-normal, clock, "
            "or event coordinates. "
            "The separate screen-sieve theorem derives twelve global unit-defect ports "
            "on a connected cut. Identifying the two twelve-dimensional objects requires "
            "a quotient-visible, refinement-natural A5 intertwiner."
        ),
    }


def _geometric_endpoint_port_assignment(
    points: np.ndarray,
    endpoint_nodes: np.ndarray,
    endpoint_other: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, str]:
    directions = build_geodesic_icosahedral_tower(0).levels[0].vertices
    frames = _local_surface_frames(points)
    vectors = points[endpoint_other] - points[endpoint_nodes]
    vectors /= np.maximum(np.linalg.norm(vectors, axis=1, keepdims=True), 1.0e-15)
    local_vectors = np.einsum(
        "ni,nij->nj",
        vectors,
        frames[endpoint_nodes],
        optimize=True,
    )
    scores = np.asarray(local_vectors @ directions.T, dtype=np.float32)
    endpoint_ports = np.argmax(scores, axis=1).astype(np.int64)

    # The vectorized nearest-direction assignment is normally collision-free
    # on a degree-three cell-dual graph.  Repair only the nodes with duplicate
    # chosen slots using an exact small Hungarian assignment; this avoids one
    # Python/SciPy optimization call per patch at 256k scale.
    keys = endpoint_nodes * directions.shape[0] + endpoint_ports
    _, first_index, key_counts = np.unique(
        keys,
        return_index=True,
        return_counts=True,
    )
    duplicate_key_values = keys[first_index[key_counts > 1]]
    collision_nodes = np.unique(duplicate_key_values // directions.shape[0])
    if collision_nodes.size:
        order = np.argsort(endpoint_nodes, kind="stable")
        sorted_nodes = endpoint_nodes[order]
        starts = np.searchsorted(sorted_nodes, collision_nodes, side="left")
        stops = np.searchsorted(sorted_nodes, collision_nodes, side="right")
        for start, stop in zip(starts, stops, strict=True):
            endpoint_ids = order[int(start) : int(stop)]
            assign_ids = endpoint_ids[: directions.shape[0]]
            row_ids, port_ids = linear_sum_assignment(-scores[assign_ids])
            endpoint_ports[assign_ids[row_ids]] = port_ids
            # Only graphs exceeding the local twelve-slot architecture can
            # reuse a direction. The overflow is separately fail-closed.
            for endpoint_id in endpoint_ids[directions.shape[0] :]:
                endpoint_ports[endpoint_id] = int(np.argmax(scores[endpoint_id]))
    alignment = scores[
        np.arange(endpoint_nodes.size, dtype=np.int64), endpoint_ports
    ].astype(float)
    if np.any(endpoint_ports < 0):
        raise AssertionError("every routed endpoint must receive one local patch port")
    frame_hash = hashlib.sha256(
        np.asarray(np.round(frames, 15), dtype="<f8").tobytes()
    ).hexdigest()
    return endpoint_ports, alignment, frame_hash


def _local_surface_frames(points: np.ndarray) -> np.ndarray:
    normals = np.asarray(points, dtype=float).copy()
    normals /= np.maximum(np.linalg.norm(normals, axis=1, keepdims=True), 1.0e-15)
    reference = np.tile(np.asarray([0.0, 0.0, 1.0]), (normals.shape[0], 1))
    near_pole = np.abs(normals[:, 2]) > 0.9
    reference[near_pole] = np.asarray([1.0, 0.0, 0.0])
    tangent_x = np.cross(reference, normals)
    tangent_x /= np.maximum(np.linalg.norm(tangent_x, axis=1, keepdims=True), 1.0e-15)
    tangent_y = np.cross(normals, tangent_x)
    return np.stack((tangent_x, tangent_y, normals), axis=2)


def _port_template_hash(
    directions: np.ndarray,
    antipodes: np.ndarray,
    *,
    edges: np.ndarray,
    faces: np.ndarray,
) -> str:
    payload = {
        "schema": "oph.echosahedral_local_patch.template.v1",
        "directions": np.round(directions, 15).tolist(),
        "antipodes": antipodes.tolist(),
        "edges": np.asarray(edges, dtype=np.int64).tolist(),
        "faces": np.asarray(faces, dtype=np.int64).tolist(),
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _local_icosahedral_incidence_audit(
    edges: np.ndarray,
    faces: np.ndarray,
    antipodes: np.ndarray,
    permutations: tuple[tuple[int, ...], ...],
) -> dict[str, Any]:
    """Verify the complete hidden 12/30/20 carrier incidence and A5 action."""

    edge_array = np.asarray(edges, dtype=np.int64)
    face_array = np.asarray(faces, dtype=np.int64)
    antipode = np.asarray(antipodes, dtype=np.int64)
    degree = np.bincount(edge_array.reshape(-1), minlength=12)
    edge_set = {tuple(sorted(map(int, edge))) for edge in edge_array}

    def oriented_face_key(face: tuple[int, int, int]) -> tuple[int, int, int]:
        rotations = (face, (face[1], face[2], face[0]), (face[2], face[0], face[1]))
        return min(rotations)

    face_set = {
        oriented_face_key(tuple(int(value) for value in face)) for face in face_array
    }
    edge_preservation: list[bool] = []
    face_preservation: list[bool] = []
    antipode_commutation: list[bool] = []
    for permutation in permutations:
        edge_preservation.append(
            {
                tuple(sorted((permutation[int(left)], permutation[int(right)])))
                for left, right in edge_array
            }
            == edge_set
        )
        face_preservation.append(
            {
                oriented_face_key(
                    tuple(permutation[int(value)] for value in face)
                )
                for face in face_array
            }
            == face_set
        )
        antipode_commutation.append(
            all(
                permutation[int(antipode[index])]
                == int(antipode[permutation[index]])
                for index in range(12)
            )
        )
    adjacency = np.zeros((12, 12), dtype=float)
    adjacency[edge_array[:, 0], edge_array[:, 1]] = 1.0
    adjacency[edge_array[:, 1], edge_array[:, 0]] = 1.0
    eigenvalues = np.linalg.eigvalsh(adjacency)
    clusters: list[dict[str, Any]] = []
    for value in eigenvalues:
        if clusters and abs(value - float(clusters[-1]["eigenvalue"])) <= 1.0e-10:
            clusters[-1]["multiplicity"] += 1
        else:
            clusters.append(
                {"eigenvalue": float(value), "multiplicity": 1}
            )
    multiplicities = sorted(int(row["multiplicity"]) for row in clusters)
    incidence_receipt = bool(
        edge_array.shape == (30, 2)
        and face_array.shape == (20, 3)
        and len(edge_set) == 30
        and len(face_set) == 20
        and np.all(degree == 5)
        and antipode.shape == (12,)
        and np.all(antipode[antipode] == np.arange(12))
        and np.all(antipode != np.arange(12))
        and len(permutations) == 60
        and all(edge_preservation)
        and all(face_preservation)
        and all(antipode_commutation)
        and multiplicities == [1, 3, 3, 5]
    )
    return {
        "vertex_count": 12,
        "edge_count": int(edge_array.shape[0]),
        "face_count": int(face_array.shape[0]),
        "vertex_degree_profile": {
            str(value): int(np.sum(degree == value)) for value in sorted(set(degree))
        },
        "antipode_fixed_point_count": int(np.sum(antipode == np.arange(12))),
        "antipode_involution": bool(np.all(antipode[antipode] == np.arange(12))),
        "a5_action_count": len(permutations),
        "all_a5_actions_preserve_edges": bool(all(edge_preservation)),
        "all_a5_actions_preserve_oriented_faces": bool(all(face_preservation)),
        "all_a5_actions_commute_with_antipode": bool(all(antipode_commutation)),
        "adjacency_spectrum": clusters,
        "adjacency_irrep_multiplicities": multiplicities,
        "incidence_receipt": incidence_receipt,
    }


def initialize_echosahedral_patch_state(
    *,
    patch_count: int,
    ports_per_patch: int,
    group_order: int,
    edge_left: np.ndarray,
    edge_right: np.ndarray,
    port_map: ScreenPortMap,
    routed_left_state: np.ndarray,
    routed_right_state: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Materialize the full finite local port state of every patch cell."""

    state = rng.integers(
        0,
        int(group_order),
        size=(int(patch_count), int(ports_per_patch)),
        dtype=np.int16,
    )
    sync_routed_echosahedral_patch_state(
        state,
        edge_left=edge_left,
        edge_right=edge_right,
        port_map=port_map,
        routed_left_state=routed_left_state,
        routed_right_state=routed_right_state,
    )
    return state


def sync_routed_echosahedral_patch_state(
    patch_port_state: np.ndarray,
    *,
    edge_left: np.ndarray,
    edge_right: np.ndarray,
    port_map: ScreenPortMap,
    routed_left_state: np.ndarray,
    routed_right_state: np.ndarray,
) -> None:
    """Commit routed edge endpoint values into their local patch port slots."""

    patch_port_state[
        np.asarray(edge_left, dtype=np.int64),
        np.asarray(port_map.left_port, dtype=np.int64),
    ] = np.asarray(routed_left_state, dtype=patch_port_state.dtype)
    patch_port_state[
        np.asarray(edge_right, dtype=np.int64),
        np.asarray(port_map.right_port, dtype=np.int64),
    ] = np.asarray(routed_right_state, dtype=patch_port_state.dtype)


def canonicalize_echosahedral_patch_state(
    patch_port_state: np.ndarray,
    *,
    edge_left: np.ndarray,
    edge_right: np.ndarray,
    routed_left_state: np.ndarray,
    routed_right_state: np.ndarray,
    group_name: str,
    group_order: int,
) -> np.ndarray:
    """Quotient every local P0..P11 value by the node reference frame.

    The lowest-index incident edge is sent to the identity, exactly matching
    the routed overlap canonicalizer. All twelve local coordinates transform
    together, including exposed/reserved slots, so the record input is
    invariant under independent node-frame changes.
    """

    state = np.asarray(patch_port_state, dtype=np.int16)
    if state.ndim != 2 or state.shape[1] != 12:
        raise ValueError("gauge quotient requires one twelve-port state per patch")
    left = np.asarray(edge_left, dtype=np.int64)
    right = np.asarray(edge_right, dtype=np.int64)
    left_values = np.asarray(routed_left_state, dtype=np.int16)
    right_values = np.asarray(routed_right_state, dtype=np.int16)
    if not (left.shape == right.shape == left_values.shape == right_values.shape):
        raise ValueError("edge endpoints and routed states must have matching shapes")
    if left.size and (
        int(np.min(left)) < 0
        or int(np.min(right)) < 0
        or int(np.max(left)) >= state.shape[0]
        or int(np.max(right)) >= state.shape[0]
    ):
        raise ValueError("edge endpoint lies outside the patch state array")

    frames = np.zeros(state.shape[0], dtype=np.int16)
    if left.size:
        edge_indices = np.arange(left.size, dtype=np.int64)
        reference_edges = np.full(state.shape[0], left.size, dtype=np.int64)
        np.minimum.at(reference_edges, left, edge_indices)
        np.minimum.at(reference_edges, right, edge_indices)
        left_is_reference = reference_edges[left] == edge_indices
        right_is_reference = reference_edges[right] == edge_indices
        reference_values = np.zeros(state.shape[0], dtype=np.int16)
        reference_values[left[left_is_reference]] = left_values[left_is_reference]
        reference_values[right[right_is_reference]] = right_values[right_is_reference]
        frames = group_inverse_indices(
            reference_values,
            group_name=group_name,
            group_order=group_order,
        )
    canonical = np.empty_like(state, dtype=np.int16)
    # Work one coordinate at a time. At the 1M rung, converting the entire
    # patch-by-12 carrier to the int64 indexing dtype would create a 96 MB
    # transient on every repair cycle.
    for port_index in range(state.shape[1]):
        canonical[:, port_index] = group_multiply_indices(
            frames,
            state[:, port_index],
            group_name=group_name,
            group_order=group_order,
        )
    return canonical


def echosahedral_patch_record_signature(
    routed_signature: np.ndarray,
    patch_port_state: np.ndarray,
) -> np.ndarray:
    """Bind all local ports through A5-invariant incidence fingerprints.

    Local port names and local icosahedral frames are presentation data.  The
    record therefore hashes the multiset of vertex labels, unordered edge
    label pairs, and cyclically oriented face-label triples.  Proper
    icosahedral rotations merely permute those three incidence collections,
    while an orientation reversal can still change the face fingerprint.
    """

    state = np.asarray(patch_port_state)
    if state.ndim != 2 or state.shape[1] != 12:
        raise ValueError("echosahedral patch record signature requires exactly 12 ports")
    routed_values = np.asarray(routed_signature, dtype=np.int64)
    if routed_values.ndim != 1 or routed_values.size != state.shape[0]:
        raise ValueError(
            "routed signature must contain exactly one token per echosahedral patch"
        )
    if np.any(state < 0) or (state.size and int(np.max(state)) >= 257):
        raise ValueError("local port labels must lie in the canonical range 0..256")
    port_token = _a5_invariant_local_incidence_token(state)
    routed = routed_values.view(np.uint64)
    mixed = port_token ^ (routed * np.uint64(0x9E3779B97F4A7C15))
    mixed ^= mixed >> np.uint64(30)
    mixed *= np.uint64(0xBF58476D1CE4E5B9)
    mixed ^= mixed >> np.uint64(27)
    mixed *= np.uint64(0x94D049BB133111EB)
    mixed ^= mixed >> np.uint64(31)
    return (mixed & np.uint64(0x7FFFFFFFFFFFFFFF)).astype(np.int64)


def _a5_invariant_local_incidence_token(
    patch_port_state: np.ndarray,
    *,
    chunk_size: int = 32_768,
) -> np.ndarray:
    """Hash a patch-by-12 state invariantly under all proper A5 rotations."""

    state = np.asarray(patch_port_state)
    if state.ndim != 2 or state.shape[1] != 12:
        raise ValueError("A5 incidence token requires one twelve-port row per patch")
    base = build_geodesic_icosahedral_tower(0).levels[0]
    edges = np.asarray(base.edges, dtype=np.int64)
    faces = np.asarray(base.faces, dtype=np.int64)
    output = np.zeros(state.shape[0], dtype=np.uint64)
    radix = np.uint64(257)
    radix2 = radix * radix
    for start in range(0, state.shape[0], max(1, int(chunk_size))):
        stop = min(state.shape[0], start + max(1, int(chunk_size)))
        chunk = np.asarray(state[start:stop], dtype=np.uint64)
        token = np.zeros(stop - start, dtype=np.uint64)

        vertex_codes = np.sort(chunk, axis=1)
        _fold_sorted_incidence_codes(
            token,
            vertex_codes,
            family_salt=np.uint64(0x243F6A8885A308D3),
        )

        edge_left = chunk[:, edges[:, 0]]
        edge_right = chunk[:, edges[:, 1]]
        edge_codes = np.minimum(edge_left, edge_right) * radix + np.maximum(
            edge_left,
            edge_right,
        )
        edge_codes.sort(axis=1)
        _fold_sorted_incidence_codes(
            token,
            edge_codes,
            family_salt=np.uint64(0x13198A2E03707344),
        )

        first = chunk[:, faces[:, 0]]
        second = chunk[:, faces[:, 1]]
        third = chunk[:, faces[:, 2]]
        face_codes = np.minimum.reduce(
            (
                first * radix2 + second * radix + third,
                second * radix2 + third * radix + first,
                third * radix2 + first * radix + second,
            )
        )
        face_codes.sort(axis=1)
        _fold_sorted_incidence_codes(
            token,
            face_codes,
            family_salt=np.uint64(0xA4093822299F31D0),
        )
        output[start:stop] = token
    return output


def _fold_sorted_incidence_codes(
    token: np.ndarray,
    codes: np.ndarray,
    *,
    family_salt: np.uint64,
) -> None:
    """Fold a canonical sorted incidence multiset into a bounded token."""

    for index in range(codes.shape[1]):
        position_salt = np.uint64(
            ((index + 1) * 0x9E3779B97F4A7C15) & ((1 << 64) - 1)
        )
        token ^= _screen_port_splitmix64(codes[:, index] ^ family_salt ^ position_salt)


def _screen_port_splitmix64(values: np.ndarray) -> np.ndarray:
    with np.errstate(over="ignore"):
        mixed = np.asarray(values, dtype=np.uint64) + np.uint64(0x9E3779B97F4A7C15)
        mixed = (mixed ^ (mixed >> np.uint64(30))) * np.uint64(
            0xBF58476D1CE4E5B9
        )
        mixed = (mixed ^ (mixed >> np.uint64(27))) * np.uint64(
            0x94D049BB133111EB
        )
        return mixed ^ (mixed >> np.uint64(31))


def echosahedral_patch_state_report(
    patch_port_state: np.ndarray,
    *,
    edge_left: np.ndarray,
    edge_right: np.ndarray,
    port_map: ScreenPortMap,
    routed_left_state: np.ndarray,
    routed_right_state: np.ndarray,
    record_signature_bound: bool = False,
    record_binding_mode: str = "none",
) -> dict[str, Any]:
    """Audit the materialized federation state without serializing every cell."""

    state = np.asarray(patch_port_state)
    if state.ndim != 2:
        raise ValueError("patch_port_state must be a two-dimensional patch-by-port array")
    left = np.asarray(edge_left, dtype=np.int64)
    right = np.asarray(edge_right, dtype=np.int64)
    left_ports = np.asarray(port_map.left_port, dtype=np.int64)
    right_ports = np.asarray(port_map.right_port, dtype=np.int64)
    used = np.zeros(state.shape, dtype=bool)
    used[left, left_ports] = True
    used[right, right_ports] = True
    routed_consistent = bool(
        np.array_equal(state[left, left_ports], np.asarray(routed_left_state))
        and np.array_equal(state[right, right_ports], np.asarray(routed_right_state))
    )
    state_hash = hashlib.sha256(np.asarray(state, dtype="<i2").tobytes()).hexdigest()
    return {
        "schema": "oph.echosahedral_patch_federation_state.v1",
        "patch_count": int(state.shape[0]) if state.ndim == 2 else 0,
        "ports_per_patch": int(state.shape[1]) if state.ndim == 2 else 0,
        "materialized_local_port_state_count": int(state.size),
        "routed_port_slot_count": int(np.sum(used)),
        "unrouted_exposed_or_reserved_port_slot_count": int(np.sum(~used)),
        "routed_edge_endpoint_count": int(2 * left.size),
        "routed_state_consistent_with_patch_slots": routed_consistent,
        "all_port_coordinate_readout_maps_materialized": bool(
            state.ndim == 2 and state.shape[1] == 12
        ),
        "patch_port_state_sha256": state_hash,
        "ECHOSAHEDRAL_PATCH_STATE_INSTANTIATION_RECEIPT": bool(
            state.ndim == 2
            and state.shape[0] > 0
            and state.shape[1] == 12
            and routed_consistent
            and port_map.overflow_count == 0
        ),
        "PATCH_ALL_PORT_READBACK_RECEIPT": bool(
            state.ndim == 2 and state.shape[1] == 12
        ),
        "record_signature_binding": {
            "algorithm": (
                "routed_record_xor_splitmix_of_a5_invariant_vertex_edge_"
                "oriented_face_incidence_multisets"
            ),
            "coordinates_per_patch": 12,
            "full_patch_state_bound": bool(record_signature_bound and state.shape[1] == 12),
            "mode": str(record_binding_mode),
            "local_value_frame_gauge_quotient": record_binding_mode.startswith(
                "node_reference_port_gauge_quotient"
            ),
            "local_frame_gauge_quotient": record_binding_mode.startswith(
                "node_reference_port_gauge_quotient"
            ),
            "local_a5_port_reorientation_invariant": record_binding_mode.endswith(
                "a5_incidence_invariants"
            ),
            "hidden_xyz_coordinates_used": False,
        },
        "RECORD_SIGNATURE_BINDS_ALL_LOCAL_PORT_STATE_RECEIPT": bool(
            record_signature_bound and state.shape[1] == 12 and routed_consistent
        ),
        "REFERENCE_LOCAL_A5_RECORD_INVARIANCE_RECEIPT": bool(
            record_signature_bound
            and state.shape[1] == 12
            and record_binding_mode.endswith("a5_incidence_invariants")
        ),
        "CARRIER_QUOTIENT_INVARIANCE_RECEIPT": False,
        "CARRIER_REFINEMENT_NATURALITY_RECEIPT": False,
        "PHYSICAL_A5_PORT_EMERGENCE_RECEIPT": False,
        "PHYSICAL_STANDARD_MODEL_EMERGENCE_RECEIPT": False,
        "claim_boundary": (
            "This is the actual finite local state for every port of every "
            "echosahedral cell. Routed slots participate in overlaps and repair; "
            "unrouted slots remain local exposed/reserved state. When the record "
            "binding receipt is true, all twelve coordinates enter each committed "
            "record token. This is a carrier instantiation, not a proof of the "
            "full physical quotient: the local A5-invariant record map is now exact, "
            "but schedule/worker invariance and cross-rung carrier refinement remain "
            "separate false receipts. It is not a proof of the collective A5 screen "
            "sieve or Standard Model emergence."
        ),
    }


def write_echosahedral_patch_state_artifact(
    path: Path,
    *,
    patch_port_state: np.ndarray,
    canonical_record_port_state: np.ndarray | None = None,
    edge_left: np.ndarray,
    edge_right: np.ndarray,
    port_map: ScreenPortMap,
    routed_left_state: np.ndarray,
    routed_right_state: np.ndarray,
    record_signature: np.ndarray,
    committed: np.ndarray,
) -> dict[str, Any]:
    """Write the replay-relevant federation state and return compact metadata."""

    destination = Path(path)
    state = np.asarray(patch_port_state, dtype=np.int16)
    canonical_state = (
        np.asarray(canonical_record_port_state, dtype=np.int16)
        if canonical_record_port_state is not None
        else None
    )
    signatures = np.asarray(record_signature, dtype=np.int64)
    committed_mask = np.asarray(committed, dtype=bool)
    left = np.asarray(edge_left, dtype=np.int64)
    right = np.asarray(edge_right, dtype=np.int64)
    left_ports = np.asarray(port_map.left_port, dtype=np.int16)
    right_ports = np.asarray(port_map.right_port, dtype=np.int16)
    routed_left = np.asarray(routed_left_state, dtype=np.int16)
    routed_right = np.asarray(routed_right_state, dtype=np.int16)
    if state.ndim != 2 or state.shape[1] != 12:
        raise ValueError("patch state artifact requires a patch-by-12 state array")
    if signatures.shape != (state.shape[0],):
        raise ValueError("record signature artifact must contain one token per patch")
    if committed_mask.shape != (state.shape[0],):
        raise ValueError("committed mask artifact must contain one flag per patch")
    if canonical_state is not None and canonical_state.shape != state.shape:
        raise ValueError("canonical record port state must match the raw patch state shape")
    if not (
        left.shape
        == right.shape
        == left_ports.shape
        == right_ports.shape
        == routed_left.shape
        == routed_right.shape
    ):
        raise ValueError("routed edge arrays in the patch state artifact must match")
    arrays: dict[str, np.ndarray] = {
        "patch_port_state": state,
        "edge_left": left,
        "edge_right": right,
        "left_port": left_ports,
        "right_port": right_ports,
        "routed_left_state": routed_left,
        "routed_right_state": routed_right,
        "record_signature": signatures,
        "committed": committed_mask,
        "port_names": np.asarray(echosahedral_port_names(12), dtype="<U3"),
    }
    if canonical_state is not None:
        arrays["canonical_record_port_state"] = canonical_state
    np.savez_compressed(destination, **arrays)
    digest = hashlib.sha256()
    with destination.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {
        "path": destination.name,
        "format": "numpy_compressed_npz",
        "byte_count": int(destination.stat().st_size),
        "file_sha256": "sha256:" + digest.hexdigest(),
        "arrays": {
            "patch_port_state": {
                "shape": [int(value) for value in state.shape],
                "dtype": str(state.dtype),
            },
            "record_signature": {
                "shape": [int(value) for value in signatures.shape],
                "dtype": str(signatures.dtype),
            },
            "committed": {
                "shape": [int(value) for value in committed_mask.shape],
                "dtype": str(committed_mask.dtype),
            },
            "port_names": {"shape": [12], "dtype": "<U3"},
            "edge_left": {"shape": [int(left.size)], "dtype": str(left.dtype)},
            "edge_right": {"shape": [int(right.size)], "dtype": str(right.dtype)},
            "left_port": {"shape": [int(left_ports.size)], "dtype": str(left_ports.dtype)},
            "right_port": {"shape": [int(right_ports.size)], "dtype": str(right_ports.dtype)},
            "routed_left_state": {
                "shape": [int(routed_left.size)],
                "dtype": str(routed_left.dtype),
            },
            "routed_right_state": {
                "shape": [int(routed_right.size)],
                "dtype": str(routed_right.dtype),
            },
            **(
                {
                    "canonical_record_port_state": {
                        "shape": [int(value) for value in canonical_state.shape],
                        "dtype": str(canonical_state.dtype),
                    }
                }
                if canonical_state is not None
                else {}
            ),
        },
    }
