from __future__ import annotations

import hashlib
from typing import Any

import numpy as np

from oph_fpe.finite_groups import S3_INV, S3_MUL


GAUGE_COVARIANT_OVERLAP_SCHEMA = "left_equals_g_ij_times_right_v1"
GAUGE_QUOTIENT_CANONICALIZER = "node_reference_port_frame_v1"


def group_multiply_indices(
    left: np.ndarray | int,
    right: np.ndarray | int,
    *,
    group_name: str,
    group_order: int,
) -> np.ndarray:
    """Multiply finite-group element indices without treating S3 as a clock."""

    name, order = _normalized_group(group_name, group_order)
    left_values = np.asarray(left, dtype=np.int64)
    right_values = np.asarray(right, dtype=np.int64)
    if name == "S3":
        return S3_MUL[left_values, right_values].astype(np.int16, copy=False)
    return ((left_values + right_values) % order).astype(np.int16, copy=False)


def group_inverse_indices(
    values: np.ndarray | int,
    *,
    group_name: str,
    group_order: int,
) -> np.ndarray:
    """Invert finite-group element indices without assuming commutativity."""

    name, order = _normalized_group(group_name, group_order)
    elements = np.asarray(values, dtype=np.int64)
    if name == "S3":
        return S3_INV[elements].astype(np.int16, copy=False)
    return ((-elements) % order).astype(np.int16, copy=False)


def transport_right_to_left(
    port_right: np.ndarray,
    gauge: np.ndarray,
    *,
    group_name: str,
    group_order: int,
) -> np.ndarray:
    """Transport right-frame port labels through ``g_ij`` into the left frame.

    The convention is ``g_ij: frame j -> frame i``.  Consequently a consistent
    edge obeys ``port_left = g_ij * port_right``.  Under local frame changes,
    ``p_i -> h_i p_i`` and ``g_ij -> h_i g_ij h_j^-1``, so this equality is
    presentation independent.
    """

    right_values, gauge_values = _matching_arrays(port_right, gauge, "right ports", "gauge")
    return group_multiply_indices(
        gauge_values,
        right_values,
        group_name=group_name,
        group_order=group_order,
    )


def covariant_mismatch_mask(
    port_left: np.ndarray,
    port_right: np.ndarray,
    gauge: np.ndarray,
    *,
    group_name: str,
    group_order: int,
) -> np.ndarray:
    """Return the exact gauge-covariant mismatch predicate on every edge."""

    left_values, right_values = _matching_arrays(port_left, port_right, "left ports", "right ports")
    _, gauge_values = _matching_arrays(left_values, gauge, "ports", "gauge")
    transported = transport_right_to_left(
        right_values,
        gauge_values,
        group_name=group_name,
        group_order=group_order,
    )
    return left_values != transported


def covariant_discrepancy(
    port_left: np.ndarray,
    port_right: np.ndarray,
    gauge: np.ndarray,
    *,
    group_name: str,
    group_order: int,
) -> np.ndarray:
    """Return ``port_left * (g_ij * port_right)^-1`` in the left frame."""

    left_values, right_values = _matching_arrays(port_left, port_right, "left ports", "right ports")
    _, gauge_values = _matching_arrays(left_values, gauge, "ports", "gauge")
    transported = transport_right_to_left(
        right_values,
        gauge_values,
        group_name=group_name,
        group_order=group_order,
    )
    return group_multiply_indices(
        left_values,
        group_inverse_indices(transported, group_name=group_name, group_order=group_order),
        group_name=group_name,
        group_order=group_order,
    )


def repair_covariant_port_pairs(
    port_left: np.ndarray,
    port_right: np.ndarray,
    gauge: np.ndarray,
    edges: np.ndarray,
    repair_left: np.ndarray,
    *,
    group_name: str,
    group_order: int,
) -> None:
    """Make selected edge slots covariantly consistent while holding links fixed.

    ``repair_left`` chooses the endpoint that is rewritten.  Both branches are
    gauge covariant: either ``p_i <- g_ij p_j`` or
    ``p_j <- g_ij^-1 p_i``.
    """

    left_values, right_values = _matching_arrays(port_left, port_right, "left ports", "right ports")
    _, gauge_values = _matching_arrays(left_values, gauge, "ports", "gauge")
    selected = np.asarray(edges, dtype=np.int64)
    directions = np.asarray(repair_left, dtype=bool)
    if selected.ndim != 1 or directions.ndim != 1 or selected.shape != directions.shape:
        raise ValueError("selected edges and repair directions must be matching one-dimensional arrays")
    if selected.size == 0:
        return
    left_edges = selected[directions]
    if left_edges.size:
        left_values[left_edges] = transport_right_to_left(
            right_values[left_edges],
            gauge_values[left_edges],
            group_name=group_name,
            group_order=group_order,
        )
    right_edges = selected[~directions]
    if right_edges.size:
        inverse_gauge = group_inverse_indices(
            gauge_values[right_edges],
            group_name=group_name,
            group_order=group_order,
        )
        right_values[right_edges] = group_multiply_indices(
            inverse_gauge,
            left_values[right_edges],
            group_name=group_name,
            group_order=group_order,
        )


def absorb_discrepancy_into_gauge(
    gauge: np.ndarray,
    edges: np.ndarray,
    discrepancy: np.ndarray,
    *,
    group_name: str,
    group_order: int,
) -> None:
    """Apply ``g_ij <- delta_i * g_ij`` for selected left-frame deltas.

    For ``delta_i = p_i (g_ij p_j)^-1`` this gives
    ``g_ij' = p_i p_j^-1`` exactly, so the link absorbs the overlap defect.
    Left multiplication is required in the non-Abelian case and preserves the
    link transformation law.
    """

    gauge_values = np.asarray(gauge)
    selected = np.asarray(edges, dtype=np.int64)
    deltas = np.asarray(discrepancy, dtype=np.int16)
    if selected.ndim != 1 or deltas.ndim != 1 or selected.shape != deltas.shape:
        raise ValueError("selected edges and discrepancy must be matching one-dimensional arrays")
    if selected.size:
        gauge_values[selected] = group_multiply_indices(
            deltas,
            gauge_values[selected],
            group_name=group_name,
            group_order=group_order,
        )


def repair_production_sector_links(
    gauge: np.ndarray,
    chosen_edges: np.ndarray,
    chosen_discrepancy: np.ndarray,
    *,
    group_name: str,
    group_order: int,
    rng: np.random.Generator,
    config: dict[str, Any],
) -> int:
    """Replay the production BW sector-link mutation exactly.

    The production move samples a configured subset of already-selected
    overlap repairs, discards identity discrepancies, and mutates the link
    before the endpoint repair is applied.  Keeping this primitive outside the
    production runner lets perturb/resettle probes exercise the identical move
    instead of a look-alike implementation.

    ``repair_coupled_group_compose`` is gauge covariant because the
    left-frame discrepancy multiplies ``g_ij`` on the left.  The legacy
    ``cool_to_identity`` branch is retained so replay is exact, but callers
    must not issue a gauge-covariance receipt for that branch.
    """

    selected = np.asarray(chosen_edges, dtype=np.int64)
    discrepancy = np.asarray(chosen_discrepancy, dtype=np.int16)
    if selected.ndim != 1 or discrepancy.ndim != 1 or selected.shape != discrepancy.shape:
        raise ValueError("chosen edges and discrepancies must be matching one-dimensional arrays")
    if (
        str(group_name).upper() != "S3"
        or int(group_order) != 6
        or not bool(config.get("enabled", False))
        or selected.size == 0
    ):
        return 0
    probability = float(config.get("probability", 0.0))
    if not np.isfinite(probability) or probability <= 0.0:
        return 0
    mask = rng.random(selected.size) < min(max(probability, 0.0), 1.0)
    if not np.any(mask):
        return 0
    edges = selected[mask]
    proposals = discrepancy[mask].astype(np.int64, copy=False) % 6
    nontrivial = proposals != 0
    if not np.any(nontrivial):
        return 0
    edges = edges[nontrivial]
    proposals = proposals[nontrivial]
    if str(config.get("mode", "repair_coupled_group_compose")) == "cool_to_identity":
        np.asarray(gauge)[edges] = 0
    else:
        absorb_discrepancy_into_gauge(
            gauge,
            edges,
            proposals,
            group_name=group_name,
            group_order=group_order,
        )
    return int(edges.size)


def transform_local_frames(
    port_left: np.ndarray,
    port_right: np.ndarray,
    gauge: np.ndarray,
    edge_left: np.ndarray,
    edge_right: np.ndarray,
    frames: np.ndarray,
    *,
    group_name: str,
    group_order: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Apply independent node-frame changes to ports and directed links."""

    left_values, right_values, gauge_values, graph_left, graph_right = _state_arrays(
        port_left,
        port_right,
        gauge,
        edge_left,
        edge_right,
    )
    frame_values = np.asarray(frames, dtype=np.int16)
    node_count = _node_count(graph_left, graph_right)
    if frame_values.ndim != 1 or frame_values.size < node_count:
        raise ValueError("frame array must contain one finite-group element per represented node")
    left_frames = frame_values[graph_left]
    right_frames = frame_values[graph_right]
    transformed_left = group_multiply_indices(
        left_frames,
        left_values,
        group_name=group_name,
        group_order=group_order,
    )
    transformed_right = group_multiply_indices(
        right_frames,
        right_values,
        group_name=group_name,
        group_order=group_order,
    )
    transformed_gauge = group_multiply_indices(
        group_multiply_indices(
            left_frames,
            gauge_values,
            group_name=group_name,
            group_order=group_order,
        ),
        group_inverse_indices(right_frames, group_name=group_name, group_order=group_order),
        group_name=group_name,
        group_order=group_order,
    )
    return transformed_left, transformed_right, transformed_gauge


def canonicalize_gauge_quotient_state(
    port_left: np.ndarray,
    port_right: np.ndarray,
    gauge: np.ndarray,
    *,
    edge_left: np.ndarray | None,
    edge_right: np.ndarray | None,
    group_name: str,
    group_order: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return a deterministic representative of the local-frame gauge orbit.

    Each represented node fixes its frame by sending the local label on its
    lowest-index incident edge to the identity.  The regular left action on a
    group-valued port is free, so this fixes every local frame uniquely.  If
    graph endpoints are absent, every edge endpoint is treated as a distinct
    node, which is the exact quotient of the edge-slot carrier.
    """

    left_values, right_values = _matching_arrays(port_left, port_right, "left ports", "right ports")
    _, gauge_values = _matching_arrays(left_values, gauge, "ports", "gauge")
    graph_left, graph_right = _graph_arrays(edge_left, edge_right, left_values.size)
    node_count = _node_count(graph_left, graph_right)
    frames = np.zeros(node_count, dtype=np.int16)
    if left_values.size:
        sentinel = np.int64(left_values.size)
        reference_edges = np.full(node_count, sentinel, dtype=np.int64)
        edge_indices = np.arange(left_values.size, dtype=np.int64)
        np.minimum.at(reference_edges, graph_left, edge_indices)
        np.minimum.at(reference_edges, graph_right, edge_indices)
        left_is_reference = reference_edges[graph_left] == edge_indices
        right_is_reference = reference_edges[graph_right] == edge_indices
        if np.any(left_is_reference):
            frames[graph_left[left_is_reference]] = group_inverse_indices(
                left_values[left_is_reference],
                group_name=group_name,
                group_order=group_order,
            )
        if np.any(right_is_reference):
            frames[graph_right[right_is_reference]] = group_inverse_indices(
                right_values[right_is_reference],
                group_name=group_name,
                group_order=group_order,
            )
    canonical_left, canonical_right, canonical_gauge = transform_local_frames(
        left_values,
        right_values,
        gauge_values,
        graph_left,
        graph_right,
        frames,
        group_name=group_name,
        group_order=group_order,
    )
    return canonical_left, canonical_right, canonical_gauge, graph_left, graph_right


def gauge_quotient_state_hash(
    port_left: np.ndarray,
    port_right: np.ndarray,
    gauge: np.ndarray,
    *,
    edge_left: np.ndarray | None,
    edge_right: np.ndarray | None,
    group_name: str,
    group_order: int,
) -> str:
    """Hash the exact integer gauge-quotient representative of a coupled state."""

    canonical = canonicalize_gauge_quotient_state(
        port_left,
        port_right,
        gauge,
        edge_left=edge_left,
        edge_right=edge_right,
        group_name=group_name,
        group_order=group_order,
    )
    return _hash_state_arrays(
        *canonical,
        group_name=group_name,
        group_order=group_order,
        schema=f"gauge-quotient:{GAUGE_QUOTIENT_CANONICALIZER}",
    )


def coupled_state_hash(
    port_left: np.ndarray,
    port_right: np.ndarray,
    gauge: np.ndarray,
    *,
    edge_left: np.ndarray | None,
    edge_right: np.ndarray | None,
    group_name: str,
    group_order: int,
) -> str:
    """Hash one exact coupled representative without quotienting its frames."""

    left_values, right_values = _matching_arrays(port_left, port_right, "left ports", "right ports")
    _, gauge_values = _matching_arrays(left_values, gauge, "ports", "gauge")
    graph_left, graph_right = _graph_arrays(edge_left, edge_right, left_values.size)
    return _hash_state_arrays(
        left_values,
        right_values,
        gauge_values,
        graph_left,
        graph_right,
        group_name=group_name,
        group_order=group_order,
        schema="coupled-representative-v1",
    )


def gauge_invariant_edge_residual(
    port_left: np.ndarray,
    port_right: np.ndarray,
    gauge: np.ndarray,
    *,
    group_name: str,
    group_order: int,
) -> np.ndarray:
    """Return ``p_i^-1 g_ij p_j``, an exact local-frame invariant per edge."""

    left_values, right_values = _matching_arrays(port_left, port_right, "left ports", "right ports")
    _, gauge_values = _matching_arrays(left_values, gauge, "ports", "gauge")
    return group_multiply_indices(
        group_multiply_indices(
            group_inverse_indices(left_values, group_name=group_name, group_order=group_order),
            gauge_values,
            group_name=group_name,
            group_order=group_order,
        ),
        right_values,
        group_name=group_name,
        group_order=group_order,
    )


def _hash_state_arrays(
    port_left: np.ndarray,
    port_right: np.ndarray,
    gauge: np.ndarray,
    edge_left: np.ndarray,
    edge_right: np.ndarray,
    *,
    group_name: str,
    group_order: int,
    schema: str,
) -> str:
    hasher = hashlib.sha256()
    header = f"oph-covariant-overlap-state-v1\0{schema}\0{str(group_name).upper()}\0{int(group_order)}\0"
    hasher.update(header.encode("ascii"))
    for values, dtype in (
        (port_left, "<i2"),
        (port_right, "<i2"),
        (gauge, "<i2"),
        (edge_left, "<i8"),
        (edge_right, "<i8"),
    ):
        array = np.ascontiguousarray(np.asarray(values, dtype=dtype))
        hasher.update(int(array.size).to_bytes(8, "little", signed=False))
        hasher.update(array.tobytes())
    return "sha256:" + hasher.hexdigest()


def _normalized_group(group_name: str, group_order: int) -> tuple[str, int]:
    name = str(group_name).upper()
    order = int(group_order)
    if order < 1:
        raise ValueError("finite group order must be positive")
    if name == "S3" and order != int(S3_MUL.shape[0]):
        raise ValueError("S3 overlap operations require the six-element S3 table")
    return name, order


def _matching_arrays(
    left: np.ndarray,
    right: np.ndarray,
    left_name: str,
    right_name: str,
) -> tuple[np.ndarray, np.ndarray]:
    left_values = np.asarray(left)
    right_values = np.asarray(right)
    if left_values.ndim != 1 or right_values.ndim != 1 or left_values.shape != right_values.shape:
        raise ValueError(f"{left_name} and {right_name} must be matching one-dimensional arrays")
    return left_values, right_values


def _state_arrays(
    port_left: np.ndarray,
    port_right: np.ndarray,
    gauge: np.ndarray,
    edge_left: np.ndarray,
    edge_right: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    left_values, right_values = _matching_arrays(port_left, port_right, "left ports", "right ports")
    _, gauge_values = _matching_arrays(left_values, gauge, "ports", "gauge")
    _, graph_left = _matching_arrays(left_values, edge_left, "ports", "left graph endpoints")
    _, graph_right = _matching_arrays(left_values, edge_right, "ports", "right graph endpoints")
    graph_left = np.asarray(graph_left, dtype=np.int64)
    graph_right = np.asarray(graph_right, dtype=np.int64)
    if np.any(graph_left < 0) or np.any(graph_right < 0):
        raise ValueError("graph endpoints must be nonnegative")
    return left_values, right_values, gauge_values, graph_left, graph_right


def _graph_arrays(
    edge_left: np.ndarray | None,
    edge_right: np.ndarray | None,
    edge_count: int,
) -> tuple[np.ndarray, np.ndarray]:
    if edge_left is None and edge_right is None:
        graph_left = np.arange(edge_count, dtype=np.int64)
        graph_right = graph_left + int(edge_count)
        return graph_left, graph_right
    if edge_left is None or edge_right is None:
        raise ValueError("both graph endpoint arrays are required together")
    graph_left = np.asarray(edge_left, dtype=np.int64)
    graph_right = np.asarray(edge_right, dtype=np.int64)
    if graph_left.ndim != 1 or graph_right.ndim != 1 or graph_left.shape != (edge_count,) or graph_right.shape != (edge_count,):
        raise ValueError("graph endpoint arrays must match the edge-slot arrays")
    if np.any(graph_left < 0) or np.any(graph_right < 0):
        raise ValueError("graph endpoints must be nonnegative")
    return graph_left, graph_right


def _node_count(edge_left: np.ndarray, edge_right: np.ndarray) -> int:
    if edge_left.size == 0:
        return 0
    return int(max(int(np.max(edge_left)), int(np.max(edge_right))) + 1)


def overlap_contract_metadata() -> dict[str, Any]:
    return {
        "schema": GAUGE_COVARIANT_OVERLAP_SCHEMA,
        "consistency_equation": "port_left = g_ij * port_right",
        "link_orientation": "g_ij maps the right endpoint frame j into the left endpoint frame i",
        "local_frame_action": {
            "port_i": "h_i * port_i",
            "link_ij": "h_i * g_ij * inverse(h_j)",
        },
        "quotient_canonicalizer": GAUGE_QUOTIENT_CANONICALIZER,
    }
