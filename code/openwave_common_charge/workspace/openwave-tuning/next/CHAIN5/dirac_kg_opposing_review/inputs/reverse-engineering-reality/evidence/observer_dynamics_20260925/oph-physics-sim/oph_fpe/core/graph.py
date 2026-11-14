from __future__ import annotations

import math
from itertools import product

import networkx as nx
import numpy as np
from scipy.spatial import cKDTree

from oph_fpe.core.icosahedral import (
    geodesic_icosahedral_graph,
    resolve_icosahedral_level,
    supported_icosahedral_count,
)


def build_patch_graph(config: dict, seed: int) -> nx.Graph:
    family = config.get("family", "cycle")
    patch_count = int(config.get("patch_count", config.get("nodes", 16)))
    rng = np.random.default_rng(seed)

    if family == "cycle":
        graph = nx.cycle_graph(patch_count)
    elif family == "path":
        graph = nx.path_graph(patch_count)
    elif family == "grid_2d":
        graph = _grid_2d(patch_count)
    elif family == "lattice_3d":
        graph = _lattice_3d(patch_count)
    elif family == "random_regular":
        degree = int(config.get("degree", 4))
        graph = nx.random_regular_graph(degree, patch_count, seed=seed)
    elif family == "small_world_screen":
        degree = int(config.get("degree", 6))
        rewire = float(config.get("rewire_probability", 0.02))
        graph = nx.watts_strogatz_graph(patch_count, degree, rewire, seed=seed)
    elif family in {
        "subdivided_icosahedral_screen",
        "nested_geodesic_icosahedral",
        "geodesic_icosahedral_refinement",
    }:
        patch_basis = str(config.get("patch_basis", config.get("node_basis", "cells")))
        if patch_basis not in {"cells", "vertices"}:
            raise ValueError("icosahedral patch_basis must be 'cells' or 'vertices'")
        nominal_patch_count = config.get("nominal_patch_count")
        if "refinement_level" in config:
            level = int(config["refinement_level"])
            expected_count = supported_icosahedral_count(level, patch_basis)
            if "patch_count" in config and int(config["patch_count"]) != expected_count:
                raise ValueError(
                    f"refinement_level={level} with patch_basis={patch_basis!r} has "
                    f"exactly {expected_count} patches, not {int(config['patch_count'])}; "
                    "use nominal_patch_count for a campaign budget label"
                )
        elif "patch_count" in config or "nodes" in config:
            policy = str(config.get("patch_count_policy", "exact"))
            level, count_mapping = resolve_icosahedral_level(
                patch_count,
                patch_basis=patch_basis,
                policy=policy,
            )
            if not count_mapping["exact_supported_count"] and nominal_patch_count is None:
                nominal_patch_count = patch_count
        else:
            level = 0
        graph = geodesic_icosahedral_graph(
            level,
            patch_basis=patch_basis,
            nominal_patch_count=(
                None if nominal_patch_count is None else int(nominal_patch_count)
            ),
        )
    elif family in {"fibonacci_sphere", "spherical_knn_screen"}:
        neighbors = int(config.get("neighbors", config.get("degree", 8)))
        graph = _fibonacci_sphere_graph(patch_count, neighbors)
    else:
        raise ValueError(f"unknown graph family: {family}")

    for left, right in graph.edges:
        graph.edges[left, right]["weight"] = float(config.get("edge_weight", 1.0))
    return graph


def fibonacci_sphere_points(count: int) -> np.ndarray:
    if count < 1:
        return np.zeros((0, 3), dtype=float)
    points = np.zeros((count, 3), dtype=float)
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    for index in range(count):
        y = 1.0 - (2.0 * index + 1.0) / count
        radius = math.sqrt(max(0.0, 1.0 - y * y))
        theta = golden_angle * index
        points[index] = (math.cos(theta) * radius, y, math.sin(theta) * radius)
    return points


def _fibonacci_sphere_graph(patch_count: int, neighbors: int) -> nx.Graph:
    points = fibonacci_sphere_points(patch_count)
    graph = nx.Graph()
    for node, point in enumerate(points):
        graph.add_node(node, screen_xyz=tuple(float(value) for value in point))
    if patch_count <= 1:
        return graph
    k = min(patch_count, max(2, neighbors + 1))
    tree = cKDTree(points)
    _, indices = tree.query(points, k=k)
    if indices.ndim == 1:
        indices = indices[:, None]
    neighbor_sets = {node: set(int(value) for value in row[1:]) for node, row in enumerate(indices)}
    for node, local_neighbors in neighbor_sets.items():
        for neighbor in local_neighbors:
            if node in neighbor_sets.get(neighbor, set()):
                graph.add_edge(int(node), int(neighbor))
    return graph


def _grid_2d(patch_count: int) -> nx.Graph:
    side = max(2, round(math.sqrt(patch_count)))
    graph = nx.grid_2d_graph(side, side)
    graph = nx.convert_node_labels_to_integers(graph)
    return graph.subgraph(range(min(patch_count, graph.number_of_nodes()))).copy()


def _lattice_3d(patch_count: int) -> nx.Graph:
    side = max(2, round(patch_count ** (1.0 / 3.0)))
    graph = nx.Graph()
    nodes = list(product(range(side), range(side), range(side)))
    graph.add_nodes_from(nodes)
    for x, y, z in nodes:
        for dx, dy, dz in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            nxt = (x + dx, y + dy, z + dz)
            if nxt in graph:
                graph.add_edge((x, y, z), nxt)
    graph = nx.convert_node_labels_to_integers(graph)
    return graph.subgraph(range(min(patch_count, graph.number_of_nodes()))).copy()
