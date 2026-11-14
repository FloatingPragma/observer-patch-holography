from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import networkx as nx
import numpy as np

from oph_fpe.core.patch_state import PatchState
from oph_fpe.groups.base import FiniteGroup, GroupElement


@dataclass
class EdgeMismatch:
    left: int
    right: int
    weight: float
    distance: float

    @property
    def phi(self) -> float:
        return self.weight * self.distance


class PatchNet:
    def __init__(self, graph: nx.Graph, group: FiniteGroup, states: dict[int, PatchState]):
        self.graph = graph
        self.group = group
        self.states = states

    @classmethod
    def random(cls, graph: nx.Graph, group: FiniteGroup, seed: int) -> "PatchNet":
        rng = np.random.default_rng(seed)
        states = {
            node: PatchState(
                hidden=int(rng.integers(0, len(group.elements))),
                scalar=float(rng.normal()),
                phase=float(rng.random()),
                modular_time=0.0,
                modular_depth=float(rng.random()),
                repair_load=0.0,
                capacity=graph.degree[node],
            )
            for node in graph.nodes
        }
        for left, right in graph.edges:
            left_port = _choice(rng, group.elements)
            right_port = _choice(rng, group.elements)
            gauge = _choice(rng, group.elements)
            states[left].ports[right] = left_port
            states[right].ports[left] = right_port
            states[left].gauges[right] = gauge
            states[right].gauges[left] = group.inverse(gauge)
        return cls(graph, group, states)

    @classmethod
    def synchronized(cls, graph: nx.Graph, group: FiniteGroup) -> "PatchNet":
        states = {
            node: PatchState(hidden=0, scalar=0.0, phase=0.0, capacity=graph.degree[node])
            for node in graph.nodes
        }
        for left, right in graph.edges:
            states[left].ports[right] = group.identity
            states[right].ports[left] = group.identity
            states[left].gauges[right] = group.identity
            states[right].gauges[left] = group.identity
        return cls(graph, group, states)

    def edge_mismatch(self, left: int, right: int) -> EdgeMismatch:
        weight = float(self.graph.edges[left, right].get("weight", 1.0))
        transported_right = self.group.multiply(
            self.states[left].gauges[right],
            self.states[right].ports[left],
        )
        distance = self.group.mismatch(self.states[left].ports[right], transported_right)
        return EdgeMismatch(left=left, right=right, weight=weight, distance=distance)

    def total_phi(self) -> float:
        return float(sum(self.edge_mismatch(left, right).phi for left, right in self.graph.edges))

    def touched_phi(self, node: int) -> float:
        return float(sum(self.edge_mismatch(node, neighbor).phi for neighbor in self.graph.neighbors(node)))

    def mismatch_edges(self) -> list[EdgeMismatch]:
        return [
            mismatch
            for left, right in self.graph.edges
            if (mismatch := self.edge_mismatch(left, right)).distance > 0
        ]

    def clone_states(self) -> dict[int, PatchState]:
        return {node: state.copy() for node, state in self.states.items()}

    def as_jsonable(self) -> dict[str, Any]:
        return {
            "group": self.group.name,
            "nodes": {
                str(node): state.to_jsonable(self.group.label)
                for node, state in sorted(self.states.items())
            },
            "edges": [
                {
                    "left": left,
                    "right": right,
                    "weight": float(self.graph.edges[left, right].get("weight", 1.0)),
                    "mismatch": self.edge_mismatch(left, right).distance,
                }
                for left, right in sorted(self.graph.edges)
            ],
        }

    def set_directed_gauge(self, left: int, right: int, value: GroupElement) -> None:
        if not self.graph.has_edge(left, right):
            raise ValueError(f"no edge {left}-{right}")
        self.states[left].gauges[right] = value
        self.states[right].gauges[left] = self.group.inverse(value)


def _choice(rng: np.random.Generator, values: tuple[GroupElement, ...]) -> GroupElement:
    return values[int(rng.integers(0, len(values)))]
