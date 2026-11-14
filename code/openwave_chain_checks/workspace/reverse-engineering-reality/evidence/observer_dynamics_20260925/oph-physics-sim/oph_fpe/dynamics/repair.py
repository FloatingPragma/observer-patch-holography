from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from oph_fpe.core.patchnet import PatchNet
from oph_fpe.groups.base import GroupElement


@dataclass(frozen=True)
class RepairEvent:
    cycle: int
    node: int
    beta: float
    phi_before: float
    phi_after: float
    delta_phi: float
    accepted: bool
    reason: str


class RepairKernel:
    """Local repair kernel for E0/E1 CPU calibration."""

    def __init__(self, mode: str, hot_metropolis: bool, seed: int):
        self.mode = mode
        self.hot_metropolis = hot_metropolis
        self.rng = np.random.default_rng(seed)

    def step(self, net: PatchNet, cycle: int, beta: float) -> RepairEvent:
        active = _active_nodes(net)
        node = int(self.rng.choice(active if active else list(net.graph.nodes)))
        before = net.touched_phi(node)
        original = net.states[node].copy()

        if self.mode in {"local_best", "local_best_plus_metropolis_hot_phase"}:
            self._propose_local_best(net, node)
        elif self.mode == "random_port":
            self._propose_random_port(net, node)
        else:
            raise ValueError(f"unknown repair mode: {self.mode}")

        after = net.touched_phi(node)
        delta = after - before
        accepted, reason = self._accept(delta, beta)
        if not accepted:
            net.states[node] = original
            after = before
            delta = 0.0
        else:
            improvement = max(0.0, before - after)
            net.states[node].repair_load = 0.9 * net.states[node].repair_load + 0.1 * improvement

        return RepairEvent(
            cycle=cycle,
            node=node,
            beta=beta,
            phi_before=before,
            phi_after=after,
            delta_phi=delta,
            accepted=accepted,
            reason=reason,
        )

    def _propose_local_best(self, net: PatchNet, node: int) -> None:
        state = net.states[node]
        neighbors = list(net.graph.neighbors(node))
        self.rng.shuffle(neighbors)
        for neighbor in neighbors:
            state.ports[neighbor] = net.group.multiply(
                state.gauges[neighbor],
                net.states[neighbor].ports[node],
            )
        if neighbors:
            state.scalar = float(np.mean([net.states[neighbor].scalar for neighbor in neighbors]))
            neighbor_hidden = [net.states[neighbor].hidden for neighbor in neighbors]
            state.hidden = int(max(set(neighbor_hidden), key=neighbor_hidden.count))

    def _propose_random_port(self, net: PatchNet, node: int) -> None:
        neighbors = list(net.graph.neighbors(node))
        if not neighbors:
            return
        neighbor = int(self.rng.choice(neighbors))
        state = net.states[node]
        state.ports[neighbor] = _choice(self.rng, net.group.elements)
        state.scalar = float(state.scalar + self.rng.normal(scale=0.1))

    def _accept(self, delta_phi: float, beta: float) -> tuple[bool, str]:
        if delta_phi <= 0:
            return True, "non_increase"
        if self.hot_metropolis:
            probability = min(1.0, float(np.exp(-beta * delta_phi)))
            if self.rng.random() < probability:
                return True, "metropolis_hot"
        return False, "rejected_increase"


def _active_nodes(net: PatchNet) -> list[int]:
    nodes: set[int] = set()
    for mismatch in net.mismatch_edges():
        nodes.add(mismatch.left)
        nodes.add(mismatch.right)
    return sorted(nodes)


def _choice(rng: np.random.Generator, values: tuple[GroupElement, ...]) -> GroupElement:
    return values[int(rng.integers(0, len(values)))]
