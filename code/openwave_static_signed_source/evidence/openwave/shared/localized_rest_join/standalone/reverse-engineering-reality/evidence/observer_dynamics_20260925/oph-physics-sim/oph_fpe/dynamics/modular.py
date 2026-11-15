from __future__ import annotations

import math
from typing import Any

import numpy as np

from oph_fpe.core.patchnet import PatchNet
from oph_fpe.core.graph import fibonacci_sphere_points


def apply_modular_flow(net: PatchNet, config: dict[str, Any], cycle: int) -> None:
    """Finite-cutoff modular-flow surrogate.

    This is not the continuum BW theorem. It is a regulator-side diagnostic that
    carries the paper's intended ingredients into the simulator: spherical cap
    chart, quasi-local branch motion, and the 2*pi null-dilation normalization.
    """

    if not config.get("enabled", False):
        return
    dt = float(config.get("dt", 0.05))
    damping = float(config.get("damping", 0.03))
    load_coupling = float(config.get("repair_load_coupling", 0.35))
    cap_coupling = float(config.get("cap_coupling", 0.08))
    diffusion = float(config.get("diffusion", 0.05))
    axes = _cap_axes(int(config.get("cap_axes", 6)))
    two_pi = 2.0 * math.pi

    loads = {node: net.touched_phi(node) / max(1, net.graph.degree[node]) for node in net.graph.nodes}
    mean_load = float(np.mean(list(loads.values()))) if loads else 0.0
    previous_depths = {node: net.states[node].modular_depth for node in net.graph.nodes}

    for node in net.graph.nodes:
        state = net.states[node]
        xyz = np.asarray(net.graph.nodes[node].get("screen_xyz", (0.0, 0.0, 1.0)), dtype=float)
        cap_drive = _cap_drive(xyz, axes, state.modular_time)
        neighbor_depth = _neighbor_mean(previous_depths, net, node)
        centered_load = loads[node] - mean_load
        # The e^{-2*pi t} continuum normalization appears here as the finite
        # regulator's phase advance; the depth variable is the record-visible
        # radial/modular response, not a hidden initialized bulk coordinate.
        state.modular_time += dt * (1.0 + cap_coupling * cap_drive)
        state.modular_depth += dt * (
            load_coupling * centered_load
            + cap_coupling * cap_drive
            + diffusion * (neighbor_depth - state.modular_depth)
            - damping * state.modular_depth
        )
        state.phase = math.exp(-two_pi * state.modular_time)


def collect_modular_sample(net: PatchNet) -> dict[int, float]:
    return {int(node): float(state.modular_depth) for node, state in net.states.items()}


def _cap_axes(count: int) -> np.ndarray:
    return fibonacci_sphere_points(max(1, count))


def _cap_drive(xyz: np.ndarray, axes: np.ndarray, modular_time: float) -> float:
    dots = axes @ xyz
    weights = np.cos((np.arange(len(axes)) + 1.0) * modular_time)
    return float(np.mean(np.tanh(3.0 * dots) * weights))


def _neighbor_mean(depths: dict[int, float], net: PatchNet, node: int) -> float:
    neighbors = list(net.graph.neighbors(node))
    if not neighbors:
        return depths[node]
    return float(np.mean([depths[neighbor] for neighbor in neighbors]))
