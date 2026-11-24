"""Finite programs making the period-three reference with actual cell contacts."""

import argparse
import json
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.candidate_laws.interval_echo import candidate, currents
from model.executable_model import carrier
from . import bulk_order


def pairs(sec, a, b):
    left = []
    right = []
    for i, sites in enumerate(sec.states):
        if a in sites and b not in sites:
            left.append(i)
            right.append(sec.index[tuple(sorted((set(sites) - {a}) | {b}))])
    return np.array(left), np.array(right)


def pulse(state, ids, time, k):
    a, b = ids
    if time == np.pi and k in (0, 8):
        state[a], state[b] = state[b], state[a]
        return
    phase = np.exp(-1j * time)
    common = np.exp(-1j * k * time)
    x = common * (1 + phase) / 2
    y = common * (1 - phase) / 2
    state[a], state[b] = x * state[a] + y * state[b], y * state[a] + x * state[b]


def run():
    _, ring, _, graphs = bulk_order.geometry()
    cell = carrier.cell()
    graph = nx.Graph(cell.seams)
    square = carrier.Carrier("comparison", 12, graphs["completed_square"])
    sector = BinarySector(square, 4)
    tools = {e: pairs(sector, *e) for e in cell.seams}
    program = []
    for a, b in sorted(square.seams):
        route = min(nx.all_shortest_paths(graph, a, b))
        keys = [tuple(sorted(e)) for e in zip(route, route[1:])]
        program.append(keys)
    swaps_per_step = 4 * sum(len(keys) - 1 for keys in program)
    occupation = np.array([[int(p in s) for p in ring] for s in sector.states])
    launch = [
        sector.index[tuple(sorted(ring[(j + i) % 12] for j in (0, 3, 6, 9)))]
        for i in (0, 4, 8)
    ]
    tx = np.exp(2j * np.pi * (occupation @ np.arange(12)) / 12)
    start = np.eye(len(sector.states))[:, launch].astype(complex)
    data = []
    for k in (int(candidate.KAPPA), 0):
        # Check the complete native pulse on the nine blank input columns.
        h0, h1 = law.local_term_parts()
        u = expm(-1j * np.pi * (h0 + k * h1))
        for j, (a, b) in enumerate(law.SRC):
            np.testing.assert_allclose(
                u[:, j], np.eye(99)[:, law.SRC.index((b, a))], atol=3e-13
            )
        e, v = np.linalg.eigh(sector.hamiltonian(k))
        total = 4
        exact = v @ (np.exp(-1j * total * e)[:, None] * (v.T @ start))
        for steps in (4, 16, 64):
            state = start.copy()
            for _ in range(steps):
                for keys in (*program, *reversed(program)):
                    for edge in keys[:-1]:
                        pulse(state, tools[edge], np.pi, k)
                    pulse(state, tools[keys[-1]], total / (2 * steps), k)
                    for edge in reversed(keys[:-1]):
                        pulse(state, tools[edge], np.pi, k)
            # The shifted square Hamiltonian differs by a constant in this sector.
            prob = abs(state) ** 2
            turns = prob.T @ tx
            data.append(
                dict(
                    kappa=k,
                    steps=steps,
                    signal_time=total,
                    serial_native_exposure=len(program) * total
                    + np.pi * swaps_per_step * steps,
                    overlaps=abs(np.sum(state.conj() * exact, axis=0)) ** 2,
                    crystalline=prob[launch].sum(axis=0),
                    twist_visibility=abs(turns),
                    phase_differences_turns=np.angle(np.roll(turns, -1) * turns.conj())
                    / (2 * np.pi),
                    norm_errors=abs(prob.sum(axis=0) - 1),
                )
            )
    return currents._portable(
        dict(routes=program, swaps_per_step=swaps_per_step, panels=data)
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    result = run()
    with a.output.open("x") as f:
        f.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
