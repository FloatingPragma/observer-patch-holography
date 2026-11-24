"""Spectral multiplicity and a native rotation on the same two wall sources."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.candidate_laws.interval_echo import candidate, currents
from model.executable_model import carrier
from . import belt


def groups(values, tol=1e-8):
    return np.split(np.arange(len(values)), np.where(np.diff(values) > tol)[0] + 1)


def run():
    exchange, cost = law.local_term_parts()
    swap = expm(-1j * np.pi * (exchange + candidate.KAPPA * cost))[:, :9]
    perm, _ = law._load_ops()
    swap_error = np.linalg.norm(swap - np.vstack([perm, np.zeros_like(swap[9:])]), 2)
    poles, order = next(iter(belt.belts().items()))
    # Each pentagon is rotated by a five-cycle of adjacent addressed transpositions.
    sequence = [
        (order[j], order[j + 2]) for i in (0, 1) for j in range(i + 6, i - 1, -2)
    ]
    cell = carrier.cell()
    assert all(tuple(sorted(e)) in cell.seams for e in sequence)
    positions = np.arange(12)
    for a, b in sequence:
        positions[[a, b]] = positions[[b, a]]
    permuted = np.argsort(positions)
    assert all(permuted[p] == order[(j + 2) % 10] for j, p in enumerate(order))
    assert all(permuted[p] == p for p in poles)
    assert {tuple(sorted((permuted[a], permuted[b]))) for a, b in cell.seams} == set(
        cell.seams
    )
    np.testing.assert_array_equal(
        permuted[permuted[permuted[permuted[permuted]]]], np.arange(12)
    )
    experiments = {}
    for n in (3, 4):
        sec = BinarySector(belt.ring(10), n)
        h = sec.hamiltonian(candidate.KAPPA)
        energy, vec = np.linalg.eigh(h)
        top = energy[-10:]
        (source,) = [
            sec.states[i]
            for i in belt.selected(sec)
            if belt.wall(sec.states[i], 10)[0] == 0
        ]
        whole = BinarySector(cell, n)
        full, fv = np.linalg.eigh(whole.hamiltonian(candidate.KAPPA))
        addr = whole.index[tuple(sorted(order[i] for i in source))]
        weight = np.abs(fv[addr]) ** 2
        blocks = groups(full)
        start = {order[i] for i in source}
        occupied = set(start)
        amplitude = 1 + 0j
        for a, b in sequence * 5:
            pair = (int(a in occupied), int(b in occupied))
            amplitude *= swap[law.SRC.index(pair[::-1]), law.SRC.index(pair)]
            if pair[0] != pair[1]:
                occupied ^= {a, b}
        assert occupied == start
        experiments[str(n)] = dict(
            wall_band=[
                dict(energy=float(top[ix].mean()), multiplicity=len(ix))
                for ix in groups(top)
            ],
            smallest_three_level_span=min(top[2:] - top[:-2]),
            wall_eigenpair_residual=np.linalg.norm(
                h @ vec[:, -10:] - vec[:, -10:] * top, 2
            ),
            intact_cell_multiplicity_weights={
                str(m): sum(weight[g].sum() for g in blocks if len(g) == m)
                for m in sorted(set(map(len, blocks)))
            },
            full_turn_source_bright_x=(1 + amplitude.real) / 2,
            spinor_control_bright_x=0,
        )
    return currents._portable(
        dict(
            native_pi_swap_operator_error=swap_error,
            axis=poles,
            fifth_turn_swaps=sequence,
            fifth_turn_permutation=permuted.tolist(),
            cases=experiments,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    encoded = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    with args.output.open("x") as f:
        f.write(encoded)
