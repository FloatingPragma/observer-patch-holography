"""Coherent returns of a three-load ground state in a rotating scalar face well."""

import argparse
import json
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import eigsh

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.candidate_laws.interval_echo import candidate, currents
from model.executable_model import carrier


def panel(sec, diagonals, k, steps):
    base = sec.hamiltonian(k)
    # All five internal contacts are active. The original degree reference
    # retains the relative one-unit difference at outgoing gateways.
    ref = 6 - carrier.federation(0).degree[:12]
    diag = sum(ref[np.array(sec.states).T])
    base[np.diag_indices_from(base)] += k * diag
    h = csr_matrix(base)
    states = []
    gaps = []
    for i, left in enumerate(diagonals):
        right = diagonals[(i + 1) % 5]
        for j in range(steps):
            f = j / steps
            operator = h + diags((1 - f) * left + f * right)
            val, vec = eigsh(
                operator, k=2, which="SA", v0=np.linspace(1, 2, len(base)), tol=1e-12
            )
            q = vec[:, 0]
            if q.sum() < 0:
                q = -q
            assert q.min() > -1e-10 and np.linalg.norm(operator @ q - val[0] * q) < 1e-9
            states.append(q)
            gaps.append(val[1] - val[0])
    links = np.array(
        [np.vdot(q, states[(i + 1) % len(states)]) for i, q in enumerate(states)]
    )
    amplitude = float(np.prod(links))

    # An ideal nondestructive energy filter in the moving arm heralds a miss
    # into an orthogonal terminal port. The other arm holds the initial state.
    def receiver(a):
        return [(1 + a) ** 2 / 4, (1 - a) ** 2 / 4, (1 - a * a) / 2]

    return dict(
        steps_per_face=steps,
        minimum_gap=min(gaps),
        minimum_overlap=min(links),
        cyclic_probability=amplitude**2,
        geometric_phase=float(np.angle(amplitude)),
        reversed_phase=float(np.angle(np.prod(links[::-1]))),
        bright_dark_miss=receiver(amplitude),
        supplied_negative_return=receiver(-amplitude),
    )


def run():
    cell = carrier.cell()
    rim = nx.Graph(cell.seams).subgraph(cell.neighbours[0])
    (cycle,) = nx.cycle_basis(rim)
    cycle = min(
        tuple(seq[i:] + seq[:i]) for seq in (cycle, cycle[::-1]) for i in range(5)
    )
    faces = [(0, p, cycle[(i + 1) % 5]) for i, p in enumerate(cycle)]
    sector = BinarySector(cell, 3)
    strength = 4
    potentials = [
        -strength * np.array([len(set(w) & set(f)) for w in sector.states])
        for f in faces
    ]
    return currents._portable(
        dict(
            faces=faces,
            well_strength=strength,
            panels={
                str(k): [panel(sector, potentials, k, m) for m in (4, 12)]
                for k in (int(candidate.KAPPA), 0)
            },
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    result = run()
    with args.output.open("x") as f:
        f.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
