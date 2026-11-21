"""Diagnose gateway dependence of the audited low bound-pair frequencies."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import eigsh

from model.executable_model import carrier, federation
from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate
from .experiment import HERE, SEED, lifts_and_geometry, portable


def run():
    sector = PairSector(federation(0))
    cells = sector.car.cells
    lifts, _, _ = lifts_and_geometry(sector)
    start = np.random.default_rng(SEED + 2).normal(size=sector.P)
    h = sector.hamiltonian(candidate.KAPPA)
    energies, states = eigsh(h, k=cells+4, ncv=84, which="SA", v0=start, tol=1e-12)
    energies, states = energies[:cells], states[:, :cells]
    residual = np.max(np.linalg.norm(h @ states - energies*states, axis=0))
    if residual > 1e-9:
        raise ArithmeticError("inaccurate eigenspace")
    projections = states.T @ lifts["uniform_seams"]
    eigenvalues, directions = np.linalg.eigh(projections.T @ projections)
    if min(eigenvalues) < .5:
        raise ArithmeticError("cell orbitals fail to span the low band reliably")
    normalized = projections @ ((directions / np.sqrt(eigenvalues)) @ directions.T)
    np.testing.assert_allclose(normalized.T @ normalized, np.eye(cells), atol=1e-10)
    effective = normalized.T @ ((energies - energies[0])[:, None] * normalized)
    split = int(np.argmax(np.diff(energies))) + 1
    lower_weight = np.sum(normalized[:split]**2, axis=0)
    d = shortest_path(carrier.cell().adjacency(), directed=False, unweighted=True)
    by_geometry = {}
    for cell in range(cells):
        ports = sorted(p%12 for e in sector.car.external for p in e if p//12 == cell)
        distances = sorted(int(d[a,b]) for j, a in enumerate(ports) for b in ports[j+1:])
        by_geometry.setdefault(",".join(map(str, distances)), []).append(cell)
    groups = {}
    for key, indices in by_geometry.items():
        diag = np.diag(effective)[indices]
        weights = lower_weight[indices]
        groups[key] = dict(cells=indices, diagonal_range=[min(diag), max(diag)],
                           lower_cluster_weight_range=[min(weights), max(weights)])
    return portable(dict(candidate=candidate.NAME, split=split, cluster_gap=energies[split]-energies[split-1],
                         min_overlap_eigenvalue=min(eigenvalues), groups=groups,
                         maximum_off_diagonal=np.max(abs(effective-np.diag(np.diag(effective)))),
                         effective_relative_to_bottom=effective.tolist(),
                         lower_cluster_weight_by_cell=lower_weight.tolist()))


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, help="write a separate receipt to a new path")
    output = cli.parse_args().output
    record = run()
    if output:
        with output.open("x") as dst:
            json.dump(record, dst, indent=2, allow_nan=False)
            dst.write("\n")
    elif record != json.loads((HERE/"gateways.json").read_text()):
        raise SystemExit("gateway assignment differs from the archive")
    for group, values in record["groups"].items():
        print(group, values["lower_cluster_weight_range"])
