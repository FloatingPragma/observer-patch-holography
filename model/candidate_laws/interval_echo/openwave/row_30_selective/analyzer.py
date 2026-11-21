"""Physically compiled wavelength analyzers for the dressed path graph."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import laplacian

from model.candidate_laws.interval_echo import candidate
from .source import Source, STEP, LAYERS


def eigengroups(values):
    return np.split(np.arange(len(values)), np.flatnonzero(np.diff(values) > 1e-8) + 1)


def spatial_basis(source):
    adjacency = -source.h / (candidate.KAPPA / 2)
    np.fill_diagonal(adjacency, 0)
    if not np.array_equal(adjacency, (adjacency > 0).astype(float)):
        raise ArithmeticError("Expected one native hop on every path edge")
    eigenvalues, vectors = np.linalg.eigh(laplacian(adjacency))
    columns = []
    for group in eigengroups(eigenvalues):
        projector = vectors[:, group] @ vectors[:, group].T
        orthogonal = []
        for column in projector.T:
            vector = column.copy()
            for basis in orthogonal:
                vector -= basis * (basis @ vector)
            norm = np.linalg.norm(vector)
            if norm > 1e-7:
                orthogonal.append(vector / norm)
            if len(orthogonal) == len(group):
                break
        columns.extend(orthogonal)
    basis = np.column_stack(columns)
    if not np.allclose(basis.T @ basis, np.eye(len(basis)), atol=1e-9):
        raise ArithmeticError("Spatial analyzer basis is incomplete")
    return eigenvalues, basis


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit("Use a fresh path for independently compiled analyzers")
    source = Source()
    _, basis = spatial_basis(source)
    pulses = []
    for i, target in enumerate(basis.T):
        pulses.append(source.compile(target, 141412 + i))
        print(i, pulses[-1]["infidelity"], flush=True)
    with args.output.open("x") as handle:
        json.dump(dict(ports=source.ports, layers=LAYERS, step=STEP, pulses=pulses), handle,
                  separators=(",", ":"), allow_nan=False)
        handle.write("\n")
