"""Source preparation and home-port observations of a blank-federation optical mode."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import linalg, sparse
from scipy.sparse.csgraph import connected_components, shortest_path, laplacian
from scipy.sparse.linalg import expm_multiply

from model.executable_model import cell, federation
from ..row_30_formed.experiment import synthesize

HERE = Path(__file__).resolve().parent
SHOTS = 2**23


def preparations(board):
    raw = shortest_path(cell().adjacency(), directed=False, unweighted=True)
    # Opposite ports are unique at internal distance three.
    assert np.all((raw == 3).sum(axis=1) == 1)
    left, right = np.where(np.triu(raw == 3))
    axis = np.empty(12, dtype=int)
    axis[left] = axis[right] = np.arange(6)
    axes = (axis + 6 * np.arange(board.cells)[:, None]).ravel()
    edge = np.array(board.external)
    endpoints = edge // 12
    cells = sparse.coo_matrix(
        (np.ones(len(endpoints)), endpoints.T), shape=(board.cells, board.cells)
    ).tocsr()
    cells = cells + cells.T
    assert np.all(cells.getnnz(axis=1) == 3)
    lcell = laplacian(cells).toarray()
    gateway = np.zeros(6 * board.cells, dtype=bool)
    gateway[axes[edge.ravel()]] = True
    gateway = gateway.reshape(-1, 6)
    n = gateway.sum(axis=1)[:, None]
    profile = np.where(gateway, 1.0, -n / (6 - n))[:, axis]
    norm = np.sum(profile**2, axis=1)
    return axes, edge, profile, norm, lcell


def exact_stationary(board, axes, edges):
    nodes = axes[edges]
    graph = sparse.coo_matrix(
        (np.ones(len(nodes)), nodes.T), shape=(6 * board.cells, 6 * board.cells)
    ).tocsr()
    groups, labels = connected_components(graph, directed=False)
    site_labels = labels[axes]
    weights = np.bincount(site_labels)
    incidence = sparse.coo_matrix(
        (np.ones(len(labels)), (np.arange(len(labels)) // 6, labels)),
        shape=(board.cells, groups),
    ).tocsr()
    inv_weight = sparse.diags(1.0 / weights)
    # Each cell has its own non-gateway axes, so these constraints are independent.
    solve = linalg.cho_factor((incidence @ inv_weight @ incidence.T).toarray())

    def project(psi):
        avg = np.bincount(site_labels, weights=psi) / weights
        enforced = linalg.cho_solve(solve, incidence @ avg)
        return (avg - (incidence.T @ enforced) / weights)[site_labels]

    return groups - board.cells, project


def spatial_modes(lcell, norm):
    eig, modes = linalg.eigh(lcell, np.diag(norm))
    levels = []
    for x in eig:
        if not levels or x > levels[-1] + 1e-9:
            levels.append(x)
    for x in levels[:6]:
        basis = modes[:, abs(eig - x) < 1e-9]
        anchor = np.flatnonzero(np.sum(basis**2, axis=1) > 1e-10)[0]
        wave = basis @ basis[anchor]
        wave /= np.sqrt(np.dot(norm, wave**2))
        yield float(x), wave


def counts(rng, probabilities):
    # Round at a declared precision and assign the remainder to the last bin.
    # Reject unstable or nonnormalised input before drawing a detector outcome.
    p = np.round(probabilities, 12).ravel()
    assert p.min() >= 0 and abs(sum(p) - 1) < 1e-10
    return rng.multinomial(SHOTS, p).tolist()


def panel(level):
    board = federation(level)
    axes, edges, profile, norm, lcell = preparations(board)
    dim, project = exact_stationary(board, axes, edges)
    h = board.laplacian() / 2  # uniform reference removes the common 48
    local = sparse.kron(sparse.eye(board.cells), cell().laplacian() / 2, format="csr")
    root = int(np.flatnonzero(board.degree == 5)[0])
    gate = int(np.min(edges))
    uniform = np.full(board.n, 1 / np.sqrt(board.n))
    psi0 = profile.ravel() / np.sqrt(norm.sum())
    assert np.linalg.norm(h @ psi0 - 3 * psi0) < 1e-12
    rng = np.random.default_rng(303440 + level)
    readings = []
    ports = board.adjacency()
    # Tomography uses native spatial preparation up to an arbitrary common phase.
    for ratio in (1, 2):
        generator = h + (ratio - 1) * local
        for phase in (1, 1j, -1, -1j):
            source, _ = synthesize(
                (uniform + phase * psi0) / np.sqrt(2), ports, board.degree, root
            )
            for t in (0.0, 0.125, 0.25, 0.5, 0.75, 1.0):
                out = expm_multiply(-1j * t * generator, source)
                hit = abs(out[gate]) ** 2
                hand = np.sin(t / 2) ** 2
                readings.append(
                    dict(
                        internal_exchange=ratio,
                        phase=[phase.real, phase.imag],
                        appointment=t,
                        counts=counts(rng, np.outer([1 - hand, hand], [1 - hit, hit])),
                    )
                )
    modes = []
    for i, (wave_number_squared, amplitudes) in enumerate(spatial_modes(lcell, norm)):
        target = (profile * amplitudes[:, None]).ravel()
        psi, _ = synthesize(target, ports, board.degree, root)
        frozen = project(target)
        assert np.linalg.norm(h @ frozen - 3 * frozen) < 1e-12
        contrast = psi[edges[:, 0]] - psi[edges[:, 1]]
        p_odd = np.vdot(contrast, contrast).real / (2 * len(edges))
        variance = np.linalg.norm(h @ psi - np.vdot(psi, h @ psi) * psi) ** 2
        delta = wave_number_squared / 2
        assert abs(variance - delta * (1 - delta)) < 1e-12
        visits = []
        for t in (0.25, 1.0, 4.0):
            out = expm_multiply(-1j * t * h, psi)
            miss = np.linalg.norm(out - np.vdot(psi, out) * psi) ** 2
            hand = np.sin(t / 2) ** 2
            visits.append(
                [t, counts(rng, np.outer([1 - hand, hand], [1 - miss, miss]))]
            )
        modes.append(
            dict(
                index=i,
                k2=round(wave_number_squared, 10),
                odd=rng.binomial(SHOTS, p_odd),
                variance=round(variance, 10),
                stationary_weight=round(float(frozen @ frozen), 10),
                visits=visits,
            )
        )
    localized = np.zeros(board.n)
    localized[:12] = profile[0] / np.sqrt(norm[0])
    src, _ = synthesize(localized, ports, board.degree, root)
    distances = shortest_path(
        sparse.csr_matrix(np.diag(np.diag(lcell)) - lcell), unweighted=True, indices=0
    ).astype(int)
    flight = []
    for t in (0, 2, 8, 32):
        out = expm_multiply(-1j * t * h, src)
        p = np.bincount(distances, weights=(abs(out.reshape(-1, 12)) ** 2).sum(axis=1))
        flight.append([t, counts(rng, p)])
    return dict(
        cells=board.cells,
        shared_seams=len(edges),
        stationary_dimension=dim,
        source=root,
        detector=gate,
        optical_at_detector=round(float(psi0[gate]), 12),
        reference_at_detector=round(float(uniform[gate]), 12),
        spatial=modes,
        local_source_stationary=round(
            float(np.linalg.norm(project(localized)) ** 2), 10
        ),
        occupation_by_hops=flight,
        interferometry=readings,
    )


def run():
    return dict(
        shots=SHOTS,
        joint_order=["hand0_quiet", "hand0_yes", "hand1_quiet", "hand1_yes"],
        panels=[panel(level) for level in range(3)],
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--record", action="store_true")
    args = p.parse_args()
    result = run()
    path = HERE / "results.json"
    if args.record:
        with path.open("x") as f:
            json.dump(result, f, indent=1)
            f.write("\n")
    elif result != json.loads(path.read_text()):
        raise SystemExit("Optical cell measurement differs from the saved receipt")
    else:
        print(
            "Verified joint clocks, local phases, spatial returns and flights on three intact boards."
        )


if __name__ == "__main__":
    main()
