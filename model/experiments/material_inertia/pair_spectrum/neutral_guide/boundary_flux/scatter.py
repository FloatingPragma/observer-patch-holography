"""Unrestricted flight on the actual memory register, with independent local-defect controls."""

import numpy as np
from scipy.linalg import eigh

from model.candidate_laws.contrast_memory_spherical_v1 import law
from .field import mount_geometry

TIMES = (1, 4, 12, 32)
SHOTS = 2**20


def carrier_cut(board, seam, sign):
    adj = board.adjacency().toarray()
    diagonal = 8.5 * board.degree + 8 * (6 - board.degree) * ((1 + sign)**2 - 1)
    full = np.diag(diagonal) - .5 * adj
    cut = full.copy()
    a, b = seam
    cut[a, a] -= 8.5
    cut[b, b] -= 8.5
    cut[a, b] = cut[b, a] = 0
    return full, cut


def local_embedding(labels, loads):
    mat = np.zeros((law.DIM, len(labels)))
    for i, d in enumerate(labels):
        for c in (0, 1):
            address = 9 * (1 + law.LABELS.index((c, d))) + 3 * loads[0] + loads[1]
            mat[address, i] = 1 / np.sqrt(2)
    return mat


def central_blocks(local, sign, labels):
    words = [(1, 1), (1 + sign, 1), (1, 1 + sign)]
    embeddings = [local_embedding(labels, w) for w in words]
    blocks = []
    for embed in embeddings:
        block = embed.T @ local @ embed
        # This norm checks all other labels, coin directions and ternary port words.
        np.testing.assert_allclose(local @ embed, embed @ block, atol=3e-13, rtol=0)
        blocks.append(block)
    return blocks


def written_generator(cut, seam, blocks):
    n = len(cut)
    h = np.kron(np.eye(len(blocks[0])), cut)
    for port in range(n):
        rows = port + n * np.arange(len(blocks[0]))
        sector = 1 if port == seam[0] else 2 if port == seam[1] else 0
        h[np.ix_(rows, rows)] += blocks[sector]
    return h


def flights(h, n, launches, source_memory):
    energies, basis = eigh(h)
    input_vectors = np.kron(np.asarray(source_memory)[:, None], np.eye(n)[:, launches])
    amplitude = basis.T @ input_vectors
    return {t: (abs((basis @ (np.exp(-1j * t * energies)[:, None] * amplitude))
                    .reshape(len(source_memory), n, len(launches)))**2).sum(axis=0) for t in TIMES}


def measure(board, seams, local):
    rng = np.random.default_rng(191102)
    rows = []
    diagnostics = []
    for seam in seams:
        _, launches = mount_geometry(board, seam)
        for sign in (1, -1):
            full, cut = carrier_cut(board, seam, sign)
            extreme = central_blocks(local, sign, [-2, 2])
            ordinary = central_blocks(local, sign, [0])
            he = written_generator(cut, seam, extreme)
            ho = written_generator(cut, seam, ordinary)
            p = dict(host=flights(full, board.n, launches, [1]), cut=flights(cut, board.n, launches, [1]),
                     write02=flights(he, board.n, launches, [1, 0]), write20=flights(he, board.n, launches, [0, 1]),
                     write11=flights(ho, board.n, launches, [1]))
            potential = np.zeros(board.n)
            potential[list(seam)] = 1
            barriers = [flights(cut + e * np.diag(potential), board.n, launches, [1]) for e in (32, 33)]
            twin_error = max(float(np.max(abs(p['write02'][t] - (barriers[0][t] + barriers[1][t]) / 2)))
                             for t in TIMES)
            diagnostics.append(dict(seam=list(seam), sign=sign, extreme_blocks=[m.tolist() for m in extreme],
                ordinary_blocks=[m.tolist() for m in ordinary], exact_barrier_twin_max_port_error=twin_error,
                reversal_max_port_error=max(float(np.max(abs(p['write02'][t] - p['write20'][t]))) for t in TIMES),
                cut_max_port_error=max(float(np.max(abs(p['write11'][t] - p['cut'][t]))) for t in TIMES)))
            for t in TIMES:
                for col, site in enumerate(launches):
                    records = {}
                    for source, times in p.items():
                        ports = times[t][:, col]
                        np.testing.assert_allclose(ports.sum(), 1, atol=3e-12, rtol=0)
                        owners = ports.reshape(board.cells, 12).sum(axis=1)
                        hits = float(sum(ports[list(seam)]))
                        records[source] = dict(owner_counts=rng.multinomial(SHOTS, owners / owners.sum()).tolist(),
                            endpoint_hits=int(rng.binomial(SHOTS, min(1, hits))),
                            expected_endpoint_probability=hits)
                    rows.append(dict(seam=list(seam), sign=sign, time=t, launch=site, receiver_records=records))
    return dict(seed=191102, shots=SHOTS, source_order=['host', 'cut', 'write02', 'write20', 'write11'],
                diagnostics=diagnostics, banks=rows)
