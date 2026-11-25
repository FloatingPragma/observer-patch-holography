"""Independent four-tick material mounts on the same marked Interval Echo source."""

from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from .. import full_source, measure as cuts

HERE = Path(__file__).parent
MOUNTS = ((0, 19, 0), (3, 4, 1), (3, 16, 2))
REGIONAL = cuts.regional.NAMES


def apparatus():
    a = cuts.native.Atlas(1)
    board = a.board
    route = cuts.regional.bridge.native.iso.pairing.formation.region(board)
    start = cuts.regional.bridge.beam.indices(
        sum(cuts.regional.bridge.beam.screens(board, route), [])
    )
    patches = [a.patch(*m) for m in MOUNTS]
    return board, route, start, patches


def channels(board, route, patch, edges, charge, odd):
    h, ops = cuts.native.op_table(board, route, patch, charge, odd)
    rows = {name: ops[name] for name in REGIONAL}
    for u, v in sorted(edges):
        if (v, u) not in edges or u < v:
            rows[f"{u}/{v}"] = edges[u, v], 4
    return h, rows


def source_components(blocks):
    j = cuts.regional.cert.iso.basis(True)
    m = np.array(full_source.local_basis(), dtype=complex)
    return np.r_[
        np.einsum("aij,bji->ba", j, blocks).real.ravel(),
        np.einsum("aij,bji->ba", m, blocks).real.ravel(),
        (np.trace(blocks[0]) - np.trace(blocks[1])).real / np.sqrt(8),
    ]


def float_effects(variations=None, setup=apparatus, regional_only=False):
    board, route, launch, patches = setup()
    edge_sets = [cuts.native.flux_on_edges(board, p) for p in patches]
    by = [{} for _ in patches]
    for q in (-1, 1):
        for odd in (False, True):
            h = cuts.regional.bridge.native.generator(board, route, q, odd)
            g = -2j * h
            starts = np.zeros((h.shape[0], len(launch)))
            starts[launch] = np.eye(len(launch))
            w = expm_multiply(g, starts, traceA=g.diagonal().sum())
            for entries, patch, edges in zip(by, patches, edge_sets):
                mat, ops = channels(board, route, patch, edges, q, odd)
                assert (mat - h).count_nonzero() == 0
                for name, (op, unit) in ops.items():
                    if regional_only and name not in REGIONAL:
                        continue
                    ids = np.unique(op.nonzero()[0])
                    z = w[ids]
                    acted = op[ids, :][:, ids] @ z
                    effect = sum(
                        x.conj().T @ y
                        for x, y in zip(np.hsplit(z, 10), np.hsplit(acted, 10))
                    ) / (10 * unit)
                    entries.setdefault((q, name), []).append(effect)
    if variations is None:
        return [{k: source_components(v) for k, v in entries.items()} for entries in by]
    return [
        {
            k: np.einsum("bij,xbji->x", pair, variations).real
            for k, pair in entries.items()
        }
        for entries in by
    ]
