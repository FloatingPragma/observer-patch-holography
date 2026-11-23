"""Exact shortest crossings between all internal-seam pair addresses of the 20-cell board."""

from collections import Counter

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components

from ..submission.row_04_check import screen
from .pair_phase_check import BOARD


def crossing():
    edges, lookup, kinetic, adjacent = screen(BOARD)
    # Remove the diagonal and recover the unsigned, integer, single-swap graph.
    off = kinetic - sparse.diags(kinetic.diagonal())
    if not np.all(off.data == -0.5):
        raise ValueError("Unexpected elementary move")
    moves = (-2 * off).astype(np.int64).tocsr()
    seam_ids = np.flatnonzero(adjacent)
    _, blocks = connected_components(moves[seam_ids][:, seam_ids], directed=False)
    home = edges[edges[:, 0] // 12 == edges[:, 1] // 12]
    ids = lookup[tuple(home.T)]
    owners = home[:, 0] // 12
    bridges = edges[edges[:, 0] // 12 != edges[:, 1] // 12]
    block_id = dict(zip(seam_ids, blocks, strict=True))
    home_blocks = [{block_id[i] for i in ids[owners == cell]} for cell in range(20)]
    bridge_blocks = [block_id[lookup[p, q]] for p, q in bridges]
    if (
        any(len(block) != 1 for block in home_blocks)
        or len(set.union(*home_blocks)) != 20
        or len(set(bridge_blocks)) != 30
        or set.union(*home_blocks).intersection(bridge_blocks)
    ):
        raise ArithmeticError("An adjacent-only walk crosses the boundary")

    # A shortest crossing needs both loads to use the bridge. In its middle,
    # either the trailing load enters the bridge or the leading load leaves it first.
    boundary = np.zeros((len(ids), len(ids)), dtype=np.int64)
    for p, q in bridges:
        near = np.flatnonzero(np.any(home == p, axis=1))
        far = np.flatnonzero(np.any(home == q, axis=1))
        boundary[np.ix_(near, far)] = boundary[np.ix_(far, near)] = 2
    different = owners[:, None] != owners
    paths = sparse.csc_matrix(
        (np.ones(len(ids), dtype=np.int64), (ids, np.arange(len(ids)))),
        shape=(len(adjacent), len(ids)),
    )
    for step in range(1, 5):
        paths = moves @ paths
        across = paths[ids].toarray()[different]
        np.testing.assert_array_equal(across, boundary[different] if step == 4 else 0)
    return dict(
        adjacent_components=dict(Counter(np.bincount(blocks).tolist())),
        ordered_four_move_crossings=int(np.count_nonzero(boundary)),
        paths_per_crossing=int(boundary.max()),
    )


if __name__ == "__main__":
    print(crossing())
