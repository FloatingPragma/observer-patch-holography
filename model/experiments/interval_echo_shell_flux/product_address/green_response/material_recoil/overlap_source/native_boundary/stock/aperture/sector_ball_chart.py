"""Fixed-volume marker balls read by reversible ordinary round trips."""

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import shortest_path

from . import closed_balls as balls
from . import landmark_arrival as arrival


def charts():
    words = arrival.source.records()[0]
    board, homes, seams, *_ = balls.energy.returned(words)
    _, marks, tables, _ = arrival.pulse_inputs()
    assert marks == seams
    n = board.n
    base, hops, _ = balls.radar.chart(board, seams)
    edges = np.array(board.seams)
    geometry = {}
    for load, table in tables.items():
        # A path and its reversal both have nonzero native ticket weight.
        # The shortest return through a target therefore uses twice the
        # shortest lifted outward distance while recovering the initial head.
        tickets = sparse.coo_matrix(
            (
                np.ones(table.size, dtype=np.int16),
                (np.repeat(np.arange(len(table)), table.shape[1]), table.ravel()),
            ),
            shape=(len(table), len(table)),
        ).tocsr()
        assert (tickets != tickets.T).nnz == 0
        panels = []
        for h in range(11):
            paths = shortest_path(
                tickets,
                directed=True,
                unweighted=True,
                indices=np.arange(h * n, (h + 1) * n),
            )
            dist = paths.reshape(n, 11, n).min(1)
            assert np.isfinite(dist).all()
            hops_h = dist.astype(np.int16)
            adj = hops_h == 1
            # An initial head may open a directed passage. Include both
            # orientations even when the material echo retraces its path.
            assert np.array_equal(
                shortest_path(sparse.csr_matrix(adj), directed=True, unweighted=True),
                hops_h,
            )
            accessible = np.column_stack(
                (adj[edges[:, 0], edges[:, 1]], adj[edges[:, 1], edges[:, 0]])
            )
            panels.append((hops_h, accessible))
        geometry[load] = panels
    return words, hops, homes, seams, edges, geometry


def regional_masks(hops, homes, seams):
    distances = np.column_stack(
        (hops[:, homes], *[hops[:, [a, b]].max(1) for a, b in seams])
    )
    return {r: (distances < r) @ (2 ** np.arange(7)) for r in balls.RADII}


def fixed_volumes(hops, homes, seams, edges, geometry):
    """Keep all control volumes and both tie orders, including insensitive cuts."""
    markers = np.arange(len(hops))
    u, v = edges.T
    masks = 2 ** np.arange(7)
    result = {}
    for load, panels in geometry.items():
        # Identical reachable-range tables may serve several measured heads.
        classes = []
        gids = []
        for d, e in panels:
            j = next(
                (
                    i
                    for i, (old, oe) in enumerate(classes)
                    if np.array_equal(d, old) and np.array_equal(e, oe)
                ),
                len(classes),
            )
            if j == len(classes):
                classes.append((d, e))
            gids.append(j)
        for descending in (False, True):
            stable = markers[::-1] if descending else markers
            ranks = []
            for d, _ in classes:
                order = stable[np.argsort(d[:, stable], axis=1, kind="stable")]
                ranks.append(np.argsort(order, axis=1))
            for r in range(2, int(hops.max()) + 2):
                size = (hops < r).sum(1)
                interiors = [rank < size[:, None] for rank in ranks]
                cuts = []
                for mask, (_, e) in zip(interiors, classes, strict=True):
                    forward = mask[:, u] & ~mask[:, v]
                    backward = mask[:, v] & ~mask[:, u]
                    out = (forward & e[:, 0]).sum(1) + (backward & e[:, 1]).sum(1)
                    inward = (forward & e[:, 1]).sum(1) + (backward & e[:, 0]).sum(1)
                    cuts.append(np.column_stack((out, inward)))
                cut = np.array(cuts)
                records = np.array(
                    [
                        np.column_stack(
                            (
                                mask[:, homes],
                                *[mask[:, a] & mask[:, b] for a, b in seams],
                            )
                        )
                        @ masks
                        for mask in interiors
                    ]
                )
                public = np.array(
                    [
                        np.column_stack((np.packbits(mask, axis=1), c))
                        for mask, c in zip(interiors, cut, strict=True)
                    ]
                )
                labels = np.column_stack(
                    [
                        np.unique(public[:, i], axis=0, return_inverse=True)[1]
                        for i in range(len(hops))
                    ]
                )
                same = np.logical_and.reduce(
                    [(mask == interiors[0]).all(1) for mask in interiors]
                )
                result[load, descending, r] = dict(
                    volumes=size,
                    # Rows run over initial head and columns over marker centers.
                    areas=cut[gids],
                    records=records[gids],
                    geometry_labels=labels[gids],
                    fixed_members=same,
                    fixed_boundary_count=(np.ptp(cut, axis=0) == 0).all(1),
                )
    return result
