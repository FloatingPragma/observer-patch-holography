"""One fundamental in each of three disjoint receiver domains."""

from itertools import combinations

import numpy as np

from model.executable_model import carrier


def groups_on_cell():
    pool = set(range(12))
    for na, nb, nc in ((3, 3, 6), (3, 4, 5), (4, 4, 4)):
        for a in combinations(range(12), na):
            rest = pool.difference(a)
            for b in combinations(sorted(rest), nb):
                c = tuple(sorted(rest.difference(b)))
                if na == nb and a > b or nb == nc and b > c:
                    continue
                yield (a, b, c)


def canonical(groups):
    return tuple(sorted(groups, key=lambda group: (len(group), group)))


def group_orbits(cell):
    unseen = set(groups_on_cell())
    symmetries = carrier.automorphisms(cell.seams)
    while unseen:
        g = next(iter(unseen))
        orbit = {
            canonical([tuple(sorted(perm[x] for x in row)) for row in g])
            for perm in symmetries
        }
        unseen -= orbit
        yield min(orbit), len(orbit)


def maxima(cell, sector, lines):
    states = np.array(sector.states)
    peaks = np.zeros(len(lines))
    peak_groups = [None] * len(lines)
    sizes = {}
    representatives = sorted(group_orbits(cell))
    for groups, size in representatives:
        key = ",".join(str(len(g)) for g in groups)
        sizes[key] = sizes.get(key, 0) + size
        owner = np.zeros(12, int)
        for j, g in enumerate(groups):
            owner[list(g)] = j
        arr = owner[states]
        keep = (
            (arr[:, 0] != arr[:, 1])
            & (arr[:, 1] != arr[:, 2])
            & (arr[:, 0] != arr[:, 2])
        )
        for j, (_, basis) in enumerate(lines):
            cut = basis[keep]
            score = np.linalg.eigvalsh(cut.T @ cut)[-1]
            if score > peaks[j] + 1e-12:
                peaks[j], peak_groups[j] = score, groups
    return {
        "counts_by_size": sizes,
        "geometry_orbits": len(representatives),
        "energies": [e for e, _ in lines],
        "dimensions": [v.shape[1] for _, v in lines],
        "one_per_group_maxima": peaks.tolist(),
        "maximizing_groups": peak_groups,
    }


def modular_rank(matrix, prime=1009):
    a = np.asarray(matrix, dtype=np.int64).copy() % prime
    nrows, ncols = a.shape
    pivot = 0
    for col in range(ncols):
        candidates = np.flatnonzero(a[pivot:, col])
        if not len(candidates):
            continue
        idx = pivot + candidates[0]
        a[[pivot, idx]] = a[[idx, pivot]]
        a[pivot] = a[pivot] * pow(int(a[pivot, col]), -1, prime) % prime
        a[pivot + 1 :] = (a[pivot + 1 :] - a[pivot + 1 :, col, None] * a[pivot]) % prime
        pivot += 1
        if pivot == nrows:
            break
    return pivot


def transfer_certificate(cell, sector, h):
    twice = np.rint(2 * h).astype(np.int64)
    np.testing.assert_array_equal(twice, 2 * h)
    states = np.asarray(sector.states)
    records = []
    for groups, size in sorted(group_orbits(cell)):
        owner = np.zeros(12, int)
        for j, g in enumerate(groups):
            owner[list(g)] = j
        x = owner[states]
        mask = (x[:, 0] != x[:, 1]) & (x[:, 1] != x[:, 2]) & (x[:, 0] != x[:, 2])
        a = twice[np.ix_(mask, mask)]
        b = twice[np.ix_(~mask, mask)]
        b = b[b.any(axis=1)]
        rank = modular_rank(b)
        final = rank if rank == len(a) else modular_rank(np.vstack([b, b @ a]))
        records.append(
            dict(
                groups=groups, orbit=size, columns=len(a), rank_b=rank, rank_b_ba=final
            )
        )
    return records
