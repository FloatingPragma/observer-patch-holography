"""Integer graph and block bounds for subthreshold material on the registered boards."""
from fractions import Fraction
from itertools import combinations

import numpy as np

from model.executable_model import carrier


def partitions(total, largest=None):
    if not total:
        yield ()
        return
    for n in range(1, min(total, largest or total) + 1):
        for tail in partitions(total-n, n):
            yield (n, *tail)


def compactness():
    cell = carrier.cell()
    seams = set(cell.seams)
    inside = {n: max(sum(pair in seams for pair in combinations(word, 2))
                     for word in combinations(range(12), n)) for n in range(1, 7)}
    # Each crossing edge consumes two distinct endpoints and at least one
    # endpoint outside a most populated cell.
    global_max = {}
    for n in range(1, 7):
        upper = [sum(inside[k] for k in p) + min(n//2, n-max(p))
                 for p in partitions(n) if len(p) > 1]
        global_max[n] = max([inside[n], *upper])
    disconnected = max(sum(global_max[k] for k in p)
                       for p in partitions(6) if len(p) > 1)
    min_energy = 8*(36-2*global_max[6])
    dis_energy = 8*(36-2*disconnected)
    return dict(cell_max_edges=inside, federation_max_edges=global_max,
                disconnected_max_edges=disconnected, potential_lower=min_energy,
                disconnected_lower=dis_energy)


def face_band(sector, shape):
    cap, faces = shape.caps, shape.pairs
    other = ~(cap | faces)
    moves = sector.K < 0
    counts = []
    for source, targets in ((cap, ~cap), (faces, other)):
        block = moves[np.ix_(source, targets)]
        counts.append([int(block.sum(axis=1).max()), int(block.sum(axis=0).max())])
    if np.any(moves[np.ix_(cap, faces)]):
        raise ValueError("Detected direct cap to two-face exchange")
    if not np.array_equal(sector.hamiltonian(8, False)[np.ix_(cap, cap)], 133*np.eye(cap.sum())):
        raise ValueError("Cap compression changed")
    gap = int(min(8*(36-2*sector.agreed[~cap])) - 133)
    squared = Fraction(counts[0][0]*counts[0][1]*counts[1][0]*counts[1][1],
                       16 * gap**4)
    return dict(caps=int(cap.sum()), face_pairs=int(faces.sum()), block_degrees=counts,
                complement_gap=gap, two_face_weight_upper=str(squared))
