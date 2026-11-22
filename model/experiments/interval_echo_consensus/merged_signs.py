"""Remove parallel sign-flip destinations on the compact bank, keeping all material moves."""
import json
from collections import defaultdict
from itertools import combinations

import numpy as np
from scipy import sparse

from model.executable_model import carrier, federation
from .record_release import pair, source
from .sign_graph import questions


def compare(arm):
    board = federation(0)
    offers, _ = questions(board, arm)
    by_marker = defaultdict(list)
    for a, b, sign in offers:
        by_marker[sign].append((a, b))
    positions = np.array(list(combinations(range(board.n), 2)))
    i, j = positions.T
    _, index, exchange, adjacency = pair.screen(carrier.HERE / "carriers/FEDERATION_LEVEL0.npz")
    _, launch, patch, _ = source(arm)
    penalties = np.zeros(len(i), dtype=int)
    launch_multiplicity = 0
    for checks in by_marker.values():
        active = np.zeros(len(i), dtype=int)
        for a, b in checks:
            at_b = False if b is None else (i == b) | (j == b)
            active += ((i == a) | (j == a)) ^ at_b
        penalties += active != 0
        launch_multiplicity = max(launch_multiplicity, int(active[launch]))
    if not np.array_equal(index[i, j], np.arange(len(i))):
        raise ValueError("Wrong input order")
    h = exchange + sparse.diags(penalties)
    readings = []
    for t, state in zip(pair.EXPOSURES, pair.evolve(h, launch, pair.EXPOSURES), strict=True):
        p = pair.local_law(state, patch)
        readings.append(dict(exposure=t, probabilities=p.tolist(),
                             adjacent=float((abs(state[adjacency.astype(bool)])**2).sum())))
    return dict(arm=arm, largest_launch_multiplicity=launch_multiplicity, settings=readings)


if __name__ == "__main__":
    print(json.dumps([compare(arm) for arm in ("primitive", "assembled")], allow_nan=False))
