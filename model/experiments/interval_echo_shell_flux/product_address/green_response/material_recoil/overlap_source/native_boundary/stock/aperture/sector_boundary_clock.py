"""Ordinary calibration by transporting a double and timing boundary exchanges."""

from itertools import product

import numpy as np
from scipy.sparse.csgraph import shortest_path

from . import exchange_geometry


def move_double(graph, homes, words, acquired):
    before = 18
    after = 19
    state = np.zeros(graph.shape[0], int)
    assert np.array_equal(acquired, words[before])
    state[homes] = acquired[:4]
    start, end = homes[0], homes[1]
    assert (
        state[start] == 2
        and state[end] == 0
        and np.all(words[before, 4:] == words[after, 4:])
    )
    excluded = state > 0
    excluded[start] = False
    adj = graph.toarray().astype(bool)
    adj[excluded] = False
    adj[:, excluded] = False
    _, pred = shortest_path(
        adj, indices=start, unweighted=True, return_predecessors=True
    )
    path = [end]
    while path[-1] != start:
        path.append(int(pred[path[-1]]))
    path = path[::-1]
    swap = exchange_geometry.contact_swap().argmax(0)
    for u, v in zip(path[:-1], path[1:], strict=True):
        state[u], state[v] = divmod(int(swap[3 * state[u] + state[v]]), 3)
    assert (
        np.array_equal(state[homes], words[after, :4])
        and state.sum() == words[after, :4].sum()
    )
    return before, after, path


def rates(graph, materials, homes, u, v):
    # Measure each link on its own fresh input. Protect every occupied
    # endpoint in a distinct empty native neighbor before inserting the
    # ordinary pulse. Restore both residents after acknowledgement.
    adj = graph.toarray().astype(bool)
    swap = exchange_geometry.contact_swap().argmax(0)
    durations = []
    for material in materials:
        state = np.zeros(len(adj), int)
        state[homes] = material
        occupied = set(np.flatnonzero(state))
        row = []
        for x, y in zip(u, v, strict=True):
            on = sorted({x, y} & occupied)
            choices = [
                sorted(set(np.flatnonzero(adj[i])) - occupied - {x, y}) for i in on
            ]
            rails = next(
                parking
                for parking in product(*choices)
                if len(set(parking)) == len(parking)
            )
            # Unoccupied links have the same exact two-pulse return regardless
            # of loads outside their endpoints. Test both loads on each bit.
            for load in (1, 2):
                trial = state.copy()

                def pulse(i, j):
                    assert adj[i, j]
                    trial[i], trial[j] = divmod(int(swap[3 * trial[i] + trial[j]]), 3)

                for i, j in zip(on, rails, strict=True):
                    pulse(i, j)
                trial[x] = load
                pulse(x, y)
                assert trial[y] == load
                pulse(y, x)
                assert trial[x] == load
                trial[x] = 0
                for i, j in zip(on, rails, strict=True):
                    pulse(i, j)
                assert np.array_equal(trial, state)
            row.append(1 + len(on))
        durations.append(row)
    return np.array(durations)
