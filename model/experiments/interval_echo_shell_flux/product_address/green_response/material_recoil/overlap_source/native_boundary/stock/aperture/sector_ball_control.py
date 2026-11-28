"""Matched native preparations and an addressed material round trip."""

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path

from . import landmark_arrival
from . import peer_acquisition
from . import sector_thermal as thermal


def contrasts(words, meter, fields, prices, subset, full, mats):
    # Explicit examples chosen by searching sharp records at fixed ordinary
    # totals. Validate every equality and acquire the chosen words afresh.
    low = np.array([1, 164])
    high = np.array([43, 198])
    open_word, closed_word = 0, 160
    ids = np.array(meter["word_to_spectral_class"])
    home = words[:, :4]
    head = words[:, -1]
    H = thermal.apparatus()[-1][0].diagonal(axis1=1, axis2=2).T
    np.testing.assert_array_equal(home[high].sum(0), home[low].sum(0))
    np.testing.assert_array_equal(H[high].sum(0), H[low].sum(0))
    np.testing.assert_array_equal(np.sort(head[high]), np.sort(head[low]))
    np.testing.assert_array_equal(home[open_word], home[closed_word])
    np.testing.assert_array_equal(H[open_word], H[closed_word])
    assert ids[open_word] == ids[closed_word]
    assert np.all(words[open_word, :6] == words[closed_word, :6])
    needed = np.unique([*low, *high, open_word, closed_word, 18, 19])
    native = np.full(len(full), -1)
    native[subset] = np.arange(len(subset))
    rng = np.random.default_rng(934856)
    tickets = peer_acquisition.returns.ticket_tables(mats)
    found = np.zeros(len(words), int)
    tries = 0
    held = {i: [] for i in needed}
    while min(found[needed]) < 512:
        ends, _ = peer_acquisition.acquire(rng, 65536, subset, native >= 0, tickets)
        for i in needed:
            held[i].extend(ends[native[ends] == i][: max(0, 512 - len(held[i]))])
        found += np.bincount(native[ends], minlength=len(words))
        tries += len(ends)

    # The full price uses the independent ordinary-clock intervals for each
    # measured group. Common additive reference normalization cancels.
    def empirical(pair):
        slots = native[np.concatenate([held[i] for i in pair])]
        return np.bincount(slots, minlength=len(words)) / len(slots)

    weights = empirical(high) - empirical(low)

    def bounds(d):
        pos = np.maximum(d, 0)
        neg = np.minimum(d, 0)
        return (pos @ prices + neg @ prices[:, ::-1]).tolist()

    board, marks, tables, _ = landmark_arrival.pulse_inputs()
    n = board.n
    ordinary = fields[2][head[closed_word]][0]
    opened = fields[2][head[open_word]][0]
    dist_to_contact = ordinary[:, np.ravel(marks)].min(1)
    valid = (dist_to_contact[:, None] > 0) & (dist_to_contact[None, :] > 0)
    # Largest shortening, then shortest outward route, then station index.
    change = ordinary - opened
    pairs = np.argwhere(valid & (change == change[valid].max()))
    start, target = min(pairs, key=lambda v: (opened[tuple(v)], *v))
    _, pred = shortest_path(
        csr_matrix(opened == 1),
        indices=start,
        return_predecessors=True,
        unweighted=True,
    )
    route = [int(target)]
    while route[-1] != start:
        route.append(int(pred[route[-1]]))
    route = route[::-1] + route[1:]
    assert len(route) - 1 == 2 * opened[start, target]
    slots = landmark_arrival.passage.contacts(board, [])[1]
    checks = {}
    for load, table in tables.items():
        result = {}
        for word in (open_word, closed_word):
            trials = 512
            launched = int(start) + n * full[held[word], -1]
            site = launched.copy()
            reached = np.zeros(trials, bool)
            for j, (u, v) in enumerate(zip(route[:-1], route[1:], strict=True)):
                slot = int(np.flatnonzero(slots[u] == v)[0])
                site = table[site, 4 * slot + rng.integers(4, size=trials)]
                if j + 1 == (len(route) - 1) // 2:
                    reached = site % n == target
            wins = int(sum(reached & (site == launched)))
            assert wins == (trials if word == open_word else 0)
            d = fields[load][head[word]][0]
            result[str(word)] = dict(
                returned_via_target=wins,
                trials=trials,
                minimum_round_trip=int(2 * d[start, target]),
            )
        checks[str(load)] = result
    return dict(
        conditioned_native_returns=tries,
        ordinary_move_input=full[held[18][0]].tolist(),
        accepted_by_word={str(i): int(found[i]) for i in needed},
        equal_geometry_pair=dict(
            lower=low.tolist(),
            higher=high.tolist(),
            mean_charge_change_band=bounds(weights),
            sum_four_home_loads=home[high].sum(0).tolist(),
            sum_four_times_contact_means=H[high].sum(0).tolist(),
            heads=head[high].tolist(),
        ),
        equal_charge_pair=dict(
            open_word=open_word,
            closed_word=closed_word,
            spectral_id=int(ids[open_word]),
            four_times_contact_means=H[open_word].tolist(),
            route=route,
            checks=checks,
        ),
    )
