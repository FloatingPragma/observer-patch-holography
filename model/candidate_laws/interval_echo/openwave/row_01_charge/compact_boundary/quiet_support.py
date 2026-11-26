"""Exact support of one bridge insertion after an arbitrary parked flight.

For every exterior input, E embeds the compact source (or its empty host).
With P the projector onto the starting interior defect count, the three
integer identities H0 E = E K, [P,H0] = 0 and P V E = 0 put the parked
packet and its first bridge variation in orthogonal number sectors for all
times. QuietBoundary.lean turns that disjointness into zero first derivative
of each sharp number record and each fixed coarsening.
"""

import json
from itertools import product

import numpy as np

from . import measure


def certify(source_sign, probe_sign, strength):
    outside = measure.PORTS[12:]
    signs = tuple(sorted(s for s in (source_sign, probe_sign) if s))
    index, parked = measure.sector(signs, strength, False)
    open_index, opened = measure.sector(signs, strength, True)
    if index != open_index:
        raise ValueError("Incompatible bridge sectors")
    if any(not np.issubdtype(h.dtype, np.integer) for h in (parked, opened)):
        raise ValueError("The source certificates require integer eighth-energy units")

    _, _, source = measure.geometry()
    singles, hsingle = measure.sector((probe_sign,), strength, False)
    ext = [singles[measure.word(((p, probe_sign),))] for p in outside]
    k = hsingle[ext, :][:, ext].toarray()
    e = np.zeros((len(index), len(outside)), dtype=int)
    for j, port in enumerate(outside):
        probe = (port, probe_sign)
        if source_sign:
            for x in np.flatnonzero(source):
                e[index[measure.word(((x, source_sign), probe))], j] = source[x]
        else:
            e[index[measure.word((probe,))], j] = 1
    if source_sign:
        k += measure.compact_energy(source_sign, strength) * np.eye(len(ext), dtype=int)
    np.testing.assert_array_equal(parked @ e, e @ k)
    gram = e.T @ e
    assert gram[0, 0] > 0
    np.testing.assert_array_equal(gram, gram[0, 0] * np.eye(len(ext), dtype=int))

    # Count complete physical number words. Keep the emitter's host case at zero.
    inner_count = np.empty(len(index), dtype=int)
    for word, i in index.items():
        inner_count[i] = sum(p < 12 for p, _ in word)
    p = inner_count == int(bool(source_sign))
    assert not np.any(e[~p])
    transitions = parked.tocoo()
    nonzero = transitions.data != 0
    np.testing.assert_array_equal(
        inner_count[transitions.row[nonzero]], inner_count[transitions.col[nonzero]]
    )
    bridge = opened - parked
    injected = bridge @ e
    assert not np.any(injected[p])
    used = np.flatnonzero(np.any(injected, axis=1))
    np.testing.assert_array_equal(
        np.unique(inner_count[used]), [int(bool(source_sign)) + 1]
    )
    # A sharp emitter sitting at the open bridge gives a same-count first
    # variation when the two opposite physical loads swap across that bridge.
    mimic = None
    if source_sign == -probe_sign:
        inside, outer = measure.BRIDGE
        word = measure.word(((inside, source_sign), (outer, probe_sign)))
        exchanged = measure.word(((inside, probe_sign), (outer, source_sign)))
        i, j = index[exchanged], index[word]
        mimic = int(bridge[i, j])
        assert p[i] and p[j] and mimic != 0
    return dict(
        source=source_sign,
        probe=probe_sign,
        strength=strength,
        number_words=len(index),
        exterior_columns=len(outside),
        first_variation_seed_rows=len(used),
        inner_count_after_one_insertion=int(inner_count[used[0]]),
        sharp_bridge_same_count_exchange_eighths=mimic,
    )


if __name__ == "__main__":
    rows = [certify(ss, ps, k) for k, ss, ps in product((8, 0), (-1, 0, 1), (-1, 1))]
    print(json.dumps(rows, indent=2))
