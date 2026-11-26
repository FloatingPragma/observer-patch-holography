"""Meter the one-vacancy inlet before the quiet collar writes seal the fork."""

from collections import deque

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from .material import edge, mount


def inlet():
    board, dimers, heads, _ = mount()
    p = dimers[0][0]
    involved = set(sum(dimers, ()))
    # Leave the bank outside both measured owners and away from each quiet write.
    near = involved | {q for site in involved for q in board.neighbours[site]}
    owners = {p // 12 for p in involved}
    pending = deque([[p]])
    visited = {p}
    while pending:
        path = pending.popleft()
        tip = path[-1]
        if tip not in near and tip // 12 not in owners:
            break
        for q in sorted(board.neighbours[tip]):
            if q not in visited and edge(tip, q) not in heads:
                pending.append(path + [q])
                visited.add(q)
    else:
        raise ValueError("No accessible vacancy bank")
    donor = path[-1]
    loaded = np.ones(board.n, dtype=int)
    loaded[donor] = 0
    loaded[[x[0] for x in dimers[1:]]] = 0
    sham = loaded.copy()
    forward = list(zip(path, path[1:]))
    for a, b in forward + forward[-2::-1]:
        loaded[a], loaded[b] = loaded[b], loaded[a]
    np.testing.assert_array_equal(np.flatnonzero(loaded - sham), sorted([p, donor]))
    # On all binary load words the full native blank pulse is an exact swap.
    law = candidate.native_law
    ex, ct = law.local_term_parts()
    u = expm(-1j * np.pi * (ex + candidate.KAPPA * ct))
    for a in (0, 1):
        for b in (0, 1):
            col = u[:, law.SRC.index((a, b))]
            np.testing.assert_allclose(
                col, np.eye(law.DIM)[:, law.SRC.index((b, a))], atol=2e-12
            )

    def read(x):
        return [*(int(x[list(pair)].sum()) for pair in dimers), int(x[donor])]

    return dict(
        path=path,
        words_before=read(sham),
        words_after=read(loaded),
        total=int(sum(loaded)),
    )


if __name__ == "__main__":
    print(inlet())
