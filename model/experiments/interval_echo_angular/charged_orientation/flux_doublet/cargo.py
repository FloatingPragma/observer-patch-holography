"""Invariant cargo sectors for a blank-cell vacancy carrying any binary load pattern."""

from collections import Counter
from math import factorial, prod

import numpy as np

from model.executable_model import cell
from model.experiments.interval_echo_gauge.loop_transport import holonomy


def conjugacy_types(n, minimum=1):
    if n == 0:
        yield ()
    for part in range(minimum, n + 1):
        for rest in conjugacy_types(n - part, part):
            yield (part, *rest)


def characters(n, cycles):
    # A subset is fixed exactly when it is a union of entire permutation cycles.
    fixed = np.array([1], dtype=int)
    for length in cycles:
        factor = np.zeros(length + 1, dtype=int)
        factor[[0, -1]] = 1
        fixed = np.convolve(fixed, factor)
    return np.diff(fixed[: n // 2 + 1], prepend=0)


def census():
    board = cell()
    root = 0
    seen = {1}
    transpositions = []
    for a, b in board.seams:
        if root in (a, b):
            continue
        tip = min(set(board.neighbours[a]) & set(board.neighbours[b]))
        route = holonomy.path(board.neighbours, root, tip, (a, b))
        word = (
            [*route, a, b, tip, *route[-2::-1]]
            if len(route) > 1
            else [root, a, b, root]
        )
        perm = holonomy.move(board, word)
        wanted = np.arange(board.n)
        wanted[[a, b]] = wanted[[b, a]]
        np.testing.assert_equal(perm, wanted)
        transpositions.append((a, b))
    for _ in range(board.n):
        for a, b in transpositions:
            if a in seen or b in seen:
                seen.update((a, b))
    assert seen == set(range(1, board.n))
    n = board.n - 1
    gram = np.zeros((n // 2 + 1, n // 2 + 1), dtype=np.int64)
    count = 0
    for kind in conjugacy_types(n):
        centralizer = prod(
            length**k * factorial(k) for length, k in Counter(kind).items()
        )
        size = factorial(n) // centralizer
        diff = characters(n, kind)
        gram += size * np.outer(diff, diff)
        count += 1
    np.testing.assert_equal(gram, factorial(n) * np.eye(n // 2 + 1, dtype=np.int64))
    dims = characters(n, (1,) * n).tolist()
    return dict(
        transposition_edges=transpositions,
        conjugacy_classes=count,
        exact_gram_numerator=gram.tolist(),
        group_order=factorial(n),
        constituent_dimensions=dims,
        by_double_count=[dims[: min(k, n - k) + 1] for k in range(n + 1)],
    )
