"""Pull the coarse mark effect through the installed total-two priority sweep."""

from collections import Counter
from math import comb
import json

import numpy as np

from model.candidate_laws.interval_echo.candidate import native_law
from model.executable_model import federation
from ..asynchronous import priority


def native_mark_gram():
    """Integer 2 W* M W, using the actual signed label and mixed-coin columns."""
    w = native_law._write_pattern()
    selected = np.array(
        [
            abs(d) == 2 and out == (1, 1)
            for coin, d in native_law.LABELS
            for out in native_law.SRC
        ]
    )
    gram = w.T @ (selected[:, None] * w)
    expected = np.diag([2 * int(abs(a - b) == 2) for a, b in native_law.SRC])
    np.testing.assert_array_equal(gram, expected)
    marked = {}
    for key in ((2, 0), (0, 2)):
        column = w[:, native_law.SRC.index(key)]
        support = np.flatnonzero(column)
        if (
            not all(selected[s] and column[s] == 1 for s in support)
            or len(support) != 2
        ):
            raise ArithmeticError(
                "A strict write lost a native coin or changed its material output"
            )
        marked[key] = tuple(native_law.LABELS[s // native_law.N_SRC] for s in support)
    return marked


def certify(board, marked):
    # Equal priority on disjoint edges is harmless. Visit whole priority
    # layers and assign a previously unwritten double to its first hit.
    first = {}
    for edge in sorted(board.seams, key=priority):
        for port in edge:
            first.setdefault(port, edge)
    if len(first) != board.n:
        raise ArithmeticError("Complete sweep missed an isolated double")
    rows = set()
    for x, edge in first.items():
        others = [tuple(sorted((x, y))) for y in board.neighbours[x]]
        ranks = list(map(priority, others))
        if len(set(ranks)) != len(ranks) or priority(edge) != min(ranks):
            raise ArithmeticError(
                "Published address disagrees with installed strict priority"
            )
        # A final marked-output row contains the winning edge, two unit
        # loads at its endpoints and the native coin and signed label.
        key = (2, 0) if x == edge[0] else (0, 2)
        for coin, difference in marked[key]:
            row = (edge, coin, difference)
            if row in rows:
                raise ArithmeticError(
                    "Two origin columns interfere in the full mark effect"
                )
            rows.add(row)
    # All doubled columns have two amplitudes 1/sqrt(2), disjoint across
    # sources. Every distinct-site binary input is quiet on every edge
    # by native_mark_gram and the native inverse-echo quiet guard.
    sizes = Counter(Counter(first.values()).values())
    return dict(
        cells=board.cells,
        blank_binary_columns=comb(board.n, 2),
        doubled_columns=len(first),
        signed_mark_rows=len(rows),
        one_origin_winners=sizes[1],
        two_origin_winners=sizes[2],
    )


def replay():
    marked = native_mark_gram()
    return [certify(federation(level), marked) for level in range(3)]


if __name__ == "__main__":
    print(json.dumps(replay(), indent=2))
