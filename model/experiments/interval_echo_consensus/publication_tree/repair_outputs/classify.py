"""Exactly readable source partitions on every enrolled subset of the nine inputs."""

from itertools import product

from .apparatus import INPUTS, possible


def partitions():
    """The union-find partition is tested independently against all binary questions."""
    observations = [set(possible(pair)) for pair in INPUTS]
    choices = list(product((0, 1), repeat=3))
    result = []
    for mask in range(1, 1 << 9):
        labels = [i for i in range(9) if mask & (1 << i)]
        parts = []
        for i in labels:
            touched = [
                p for p in parts if any(observations[i] & observations[j] for j in p)
            ]
            parts = [p for p in parts if p not in touched]
            parts.append({i}.union(*touched))
        parts.sort(key=min)
        # An arbitrary deterministic binary function of each single-site ternary output.
        # Its response must be sharp on every included source. Count resulting distinct
        # source functions, discarding changes on impossible output symbols.
        answers = set()
        for decoder in choices:
            if all(len({decoder[o] for o in observations[i]}) == 1 for i in labels):
                answers.add(tuple(decoder[min(observations[i])] for i in labels))
        if len(answers) != 2 ** len(parts):
            raise ArithmeticError(
                "component partition differs from the source-question enumeration"
            )
        result.append([mask, [sum(1 << j for j in p) for p in parts]])
    return result
