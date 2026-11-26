"""Compile truth tables at enrollment; the physical interpreter only receives contact instructions."""

from itertools import product


def minterm(literals, target, scratch):
    """Return signed NOT, controlled-NOT and doubly controlled-NOT instructions."""
    if len(literals) <= 2:
        return [(tuple(literals), target)]
    if len(scratch) < len(literals) - 2:
        raise ValueError("Missing clean scratch")
    first, *rest = literals
    prefix = first
    ladder = []
    for j, literal in enumerate(rest[:-1]):
        ladder.append(((prefix, literal), scratch[j]))
        prefix = (scratch[j], 1)
    return ladder + [((prefix, rest[-1]), target)] + ladder[::-1]


def compile_table(sources, targets, table, scratch):
    """XOR independent outputs with the tabulated finite vector map."""
    layout = (*sources, *targets, *scratch)
    if any(type(r) is not int or r < 0 for r in layout):
        raise ValueError("Expected nonnegative register addresses")
    if len(set(layout)) != len(layout) or len(scratch) < max(0, len(sources) - 2):
        raise ValueError("Distinct public words and clean scratch are required")
    if len(table) != 1 << len(sources):
        raise ValueError("Expected every input row")
    code = []
    for mask, row in enumerate(table):
        if not 0 <= row < 1 << len(targets):
            raise ValueError("Output is outside the registered alphabet")
        literals = [(reg, (mask >> i) & 1) for i, reg in enumerate(sources)]
        for j, target in enumerate(targets):
            if (row >> j) & 1:
                code.extend(minterm(literals, target, scratch))
    return code


def support_inputs():
    return tuple(product(range(8), range(2)))
