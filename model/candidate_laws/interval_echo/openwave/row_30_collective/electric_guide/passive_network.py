"""Check the short-time dimer inventory read on a network of native contacts."""

import json
from itertools import product

from .connector_guard import exact_columns


def check_network(dimers):
    # Every physical crossing between neighboring dimers appears once.
    seams = [(2 * i + a, 2 * i + 2 + b)
             for i in range(dimers - 1) for a, b in product(range(2), repeat=2)]
    words = [tuple(load for sign in signs for load in (sign, 2 - sign))
             for signs in product((0, 2), repeat=dimers)]
    weights = range(1, len(seams) + 1)
    panels = 0
    for blanks in product((False, True), repeat=len(seams)):
        total = 0
        for word in words:
            # Twice the native swap amplitude is an integer. Add amplitudes
            # onto their actual terminal word before squaring.
            exits = {}
            for (i, j), weight, blank in zip(seams, weights, blanks, strict=True):
                if blank and word[i] != word[j]:
                    target = list(word)
                    target[i], target[j] = target[j], target[i]
                    if all(target[k] + target[k + 1] == 2 for k in range(0, len(target), 2)):
                        raise ArithmeticError("Crossing failed to trip the inventory meter")
                    target = tuple(target)
                    exits[target] = exits.get(target, 0) - weight
            total += sum(amplitude**2 for amplitude in exits.values())
        expected = sum(w**2 for w, blank in zip(weights, blanks, strict=True) if blank)
        if 2 * total != len(words) * expected:
            raise ArithmeticError("The terminal inventory coefficient disagrees")
        panels += 1
    return {"dimers": dimers, "seams": len(seams), "sharp_inputs": len(words),
            "memory_patterns": panels}


if __name__ == "__main__":
    print(json.dumps({"native_columns": exact_columns(),
                      "networks": [check_network(2), check_network(3)]}, indent=2))
