"""Read coincidence frequencies using the two observers' complete count banks."""

import json
from functools import cache
from itertools import product
from pathlib import Path

from scipy.stats import beta

from ..measure import OWNERS, TIMES
from .measure import MODES


@cache
def interval(success, n):
    tail = 0.01 / (2 * 4096)
    low = beta.ppf(tail, success, n - success + 1) if success else 0
    high = beta.isf(tail, success + 1, n - success) if success < n else 1
    return (float(low), float(high))


def difference(first, second):
    x, a = first
    y, b = second
    lx, hx = interval(x, a)
    ly, hy = interval(y, b)
    return dict(
        central=x / a - y / b,
        limits=[lx - hy, hx - ly],
        resolved=(lx > hy or hx < ly),
        counts=[first, second],
    )


def read(record):
    size = record["starts"]
    banks = {tuple(r["key"]): r["counts"] for r in record["held"]}
    if len(banks) != 96 or not all(
        sum(v) == size and min(v) >= 0 for v in banks.values()
    ):
        raise ValueError("missing original starts, negative counts, or duplicate bank")
    results = []
    for mode, owner, time, panel in product(MODES, OWNERS, TIMES, (-1, 1)):
        signs = [banks[owner, time, j, panel, mode] for j in (1, 2)]
        conditional = [(counts[0], sum(counts[:3])) for counts in signs]
        far = [(n, size) for _, n in conditional]
        results.append(
            dict(
                mode=mode,
                owner=owner,
                time=time,
                panel=panel,
                joint_given_far=difference(*conditional),
                far_marginal=difference(*far),
            )
        )
    return results


if __name__ == "__main__":
    data = json.loads(Path(__file__).with_name("counts.json").read_text())
    print(json.dumps(read(data), indent=1))
