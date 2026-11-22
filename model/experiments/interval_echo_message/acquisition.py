"""Trial enrolment, detector sampling and scoring outside the receiver."""
import math

import numpy as np

from .channel import SPEC, require, validate_law
from .decoder import decode


def stream(role, level, grid, cohort=0):
    return np.random.Generator(np.random.PCG64([SPEC["seeds"][role], level, grid, cohort]))


def outcomes(p, bits, rng):
    validate_law(p)
    u = rng.random(len(bits))
    records = np.empty(len(bits), dtype=np.uint8)
    for bit in (0, 1):
        mask = bits == bit
        # Half-open port intervals. The final interval is the complete quiet branch.
        records[mask] = np.searchsorted(np.cumsum(p[bit, :-1]), u[mask], side="right")
    return records


def counts(bits, records, per_bit):
    require(bits.ndim == records.ndim == 1 and len(bits) == len(records) == 2 * per_bit
            and bits.dtype.kind in "iu" and records.dtype.kind in "iu"
            and bits.min() >= 0 and bits.max() <= 1 and records.min() >= 0
            and records.max() <= 12, "missing or invalid enrolled trial")
    table = np.bincount(bits.astype(int) * 13 + records, minlength=26).reshape(2, 13)
    require(np.all(table.sum(axis=1) == per_bit), "lost balanced enrolment")
    return table


def train(p, level, grid):
    n = SPEC["calibration_per_bit"]
    bits = stream("setup", level, grid).permutation(np.repeat([0, 1], n))
    local = outcomes(p, bits, stream("calibration", level, grid))
    c = counts(bits, local, n)
    return c.tolist(), (c[1] > c[0]).astype(np.uint8).tobytes()


def accuracy(correct, n, permutation=False):
    radius = math.sqrt(math.log(2 * SPEC["comparisons"] / SPEC["alpha"])
                       / (n if permutation else 2 * n))
    rate = correct / n
    return dict(correct=correct, total=n, accuracy=rate,
                interval=[max(0., rate - radius), min(1., rate + radius)])


def evaluate(p, lookup, level, grid, cohort):
    n = SPEC["trials"]
    bits = stream("bits", level, grid, cohort).permutation(np.repeat([0, 1], n // 2))
    records = outcomes(p, bits, stream("detector", level, grid, cohort))
    guesses = np.fromiter((decode(int(y), lookup) for y in records), dtype=np.uint8, count=n)
    c = counts(bits, records, n // 2)
    right = guesses == bits
    quiet = records == 12
    row = dict(counts=c.tolist(), quiet=int(quiet.sum()), quiet_correct=int(c[lookup[-1], -1]), errors=int((~right).sum()), abstentions=0,
               **accuracy(int(right.sum()), n))
    if cohort == 0:
        shuffled = stream("shuffle", level, grid).permutation(bits)
        row["permuted"] = accuracy(int((guesses == shuffled).sum()), n, permutation=True)
    return row
