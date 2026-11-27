"""Acquire and independently weigh the complete native sharp branch banks."""

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import gzip
import hashlib
from itertools import product
import json
from math import prod
from pathlib import Path

import numpy as np

from .instrument import BOUNDARY, MASKS, NOMINAL, acquire, pulse_calibration

HERE = Path(__file__).resolve().parent
ALL = (1 << 6) - 1
WEIGHTS = (12, 3, 1)


def flags(side, m):
    return sum(1 << (3 * side + i) for i in MASKS[m])


def one(q, mask, mode):
    return (q & ~mask, q, q | mask)[mode]


def joint(q, order, held=False, reverse=False):
    if held:
        out = Counter()
        for mode, weight in enumerate(WEIGHTS):
            state = q
            for side, m in order:
                state = one(state, flags(side, m), mode)
            out[state] += weight
        return out, 16
    states = {q: 1}
    for side, m in order:
        row = defaultdict(int)
        for old, w in states.items():
            for mode, n in enumerate(WEIGHTS):
                key = 3 - mode if reverse and mode else mode
                row[one(old, flags(side, m), key)] += w * n
        states = row
    return states, 16 ** len(order)


def weights(length):
    return np.array(
        [prod(WEIGHTS[m] for m in modes) for modes in product(range(3), repeat=length)],
        dtype=np.int64,
    )


def unpack(value, n):
    width = (n + 7) // 8
    return np.unpackbits(
        np.frombuffer(int(value, 16).to_bytes(width, "little"), dtype=np.uint8),
        bitorder="little",
    )[:n].astype(np.int64)


def extract_state(snapshot, n):
    state = np.zeros(n, dtype=np.int64)
    for side, bank in enumerate(snapshot):
        b = unpack(bank[0], n)
        if not np.all(b == BOUNDARY[side]):
            raise ArithmeticError("Protected boundary changed")
        decoded = [unpack(v, n) ^ b for v in bank[1:4]]
        rank = sum(decoded)
        for k, selector in enumerate(bank[4:]):
            if not np.array_equal(unpack(selector, n), rank == k):
                raise ArithmeticError("Wrong measured rank address")
        for k, mismatch in enumerate(decoded):
            state += mismatch << (3 * side + k)
    return state


def census(bank):
    n = bank["leaves"]
    w = weights(1 if bank["held"] else len(bank["order"]))
    rows = []
    for snap in bank["state"]:
        states = extract_state(snap, n)
        hist = [int(w[states == q].sum()) for q in range(64)]
        rows.append(hist)
    q = next(i for i, m in enumerate(rows[0]) if m)
    total = int(w.sum())
    for t, row in enumerate(rows):
        if not t:
            continue
        expected, den = joint(
            q, bank["order"][:t], bank["held"], bank["reverse_faults"]
        )
        if any(row[v] * den != expected.get(v, 0) * total for v in range(64)):
            raise ArithmeticError(
                "Routed joint record differs from independently composed kernel"
            )
    return dict(
        denominator=total,
        joint_rows=rows,
        expected_distances=[
            str(Fraction(sum(q.bit_count() * v for q, v in enumerate(row)), total))
            for row in rows
        ],
    )


def certificate():
    # Uniform over every order. At fixed block start a prefix records a subset
    # of eight labelled offers; only the number touching each address affects its mean.
    offers = [flags(s, m) for s, m in NOMINAL]
    rows = []
    largest = 0
    for sub in range(256):
        hits = [
            sum(
                bool(mask & (1 << a)) for j, mask in enumerate(offers) if sub & (1 << j)
            )
            for a in range(6)
        ]
        with_error = []
        for q in range(64):
            distance = Fraction()
            for i, k in enumerate(hits):
                f = Fraction((q >> i) & 1)
                for _ in range(k):
                    f = (3 * f + 1) / 16
                distance += f
            scaled = distance * 4096
            if scaled.denominator != 1:
                raise ArithmeticError("Distance grid failed")
            value = int(scaled)
            with_error.append(value)
            if sub != 255:
                largest = max(largest, value - 4096 * q.bit_count())
            elif value != 27 * q.bit_count() + 1878:
                raise ArithmeticError("Block inequality failed")
        rows.append(dict(offers=sub, hits=hits, expected_4096=with_error))
    if largest != 1860:
        raise ArithmeticError("Intermediate sharp constant failed")
    return dict(
        active_states=256,
        normal_states=4,
        block_orders=40320,
        subsets=rows,
        distance_denominator=4096,
        largest_strict_prefix_excess=largest,
    )


def report():
    wrong = ((1, 1, 1), (0, 0, 0))
    normal = ((0, 0, 0), (1, 1, 1))
    rev = tuple(reversed(NOMINAL))
    banks = [acquire(o, q) for o in (NOMINAL, rev) for q in (normal, wrong)]
    banks.extend(
        [
            acquire(NOMINAL, wrong, held=True),
            acquire(tuple((0, m) for _ in (0, 1) for m in range(4)), wrong),
            acquire(NOMINAL, normal, reverse_faults=True),
        ]
    )
    for bank in banks:
        bank["receiver_census"] = census(bank)
    margin = [b["receiver_census"]["expected_distances"][-1] for b in banks]
    if margin[:5] != ["939/2048", "255/512", "939/2048", "255/512", "3/2"]:
        raise ArithmeticError("Prospective terminal expectations failed")
    from .. import apparatus as triangle
    from ....seam_word_meter import apparatus as seam
    from ....prepared_quorum import transport
    from ....boolean_nodes import material

    root = HERE.parents[5]
    paths = [
        HERE / "PROTOCOL.md",
        HERE / "instrument.py",
        Path(__file__).resolve(),
        *[
            Path(m.__file__).resolve()
            for m in (triangle, seam, transport, material, seam.LAW)
        ],
    ]
    return dict(
        provenance="Complete conditional native outcome branches and exact weights, without sampled launches.",
        branch_order="lexicographic, correct=0, preserve=1, wrong=2",
        selector=pulse_calibration(),
        source_sha256={
            str(f.relative_to(root)): hashlib.sha256(f.read_bytes()).hexdigest()
            for f in paths
        },
        certificate=certificate(),
        banks=banks,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    options = parser.parse_args()
    data = (json.dumps(report(), separators=(",", ":")) + "\n").encode()
    if options.output:
        with options.output.open("xb") as f:
            f.write(gzip.compress(data, mtime=0))
    elif gzip.decompress((HERE / "records.json.gz").read_bytes()) != data:
        raise ValueError("Endpoint or intermediate record differs")
