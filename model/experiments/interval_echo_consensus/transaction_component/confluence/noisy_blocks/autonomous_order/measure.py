"""Regenerate the exact sharp receiver tapes and compare their joint dynamics."""

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import gzip
import json

from ..instrument import NOMINAL
from ..measure import flags, one, extract_state
from . import dispatcher, scalar

HERE = Path(__file__).parent
WEIGHTS = (12, 3, 1)


def expected_label(q, mask):
    active = [j for j in range(8) if (mask >> j) & 1]
    return (max if q.bit_count() % 2 else min)(active)


@lru_cache(maxsize=1024)
def column(h):
    return int(h, 16)


def word(columns, leaf):
    return sum(((column(h) >> leaf) & 1) << j for j, h in enumerate(columns))


def decode(bank):
    flags_by_step = [extract_state(s["state"], bank["leaves"]) for s in bank["steps"]]
    terminals = Counter()
    routes = Counter()
    sums = [0] * 8
    for leaf, noise in enumerate(product(range(3), repeat=8)):
        q, mask = 8, 255
        weight = prod(WEIGHTS[m] for m in noise)
        route = []
        for tick, step in enumerate(bank["steps"]):
            j = expected_label(q, mask)
            assert word(step["flags"], leaf) == q
            assert (
                word(step["before"], leaf) == mask
                and word(step["groups"], leaf) == 1 << j
            )
            assert word(step["fault"], leaf) == 1 << noise[tick]
            mask ^= 1 << j
            assert word(step["retired"], leaf) == mask
            assert word(step["final"], leaf) == (mask if mask else 511)
            q = one(q, flags(*NOMINAL[j]), noise[tick])
            assert flags_by_step[tick][leaf] == q
            route.append(j)
            sums[tick] += weight * q.bit_count()
        routes[tuple(route)] += weight
        terminals[q] += weight
    assert sum(terminals.values()) == 16**8
    return {
        "branch_total": 3**8,
        "integer_weight": 16**8,
        "distinct_orders": len(routes),
        "terminal_weights": [terminals[i] for i in range(64)],
        "prefix_means": [str(Fraction(v, 16**8)) for v in sums],
    }


def validate(source):
    for odd, table in enumerate(source["mask_probes"]):
        row = table["proposal"]
        for m in range(256):
            assert word(row["before"], m) == m
            j = expected_label(odd, m) if m else None
            selection = 0 if j is None else 1 << j
            assert word(row["groups"], m) == selection
            assert word(row["retired"], m) == m ^ selection
    s = source["parked"]["steps"]
    assert len(s) == 1 and s[0]["status"] == "refused" and "fault" not in s[0]
    assert s[0]["before"] == s[0]["retired"]
    assert source["parked"]["initial"] == source["parked"]["ending"]
    scalar_bank = source["scalar"]["stages"]
    single = source["two_blocks"]["steps"]
    for tick, (s, row) in enumerate(zip(scalar_bank, single, strict=True)):
        for key, raw_key in (
            ("mask", "before"),
            ("retired", "retired"),
            ("final", "final"),
            ("flags", "flags"),
        ):
            assert sum(v << i for i, v in enumerate(s[key])) == word(row[raw_key], 0)
        assert 1 << s["label"] == word(row["groups"], 0) and s["mode"] == tick % 3
        assert (
            sum(v << i for i, v in enumerate(s["flags_after"]))
            == extract_state(row["state"], 1)[0]
        )
    result = decode(source["tree"])
    assert result["prefix_means"][-1] == str(Fraction(27 + 1878, 4096))
    assert result["distinct_orders"] >= 2
    return result


def report():
    raw = dispatcher.report()
    raw["scalar"] = scalar.acquire()
    raw["summary"] = validate(raw)
    raw["origin"] = (
        "All sharp branches of the calibrated projected native source and their exact joint weights; scalar pointer records on a two-block held trajectory; no sampled shots."
    )
    return raw


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--write", action="store_true")
    write = cli.parse_args().write
    data = report()
    encoded = (json.dumps(data, separators=(",", ":")) + "\n").encode()
    target = HERE / "records.json.gz"
    if write:
        target.write_bytes(gzip.compress(encoded, compresslevel=6, mtime=0))
    elif gzip.decompress(target.read_bytes()) != encoded:
        raise ValueError("Regenerated records differ from the archived receiver bank")
    print(json.dumps(data["summary"]))
