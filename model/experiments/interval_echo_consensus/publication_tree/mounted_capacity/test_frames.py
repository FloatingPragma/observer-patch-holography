"""Independent checks of all addresses and of the public frame comparisons."""

import gzip
import hashlib
import json
from collections import defaultdict


from .apparatus import AUDITORS, MARKS, WORDS, Board
from .measure import HERE, ROOT


def archive(name):
    return json.loads(gzip.decompress((HERE / name).read_bytes()))


def clicks(digits, launch):
    return [tuple((int(v, 16) >> launch) & 1 for v in pair) for pair in digits]


def test_every_rotation_and_record_at_both_observers():
    d = archive("records.json.gz")
    for path, h in d["sources"].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == h
    assert len({tuple(r["mount"]) for r in d["receipts"]}) == 60
    for r in d["receipts"]:
        g = r["mount"]
        for readings in r["raw"]:
            for i, j in enumerate(g):
                a, b = [int(x, 16) for x in readings[j]]
                # The truth table of input bit i is alternating zero/one runs.
                word = sum(
                    (((1 << (1 << i)) - 1) << (k + (1 << i)))
                    for k in range(0, WORDS, 2 << i)
                )
                assert b == word and a ^ b == (1 << WORDS) - 1
                for column, mark in enumerate(MARKS, 1):
                    assert [
                        int(x, 16).bit_count() for x in readings[12 * column + j]
                    ] == ([0, WORDS] if i == mark else [WORDS, 0])


def test_original_scalar_echo_and_route_on_held_launch():
    d = archive("records.json.gz")
    board = Board()
    for row in (d["receipts"][7], d["receipts"][38]):
        g = row["mount"]
        word = 0b101001011010
        loads = [0] * board.device.car.n
        for dim in board.device.rows:
            loads[dim[2][0]] = 2
        for offset, value in enumerate((word, 1 << MARKS[0], 1 << MARKS[1])):
            for i in range(12):
                a, b = board.slots[12 * offset + i]
                bit = (value >> i) & 1
                loads[a], loads[b] = 2 * bit, 2 * (1 - bit)
        for a, b in board.program(tuple(g), False):
            loads[a], loads[b] = loads[b], loads[a]
        for k, station in enumerate(AUDITORS):
            row_s = board.device.rows[station]
            read = [
                tuple(
                    board.device.probe(loads, (rail, row_s[2][0]), row_s[3], "native")
                    for rail in pair
                )
                for pair in board.slots
            ]
            assert read == clicks(row["raw"][k], word)


def test_veto_on_shared_records_and_all_ambiguity_counts():
    d = archive("records.json.gz")
    v = archive("veto.json.gz")
    for name, h in v["sources"].items():
        assert hashlib.sha256((HERE / name).read_bytes()).hexdigest() == h
    source_bit = sum(1 << i for i in range(WORDS) if i & (1 << MARKS[0]))
    for row in v["rows"]:
        assert all(
            int(r[0][1], 16) == source_bit
            and int(r[1][1], 16) == ((1 << WORDS) - 1) ^ source_bit
            for r in row["output"]
        )
    assert v["joint_vetoes"] == WORDS * 30
    assert v["controls"]["changed"] == 2048
    assert all(int(x, 16) == 0 for p in v["controls"]["parked"]["survey"] for x in p)
    # Independently form stabilizer orbits; classify the physical host values in each orbit.
    for menu, m in enumerate(("unmarked", "host", "edge", "ordered")):
        counts = defaultdict(int)
        gset = []
        for r in d["receipts"]:
            g = r["mount"]
            a, b = MARKS
            if [True, g[a] == a, {g[a], g[b]} == {a, b}, (g[a], g[b]) == (a, b)][menu]:
                gset.append(g)
        remaining = set(range(WORDS))
        while remaining:
            w = min(remaining)
            orbit = {sum(((w >> i) & 1) << j for i, j in enumerate(g)) for g in gset}
            remaining.difference_update(orbit)
            bs = {bool(other & (1 << MARKS[0])) for other in orbit}
            key = "ambiguous" if len(bs) == 2 else ("veto" if True in bs else "allow")
            counts[key] += 1
        assert {k: counts[k] for k in v["census"][m]} == v["census"][m]
        assert sum(counts.values()) == d["sizes"][m]["code_size"]
