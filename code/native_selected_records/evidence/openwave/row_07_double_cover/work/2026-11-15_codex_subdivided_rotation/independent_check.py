"""Independent full-ledger audit. No imports from any producer module."""
import json
import math
import sys
from fractions import Fraction as F
from pathlib import Path


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def mv(a, v):
    return [dot(row, v) for row in a]


def branches(edges, law):
    out = []
    for x, y in edges:
        for coin in range(2):
            flags = [[[0] * 3 for _ in range(3)] for _ in range(3)]
            for source in range(3):
                old = [int(k == source) for k in range(3)]
                new = old[:]
                if law == "native":
                    low, high = (old[x] + old[y]) // 2, (old[x] + old[y] + 1) // 2
                    new[x], new[y] = (low, high) if coin == 0 else (high, low)
                else:
                    new[x], new[y] = old[y], old[x]
                flags[old[x] - new[x] + 1][new.index(1)][source] = 1
            out += flags
    return out


def expected(slot, m, law):
    v, k = divmod(slot, m)
    edges = [(0, 1), (0, 2), (1, 2)]
    # Order matters for individual flags: producer EDGES=(01,12,02).
    ordered = [(0, 1), (1, 2), (0, 2)]
    addresses = [e for e in ordered if v in e] * (m - k)
    addresses += [e for e in ordered if (v + 1) % 3 in e] * k
    fine = branches(addresses, law)
    prefix = [[0] * 3 for _ in range(3)]
    rows = []
    for i, a in enumerate(fine):
        if i:
            rows.append(([[prefix[r][c] - i * a[r][c] for c in range(3)] for r in range(3)],
                         4 * m * i * (i + 1)))
        prefix = [[prefix[r][c] + a[r][c] for c in range(3)] for r in range(3)]
    rows.insert(0, (prefix, 48 * m * m))
    mean = [[F(x, 4 * m) for x in row] for row in prefix]
    return mean, rows


def read_ops(rows):
    return [([list(map(int, row)) for row in op["numerator"]], op["squared_denominator"]) for op in rows]


base = Path(__file__).resolve().parent
r = json.loads((base / (sys.argv[1] if len(sys.argv) > 1 else "CONFIRM_RECEIPT.json")).read_text())
assert not r["gate_failures"]
tables, instruments = {}, 0
for key, catalog in r["instrument_catalog"].items():
    m, law = int(key.split(":")[0]), key.split(":")[1]
    means, ops = [], []
    for slot, encoded in enumerate(catalog["target"]):
        mean, predicted = expected(slot, m, law)
        actual = read_ops(encoded)
        assert predicted == actual
        means.append(mean)
        ops.append(actual)
    plane = read_ops(catalog["plane_filter"])
    for outputs in [*ops, plane]:
        for i in range(3):
            for j in range(3):
                assert sum(F(sum(a[k][i] * a[k][j] for k in range(3)), d)
                           for a, d in outputs) == int(i == j)
        instruments += 1
    tables[key] = means, ops, plane
checked_flags, checked_histories = 0, 0
for case in r["cases"]:
    m, law = case["M"], case["law"]
    means, ops, plane = tables[f"{m}:{law}"]
    a = b = [F(1), F(-1), F(0)]
    preparation = F(5, 144) if law == "native" else F(1, 9)
    assert preparation == F(case["preparation"]["accepted_probability"])
    reach, rejected = preparation, 1 - preparation
    positions = {"rotating": 2 * m, "reference": 2 * m}
    for row in case["records"]:
        assert F(row["reach"]) == reach
        if row["kind"] == "marker":
            arm, old = row["arm"], positions[row["arm"]]
            x, y = row["address"]
            success = F(0)
            for history in row["histories"]:
                position = old
                for coin, current in zip(history["coins"], history["currents"]):
                    before = int(position == x)
                    after = (0 if coin == 0 else 1) if law == "native" else 1 - before
                    assert before - after == current
                    if current:
                        position = y if current == 1 else x
                assert position == history["decoded_position"]
                probability = F(1, 2 ** len(history["coins"]))
                assert probability == F(history["probability"])
                if history["outcome"] == "hop":
                    success += probability
                    positions[arm] = position
                checked_histories += 1
            assert F(row["timeout"]) == reach * (1 - success)
            rejected += reach * (1 - success)
            reach *= success
            continue
        incoming = [list(map(F, v)) for v in row["audit_input_vectors"]]
        assert incoming == [a, b]
        common = math.lcm(*(x.denominator for x in a + b))
        aa, bb = [int(x * common) for x in a], [int(x * common) for x in b]
        norm = dot(aa, aa) + dot(bb, bb)
        if row["kind"] == "target":
            assert row["rotating_slot"] == positions["rotating"]
            assert row["reference_slot"] == positions["reference"]
            left, right = ops[positions["rotating"]], ops[positions["reference"]]
            accept = 0
        else:
            left = right = plane
            accept = 1
        probabilities = []
        for (g, d), (h, dd) in zip(left, right):
            av, bv = mv(g, aa), mv(h, bb)
            probabilities.append(F(dot(av, av) + dot(bv, bv), d * norm))
        assert probabilities == list(map(F, row["flag_probabilities"]))
        assert sum(probabilities) == 1
        checked_flags += len(probabilities)
        failure = reach * (1 - probabilities[accept])
        assert failure == F(row["first_failure"])
        rejected += failure
        reach *= probabilities[accept]
        if row["kind"] == "target":
            a, b = mv(means[positions["rotating"]], a), mv(means[positions["reference"]], b)
    assert reach == F(case["success"]) and rejected + reach == 1
    assert 2 * dot(a, b) / (dot(a, a) + dot(b, b)) == F(case["exact_visibility"])
print(json.dumps({"passed": True, "producer_imports": False, "exact_cases": len(r["cases"]),
                  "complete_instruments": instruments, "individually_checked_outcome_flags": checked_flags,
                  "individually_replayed_marker_histories": checked_histories,
                  "invalid_extension": "all3x3valid instruments plus orthogonal invalid projector give I4"}, indent=2))
