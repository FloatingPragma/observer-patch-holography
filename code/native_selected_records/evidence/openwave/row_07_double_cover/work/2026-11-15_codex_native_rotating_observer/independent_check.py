"""No producer imports. Independently checks rational physics and full catalog."""
import json
from fractions import Fraction as Q
from pathlib import Path


def mul(a, b):
    return [[sum(x * y for x, y in zip(r, c)) for c in zip(*b)] for r in a]


def identity(n):
    return [[Q(i == j) for j in range(n)] for i in range(n)]


def mean(vertex, law):
    p = identity(3)
    for leaf in range(3):
        if leaf != vertex:
            c = Q(1, 4 if law == "native" else 2)
            p[vertex][vertex] -= c
            p[leaf][leaf] -= c
            p[vertex][leaf] += c
            p[leaf][vertex] += c
    return p


base = Path(__file__).resolve().parent
r = json.loads((base / "CONFIRM_RECEIPT.json").read_text())
for case in r["cases"]:
    v, h, turns, law = (case[k] for k in ("start", "h", "turns", "law"))
    leaves = [i for i in range(3) if i != v]
    a = [[Q(int(i == leaves[0]) - int(i == leaves[1]))] for i in range(3)]
    b = [row[:] for row in a]
    fixed = mean(v, law)
    for hop in range(3 * turns):
        v = (v + case["direction"]) % 3
        for _ in range(h):
            a, b = mul(mean(v, law), a), mul(fixed, b)
    a, b = [row[0] for row in a], [row[0] / 8 ** turns for row in b]
    dot = lambda x, y: sum(u * w for u, w in zip(x, y))
    visibility = 2 * dot(a, b) / (dot(a, a) + dot(b, b))
    assert visibility == Q(case["visibility"]["exact"])
    assert a == list(map(Q, case["final_arm_vectors"][0]))
    assert b == list(map(Q, case["final_arm_vectors"][1]))
    prep = Q(case["preparation"]["success_probability"])
    marker = (1 - Q(1, 2 ** case["cap"])) ** (3 * turns) if law == "native" else Q(1)
    success = prep * marker * (dot(a, a) + dot(b, b)) / (4 * 3 ** (3 * turns * h))
    assert success == Q(case["success"]["exact"])
    rejected = 1 - prep
    for row in case["records"]:
        if row["kind"] == "marker":
            assert sum(Q(t["probability"]) for t in row["histories"]) == 1
            for history in row["histories"]:
                assert history["outcome"] == "timeout" or history["currents"][-1] != 0
            rejected += Q(row["timeout_probability"])
        elif row["kind"] == "source":
            assert sum(map(Q, row["conditional_flags"])) == 1
            rejected += Q(row["failure_probability"])
        else:
            rejected += Q(row["first_failure_probability"])
    assert rejected + success == 1
catalog = json.loads((base / "INSTRUMENTS.json").read_text())["catalog"]
checked = 0
for law in catalog.values():
    for frame in law.values():
        for kind, ops in frame.items():
            total = [[Q(0) for _ in range(4)] for _ in range(4)]
            for op in ops:
                a = [list(map(Q, row)) for row in op["matrix"]]
                product = mul(list(map(list, zip(*a))), a)
                for i in range(4):
                    for j in range(4):
                        total[i][j] += product[i][j] / op["squared_denominator"]
            assert total == identity(4)
            checked += 1
assert not r["gate_failures"]
print(json.dumps({"independence": "standalone rational endpoint-average derivation; imports no producer code",
                  "exact_case_products_and_probabilities": len(r["cases"]),
                  "complete_four_dimensional_instruments": checked,
                  "all_first_rejection_sums_exact": True, "passed": True}, indent=2))
