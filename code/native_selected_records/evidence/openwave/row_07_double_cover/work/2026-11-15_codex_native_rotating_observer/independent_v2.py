"""Independent v2 auditor: no producer imports, analytic error mixture."""
import json
import math
from fractions import Fraction as F
from pathlib import Path


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def mv(a, v):
    return [dot(row, v) for row in a]


def mm(a, b):
    return [[dot(row, col) for col in zip(*b)] for row in a]


def eye():
    return [[F(i == j) for j in range(3)] for i in range(3)]


def star(v):
    a = eye()
    for w in range(3):
        if w != v:
            a[v][v] -= F(1, 4)
            a[w][w] -= F(1, 4)
            a[v][w] += F(1, 4)
            a[w][v] += F(1, 4)
    return a


base = Path(__file__).resolve().parent
original = json.loads((base / "CONFIRM_RECEIPT.json").read_text())["cases"]
filtered = json.loads((base / "FILTERED_RECEIPT.json").read_text())["cases"]
checked = 0
for old, new in zip(original, filtered):
    ratio = F(1, (48 if new["law"] == "native" else 12) ** 3)
    assert new["visibility"] == old["visibility"]
    assert new["final_arm_vectors"] == old["final_arm_vectors"]
    assert F(new["success"]["exact"]) == ratio * F(old["success"]["exact"])
    total = [[F(0) for _ in range(3)] for _ in range(3)]
    for op in new["optical_filter_kraus"]:
        a = [list(map(F, row)) for row in op["matrix"]]
        b = mm(list(map(list, zip(*a))), a)
        for i in range(3):
            for j in range(3):
                total[i][j] += b[i][j] / op["squared_denominator"]
    assert total == eye()
    rejected = 1 - F(new["preparation"]["success_probability"])
    for row in new["optical_filters"]:
        assert sum(map(F, row["conditional_probabilities"])) == 1
        rejected += F(row["first_failure_probability"])
    for row in new["records"]:
        key = {"marker": "timeout_probability", "source": "failure_probability",
               "reference_coin": "first_failure_probability"}[row["kind"]]
        rejected += F(row[key])
    assert rejected + F(new["success"]["exact"]) == 1
    checked += 1
errors = json.loads((base / "FILTERED_ERROR.json").read_text())
assert not errors["gate_failures"]
max_disagreement = 0.0
for row in errors["cases"]:
    delta, turns = row["delta_radians"], row["turns"]
    t = 4 * math.sin(delta / 2) ** 2
    optical_filter_squared = (.25 + t / 2) / 12
    weights = [optical_filter_squared ** 3 / 54,
               t * optical_filter_squared ** 3 / 54,
               5 * t / 108 * (t / 12) ** 3]
    basis = [[1 / 2 ** .5, -1 / 2 ** .5, 0],
             [1 / 6 ** .5, 1 / 6 ** .5, -2 / 6 ** .5],
             [1 / 3 ** .5] * 3]
    a = b = eye()
    vertex = 2
    for hop in range(3 * turns):
        vertex = (vertex + 1) % 3
        for _ in range(4):
            a, b = mm(star(vertex), a), mm(star(2), b)
    left, right = [mv(a, v) for v in basis], [[x / 8 ** turns for x in mv(b, v)] for v in basis]
    na = sum(w * dot(v, v) for w, v in zip(weights, left))
    nb = sum(w * dot(v, v) for w, v in zip(weights, right))
    cross = sum(w * dot(v, u) for w, v, u in zip(weights, left, right))
    vis = 2 * cross / (na + nb)
    fidelity = sum(w * dot(v, basis[0]) ** 2 for w, v in zip(weights, left)) / na
    success = (na + nb) / (2 * 3 ** (12 * turns)) * (255 / 256) ** (3 * turns)
    differences = [abs(vis - row["visibility"]), abs(fidelity - row["rotated_contrast_fidelity"]),
                   abs(success / row["success"] - 1)]
    max_disagreement = max(max_disagreement, *differences)
    assert max(differences) < 2e-10
print(json.dumps({"passed": True, "exact_filtered_cases": checked,
                  "full_filter_completeness_and_all_rejections": True,
                  "independent_spectral_mixture_error_cases": len(errors["cases"]),
                  "maximum_disagreement": max_disagreement,
                  "producer_imports": False}, indent=2))
