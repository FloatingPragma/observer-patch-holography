"""Independent analytic three-mode audit; imports no producer source."""
import json
import math
from fractions import Fraction as F
from pathlib import Path


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def mv(a, v):
    return [dot(row, v) for row in a]


def mean(slot, m):
    v, k = divmod(slot, m)
    result = [[F(i == j) for j in range(3)] for i in range(3)]
    for center, weight in ((v, F(m - k, m)), ((v + 1) % 3, F(k, m))):
        for leaf in range(3):
            if leaf != center:
                c = weight / 4
                result[center][center] -= c
                result[leaf][leaf] -= c
                result[center][leaf] += c
                result[leaf][center] += c
    return result


base = Path(__file__).resolve().parent
r = json.loads((base / "ERROR_RECEIPT.json").read_text())
maximum, probabilities = 0.0, 0
for case in r["cases"]:
    m, turns, delta = case["M"], case["turns"], case["delta_radians"]
    frames = 3 * m * turns
    t = 4 * math.sin(delta / 2) ** 2
    optical = (.25 + t / 2) / 12
    weights = [optical ** frames / 54, t * optical ** frames / 54,
               5 * t / 108 * (t / 12) ** frames]
    left = [[F(1), F(-1), F(0)], [F(1), F(1), F(-2)], [F(1)] * 3]
    right = [v[:] for v in left]
    for step in range(1, frames + 1):
        k = step % m
        rot, ref = (2 * m + step) % (3 * m), (2 * m + min(k, m - k)) % (3 * m)
        p, q = mean(rot, m), mean(ref, m)
        for _ in range(4):
            left, right = [mv(p, v) for v in left], [mv(q, v) for v in right]
    left = [[float(x) / math.sqrt(n) for x in v] for v, n in zip(left, (2, 6, 3))]
    right = [[float(x) / math.sqrt(n) for x in v] for v, n in zip(right, (2, 6, 3))]
    na = sum(w * dot(v, v) for w, v in zip(weights, left))
    nb = sum(w * dot(v, v) for w, v in zip(weights, right))
    cross = sum(w * dot(a, b) for w, a, b in zip(weights, left, right))
    visibility = 2 * cross / (na + nb)
    fidelity = sum(w * (v[0] - v[1]) ** 2 / 2 for w, v in zip(weights, left)) / na
    success = (na + nb) / (2 * 3 ** (4 * frames)) * (255 / 256) ** (2 * frames)
    errors = [abs(visibility - case["visibility"]), abs(fidelity - case["rotated_contrast_fidelity"]),
              abs(success / case["success"] - 1)]
    maximum = max(maximum, *errors)
    assert max(errors) < 2e-10
    for row in case["records"]:
        if "flag_probabilities" in row:
            ps = row["flag_probabilities"]
            assert min(ps) >= 0
            assert abs(sum(ps) / row["reach"] - 1) < 2e-11
            probabilities += len(ps)
assert not r["gate_failures"]
print(json.dumps({"passed": True, "producer_imports": False,
                  "independent_analytic_error_cases": len(r["cases"]),
                  "recorded_nonnegative_outcomes_and_stage_sums_checked": probabilities,
                  "maximum_spectral_prediction_disagreement": maximum,
                  "scope": "independent entire success path and every recorded probability sum; not a second flag-output-state replay"}, indent=2))
