"""Standalone inspectable core: native records -> preparation -> turn fringes.

Stdlib only; no file, network, process or existing-package access. Full outcome
ledgers and phase-error replay are separate, explicitly larger payloads.
"""
from fractions import Fraction as F
from math import log10, sqrt
import json

CALLS = 0
EDGES = ((0, 1), (1, 2), (0, 2))
I = [[F(i == j) for j in range(3)] for i in range(3)]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def mv(a, v):
    return [dot(row, v) for row in a]


def fine(edges, law):
    global CALLS
    result = []
    for a, b in edges:
        for coin in (0, 1):
            flags = [[[0] * 3 for _ in range(3)] for _ in range(3)]
            for x in range(3):
                s = [int(i == x) for i in range(3)]
                before = s[a]
                CALLS += 1
                if law == "native":
                    lo = (s[a] + s[b]) // 2
                    hi = s[a] + s[b] - lo
                    s[a], s[b] = (lo, hi) if coin == 0 else (hi, lo)
                else:
                    s[a], s[b] = s[b], s[a]
                flags[before - s[a] + 1][s.index(1)][x] = 1
            result.extend(flags)
    return result


def instrument(branches):
    count = len(branches) // 3
    prefix = [[0] * 3 for _ in range(3)]
    ops = []
    for k, a in enumerate(branches):
        if k:
            ops.append(([[prefix[i][j] - k * a[i][j] for j in range(3)] for i in range(3)],
                        count * k * (k + 1)))
        prefix = [[prefix[i][j] + a[i][j] for j in range(3)] for i in range(3)]
    ops.insert(0, (prefix, count * len(branches)))
    total = [[sum(F(sum(a[k][i] * a[k][j] for k in range(3)), d) for a, d in ops)
              for j in range(3)] for i in range(3)]
    assert total == I
    return ops


def mean(edges, law):
    ops = instrument(fine(edges, law))
    # Uniform record-erasure success is this mean divided by sqrt(3).
    return [[F(x, 2 * len(edges)) for x in row] for row in ops[0][0]]


results = []
for law in ("native", "swap"):
    stars = [mean(tuple(e for e in EDGES if v in e), law) for v in range(3)]
    full = mean(EDGES, law)
    edge = mean(((0, 1),), law)
    f = [F(1), F(-1), F(0)]
    # Native edge/vacuum path-minus preparation: (I-Pedge)/(2sqrt3).
    contrast = [[I[i][j] - edge[i][j] for j in range(3)] for i in range(3)]
    assert mv(contrast, [1, 1, 1]) == [0, 0, 0]
    assert all(contrast[0][j] == -contrast[1][j] and contrast[2][j] == 0 for j in range(3))
    # Full-triangle warmup + filter effect, checked on the whole input matrix.
    warm = fine(EDGES, law)
    effect = [[sum(F(dot(mv(contrast, [m[k][i] for k in range(3)]),
                            mv(contrast, [m[k][j] for k in range(3)])), 72)
                   for m in warm) for j in range(3)] for i in range(3)]
    if law == "native":
        assert effect == [[F(5, 144), -F(1, 72), 0], [-F(1, 72), F(5, 144), 0], [0, 0, F(1, 72)]]
    else:
        assert effect == [[F(i == j, 9) for j in range(3)] for i in range(3)]
    filter_vector = [f[i] - mv(full, f)[i] for i in range(3)]
    filter_probability = F(dot(filter_vector, filter_vector), 12 * dot(f, f))
    assert filter_probability == F(1, 48 if law == "native" else 12)
    for h in (4, 6):
        for turns in (1, 2):
            a = b = f
            vertex = 2
            for _ in range(3 * turns):
                vertex = (vertex + 1) % 3
                for event in range(h):
                    a, b = mv(stars[vertex], a), mv(stars[2], b)
            b = [x / 8 ** turns for x in b]
            visibility = 2 * dot(a, b) / (dot(a, a) + dot(b, b))
            overlap = float(dot(a, b)) / sqrt(float(dot(a, a) * dot(b, b)))
            # Basis0 input; each marker hop is current-confirmed, cap8.
            preparation = effect[0][0]
            marker = F(255, 256) ** (3 * turns) if law == "native" else F(1)
            success = preparation * filter_probability ** 3 * marker
            success *= F(dot(a, a) + dot(b, b), 4 * 3 ** (3 * turns * h))
            assert visibility < 0 if law == "native" and turns == 1 else visibility > 0
            results.append({"law": law, "events_per_frame": h, "turns": turns,
                            "visibility": float(visibility), "exact_visibility": str(visibility),
                            "ray_overlap": overlap,
                            "log10_total_success": log10(success.numerator) - log10(success.denominator)})
print(json.dumps({"candidate": "A-native-incidence-rotation-filtered-v2",
                  "native_endpoint_evaluations_in_this_check": CALLS,
                  "complete_finite_instruments": True, "results": results,
                  "supplied_architecture": ["coherent path/current memory and analyzers",
                                            "equilateral observer frame and local incidence routing"],
                  "scope": "one-axis native observer holonomy; no exchange or general SO3 claim",
                  "official_row_awarded": False}, indent=2))
