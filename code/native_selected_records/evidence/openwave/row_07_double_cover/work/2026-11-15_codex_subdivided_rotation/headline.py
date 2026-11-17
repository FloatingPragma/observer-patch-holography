"""Standalone finite-frame check. Stdlib only; no file/network/process access.

Declared apparatus: coherent records, path erasers, triangular routing.
Native law supplies the source transfers and the observed minus/plus pattern.
Full eraser-flag ledgers are checked separately by independent_check.py.
"""
from fractions import Fraction as F
from math import atan2, degrees, log10, sqrt
import json

EDGES = ((0, 1), (1, 2), (0, 2))
CALLS = 0


def body(a, b, coin, law):
    global CALLS
    CALLS += 1
    if law == "swap":
        return b, a
    lo, hi = (a + b) // 2, (a + b + 1) // 2
    return (lo, hi) if coin == 0 else (hi, lo)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def mv(a, v):
    return [dot(row, v) for row in a]


def mean(addresses, law):
    fine = []
    for a, b in addresses:
        for coin in range(2):
            flags = [[[0] * 3 for _ in range(3)] for _ in range(3)]
            for x in range(3):
                s = [int(i == x) for i in range(3)]
                before = s[a]
                s[a], s[b] = body(s[a], s[b], coin, law)
                flags[before - s[a] + 1][s.index(1)][x] = 1
            fine += flags
    choices = 2 * len(addresses)
    for i in range(3):
        for j in range(3):
            assert sum(sum(a[k][i] * a[k][j] for k in range(3)) for a in fine) == choices * (i == j)
    return [[F(sum(a[i][j] for a in fine), choices) for j in range(3)] for i in range(3)]


def routing(slot, m, law):
    v, k = divmod(slot, m)
    addresses = []
    for word in range(2 * m):
        for shift in reversed(range(m.bit_length())):
            assert body(0, 0, (word >> shift) & 1, law) == (0, 0)
        center = v if word // 2 < m - k else (v + 1) % 3
        addresses.append([e for e in EDGES if center in e][word % 2])
    return addresses


def hop(old, new, law):
    lo, hi = sorted((old, new))
    hold, success = F(1), F(0)
    for _ in range(8):
        next_hold = F(0)
        for coin in range(2):
            before = int(old == lo)
            after, other = body(before, 1 - before, coin, law)
            current = before - after
            if current:
                assert (hi if current == 1 else lo) == new
                success += hold / 2
            else:
                next_hold += hold / 2
        hold = next_hold
    assert success + hold == 1
    return success


def xy(slot, m):
    points = ((1, 0), (-.5, sqrt(3) / 2), (-.5, -sqrt(3) / 2))
    v, k = divmod(slot, m)
    return [(1 - k / m) * points[v][i] + k / m * points[(v + 1) % 3][i] for i in range(2)]


def angle(a, b):
    return degrees(atan2(a[0] * b[1] - a[1] * b[0], dot(a, b)))


results = []
for m in (4, 8):
    for law in ("native", "swap"):
        means = [mean(routing(slot, m, law), law) for slot in range(3 * m)]
        all_mean = mean(EDGES, law)
        f = [F(1), F(-1), F(0)]
        filtered = [x - y for x, y in zip(f, mv(all_mean, f))]
        plane_probability = F(dot(filtered, filtered), 12 * dot(f, f))
        assert plane_probability == F(1, 48 if law == "native" else 12)
        for turns, reference in ((1, "retraced"), (2, "retraced"), (1, "stationary"), (2, "stationary")):
            a = b = f
            rot = ref = 2 * m
            marker = F(1)
            angles, reference_angles = [], []
            for step in range(1, 3 * m * turns + 1):
                k = step % m
                new_rot = (2 * m + step) % (3 * m)
                new_ref = (2 * m + min(k, m - k)) % (3 * m) if reference == "retraced" else 2 * m
                marker *= hop(rot, new_rot, law)
                if new_ref != ref:
                    marker *= hop(ref, new_ref, law)
                angles.append(angle(xy(rot, m), xy(new_rot, m)))
                reference_angles.append(angle(xy(ref, m), xy(new_ref, m)))
                rot, ref = new_rot, new_ref
                for _ in range(4):
                    a, b = mv(means[rot], a), mv(means[ref], b)
            ab, norm = dot(a, b), dot(a, a) + dot(b, b)
            visibility = 2 * ab / norm
            overlap = (1 if ab >= 0 else -1) * sqrt(float(ab * ab / (dot(a, a) * dot(b, b))))
            prep = F(5, 144) if law == "native" else F(1, 9)  # independently proved all-input instrument
            success = prep * marker * plane_probability ** (3 * m * turns) * norm / (4 * 3 ** (12 * m * turns))
            assert abs(sum(angles) - 360 * turns) < 1e-10 and abs(sum(reference_angles)) < 1e-10
            assert max(map(abs, angles)) < 90
            assert visibility < 0 if law == "native" and turns == 1 else visibility > 0
            if law == "swap" and reference == "retraced":
                assert visibility == 1
            results.append({"M": m, "law": law, "turns": turns, "reference": reference, "visibility": float(visibility),
                            "exact_visibility": str(visibility), "ray_overlap": overlap,
                            "max_frame_step_degrees": max(map(abs, angles)),
                            "frame_winding_degrees": sum(angles), "reference_winding_degrees": sum(reference_angles),
                            "log10_total_success": log10(success.numerator) - log10(success.denominator)})
print(json.dumps({"candidate": "A-native-subdivided-observer-rotation-v3", "results": results,
                  "literal_endpoint_evaluations": CALLS, "native_fine_isometries_checked": True,
                  "scope": "recorded planar observer rotation with local routing; no arbitrary-axis or exchange claim",
                  "official_row_awarded": False}, indent=2))
