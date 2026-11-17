"""Prospective fixed exploratory grid; every case retained."""
import json
from fractions import Fraction as F
from math import sqrt
from native import CALLS, dot, logprob, mv, source

cases = []
for m in (2, 4, 8, 16):
    for law in ("native", "swap"):
        means = [source(i, m, law)[0] for i in range(3 * m)]
        for h in (2, 4):
            for turns in (1, 2):
                for direction in (1, -1):
                    for reference in ("fixed", "retraced"):
                        a = b = (F(1), F(-1), F(0))
                        start = 2 * m
                        for step in range(1, 3 * m * turns + 1):
                            rot = (start + direction * step) % (3 * m)
                            k = step % m
                            offset = min(k, m - k) if reference == "retraced" else 0
                            ref = (start + direction * offset) % (3 * m)
                            for event in range(h):
                                a, b = mv(means[rot], a), mv(means[ref], b)
                        norm = dot(a, a) + dot(b, b)
                        visibility = 2 * dot(a, b) / norm
                        ab = dot(a, b)
                        overlap = (1 if ab >= 0 else -1) * sqrt(float(ab * ab / (dot(a, a) * dot(b, b))))
                        frames, events = 3 * m * turns, 3 * m * turns * h
                        prep = F(5, 144) if law == "native" else F(1, 9)
                        marker = F(255, 256) ** (frames * (2 if reference == "retraced" else 1)) if law == "native" else F(1)
                        filters = F(1, 48 if law == "native" else 12) ** frames
                        success = prep * marker * filters * norm / (4 * 3 ** events)
                        cases.append({"M": m, "h": h, "turns": turns, "law": law,
                                      "direction": direction, "reference": reference,
                                      "visibility": float(visibility), "exact_visibility": str(visibility),
                                      "ray_overlap": overlap, "log10_success": logprob(success)})
print(json.dumps({"candidate": "A-native-subdivided-observer-rotation-v3", "cases": cases,
                  "literal_endpoint_calls": CALLS, "kind": "exploratory full fixed grid"}, indent=2))
