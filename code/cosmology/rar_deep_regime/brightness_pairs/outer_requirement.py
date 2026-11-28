"""Nuisance-bounded inward acceleration and a free brightness ordering on SPARC outer rings."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
import math

from mpmath import iv

from . import measure

HERE = measure.HERE
RECT = [(F("0.1"), F("0.1")), (F(1), F("0.1")), (F(1), F("1.4")), (F("0.1"), F("1.4"))]


def cut(poly, a, b, h):
    out = []
    for p, q in zip(poly, poly[1:] + poly[:1], strict=True):
        u, v = (a * z[0] + b * z[1] - h for z in (p, q))
        if u <= 0:
            out.append(p)
        if u * v < 0:
            out.append(
                tuple(x + u * (y - x) / (u - v) for x, y in zip(p, q, strict=True))
            )
    return out


def box(value):
    return iv.mpf(value.numerator) / value.denominator


def lower(v):
    return math.nextafter(float(v.a), -math.inf)


def upper(v):
    return math.nextafter(float(v.b), math.inf)


def slope_box(row):
    d, e, i, di = [F(row[a:b]) for a, b in ((15, 21), (22, 27), (30, 34), (35, 39))]

    def sini(z):
        return iv.sin(box(z) * iv.pi / 180) ** 2

    return (
        box(max(0, d - 2 * e) / d) * sini(max(0, i - 2 * di)) / sini(i),
        box((d + 2 * e) / d) * sini(min(90, i + 2 * di)) / sini(i),
    )


def run():
    iv.dps = 36
    files = [measure.DATA / n for n in ("table1.dat", "table2.dat")]
    meta = {
        t[:11].strip(): t
        for t in files[0].read_text().splitlines()
        if int(t[112:115]) in (1, 2) and F(t[30:34]) >= 30
    }
    all_rings = {name: [] for name in meta}
    for line, t in enumerate(files[1].read_text().splitlines(), 1):
        name = t[:11].strip()
        if name not in meta:
            continue
        r = F(t[19:25])
        if r < 3 * F(meta[name][71:76]):
            continue
        v, e = F(t[26:32]), F(t[33:38])
        gas, disk, bulge = [F(t[a:b]) for a, b in ((39, 45), (46, 52), (53, 59))]
        terms = [z * abs(z) for z in (gas, disk, bulge)]
        all_rings[name].append(
            (
                line,
                r,
                F(t[60:67]) + F(t[68:76]),
                max(0, v - 2 * e) ** 2,
                (v + 2 * e) ** 2,
                terms,
            )
        )
    results = {}
    for name, rings in all_rings.items():
        if not rings:
            continue
        hmin, hmax = slope_box(meta[name])
        polygon = RECT.copy()
        inner = RECT.copy()
        clean = []
        outputs = []
        skipped = []
        cuts = 0
        for line, r, brightness, slow, fast, terms in rings:
            g, d, b = terms
            values = [g + d * x + b * y for x, y in RECT]
            if min(values) <= 0:
                skipped.append(line)
                continue
            minimum = lower(box(slow) / (hmax * box(max(values))))
            maximum = (
                upper(box(fast) / (hmin * box(min(values)))) if hmin.a > 0 else None
            )
            outputs.append(
                {
                    "line": line,
                    "radius_kpc": float(r),
                    "brightness": float(brightness),
                    "total_over_baryonic": [minimum, maximum],
                }
            )
            clean.append((line, brightness, slow, fast, terms))
            # A nonnegative excess can always choose the lowest allowable h.
            # Round h downward for the feasibility test to favor the candidate.
            if hmin.a > 0:
                floor = F(lower(hmin))
                polygon = cut(polygon, d, b, fast / floor - g)
                inner = cut(inner, d, b, fast / F(upper(hmin)) - g)
        assert bool(inner) == bool(polygon)
        inward = bool(polygon)
        disallowed_pair = None
        for left, right in combinations(clean, 2):
            bright, faint = sorted((left, right), key=lambda z: z[1], reverse=True)
            if faint[1] <= 0 or bright[1] < 2 * faint[1]:
                continue
            cuts += 1
            # An increasing excess towards the fainter ring requires
            # slow_bright * B_faint <= fast_faint * B_bright.
            a, d, b = (
                bright[2] * f - faint[3] * g
                for f, g in zip(faint[4], bright[4], strict=True)
            )
            inner = cut(inner, d, b, -a)
            if polygon:
                polygon = cut(polygon, d, b, -a)
                if not polygon:
                    disallowed_pair = [bright[0], faint[0]]
        assert bool(inner) == bool(polygon)
        results[name] = {
            "rings": outputs,
            "nonpositive_template_possible_lines": skipped,
            "inward_allowed": inward,
            "brightness_constraints": cuts,
            "inward_and_ordering_allowed": bool(polygon),
            "last_constraint_in_empty_intersection": disallowed_pair,
        }
    rows = [(n, row) for n, d in results.items() for row in d["rings"]]
    strongest = max(rows, key=lambda v: v[1]["total_over_baryonic"][0])
    summaries = {
        "selected_galaxies": len(results),
        "enrolled_outer_rings": sum(map(len, all_rings.values())),
        "positive_template_rings": len(rows),
        "points_requiring_extra": sum(r["total_over_baryonic"][0] > 1 for _, r in rows),
        "galaxies_requiring_extra": sum(
            any(row["total_over_baryonic"][0] > 1 for row in d["rings"])
            for d in results.values()
        ),
        "galaxies_with_ordering_tests": sum(
            d["brightness_constraints"] > 0 for d in results.values()
        ),
        "inward_incompatible": [
            n for n, d in results.items() if not d["inward_allowed"]
        ],
        "ordering_incompatible_with_inward_allowed": [
            n
            for n, d in results.items()
            if d["inward_allowed"] and not d["inward_and_ordering_allowed"]
        ],
        "strongest_required_ratio": {"galaxy": strongest[0], **strongest[1]},
    }
    return {
        "summary": summaries,
        "galaxies": results,
        "inputs": {
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in [*files, HERE / "OUTER_REQUIREMENT.md"]
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    d = run()
    path = HERE / "outer_requirement.json"
    if cli.parse_args().write:
        path.write_text(json.dumps(d, indent=1, allow_nan=False) + "\n")
    else:
        assert json.loads(path.read_bytes()) == d
    print(d["summary"])
