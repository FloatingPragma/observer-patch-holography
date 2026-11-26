"""Exact interval audit of all tetrahedra in the deposited six-station record."""

import hashlib
import json
from fractions import Fraction
from itertools import combinations, permutations
from math import isqrt, prod
from pathlib import Path

from ..observer_banks import audit

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[3] / "full_geometry" / "curved.json"
SCALE = 10**9
EDGES = list(combinations(range(4), 2))


def interval(text):
    return tuple(map(Fraction, text.strip("[]").split(",")))


def multiply(*positive):
    return prod(i[0] for i in positive), prod(i[1] for i in positive)


def divide(a, b):
    assert b[0] > 0
    return a[0] / b[1], a[1] / b[0]


def grid(i):
    return int(i[0] * SCALE), -((-i[1] * SCALE) // 1)


def fourth_root_grid(x):
    scaled = x * SCALE**4
    return isqrt(isqrt(scaled.numerator // scaled.denominator))


def fitted(quad, ranges, pairs):
    # N_i^4 = product_{j != i} r_ij/r_ji under the original four-clock fit.
    gains4 = {
        i: multiply(*(divide(ranges[i, j], ranges[j, i]) for j in quad if i != j))
        for i in quad
    }
    out = {}
    for i, j in combinations(quad, 2):
        squared4 = divide(
            multiply(*([pairs[i, j]] * 4)), multiply(gains4[i], gains4[j])
        )
        out[i, j] = fourth_root_grid(squared4[0]), fourth_root_grid(squared4[1]) + 1
    return out


def fingerprint(quad, pairs):
    a, b, c, d = quad
    opposite = (
        multiply(pairs[a, b], pairs[c, d]),
        multiply(pairs[a, c], pairs[b, d]),
        multiply(pairs[a, d], pairs[b, c]),
    )
    # Each order statistic is monotone in every entry, even for tied intervals.
    ordered = list(
        zip(sorted(x[0] for x in opposite), sorted(x[1] for x in opposite), strict=True)
    )
    return [grid(divide(x, ordered[-1])) for x in ordered[:2]]


def gap(a, b):
    return max(0, a[0] - b[1], b[0] - a[1])


def metric_gap(qa, qb, a, b):
    def distance(v):
        return max(
            gap(a[qa[i], qa[j]], b[tuple(sorted((v[i], v[j])))]) for i, j in EDGES
        )

    return min(distance(v) for v in permutations(qb))


def endpoint_trace_bounds():
    """Compare each independently observed head density with the exact-source enclosure."""
    bank_path = HERE.parent / "observer_banks" / "records.json"
    bank = json.loads(bank_path.read_bytes())
    audit.marginal_check(bank)
    head, _ = audit.instrument_widths()
    assert head.b < audit.iv.mpf("0.0000000008")
    full = HERE.parents[2]
    certificate = json.loads(
        (full / "certificate.json").read_text(), parse_float=Fraction
    )
    record = json.loads((full / "results.json").read_text())
    # Match the exact, Hermitian, trace-normalized binary64 center of the root proof.
    real, imag = record["x_real"], record["x_imag"]
    norm = sum(Fraction(real[b][i][i]) for b in range(2) for i in range(4))
    bounds = (Fraction(".000852"), Fraction(".001156"))
    # Sqrt(8) < 2.829; pointwise head moment uncertainty gives trace error
    # at most 8*head by orthogonality of the 32 Pauli block coordinates.
    sqrt8_upper = Fraction("2.829")
    slack = sqrt8_upper * certificate["state_HS_error_upper"] + 8 * Fraction(
        ".0000000008"
    )
    for arm, trace_bound in zip(bank["counts"], bounds, strict=True):
        hs_squared = 0
        for b in range(2):
            observed = audit.head_block(arm["pauli"], b)
            for i in range(4):
                for j in range(4):
                    x = (Fraction(real[b][i][j]) + Fraction(real[b][j][i])) / (2 * norm)
                    y = (Fraction(imag[b][i][j]) - Fraction(imag[b][j][i])) / (2 * norm)
                    rr, ii = observed[i][j]
                    hs_squared += (Fraction(rr, 4 * audit.N) - x) ** 2
                    hs_squared += (Fraction(ii, 4 * audit.N) - y) ** 2
        assert slack < trace_bound and 8 * hs_squared < (trace_bound - slack) ** 2
    return bounds


def report():
    raw = SOURCE.read_bytes()
    data = json.loads(raw)
    stations = data["stations"]
    indices = tuple(range(len(stations)))
    assert len(stations) == 6 and len(set(stations)) == 6
    ranges = {
        (i, j): interval(data["directed_mean_range_bounds"][i][j])
        for i, j in permutations(indices, 2)
    }
    pairs = {
        (i, j): multiply(ranges[i, j], ranges[j, i])
        for i, j in combinations(indices, 2)
    }
    held = {
        edge: grid(
            interval(data["projected_squared_distance_bounds"][edge[0]][edge[1]])
        )
        for edge in pairs
    }
    quads = list(combinations(indices, 4))
    local = {q: fitted(q, ranges, pairs) for q in quads}
    keys = {q: fingerprint(q, pairs) for q in quads}
    gaps = {
        "published_five_plus_one": [],
        "separate_four_clock": [],
        "arbitrary_clock_gains": [],
    }
    for qa, qb in combinations(quads, 2):
        gaps["published_five_plus_one"].append((metric_gap(qa, qb, held, held), qa, qb))
        gaps["separate_four_clock"].append(
            (metric_gap(qa, qb, local[qa], local[qb]), qa, qb)
        )
        invariant = max(gap(x, y) for x, y in zip(keys[qa], keys[qb], strict=True))
        gaps["arbitrary_clock_gains"].append((invariant, qa, qb))

    def closest(rows):
        lower, a, b = min(rows)
        return {
            "lower_billionths": lower,
            "station_sets": [[stations[i] for i in q] for q in (a, b)],
        }

    endpoint = []
    for trace_bound in endpoint_trace_bounds():
        # Every true native response effect has spectrum in [0,3]. The
        # source difference has trace zero, giving |delta r| <= 1.5 ||D||_1.
        shift = 3 * trace_bound / 2
        shifted = {k: (v[0] - shift, v[1] + shift) for k, v in ranges.items()}
        products = {
            edge: multiply(shifted[edge], shifted[edge[::-1]]) for edge in pairs
        }
        f = {q: fingerprint(q, products) for q in quads}
        clearances = [
            (max(gap(x, y) for x, y in zip(f[a], f[b], strict=True)), a, b)
            for a, b in combinations(quads, 2)
        ]
        endpoint.append(
            {
                "source_trace_norm_upper_millionths": int(trace_bound * 10**6),
                "excluded_of_105": sum(row[0] > 0 for row in clearances),
                "reference_vs_any": closest(
                    [v for v in clearances if v[1] == quads[0]]
                ),
            }
        )
    return {
        "curve_record_sha256": hashlib.sha256(raw).hexdigest(),
        "reference": [stations[i] for i in quads[0]],
        "subsets": len(quads),
        "pairs_per_comparison": len(gaps["arbitrary_clock_gains"]),
        "permutations_per_metric_pair": 24,
        "transferred_to_two_counted_endpoints": endpoint,
        "lower_bounds": {
            mode: {
                "any_pair": closest(rows),
                "reference_vs_any": closest([v for v in rows if v[1] == quads[0]]),
            }
            for mode, rows in gaps.items()
        },
    }


if __name__ == "__main__":
    print(json.dumps(report(), indent=2))
