"""Independent recount of the seam-count area-law receipt; the limit is analytic.

The producer and the evidence-package builder are never imported.  Every
level is rebuilt here from the definitions alone:

* the golden orbit by the Fibonacci residue trick of the source paper,
  ``floor(b phi) = floor(b p / q)`` for ``0 <= b < q`` with ``p`` the next
  Fibonacci number, each floor then validated exactly in ``Q(phi)``;
* exact signs in ``Q(sqrt 5)`` by the same-sign shortcut and the comparison
  of ``c^2`` with ``5 b^2`` for mixed signs;
* the read relation from per-axis near pairs (each axis term of the squared
  distance is at most ``1/q``) combined over the three axes and decided
  exactly on the sum, with the site index ``x1 q^2 + x2 q + x3``, streamed
  one ``x1`` chunk at a time so that no level needs the whole relation in
  memory;
* the crossing counts of the five planes from strictly opposite exact sides;
* the predictions from their formulas, the box-corrected value by the same
  fixed Gauss-Legendre rule and, as a cross-check, by a finer rule.

Counts, digests and exact expressions are compared byte for byte after
canonicalization; floating fields must be printed at twelve significant
digits and must agree with the recomputation within ``1e-9`` relative
tolerance.  The verifier checks the pins of the parent files as well.
"""
from fractions import Fraction
from hashlib import sha256
import json
from math import asin, pi, sqrt
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "crossing_read_area_law_receipt.json"
EVIDENCE_RECEIPT = "evidence/source_net_causal_poset/source_net_causal_limit_receipt.json"
PIN_PATHS = (
    "code/causal_refinement/crossing_read_area_law.py",
    "code/causal_refinement/verify_crossing_read_area_law.py",
    "code/causal_refinement/test_crossing_read_area_law.py",
    "Lean/Geometry/CrossingReadAreaLaw.lean",
    "evidence/source_net_causal_poset/build_causal_poset.py",
)
LEVELS = (5, 8, 13, 21, 34, 55)
FIBONACCI_INDEX = {5: 5, 8: 6, 13: 7, 21: 8, 34: 9, 55: 10}
PLANES = ((1, Fraction(1, 2)), (1, Fraction(1, 4)), (1, Fraction(3, 4)),
          (2, Fraction(1, 2)), (3, Fraction(1, 2)))
GAUSS_NODES = 64
CROSS_CHECK_NODES = 96
FLOAT_TOLERANCE = 1e-9


def require(condition, message):
    if not condition:
        raise ValueError(message)


# --------------------------------------------------------------------------
# Strict JSON
# --------------------------------------------------------------------------


def pairs(items):
    result = {}
    for key, value in items:
        require(key not in result, "duplicate JSON key")
        result[key] = value
    return result


def forbidden(token):
    raise ValueError("non-finite JSON token")


def load(path=OUTPUT):
    data = Path(path).read_bytes()
    require(len(data) <= 200_000, "receipt size")
    return json.loads(data.decode("ascii"), object_pairs_hook=pairs, parse_constant=forbidden)


def canonical(material):
    return (json.dumps(material, sort_keys=True, separators=(",", ":"),
                       ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")


def twelve_digits(value):
    require(type(value) is float, "float field required")
    require(value == float(f"{value:.12g}"), "float not printed at twelve significant digits")
    return value


def close(actual, expected):
    twelve_digits(actual)
    require(abs(actual - expected) <= FLOAT_TOLERANCE * max(1.0, abs(expected)),
            f"float {actual} differs from the recomputed {expected}")
    return actual


# --------------------------------------------------------------------------
# Exact arithmetic in Q(sqrt 5)
# --------------------------------------------------------------------------


def sqrt5_sign(c, b):
    """Exact sign of ``c + b sqrt5`` for int64 arrays.

    Same signs (zeros included) decide directly; for mixed signs ``c > 0 > b``
    the sign is that of ``c^2 - 5 b^2`` and for ``c < 0 < b`` its negative.
    """
    c = np.asarray(c, dtype=np.int64)
    b = np.asarray(b, dtype=np.int64)
    sc = np.sign(c)
    sb = np.sign(b)
    cmp = np.sign(c * c - 5 * b * b)
    same = np.where(sc != 0, sc, sb)
    mixed = np.where(sc > 0, cmp, -cmp)
    return np.where(sc * sb >= 0, same, mixed).astype(np.int8)


def phi_sign(a, b):
    """Exact sign of ``a + b phi`` (int64 arrays): ``2a + b + b sqrt5`` over two."""
    a = np.asarray(a, dtype=np.int64)
    b = np.asarray(b, dtype=np.int64)
    return sqrt5_sign(2 * a + b, b)


def fibonacci_successor(q):
    a, b = 1, 1
    while b < q:
        a, b = b, a + b
    require(b == q, "q is not a Fibonacci number")
    return a + b


def orbit(q):
    """``(m_b, b)`` with ``xi_b = m_b + b phi``, ``m_b = -floor(b p / q)``, validated exactly."""
    p = fibonacci_successor(q)
    m = np.array([-((b * p) // q) for b in range(q)], dtype=np.int64)
    b = np.arange(q, dtype=np.int64)
    require(bool(np.all(phi_sign(m, b) >= 0)), "orbit value below zero")
    require(bool(np.all(phi_sign(m - 1, b) < 0)), "orbit value not below one")
    require(len({(int(x), int(y)) for x, y in zip(m, b)}) == q, "orbit labels not distinct")
    return m, b


def axis_terms(m, b):
    """``(xi_x - xi_y)^2 = A + B phi`` for every label pair, in units of ``L^2``."""
    dm = m[:, None] - m[None, :]
    db = b[:, None] - b[None, :]
    return dm * dm + db * db, 2 * dm * db + db * db


# --------------------------------------------------------------------------
# Read relation and crossing counts
# --------------------------------------------------------------------------


def read_relation(q, m, b):
    """Ordered neighbour pairs ``(s, t)`` including ``t = s``, one chunk per ``x1``.

    A pair is a neighbour pair iff ``q (A + B phi) <= 1`` for the summed axis
    terms.  Each axis term is nonnegative, so a neighbour pair has every axis
    term at most ``1/q``; the near pairs of one axis are enumerated exactly and
    the three-axis sums are decided exactly.  The chunk of ``x1`` holds the
    complete neighbour lists of the ``q^2`` sites ``x1 q^2 + x2 q + x3``,
    sorted by ``(s, t)``; the generator yields ``(s, t)`` as int64 arrays.
    """
    A, B = axis_terms(m, b)
    near = phi_sign(q * A - 1, q * B) <= 0
    x, y = np.nonzero(near)
    xa, ya = x.astype(np.int64), y.astype(np.int64)
    Aa, Ba = A[x, y], B[x, y]
    for x1 in range(q):
        sel = xa == x1
        y1, A1, B1 = ya[sel], Aa[sel], Ba[sel]
        Asum = A1[:, None, None] + Aa[None, :, None] + Aa[None, None, :]
        Bsum = B1[:, None, None] + Ba[None, :, None] + Ba[None, None, :]
        ok = phi_sign(q * Asum - 1, q * Bsum) <= 0
        i1, i2, i3 = np.nonzero(ok)
        s = x1 * q * q + xa[i2] * q + xa[i3]
        t = y1[i1] * q * q + ya[i2] * q + ya[i3]
        order = np.lexsort((t, s))
        yield s[order], t[order]


def sides(q, m, b, axis, position):
    """Exact side of every site (index ``x1 q^2 + x2 q + x3``) relative to ``x_axis = position L``."""
    num, den = position.numerator, position.denominator
    label_side = phi_sign(den * m - num, den * b)
    idx = np.arange(q ** 3)
    label = (idx // (q * q), (idx // q) % q, idx % q)[axis - 1]
    return label_side[label]


def crossing_count(s, t, side):
    """Unordered neighbour pairs ``{s, t}``, ``s < t``, with strictly opposite sides, in one chunk."""
    keep = s < t
    return int(np.count_nonzero(side[s[keep]] * side[t[keep]] == -1))


def stream_relation(q, m, b, side_tables):
    """Digest, degrees, edge count and crossing counts, accumulated chunk by chunk.

    The digest is the SHA-256 of ``[[...row 0...],[...row 1...],...]\\n`` with
    rows ascending by site and entries ascending inside a row.
    """
    n = q ** 3
    h = sha256()
    h.update(b"[")
    degrees = np.zeros(n, dtype=np.int64)
    self_reads = 0
    entries = 0
    counts = [0] * len(side_tables)
    for x1, (s, t) in enumerate(read_relation(q, m, b)):
        lo = x1 * q * q
        starts = np.searchsorted(s, np.arange(lo, lo + q * q + 1))
        require(starts[0] == 0 and starts[-1] == len(s), "chunk holds exactly its own sites")
        for i in range(q * q):
            if lo + i:
                h.update(b",")
            h.update(("[" + ",".join(map(str, t[starts[i]:starts[i + 1]].tolist())) + "]").encode("ascii"))
        degrees[lo:lo + q * q] = np.diff(starts)
        self_reads += int(np.count_nonzero(s == t))
        entries += len(s)
        for k, side in enumerate(side_tables):
            counts[k] += crossing_count(s, t, side)
    h.update(b"]\n")
    require(self_reads == n, "every site reads itself exactly once")
    require(bool(np.all(degrees >= 1)), "degree table")
    require((entries - n) % 2 == 0, "the relation is symmetric")
    return h.hexdigest(), degrees, (entries - n) // 2, counts


# --------------------------------------------------------------------------
# Predictions
# --------------------------------------------------------------------------


def infinite_plane_prediction(q):
    return pi * q ** 4 / 4


def cap_truncated(q, position):
    margin = min(position, 1 - position)
    return Fraction(1, q) > margin * margin


def closed_form(q):
    return q ** 4 * (pi / 4 - 8 / (15 * sqrt(q)) + 1 / (12 * q))


def box_corrected_prediction(q, position, nodes):
    """``q^6 int_0^a len_P(d) G(sqrt(a^2 - d^2)) dd`` with ``L = 1``, ``a = 1/sqrt q``."""
    a = 1.0 / sqrt(q)
    P = float(position)
    breaks = sorted({x for x in (P, 1.0 - P) if 0.0 < x < a})
    edges = [0.0] + breaks + [a]
    x, w = np.polynomial.legendre.leggauss(nodes)
    total = 0.0
    for lo, hi in zip(edges, edges[1:]):
        tlo, thi = asin(lo / a), asin(hi / a)
        theta = 0.5 * (thi - tlo) * x + 0.5 * (thi + tlo)
        d = a * np.sin(theta)
        rho = a * np.cos(theta)
        length = np.maximum(0.0, np.minimum(1.0 - d, P) - np.maximum(0.0, P - d))
        G = pi * rho ** 2 - (8.0 / 3.0) * rho ** 3 + rho ** 4 / 2.0
        total += 0.5 * (thi - tlo) * float(np.dot(w, length * G * rho))
    return float(q) ** 6 * total


def plane_label(axis, position):
    return f"x{axis} = {position.numerator}L/{position.denominator}" if position.numerator != 1 \
        else f"x{axis} = L/{position.denominator}"


def box_corrected_ball_neighbours(q):
    """``(4 pi/3) q^{3/2} - (3 pi/2) q + (8/5) q^{1/2} - 1/6``: the cube mean of ``n vol(B(s, a) cap Omega)``."""
    return 4 * pi / 3 * q ** 1.5 - 3 * pi / 2 * q + 8 / 5 * sqrt(q) - 1 / 6


# --------------------------------------------------------------------------
# Level verification
# --------------------------------------------------------------------------


def rebuild_level(q, evidence_family):
    """Independent rebuild of one level; returns the expected level row with recomputed floats."""
    m, b = orbit(q)
    n = q ** 3
    side_tables = [sides(q, m, b, axis, position) for axis, position in PLANES]
    digest, degrees, edges, counts = stream_relation(q, m, b, side_tables)
    prediction = float(f"{infinite_plane_prediction(q):.12g}")
    crossings = []
    for (axis, position), side, count in zip(PLANES, side_tables, counts):
        truncated = cap_truncated(q, position)
        box = box_corrected_prediction(q, position, GAUSS_NODES)
        finer = box_corrected_prediction(q, position, CROSS_CHECK_NODES)
        require(abs(box - finer) <= 1e-11 * box, "quadrature rules disagree")
        if not truncated:
            require(abs(box - closed_form(q)) <= 1e-11 * box, "closed form disagrees with the quadrature")
        box = float(f"{box:.12g}")
        crossings.append({
            "plane": plane_label(axis, position),
            "axis": axis,
            "position_over_L": str(position),
            "sites_on_plane": int(np.count_nonzero(side == 0)),
            "sites_below": int(np.count_nonzero(side < 0)),
            "sites_above": int(np.count_nonzero(side > 0)),
            "count": count,
            "ratio_to_infinite_plane": float(f"{count / prediction:.12g}"),
            "cap_truncated_by_far_face": truncated,
            "box_corrected_prediction": box,
            "box_corrected_closed_form": None if truncated else "q^4*(pi/4 - 8/(15*sqrt(q)) + 1/(12*q))",
            "box_corrected_closed_form_value": None if truncated else float(f"{closed_form(q):.12g}"),
            "ratio_to_box_corrected": float(f"{count / box:.12g}"),
        })
    mid = [r["count"] for r in crossings if r["position_over_L"] == "1/2"]
    return {
        "q": q,
        "fibonacci_index": FIBONACCI_INDEX[q],
        "site_count": n,
        "undirected_edges": edges,
        "self_reads": n,
        "minimum_degree_including_self": int(degrees.min()),
        "maximum_degree_including_self": int(degrees.max()),
        "neighbour_digest_sha256": digest,
        "equals_evidence_receipt_digest": digest == evidence_family["neighbors_including_wait_sha256"],
        "equals_evidence_receipt_edges": edges == evidence_family["undirected_spatial_edges"],
        "read_radius_over_L_exact": "1/sqrt(q)",
        "read_radius_over_L": float(f"{1.0 / sqrt(q):.12g}"),
        "density_times_L_cubed": q ** 3,
        "n2_a4_L2_exact": q ** 4,
        "mean_neighbours_excluding_self": float(f"{2 * edges / n:.12g}"),
        "continuum_ball_neighbours": float(f"{4 * pi / 3 * q ** 1.5:.12g}"),
        "continuum_ball_neighbours_exact": "(4*pi/3)*q^(3/2)",
        "box_corrected_ball_neighbours": float(f"{box_corrected_ball_neighbours(q):.12g}"),
        "box_corrected_ball_neighbours_exact": "(4*pi/3)*q^(3/2) - (3*pi/2)*q + (8/5)*q^(1/2) - 1/6",
        "infinite_plane_prediction_exact": "pi*q^4/4",
        "infinite_plane_prediction": prediction,
        "planck_length_over_L_exact": "1/(sqrt(pi)*q^2)",
        "planck_length_over_L": float(f"{1.0 / (sqrt(pi) * q ** 2):.12g}"),
        "crossings": crossings,
        "midplane_counts_agree_across_axes": len(set(mid)) == 1,
    }


FLOAT_FIELDS = ("read_radius_over_L", "mean_neighbours_excluding_self", "continuum_ball_neighbours",
                "box_corrected_ball_neighbours", "infinite_plane_prediction", "planck_length_over_L")
CROSSING_FLOAT_FIELDS = ("ratio_to_infinite_plane", "box_corrected_prediction",
                         "box_corrected_closed_form_value", "ratio_to_box_corrected")


def accept_floats(row, expected):
    """Check the receipt's floats against the recomputation and return the expected row
    carrying the receipt's accepted float values, for the byte-for-byte comparison."""
    require(type(row) is dict and row.keys() == expected.keys(), "level fields")
    merged = dict(expected)
    for key in FLOAT_FIELDS:
        merged[key] = close(row[key], expected[key])
    require(type(row["crossings"]) is list and len(row["crossings"]) == len(expected["crossings"]),
            "crossing rows")
    merged["crossings"] = []
    for actual, exp in zip(row["crossings"], expected["crossings"]):
        require(type(actual) is dict and actual.keys() == exp.keys(), "crossing fields")
        item = dict(exp)
        for key in CROSSING_FLOAT_FIELDS:
            if exp[key] is None:
                require(actual[key] is None, "closed form present for a truncated cap")
            else:
                item[key] = close(actual[key], exp[key])
        merged["crossings"].append(item)
    return merged


def verify_level(row, evidence_family):
    require(type(row) is dict and type(row.get("q")) is int and row["q"] in LEVELS, "level q")
    expected = rebuild_level(row["q"], evidence_family)
    merged = accept_floats(row, expected)
    require(canonical(row) == canonical(merged), f"level q={row['q']} differs from the independent rebuild")
    return expected


def evidence_families(receipt_bytes):
    receipt = json.loads(receipt_bytes.decode("utf-8"), object_pairs_hook=pairs)
    families = {}
    for level in receipt["levels"]:
        for fam in level["families"]:
            if fam["dimension"] == 3:
                families[level["q"]] = fam
    return families


def verify(receipt, levels=None):
    """Verify the receipt; ``levels`` restricts the expensive rebuild (default: every level)."""
    require(type(receipt) is dict, "receipt object")
    require(receipt.keys() == {"schema", "source_pins", "evidence_receipt", "evidence_receipt_sha256",
                               "family", "quadrature", "planes", "levels", "scope"}, "receipt fields")
    require(receipt["schema"] == "crossing-read-area-law-v1", "schema")
    require(receipt["source_pins"] == {p: sha256((ROOT / p).read_bytes()).hexdigest() for p in PIN_PATHS},
            "stale parent source pin")
    require(receipt["evidence_receipt"] == EVIDENCE_RECEIPT, "evidence receipt path")
    evidence_bytes = (ROOT / EVIDENCE_RECEIPT).read_bytes()
    require(receipt["evidence_receipt_sha256"] == sha256(evidence_bytes).hexdigest(), "evidence receipt pin")
    families = evidence_families(evidence_bytes)
    family = receipt["family"]
    require(type(family) is dict and family.keys() == {
        "L_exact", "L", "population", "read_law", "density", "pair_convention", "tie_rule",
        "infinite_plane_prediction", "box_corrected_prediction", "dictionary"}, "family fields")
    require(family["L_exact"] == "2/sqrt(phi+2)", "L expression")
    close(family["L"], 2.0 / sqrt((1.0 + sqrt(5.0)) / 2.0 + 2.0))
    for key in ("population", "read_law", "density", "pair_convention", "tie_rule",
                "infinite_plane_prediction", "box_corrected_prediction", "dictionary"):
        require(type(family[key]) is str and 0 < len(family[key]) < 400, "family text")
    require(receipt["quadrature"] == {
        "rule": "Gauss-Legendre",
        "nodes_per_piece": GAUSS_NODES,
        "substitution": "d = a sin(theta) on every piece between the breakpoints of len_P",
        "transverse_disc_integral": "G(rho) = pi L^2 rho^2 - (8 L/3) rho^3 + rho^4/2",
    }, "quadrature declaration")
    require(receipt["planes"] == [{"axis": axis, "position_over_L": str(position),
                                   "plane": plane_label(axis, position)} for axis, position in PLANES],
            "plane list")
    require(receipt["scope"] == {
        "pair_convention_unordered_per_layer": True,
        "ordered_reads_are_twice_the_unordered_count": True,
        "self_reads_excluded_and_never_crossing": True,
        "side_decisions_exact_in_Q_phi": True,
        "neighbour_relation_from_evidence_package_builder": True,
        "continuum_limit_is_analytic": True,
        "finite_counts_demonstrate_the_limit": False,
        "horizon_constructed": False,
        "temperature_or_hawking_flux_claimed": False,
        "planck_length_physically_identified": False,
        "nat_per_crossing_read_is_declared_dictionary": True,
        "population_and_read_law_supplied": True,
        "native_read_law_selected": False,
    }, "scope flags")
    rows = receipt["levels"]
    require(type(rows) is list and [r.get("q") for r in rows] == list(LEVELS[:len(rows)]) and rows,
            "level list")
    summary = []
    for row in rows:
        require(type(row) is dict and row.get("q") in families, "evidence family for the level")
        if levels is not None and row["q"] not in levels:
            require(row["equals_evidence_receipt_digest"] is True and row["equals_evidence_receipt_edges"] is True,
                    "unverified level must at least claim the evidence digest")
            for c in row["crossings"]:
                require(type(c["count"]) is int and c["count"] > 0 and c["sites_on_plane"] == 0, "crossing row")
            continue
        expected = verify_level(row, families[row["q"]])
        require(expected["equals_evidence_receipt_digest"] and expected["equals_evidence_receipt_edges"],
                "rebuilt relation differs from the evidence package receipt")
        for c in expected["crossings"]:
            require(c["sites_on_plane"] == 0, "a site lies on a plane")
        summary.append({"q": row["q"], "undirected_edges": expected["undirected_edges"],
                        "midplane_count": expected["crossings"][0]["count"],
                        "ratio_to_infinite_plane": expected["crossings"][0]["ratio_to_infinite_plane"],
                        "ratio_to_box_corrected": expected["crossings"][0]["ratio_to_box_corrected"]})
    return {"schema": "crossing-read-area-law-v1", "levels_rebuilt": summary,
            "levels_in_receipt": [r["q"] for r in rows],
            "continuum_limit_is_analytic": True, "planck_length_physically_identified": False}


if __name__ == "__main__":
    selected = tuple(int(x) for x in sys.argv[1:]) or None
    print(json.dumps(verify(load(), selected), sort_keys=True))
