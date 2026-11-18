"""Exact seam-count area-law receipt on the golden source family.

For every Fibonacci level ``q`` the population ``S_q`` of ``q^3`` sites in the
cube ``Omega = [0, L]^3`` (``L = 2/sqrt(phi + 2)``) and the complete-neighbour
read relation inside the radius ``a_q = L/sqrt(q)`` are rebuilt exactly as the
evidence package ``evidence/source_net_causal_poset/build_causal_poset.py``
builds them; that module is imported (its ``__main__`` guard keeps the import
free of side effects) and its functions ``orbit``, ``axis_tables``,
``site_coordinates``, ``site_graph``, ``phi_sign_array`` and
``neighbour_digest`` are used unchanged.

Crossing count.  For a plane ``x_i = P`` the receipt counts the unordered
site pairs ``{s, t}`` with ``s != t``, ``|t - s| <= a_q`` and the two ``i``-th
coordinates on strictly opposite sides of the plane.  One unordered crossing
pair is one crossing read per layer from the lower layer to the upper layer;
the ordered read count (both directions) is twice the unordered count.  A
same-site read never crosses a plane.  Tie rule: a site whose coordinate lies
exactly on the plane would belong to neither side; every side decision is
exact in ``Q(phi)`` and the receipt records that no site lies on any of the
five planes ``x_1 = L/2``, ``x_1 = L/4``, ``x_1 = 3L/4``, ``x_2 = L/2``,
``x_3 = L/2``.

Predictions.  The infinite-plane continuum value per layer through the full
cross-section of area ``L^2`` is ``L^2 (pi/4) n^2 a^4 = pi q^4 / 4`` with
``n = q^3/L^3``.  The box-corrected continuum value is
``n^2 * int_{s in Omega, s_i < P} vol{t in Omega : t_i > P, |t - s| <= a} d^3 s``.
By Fubini in the difference vector ``d = t - s`` this equals
``n^2 * int_{|d| <= a, d_i > 0} len_P(d_i) prod_{j != i} (L - |d_j|) d^3 d`` with
``len_P(d_i) = |[0, L - d_i] cap (P - d_i, P)|``, and the two transverse
coordinates integrate over the disc of radius ``rho = sqrt(a^2 - d_i^2)`` to
``G(rho) = pi L^2 rho^2 - (8L/3) rho^3 + rho^4 / 2``.  The remaining
one-dimensional integral is evaluated with a fixed Gauss-Legendre rule after
the substitution ``d_i = a sin(theta)``, piecewise between the breakpoints of
``len_P``.  When the cap is not truncated by a far face
(``a <= min(P, L - P)``) the value has the closed form
``q^4 (pi/4 - 8/(15 sqrt q) + 1/(12 q))``.

Dictionary.  Reading one nat per crossing read gives the cut entropy per unit
area ``(pi/4) n^2 a^4``; identifying it with ``1/(4 l^2)`` fixes
``l = 1/(sqrt(pi) n a^2) = L/(sqrt(pi) q^2)``.  This is a declared calibration.

The receipt carries no timestamps and no absolute paths.  Derived floats are
printed at twelve significant digits; counts and exact expressions are exact.
"""
from fractions import Fraction
from hashlib import sha256
import json
from math import asin, pi, sqrt
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
RER = HERE.parents[1]
OUTPUT = HERE / "crossing_read_area_law_receipt.json"
EVIDENCE = RER / "evidence" / "source_net_causal_poset"
EVIDENCE_RECEIPT = "evidence/source_net_causal_poset/source_net_causal_limit_receipt.json"
PIN_PATHS = (
    "code/causal_refinement/crossing_read_area_law.py",
    "code/causal_refinement/verify_crossing_read_area_law.py",
    "code/causal_refinement/test_crossing_read_area_law.py",
    "Lean/Geometry/CrossingReadAreaLaw.lean",
    "evidence/source_net_causal_poset/build_causal_poset.py",
)
LEVELS = (5, 8, 13, 21, 34, 55)
# (axis, position over L): the mid-plane and quarter planes in x1, the mid-planes in x2 and x3.
PLANES = ((1, Fraction(1, 2)), (1, Fraction(1, 4)), (1, Fraction(3, 4)),
          (2, Fraction(1, 2)), (3, Fraction(1, 2)))
GAUSS_NODES = 64
ROW_CHUNK = 4096

if str(EVIDENCE) not in sys.path:
    sys.path.insert(0, str(EVIDENCE))
import build_causal_poset as bcp  # noqa: E402


def rounded(x):
    """Twelve significant digits, the receipt convention for derived floats."""
    return float(f"{float(x):.12g}")


def canonical(material):
    return (json.dumps(material, sort_keys=True, separators=(",", ":"),
                       ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")


def plane_label(axis, position):
    return f"x{axis} = {position.numerator}L/{position.denominator}" if position.numerator != 1 \
        else f"x{axis} = L/{position.denominator}"


# --------------------------------------------------------------------------
# Population and read relation (the evidence package's construction)
# --------------------------------------------------------------------------


def population(q):
    """Golden orbit as (m, b) pairs, the per-axis Q(phi) tables and the site table."""
    values = bcp.orbit(q)
    A, B = bcp.axis_tables(values)
    sites = bcp.site_coordinates(q, 3)
    return values, A, B, sites


def read_relation(q, values, A, B, sites):
    """CSR neighbour lists including the same-site read, exactly as the evidence package builds them."""
    return bcp.site_graph(q, 3, values, A, B, sites)


def scaled_read_relation(q, values, A, B, sites, radius_sq_multiplier=Fraction(1)):
    """The same construction with the decision ``q (A + B phi) <= multiplier``.

    Generalizes ``build_causal_poset.site_graph`` (same per-axis candidate
    lists, same exact sign rule, same row order) by an exact rational factor on
    the squared radius; at multiplier one it reproduces ``site_graph`` edge for
    edge.  Used by the mutation tests (halved radius).
    """
    num, den = radius_sq_multiplier.numerator, radius_sq_multiplier.denominator
    if num <= 0 or den <= 0:
        raise ValueError("positive radius multiplier required")
    within = bcp.phi_sign_array(den * q * A - num, den * q * B) <= 0
    cand = [np.flatnonzero(within[x]) for x in range(q)]
    n = q ** 3
    powers = [q * q, q, 1]
    rows_out = []
    for s in range(n):
        coords = sites[s]
        lists = [cand[int(coords[axis])] for axis in range(3)]
        grids = np.meshgrid(*lists, indexing="ij")
        Asum = np.zeros(grids[0].shape, dtype=np.int64)
        Bsum = np.zeros(grids[0].shape, dtype=np.int64)
        nb = np.zeros(grids[0].shape, dtype=np.int64)
        for axis in range(3):
            Asum += A[int(coords[axis])][grids[axis]]
            Bsum += B[int(coords[axis])][grids[axis]]
            nb += grids[axis] * powers[axis]
        ok = bcp.phi_sign_array(den * q * Asum - num, den * q * Bsum) <= 0
        rows_out.append(np.sort(nb[ok]))
    indptr = np.zeros(n + 1, dtype=np.int64)
    indptr[1:] = np.cumsum([len(r) for r in rows_out])
    indices = np.concatenate(rows_out).astype(np.int32)
    return indptr, indices


# --------------------------------------------------------------------------
# Sides of a plane and the crossing count
# --------------------------------------------------------------------------


def axis_sides(values, position):
    """Exact side (-1, 0, +1) of every orbit label relative to ``xi = position``."""
    num, den = position.numerator, position.denominator
    a = np.array([den * m - num for m, _b in values], dtype=np.int64)
    b = np.array([den * b for _m, b in values], dtype=np.int64)
    return bcp.phi_sign_array(a, b)


def site_sides(values, sites, axis, position):
    """Side of every site relative to the plane ``x_axis = position * L`` (axis is 1, 2 or 3)."""
    return axis_sides(values, position)[sites[:, axis - 1]]


def crossing_count(indptr, indices, side, include_self=False, ordered=False):
    """Number of neighbour pairs with strictly opposite sides.

    Default: unordered pairs ``{s, t}`` with ``t > s``.  ``include_self`` admits
    ``t = s`` (never a crossing); ``ordered`` counts both directions.
    """
    n = len(indptr) - 1
    side = np.asarray(side, dtype=np.int8)
    total = 0
    for lo in range(0, n, ROW_CHUNK):
        hi = min(n, lo + ROW_CHUNK)
        cols = indices[indptr[lo]:indptr[hi]].astype(np.int64)
        rws = np.repeat(np.arange(lo, hi, dtype=np.int64), np.diff(indptr[lo:hi + 1]))
        if ordered:
            keep = cols != rws if not include_self else np.ones(len(cols), dtype=bool)
        else:
            keep = cols >= rws if include_self else cols > rws
        total += int(np.count_nonzero(side[rws[keep]] * side[cols[keep]] == -1))
    return total


def pair_count(indptr, include_self=False):
    """Unordered neighbour pairs; with ``include_self`` the ``n`` same-site reads are added."""
    n = len(indptr) - 1
    edges = (int(indptr[-1]) - n) // 2
    return edges + n if include_self else edges


# --------------------------------------------------------------------------
# Continuum predictions
# --------------------------------------------------------------------------


def infinite_plane_prediction(q):
    """``L^2 (pi/4) n^2 a^4 = pi q^4 / 4`` unordered crossing pairs per layer."""
    return pi * q ** 4 / 4


def cap_truncated(q, position):
    """Whether the read ball of a site next to the plane reaches a far face."""
    a = Fraction(1, 1)  # compare a^2 = 1/q with min(P, 1 - P)^2 exactly
    margin = min(position, 1 - position)
    return Fraction(1, q) > margin * margin * a


def box_corrected_closed_form(q):
    """``q^4 (pi/4 - 8/(15 sqrt q) + 1/(12 q))``, valid when the cap is not truncated."""
    return q ** 4 * (pi / 4 - 8 / (15 * sqrt(q)) + 1 / (12 * q))


def box_corrected_prediction(q, position, nodes=GAUSS_NODES):
    """Box-corrected continuum count, units ``L = 1``: ``n = q^3``, ``a = 1/sqrt q``.

    ``n^2 int_0^a len_P(d) G(sqrt(a^2 - d^2)) dd`` with
    ``G(rho) = pi rho^2 - (8/3) rho^3 + rho^4/2`` and
    ``len_P(d) = max(0, min(1 - d, P) - max(0, P - d))``; substitution
    ``d = a sin(theta)``, Gauss-Legendre with ``nodes`` points on every piece
    between the breakpoints of ``len_P``.
    """
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


def planck_length_over_L(q):
    """The dictionary value ``l/L = 1/(sqrt(pi) q^2)``."""
    return 1.0 / (sqrt(pi) * q ** 2)


def box_corrected_ball_neighbours(q):
    """Mean over the cube of ``n vol(B(s, a) cap Omega)``, the ball moments with the box factors:
    ``(4 pi/3) q^{3/2} - (3 pi/2) q + (8/5) q^{1/2} - 1/6``."""
    return 4 * pi / 3 * q ** 1.5 - 3 * pi / 2 * q + 8 / 5 * sqrt(q) - 1 / 6


# --------------------------------------------------------------------------
# One level
# --------------------------------------------------------------------------


def level_report(q, evidence_family=None):
    if q not in bcp.FIBONACCI_INDEX or q < 5:
        raise ValueError("q must be a Fibonacci number of the receipt family, at least 5")
    values, A, B, sites = population(q)
    indptr, indices = read_relation(q, values, A, B, sites)
    n = q ** 3
    degrees = np.diff(indptr)
    edges = pair_count(indptr)
    digest = bcp.neighbour_digest(indptr, indices)
    prediction = rounded(infinite_plane_prediction(q))
    crossings = []
    for axis, position in PLANES:
        side = site_sides(values, sites, axis, position)
        count = crossing_count(indptr, indices, side)
        truncated = cap_truncated(q, position)
        box = rounded(box_corrected_prediction(q, position))
        row = {
            "plane": plane_label(axis, position),
            "axis": axis,
            "position_over_L": str(position),
            "sites_on_plane": int(np.count_nonzero(side == 0)),
            "sites_below": int(np.count_nonzero(side < 0)),
            "sites_above": int(np.count_nonzero(side > 0)),
            "count": count,
            "ratio_to_infinite_plane": rounded(count / prediction),
            "cap_truncated_by_far_face": truncated,
            "box_corrected_prediction": box,
            "box_corrected_closed_form": None if truncated else "q^4*(pi/4 - 8/(15*sqrt(q)) + 1/(12*q))",
            "box_corrected_closed_form_value": None if truncated else rounded(box_corrected_closed_form(q)),
            "ratio_to_box_corrected": rounded(count / box),
        }
        crossings.append(row)
    mid = [r["count"] for r in crossings if r["position_over_L"] == "1/2"]
    report = {
        "q": q,
        "fibonacci_index": bcp.FIBONACCI_INDEX[q],
        "site_count": n,
        "undirected_edges": edges,
        "self_reads": n,
        "minimum_degree_including_self": int(degrees.min()),
        "maximum_degree_including_self": int(degrees.max()),
        "neighbour_digest_sha256": digest,
        "equals_evidence_receipt_digest": None if evidence_family is None
        else digest == evidence_family["neighbors_including_wait_sha256"],
        "equals_evidence_receipt_edges": None if evidence_family is None
        else edges == evidence_family["undirected_spatial_edges"],
        "read_radius_over_L_exact": "1/sqrt(q)",
        "read_radius_over_L": rounded(1.0 / sqrt(q)),
        "density_times_L_cubed": q ** 3,
        "n2_a4_L2_exact": q ** 4,
        "mean_neighbours_excluding_self": rounded(2 * edges / n),
        "continuum_ball_neighbours": rounded(4 * pi / 3 * q ** 1.5),
        "continuum_ball_neighbours_exact": "(4*pi/3)*q^(3/2)",
        "box_corrected_ball_neighbours": rounded(box_corrected_ball_neighbours(q)),
        "box_corrected_ball_neighbours_exact": "(4*pi/3)*q^(3/2) - (3*pi/2)*q + (8/5)*q^(1/2) - 1/6",
        "infinite_plane_prediction_exact": "pi*q^4/4",
        "infinite_plane_prediction": prediction,
        "planck_length_over_L_exact": "1/(sqrt(pi)*q^2)",
        "planck_length_over_L": rounded(planck_length_over_L(q)),
        "crossings": crossings,
        "midplane_counts_agree_across_axes": len(set(mid)) == 1,
    }
    return report


def evidence_families():
    receipt = json.loads((RER / EVIDENCE_RECEIPT).read_bytes())
    families = {}
    for level in receipt["levels"]:
        for fam in level["families"]:
            if fam["dimension"] == 3:
                families[level["q"]] = fam
    return families


def report(levels=LEVELS, log=None):
    families = evidence_families()
    rows = []
    for q in levels:
        if log is not None:
            print(f"q={q}: rebuilding", file=log, flush=True)
        rows.append(level_report(q, families.get(q)))
        if log is not None:
            r = rows[-1]
            print(f"q={q}: edges {r['undirected_edges']} midplane {r['crossings'][0]['count']} "
                  f"ratio_inf {r['crossings'][0]['ratio_to_infinite_plane']} "
                  f"ratio_box {r['crossings'][0]['ratio_to_box_corrected']}", file=log, flush=True)
    return {
        "schema": "crossing-read-area-law-v1",
        "source_pins": {name: sha256((RER / name).read_bytes()).hexdigest() for name in PIN_PATHS},
        "evidence_receipt": EVIDENCE_RECEIPT,
        "evidence_receipt_sha256": sha256((RER / EVIDENCE_RECEIPT).read_bytes()).hexdigest(),
        "family": {
            "L_exact": "2/sqrt(phi+2)",
            "L": rounded(2.0 / sqrt((1.0 + sqrt(5.0)) / 2.0 + 2.0)),
            "population": "s(b) = L (xi_b1, xi_b2, xi_b3), xi_b = b phi - floor(b phi), 0 <= b_i < q, q^3 sites in [0, L]^3",
            "read_law": "an event at layer j+1 reads every event at layer j within distance a_q = L/sqrt(q), the same site included",
            "density": "n_q = q^3 / L^3 sites per unit volume",
            "pair_convention": "unordered site pairs {s, t}, s != t, |t - s| <= a_q, sides strictly opposite; one such pair is one crossing read per layer; the ordered read count is twice this",
            "tie_rule": "a site with its coordinate exactly on the plane belongs to neither side; side decisions are exact in Q(phi); sites_on_plane records that no site lies on any plane",
            "infinite_plane_prediction": "L^2 (pi/4) n^2 a^4 = pi q^4 / 4 unordered crossing pairs per layer through the full cross-section",
            "box_corrected_prediction": "n^2 int_{s in Omega, s_i < P} vol{t in Omega : t_i > P, |t - s| <= a} d^3 s, the cap truncated by the box",
            "dictionary": "one nat per crossing read; (pi/4) n^2 a^4 = 1/(4 l^2) fixes l = 1/(sqrt(pi) n a^2) = L/(sqrt(pi) q^2)",
        },
        "quadrature": {
            "rule": "Gauss-Legendre",
            "nodes_per_piece": GAUSS_NODES,
            "substitution": "d = a sin(theta) on every piece between the breakpoints of len_P",
            "transverse_disc_integral": "G(rho) = pi L^2 rho^2 - (8 L/3) rho^3 + rho^4/2",
        },
        "planes": [{"axis": axis, "position_over_L": str(position), "plane": plane_label(axis, position)}
                   for axis, position in PLANES],
        "levels": rows,
        "scope": {
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
        },
    }


if __name__ == "__main__":
    levels = tuple(int(x) for x in sys.argv[1:]) or LEVELS
    OUTPUT.write_bytes((json.dumps(report(levels, log=sys.stderr), indent=2, sort_keys=True,
                                   ensure_ascii=True, allow_nan=False) + "\n").encode("ascii"))
    print(OUTPUT)
