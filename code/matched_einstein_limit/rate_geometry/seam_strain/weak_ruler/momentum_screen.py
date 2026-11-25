"""Normal and tangential Codazzi coefficients of the evolving five-station radar star.

Identically prepared copies of the vacancy pump supply both measurements.
Read the energy current before discarding the spent mobile. For radar,
retain its two heads and bridge strain and send a fresh probe; pool both
probe charges and orientations, keeping failures at their deadline. The four
tetrahedra use the same ten pooled weak round-trip lengths at every pump age. Barycentric worldlines are comoving, the lapse is one native tick
and the shift vanishes. Write B=dot(g)/2, R_0i=D_j B^j_i-D_i tr B. On an
open shared face orient n from left to right. The face-delta coefficient is
[B(n,i)]-n_i[tr B]. Its two tangential components can differ from zero.
Its normal component is minus [tr_face B]. Both sides inherit the same three
time-dependent lengths and hence the same induced Gram matrix and its time
derivative: the normal coefficient vanishes exactly. The bulk mixed
curvature vanishes in either homogeneous tetrahedron. Thus the odd null
contraction difference for the unit frame k±=u±n, namely 4 R(u,n), also
has zero face coefficient. This local statement excludes hinge edges and
outer boundaries. Computing the full even null response requires time-time
and space-space curvature as well.

For an independent material comparison take the face separating the native
head ports 0 and 39. Use graph-hop distance to its five named stations to
enroll two sharp side regions. All shared-zone ports enter a symmetrized
charge with half weight. Its signed current is the average of the sharp-cut
continuity currents. A convex lottery of an ordinary rail switch and an
independent re-pairing is calibrated solely on the number meter to match
the paired reference. Its energy gap is invariant under H -> H + cI. The graph partition does not supply a proper normal,
lapse, area, or a pullback of that cut to a shrinking normal prism.
"""

import argparse
from collections import deque
from functools import lru_cache
from itertools import combinations
import json

import numpy as np
import sympy as sp
from mpmath import iv

from . import current, flight, hinge
from .. import live_bounds, live_source
from ..source import certificate as cert

STATIONS = flight.STATIONS[:5]
# Poles are the actual stationary heads at the ends of the material route.
POLES = (0, 39)
FACE = (1, 38, 5)
SOURCES = ("reference", "ordinary", "decorrelated")


def distances(board, root):
    d = np.full(board.n, board.n, dtype=int)
    d[root] = 0
    queue = deque([root])
    while queue:
        p = queue.popleft()
        for j in board.neighbours[p]:
            if d[j] == board.n:
                d[j] = d[p] + 1
                queue.append(j)
    assert d.max() < board.n
    return d


def zones(board):
    """The shared zone contains the three face cells and all boundary ties."""
    assert set(STATIONS) == set((*POLES, *FACE))
    ds = [distances(board, p) for p in (*POLES, *FACE)]
    left = ds[0] < np.minimum.reduce(ds[1:])
    right = ds[1] < np.minimum.reduce([ds[0], *ds[2:]])
    zone = ~(left | right)
    assert left[POLES[0]] and right[POLES[1]] and all(zone[p] for p in FACE)
    return left, zone, right


@lru_cache(maxsize=1)
def source_states():
    (formed, rerolled), ferr, _ = live_source.formation_center()
    # Replacing the two independent marginals costs <=2*ferr in trace norm.
    return {
        "reference": (formed, (2, 1), ferr),
        "ordinary": (formed, (1, 2), ferr),
        "decorrelated": (rerolled, (2, 1), 2 * ferr),
    }


def source_currents(ports, sign):
    out = {}
    for tag, (effect, limit, pump_error) in current.current_effects(ports).items():
        readings = {}
        for source, (head, ratios, err) in source_states().items():
            full = sp.kronecker_product(sp.diag(*ratios) / 3, head)
            value = sign * sp.re(sp.trace(effect * full))
            readings[source] = (value, limit * err + (1 + err) * pump_error)
        out[tag] = readings
    return out


def material():
    src = live_source.source
    board = src.formation.federation(1)
    assert (src.formation.region(board)[0], src.formation.region(board)[-1]) == POLES
    left, shared, right = zones(board)
    ls = source_currents(np.flatnonzero(left).tolist(), 1)
    rs = source_currents(np.flatnonzero(right).tolist(), -1)
    # Outward from left and inward to right use one common orientation.
    mean = {
        q: {
            s: tuple((a + b) / 2 for a, b in zip(ls[q][s], rs[q][s], strict=True))
            for s in SOURCES
        }
        for q in ls
    }
    out = {}
    for k, meter in (
        ("left_boundary", ls),
        ("right_boundary", rs),
        ("half_shared_zone", mean),
    ):
        panel = {}
        for q, reads in meter.items():
            v, e = reads["reference"]
            panel[q] = {s: cert.enclosed(a, b) for s, (a, b) in reads.items()}
            panel[q].update(
                {
                    s + "_minus_reference": cert.enclosed(a - v, b + e)
                    for s, (a, b) in reads.items()
                    if s != "reference"
                }
            )
        # E -> E+cN. A common additive constant removing BOTH source
        # variations must be a zero of this two-row energy/number determinant.
        a, b = [x + "_minus_reference" for x in SOURCES[1:]]
        e, n = panel["energy"], panel["number"]
        ivv = iv.mpf
        det = ivv(e[a]) * ivv(n[b]) - ivv(e[b]) * ivv(n[a])
        panel["energy_number_contrast_determinant"] = live_bounds.endpoints(det)
        panel["zeroing_energy_origin"] = {
            s: live_bounds.endpoints(-ivv(e[s]) / ivv(n[s]))
            for s in (a, b)
            if ivv(n[s]).a * ivv(n[s]).b > 0
        }
        no, nd = ivv(n[a]), ivv(n[b])
        if no.b < 0 < nd.a:
            # Calibrate ONLY the ordinary admixture of the repaired state on
            # the number meter. True endpoint currents are affine in this mix.
            fraction = nd / (nd - no)
            energy = fraction * ivv(e[a]) + (1 - fraction) * ivv(e[b])
            assert 0 < fraction.a and fraction.b < 1
            panel["matched_number_lottery"] = dict(
                ordinary_probability=live_bounds.endpoints(fraction),
                decorrelated_probability=live_bounds.endpoints(1 - fraction),
                number_minus_reference=[0, 0],
                energy_minus_reference=live_bounds.endpoints(energy),
            )
        minimax = abs(det) / (abs(no) + abs(nd))
        panel["any_energy_origin_max_gap_lower"] = live_bounds.endpoints(minimax)[0]
        out[k] = panel
    return dict(
        ports=[int(z.sum()) for z in (left, shared, right)],
        boundary_bonds=[
            sum(int(v[i] != v[j]) for i, j in board.seams) for v in (left, right)
        ],
        oriented_cuts=out,
    )


def derivative_grams(h, hi, lo):
    """Differentiate the pump isometry at its 82-bit rational center."""
    peak = int(abs(h).sum(axis=1).max())
    assert peak * max(int(abs(hi).max()), int(abs(lo).max())) < 2**62
    base = hi.astype(object) * (1 << 36) + lo
    action = (h @ hi).astype(object) * (1 << 36) + h @ lo
    re, im = np.hsplit(action, 2)
    derivative = np.hstack([im, -re])  # -i (40H) V; divide by 40 after read

    def columns(w, rail):
        return np.split(
            w.reshape(2, -1, 4, 2, 2, 4)[:, :, :, :, rail, :]
            .transpose(1, 3, 0, 2, 4)
            .reshape(-1, 64),
            2,
            axis=1,
        )

    result = []
    for rail in (0, 1):
        x, y = columns(base, rail)
        a, b = columns(derivative, rail)
        r, i = x.T @ a + y.T @ b, x.T @ b - y.T @ a
        result.append(((r + r.T).reshape(8, 4, 8, 4), (i - i.T).reshape(8, 4, 8, 4)))
    return result


def radar_states():
    hi, lo, col, gram = live_source.pump_center(True)
    h = flight.hardware(-1)[0]
    dg = derivative_grams(h, hi, lo)
    hnorm = sp.Rational(int(abs(h).sum(axis=1).max()), 40)
    states = {}
    for name, (head, weights, ferr) in source_states().items():
        st = live_source.rational_read(gram, head, weights)
        sd = live_source.rational_read(dg, head, weights) / 40
        error = ferr + (1 + ferr) * (6 * col + 9 * col**2)
        states["dot_" + name] = (sd, 2 * hnorm * error)
        if name == "reference":
            ref, re = st, error
            states[name] = (st, error)
        else:
            states[name] = (st - ref, error + re)
    return states


def face_momentum(distances, velocities):
    """Three face-edge lengths fix B on the common tangent plane.

    Components are coefficients of the unit-normal face delta paired with
    the two coordinate vectors toward the other shared-face stations.
    """
    edges = list(combinations(STATIONS, 2))
    lengths = {frozenset(e): v for e, v in zip(edges, distances, strict=True)}
    rates = {frozenset(e): v for e, v in zip(edges, velocities, strict=True)}
    origin, first, second = FACE

    def entry(i, j, values):
        return 0 if i == j else values[frozenset((i, j))]

    def g(i, j):
        return (
            entry(origin, i, lengths) ** 2
            + entry(origin, j, lengths) ** 2
            - entry(i, j, lengths) ** 2
        ) / 2

    def half_dot_g(i, j):
        def term(a, b):
            return entry(a, b, lengths) * entry(a, b, rates)

        return (term(origin, i) + term(origin, j) - term(i, j)) / 2

    x, y, z = g(first, first), g(first, second), g(second, second)
    det = x * z - y * y
    assert det > 0
    b11, b12, b22 = (
        half_dot_g(first, first),
        half_dot_g(first, second),
        half_dot_g(second, second),
    )
    face_trace = (z * b11 - 2 * y * b12 + x * b22) / det
    out = []
    for apex in POLES:
        c, d = g(first, apex), g(second, apex)
        p, q = (z * c - y * d) / det, (x * d - y * c) / det
        height2 = g(apex, apex) - c * p - d * q
        assert height2 > 0
        height = height2**0.5
        out.append(
            [
                (half_dot_g(first, apex) - b11 * p - b12 * q) / height,
                (half_dot_g(second, apex) - b12 * p - b22 * q) / height,
            ]
        )
    # Right inward normal is +n, left inward normal is -n.
    return face_trace, [a + b for a, b in zip(*out, strict=True)]


@lru_cache(maxsize=1)
def geometry():
    edges, den, err = hinge.exact_ranges()
    limits = hinge.edge_errors(edges, den, err)
    data = hinge.paired_inputs(edges, den, limits, radar_states())
    out, shear = {}, {}
    for name in SOURCES:
        d = (
            data["reference"]
            if name == "reference"
            else [a + b for a, b in zip(data["reference"], data[name], strict=True)]
        )
        tr, coeff = face_momentum(d, data["dot_" + name])
        shear[name] = coeff
        out[name] = dict(
            common_tangential_expansion=live_bounds.endpoints(tr),
            tangential_face_coefficients=[live_bounds.endpoints(x) for x in coeff],
            normal_face_coefficient=[0, 0],
        )
    for name in SOURCES[1:]:
        out[name]["tangential_minus_reference"] = [
            live_bounds.endpoints(a - b)
            for a, b in zip(shear[name], shear["reference"], strict=True)
        ]
    return out


@lru_cache(maxsize=1)
def run():
    old = iv.dps
    try:
        iv.dps = 35
        return dict(
            stations=STATIONS,
            face=FACE,
            oriented_poles=POLES,
            pump_age=4,
            radar_divisor=hinge.DIVISOR,
            geometry=geometry(),
            material=material(),
        )
    finally:
        iv.dps = old


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    payload = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (hinge.HERE / "momentum_screen.json").write_text(payload)
    else:
        print(payload, end="")
