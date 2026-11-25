"""Clock length, geodesic and gain checks throughout the certified source enclosure."""

import argparse
from itertools import combinations
import json
import numpy as np
import sympy as sp
from mpmath import iv
from scipy.sparse.csgraph import shortest_path

from ..full_algebra import certificate as full
from ..isometric_collar import certify as PRIOR
from .measure import HERE, Ruler, graph, optics


def bounds(x):
    return [float(x.a), float(x.b)]


def inputs():
    # Only the full station is consumed. The separate even source does not enter.
    mid = [PRIOR.sym_matrix(p) for p in optics.source()["state"]]
    mid = [(p + p.H) / 2 for p in mid]
    mid = [p / sum(q.trace() for q in mid) for p in mid]
    assert all(
        all((p - sp.Rational(".02") * sp.eye(4))[:i, :i].det() > 0 for i in range(1, 5))
        for p in mid
    )
    total = sum(mid, sp.zeros(4))
    right = sp.Matrix(
        2, 2, lambda i, j: sum(total[2 * k + i, 2 * k + j] for k in range(2))
    )
    left = [
        sp.Matrix(2, 2, lambda i, j: sum(m[2 * i + k, 2 * j + k] for k in range(2)))
        for m in mid
    ]
    targets = {
        "ordinary": [total / 2 - m for m in mid],
        "unpair": [sp.kronecker_product(part, right) - m for part, m in zip(left, mid)],
    }
    paulis = [
        sp.Matrix([[0, 1], [1, 0]]),
        sp.Matrix([[0, -sp.I], [sp.I, 0]]),
        sp.diag(1, -1),
    ]
    gen = [sp.kronecker_product(sp.eye(2), p) for p in paulis]
    r = sp.Rational("2e-7")
    b = [full.ri_band(sum(full.inner(m, g) for m in mid), 3 * r) for g in gen]
    length = iv.sqrt(sum(a**2 for a in b))
    rotations = [
        [(-sp.I * (g * m - m * g) / 2).applyfunc(sp.expand) for m in mid] for g in gen
    ]
    return mid, targets, [a / length for a in b], rotations, r


def variable_coefficients(e, source):
    # Both unpairing and the right-marginal rotation fix both parity-block left
    # marginals and the summed right marginal. Discard those exact annihilators.
    mid, targets, axis, rots, r = source
    identity = sp.eye(2)
    labels = [
        identity,
        sp.Matrix([[0, 1], [1, 0]]),
        sp.Matrix([[0, -sp.I], [sp.I, 0]]),
        sp.diag(1, -1),
    ]
    f = list(e)
    for i in range(4):
        for j in range(4):
            p = sp.kronecker_product(labels[i], labels[j])
            coeff = [full.inner(a, p) / 4 for a in e]
            if j == 0:
                subtract = coeff
            elif i == 0:
                subtract = [sum(coeff) / 2] * 2
            else:
                continue
            f = [a - c * p for a, c in zip(f, subtract)]
    norm = iv.sqrt(full.ri(sum(full.inner(a, a) for a in f)))

    def mean(d):
        return full.ri(sum(full.inner(a, b) for a, b in zip(f, d)))

    u = mean(targets["unpair"]) + iv.mpf([-1, 1]) * 3 * iv.sqrt(8) * full.ri(r) * norm
    v = sum(
        a * (mean(d) + iv.mpf([-1, 1]) * full.ri(r) * norm) for a, d in zip(axis, rots)
    )
    return u, v


def weights(ruler, source):
    mid, targets, _, _, r = source
    scale = 2**35
    n = 2 * ruler.order
    lo = []
    hi = []
    physical = []
    for blocks in ruler.menu:
        e = [sp.Matrix(a) / 2 for a in blocks]

        def read(target):
            return sum(full.inner(h, x) for h, x in zip(e, target))

        c = read(mid)
        mu = sum(h.trace() for h in e) / 8
        norm = iv.sqrt(
            full.ri(sum(full.inner(h - mu * sp.eye(4), h - mu * sp.eye(4)) for h in e))
        )
        # Each e is PSD. X >= (1-1e-5) Xmid since Xmid >= .02 I.
        dc = sp.Rational("1e-5") * c
        for end, array, up in ((c + dc, lo, 0), (c - dc, hi, 1)):
            t = sp.Rational(scale) ** n / end
            root, _ = sp.integer_nthroot(t.p // t.q, n)
            array.append(int(root) + up)
        coeff = full.ri(c) + iv.mpf([-1, 1]) * full.ri(r) * norm
        w = coeff ** (-iv.mpf(1) / n)
        # The ordinary tangent is linear in the station with trace-norm
        # change bounded by twice the trace-norm change of the station.
        ordinary = full.ri(read(targets["ordinary"])) + iv.mpf(
            [-1, 1]
        ) * norm * 2 * iv.sqrt(8) * full.ri(r)
        updates = {"ordinary": -w * ordinary / (n * coeff)}
        u, v = variable_coefficients(e, source)
        updates["rotation"] = -w * v / (n * coeff)
        updates["unpair"] = -w * u / (n * coeff)
        physical.append(dict(weight=w, responses=updates, factor=-w / (n * coeff), e=e))
    return np.array(lo, dtype=np.int64), np.array(hi, dtype=np.int64), physical, scale


def path_classes(ruler, lg, hi_matrix, station_group):
    # Integer lower-cost graph; all path sums fit exactly in a binary64 integer.
    assert ruler.n * int(max(lg.data)) < 2**52
    dist = shortest_path(lg, directed=False, indices=station_group)
    rows = {}
    for i, j in combinations(range(len(station_group)), 2):
        origin, target = station_group[i], station_group[j]
        possible = set()

        def visit(p, total, types):
            if p == target:
                possible.add(tuple(sorted(types)))
                return
            for k in range(lg.indptr[p], lg.indptr[p + 1]):
                dest = int(lg.indices[k])
                next_sum = total + int(lg.data[k])
                if next_sum + dist[j, dest] <= hi_matrix[i, j]:
                    kind = int(ruler.index[ruler.ids[tuple(sorted((p, dest)))]])
                    visit(dest, next_sum, types + [kind])

        visit(origin, 0, [])
        assert possible
        rows[i, j] = sorted(possible)
    return rows


def exact_gram(lo, hi, scale, count):
    # Anchor the last point; for a Euclidean fit this Gram must be positive semidefinite.
    r = sp.Matrix((lo[:count, :count] + hi[:count, :count]).astype("int64")) / (
        2 * scale
    )
    d2 = r.applyfunc(lambda x: x * x)
    g = sp.Matrix(
        count - 1,
        count - 1,
        lambda i, j: (d2[i, count - 1] + d2[j, count - 1] - d2[i, j]) / 2,
    )
    # Uniform bound on the perturbation of the squared distances and Gram.
    sqerr = max(
        sp.Rational(int(hi[i, j]) ** 2 - int(lo[i, j]) ** 2, scale**2)
        for i, j in combinations(range(count), 2)
    )
    err = (count - 1) * sp.Rational(3, 2) * sqerr
    safe_lower = g - err * sp.eye(count - 1)
    if all(safe_lower[:k, :k].det() > 0 for k in range(1, count)):
        return "positive"
    _, basis = np.linalg.eigh(np.array(g, dtype=float))
    vec = sp.Matrix(np.rint(100000 * basis[:, 0]).astype(int).tolist())
    if (vec.T * (g + err * sp.eye(count - 1)) * vec)[0] < 0:
        return "negative_direction"
    return "unresolved"


def certify():
    iv.dps = 32
    s = inputs()
    board = optics.pairing.formation.federation(1)
    route = optics.pairing.formation.region(board)
    stations = (*route, 5)
    out = {}
    for order in (3, 4):
        ruler = Ruler(board, route, order)
        low, high, phys, scale = weights(ruler, s)
        dl, _, lg = graph(ruler.n, stations, low[ruler.index], ruler.edges)
        dh, _, _ = graph(ruler.n, stations, high[ruler.index], ruler.edges)
        row = dict(
            four_gram=exact_gram(dl, dh, scale, 4),
            five_gram=exact_gram(dl, dh, scale, 5),
        )
        if order == 3:
            direct = ruler.index[ruler.ids[(0, 39)]]
            assert high[direct] < 2 * low.min()
            # If every minimizing path contains the given intermediate station,
            # the corresponding triangle equality is exact throughout the enclosure.
            for i, j, cut in ((0, 1, 39), (3, 2, 0)):
                keep = np.arange(ruler.n) != cut
                ids = np.flatnonzero(keep)
                alt = shortest_path(
                    lg[keep][:, keep], indices=np.searchsorted(ids, route[i])
                )[np.searchsorted(ids, route[j])]
                assert alt > dh[i, j]
            order_gap = dl[0, 2] + dl[0, 3] + dl[1, 3] - dh[1, 2]
            assert order_gap > scale / 2
            panel = phys[direct]
            known = json.loads(
                (HERE.parent / "null_optics/certificate.json").read_text()
            )["matched_pairing"]
            angles = {"unpair": (-297, -291), "ordinary": (-997, -983)}
            for key, identifier in (
                ("unpair", "angle_interval"),
                ("ordinary", "ordinary_angle_interval"),
            ):
                lo, hi = known[identifier]
                assert angles[key][0] < lo < hi < angles[key][1]
            row.update(
                direct=bounds(panel["weight"]),
                two_link_lower=float(2 * low.min() / scale),
                convex_line_order_gap=float(order_gap / scale),
                matched_slopes={
                    k: bounds(
                        panel["responses"][k]
                        + iv.mpf(v) * panel["responses"]["rotation"]
                    )
                    for k, v in angles.items()
                },
            )
        else:
            candidates = path_classes(ruler, lg, dh, stations[:4])
            entries = {}
            for edge, ways in candidates.items():

                def cost(p):
                    return sum(phys[k]["weight"] for k in p)

                best = min(ways, key=lambda w: float(cost(w).mid))
                for alternate in ways:
                    if alternate == best:
                        continue
                    # Cancel shared weights before comparing intervals.
                    coeff = np.bincount(alternate, minlength=len(phys)) - np.bincount(
                        best, minlength=len(phys)
                    )
                    assert sum(int(n) * p["weight"] for n, p in zip(coeff, phys)).a > 0
                entries[edge] = best

            def path(i, j):
                return entries[tuple(sorted((int(i), int(j))))]

            order_pairs = list(zip(*full.ORDER))
            d = [sum(phys[k]["weight"] for k in path(i, j)) for i, j in order_pairs]
            ga, gv = full.enc.geometry_jet(d)
            derivatives = {}
            for name in ("rotation", "unpair", "ordinary"):
                direction = [
                    sum(phys[k]["responses"][name] for k in path(i, j))
                    for i, j in order_pairs
                ]
                derivatives[name] = [
                    sum(g * t for g, t in zip(jet.g, direction)) for jet in (ga, gv)
                ]
            ra, rv = derivatives["rotation"]
            assert rv.b < 0
            fixed = {}
            for name in ("ordinary", "unpair"):
                a, v = derivatives[name]
                fixed[name] = a - ra * v / rv
            changing = [
                k
                for k, p in enumerate(phys)
                if p["responses"]["unpair"] != 0 or p["responses"]["rotation"] != 0
            ]
            assert len(changing) == 2
            i, j = changing
            pi, pj = phys[i], phys[j]
            u, v = variable_coefficients(pi["e"], s)
            du, dv = variable_coefficients([b - a for a, b in zip(pi["e"], pj["e"])], s)
            directional_det = pi["factor"] * pj["factor"] * (u * dv - du * v)

            def coeffs(jet, k):
                return sum(
                    g * path(a, b).count(k) for g, (a, b) in zip(jet.g, order_pairs)
                )

            geometry_det = coeffs(ga, i) * coeffs(gv, j) - coeffs(ga, j) * coeffs(gv, i)
            fixed["unpair"] = geometry_det * directional_det / rv
            row.update(
                four_volume=bounds(gv.v),
                four_area=bounds(ga.v),
                fixed_area={k: bounds(x) for k, x in fixed.items()},
                unique_path_weights={f"{i},{j}": v for (i, j), v in entries.items()},
            )
            so = iv.mpf([0.4094, 0.4095])
            sd = iv.mpf([0.2310, 0.2312])
            upstream = json.loads(
                (HERE.parent / "full_algebra/certificate.json").read_text()
            )["signed_transfers"]
            for outer, name in ((so, "ordinary_entropy"), (sd, "pairing_entropy")):
                # Preserve decimal endpoints, converting to outward-rounded intervals.
                endpoints = upstream[name].strip("[]").split(",")
                inner = iv.mpf([v.strip() for v in endpoints])
                assert outer.a < inner.a < inner.b < outer.b
            row["held_fraction"] = bounds(
                so * fixed["unpair"] / (sd * fixed["ordinary"])
            )
            assert row["held_fraction"][1] < 1
        out[str(order)] = row
    return out


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    body = json.dumps(certify(), indent=2, allow_nan=False) + "\n"
    if cli.parse_args().write:
        (HERE / "certificate.json").write_text(body)
    else:
        print(body, end="")
