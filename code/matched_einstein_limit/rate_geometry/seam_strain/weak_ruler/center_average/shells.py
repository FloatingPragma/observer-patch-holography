"""Radial sources for the individually volume-matched area coefficient.

The measured response is dC in W(4 pi ell^2 - <A>) = C ell^4, with the
observer sampling volume W fixed and every ball at volume 4 pi ell^3/3.
C has radar length units. Number is the spent pump projection P, and energy
is (PH+HP)/2 in exchange units. A common energy-origin shift adds cP to
each energy meter. The ordinary rail direction and nine Pauli-product/4
correlations span the ten calibrated source directions.

With four physical variations, four independently signed coefficients fit
if the meter matrix is invertible. A common shifted-energy profile further
requires N0*E1 - N1*E0 = 0 for these four fit coefficients. Over the nine
fixed-marginal inputs the meter rows have a kernel. An outward inverse
bound constructs an exact direction in this kernel for which dC is positive.
It excludes even a differentiable function of the enrolled shell readings.
The direction is normalized to maximum input coefficient at most one.

Two positive lotteries (quarter turns about the actual two marginal axes
and independent re-pairing) implement a scaled copy of that difference.
At common admixture t into the formed reference their shell energies and
numbers agree exactly for every t by linearity. Their area difference per
unit admixture has the reported positive derivative at zero. Thus for all
sufficiently small positive t the same material inputs give different C.
This last conclusion excludes even nondifferentiable single-valued maps
from these particular shell readings to C. The construction
concerns the enrolled graph sets and supplied radar filling. Identifying
their shells with proper-radius shells of a spacetime ball is a separate
physical condition.
"""

import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np
import sympy as sp
from mpmath import iv
from scipy.linalg import qr

from . import native
from .. import aperture, hinge, jet_lottery
from ... import live_bounds, live_modular, live_radar, live_source
from ...source import certificate as cert, measure as source

HERE = Path(__file__).parent
AXES = tuple(a + b for a, b in product("XYZ", repeat=2))
INPUTS = ("ordinary", *AXES)
Q = sp.Rational


def axes_and_preparations():
    axes = dict(zip(AXES, live_source.qubit_hidden.correlation_axes(), strict=True))
    (ref, prod), err, _ = live_source.formation_center()
    arms = {
        "reference": ((2, 1), ref, err, 1),
        "ordinary": ((-1, 1), ref, err, 1),
        "decorrelated": ((2, 1), prod - ref, 3 * err, 2),
    }
    arms.update(
        {
            t: ((2, 1), rho, e, 1)
            for t, (rho, e) in live_modular.rotation_centers(ref, err).items()
        }
    )
    arms.update({t: ((2, 1), rho, 0, 1) for t, rho in axes.items()})
    return arms


def shell_groups(board):
    cut = aperture.enroll(board, 1)[0]
    distance, front = {0: 0}, {0}
    while not set(cut) <= distance.keys():
        front = {n for p in front for n in board.neighbours[p]} - distance.keys()
        for n in front:
            distance[n] = 1 + min(
                distance[v] for v in board.neighbours[n] if v in distance
            )
    central = set(source.formation.region(board)[1:3])
    return {
        "zero_and_aperture_remainder": [[0], sorted(set(cut) - {0})],
        "bridge_and_neighbors": [sorted(central), sorted(set(cut) - central)],
        "distance_from_zero_within_aperture": [
            [p for p in cut if distance[p] == r]
            for r in range(1 + max(distance[p] for p in cut))
        ],
        "graph_ball_about_zero": [[0], sorted(board.neighbours[0])],
    }


def meters(board, groups, preparations):
    hi, lo, colerr, _ = live_source.pump_center(True)
    re, im = np.hsplit(hi.astype(object) * (1 << 36) + lo, 2)
    h = source.system(board, source.formation.region(board), -1)[0]
    norm = Q(int(abs(h).sum(axis=1).max()), 40)
    out = {}
    for name, cells in groups.items():
        rows = []
        for ports in cells:
            inside = live_radar.station_rows(board, ports).ravel()
            sub = h[inside].tocoo()
            hr, hi_part = [np.zeros((len(inside), 8), dtype=object) for _ in range(2)]
            for part, w in ((hr, re), (hi_part, im)):
                np.add.at(part, sub.row, sub.data[:, None] * w[sub.col])

            def pair(br, bi):
                return (
                    sp.Matrix(re[inside].T @ br + im[inside].T @ bi)
                    + sp.I * sp.Matrix(re[inside].T @ bi - im[inside].T @ br)
                ) / (1 << 164)

            number = pair(re[inside], im[inside])
            ph = pair(hr, hi_part)
            for label, op, bound in (
                ("number", number, 1),
                ("energy", (ph + ph.H) / 80, norm),
            ):
                inputs = {}
                for tag, (weights, rho, err, tn) in preparations.items():
                    state = sp.kronecker_product(sp.diag(*weights) / 3, rho)
                    mid = sp.re(
                        sum(op[i, j] * state[j, i] for i in range(8) for j in range(8))
                    )
                    uncertainty = bound * (
                        err + (tn + err) * (6 * colerr + 9 * colerr**2)
                    )
                    inputs[tag] = cert.enclosed(mid, uncertainty)
                rows.append(dict(ports=ports, meter=label, inputs=inputs))
        out[name] = rows
    return out


def area_axes(preparations, archive):
    edges, den, err = hinge.exact_ranges()
    _, _, col, g = live_source.pump_center(True)
    e = 4 * col + 4 * col**2
    states = {
        t: (live_source.rational_read(g, rho, weights), inp + (1 + inp) * e)
        for t, (weights, rho, inp, tn) in preparations.items()
        if t in (*AXES, "reference")
    }
    inputs = hinge.paired_inputs(edges, den, hinge.edge_errors(edges, den, err), states)
    lengths = inputs.pop("reference")
    out = dict(archive["individual_area_slopes"])
    for batch in [AXES[:5], AXES[5:]]:
        names = [*batch, *[batch[-1]] * (5 - len(batch))]
        radial = hinge.radial_jets(
            lengths,
            {key: inputs[t] for t, key in zip(names, live_bounds.KEYS, strict=True)},
        )
        for k, tag in enumerate(batch):
            total = iv.mpf(0)
            for (length, angle), row in zip(
                radial, archive["radial_edges"], strict=True
            ):
                coef = iv.mpf(row["individually_fixed_volume_area_coefficient"])
                deriv = iv.mpf(row["individually_fixed_volume_angle_derivative"])
                total += (coef * length.g[k] + length.v * deriv * angle.g[k]) / 3
            out[tag] = live_bounds.endpoints(total)
    return out


def table(rows, tags):
    m = sp.Matrix([[sum(map(Q, r[t])) / 2 for t in tags] for r in rows])
    e = sp.Matrix([[(Q(r[t][1]) - Q(r[t][0])) / 2 for t in tags] for r in rows])
    return jet_lottery.round_center(m, e, bits=68)


def certificate(rows, area, tags=INPUTS):
    # Choose a well-conditioned material minor, then build an EXACT kernel
    # vector using the true physical minor. All uncertainty arithmetic below
    # is rational and encloses the change from the computed center.
    M, E = table([r["inputs"] for r in rows], tags)
    A, F = table([area], tags)
    v = np.array(M, float)
    piv = qr(v / np.max(abs(v), axis=1)[:, None], pivoting=True, mode="economic")[2][
        : M.rows
    ].tolist()
    B = M[:, piv]
    if B.det() == 0:
        return dict(status="unresolved rank", columns=list(tags))
    inv = B.inv()
    eta = jet_lottery.inf_norm(inv.applyfunc(abs) * E[:, piv])
    if eta >= 1:
        return dict(
            status="unresolved inverse", columns=list(tags), inverse_defect=float(eta)
        )
    candidates = []
    for j in range(len(tags)):
        if j in piv:
            continue
        c = -inv * M[:, j]
        dc = jet_lottery.inf_norm(
            inv.applyfunc(abs) * (E[:, j] + E[:, piv] * c.applyfunc(abs))
        ) / (1 - eta)
        # Normalize by a common rational upper bound on max |head/rail weight|.
        scale = max(sp.S.One, max(map(abs, c)) + dc)
        w = sp.zeros(len(tags), 1)
        w[j] = 1
        for i, p in enumerate(piv):
            w[p] = c[i]
        err = (F * w.applyfunc(abs))[0] + sum(abs(A[0, p]) + F[0, p] for p in piv) * dc
        val = (A * w)[0]
        if val < 0:
            w = -w
            val = -val
        candidates.append(
            (
                float((val - err) / scale),
                j,
                w / scale,
                dc / scale,
                val / scale,
                err / scale,
            )
        )
    lower, j, w, dc, val, err = max(candidates, key=lambda x: x[0])
    return dict(
        status="separated" if lower > 0 else "unresolved area",
        columns=list(tags),
        material_rows=len(rows),
        pivot_columns=[tags[i] for i in piv],
        inverse_defect_upper=cert.enclosed(eta, 0)[1],
        free_column=tags[j],
        direction={
            t: cert.enclosed(w[i], dc if i in piv else 0) for i, t in enumerate(tags)
        },
        individually_fixed_volume_area_slope=cert.enclosed(val, err),
    )


def physical_arms(rows, area):
    tags = native.ARMS
    M, E = table([r["inputs"] for r in rows], tags)
    A, F = table([area], tags)
    T = M.T
    inv = T.inv()
    bound = inv.applyfunc(abs)
    defect = bound * E.T
    if jet_lottery.inf_norm(defect) >= 1:
        return dict(status="unresolved inverse")
    x = inv * A.T
    rad = (sp.eye(4) - defect).inv() * bound * (F.T + E.T * x.applyfunc(abs))
    coefficients = [cert.enclosed(u, e) for u, e in zip(x, rad, strict=True)]
    n0, h0, n1, h1 = map(iv.mpf, coefficients)
    return dict(
        status="four free signed coefficients fit four inputs",
        arms=list(tags),
        number_energy_coefficients=coefficients,
        common_origin_minor=live_bounds.endpoints(n0 * h1 - n1 * h0),
    )


def native_lotteries(witness):
    names, bank, entry_error = jet_lottery.material()
    entry_error = Q(sp.ceiling(entry_error * 2**75), 2**75)
    bank, errors = jet_lottery.round_center(bank, sp.ones(9) * entry_error, bits=75)
    inv = bank.inv()
    raw = sp.Matrix([list(map(Q, witness["direction"][a])) for a in AXES])
    center = (raw[:, 0] + raw[:, 1]) / 2
    width = (raw[:, 1] - raw[:, 0]) / 2
    w = inv * center
    defect = inv.applyfunc(abs) * errors
    assert jet_lottery.inf_norm(defect) < 1
    radius = (
        (sp.eye(9) - defect).inv()
        * inv.applyfunc(abs)
        * (width + errors * w.applyfunc(abs))
    )
    scale = (
        sp.ceiling(
            max(sum(max(0, sg * v + e) for v, e in zip(w, radius)) for sg in (1, -1))
        )
        + 1
    )
    weights = {}
    for sign, name in ((1, "plus"), (-1, "minus")):
        low = [max(0, sign * v - e) / scale for v, e in zip(w, radius)]
        high = [max(0, sign * v + e) / scale for v, e in zip(w, radius)]
        parts = [*zip(low, high), (1 - sum(high), 1 - sum(low))]
        weights[name] = {
            k: [max(0, cert.enclosed(a, 0)[0]), min(1, cert.enclosed(b, 0)[1])]
            for k, (a, b) in zip([*names, "reference"], parts, strict=True)
        }
    slope = iv.mpf(witness["individually_fixed_volume_area_slope"]) / int(scale)
    return dict(
        common_scale=int(scale),
        probabilities=weights,
        area_slope_difference=live_bounds.endpoints(slope),
    )


def run():
    archive = json.loads((HERE / "native.json").read_text())
    board = source.formation.federation(1)
    groups = shell_groups(board)
    p = axes_and_preparations()
    area = area_axes(p, archive)
    material = meters(board, groups, p)
    comparisons = {}
    for name, rows in material.items():
        comparisons[name] = {
            "meters": rows,
            "ten_input_certificate": certificate(rows, area),
            "fixed_marginals_certificate": certificate(rows, area, AXES),
        }
        if len(rows) == 4:
            comparisons[name]["four_physical_arms"] = physical_arms(rows, area)
        witness = comparisons[name]["fixed_marginals_certificate"]
        if witness["status"] == "separated":
            comparisons[name]["same_marginal_lotteries"] = native_lotteries(witness)
    return dict(
        preparation="Same four-tick, coupled 80-cell vacancy pump and formed heads.",
        normalization="C has radar length units in W times (Euclidean area minus individually matched mean area) = C times radius^4. Raw kernel directions have maximum Pauli-product/4 or rail coefficient at most one. Lottery responses are per unit admixture around the formed reference.",
        area_slopes=area,
        groups=comparisons,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if parser.parse_args().write:
        (HERE / "shells.json").write_text(data)
    else:
        print(data, end="")
