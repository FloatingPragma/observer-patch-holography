"""Does sweeping the shared zone reconcile the four-tick native currents?

P(t) = P_left + f(t) P_zone, with f(4)=1/2 and v=f'(4). For H_R={P,H}/2,
the outward energy rate is J_E-v E_zone at that epoch. Adding cI to H gives
J_E+c J_N-v E_zone-vc N_zone. All four meters use the original pump and the
490-port shared zone, with no change to the source or any radar reading.
Finite quarter-turn preparations preserve the original two head marginals.
The two additional derivative rows differentiate these physical rotations
at zero angle in the existing native angular unit.

A nuisance-only cofactor identity constructs two physical source mixtures
with all three sweep/origin terms exactly matched in population. The
probabilities use the true three meters. They are computable reals because
the finite Hamiltonian, its exponential and the gapped marginal axes are
computable and the normalization is dyadic. A fair-bit threshold lottery
which refines the probability enclosures terminates almost surely. This
argument applies to one spatially uniform, preparation-independent zone
speed. An independently measured velocity field would need its own test.
"""

import argparse
from functools import lru_cache
import json
import math
from itertools import permutations

import numpy as np
import sympy as sp
from mpmath.ctx_iv import MPIntervalContext

from . import current, momentum_screen as face, observability
from .center_average import shells
from .. import live_modular, live_source
from ..source import certificate as cert

COLUMN = ("energy_current", "number_current", "zone_energy", "zone_number")
FINITE = ("ordinary", "decorrelated", "left_quarter", "right_quarter")
LOCAL = ("ordinary", "decorrelated", "left", "right")


def preparations():
    (ref, prod), ferr, _ = live_source.formation_center()
    out = {
        "ordinary": ((-1, 1), ref, ferr, 1),
        "decorrelated": ((2, 1), prod - ref, 3 * ferr, 2),
    }
    for side, (name, (drot, err)) in enumerate(
        live_modular.rotation_centers(ref, ferr).items()
    ):
        assert name == ("left", "right")[side]
        out[name] = ((2, 1), drot, err, 1)
        g, ge = observability.marginal_generator(ref, ferr, side)
        alt = observability.quarter_turn(ref, g, 1)
        # U=cos(pi/4)I-i sin(pi/4)2G. The norm error in U is <=2 ge.
        # quarter_turn uses G in the expanded conjugation, preserving this
        # bound even when the center itself is not an exact involution.
        ue = 2 * ge
        err = 2 * ferr + (1 + ferr) * (2 * ue + ue**2)
        out[name + "_quarter"] = ((2, 1), alt - ref, err, 2)
    return out


def controls():
    arms = preparations()
    src = live_source.source
    board = src.formation.federation(1)
    left, shared, right = face.zones(board)
    ls = current.current_effects(np.flatnonzero(left).tolist())
    rs = current.current_effects(np.flatnonzero(right).tolist())
    out = {s: {} for s in arms}
    for k in ("energy", "number"):
        op = (ls[k][0] - rs[k][0]) / 2
        limit = (ls[k][1] + rs[k][1]) / 2
        pump_error = (ls[k][2] + rs[k][2]) / 2
        for s, (rail, head, err, norm) in arms.items():
            sigma = sp.kronecker_product(sp.diag(*rail) / 3, head)
            center = sp.re(sp.trace(op * sigma))
            out[s][k + "_current"] = cert.enclosed(
                center, limit * err + (norm + err) * pump_error
            )
    met = shells.meters(board, {"zone": [np.flatnonzero(shared).tolist()]}, arms)[
        "zone"
    ]
    for row in met:
        for s in arms:
            out[s]["zone_" + row["meter"]] = row["inputs"][s]
    return out


def determinant(rows):
    n = len(rows)
    value = 0
    for order in permutations(range(n)):
        sign = (-1) ** sum(
            order[a] > order[b] for a in range(n) for b in range(a + 1, n)
        )
        term = sign
        for i, j in enumerate(order):
            term *= rows[i][j]
        value += term
    return value


def ends(x):
    return [
        float(np.nextafter(float(x.a), -np.inf)),
        float(np.nextafter(float(x.b), np.inf)),
    ]


def solve(data, arms):
    ctx = MPIntervalContext()
    ctx.dps = 50
    a = [[ctx.mpf(data[r][m]) for m in COLUMN] for r in arms]
    # Signed minors of the three nuisance columns depend on no energy current.
    w = [
        (-1) ** i * determinant([row[1:] for j, row in enumerate(a) if j != i])
        for i in range(4)
    ]
    det = sum(x * row[0] for x, row in zip(w, a, strict=True))
    l1 = sum(abs(x) for x in w)
    if abs(det).a <= 0 or l1.a <= 0:
        return dict(four_row_determinant=ends(det), rank_four=False)
    out = dict(
        four_row_determinant=ends(det),
        rank_four=True,
        any_origin_and_sweep_max_absolute_residual_lower=ends(abs(det) / l1)[0],
    )
    if arms != FINITE:
        return out
    # The exact weights are the corresponding true minors, affinely normalized.
    # Use an outward dyadic upper bound for either mass; unused weight prepares
    # the same reference in each source. Sign orientation uses the first
    # nuisance minor; the energy read is held out of the construction.
    assert all(x.a * x.b > 0 for x in w)
    if w[0].b < 0:
        w = [-x for x in w]
        det = -det
    splus = sum(x for x in w if x.a > 0)
    sminus = sum(-x for x in w if x.b < 0)
    # A 44-bit mesh is fine here because the minors are physical/source
    # rather than arbitrarily small arithmetic residuals.
    denominator = 1 << 44
    upper = max(ends(splus)[1], ends(sminus)[1])
    numerator = math.ceil(upper * denominator)
    norm = ctx.mpf(numerator) / denominator
    assert splus.b <= norm.a and sminus.b <= norm.a
    lotteries = {}
    for tag, sign in (("positive", 1), ("negative", -1)):
        p = {
            s: ends(sign * x / norm) if (sign * x).a > 0 else [0, 0]
            for s, x in zip(arms, w, strict=True)
        }
        mass = splus if sign == 1 else sminus
        p["reference"] = [max(0.0, ends(1 - mass / norm)[0]), ends(1 - mass / norm)[1]]
        lotteries[tag] = p
    out.update(
        lotteries=lotteries,
        cofactor_normalization=[numerator, denominator],
        exactly_matched=["number_current", "zone_energy", "zone_number"],
        signed_energy_current_gap=ends(det / norm),
    )
    return out


@lru_cache(maxsize=1)
def run():
    data = controls()
    return dict(
        pump_age=4,
        shared_weight=0.5,
        shared_ports=490,
        columns=COLUMN,
        contrasts=data,
        finite_preparations=solve(data, FINITE),
        infinitesimal_rotations=solve(data, LOCAL),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    result = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (face.hinge.HERE / "moving_cut.json").write_text(result)
    else:
        print(result, end="")
