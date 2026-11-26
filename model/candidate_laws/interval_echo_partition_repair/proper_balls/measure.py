"""Compatible boundary words on volume-matched balls in ten supplied rest frames."""
import argparse
from decimal import Decimal, localcontext
from functools import cache
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

HERE = Path(__file__).resolve().parent
PI = Decimal("3.14159265358979323846264338327950288419716939937510582097494")
ROUND_RADIUS = 0.1 / (math.pi / 2 - math.pi / 20)
SCALES = (256, 1024, 4096)
AXES = dict(rest=(0, 0, 0), x=(1, 0, 0), y=(0, 1, 0), z=(0, 0, 1),
            xy=(1, 1, 0), xz=(1, 0, 1), yz=(0, 1, 1))
FRAMES = dict(rest="rest", **{s+k: k for k in "xyz" for s in ("+", "-")},
              xy="xy", xz="xz", yz="yz")


def lift(name):
    axis = np.array(AXES[name], dtype=float)
    gamma = {0: 1., 1: 5/3, 2: 3.}[int(axis.sum())]
    if gamma > 1:
        axis /= np.linalg.norm(axis)
    return gamma, np.eye(3) + (gamma-1)*np.outer(axis, axis)


@cache
def sampling(order):
    z, w = leggauss(order)
    phis = np.arange(2*order)*np.pi/order
    z, phi = (v.ravel() for v in np.broadcast_arrays(z[:, None], phis[None, :]))
    s, c, h = np.sqrt(1-z*z), np.cos(phi), np.sin(phi)
    n = np.column_stack((z, s*c, s*h))
    a = np.column_stack((-s, z*c, z*h))
    b = np.column_stack((np.zeros_like(z), -h, c))
    # We integrate the correction to the constant, normalized surface density.
    return n, a, b, np.repeat(w, 2*order)/(4*order), (leggauss(order))


def inner(a, b):
    return (a*b).sum(axis=1)


class Survey:
    """Tangential form of the chamber metric; subtract the flat answer analytically."""
    def __init__(self, name, sign, order):
        n, a, b, self.angular, (u, w) = sampling(order)
        self.nodes, self.weights = (u+1)/2, w/2
        _, L = lift(name)
        X, A, B = n@L, a@L, b@L
        q2 = inner(X, X)
        self.q2 = q2
        C, D = np.cross(X, A), np.cross(X, B)
        diag = np.array((2+sign, 2-sign, 2))
        self.round = np.array([inner(U, V)-inner(U, X)*inner(V, X)/q2
                               for U, V in ((A, A), (B, B), (A, B))])
        self.chamber = np.array([inner(U*diag, V) for U, V in ((C, C), (D, D), (C, D))])

    def surface_excess(self, radius):
        if radius*math.sqrt(max(self.q2)) >= 0.02:
            raise ValueError("outside the constant-cutoff chamber")
        u2 = radius**2 * self.q2 / ROUND_RADIUS**2
        # sinc(sqrt(u2))^2-1. The next Taylor term is smaller than 1e-18
        # throughout |x|<.02; the largest actual |x| here is below .012.
        correction = sum((-1)**j*2**(2*j+1)*u2**j/math.factorial(2*j+2)
                         for j in range(1, 8))
        aa, bb, ab = self.round*correction + self.chamber*radius**2
        det_excess = aa+bb+aa*bb-ab*ab
        return float(self.angular @ (det_excess/(1+np.sqrt(1+det_excess))))

    def ball(self, denominator):
        ell = 1/denominator
        # Write r/ell=1+epsilon and solve the small correction itself. Gauss'
        # lemma gives dV/dr=A, including the chamber term.
        def vol_excess(eps):
            r = ell*(1+eps)
            excess = sum(w*u*u*self.surface_excess(r*u)
                         for u, w in zip(self.nodes, self.weights, strict=True))
            return eps*(3+eps*(3+eps)) + 3*(1+eps)**3*excess

        eps = 0.
        for _ in range(4):
            eps -= vol_excess(eps)/(3*(1+eps)**2*(1+self.surface_excess(ell*(1+eps))))
        radius = (1+Decimal.from_float(eps))/denominator
        area = 4*PI*radius**2*(1+Decimal.from_float(self.surface_excess(ell*(1+eps))))
        return radius, area, abs(vol_excess(eps))


def multiplicities(denominator, areas):
    M = denominator**5
    plus, minus = (int(M*a.exp()) for a in areas)
    # Each metric has N equiprobable legal pairs; d distinct completions
    # per endpoint. No arrays of size N are needed.
    N = math.lcm(plus, minus)*max(plus, minus)**2
    return M, plus, minus, N


def word_from_rank(rank, d, N):
    """Read a uniform legal pair using a uniform integer rank, with either metric."""
    if N % d or N < d*d or not 0 <= rank < N:
        raise ValueError("invalid compatible-word parameters")
    i, k = divmod(rank, d)
    return i, (i+k) % (N//d)


def direction_bank():
    j = np.arange(384)
    z = (2*j+1)/len(j)-1
    phi = math.pi*(3-math.sqrt(5))*j
    return np.column_stack((np.sqrt(1-z*z)*np.cos(phi), np.sqrt(1-z*z)*np.sin(phi), z))


def shutters(rows):
    directions = direction_bank()
    bands = []
    for row in rows:
        inverse = np.linalg.inv(lift(row["shape"])[1])
        norm = np.linalg.norm(directions @ inverse, axis=1)
        rr = [float(row[s]["radius"]) for s in ("plus", "minus")]
        bands.append(np.stack((min(rr)/norm, max(rr)/norm), axis=1))
    for i, row in enumerate(rows):
        low, high = bands[i].T
        clearance = low.copy()
        for j, band in enumerate(bands):
            if i != j:
                clearance = np.minimum(clearance, np.maximum(low-band[:, 1], band[:, 0]-high))
        k = int(np.argmax(clearance/low))
        if clearance[k] <= 0:
            raise ValueError("no isolated crossing")
        # Round coordinates much less than the verified clearance. The frozen
        # positions are used directly when the two memberships are replayed.
        endpoints = np.array((low[k]-clearance[k]/4, high[k]+clearance[k]/4))
        row["worldlines"] = (endpoints[:, None]*directions[k]).round(14).tolist()
        row["relative_radial_clearance"] = float(clearance[k]/low[k])
    return rows


def memberships(rows):
    whole = []
    x = np.array([r["worldlines"] for r in rows])
    for i, row in enumerate(rows):
        L = lift(row["shape"])[1]
        distance = np.linalg.norm(x @ np.linalg.inv(L), axis=2)
        plus, minus = (distance < float(row[label]["radius"]) for label in ("plus", "minus"))
        assert np.array_equal(plus, minus)
        assert np.array_equal(np.where(plus[:, 0] != plus[:, 1])[0], [i])
        assert plus[i, 0] and not plus[i, 1]
        whole.append(np.flatnonzero(plus[:, 0] & plus[:, 1]).tolist())
    return whole


def measure():
    with localcontext() as ctx:
        ctx.prec = 55
        rows = []
        for denominator in SCALES:
            for name in AXES:
                row = {"shape": name, "inverse_ell": denominator}
                areas = []
                comparisons = []
                residuals = []
                for sign, label in ((1, "plus"), (-1, "minus")):
                    radius, area, residual = Survey(name, sign, 14).ball(denominator)
                    _, check_area, _ = Survey(name, sign, 9).ball(denominator)
                    areas.append(area)
                    comparisons.append(check_area)
                    residuals.append(residual)
                    row[label] = {"radius": str(radius), "area": str(area)}
                M, dp, dm, N = multiplicities(denominator, areas)
                assert (M, dp, dm, N) == multiplicities(denominator, comparisons)
                # Antisymmetric area vanishes under (x,y) exchange for these planes.
                if name in ("rest", "z", "xy"):
                    assert dp == dm
                logratio = (Decimal(dp)/dm).ln()
                da = areas[0]-areas[1]
                rp, rm = (Decimal(row[k]["radius"]) for k in ("plus", "minus"))
                # A fixed homogeneous line flux on convex projected balls scales
                # as r^2. Give it a separate correct flat calibration per frame.
                crofton_area_change = 4*PI*(rp-rm)*(rp+rm)
                norm = denominator**4
                floor_margins = []
                for area in areas:
                    v = M*area.exp()
                    f = v-int(v)
                    floor_margins.append(min(f, 1-f))
                row.update(M=M, legal_pairs=N, plus_d=dp, minus_d=dm,
                           min_word_floor_margin=float(min(floor_margins)),
                           biased_reference_residual_nats=float(da-logratio-Decimal(2).ln()),
                           area_slope_over_ell4=float(da*norm),
                           homogeneous_line_area_slope_over_ell4=float(crofton_area_change*norm),
                           entropy_slope_over_ell4=float(-logratio*norm),
                           balance_over_ell4=float((da-logratio)*norm),
                           swapped_balance_over_ell4=float((da+logratio)*norm),
                           M_times_max_area_order_change=float(max(abs(a-b)*M for a, b in zip(areas, comparisons, strict=True))),
                           max_relative_volume_residual=max(residuals))
                rows.append(row)
        shutters(rows)
        for row, counts in zip(rows, memberships(rows), strict=True):
            row["whole_bond_indices"] = counts
        return {"inverse_scales": list(SCALES), "frames": FRAMES, "gain": 1, "rows": rows}


def verify_saved(output):
    """Replay archived counts, the geometry-blind control and shared memberships."""
    archived = json.loads((HERE / "results.json").read_text())
    assert archived["frames"] == output["frames"] and archived["inverse_scales"] == output["inverse_scales"]
    assert memberships(archived["rows"]) == [r["whole_bond_indices"] for r in archived["rows"]]
    for old, new in zip(archived["rows"], output["rows"], strict=True):
        for field in ("shape", "inverse_ell", "plus_d", "minus_d", "legal_pairs"):
            assert old[field] == new[field]
        for field in ("area_slope_over_ell4", "balance_over_ell4", "homogeneous_line_area_slope_over_ell4"):
            assert abs(old[field]-new[field]) < 1e-9


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = measure()
    if args.write:
        (HERE / "results.json").write_text(json.dumps(output, indent=2, allow_nan=False)+"\n")
    else:
        verify_saved(output)
        print("Verified 30 frame/radius surveys and their shared compatible-word branches.")
