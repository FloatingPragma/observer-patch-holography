"""Certify the registered home-cell release over all nonzero quarter-contact aliases."""
import argparse
from fractions import Fraction
import json
from math import acos, cos, isqrt, pi, sqrt

import numpy as np
from scipy import sparse

from model.candidate_laws.interval_echo.openwave.submission import row_04_check as pair
from model.executable_model import carrier
from model.experiments.interval_echo_encounter.check_transport import checked_steps


def setup():
    edges, indices, exchange, adjacent = pair.screen(carrier.HERE / "carriers/FEDERATION_LEVEL0.npz")
    home = edges[np.all(edges // 12 == 13, axis=1)]
    gateways = np.unique(edges[edges[:, 0] // 12 != edges[:, 1] // 12])
    launch = home[np.isin(home, gateways).any(axis=1)][0]
    if tuple(launch) != (156, 157) or len(home) != 30:
        raise ValueError("The registered preparation or terminal record set changed")
    return exchange, adjacent, indices[home[:, 0], home[:, 1]], indices[tuple(launch)]


def gram_bound(matrix, denominator):
    """Exact Schur bound on the squared operator norm of a dyadic matrix."""
    gram = matrix.T @ matrix
    return Fraction(int(np.asarray(abs(gram).sum(axis=1)).max()), denominator**2)


def upper_sqrt(value):
    scale = 2**24
    return Fraction(isqrt(int(value * scale**2)) + 1, scale)


def uniform_tail(exchange, adjacent, patch):
    twice = (2 * exchange).astype(np.int64)
    if (twice != 2 * exchange).nnz:
        raise ArithmeticError("Exchange must be an integer half-matrix")
    embed = sparse.csr_matrix((np.ones(len(patch), dtype=np.int64), (patch, np.arange(len(patch)))),
                              shape=(exchange.shape[0], len(patch)))
    q = sparse.diags((1-adjacent).astype(np.int64), format="csr", dtype=np.int64)
    a = embed.T @ twice @ embed                       # 2 A
    x = -q @ twice @ embed                            # 2 X
    c = embed.T @ twice @ x                           # 4 C
    y = x @ a - q @ twice @ x                         # 4 Y
    if (twice @ embed + x - embed @ a).nnz or (c != c.T).nnz:
        raise ArithmeticError("Home exchange leaves its adjacent block or the correction is not Hermitian")
    r1 = twice @ x - q @ twice @ x - embed @ c        # 4 R1
    r2 = twice @ y - y @ a - x @ c                    # 8 R2
    r3 = -y @ c                                       # 16 R3
    squares = [gram_bound(m, d) for m, d in ((x,2),(y,4),(r1,4),(r2,8),(r3,16))]
    nx, ny, n1, n2, n3 = map(upper_sqrt, squares)
    u = Fraction(1, 2*64)
    error = 2*(u*nx + u**2*ny) + 37*(u*n1 + u**2*n2 + u**3*n3)
    if error >= 1 or 1-error**2 <= Fraction(94, 100):
        raise ArithmeticError("The large-strength home bound missed its probability gate")
    return dict(squared_norm_bounds=[str(s) for s in squares], lower_bound=float(1-error**2))


def finite_window(exchange, adjacent, patch, launch, strength):
    h = exchange - 2*strength*sparse.diags(adjacent, format="csr")
    column = h.getcol(launch)
    variance = float(column.multiply(column).sum() - h[launch, launch]**2)
    if variance != 2.25:
        raise ArithmeticError("The registered launch energy variance changed")
    start = np.zeros(h.shape[0], complex)
    start[launch] = 1
    grid_floor = 1.
    readings = {}
    # One time unit of the scaled matrix is 1/8 of a native exposure.
    for tick, state, error in checked_steps(h / 8, start, range(1, 297), degree=52, steps_per_unit=1):
        if not 0 <= error < 1e-5:
            raise ArithmeticError("The Taylor probability enclosure exceeded its budget")
        local = abs(state[patch])**2
        home = float(local.sum())
        grid_floor = min(grid_floor, home-error)
        if tick in (8, 56, 296):
            readings[tick//8] = dict(seams=local.tolist(), quiet=float(np.vdot(state, state).real-home),
                                     probability_error=error)
    angle = acos(sqrt(grid_floor)) + sqrt(variance)/16
    lower = cos(min(pi/2, angle))**2
    if strength and lower <= 0.940001:
        raise ArithmeticError(f"Strength {strength} missed the continuous home-cell gate: {lower}")
    if strength == 0 and 1-readings[37]["quiet"] + readings[37]["probability_error"] >= 0.65:
        raise ArithmeticError("Exchange alone did not clear the rejection control")
    return dict(strength=strength, continuous_lower=lower, readings=readings)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true", help="also propagate all fourteen finite aliases and zero")
    args = parser.parse_args()
    exchange, adjacent, patch, launch = setup()
    print(json.dumps(dict(tail=uniform_tail(exchange, adjacent, patch))), flush=True)
    if args.full:
        for k in range(-56, 57, 8):
            print(json.dumps(finite_window(exchange, adjacent, patch, launch, k), allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
