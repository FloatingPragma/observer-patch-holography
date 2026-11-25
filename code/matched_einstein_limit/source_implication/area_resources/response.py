"""Area response resources for the unflagged formed-collar comparison bit."""

from math import atanh, exp, log1p, pi, sqrt

import numpy as np
from scipy.integrate import quad
from scipy.stats import binom

from ...collar_identifiability.refinement import admission
from ...collar_identifiability.refinement.boundary_sources import absolute_flux


def native_input():
    source = admission.source
    board = source.formation.federation(1)
    route = source.formation.region(board)
    joint = source.prepared_states(board, route, 1)[0]
    joint = joint / np.trace(joint)
    sigma = np.kron(*source.formation.marginals(joint))
    dmatrix = joint - sigma
    values, axes = np.linalg.eigh(sigma)
    relative = (axes.conj().T @ dmatrix @ axes) / np.sqrt(values[:, None]*values[None, :])
    half_width = 1 / np.max(np.abs(np.linalg.eigvalsh(relative)))
    outputs = {}
    for strength in (8, 0):
        record = admission.material_panel(1, 1, strength)["records"][-1]
        z = record["d"][-1]
        fresh = record["q"][-1]
        p = admission.compare_returns(fresh + z, fresh)
        np.testing.assert_allclose(2*p - 1, z, atol=1e-15, rtol=0)
        outputs[strength] = z
    return float(half_width), outputs


def curvature_bound(copies, half_width, x, limit=1):
    if copies < 1 or abs(x) >= half_width or limit < 0:
        raise ValueError("invalid copy count, output bound or position on the admissible chord")
    return limit * sqrt(2 * copies * (copies-1)) / (half_width**2 - x**2)


def latent_extremum(copies, half_width, x):
    """Tight bound if the experiment could read all hidden endpoint choices."""
    k = np.arange(copies+1)
    plus, minus = 1/(half_width+x), -1/(half_width-x)
    kernel = (k*plus + (copies-k)*minus)**2 - (k*plus**2 + (copies-k)*minus**2)
    mass = binom.pmf(k, copies, (1 + x/half_width)/2)
    return float(mass @ np.abs(kernel)), float(mass @ kernel**2)


def majority(bank, y):
    if bank < 1 or bank % 2 == 0 or abs(y) >= 1:
        raise ValueError("the majority bank needs a positive odd size and interior bit bias")
    b = absolute_flux(bank)
    k = (bank-1)//2
    derivative = b * exp(k * log1p(-y*y))
    mean = b * y * quad(lambda v: exp(k * log1p(-(y*v)**2)), 0, 1)[0]
    return mean, derivative


def majority_stationarity(inverse_radius, z, s, eta=1):
    """A fixed positive gain is chosen to match one scalar on the complete source interval."""
    n = inverse_radius
    ell = 1/n
    bank = n*n if n%2 else n*n + 1
    f = admission.any_hit(admission.sphere_count(n)/n**6, n)
    y = s * f * z
    entropy_slope = -f*z*atanh(y)
    mu, prime = majority(bank, y)
    gain = pi / (4*eta)
    area = gain * ell**2 * (1 + mu**2)
    calibrated_area_slope = eta * gain * ell**2 * 2*mu*prime*f*z
    return dict(ell=ell, f=f, bank=bank, area=area, entropy=entropy_slope,
                calibrated_area=calibrated_area_slope,
                residual=(entropy_slope+calibrated_area_slope)/ell**4)


if __name__ == "__main__":
    width, signals = native_input()
    print(f"symmetric input chord {width:.9g}, measured contrasts {signals}")
    for n in (8, 16, 32, 64):
        rows = [majority_stationarity(n, signals[8], s) for s in (.25,.75)]
        ell = 1/n
        spread = (rows[0]["entropy"] - rows[1]["entropy"]) / ell**4
        residuals = tuple(row["residual"] for row in rows)
        print(f"ell=1/{n}, copies {2*rows[0]['bank']}, entropy slope spread/ell^4 "
              f"{spread:.9g}, residuals/ell^4 {residuals}")
