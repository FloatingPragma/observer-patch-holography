import hashlib
import json
import math

import numpy as np
from pathlib import Path

import sympy as sp

from . import static_directions as directions


def test_apertures_and_static_product_curvature():
    rec = directions.report()
    saved = json.loads(Path(directions.__file__).with_suffix(".json").read_text())
    original = saved["native_single_cell"]
    payload = json.dumps(original, sort_keys=True, separators=(",", ":")).encode()
    assert hashlib.sha256(payload).hexdigest() == "5718c48f49c1202cc988dcc8bcd50e60aced56ec3867f9e31430966bc4ad1e6b"
    vectors, _, raw, weighted = directions.design()
    axial = np.ones(12, dtype=int)
    axial[[0, (vectors @ vectors[0]).argmin()]] = 5
    probability = [axial/20, (axial == 1)/10]
    # Simultaneous Bernstein gate on the 48 measured/frozen one-port frequencies.
    log_tail = math.log(2 * 48 / .01)
    for i, expectation in enumerate((0.375, 0.25)):
        row = rec["native_single_cell"][i]
        for sample in (row, original[i]):
            n = sample["shots"]
            bins = np.array(sample["counts_by_port"])
            p = probability[i]
            assert bins.sum() == n == 2**20
            assert np.all(bins[p == 0] == 0)
            assert np.all(abs(bins-n*p) <= np.sqrt(2*n*p*(1-p)*log_tail) + 2*log_tail/3)
            np.testing.assert_array_equal(raw @ bins / n, sample["equal_apertures"])
            np.testing.assert_array_equal(weighted @ bins / (2*n), sample["pole_weighted_apertures"])
        assert abs(row["reference_pole_exact"] - expectation) < 1e-13
    for control in rec["static_product_comparisons"][:2]:
        errors = [r["max_absolute_error"] for r in control["radii"]]
        assert errors[-1] < 2.1e-4
        assert all(a > 3.9 * b for a, b in zip(errors, errors[1:]))
    assert rec["static_product_comparisons"][2]["radii"][-1]["max_absolute_error"] == 0


def test_static_conformal_clock_witness_and_even_decoder():
    x, y, z, s = sp.symbols("x y z s")
    space = x, y, z
    f = s * (x * x - y * y) / 2
    assert sp.expand(f.subs({x: y, y: x}, simultaneous=True) + f) == 0
    hess = sp.hessian(f, space)
    # Local conformal spatial Ricci in dimension three and static lapse Hessian.
    spatial = -hess - sp.trace(hess) * sp.eye(3)
    ric = spatial + sp.trace(hess) * sp.eye(3) - hess
    assert ric == sp.diag(-2 * s, 2 * s, 0)
    ones = sp.ones(6)
    encode = (ones + sp.eye(6)) / 7
    decode = 7 * sp.eye(6) - ones
    assert encode * decode == sp.eye(6)
    assert ones.rank() == 1
    clock_weights = ones / 2 - decode
    assert sum(abs(x) for x in decode[0, :]) == 11
    assert sum(abs(x) for x in clock_weights[0, :]) == 13


def test_first_variation_all_rest_balls_inside_harmonic_core():
    x, y, s, radius = sp.symbols("x y s r", real=True)
    nx, ny = sp.symbols("nx ny")
    field = ((x + s * nx) ** 2 - (y + s * ny) ** 2) / 2
    # Odd monomials average to zero; each spatial coordinate has mean square 1/3.
    spherical_mean = sp.Poly(field, nx, ny)
    moments = {(0, 0): 1, (2, 0): sp.Rational(1, 3), (0, 2): sp.Rational(1, 3)}
    mean = sum(c * moments.get(power, 0) for power, c in spherical_mean.terms())
    volume = 3 * sp.integrate(s**2 * mean, (s, 0, radius))
    moving_boundary = -(radius**2) * sp.integrate(mean, (s, 0, radius))
    assert sp.expand(volume + moving_boundary) == 0
    assert sp.diff(sp.expand(volume + moving_boundary), radius) == 0
    isotropic = s**2
    control = 3 * sp.integrate(s**2 * isotropic, (s, 0, radius))
    control -= radius**2 * sp.integrate(isotropic, (s, 0, radius))
    assert sp.expand(control) == 4 * radius**5 / 15
