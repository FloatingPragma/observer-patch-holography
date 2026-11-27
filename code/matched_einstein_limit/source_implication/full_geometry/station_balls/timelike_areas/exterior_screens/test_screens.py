"""Checks from the native banks and a direct projected-metric area integration."""
import json
from math import cos, pi, sin, sqrt

import numpy as np
from numpy.polynomial.legendre import leggauss

from ..measure import ROUND_RADIUS
from . import measure as m


def test_pointer_statistics_and_interference_identity():
    archive = json.loads((m.HERE/"results.json").read_text())
    m.verify(archive)
    t = m.probabilities()
    blank, full, half = t.sum(axis=1)
    np.testing.assert_allclose(half, (full+blank)/2, atol=5e-14, rtol=0)
    np.testing.assert_allclose(t[2, 1], .75*t[1, 1], atol=5e-14, rtol=0)
    assert archive["summary"]["simultaneous_interval"][0] > .01217
    assert archive["summary"]["resolved_indices"] == [0, 1]


def pixels(r, order):
    # Direct Gram determinant of the two lifted orthonormal angular vectors.
    # Test the interior cap and the entire exterior Jacobi map on many bearings.
    c = sqrt(3)/2
    u, w = leggauss(order)
    mu = c+(1+u)*(1-c)/2
    a = 2*pi*np.arange(2*order)/(2*order)
    mu, a = np.broadcast_arrays(mu[:, None], a[None, :])
    mu, a = mu.ravel(), a.ravel()
    h = np.sqrt(1-mu*mu)
    n = np.column_stack((mu, h*np.cos(a), h*np.sin(a)))
    E = np.stack((np.column_stack((-h, mu*np.cos(a), mu*np.sin(a))),
                  np.column_stack((0*a, -np.sin(a), np.cos(a)))), axis=1)
    X = np.cross(n[:, None, :], E)
    t = (sin(r/ROUND_RADIUS)/(r/ROUND_RADIUS))**2
    if r <= .02:
        cutoff = 1
    elif r >= .03:
        cutoff = 0
    else:
        raise ValueError("transition band is outside this test")
    plus, minus = (t*np.eye(2) + cutoff*r*r*(X*np.array(B)) @ X.transpose(0, 2, 1)
                   for B in ((3, 1, 2), (1, 3, 2)))
    if not cutoff:
        np.testing.assert_array_equal(plus, minus)
    weights = np.repeat(w, 2*order)*(1-c)*pi/(2*order)
    # Rationalize the difference by subtracting its exact determinant formula.
    # Compare as well with the two independently evaluated direct Gram matrices.
    denominator = np.sqrt(np.linalg.det(plus))+np.sqrt(np.linalg.det(minus))
    exact = -2*cutoff*(t+2*cutoff*r*r)*(n[:, 0]**2-n[:, 1]**2)/denominator
    direct = (np.sqrt(np.linalg.det(plus))-np.sqrt(np.linalg.det(minus)))/r**2
    np.testing.assert_allclose(exact, direct, rtol=0, atol=1e-9)
    return weights @ exact


def test_separating_interior_cap_and_exterior_geometry():
    angular_factor = pi*sqrt(3)/8
    for r in (.02, 1/128, 1/256, 1/1024):
        val = pixels(r, 12)
        assert -angular_factor*(1+r*r/.97) < val < -angular_factor*(1-r*r/.97)
        assert abs(val-pixels(r, 18)) < 1e-13
    assert pixels(.031, 12) == 0 and pixels(.04, 12) == 0
    survey = json.loads((m.HERE.parent/"results.json").read_text())
    axial = survey["tables"]["axial"]
    assert all(abs(a["signed_half_difference_per_s_ell4"]) > 64*pi/45 for a in axial)
    assert abs(axial[-1]["signed_half_difference_per_s_ell4"]+64*pi/45) < .0012
