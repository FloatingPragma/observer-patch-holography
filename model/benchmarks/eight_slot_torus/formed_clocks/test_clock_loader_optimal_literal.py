"""Prospective full-coordinate audit of a 26-service clock preparation."""
import math
from fractions import Fraction

import numpy as np
import pytest

from . import clock_loader_optimal
from .test_clock_loader_literal import PANELS, _loop, _run, _shift


def test_twenty_five_native_hops_have_insufficient_angular_budget():
    # Sum theta/pi <= (22/7)^2 sum(k^2)/(8*128^2), using 1-cos(x)<=x^2/2.
    upper = Fraction(22, 7)**2 * sum(k*k for k in range(1, 26)) / (8 * 128**2)
    assert upper == Fraction(668525, 1605632) < Fraction(5, 12)


def test_steering_angles_have_an_exact_root_inside_the_enclosure():
    assert clock_loader_optimal.certify_angles()['fixed_point_interior']


@pytest.mark.parametrize('shape,root', PANELS)
def test_literal_two_tag_26_service_channel_and_local_detector(shape, root):
    gates = clock_loader_optimal.schedule(shape, root)
    assert set(gates) <= set(range(1, 27)) and 26 in gates
    loaded = _run(shape, root, gates, 26)
    mu = np.array([math.pi/32, math.pi/32 + math.pi/4096])
    integral = sum(math.pi*(1-math.cos(k*math.pi/128))/4 for k in range(1, 27))
    loop = _loop(shape, root)
    # Both full rail columns and their complex phases fix the action on any reference.
    expected = np.exp(-1j*(integral+mu))[:, None] * loop[None, :] / math.sqrt(2)
    np.testing.assert_allclose(loaded, expected, atol=2e-12, rtol=0)
    ends = (root, _shift(shape, _shift(shape, root, 0, 1), 1, 1))
    for age in (26, 133):
        pair = loaded if age == 26 else _run(shape, root, gates, age)
        # Fresh terminal cohort, coherent input |+>, original birth probability 1/2.
        for site in ends:
            i = 8*np.ravel_multi_index(site, shape)
            a, b = pair[:, i:i+8]/math.sqrt(2)
            words = [np.vdot(a+h*b, a+h*b).real/4 for h in (1, 1j, -1, -1j)]
            assert abs(sum(words) - 1/8) < 2e-11
            hand = complex(words[0]-words[2], words[1]-words[3])
            predicted = np.exp(1j*(1-(age-26))*(mu[1]-mu[0]))/16
            assert abs(hand - predicted) < 2e-11
