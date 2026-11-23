from fractions import Fraction
from math import pi, sqrt

import numpy as np
import pytest

from . import corridor_release, open_pair_release


def test_connected_receiver_counts_the_entire_preparation():
    rho, record = open_pair_release.source()
    wave = np.array([-13, 0, 14])/27
    f = Fraction(364, 729)
    reaches = [f**i for i in range(8)]
    floor, exhaustion = open_pair_release.source_bound()
    assert exhaustion == f**8 and floor == (1-f**8)*729/730
    assert Fraction(9947, 10000) < floor < Fraction(995, 1000)
    np.testing.assert_allclose(record["first_vacancy_in_one_attempt"], [4/9, 4/81, 4/729], atol=3e-15, rtol=0)
    np.testing.assert_allclose(record["attempt_reach"], np.array(reaches, float), atol=3e-15, rtol=0)
    np.testing.assert_allclose(rho, float(sum(reaches))*np.outer(wave, wave), rtol=0, atol=3e-15)
    arms = open_pair_release.terminal(rho, record["exhaustion"], 0., 1.)
    for bins in arms.values():
        assert len(bins) == 31 and min(bins) >= 0 and abs(sum(bins)-1) < 2e-11
    assert float(exhaustion) <= arms["coherent"][-1] <= 1-float(floor)
    assert arms["stochastic_swap"][-1] > .20 and arms["dephased"][-1] > .10
    with pytest.raises(ValueError):
        open_pair_release.terminal(rho, 0., 0., 1.)


@pytest.mark.parametrize("budgets", [(0, 8), (3, -1), (True, 8), (3.1, 8)])
def test_source_budget_rejects_invalid_counts(budgets):
    with pytest.raises(ValueError):
        open_pair_release.source_bound(*budgets)


@pytest.mark.parametrize("rounds,attempts", [(1, 1), (3, 8), (12, 23)])
def test_larger_finite_sources_can_meet_arbitrary_tolerance(rounds, attempts):
    floor, exhausted = open_pair_release.source_bound(rounds, attempts)
    assert 1-floor <= Fraction(1, 9)**rounds + Fraction(1, 2)**attempts
    assert 0 < exhausted < 1-floor


@pytest.mark.parametrize("strength", [-8, -.25, 0, .5, 8, sqrt((2*pi)**2-.5)-.25])
def test_closed_control_including_pulse_alias(strength):
    source = corridor_release.acquire(strength)
    rho, exhausted, reaches, first = source
    assert len(first) == 128 and len(reaches) == 8
    assert abs(sum(first)*reaches[-1]-exhausted) < 1e-13 and abs(rho.trace()+exhausted-1) < 5e-12
    assert corridor_release.uniform_floor() == Fraction(2975, 3072)
    for t in (0., 1., 37.):
        arms = corridor_release.terminal(strength, t, source)
        for bins in arms.values():
            assert len(bins) == 31 and min(bins) >= 0 and abs(sum(bins)-1) < 5e-12
        assert 1-arms["coherent"][-1] >= float(corridor_release.uniform_floor())
        assert abs(1-arms["coherent"][-1]-corridor_release.analytic_home(strength, t)) < 5e-12
