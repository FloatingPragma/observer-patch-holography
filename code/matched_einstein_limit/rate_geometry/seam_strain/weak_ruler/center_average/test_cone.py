"""Direct wedge integrals, source-angle bands and volume-matching checks."""

import json
import math
from pathlib import Path

import pytest
from scipy.integrate import quad

from . import cone


def direct_fractions(alpha, u):
    eta = math.pi - alpha / 2
    y = u * math.sin(eta)
    start = math.asin(y)

    def missing(power):
        def ray(theta):
            return max(0, 1 - (y / math.sin(theta))**2)**(power / 2)

        return quad(ray, start, eta, epsabs=1e-13, epsrel=1e-12)[0] / math.pi

    return 1 - missing(3), 1 - missing(1)


def direct_individual(alpha):
    def integrand(u):
        b, p = direct_fractions(alpha, u)
        return u * p * b**(-5/3) * (1 - p * b**(-2/3))

    return 4 * math.pi * alpha * quad(integrand, 0, 1, epsabs=1e-10)[0]


def test_area_volume_profiles_and_angular_derivatives():
    alpha = 3.66
    step = 0.0001
    profile = cone._profile(cone._half_gap(alpha), derivatives=True)
    for u in (0.03125, 0.1875, 0.5, 0.875, 0.96875):
        b, p, db, dp = profile(cone._iv.mpf(u))
        assert (float(b.mid), float(p.mid)) == pytest.approx(
            direct_fractions(alpha, u), rel=0, abs=1e-12
        )
        plus, minus = direct_fractions(alpha + step, u), direct_fractions(alpha - step, u)
        numeric = [(x - y) / (2*step) for x, y in zip(plus, minus, strict=True)]
        assert (float(db.mid), float(dp.mid)) == pytest.approx(numeric, rel=0, abs=2e-10)


def test_common_volume_formula_and_weak_deficit_limit():
    alpha = 3.66
    integral = quad(lambda u: u * (1 - direct_fractions(alpha, u)[0]), 0, 1,
                    epsabs=1e-13, epsrel=1e-12)[0]
    assert cone.mean_volume_coefficient(alpha) == pytest.approx(
        alpha * 4*math.pi/3 * integral, abs=1e-12
    )
    delta = 1e-4
    linear = 4 * math.pi * delta / 15
    assert abs(3 * cone.mean_volume_coefficient(math.tau - delta) - linear) < delta**2
    assert abs(direct_individual(math.tau - delta) - linear) < delta**2


def test_each_observer_fixes_the_ball_volume():
    alpha = math.tau - 2.6225
    for offset in (0, 0.5, 0.95, 1, 2):
        box = cone.individual_ball(alpha, offset)
        radius = sum(box["radius"]) / 2
        assert box["volume"][0] <= 4*math.pi/3 <= box["volume"][1]
        if offset == 0:
            b = p = alpha / math.tau
        elif offset >= 1:
            b = p = 1
        else:
            b, p = direct_fractions(alpha, offset / radius)
        assert radius**3 * b == pytest.approx(1, abs=2e-13)
        assert 4*math.pi*radius**2*p == pytest.approx(sum(box["area"])/2, abs=2e-12)
        assert box["area"][1] - box["area"][0] < 2e-13
    assert 11.0137 < cone.individual_ball(alpha, 0.5)["area"][0] < 11.0138


def test_native_deficit_bands_distinguish_both_protocols():
    archive = Path(__file__).parents[1] / "hinge.json"
    native = json.loads(archive.read_text())["panels"][0]
    assert native["divisor"] == 16 and native["continuous"] is False
    bands = [(2.575, 2.580, 1.7648, 1.8152), (2.620, 2.625, 1.7956, 1.8469)]
    for delta in native["interior_deficits"]:
        assert any(low < delta < high for low, high, _, _ in bands)
    certified = json.loads(Path(__file__).with_name("native.json").read_text())["radial_edges"]
    assert len(certified) == 4
    for edge in certified:
        lower_angle, upper_angle = edge["angle"]
        assert any(math.tau-high < lower_angle <= upper_angle < math.tau-low
                   for low, high, _, _ in bands)
    for low, high, target_lo, target_hi in bands:
        alpha_lo = math.nextafter(math.tau - high, -math.inf)
        alpha_hi = math.nextafter(math.tau - low, math.inf)
        f_lo, f_hi = cone.mean_volume_interval(alpha_lo, alpha_hi)
        i_lo, i_hi = cone.individual_interval(alpha_lo, alpha_hi)
        assert 0 < 3*f_lo < 3*f_hi < target_lo < i_lo < i_hi < target_hi
        assert target_hi < 4*math.pi*low/15
        assert i_lo < direct_individual((alpha_lo + alpha_hi)/2) < i_hi


def test_individual_coefficient_derivative_encloses_independent_held_secant():
    alpha, step = math.tau - 2.6215998703, 1e-4
    lo, hi = cone.individual_slope_interval(alpha - step, alpha + step, bins=1024)
    difference = (direct_individual(alpha + step) - direct_individual(alpha - step)) / (2*step)
    assert lo < difference < hi < 0
    assert hi - lo < 0.07
