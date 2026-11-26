from fractions import Fraction as F
import json
import math

import pytest

from .audit import (COEFFICIENTS, HERE, RATIOS, SEPARATION, annular_kernel_over_pi,
                    finite_bands, kernel_over_pi, native_transfer, normalized_lower_bound, radius)


@pytest.mark.parametrize("low,side", [(0, "right"), (0, "left"), (1, "left"), (1, "right")])
def test_native_word_transfer(low, side):
    row = native_transfer(low, side)
    assert row["energy_laws_in_sixteenths"] == [[4, 6, 6]]*2
    assert row["transfer_load"] == sum(row["eigenwords"][1]) == 2*low+1
    assert sorted(row["eigenwords"][0]) == sorted(row["eigenwords"][1])
    assert row["undo_laws_in_sixteenths"] == [[10, 2, 2, 2], [9, 2, 2, 3]]
    assert abs(row["entropy_slope"]-math.log(5)) < 1e-12


def test_signed_interior_and_exterior_cancel_and_annular_response_has_fixed_sign():
    r = F(1, 512)
    atoms = [(F(0), -7), (r/3, 9), (r*15/16, -3), (16*r, -49), (19*r, 80)]
    b = [sum(weight*kernel_over_pi(q*r, x) for x, weight in atoms) for q in RATIOS]
    assert sum(c*value for c, value in zip(COEFFICIENTS, b, strict=True)) == 0
    for q, expected in [(F(0), 0), (F(1), 0), (F(2), 12), (F(4), 60),
                        (F(8), 48), (F(16), 0), (F(18), 0)]:
        assert annular_kernel_over_pi(r, q*r) == expected*r
    # Uniform density on all balls acts in the annuli: its diamond integral scales as r**4.
    assert sum(c*q**4 for c, q in zip(COEFFICIENTS, RATIOS, strict=True)) == 257796


def test_best_signed_two_moment_fit_on_one_band():
    c = normalized_lower_bound()
    assert c == F(9, 266500)
    # Normalize R and the common slope to one. Fit the optimum using two moments.
    y1, y4 = 1+c, 1-4**4*c
    e = (4*y4-y1)/15
    m = e-y1
    errors = [(1-(e*q-m/q))/q**4 for q in RATIOS]
    assert errors == [-c, c, -c]


def test_prefixes_use_exactly_the_same_sources_and_the_infinite_tail_is_bounded():
    short, long = finite_bands(3), finite_bands(5)
    for k in range(3):
        assert short[k]["radius"] == long[k]["radius"] == str(radius(k))
        assert short[k]["contained_active_pairs"] == list(range(k, 3))
        assert long[k]["contained_active_pairs"] == list(range(k, 5))
        extra = sum((radius(i)**4 for i in (3,4)), F())
        entire_tail = radius(k)**4/(1-F(1, SEPARATION)**4)
        short_weight = sum((radius(i)**4 for i in range(k, 3)), F())
        assert radius(k)**4 <= short_weight < short_weight+extra < entire_tail
    assert short[-1]["minimum_normalized_error"] == float(normalized_lower_bound())*math.log(5)
    assert long[-1]["minimum_normalized_error"] == short[-1]["minimum_normalized_error"]


def test_frozen_prefix_and_bank_records():
    archive = json.loads((HERE/"results.json").read_text())
    assert archive["finite_prefix"] == finite_bands(4)
    assert archive["native_dimension"] == 99
    for row in archive["native_transfers"]:
        assert native_transfer(row["low_load"], row["reversed_port"]) == row


def test_regular_density_has_one_limit_and_the_recorded_entropy_has_two():
    # Radially integrate a density a + b*r^2 with the flat diamond kernel.
    # The common geometric factor 4*pi^2 has been divided out.
    a, b = F(7,2), F(-9, 5)
    expected_limit = a*(F(1,3)-F(1,5))
    for q in RATIOS:
        t = q*radius(2)
        integral = a*(t**4/3 - t**4/5) + b*(t**6/5 - t**6/7)
        assert (integral/t**4-expected_limit)/t**2 == b*(F(1,5)-F(1,7))
    archive = json.loads((HERE/"results.json").read_text())
    large, _, small = archive["entropy_slopes_per_radius_fourth_infinite"]
    assert large == 16**4 * small
    assert abs((large-small)/2-archive["minimum_limsup_error_for_regular_stress"]) < 1e-15
