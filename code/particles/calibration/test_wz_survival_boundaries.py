from __future__ import annotations

import math
from fractions import Fraction

import pytest

from source_completion_no_go import CompletionWitnessError, two_completion_separation
from wz_strict_series import (
    KAPPA,
    StrictSeriesError,
    charged_two_loop_coefficient,
    display_energy_coordinates,
    neutral_two_loop_coefficient,
    strict_one_loop_pole,
    tree_mass_squares,
)


def test_tree_mass_map_and_one_loop_sign() -> None:
    masses = tree_mass_squares(0.65, 0.35, 246.0)
    assert masses.z > masses.w > 0
    pole = strict_one_loop_pole(masses.w, complex(-100.0, 50.0))
    assert pole.coefficient_s1 == complex(100.0, -50.0)
    assert pole.truncated_s == complex(masses.w, 0.0) + KAPPA * pole.coefficient_s1
    assert pole.strict_width > 0


def test_neutral_mixing_is_two_loop_only() -> None:
    one_loop = strict_one_loop_pole(100.0, complex(4.0, 2.0))
    unchanged = strict_one_loop_pole(100.0, complex(4.0, 2.0))
    assert one_loop == unchanged
    s2 = neutral_two_loop_coefficient(
        z=100.0,
        pi2_zz_at_z=1.0,
        s1=one_loop.coefficient_s1,
        pi1_zz_prime_at_z=2.0,
        pi1_za_at_z=3.0,
        pi1_az_at_z=5.0,
    )
    without_mixing = -1.0 - one_loop.coefficient_s1 * 2.0
    assert s2 - without_mixing == pytest.approx(15.0 / 100.0)


def test_charged_two_loop_derivative_term() -> None:
    assert charged_two_loop_coefficient(pi2_at_w=7, s1=3, pi1_prime_at_w=2) == -13


def test_strict_and_display_coordinates_are_distinct() -> None:
    pole = strict_one_loop_pole(6400.0, complex(-150.0, 80.0))
    mass, width = display_energy_coordinates(pole.truncated_s)
    strict_mass = math.sqrt(6400.0) + pole.strict_mass_shift
    assert abs(mass - strict_mass) > 0
    assert abs(width - pole.strict_width) > 0


def test_growing_sheet_is_rejected() -> None:
    with pytest.raises(StrictSeriesError):
        strict_one_loop_pole(100.0, complex(1.0, -1.0))


def test_two_completion_exact_separation() -> None:
    result = two_completion_separation(Fraction(3, 2))
    assert result.higgs_separation == 2 * result.u * (1 - result.u)
    assert result.top_separation == result.u * (1 - result.u) / 2


@pytest.mark.parametrize("pixel", [Fraction(0), Fraction(24), Fraction(25)])
def test_two_completion_domain_guard(pixel: Fraction) -> None:
    with pytest.raises(CompletionWitnessError):
        two_completion_separation(pixel)
