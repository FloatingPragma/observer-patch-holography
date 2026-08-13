from fractions import Fraction
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from modal_maxwell_factorization import (  # noqa: E402
    GZERO,
    VZERO,
    exact_replay,
    fourier_curl,
    maxwell_shaped_generator,
    momentum_dot,
    same_sign_mutation,
    second_order_residual,
)


def off_axis_witness():
    momentum = (Fraction(3), Fraction(4), Fraction(0))
    state = (
        (
            (Fraction(4), Fraction(0)),
            (Fraction(-3), Fraction(0)),
            (Fraction(2), Fraction(1)),
        ),
        (
            (Fraction(0), Fraction(4)),
            (Fraction(0), Fraction(-3)),
            (Fraction(5), Fraction(0)),
        ),
    )
    return momentum, Fraction(7), Fraction(7, 5), state


def test_exact_replay_carries_both_mutation_controls() -> None:
    replay = exact_replay()
    assert replay["arithmetic"] == "exact_gaussian_fraction"
    assert replay["same_sign_mutation_rejected"] is True
    assert replay["same_sign_residual_exactly_twice_spatial_action"] is True
    assert replay["wrong_scale_mutation_rejected"] is True
    assert replay["computes_fz12_symbol"] is False
    assert replay["physical_claim"] is False


def test_off_axis_div_curl_is_exact() -> None:
    momentum, _, curl_scale, state = off_axis_witness()
    for amplitude in state:
        assert momentum_dot(momentum, amplitude) == GZERO
        assert momentum_dot(
            momentum, fourier_curl(curl_scale, momentum, amplitude)
        ) == GZERO


def test_opposite_sign_maxwell_shape_obeys_wave_equation() -> None:
    momentum, frequency, curl_scale, state = off_axis_witness()
    assert second_order_residual(
        maxwell_shaped_generator, curl_scale, momentum, frequency, state
    ) == (VZERO, VZERO)


def test_same_sign_mutation_has_wrong_second_order_sign() -> None:
    momentum, frequency, curl_scale, state = off_axis_witness()
    assert second_order_residual(
        same_sign_mutation, curl_scale, momentum, frequency, state
    ) != (VZERO, VZERO)


def test_wrong_normalization_fails_the_wave_equation() -> None:
    momentum, frequency, _, state = off_axis_witness()
    assert second_order_residual(
        maxwell_shaped_generator, Fraction(2), momentum, frequency, state
    ) != (VZERO, VZERO)


def test_zero_frequency_zero_curl_boundary() -> None:
    momentum = (Fraction(3), Fraction(4), Fraction(0))
    state = (
        (
            (Fraction(4), Fraction(0)),
            (Fraction(-3), Fraction(0)),
            (Fraction(0), Fraction(0)),
        ),
        (
            (Fraction(0), Fraction(0)),
            (Fraction(0), Fraction(0)),
            (Fraction(1), Fraction(0)),
        ),
    )
    assert second_order_residual(
        maxwell_shaped_generator, Fraction(0), momentum, Fraction(0), state
    ) == (VZERO, VZERO)
