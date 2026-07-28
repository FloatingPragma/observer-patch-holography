"""Exact-case tests for the ball-arithmetic loop layer and contours."""

from __future__ import annotations

import sys
from pathlib import Path

import mpmath as mp
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "producers"))

import ball_arithmetic as ba  # noqa: E402


def test_a0_exact_value_and_enclosure() -> None:
    ball = ba.a0_fin(4, 1, precision=192)
    with mp.workprec(192):
        exact = mp.mpf(4) * (1 - mp.log(4))
        assert ball.contains(exact, mp.mpf(0))
        assert ball.width_bound() < mp.mpf(2) ** -150


def test_a0_scaleless_zero() -> None:
    ball = ba.a0_fin(0, 1)
    assert ball.mid_re == 0 and ball.rad == 0


def test_b0_zero_momentum_equal_masses() -> None:
    ball = ba.b0_fin(0, 9, 9, 1, precision=192)
    with mp.workprec(192):
        exact = -mp.log(mp.mpf(9))
        assert ball.contains(exact, mp.mpf(0))


def test_b0_below_threshold_real() -> None:
    ball = ba.b0_fin(1, 4, 4, 1, precision=192)
    assert abs(ball.mid_im) < mp.mpf(2) ** -60


def test_b0_absorptive_part_above_threshold() -> None:
    # equal masses m^2 = 1, p2 = 8 > 4: Im B0 = pi sqrt(1 - 4/p2)
    ball = ba.b0_fin(8, 1, 1, 1, precision=192)
    with mp.workprec(192):
        beta = mp.sqrt(1 - mp.mpf(4) / 8)
        assert abs(ball.mid_im - mp.pi * beta) < mp.mpf(10) ** -30


def test_b0_precision_presets_only() -> None:
    with pytest.raises(ValueError):
        ba.b0_fin(1, 1, 1, 1, precision=64)


def test_winding_counts_zeros() -> None:
    def f(z: complex) -> complex:
        return (z - (1 + 1j)) * (z + 2)

    assert ba.certify_winding(f, (0.5 + 0.5j, 1.5 + 1.5j), subdivisions=256) == 1
    assert ba.certify_winding(f, (-3 - 1j, 2 + 2j), subdivisions=512) == 2
    assert ba.certify_winding(f, (3 + 3j, 4 + 4j), subdivisions=128) == 0


def test_winding_refuses_boundary_zero() -> None:
    def f(z: complex) -> complex:
        return z - 1

    with pytest.raises(RuntimeError):
        ba.certify_winding(f, (1 - 1j, 3 + 1j), subdivisions=8)
