#!/usr/bin/env python3
"""Strict finite-order W/Z complex-pole algebra.

The convention is Gamma(s)=s-m0^2+kappa*Pi(s), with
kappa=(16*pi^2)^-1.  The module computes conditional coefficients only.  It
does not supply an OPH source law, an EFT matching receipt, or gauge/BRST
verification.
"""

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass


KAPPA = 1.0 / (16.0 * math.pi**2)


class StrictSeriesError(ValueError):
    """Raised when a finite-order pole packet is inconsistent."""


@dataclass(frozen=True)
class TreeMassSquares:
    w: float
    z: float


@dataclass(frozen=True)
class StrictOneLoopPole:
    tree_mass_square: float
    coefficient_s1: complex
    truncated_s: complex
    strict_mass_shift: float
    strict_width: float


def _positive(name: str, value: float) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise StrictSeriesError(f"{name} must be positive and finite")
    return value


def tree_mass_squares(g: float, gprime: float, v_f: float) -> TreeMassSquares:
    g = _positive("g", g)
    gprime = _positive("gprime", gprime)
    v_f = _positive("v_f", v_f)
    w = g * g * v_f * v_f / 4.0
    z = (g * g + gprime * gprime) * v_f * v_f / 4.0
    if not z > w:
        raise StrictSeriesError("tree masses must obey z>w>0")
    return TreeMassSquares(w=w, z=z)


def strict_one_loop_pole(tree_mass_square: float, pi_at_tree: complex) -> StrictOneLoopPole:
    tree_mass_square = _positive("tree_mass_square", tree_mass_square)
    pi_at_tree = complex(pi_at_tree)
    if not all(math.isfinite(value) for value in (pi_at_tree.real, pi_at_tree.imag)):
        raise StrictSeriesError("self energy must be finite")
    s1 = -pi_at_tree
    truncated = complex(tree_mass_square, 0.0) + KAPPA * s1
    m0 = math.sqrt(tree_mass_square)
    mass_shift = KAPPA * s1.real / (2.0 * m0)
    width = -KAPPA * s1.imag / m0
    if width < -1e-14:
        raise StrictSeriesError("pole is on the growing-state sheet")
    return StrictOneLoopPole(
        tree_mass_square=tree_mass_square,
        coefficient_s1=s1,
        truncated_s=truncated,
        strict_mass_shift=mass_shift,
        strict_width=max(0.0, width),
    )


def charged_two_loop_coefficient(
    *, pi2_at_w: complex, s1: complex, pi1_prime_at_w: complex
) -> complex:
    """Return s_W,2=-Pi_WW,2-s_W,1 Pi'_WW,1."""

    return -complex(pi2_at_w) - complex(s1) * complex(pi1_prime_at_w)


def neutral_two_loop_coefficient(
    *,
    z: float,
    pi2_zz_at_z: complex,
    s1: complex,
    pi1_zz_prime_at_z: complex,
    pi1_za_at_z: complex,
    pi1_az_at_z: complex,
) -> complex:
    """Return the strict neutral two-loop coefficient, including ZA mixing."""

    z = _positive("z", z)
    return (
        -complex(pi2_zz_at_z)
        - complex(s1) * complex(pi1_zz_prime_at_z)
        + complex(pi1_za_at_z) * complex(pi1_az_at_z) / z
    )


def display_energy_coordinates(s_pole: complex) -> tuple[float, float]:
    """Return M,Gamma from an exact square root of a truncated s pole.

    This is a display coordinate.  It is not the strict finite-order mass and
    width coefficient because the square root resums kinematic powers.
    """

    root = cmath.sqrt(complex(s_pole))
    if root.real < 0 or (root.real == 0 and root.imag > 0):
        root = -root
    if root.imag > 1e-14:
        raise StrictSeriesError("pole is on the growing-state sheet")
    return float(root.real), float(max(0.0, -2.0 * root.imag))
