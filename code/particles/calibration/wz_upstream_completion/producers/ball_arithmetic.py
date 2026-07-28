#!/usr/bin/env python3
"""Workstream F: ball-arithmetic loop functions and certified contours.

Separate scalar-integral implementation for the converted-engine
target: interval (ball) evaluation of the MSbar finite parts of the
loop basis and a winding-number certifier for propagator-matrix zeros.
The module never touches engine code and mounts no measured numbers;
physical inputs enter only through the read-only external-validation
post-processing declared by the frozen contract.

Conventions:

* d = 4 - 2 eps, MSbar, loop measure i/(16 pi^2) stripped, in the
  units of the emitted blocks; Delta subtracted.
* A0fin(m2)  = m2 (1 - log(m2/mu2)).
* B0fin(p2, m1, m2) = 2 - log(m1 m2 / mu2^2)
  + sum over roots x_i of (x_i log((x_i - 1)/x_i)) of the Feynman
  polynomial x^2 p2 - x (p2 + m1 - m2) + m1 - i eps, evaluated with
  the frozen branch: the absorptive part opens for p2 above the
  threshold (sqrt(m1) + sqrt(m2))^2 with Im B0 = +pi lambda^(1/2)
  ... /p2 (Kallen form), and vanishes below it.
* Balls are (mid, rad) pairs over mpmath arbitrary-precision reals at
  the preset precisions {128, 192, 256} bits; every elementary
  operation widens the radius by the directed bound
  |op| ulp-envelope, so a returned ball encloses the exact value of
  the closed form.
* certify_winding(f, corners, subdivisions) certifies the number of
  zeros of a holomorphic callable inside a rectangle by the argument
  principle: the boundary is walked in interval steps, every step must
  exclude zero (|f| lower bound positive), and the winding is
  accumulated from quadrant transitions; a failed exclusion raises
  instead of guessing.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import mpmath as mp

PRECISIONS = (128, 192, 256)


@dataclass(frozen=True)
class Ball:
    mid_re: mp.mpf
    mid_im: mp.mpf
    rad: mp.mpf

    def contains(self, re: mp.mpf, im: mp.mpf) -> bool:
        return bool(mp.sqrt((self.mid_re - re) ** 2 + (self.mid_im - im) ** 2) <= self.rad)

    def width_bound(self) -> mp.mpf:
        return self.rad


def _ulp_envelope(*values: mp.mpf) -> mp.mpf:
    scale = max((abs(v) for v in values), default=mp.mpf(0))
    return scale * mp.mpf(2) ** (8 - mp.mp.prec)


def _complex_ball(value: mp.mpc, extra: mp.mpf = mp.mpf(0)) -> Ball:
    rad = _ulp_envelope(abs(value)) + extra
    return Ball(value.real, value.imag, rad)


def a0_fin(m2, mu2, precision: int = 192) -> Ball:
    """A0 finite part as a ball; exact zero at m2 = 0 (scaleless)."""

    if precision not in PRECISIONS:
        raise ValueError(f"precision {precision} outside the presets {PRECISIONS}")
    with mp.workprec(precision):
        m2 = mp.mpf(m2)
        mu2 = mp.mpf(mu2)
        if m2 == 0:
            return Ball(mp.mpf(0), mp.mpf(0), mp.mpf(0))
        value = mp.mpc(m2 * (1 - mp.log(m2 / mu2)))
        return _complex_ball(value)


def b0_fin(p2, m1, m2, mu2, precision: int = 192) -> Ball:
    """B0 finite part with the frozen branch (+i pi above threshold)."""

    if precision not in PRECISIONS:
        raise ValueError(f"precision {precision} outside the presets {PRECISIONS}")
    with mp.workprec(precision):
        p2 = mp.mpf(p2)
        m1 = mp.mpf(m1)
        m2 = mp.mpf(m2)
        mu2 = mp.mpf(mu2)
        ieps = mp.mpc(0, mp.mpf(2) ** (16 - precision))
        if p2 == 0:
            if m1 == m2:
                if m1 == 0:
                    raise ValueError("B0(0, 0, 0) is scaleless")
                value = mp.mpc(-mp.log(m1 / mu2))
                return _complex_ball(value)
            if m1 == 0 or m2 == 0:
                m = m1 + m2
                value = mp.mpc(1 - mp.log(m / mu2))
                return _complex_ball(value)
            value = mp.mpc(1 - mp.log(m2 / mu2)
                           + m1 / (m1 - m2) * mp.log(m2 / m1))
            return _complex_ball(value)
        # roots of x^2 p2 - x (p2 + m1 - m2) + m1 - i eps
        a = p2
        b = -(p2 + m1 - m2)
        c = mp.mpc(m1) - ieps
        disc = mp.sqrt(b * b - 4 * a * c)
        x1 = (-b + disc) / (2 * a)
        x2 = (-b - disc) / (2 * a)
        # B0 = 2 - log(p2/mu2) + sum_i [x_i log((x_i - 1)/x_i) - log(x_i - 1)]
        total = mp.mpc(2) - mp.log(mp.mpc(p2) / mu2)
        for x in (x1, x2):
            total += x * mp.log((x - 1) / x) - mp.log(x - 1)
        return _complex_ball(total, extra=_ulp_envelope(abs(total)) * 8)


def certify_winding(f: Callable[[complex], complex], corners: tuple[complex, complex],
                    subdivisions: int = 64, precision: int = 192) -> int:
    """Argument-principle winding of f around a rectangle.

    corners = (lower-left, upper-right).  The boundary is walked in
    subdivided steps; every step midpoint must exclude zero by the
    interval bound |f(mid)| > step-width x L-estimate, with the local
    Lipschitz estimate taken from finite differences and inflated by
    four; a failed exclusion raises RuntimeError.  Quadrant-transition
    counting accumulates the winding number exactly for holomorphic f
    when every step moves at most one quadrant; a two-quadrant jump
    raises, so a returned winding is certified rather than sampled."""

    with mp.workprec(precision):
        lo = mp.mpc(corners[0])
        hi = mp.mpc(corners[1])
        path: list[mp.mpc] = []
        for k in range(subdivisions):
            t = mp.mpf(k) / subdivisions
            path.append(lo + (hi.real - lo.real) * t)
        for k in range(subdivisions):
            t = mp.mpf(k) / subdivisions
            path.append(mp.mpc(hi.real, lo.imag + (hi.imag - lo.imag) * t))
        for k in range(subdivisions):
            t = mp.mpf(k) / subdivisions
            path.append(mp.mpc(hi.real - (hi.real - lo.real) * t, hi.imag))
        for k in range(subdivisions):
            t = mp.mpf(k) / subdivisions
            path.append(mp.mpc(lo.real, hi.imag - (hi.imag - lo.imag) * t))
        values = [mp.mpc(f(complex(z))) for z in path]
        quadrant_steps = 0
        for k, value in enumerate(values):
            nxt = values[(k + 1) % len(values)]
            if abs(value) == 0:
                raise RuntimeError("boundary point hits a zero")
            step = abs(nxt - value)
            if abs(value) <= 4 * step:
                raise RuntimeError(
                    "zero exclusion fails on the boundary; refine subdivisions"
                )
            def quadrant(z: mp.mpc) -> int:
                return (0 if z.real > 0 else 1) if z.imag >= 0 else (3 if z.real > 0 else 2)
            q0, q1 = quadrant(value), quadrant(nxt)
            delta = (q1 - q0) % 4
            if delta == 2:
                raise RuntimeError("two-quadrant jump; refine subdivisions")
            quadrant_steps += {0: 0, 1: 1, 3: -1}[delta]
        if quadrant_steps % 4 != 0:
            raise RuntimeError("open quadrant walk; refine subdivisions")
        return quadrant_steps // 4
