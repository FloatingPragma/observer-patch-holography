#!/usr/bin/env python3
"""Workstream F: high-precision loop functions and sampled contours.

Separate scalar-integral implementation for the converted-engine
target: arbitrary-precision evaluation of the MSbar finite parts of the
loop basis and a sampled winding diagnostic for propagator-matrix zeros.
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
* ``DiagnosticEstimate`` values carry a midpoint and a roundoff
  heuristic.  mpmath does not provide directed rounding, and the
  radius does not enclose removal of the finite ``-i eps`` regulator.
  These objects are therefore diagnostics, not proof-bearing complex
  balls.
* ``sample_winding`` walks point samples around a rectangle and
  accumulates quadrant transitions.  Its finite-difference guard is a
  refinement heuristic, not an analytic derivative bound.  A returned
  value is not an argument-principle certificate.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import mpmath as mp

PRECISIONS = (128, 192, 256)


@dataclass(frozen=True)
class DiagnosticEstimate:
    mid_re: mp.mpf
    mid_im: mp.mpf
    rad: mp.mpf

    def within_reported_radius(self, re: mp.mpf, im: mp.mpf) -> bool:
        return bool(mp.sqrt((self.mid_re - re) ** 2 + (self.mid_im - im) ** 2) <= self.rad)

    def reported_radius(self) -> mp.mpf:
        return self.rad


def _ulp_envelope(*values: mp.mpf) -> mp.mpf:
    scale = max((abs(v) for v in values), default=mp.mpf(0))
    return scale * mp.mpf(2) ** (8 - mp.mp.prec)


def _complex_estimate(
    value: mp.mpc,
    extra: mp.mpf = mp.mpf(0),
) -> DiagnosticEstimate:
    rad = _ulp_envelope(abs(value)) + extra
    return DiagnosticEstimate(value.real, value.imag, rad)


def a0_fin(m2, mu2, precision: int = 192) -> DiagnosticEstimate:
    """A0 finite-part estimate; exact zero at m2 = 0 (scaleless)."""

    if precision not in PRECISIONS:
        raise ValueError(f"precision {precision} outside the presets {PRECISIONS}")
    with mp.workprec(precision):
        m2 = mp.mpf(m2)
        mu2 = mp.mpf(mu2)
        if m2 == 0:
            return DiagnosticEstimate(mp.mpf(0), mp.mpf(0), mp.mpf(0))
        value = mp.mpc(m2 * (1 - mp.log(m2 / mu2)))
        return _complex_estimate(value)


def b0_fin(p2, m1, m2, mu2, precision: int = 192) -> DiagnosticEstimate:
    """B0 diagnostic with a finite regulator (+i pi above threshold).

    Accepts complex p2: the root formula continues analytically off the
    real axis using the same precision-dependent ``-i eps`` prescription.
    The returned radius covers the evaluated closed form only; it does not
    certify the ``eps -> 0`` limit or a global continuation sheet."""

    if precision not in PRECISIONS:
        raise ValueError(f"precision {precision} outside the presets {PRECISIONS}")
    with mp.workprec(precision):
        if isinstance(p2, complex) or (hasattr(p2, "imag") and getattr(p2, "imag", 0) != 0):
            return _b0_fin_complex(mp.mpc(p2), mp.mpf(m1), mp.mpf(m2), mp.mpf(mu2))
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
                return _complex_estimate(value)
            if m1 == 0 or m2 == 0:
                m = m1 + m2
                value = mp.mpc(1 - mp.log(m / mu2))
                return _complex_estimate(value)
            value = mp.mpc(1 - mp.log(m2 / mu2)
                           + m1 / (m1 - m2) * mp.log(m2 / m1))
            return _complex_estimate(value)
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
        return _complex_estimate(total, extra=_ulp_envelope(abs(total)) * 8)


def sample_winding(f: Callable[[complex], complex], corners: tuple[complex, complex],
                   subdivisions: int = 64, precision: int = 192) -> int:
    """Sampled quadrant winding of ``f`` around a rectangle.

    corners = (lower-left, upper-right).  The boundary is walked in
    subdivided point steps.  The local finite-difference check refuses
    coarse walks, but it is not an upper bound on the derivative between
    samples.  A two-quadrant jump also raises.  The result is useful as a
    deterministic numerical diagnostic and carries no root-count proof."""

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


def _b0_fin_complex(
    p2: mp.mpc,
    m1: mp.mpf,
    m2: mp.mpf,
    mu2: mp.mpf,
) -> DiagnosticEstimate:
    ieps = mp.mpc(0, mp.mpf(2) ** (16 - mp.mp.prec))
    a = p2
    b = -(p2 + m1 - m2)
    c = mp.mpc(m1) - ieps
    disc = mp.sqrt(b * b - 4 * a * c)
    x1 = (-b + disc) / (2 * a)
    x2 = (-b - disc) / (2 * a)
    total = mp.mpc(2) - mp.log(p2 / mu2)
    for x in (x1, x2):
        total += x * mp.log((x - 1) / x) - mp.log(x - 1)
    return _complex_estimate(total, extra=_ulp_envelope(abs(total)) * 8)
