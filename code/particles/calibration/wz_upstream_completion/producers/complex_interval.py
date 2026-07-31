"""Directed rectangular complex interval arithmetic with sheet gates.

Complex values are enclosed by axis-aligned rectangles whose real and
imaginary parts are mpmath directed real intervals, so every
elementary operation propagates outward-rounded enclosures.  The
multivalued operations carry explicit sheet gates: a logarithm or
square root refuses, by raising SheetError, whenever its argument
rectangle meets zero or the principal branch cut cannot be excluded,
and a division refuses whenever the denominator rectangle meets zero.
A refused operation is a subdivision signal for the caller, never a
silent widening.

The argument enclosure uses the interval two-argument arctangent,
which returns the full circle when the rectangle meets the closed
negative real axis, so cut exclusion is checked directly on the
rectangle geometry before the tight enclosure is trusted.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any

from mpmath import iv


class SheetError(ValueError):
    """A branch, cut, or zero condition failed on an interval operation."""


def _iv_fraction(value: Fraction | int) -> Any:
    fraction = Fraction(value)
    return iv.mpf(fraction.numerator) / iv.mpf(fraction.denominator)


def _iv_square(x: Any) -> Any:
    """Directed interval square without the dependency defect.

    Generic interval multiplication of an interval with itself loses
    the information that both factors are one variable, so a
    zero-straddling interval would receive a negative lower bound.
    The square of an interval is computed from its endpoint moduli
    with outward rounding and is nonnegative by construction."""

    lo_abs = abs(x.a) if abs(x.a) <= abs(x.b) else abs(x.b)
    hi_abs = abs(x.a) if abs(x.a) >= abs(x.b) else abs(x.b)
    hi = (iv.mpf([hi_abs, hi_abs]) * iv.mpf([hi_abs, hi_abs])).b
    if 0 in x:
        return iv.mpf([0, hi])
    lo = (iv.mpf([lo_abs, lo_abs]) * iv.mpf([lo_abs, lo_abs])).a
    return iv.mpf([lo, hi])


class CInterval:
    __slots__ = ("re", "im")

    def __init__(self, re: Any, im: Any) -> None:
        self.re = re
        self.im = im

    @classmethod
    def from_fraction(cls, re: Fraction | int, im: Fraction | int = 0) -> "CInterval":
        return cls(_iv_fraction(re), _iv_fraction(im))

    @classmethod
    def box(
        cls,
        re_lo: Fraction,
        re_hi: Fraction,
        im_lo: Fraction,
        im_hi: Fraction,
    ) -> "CInterval":
        lo_re = _iv_fraction(re_lo)
        hi_re = _iv_fraction(re_hi)
        lo_im = _iv_fraction(im_lo)
        hi_im = _iv_fraction(im_hi)
        return cls(
            iv.mpf([lo_re.a, hi_re.b]),
            iv.mpf([lo_im.a, hi_im.b]),
        )

    def __add__(self, other: "CInterval") -> "CInterval":
        return CInterval(self.re + other.re, self.im + other.im)

    def __sub__(self, other: "CInterval") -> "CInterval":
        return CInterval(self.re - other.re, self.im - other.im)

    def __neg__(self) -> "CInterval":
        return CInterval(-self.re, -self.im)

    def __mul__(self, other: "CInterval") -> "CInterval":
        return CInterval(
            self.re * other.re - self.im * other.im,
            self.re * other.im + self.im * other.re,
        )

    def abs2(self) -> Any:
        return _iv_square(self.re) + _iv_square(self.im)

    def contains_zero(self) -> bool:
        return bool(0 in self.re and 0 in self.im)

    def __truediv__(self, other: "CInterval") -> "CInterval":
        denominator = other.abs2()
        if 0 in denominator:
            raise SheetError("division by an interval containing zero")
        numerator = self * CInterval(other.re, -other.im)
        return CInterval(numerator.re / denominator, numerator.im / denominator)

    def pow_int(self, exponent: int) -> "CInterval":
        if exponent < 0:
            raise SheetError("negative interval powers are formed by division")
        result = CInterval.from_fraction(1)
        base = self
        remaining = exponent
        while remaining:
            if remaining & 1:
                result = result * base
            base = base * base
            remaining >>= 1
        return result

    def excludes_principal_cut(self) -> bool:
        """The rectangle avoids zero and the closed negative real axis."""

        if self.contains_zero():
            return False
        if self.re.a > 0:
            return True
        if self.im.a > 0 or self.im.b < 0:
            return True
        return False

    def arg(self) -> Any:
        if not self.excludes_principal_cut():
            raise SheetError("argument rectangle meets the principal cut")
        return iv.atan2(self.im, self.re)

    def log(self, arg_width_gate: Any) -> "CInterval":
        angle = self.arg()
        if angle.delta > arg_width_gate:
            raise SheetError("argument enclosure wider than the declared gate")
        return CInterval(iv.log(self.abs2()) / iv.mpf(2), angle)

    def sqrt(self, arg_width_gate: Any) -> "CInterval":
        angle = self.arg()
        if angle.delta > arg_width_gate:
            raise SheetError("argument enclosure wider than the declared gate")
        radius = iv.sqrt(iv.sqrt(self.abs2()))
        half = angle / iv.mpf(2)
        return CInterval(radius * iv.cos(half), radius * iv.sin(half))

    def width(self) -> Any:
        return self.re.delta if self.re.delta > self.im.delta else self.im.delta

    def midpoint(self) -> complex:
        return complex(
            (float(self.re.a) + float(self.re.b)) / 2,
            (float(self.im.a) + float(self.im.b)) / 2,
        )

    def hull(self, other: "CInterval") -> "CInterval":
        return CInterval(
            iv.mpf([min(self.re.a, other.re.a), max(self.re.b, other.re.b)]),
            iv.mpf([min(self.im.a, other.im.a), max(self.im.b, other.im.b)]),
        )

    def __repr__(self) -> str:
        return f"CInterval(re={self.re}, im={self.im})"
