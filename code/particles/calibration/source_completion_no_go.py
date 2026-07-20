#!/usr/bin/env python3
"""Exact two-completion witness for the current Higgs/top source reduct."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


class CompletionWitnessError(ValueError):
    """Raised when the source coordinate is outside the witness domain."""


@dataclass(frozen=True)
class CompletionSeparation:
    u: Fraction
    linear_higgs_ratio: Fraction
    born_higgs_ratio: Fraction
    linear_top_ratio: Fraction
    born_top_ratio: Fraction
    higgs_separation: Fraction
    top_separation: Fraction


def two_completion_separation(pixel: Fraction) -> CompletionSeparation:
    """Evaluate a_lin=u and a_Born=sqrt(u) through their squared outputs."""

    pixel = Fraction(pixel)
    u = Fraction(1) - pixel / 24
    if not Fraction(0) < u < Fraction(1):
        raise CompletionWitnessError("the witness requires 0<P<24")

    linear_a2 = u * u
    born_a2 = u
    linear_higgs = 2 * (1 - linear_a2)
    born_higgs = 2 * (1 - born_a2)
    linear_top = linear_a2 / 2
    born_top = born_a2 / 2
    higgs_separation = linear_higgs - born_higgs
    top_separation = born_top - linear_top
    expected_higgs = 2 * u * (1 - u)
    expected_top = u * (1 - u) / 2
    if higgs_separation != expected_higgs or top_separation != expected_top:
        raise CompletionWitnessError("exact separation identity failed")
    if higgs_separation <= 0 or top_separation <= 0:
        raise CompletionWitnessError("the two completions are not strictly separated")
    return CompletionSeparation(
        u=u,
        linear_higgs_ratio=linear_higgs,
        born_higgs_ratio=born_higgs,
        linear_top_ratio=linear_top,
        born_top_ratio=born_top,
        higgs_separation=higgs_separation,
        top_separation=top_separation,
    )
