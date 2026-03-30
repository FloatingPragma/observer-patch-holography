#!/usr/bin/env python3
"""Concrete scaffold for the missing same-label left-handed quark orbit provider."""

from __future__ import annotations

from typing import Sequence

from sigma_ud_orbit_provider_interface import CanonicalToken, OrbitElement, SigmaUDOrbitProvider


class MissingSigmaUDOrbitProvider(SigmaUDOrbitProvider):
    """Explicit placeholder until the live solver can emit same-label left-handed orbit elements."""

    def enumerate_relative_sheets_d12(self) -> Sequence[CanonicalToken]:
        raise NotImplementedError(
            "Current local corpus still lacks an emitted finite same-label left-handed Sigma_ud orbit."
        )

    def evaluate_relative_sheet(self, token: CanonicalToken) -> OrbitElement:
        raise NotImplementedError(
            "Current local corpus still lacks a concrete same-label left-handed sigma -> CKM evaluator."
        )
