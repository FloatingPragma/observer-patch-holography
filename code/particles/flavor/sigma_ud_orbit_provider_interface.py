#!/usr/bin/env python3
"""Interface skeleton for the true missing quark relative-sheet orbit provider."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, Sequence


@dataclass(frozen=True)
class CanonicalToken:
    token: str


@dataclass(frozen=True)
class CKMTuple:
    theta_12: float
    theta_23: float
    theta_13: float
    delta_ckm: float
    jarlskog: float


@dataclass(frozen=True)
class OrbitElement:
    sigma_id: str
    canonical_token: CanonicalToken
    U_u_left: object
    U_d_left: object
    V_CKM: object
    ckm_invariants: CKMTuple


class SigmaUDOrbitProvider(Protocol):
    def enumerate_relative_sheets_d12(self) -> Sequence[CanonicalToken]:
        raise NotImplementedError(
            "Current local corpus does not expose finite Sigma_ud representatives."
        )

    def evaluate_relative_sheet(self, token: CanonicalToken) -> OrbitElement:
        raise NotImplementedError(
            "Current local corpus does not expose a same-label left-handed sigma -> CKM evaluator."
        )


def build_sigma_ud_orbit(provider: SigmaUDOrbitProvider) -> list[OrbitElement]:
    orbit: list[OrbitElement] = []
    for idx, token in enumerate(provider.enumerate_relative_sheets_d12()):
        evaluation = provider.evaluate_relative_sheet(token)
        orbit.append(
            OrbitElement(
                sigma_id=evaluation.sigma_id or f"sigma_{idx}",
                canonical_token=evaluation.canonical_token,
                U_u_left=evaluation.U_u_left,
                U_d_left=evaluation.U_d_left,
                V_CKM=evaluation.V_CKM,
                ckm_invariants=evaluation.ckm_invariants,
            )
        )
    return orbit
