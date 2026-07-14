#!/usr/bin/env python3
"""Tests for the Ward-projected payload harness (generator G1).

LABEL: development bracket, non-blind environment; the protocol pass requires
an isolated re-run.
"""

from __future__ import annotations

import json
import math
import sys
from decimal import Decimal
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import payload_harness as ph
import spectral_modules as sm


@pytest.fixture(scope="module")
def ep() -> ph.EvaluationPoint:
    return ph.build_evaluation_point()


def _canonical(payload: dict) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"))


def test_harness_determinism(ep: ph.EvaluationPoint) -> None:
    module = sm.make_pqcd("lane_central", 4.0, "zero", 3)
    first = ph.emit_delta_source(module, ep)
    second = ph.emit_delta_source(module, ep)
    assert first["content_sha256"] == second["content_sha256"]
    assert _canonical(first) == _canonical(second)
    assert first["label"] == ph.NON_BLIND_LABEL


def test_synthetic_atom_round_trip(ep: ph.EvaluationPoint) -> None:
    # The no-go countermodel atoms: unit weight at y=2 and y=3 have kernel
    # moments 1/6 and 1/12 before the 1/(3*pi) prefactor.
    module = sm.SyntheticAtomModule([(2.0, 1.0), (3.0, 1.0)])
    payload = ph.emit_delta_source(module, ep)
    expected = (1.0 / 6.0 + 1.0 / 12.0) / (3.0 * math.pi)
    assert payload["components_alpha_inv"]["delta_had"] == pytest.approx(
        expected, rel=1e-14
    )
    assert payload["diagnostics"]["positivity_ok"] is True

    negative = sm.SyntheticAtomModule([(2.0, -1.0)])
    flagged = ph.emit_delta_source(negative, ep)
    assert flagged["diagnostics"]["positivity_ok"] is False


def test_parton_moment_matches_paper_math_kernel(ep: ph.EvaluationPoint) -> None:
    sys.path.insert(0, str(ph.P_DERIVATION))
    from paper_math import PaperMathContext

    ctx = PaperMathContext(precision=30, su2_cutoff=10, su3_cutoff=8)
    for name, q2, nc in (("c", 4.0 / 9.0, 3), ("b", 1.0 / 9.0, 3), ("e", 1.0, 1)):
        mass = ep.quark_masses[name] if name in ep.quark_masses else ep.lepton_masses[name]
        reference = float(
            ctx.fermion_transport_kernel_exact(
                Decimal(ep.mz_run), Decimal(mass), Decimal(str(q2)), Decimal(nc)
            )
        )
        y_thr = 4.0 * (mass / ep.mz_run) ** 2
        assert ph.parton_moment(y_thr, nc * q2) == pytest.approx(reference, rel=1e-11)


def test_density_integrator_matches_closed_form(ep: ph.EvaluationPoint) -> None:
    mass = ep.quark_masses["c"]
    y_thr = 4.0 * (mass / ep.mz_run) ** 2
    ncq2 = 3.0 * 4.0 / 9.0

    def rho(y: float) -> float:
        beta = ph.beta_of_y(y_thr, y)
        return ncq2 * beta * (3.0 - beta * beta) / 2.0

    numeric = ph.density_moment(
        rho, y_thr, 1.0e12, sqrt_left=True, gauss_n=48, splits_per_decade=4
    ) + ph.tail_moment(1.0e12, ncq2)
    closed = ph.parton_moment(y_thr, ncq2)
    assert numeric == pytest.approx(closed, rel=1e-10)


def test_chain_baseline_reproduced(ep: ph.EvaluationPoint) -> None:
    # Values from code/P_derivation/runtime/full_p_alpha_report_current.json
    # (Stage-5 internal chain, mode thomson_structured_running).
    assert ph.lepton_transport(ep) == pytest.approx(4.309397144816788, rel=1e-12)
    assert ph.quark_naive_transport(ep) == pytest.approx(4.934816164436291, rel=1e-12)
    assert ph.implemented_screen(ep) == pytest.approx(0.886997575746281, rel=1e-12)
    free = ph.emit_delta_source(sm.make_parton_free(), ep)
    assert free["diagnostics"]["s_effective"] == pytest.approx(1.0, rel=1e-14)


def test_bracket_reproducibility(ep: ph.EvaluationPoint) -> None:
    import run_bracket

    first = run_bracket.build_bracket(ep, fast=True)
    second = run_bracket.build_bracket(ep, fast=True)
    assert first["content_sha256"] == second["content_sha256"]
    bracket = first["bracket"]["delta_source_alpha_inv"]
    assert bracket["lo"] < bracket["hi"]
    assert first["label"] == ph.NON_BLIND_LABEL
    assert first["precision_wall"]["width_over_tolerance"] > 1.0e6
