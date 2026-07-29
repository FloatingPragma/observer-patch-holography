"""Tests for the counterterm pole solution and its certificates."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "producers"))

PAYLOAD = json.loads((ROOT / "outputs" / "counterterm_pole_solution.json").read_text(encoding="utf-8"))


def test_status_and_certificates() -> None:
    assert PAYLOAD["status"] == "PARTIAL_SOLVABLE_SLICE__SCALAR_XI_OBSTRUCTION_OPEN"
    assert PAYLOAD["slice_checks_passed"] is True
    assert PAYLOAD["acceptance_complete"] is False
    assert "solvable" in PAYLOAD["uv_cancellation"]
    assert "non-abelian binding" in PAYLOAD["census_binding"]


def test_census_gauge_poles() -> None:
    g1, g2 = sp.symbols("g1 g2", positive=True)
    dg1 = sp.sympify(PAYLOAD["solution_on_slice"]["dg1"], locals={"g1": g1, "g2": g2})
    dg2 = sp.sympify(PAYLOAD["solution_on_slice"]["dg2"], locals={"g1": g1, "g2": g2})
    assert sp.simplify(dg1 - sp.Rational(41, 12) * g1 ** 3) == 0
    assert sp.simplify(dg2 + sp.Rational(19, 12) * g2 ** 3) == 0


def test_abelian_ward_and_dmu2() -> None:
    assert PAYLOAD["checks"]["abelian_ward_ZB"]["passed"]
    assert PAYLOAD["checks"]["dmu2_xi_independent"]["passed"]


def test_scheme_amendment_recorded() -> None:
    assert "dxiW" in PAYLOAD["scheme"]["gauge_fixing"]
    assert "machine finding" in PAYLOAD["scheme"]["gauge_fixing"]


def test_residual_obstruction_recorded_exactly() -> None:
    g1, g2, v, lam, mu2 = sp.symbols("g1 g2 v lam mu2", positive=True)
    invariant = sp.sympify(PAYLOAD["residual_obstruction"]["invariant"],
                           locals={"g1": g1, "g2": g2, "v": v, "lam": lam, "mu2": mu2})
    expected = 24 * v ** 2 * (g1 ** 2 + 3 * g2 ** 2) * (2 * lam * v ** 2 - mu2)
    assert sp.simplify(invariant - expected) == 0
