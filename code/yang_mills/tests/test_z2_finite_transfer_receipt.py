"""Tests for the Z2 finite ground-state-transform diagnostic."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

import z2_finite_transfer_receipt as z2  # noqa: E402


@pytest.fixture(scope="module")
def orbits() -> z2.Z2GaugeOrbits:
    return z2.Z2GaugeOrbits(2)


def test_orbit_count_and_free_action(orbits: z2.Z2GaugeOrbits) -> None:
    # 8 links, gauge group Z2^4 / global = 8 elements, free action
    assert orbits.n_links == 8
    assert len(orbits.gauge) == 8
    assert orbits.n_orbits == 32
    assert np.all(orbits.orbit_size == 8)


def test_free_control_receipt_is_exact(orbits: z2.Z2GaugeOrbits) -> None:
    """At beta_s = 0 the kinetic kernel factorises as
    prod_l (e^{beta_t} I + e^{-beta_t} X_l) = const * exp(b sum_l X_l) with
    b = log(coth beta_t) / 2, so the receipt holds exactly with every rate
    equal to log(coth beta_t)."""
    r = z2.evaluate(orbits, "wilson", beta_s=0.0, beta_t=0.5)
    fit = r["constant_rate_fit"]
    rate = math.log(1.0 / math.tanh(0.5))
    assert fit["relative_frobenius_residual"] < 1e-10
    assert fit["rate_min"] == pytest.approx(rate, abs=1e-9)
    assert fit["rate_max"] == pytest.approx(rate, abs=1e-9)
    assert r["spectral"]["gap_H"] == pytest.approx(2 * rate, abs=1e-9)
    assert r["fiber_dependent_rates"]["spread_max_over_min"] == pytest.approx(1.0, abs=1e-9)
    assert r["dobrushin"]["eta_star"] < 1e-9
    assert r["spectral"]["gap_unit_rate_heat_bath"] == pytest.approx(2.0, abs=1e-9)


def test_interacting_wilson_receipt_fails(orbits: z2.Z2GaugeOrbits) -> None:
    r = z2.evaluate(orbits, "wilson", beta_s=0.5, beta_t=0.5)
    assert r["doob_generator_is_markov"]
    assert r["constant_rate_fit"]["relative_frobenius_residual"] > 1e-2
    assert r["fiber_dependent_rates"]["offdiagonal_mass_outside_single_flip"] > 1.0
    assert r["fiber_dependent_rates"]["spread_max_over_min"] > 1.01


def test_kogut_susskind_single_flip_exact_but_fiber_dependent(orbits: z2.Z2GaugeOrbits) -> None:
    r = z2.evaluate(orbits, "kogut_susskind", lam=1.0)
    fib = r["fiber_dependent_rates"]
    assert r["doob_generator_is_markov"]
    assert r["doob_generator_offdiagonal_nonpositive"]
    assert fib["offdiagonal_mass_outside_single_flip"] < 1e-9
    assert fib["all_rates_positive"]
    assert fib["spread_max_over_min"] > 1.01
    # exact rate formula c = lam (r + 1/r) with r = Omega(o)/Omega(o') >= 2 lam
    assert fib["rate_min"] >= 2.0 - 1e-9


def test_committed_receipt_matches_code() -> None:
    path = HERE / "receipts" / "z2_finite_transfer_receipt.json"
    receipt = json.loads(path.read_text(encoding="utf-8"))
    assert receipt["schema"] == z2.SCHEMA
    assert receipt["physical_clay_receipt"] is False
    controls = [r for r in receipt["runs"] if r.get("control") == "beta_s_zero"]
    assert controls and all(
        r["constant_rate_fit"]["relative_frobenius_residual"] < 1e-10 for r in controls
    )
    interacting = [r for r in receipt["runs"] if r["transfer"] == "wilson" and "control" not in r]
    assert interacting and all(
        r["constant_rate_fit"]["relative_frobenius_residual"] > 1e-3 for r in interacting
    )
    # reproduce one L=2 row exactly
    orbits = z2.Z2GaugeOrbits(2)
    row = next(r for r in interacting if r["L"] == 2)
    fresh = z2.evaluate(orbits, "wilson", **row["parameters"])
    assert fresh["constant_rate_fit"]["relative_frobenius_residual"] == pytest.approx(
        row["constant_rate_fit"]["relative_frobenius_residual"], rel=1e-9
    )
