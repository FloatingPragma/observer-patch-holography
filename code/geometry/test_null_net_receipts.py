#!/usr/bin/env python3
"""Tests for the null-net receipt instrumentation (#503, #524)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from null_net_receipts import (  # noqa: E402
    hsm_compression_receipt,
    instrument_null_net,
    lie_closure_receipt,
    mixed_gns_cauchy_receipt,
    nti_receipt,
    separating_modulus_receipt,
    weak_additivity_receipt,
)


def test_nti_and_weak_additivity_are_witnessed():
    for n in (16, 32):
        assert nti_receipt(n)["nontrivial"]
        assert weak_additivity_receipt(n)["covers_ring"]


def test_separating_modulus_positive_at_every_stage():
    r = separating_modulus_receipt(32)
    assert r["strictly_faithful"]
    # the gap is astronomically small but strictly positive: this is why
    # the instrumentation needs extended precision
    assert r["occupation_gap_log10"] < -10


def test_mixed_gns_cauchy():
    r = mixed_gns_cauchy_receipt((16, 32, 64))
    assert r["cauchy_decreasing"]


def test_hsm_compression_has_consistent_sign_asymmetry():
    r32 = hsm_compression_receipt(32)
    r64 = hsm_compression_receipt(64)
    assert r32["asymmetry_ratio"] > 3.0
    assert r64["asymmetry_ratio"] > 3.0
    assert r32["compressing_sign"] == r64["compressing_sign"]


def test_hsm_real_packet_is_the_null_control():
    # a real (non-chiral) packet is exactly time-reversal symmetric: the
    # asymmetry must vanish, confirming the chiral carrier is what detects
    # the half-sided compression
    import numpy as np
    from scipy.linalg import expm
    from modular_clock_instrumentation import arc_entanglement_hamiltonian
    from null_net_receipts import embed
    n, m_a, m_b = 32, 16, 8
    h_a = arc_entanglement_hamiltonian(n, m_a)
    h_b = embed(arc_entanglement_hamiltonian(n, m_b), m_a)
    xs = np.arange(m_a)
    psi = np.exp(-((xs - m_b / 2.0) ** 2) / (2 * (m_b / 6.0) ** 2)).astype(complex)
    psi[m_b:] = 0.0
    psi /= np.linalg.norm(psi)

    def leakage(t):
        out = expm(1j * h_a * t) @ expm(-1j * h_b * t) @ psi
        return float(np.sum(np.abs(out[m_b:]) ** 2))

    assert abs(leakage(0.12) - leakage(-0.12)) < 1e-12


def test_lie_closure_percent_level_and_resummation_control():
    r = lie_closure_receipt(32)
    assert r["relative_residual"] < 0.02
    # the unresummed range-2 truncation is several times worse: the
    # resummation is load-bearing
    assert (r["relative_residual_unresummed_control"]
            > 5.0 * r["relative_residual"])


def test_instrumented_null_net_verdicts():
    report = instrument_null_net(rings=(16, 32))
    w = report["receipts_witnessed"]
    assert w["nti"] and w["weak_additivity"]
    assert w["separating_faithfulness"] and w["mixed_gns_cauchy"]
    assert w["hsm_compression_one_particle"]
    assert w["modular_lie_closure_percent_level"]
    # the rate clause must be recorded as open in the report
    assert any("convergence rate" in p or "rate" in p
               for p in report["receipts_pending"])
