#!/usr/bin/env python3
"""Tests for the null-net standardness receipts (GitHub #524)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from null_net_standardness import (  # noqa: E402
    assembly_receipts,
    commutant_cyclicity_implies_separating,
    full_algebra_standard,
    gibbs_nonlocality_witness,
    markov_chain_state,
    markov_modular_split_receipt,
    modular_hamiltonian,
    one_sided_split_defect,
    partial_trace,
    subalgebra_separating,
)


# ---------------------------------------------------------------------------
# stagewise standardness (finite clauses of thm:null-net-standardness)
# ---------------------------------------------------------------------------

def test_faithful_state_gives_cyclic_separating_gns_vector():
    rng = np.random.default_rng(0)
    m = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    rho = m @ m.conj().T + 0.2 * np.eye(4)
    rho = rho / np.trace(rho)
    cyclic, separating = full_algebra_standard(rho)
    assert cyclic and separating


def test_nonfaithful_state_fails_standardness():
    rho = np.diag([0.5, 0.5, 0.0, 0.0])
    cyclic, separating = full_algebra_standard(rho)
    assert not cyclic
    assert not separating


def test_separation_passes_to_halfline_subalgebras():
    rng = np.random.default_rng(1)
    dims = [2, 2, 2]
    d = 8
    m = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    rho = m @ m.conj().T + 0.2 * np.eye(d)
    rho = rho / np.trace(rho)
    assert subalgebra_separating(rho, dims, region=[0])
    assert subalgebra_separating(rho, dims, region=[0, 1])


def test_commutant_cyclicity_implies_separation_mechanism():
    b_cyclic, a_separated = commutant_cyclicity_implies_separating(d=3)
    assert b_cyclic and a_separated


# ---------------------------------------------------------------------------
# Markov modular locality (thm:markov-modular-locality)
# ---------------------------------------------------------------------------

def test_markov_split_is_exact_and_gibbs_split_is_not():
    receipts = markov_modular_split_receipt()
    assert receipts["markov_split_defect"] < 1e-8
    assert receipts["gibbs_split_defect"] > 0.05
    assert receipts["endpoint_lipschitz_growth"] < 1e-8


def test_split_defect_functional_detects_cross_terms():
    rng = np.random.default_rng(2)
    h_l = rng.normal(size=(2, 2))
    h_l = h_l + h_l.T
    h_r = rng.normal(size=(2, 2))
    h_r = h_r + h_r.T
    k_split = np.kron(h_l, np.eye(2)) + np.kron(np.eye(2), h_r)
    assert one_sided_split_defect(k_split, 2, 2) < 1e-10
    cross = np.kron(np.array([[0.0, 1.0], [1.0, 0.0]]),
                    np.array([[0.0, 1.0], [1.0, 0.0]]))
    assert one_sided_split_defect(k_split + 0.3 * cross, 2, 2) > 0.2


def test_markov_state_has_zero_collar_cmi():
    from null_net_standardness import collar_cmi
    rho = markov_chain_state(2, 2, 2, 2, weights=[1.0], seed=11)
    cmi = collar_cmi(rho, [2, 2, 2, 2], a=[0], b=[1, 2], c=[3])
    assert abs(cmi) < 1e-8


# ---------------------------------------------------------------------------
# counterexample boundary (prop:gibbs-modular-nonlocality-witness)
# ---------------------------------------------------------------------------

def test_gibbs_nonlocality_witness():
    w = gibbs_nonlocality_witness(beta=1.0, g=1.0)
    # nearest-neighbour Gibbs: endpoint coupling in K_I is genuinely nonzero
    assert w["endpoint_coupling_norm"] > 1e-3
    # and the state is strictly outside the Markov branch
    assert w["collar_cmi"] > 1e-4
    # exactly Markov comparison state: both functionals vanish
    assert w["markov_endpoint_coupling_norm"] < 1e-9
    assert w["markov_collar_cmi"] < 1e-9


def test_reduced_modular_hamiltonian_is_defined():
    rho = markov_chain_state(2, 2, 2, 2, weights=[0.5, 0.5], seed=3)
    rho_i = partial_trace(rho, [2, 2, 2, 2], keep=[0, 1])
    k_i = modular_hamiltonian(rho_i)
    assert np.allclose(k_i, k_i.conj().T, atol=1e-9)


# ---------------------------------------------------------------------------
# four-translation assembly (thm:four-translation-assembly)
# ---------------------------------------------------------------------------

def test_assembly_lie_algebra_receipts():
    r = assembly_receipts()
    assert r["commutators"] < 1e-12          # null translations commute
    assert r["covariance"] < 1e-9            # Lambda P(q) Lambda^-1 = P(Lambda q)
    assert r["linearity"] < 1e-9             # dependent null families assemble
    assert r["dual_cone_inside_min"] > 0.0   # future cone dual positivity
    assert r["dual_cone_outside_min"] < 0.0  # spacelike vectors excluded
