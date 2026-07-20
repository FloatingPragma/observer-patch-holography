#!/usr/bin/env python3
"""Tests for the Einstein branch closure receipts (GitHub #526-#528, #503, #578)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from einstein_closure_receipts import (  # noqa: E402
    ETA,
    baseline_countermodel_receipt,
    blockwise_state,
    bulk_entropy,
    charges_of,
    central_z,
    design_matrix,
    edge_entropy,
    entropy,
    eta_project_out,
    first_law_receipt,
    generic_null_directions,
    maxent_multiplier_receipt,
    reconstruct_from_charges,
    tomography_receipt,
)
from realized_branch_receipts import (  # noqa: E402
    build_report,
    receipt_report_for,
    tree_packet_net_records,
)


# ---------------------------------------------------------------------------
# null tomography (thm:null-tomography, prop:no-local-stress-countermodel)
# ---------------------------------------------------------------------------

def test_design_rank_is_nine_and_eta_is_the_kernel():
    dirs = generic_null_directions(12)
    a = design_matrix(dirs)
    assert np.linalg.matrix_rank(a) == 9
    # eta is in the kernel of the charge map: eta(k,k) = 0
    assert np.allclose(charges_of(ETA, dirs), 0.0, atol=1e-12)


def test_consistent_charges_reconstruct_tracefree_part_exactly():
    r = tomography_receipt()
    assert r["design_rank"] == 9
    assert r["consistent_residual"] < 1e-10
    assert r["tracefree_error"] < 1e-9


def test_inconsistent_charges_have_irreducible_residual():
    r = tomography_receipt()
    assert r["inconsistent_residual"] > 0.1


def test_eta_ambiguity_is_exactly_the_reconstruction_freedom():
    rng = np.random.default_rng(3)
    dirs = generic_null_directions(10, seed=5)
    m = rng.normal(size=(4, 4))
    t = (m + m.T) / 2.0
    for phi in (0.0, 1.7, -0.4):
        t_hat, res = reconstruct_from_charges(charges_of(t + phi * ETA, dirs), dirs)
        assert res < 1e-10
        assert np.allclose(eta_project_out(t_hat), eta_project_out(t), atol=1e-8)


# ---------------------------------------------------------------------------
# first law with edge term (thm:bulk-edge-central-first-law)
# ---------------------------------------------------------------------------

def test_direct_sum_state_and_entropy_convention_are_explicit():
    ps = [0.25, 0.75]
    bulk_states = [np.diag([0.8, 0.2]), np.array([[1.0]])]
    edge_dims = [2, 3]
    rho = blockwise_state(ps, bulk_states, edge_dims)

    expected = np.zeros((7, 7))
    expected[:4, :4] = ps[0] * np.kron(bulk_states[0], np.eye(2) / 2)
    expected[4:, 4:] = ps[1] * np.kron(bulk_states[1], np.eye(3) / 3)
    assert np.allclose(rho, expected)
    assert abs(np.trace(rho) - 1.0) < 1e-12
    assert abs(edge_entropy(ps, edge_dims)
               - (0.25 * np.log(2) + 0.75 * np.log(3))) < 1e-12
    assert abs(bulk_entropy(ps, bulk_states)
               - (-0.25 * np.log(0.25) - 0.75 * np.log(0.75)
                  + 0.25 * entropy(bulk_states[0]))) < 1e-12
    assert abs(entropy(rho) - bulk_entropy(ps, bulk_states)
               - edge_entropy(ps, edge_dims)) < 1e-12


def test_central_generator_contains_only_log_edge_dimensions():
    z = central_z([2, 1], [2, 3], [np.log(2), np.log(3)])
    expected = np.diag([np.log(2)] * 4 + [np.log(3)] * 3)
    assert np.allclose(z, expected)


def test_first_law_and_boxed_split_exact_on_declared_normalization():
    r = first_law_receipt()
    assert r["base_entropy_split_defect"] < 1e-12
    assert r["varied_entropy_split_defect"] < 1e-12
    assert r["first_law_defect"] < 1e-3            # O(eps) second-order remainder
    assert r["edge_identification_defect"] < 1e-3  # delta<Z> = delta S_edge
    assert r["split_identity_defect"] < 1e-3       # delta S = 2pi d<B> + dS_edge
    assert r["bulk_identity_defect"] < 1e-3        # dS_bulk = 2pi d<B>


def test_bulk_identity_includes_moving_sector_weights():
    # H(p) is in S_bulk and -log p_a is in B, so moving weights are covered.
    r = first_law_receipt(move_weights=True)
    assert r["bulk_identity_defect"] < 1e-3
    assert r["predicted_bulk_defect"] < 1e-12
    # Fixed sector weights remain covered too.
    r0 = first_law_receipt(move_weights=False)
    assert r0["bulk_identity_defect"] < 1e-3
    assert r0["predicted_bulk_defect"] < 1e-12


def test_mismatched_edge_normalization_breaks_edge_identification():
    wrong = [0.0, 0.0]  # z_alpha = 0 instead of log d_alpha
    r = first_law_receipt(z_weights=wrong)
    # the first law itself still holds (fixed operators)...
    assert r["first_law_defect"] < 1e-3
    # ...but the edge identification fails by exactly the predicted defect
    assert r["predicted_edge_defect"] > 0.1
    assert abs(r["edge_identification_defect"] - r["predicted_edge_defect"]) < 1e-2
    # The same mismatch makes B differ from the declared bulk generator.
    assert abs(r["bulk_identity_defect"] - r["predicted_bulk_defect"]) < 1e-2


# ---------------------------------------------------------------------------
# MaxEnt multiplier identity (thm:maxent-lagrange-stationarity)
# ---------------------------------------------------------------------------

def test_maxent_multiplier_identity():
    for lam in (0.4, 1.3, 2.2):
        r = maxent_multiplier_receipt(lam=lam)
        assert r["multiplier_defect"] < 1e-3, r


# ---------------------------------------------------------------------------
# baseline countermodel (prop:baseline-countermodels(iv))
# ---------------------------------------------------------------------------

def test_baseline_invisible_to_first_variations():
    r = baseline_countermodel_receipt(c=0.7)
    assert r["variation_data_difference"] < 1e-15
    assert r["baseline_difference"] > 1.0
    assert r["eta_alignment"] < 1e-12


# ---------------------------------------------------------------------------
# realized-branch evaluation (rem:realized-branch-status)
# ---------------------------------------------------------------------------

def test_tree_packet_net_fails_spherical_incidence():
    records = tree_packet_net_records()
    assert len(records) == 3
    rep = receipt_report_for(records)
    assert rep["n_triangles"] == 0
    assert rep["euler_characteristic"] == 1
    assert rep["spherical_incidence_receipt"] is False
    assert rep["surface_classification"] == "NOT_A_CLOSED_SURFACE"


def test_realized_branch_report_status_is_open():
    report = build_report()
    assert report["realized_geometric_branch_certified_nonempty"] is False
    assert report["status"].startswith("OPEN")
    designed = report["evaluations"]["designed_cellulated_sphere_icosahedron"]
    assert designed["spherical_incidence_receipt"] is True
    assert designed["provenance"].startswith("declared_geometry")
