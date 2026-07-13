"""Tests for the synthetic Cs-133 Feshbach scalarization fixture."""

from __future__ import annotations

import math

import numpy as np

import derive_cs133_feshbach_scalarization as lane


THEOREM_NAMES = {
    "CS133_FESHBACH_SPECTRAL_EQUIVALENCE_THEOREM",
    "CS133_FESHBACH_CHANNEL_SCALARIZATION_THEOREM",
    "CS133_FESHBACH_STRICT_MONOTONICITY_THEOREM",
    "CS133_SCALAR_SIGN_BRACKET_ROOT_THEOREM",
    "CS133_CLOCK_SCALAR_NORMAL_FORM_THEOREM",
    "CS133_FULL_OPERATOR_TO_GAP_INTERVAL_THEOREM",
    "CURRENT_OPH_CS133_GAP_NONENTAILMENT_THEOREM",
    "SOURCE_GEV_REQUIRES_OPERATIONAL_CLOCK_RATIO_THEOREM",
}


def test_projectors_have_ranks_seven_and_nine_and_are_orthogonal():
    h, p3, p4, _ = lane.make_synthetic_clock_model()
    c3 = lane.projector_checks(p3)
    c4 = lane.projector_checks(p4)
    assert c3["pass"] is True and c3["rank"] == 7
    assert c4["pass"] is True and c4["rank"] == 9
    assert float(np.linalg.norm(p3 @ p4, 2)) < 1e-12
    assert float(np.linalg.norm(h - h.T, 2)) < 1e-12


def test_hamiltonian_commutes_with_channel_rotations():
    h, _, _, _ = lane.make_synthetic_clock_model()
    rotation = lane.rotation_commutant_check(h)
    assert rotation["unitary_defect"] < 1e-12
    assert rotation["commutator_norm"] < 1e-12


def test_feshbach_operator_is_scalar_on_each_channel():
    h, p3, p4, _ = lane.make_synthetic_clock_model()
    sample = lane.channel_functions(h, p3, p4, 0.0)
    assert sample["scalar_defect_F3"] < 1e-12
    assert sample["scalar_defect_F4"] < 1e-12
    assert sample["cross_channel_defect"] < 1e-12


def test_derivatives_are_at_or_below_minus_one_on_sampled_interval():
    h, p3, p4, _ = lane.make_synthetic_clock_model()
    for z in lane.RESOLVENT_SAMPLE_GRID:
        sample = lane.channel_functions(h, p3, p4, z)
        assert sample["df3"] <= -1.0 + 1e-12
        assert sample["df4"] <= -1.0 + 1e-12
        assert sample["derivative_max_eigenvalue"] <= -1.0 + 1e-12


def test_sign_bracket_roots_match_exact_channel_eigenvalues():
    h, p3, p4, exact = lane.make_synthetic_clock_model()
    r3 = lane.bisect_channel(h, p3, p4, "F3", *lane.F3_BRACKET)
    r4 = lane.bisect_channel(h, p3, p4, "F4", *lane.F4_BRACKET)
    assert r3.f_lo > 0.0 > r3.f_hi
    assert r4.f_lo > 0.0 > r4.f_hi
    assert abs(r3.mid - exact["E3"]) < 1e-11
    assert abs(r4.mid - exact["E4"]) < 1e-11
    assert r4.lo - r3.hi > 0.0
    values = np.linalg.eigvalsh(h)
    assert int(np.sum(np.isclose(values, exact["E3"], atol=1e-11))) == 7
    assert int(np.sum(np.isclose(values, exact["E4"], atol=1e-11))) == 9


def test_gap_normal_form_identity_on_arbitrary_levels():
    form = lane.gap_normal_form_checks(1.0, 3.0)
    assert math.isclose(form["E0"], 2.125, rel_tol=0.0, abs_tol=1e-15)
    assert math.isclose(form["a_Cs"], 0.5, rel_tol=0.0, abs_tol=1e-15)
    assert math.isclose(form["epsilon"], 2.0, rel_tol=0.0, abs_tol=1e-15)
    assert abs(form["epsilon_minus_four_a"]) < 1e-15
    assert form["channel_identity_residual_F3"] < 1e-15
    assert form["channel_identity_residual_F4"] < 1e-15
    assert form["normal_form_residual"] < 1e-12
    assert form["K_eigenvalues"] == [-2.25, 1.75]


def test_artifact_passes_all_synthetic_checks_and_is_fail_closed():
    artifact = lane.build_artifact()
    assert artifact["checks_pass"] is True
    assert all(bool(value) for value in artifact["checks"].values())
    assert artifact["status"] == "SYNTHETIC_FIXTURE_THEOREM_CHECKS_ONLY"
    assert artifact["physical_source_prediction_ready"] is False
    assert artifact["public_clock_gap_promotion_allowed"] is False
    assert artifact["measured_particle_data_consumed"] is False
    assert artifact["scalarization_checks"]["max_scalar_defect"] < 1e-12
    assert artifact["scalarization_checks"]["max_cross_channel_defect"] < 1e-12
    assert artifact["scalarization_checks"]["max_derivative_value"] <= -1.0 + 1e-12


def test_artifact_registry_packets_and_reduction_block():
    artifact = lane.build_artifact()
    registry = artifact["theorem_registry"]
    assert set(registry) == THEOREM_NAMES
    assert all(isinstance(statement, str) and statement for statement in registry.values())
    assert artifact["open_source_packets"] == [
        "R_ALPHA_ATOMIC_SCHEME_PACKET",
        "R_ELECTRON_ABSOLUTE_RATIO_PACKET",
        "R_QCD_NUCLEAR_133CS_PACKET",
        "R_ATOMIC_FESHBACH_SCALARS_133CS",
        "R_CLOCK_REFINEMENT_LIMIT",
        "R_NO_G_CLOCK_DAG",
        "R_CLOCK_PROSPECTIVE_FREEZE",
    ]
    reduction = artifact["reduction"]
    assert "f3, f4" in reduction["pipeline"]
    assert "eps_Cs = E4 - E3" in reduction["pipeline"]
    assert "55-electron" in reduction["public_packet_form"]
