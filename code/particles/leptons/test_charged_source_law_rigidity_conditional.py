"""Tests for the conditional CFQ packet-weight rigidity theorem."""

from __future__ import annotations

import json
from fractions import Fraction

import mpmath as mp

import derive_charged_source_law_rigidity_conditional as lane


def artifact() -> dict:
    face_raw = lane.FACE_RECEIPT.read_bytes()
    declared_raw = lane.DECLARED_MAP_RECEIPT.read_bytes()
    return lane.build_artifact(
        json.loads(face_raw),
        lane.sha256(face_raw),
        json.loads(declared_raw),
        lane.sha256(declared_raw),
    )


def test_exact_register_graph_and_path_arithmetic():
    result = artifact()
    graphs = result["declared_finite_graph_models"]["graphs"]
    assert {
        name: (row["nodes"], row["edges"], row["connected"])
        for name, row in graphs.items()
    } == {
        "kappa_0": (50, 60, True),
        "kappa_1": (31, 30, True),
        "kappa_2": (10, 9, True),
        "chi_0": (512, 2304, True),
        "chi_1": (77, 340, True),
        "zeta_0": (21, 20, True),
        "zeta_1": (27, 99, True),
        "zeta_2": (5, 5, True),
    }
    assert [
        row["exact_trace_coefficient"]
        for row in result["conditional_cfq_packet"]["paths"]
    ] == ["1/50", "-1/31", "-1/310", "1/512", "1/77", "1/21", "1/27", "1/135"]


def test_packet_reproduces_prior_declared_kernel_only_conditionally():
    result = artifact()
    residuals = result["declared_map_agreement"][
        "coefficient_residuals_against_prior_receipt"
    ].values()
    assert max(abs(mp.mpf(value)) for value in residuals) < mp.mpf("1e-75")
    assert result["theorem_scope"]["rigidity_scope"].startswith(
        "Uniqueness holds inside the stipulated CFQ admissible class"
    )
    assert result["declared_map_agreement"]["downstream_mass_values_recomputed_here"] is False


def test_dimension_rank_and_multiplicity_are_assumptions_not_derived_freedoms():
    registers = lane.packet_registers()
    paths = lane.packet_paths()
    base = lane.all_path_weights(paths, registers)
    perturbed_registers = {
        **registers,
        "kappa_0": lane.Register("kappa_0", 51),
    }
    assert lane.all_path_weights(paths, perturbed_registers)["kappa_write"] != base[
        "kappa_write"
    ]
    perturbed_paths = (
        lane.PathRule(
            paths[0].name,
            paths[0].block,
            paths[0].monomial,
            paths[0].sign,
            paths[0].registers,
            multiplicity=2,
        ),
        *paths[1:],
    )
    assert lane.all_path_weights(perturbed_paths, registers)["kappa_write"] == Fraction(
        1, 25
    )
    assert lane.multiplier_admissible_inside_cfq(Fraction(1, 1)) is True
    assert lane.multiplier_admissible_inside_cfq(Fraction(2, 1)) is False


def test_every_physical_gate_is_open_and_promotion_flags_cannot_regress():
    result = artifact()
    assert not any(result["conditional_cfq_packet"]["physical_gates"].values())
    assert all(
        value is None
        for value in result["conditional_cfq_packet"]["evidence_hashes"].values()
    )
    assert result["historical_charged_target_informed"] is True
    assert result["global_source_only"] is False
    assert result["branch_tuple_coherent"] is False
    assert result["mass_scheme_certified"] is False
    assert result["public_prediction_promotion_allowed"] is False
    assert result["record_and_maxent_boundary"]["central_record_bridge_open"] is True
    assert result["separate_conditional_corollaries"]["closed_by_core_cfq_theorem"] is False
    assert result["checks_pass"] is True


def test_committed_receipt_is_byte_reproducible():
    result = artifact()
    expected = (json.dumps(result, indent=2, sort_keys=True) + "\n").encode()
    assert lane.DEFAULT_OUT.read_bytes() == expected
