#!/usr/bin/env python3
"""Tests for the finite-screen shock-sign and graph-spectrum receipt."""

from __future__ import annotations

import copy
import json
import pathlib
import subprocess
import sys

import numpy as np
import pytest


HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCRIPT = HERE / "derive_screen_shock_spectrum.py"
ARTIFACT = HERE / "runs" / "screen_shock_spectrum.json"
sys.path.insert(0, str(HERE))

from derive_screen_shock_spectrum import (  # noqa: E402
    ReceiptValidationError,
    _icosahedron_graphs,
    _laplacian,
    build_receipt,
    receipt_failures,
    validate_receipt,
)


def _run_to(path: pathlib.Path) -> dict:
    subprocess.run(
        [sys.executable, str(SCRIPT), "--output", str(path)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return json.loads(path.read_text(encoding="utf-8"))


def test_screen_shock_receipt_is_byte_exact_and_lf_only(
    tmp_path: pathlib.Path,
) -> None:
    rebuilt = tmp_path / ARTIFACT.name
    _run_to(rebuilt)
    assert rebuilt.read_bytes() == ARTIFACT.read_bytes()
    assert b"\r\n" not in rebuilt.read_bytes()
    subprocess.run(
        [sys.executable, str(SCRIPT), "--check"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )


def test_pure_de_sitter_mu_squared_equals_lambda_l1_for_d3_through_d6() -> None:
    rows = build_receipt()["pure_de_sitter_identity"]["checks"]
    assert [row["spacetime_dimension"] for row in rows] == [3, 4, 5, 6]
    assert [row["mu_squared_exact"] for row in rows] == ["1", "2", "3", "4"]
    assert all(row["mu_squared_exact"] == row["lambda_l1_exact"] for row in rows)
    assert all(row["kappa_times_r_c_exact"] == "1" for row in rows)


def test_entropy_maximum_and_area_hessian_have_distinct_statuses() -> None:
    capacity = build_receipt()["entropy_and_capacity"]
    entropy_maximum = capacity["finite_generalized_entropy_maximum"]
    area = capacity["area_analytic_relaxation"]

    assert entropy_maximum["identity"] == "S_gen(p,d)=log(M)-D_KL(p || d/M)"
    assert entropy_maximum["unique_simplex_maximum"] is True
    assert entropy_maximum["maximum_value"] == "log(M)"

    assert area["domain"] == "positive-real relaxation of integer sector dimensions d_i"
    assert area["logarithmic_area_coordinate"] == (
        "F(d)=sum_i d_i*log(d_i)/M"
    )
    assert area["gradient_component_exact"] == "1/(n*d)"
    assert area["gradient_component"] > 0.0
    assert area["unconstrained_stationary"] is False
    assert area["hessian_exact"] == "(1/(n*d^2))*(I-(2/n)*J)"
    assert area["fixed_M_tangent_eigenvalue"] > 0.0
    assert area["fixed_M_tangent_classification"] == "strict_local_minimum"
    assert area["homogeneous_eigenvalue"] < 0.0
    assert area["homogeneous_directional_second_derivative"] < 0.0


def test_area_gradient_and_hessian_match_finite_differences() -> None:
    n = 12
    d = 7.0
    point = np.full(n, d)

    def coordinate(dimensions: np.ndarray) -> float:
        return float(np.dot(dimensions, np.log(dimensions)) / np.sum(dimensions))

    epsilon = 1.0e-3
    basis0 = np.zeros(n)
    basis0[0] = 1.0
    numerical_gradient = (
        coordinate(point + epsilon * basis0)
        - coordinate(point - epsilon * basis0)
    ) / (2.0 * epsilon)
    assert numerical_gradient == pytest.approx(1.0 / (n * d), rel=2.0e-8)

    tangent = np.zeros(n)
    tangent[0] = 1.0
    tangent[1] = -1.0
    tangent_second = (
        coordinate(point + epsilon * tangent)
        - 2.0 * coordinate(point)
        + coordinate(point - epsilon * tangent)
    ) / epsilon**2
    assert tangent_second == pytest.approx(2.0 / (n * d**2), rel=1.0e-6)

    homogeneous = np.ones(n)
    homogeneous_second = (
        coordinate(point + epsilon * homogeneous)
        - 2.0 * coordinate(point)
        + coordinate(point - epsilon * homogeneous)
    ) / epsilon**2
    assert homogeneous_second == pytest.approx(-1.0 / d**2, rel=1.0e-6)


def test_one_sided_capacity_transfer_is_monotone_and_matches_budget_rows() -> None:
    transfer = build_receipt()["entropy_and_capacity"]["one_sided_capacity_transfer"]
    finite = transfer["finite_integer_statement"]
    analytic = transfer["positive_real_analytic_interpolation"]
    assert "positive integer" in finite["domain"]
    assert finite["f_equals_zero_boundary_maximum"] is True
    assert finite["monotone_on_admissible_depletions"] is True
    assert analytic["f_equals_zero_one_sided_boundary_maximum"] is True
    assert analytic["right_derivative_at_zero"] == -1.0
    assert analytic["monotone_coordinate_loss"] is True
    rows = finite["budget_checks"]
    assert [row["observer_capacity_fraction"] for row in rows] == [0.02, 0.05, 0.1]
    assert [row["coordinate_loss"] for row in rows] == [
        pytest.approx(0.020202707318),
        pytest.approx(0.051293294388),
        pytest.approx(0.105360515658),
    ]
    assert all(
        row["coordinate_change"] < 0.0 < row["coordinate_loss"] for row in rows
    )
    assert [row["initial_integer_sector_dimension"] for row in rows] == [100, 100, 100]
    assert [row["depleted_integer_sector_dimension"] for row in rows] == [98, 95, 90]


def test_port_edge_and_face_graph_spectra_are_frozen() -> None:
    raw = build_receipt()["graph_spectra"]["raw_laplacian_spectra"]
    assert [
        (row["raw_laplacian_eigenvalue_exact"], row["multiplicity"])
        for row in raw["ports"]
    ] == [
        ("0", 1),
        ("5 - sqrt(5)", 3),
        ("6", 5),
        ("5 + sqrt(5)", 3),
    ]
    assert [
        (row["raw_laplacian_eigenvalue_exact"], row["multiplicity"])
        for row in raw["edge_sectors"]
    ] == [
        ("0", 1),
        ("5 - sqrt(5)", 3),
        ("6", 5),
        ("5 + sqrt(5)", 3),
        ("10", 18),
    ]
    assert [
        (row["raw_laplacian_eigenvalue_exact"], row["multiplicity"])
        for row in raw["faces"]
    ] == [
        ("0", 1),
        ("3 - sqrt(5)", 3),
        ("2", 5),
        ("3", 4),
        ("5", 4),
        ("3 + sqrt(5)", 3),
    ]


def test_line_graph_theorem_is_checked_by_exact_integer_incidence_identities() -> None:
    theorem = build_receipt()["graph_spectra"]["line_graph_theorem"]
    assert theorem["port_identity_max_integer_residual"] == 0
    assert theorem["edge_identity_max_integer_residual"] == 0
    assert theorem["incidence_rank"] == 12
    assert theorem["extra_eigenvalue_exact"] == "2k=10"
    assert theorem["extra_multiplicity_exact"] == "m-n=18"
    assert theorem["port_edge_low_spectra_identical"] is True


def test_line_graph_negative_control_breaks_after_one_edge_is_deleted() -> None:
    graphs = _icosahedron_graphs()
    mutated = graphs["edge_adjacency"].copy()
    pair = np.argwhere(np.triu(mutated, k=1) == 1)[0]
    mutated[pair[0], pair[1]] = 0
    mutated[pair[1], pair[0]] = 0
    original_eigenvalues = np.linalg.eigvalsh(_laplacian(graphs["edge_adjacency"]))
    mutated_eigenvalues = np.linalg.eigvalsh(_laplacian(mutated))
    assert not np.allclose(original_eigenvalues, mutated_eigenvalues)
    assert not np.all(np.sum(mutated, axis=1) == 8)


def test_normalized_spectrum_is_conditional_on_both_named_premises() -> None:
    receipt = build_receipt()
    boundary = receipt["claim_boundary"]
    conditional = receipt["graph_spectra"][
        "conditional_normalized_port_shock_spectrum"
    ]
    assert boundary["assumption_tokens"] == ["DS-GAUGE", "DS-LAPLACIAN"]
    assert receipt["assumptions"]["DS-GAUGE"]["discharged"] is False
    assert receipt["assumptions"]["DS-LAPLACIAN"]["discharged"] is False
    assert boundary["physical_shock_spectrum_established"] is False
    assert boundary["promotion_allowed"] is False
    assert conditional["status"] == (
        "DERIVED_CONDITIONAL_ON_DS-GAUGE_AND_DS-LAPLACIAN"
    )
    assert conditional["scope_boundary"]["physical_prediction"] is False
    assert len(conditional["scope_boundary"]["escape_routes"]) == 3


def test_primary_source_and_physical_interpretation_boundary_are_explicit() -> None:
    receipt = build_receipt()
    source = receipt["primary_source"]
    interpretation = receipt["physical_interpretation"]
    assert source["arxiv_id"] == "2607.14042v1"
    assert "Eq. (C.9)" in " ".join(
        source["equation_map"]["higher_dimensional_shock_operator"]
    )
    assert interpretation["status"] == (
        "CONDITIONAL_INTERPRETATION_WITH_OPEN_PHYSICAL_PREMISES"
    )
    assert interpretation["premise_tokens"] == [
        "capacity_to_horizon_area",
        "capacity_ledger_to_observer_mass",
        "einstein_branch_realization",
        "coefficient_and_physical_scale",
    ]
    assert all(
        row["discharged"] is False
        for row in interpretation["premises"].values()
    )
    assert interpretation["physical_time_advance_explained"] is False
    assert interpretation["promotion_allowed"] is False
    assert len(interpretation["scope_boundary"]["escape_routes"]) == 4


def test_invariant_shock_mode_is_exactly_minus_two() -> None:
    conditional = build_receipt()["graph_spectra"][
        "conditional_normalized_port_shock_spectrum"
    ]
    spectrum = conditional["spectrum"]
    assert spectrum[0]["block"] == "1"
    assert spectrum[0]["shock_eigenvalue_exact"] == "-2"
    assert spectrum[0]["shock_eigenvalue"] == -2.0
    assert spectrum[1]["block"] == "3"
    assert spectrum[1]["shock_eigenvalue"] == 0.0
    assert conditional["invariant_minus_two_guard"]["passed"] is True


def test_repair_generator_domain_does_not_contain_sector_dimension_shock() -> None:
    boundary = build_receipt()["repair_generator_domain_boundary"]
    assert boundary["repair_domain"] == "K_r=L^2(X_r,pi_r) at fixed sector structure"
    assert boundary["shock_deformation_in_repair_domain"] is False
    assert boundary["repair_gap_constrains_shock_mode"] is False
    assert boundary["repair_gap_forbids_shock_mode"] is False
    assert boundary["repair_gap_protects_shock_mode"] is False
    assert boundary["repair_kernel_exactness_established_by_this_receipt"] is False
    assert boundary["shock_generator_supplied_by_repair_operator"] is False
    assert boundary["classification"] == "domain_separation_with_shock_generator_open"


@pytest.mark.parametrize(
    ("path", "value", "failure"),
    [
        (
            ("claim_boundary", "promotion_allowed"),
            True,
            "promotion_allowed_must_be_false",
        ),
        (
            ("claim_boundary", "physical_shock_spectrum_established"),
            True,
            "physical_shock_spectrum_established_must_be_false",
        ),
        (
            ("assumptions", "DS-GAUGE", "discharged"),
            True,
            "DS-GAUGE_must_remain_undischarged",
        ),
        (
            ("assumptions", "DS-LAPLACIAN", "discharged"),
            True,
            "DS-LAPLACIAN_must_remain_undischarged",
        ),
        (
            ("entropy_and_capacity", "area_analytic_relaxation", "unconstrained_stationary"),
            True,
            "symmetric_area_point_must_be_nonstationary",
        ),
        (
            (
                "entropy_and_capacity",
                "area_analytic_relaxation",
                "fixed_M_tangent_classification",
            ),
            "strict_local_maximum",
            "fixed_M_tangent_sign_mismatch",
        ),
        (
            (
                "graph_spectra",
                "conditional_normalized_port_shock_spectrum",
                "invariant_minus_two_guard",
                "observed",
            ),
            -2.76,
            "minus_two_observed_mismatch",
        ),
        (
            (
                "repair_generator_domain_boundary",
                "shock_deformation_in_repair_domain",
            ),
            True,
            "shock_deformation_in_repair_domain_must_be_false",
        ),
        (
            (
                "repair_generator_domain_boundary",
                "repair_kernel_exactness_established_by_this_receipt",
            ),
            True,
            "repair_kernel_exactness_established_by_this_receipt_must_be_false",
        ),
        (
            ("physical_interpretation", "physical_time_advance_explained"),
            True,
            "physical_interpretation_physical_time_advance_explained_must_be_false",
        ),
        (
            (
                "physical_interpretation",
                "premises",
                "capacity_to_horizon_area",
                "discharged",
            ),
            True,
            "capacity_to_horizon_area_must_remain_undischarged",
        ),
    ],
)
def test_adversarial_status_and_sign_mutations_fail_closed(
    path: tuple[str, ...],
    value: object,
    failure: str,
) -> None:
    mutated = copy.deepcopy(build_receipt())
    cursor = mutated
    for key in path[:-1]:
        cursor = cursor[key]
    cursor[path[-1]] = value
    assert failure in receipt_failures(mutated)
    with pytest.raises(ReceiptValidationError, match=failure):
        validate_receipt(mutated)


def test_adversarial_capacity_loss_and_line_graph_mutations_fail_closed() -> None:
    mutated_loss = copy.deepcopy(build_receipt())
    mutated_loss["entropy_and_capacity"]["one_sided_capacity_transfer"][
        "finite_integer_statement"
    ][
        "budget_checks"
    ][0]["coordinate_loss"] *= -1.0
    assert "capacity_budget_losses_mismatch" in receipt_failures(mutated_loss)

    mutated_line = copy.deepcopy(build_receipt())
    mutated_line["graph_spectra"]["line_graph_theorem"][
        "edge_identity_max_integer_residual"
    ] = 1
    assert (
        "edge_identity_max_integer_residual_must_be_zero"
        in receipt_failures(mutated_line)
    )
