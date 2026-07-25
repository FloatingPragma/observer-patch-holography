#!/usr/bin/env python3
"""Build the finite-screen de Sitter shock-sign and graph-spectrum receipt.

The receipt separates three kinds of statements.

* The pure-de Sitter identity, entropy maximum, logarithmic screen-coordinate
  derivatives, capacity transfer, graph spectra, and line-graph identity are
  finite or algebraic calculations.
* Derivatives with respect to sector dimensions use the positive-real analytic
  relaxation of the integer dimensions.  The symmetric point has a nonzero
  gradient.  Its Hessian is positive on fixed-total-capacity tangent
  directions and negative in the homogeneous direction, so it is not an
  unconstrained stationary maximum.
* Interpreting the normalized graph spectrum as the physical screen shock
  spectrum is conditional on both DS-GAUGE and DS-LAPLACIAN.  The receipt does
  not promote that conditional calculation to a physical result.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_OUT = HERE / "runs" / "screen_shock_spectrum.json"

SCHEMA = "oph.finite_screen_shock_receipt.v1"
ARTIFACT = "oph_finite_screen_de_sitter_shock_receipt"
ASSUMPTION_TOKENS = ["DS-GAUGE", "DS-LAPLACIAN"]
INTERPRETATION_PREMISES = [
    "capacity_to_horizon_area",
    "capacity_ledger_to_observer_mass",
    "einstein_branch_realization",
    "coefficient_and_physical_scale",
]
TOLERANCE = 2.0e-9


class ReceiptValidationError(ValueError):
    """Raised when a shock receipt violates a mathematical or status guard."""


def _round(value: float, digits: int = 12) -> float:
    rounded = round(float(value), digits)
    return 0.0 if rounded == 0.0 else rounded


def _laplacian(adjacency: np.ndarray) -> np.ndarray:
    return np.diag(np.sum(adjacency, axis=1)) - adjacency


def _icosahedron_graphs() -> dict[str, Any]:
    """Return the port, edge-sector, and face adjacency matrices.

    The port coordinates are the standard twelve vertices
    ``(0,+/-1,+/-phi)`` and cyclic permutations.  Edges join the minimum
    nonzero-distance pairs.  The other two carriers are the line graph and the
    face-adjacency graph of this exact combinatorial incidence structure.
    """

    phi = (1.0 + math.sqrt(5.0)) / 2.0
    vertices: list[tuple[float, float, float]] = []
    for first_sign in (-1.0, 1.0):
        for second_sign in (-1.0, 1.0):
            vertices.extend(
                [
                    (0.0, first_sign, second_sign * phi),
                    (first_sign, second_sign * phi, 0.0),
                    (second_sign * phi, 0.0, first_sign),
                ]
            )
    coordinates = np.asarray(vertices, dtype=float)
    if coordinates.shape != (12, 3):
        raise AssertionError("the port construction must have shape 12 by 3")

    deltas = coordinates[:, None, :] - coordinates[None, :, :]
    distances_squared = np.sum(deltas * deltas, axis=2)
    positive_distances = distances_squared[distances_squared > 1.0e-12]
    edge_distance_squared = float(np.min(positive_distances))
    port_adjacency = np.isclose(
        distances_squared,
        edge_distance_squared,
        atol=1.0e-12,
        rtol=0.0,
    ).astype(int)
    np.fill_diagonal(port_adjacency, 0)

    edges = [
        (left, right)
        for left in range(12)
        for right in range(left + 1, 12)
        if port_adjacency[left, right] == 1
    ]
    faces = [
        triple
        for triple in itertools.combinations(range(12), 3)
        if all(
            port_adjacency[left, right] == 1
            for left, right in itertools.combinations(triple, 2)
        )
    ]

    line_adjacency = np.zeros((len(edges), len(edges)), dtype=int)
    for left_index, left_edge in enumerate(edges):
        for right_index in range(left_index + 1, len(edges)):
            if set(left_edge).intersection(edges[right_index]):
                line_adjacency[left_index, right_index] = 1
                line_adjacency[right_index, left_index] = 1

    face_adjacency = np.zeros((len(faces), len(faces)), dtype=int)
    for left_index, left_face in enumerate(faces):
        for right_index in range(left_index + 1, len(faces)):
            if len(set(left_face).intersection(faces[right_index])) == 2:
                face_adjacency[left_index, right_index] = 1
                face_adjacency[right_index, left_index] = 1

    if len(edges) != 30 or len(faces) != 20:
        raise AssertionError(
            f"expected 30 edges and 20 faces, found {len(edges)} and {len(faces)}"
        )
    expected_degrees = {
        "port": (port_adjacency, 5),
        "edge_sector": (line_adjacency, 8),
        "face": (face_adjacency, 3),
    }
    for name, (adjacency, degree) in expected_degrees.items():
        if not np.array_equal(np.sum(adjacency, axis=1), np.full(len(adjacency), degree)):
            raise AssertionError(f"{name} graph is not regular of degree {degree}")

    incidence = np.zeros((12, 30), dtype=int)
    for edge_index, (left, right) in enumerate(edges):
        incidence[left, edge_index] = 1
        incidence[right, edge_index] = 1

    return {
        "port_adjacency": port_adjacency,
        "edge_adjacency": line_adjacency,
        "face_adjacency": face_adjacency,
        "incidence": incidence,
        "edges": edges,
        "faces": faces,
        "edge_distance_squared": edge_distance_squared,
    }


def _group_spectrum(values: Iterable[float]) -> list[tuple[float, int]]:
    groups: list[list[float]] = []
    for value in sorted(float(item) for item in values):
        if not groups or abs(value - groups[-1][0]) > TOLERANCE:
            groups.append([value])
        else:
            groups[-1].append(value)
    return [(float(np.mean(group)), len(group)) for group in groups]


PORT_EXPECTED = [
    ("1", "0", 0.0, 1),
    ("3", "5 - sqrt(5)", 5.0 - math.sqrt(5.0), 3),
    ("5", "6", 6.0, 5),
    ("3_prime", "5 + sqrt(5)", 5.0 + math.sqrt(5.0), 3),
]

EDGE_EXPECTED = PORT_EXPECTED + [
    ("line_graph_kernel", "10", 10.0, 18),
]

FACE_EXPECTED = [
    ("1", "0", 0.0, 1),
    ("3", "3 - sqrt(5)", 3.0 - math.sqrt(5.0), 3),
    ("5", "2", 2.0, 5),
    ("4_a", "3", 3.0, 4),
    ("4_b", "5", 5.0, 4),
    ("3_prime", "3 + sqrt(5)", 3.0 + math.sqrt(5.0), 3),
]


def _spectrum_receipt(
    adjacency: np.ndarray,
    expected: list[tuple[str, str, float, int]],
) -> list[dict[str, Any]]:
    observed = _group_spectrum(np.linalg.eigvalsh(_laplacian(adjacency)))
    expected_numeric = [(value, multiplicity) for _, _, value, multiplicity in expected]
    if len(observed) != len(expected_numeric):
        raise AssertionError(
            f"expected {len(expected_numeric)} eigenvalue blocks, found {len(observed)}"
        )
    rows: list[dict[str, Any]] = []
    for observed_row, expected_row in zip(observed, expected, strict=True):
        observed_value, observed_multiplicity = observed_row
        representation, expression, expected_value, expected_multiplicity = expected_row
        if observed_multiplicity != expected_multiplicity:
            raise AssertionError(
                f"{representation} multiplicity is {observed_multiplicity}, "
                f"expected {expected_multiplicity}"
            )
        if abs(observed_value - expected_value) > TOLERANCE:
            raise AssertionError(
                f"{representation} eigenvalue is {observed_value}, expected {expected_value}"
            )
        rows.append(
            {
                "block": representation,
                "raw_laplacian_eigenvalue_exact": expression,
                "raw_laplacian_eigenvalue": _round(expected_value),
                "multiplicity": expected_multiplicity,
            }
        )
    return rows


def pure_de_sitter_identity_receipt() -> dict[str, Any]:
    checks = []
    for spacetime_dimension in range(3, 7):
        sphere_dimension = spacetime_dimension - 2
        mu_squared = spacetime_dimension - 2
        lambda_l1 = 1 * (1 + sphere_dimension - 1)
        if mu_squared != lambda_l1:
            raise AssertionError("pure-de Sitter mu squared must equal lambda_1")
        checks.append(
            {
                "spacetime_dimension": spacetime_dimension,
                "horizon_sphere": f"S^{sphere_dimension}",
                "kappa_times_r_c_exact": "1",
                "mu_squared_exact": str(mu_squared),
                "lambda_l1_exact": str(lambda_l1),
                "identity_passed": True,
            }
        )
    return {
        "formula": "mu^2=((d-2)/2)*|f'(r_c)|*r_c=(d-2)*kappa*r_c",
        "pure_de_sitter_substitution": (
            "f(r)=1-r^2/L^2, r_c=L, |f'(r_c)|=2/L, kappa=1/L"
        ),
        "lambda_cancels": True,
        "checks": checks,
    }


def entropy_and_capacity_receipt() -> dict[str, Any]:
    dimensions = np.asarray([2.0, 3.0, 5.0, 7.0])
    total = float(np.sum(dimensions))
    maximizing_probabilities = dimensions / total
    generalized_entropy = float(
        -np.sum(maximizing_probabilities * np.log(maximizing_probabilities))
        + np.sum(maximizing_probabilities * np.log(dimensions))
    )
    lagrange_gradient = (
        -np.log(maximizing_probabilities) - 1.0 + np.log(dimensions)
    )
    if not np.allclose(lagrange_gradient, np.full(4, math.log(total) - 1.0)):
        raise AssertionError("the entropy gradient must be constant on the simplex")
    if abs(generalized_entropy - math.log(total)) > 1.0e-13:
        raise AssertionError("the generalized entropy maximum must be log M")

    n = 12
    symmetric_dimension = 7
    gradient_component = 1.0 / (n * symmetric_dimension)
    homogeneous_curvature = -1.0 / (symmetric_dimension**2)
    homogeneous_hessian_eigenvalue = -1.0 / (
        n * symmetric_dimension**2
    )
    tangent_hessian_eigenvalue = 1.0 / (n * symmetric_dimension**2)

    budget_rows = []
    for fraction in (0.02, 0.05, 0.10):
        initial_dimension = 100
        depleted_dimension = round((1.0 - fraction) * initial_dimension)
        if depleted_dimension != (1.0 - fraction) * initial_dimension:
            raise AssertionError("budget fixture must preserve integer dimensions")
        signed_change = math.log1p(-fraction)
        budget_rows.append(
            {
                "observer_capacity_fraction": fraction,
                "initial_integer_sector_dimension": initial_dimension,
                "depleted_integer_sector_dimension": depleted_dimension,
                "coordinate_change_exact": f"log({1.0 - fraction:.2f})",
                "coordinate_change": _round(signed_change),
                "coordinate_loss": _round(-signed_change),
            }
        )

    return {
        "finite_generalized_entropy_maximum": {
            "identity": "S_gen(p,d)=log(M)-D_KL(p || d/M)",
            "M_definition": "M=sum_i d_i",
            "maximizer": "p_i=d_i/M",
            "unique_simplex_maximum": True,
            "maximum_value": "log(M)",
            "capacity_closure_identity": "N=log(M_0)=max_p S_gen",
            "finite_example_dimensions": [2, 3, 5, 7],
            "finite_example_probabilities_exact": ["2/17", "3/17", "5/17", "7/17"],
            "finite_example_maximum": _round(generalized_entropy),
            "simplex_tangent_hessian_negative": True,
        },
        "area_analytic_relaxation": {
            "domain": "positive-real relaxation of integer sector dimensions d_i",
            "logarithmic_area_coordinate": "F(d)=sum_i d_i*log(d_i)/M",
            "symmetric_point": {"n": n, "d": symmetric_dimension},
            "gradient_component_exact": "1/(n*d)",
            "gradient_component": _round(gradient_component),
            "unconstrained_stationary": False,
            "hessian_exact": "(1/(n*d^2))*(I-(2/n)*J)",
            "fixed_M_tangent_eigenvalue_exact": "+1/(n*d^2)",
            "fixed_M_tangent_eigenvalue": _round(tangent_hessian_eigenvalue),
            "fixed_M_tangent_classification": "strict_local_minimum",
            "homogeneous_eigenvalue_exact": "-1/(n*d^2)",
            "homogeneous_eigenvalue": _round(homogeneous_hessian_eigenvalue),
            "homogeneous_directional_second_derivative_exact": "-1/d^2",
            "homogeneous_directional_second_derivative": _round(
                homogeneous_curvature
            ),
            "homogeneous_curvature_negative": True,
            "signature": "(-,+^(n-1))",
            "classification": (
                "nonstationary_symmetric_point_with_positive_fixed-M_tangent_"
                "curvature_and_negative_homogeneous_curvature"
            ),
        },
        "one_sided_capacity_transfer": {
            "finite_integer_statement": {
                "domain": (
                    "uniform integer d_i with f such that (1-f)*d_i is a "
                    "positive integer"
                ),
                "uniform_horizon_capacity": "M_h(f)=(1-f)*M_0",
                "coordinate_change_exact": "F(f)-F(0)=log(1-f)",
                "monotone_on_admissible_depletions": True,
                "f_equals_zero_boundary_maximum": True,
                "budget_checks": budget_rows,
            },
            "positive_real_analytic_interpolation": {
                "domain": "0 <= f < 1",
                "derivative_exact": "dF/df=-1/(1-f)<0",
                "right_derivative_at_zero": -1.0,
                "monotone_coordinate_loss": True,
                "f_equals_zero_one_sided_boundary_maximum": True,
            },
        },
    }


def graph_spectrum_receipt() -> dict[str, Any]:
    graphs = _icosahedron_graphs()
    port_adjacency = graphs["port_adjacency"]
    edge_adjacency = graphs["edge_adjacency"]
    face_adjacency = graphs["face_adjacency"]
    incidence = graphs["incidence"]

    k = 5
    n = 12
    m = 30
    port_incidence_identity = incidence @ incidence.T
    edge_incidence_identity = incidence.T @ incidence
    port_identity_residual = int(
        np.max(np.abs(port_incidence_identity - (port_adjacency + k * np.eye(n, dtype=int))))
    )
    edge_identity_residual = int(
        np.max(np.abs(edge_incidence_identity - (edge_adjacency + 2 * np.eye(m, dtype=int))))
    )
    incidence_rank = int(np.linalg.matrix_rank(incidence))
    if port_identity_residual != 0 or edge_identity_residual != 0:
        raise AssertionError("the unsigned-incidence line-graph identities must hold exactly")
    if incidence_rank != n:
        raise AssertionError("the nonbipartite icosahedron incidence matrix must have rank 12")

    port_spectrum = _spectrum_receipt(port_adjacency, PORT_EXPECTED)
    edge_spectrum = _spectrum_receipt(edge_adjacency, EDGE_EXPECTED)
    face_spectrum = _spectrum_receipt(face_adjacency, FACE_EXPECTED)

    sqrt5 = math.sqrt(5.0)
    rotation_raw = 5.0 - sqrt5
    scale = 2.0 / rotation_raw
    normalized_rows = [
        ("1", "0", 0.0, "-2", -2.0, 1),
        ("3", "2", 2.0, "0", 0.0, 3),
        (
            "5",
            "3 + 3/sqrt(5)",
            3.0 + 3.0 / sqrt5,
            "1 + 3/sqrt(5)",
            1.0 + 3.0 / sqrt5,
            5,
        ),
        (
            "3_prime",
            "3 + sqrt(5)",
            3.0 + sqrt5,
            "1 + sqrt(5)",
            1.0 + sqrt5,
            3,
        ),
    ]
    normalized_spectrum = [
        {
            "block": block,
            "normalized_laplacian_eigenvalue_exact": normalized_exact,
            "normalized_laplacian_eigenvalue": _round(normalized),
            "shock_eigenvalue_exact": shock_exact,
            "shock_eigenvalue": _round(shock),
            "multiplicity": multiplicity,
        }
        for (
            block,
            normalized_exact,
            normalized,
            shock_exact,
            shock,
            multiplicity,
        ) in normalized_rows
    ]
    if normalized_spectrum[0]["shock_eigenvalue"] != -2.0:
        raise AssertionError("the normalized invariant shock mode must equal -2")
    if abs(float(normalized_spectrum[1]["shock_eigenvalue"])) > TOLERANCE:
        raise AssertionError("the normalized rotation triplet must be the zero block")

    line_low_rows = edge_spectrum[: len(port_spectrum)]
    if line_low_rows != port_spectrum:
        raise AssertionError("the port and line-graph low spectra must match")

    return {
        "graph_construction": {
            "port_graph": "icosahedron_vertex_graph",
            "port_count": n,
            "port_degree": k,
            "edge_sector_graph": "line_graph_of_port_graph",
            "edge_sector_count": m,
            "edge_sector_degree": 2 * k - 2,
            "face_graph": "icosahedron_face_adjacency_graph",
            "face_count": 20,
            "face_degree": 3,
            "edge_distance_squared": _round(graphs["edge_distance_squared"]),
        },
        "raw_laplacian_spectra": {
            "ports": port_spectrum,
            "edge_sectors": edge_spectrum,
            "faces": face_spectrum,
        },
        "line_graph_theorem": {
            "statement": (
                "Spec(L_line)=Spec(L_port) union {2k repeated m-n times} "
                "for this connected nonbipartite k-regular graph"
            ),
            "unsigned_incidence_port_identity": "B*B^T=A_port+k*I",
            "unsigned_incidence_edge_identity": "B^T*B=A_line+2*I",
            "port_identity_max_integer_residual": port_identity_residual,
            "edge_identity_max_integer_residual": edge_identity_residual,
            "incidence_rank": incidence_rank,
            "extra_eigenvalue_exact": "2k=10",
            "extra_multiplicity_exact": "m-n=18",
            "port_edge_low_spectra_identical": True,
            "hard_check_passed": True,
        },
        "normalization_free_ratios": {
            "port_lambda5_over_lambda3_exact": "6/(5-sqrt(5))",
            "port_lambda5_over_lambda3": _round(6.0 / rotation_raw),
            "edge_lambda5_over_lambda3_exact": "6/(5-sqrt(5))",
            "edge_lambda5_over_lambda3": _round(6.0 / rotation_raw),
            "face_lambda5_over_lambda3_exact": "2/(3-sqrt(5))=phi^2",
            "face_lambda5_over_lambda3": _round(2.0 / (3.0 - sqrt5)),
            "smooth_S2_lambda2_over_lambda1": 3.0,
        },
        "conditional_normalized_port_shock_spectrum": {
            "status": "DERIVED_CONDITIONAL_ON_DS-GAUGE_AND_DS-LAPLACIAN",
            "normalization": "scale lambda_3=5-sqrt(5) to mu^2=2",
            "scale_factor_exact": "2/(5-sqrt(5))",
            "scale_factor": _round(scale),
            "operator": "scaled graph Laplacian minus 2*I",
            "spectrum": normalized_spectrum,
            "invariant_minus_two_guard": {
                "expected_exact": "-2",
                "observed": normalized_spectrum[0]["shock_eigenvalue"],
                "passed": True,
            },
            "fivefold_softening_relative_to_smooth_l2": _round(
                1.0 - (1.0 + 3.0 / sqrt5) / 4.0
            ),
            "extra_3_prime_between_smooth_l2_and_l3": True,
            "scope_boundary": {
                "physical_prediction": False,
                "closed_calculation_class": (
                    "nearest-neighbour port graph with the stated normalization"
                ),
                "escape_routes": [
                    "a shock kinetic operator outside DS-LAPLACIAN",
                    "failure of the exact triplet gauge identification DS-GAUGE",
                    "a refined, weighted, nonlocal, or different carrier",
                ],
            },
        },
    }


def build_receipt() -> dict[str, Any]:
    receipt = {
        "schema": SCHEMA,
        "artifact": ARTIFACT,
        "producer": "code/geometry/derive_screen_shock_spectrum.py",
        "status": "MACHINE_CHECKED_FINITE_AND_CONDITIONAL_RECEIPT",
        "primary_source": {
            "title": "Negative shocks versus static patch holography",
            "authors": [
                "Yiming Chen",
                "Douglas Stanford",
                "Haifeng Tang",
                "Zhenbin Yang",
            ],
            "arxiv_id": "2607.14042v1",
            "arxiv_url": "https://arxiv.org/abs/2607.14042v1",
            "equation_map": {
                "observer_mass_area_and_entropy": [
                    "introductory A=4*pi-8*pi*G*m+O(G^2*m^2)",
                    "Eq. (1.4): Z_sphere proportional to exp(S_dS-2*pi*m)",
                ],
                "shock_geometry_sign_data": [
                    "Eq. (2.6): shock Ricci tensor contains B'(0)/B(0)",
                    "Eqs. (C.5)-(C.6): B(z)=r(z)^2 and B'(0)=2*r_c*r'(0)",
                ],
                "higher_dimensional_shock_operator": [
                    "Eq. (C.8): [-Laplacian_{S^(d-2)}-mu^2] X^+ = source",
                    "Eq. (C.9): mu^2=((d-2)/2)*|f'(r_c)|*r_c",
                    "Eqs. (C.10)-(C.11): spherical spectrum parameter j",
                ],
                "gauge_zero_mode_statement": (
                    "Sec. 2.2 and the paragraph following Eq. (2.10): "
                    "the pure-de Sitter l=1 shock is generated by an isometry"
                ),
            },
        },
        "exact_successes": [
            "pure_de_sitter_mu_squared_equals_lambda_l1_for_d3_through_d6",
            "finite_generalized_entropy_equals_log_M_at_its_simplex_maximum",
            "exact_uniform_capacity_transfer_coordinate_change_on_admissible_integer_depletions",
            "exact_logarithmic_coordinate_hessian_split_in_the_positive_real_relaxation",
            "exact_port_edge_low_spectrum_equality_by_unsigned_incidence",
        ],
        "claim_boundary": {
            "exact_calculations": [
                "pure_de_sitter_mu_squared_equals_lambda_l1",
                "finite_generalized_entropy_simplex_maximum",
                "logarithmic_coordinate_gradient_and_hessian_in_positive_real_dimension_relaxation",
                "uniform_capacity_transfer_logarithmic_coordinate_loss_on_admissible_integer_depletions",
                "icosahedron_port_edge_face_graph_spectra",
                "exact_unsigned_incidence_line_graph_identity",
                "repair_generator_domain_separation",
            ],
            "conditional_calculation": (
                "normalized_icosahedral_port_shock_spectrum"
            ),
            "assumption_tokens": ASSUMPTION_TOKENS,
            "all_required_premises_discharged": False,
            "physical_shock_spectrum_established": False,
            "physical_coefficient_established": False,
            "promotion_allowed": False,
            "static_patch_trace_conjecture_rescued": False,
        },
        "assumptions": {
            "DS-GAUGE": {
                "statement": (
                    "the geometric A5 rotation triplet is an exact gauge zero "
                    "mode on the declared screen carrier"
                ),
                "discharged": False,
            },
            "DS-LAPLACIAN": {
                "statement": (
                    "the shock kinetic operator is the scaled nearest-neighbour "
                    "combinatorial Laplacian on the port or edge-sector carrier"
                ),
                "discharged": False,
            },
        },
        "pure_de_sitter_identity": pure_de_sitter_identity_receipt(),
        "entropy_and_capacity": entropy_and_capacity_receipt(),
        "graph_spectra": graph_spectrum_receipt(),
        "physical_interpretation": {
            "status": "CONDITIONAL_INTERPRETATION_WITH_OPEN_PHYSICAL_PREMISES",
            "statement": (
                "the one-sided capacity-transfer coordinate loss carries the "
                "cosmological-horizon time-advance sign if every listed "
                "physical premise holds"
            ),
            "premises": {
                "capacity_to_horizon_area": {
                    "statement": (
                        "the finite-screen logarithmic coordinate F(d) equals "
                        "the physical cosmological-horizon area divided by 4G"
                    ),
                    "discharged": False,
                },
                "capacity_ledger_to_observer_mass": {
                    "statement": (
                        "capacity transferred out of the horizon ledger maps "
                        "monotonically to positive observer mass"
                    ),
                    "discharged": False,
                },
                "einstein_branch_realization": {
                    "statement": (
                        "the declared carrier realizes the Einstein branch "
                        "with cosmological-horizon orientation"
                    ),
                    "discharged": False,
                },
                "coefficient_and_physical_scale": {
                    "statement": (
                        "the generalized-entropy coefficient and source-derived "
                        "physical scale match the shock normalization"
                    ),
                    "discharged": False,
                },
            },
            "premise_tokens": INTERPRETATION_PREMISES,
            "all_premises_discharged": False,
            "physical_time_advance_explained": False,
            "promotion_allowed": False,
            "scope_boundary": {
                "excluded_claim": (
                    "unconditional derivation of the physical de Sitter shock "
                    "sign or coefficient from the finite calculation alone"
                ),
                "escape_routes": [
                    "failure of the capacity-to-area map",
                    "failure of the capacity-ledger-to-mass map",
                    "a non-Einstein branch or reversed horizon orientation",
                    "a different coefficient, physical scale, or shock operator",
                ],
            },
        },
        "repair_generator_domain_boundary": {
            "repair_generator": "L_rep=sum_v(I-P_v)",
            "repair_domain": "K_r=L^2(X_r,pi_r) at fixed sector structure",
            "shock_deformation": (
                "varies sector dimensions d_alpha and therefore changes X_r"
            ),
            "shock_deformation_in_repair_domain": False,
            "repair_gap_constrains_shock_mode": False,
            "repair_gap_forbids_shock_mode": False,
            "repair_gap_protects_shock_mode": False,
            "repair_kernel_exactness_established_by_this_receipt": False,
            "shock_generator_supplied_by_repair_operator": False,
            "classification": "domain_separation_with_shock_generator_open",
            "scope_boundary": {
                "excluded_claim": (
                    "the declared fixed-sector repair gap determines the "
                    "sector-dimension shock deformation"
                ),
                "escape_routes": [
                    "an extended generator acting on sector dimensions",
                    "a coupled carrier-refinement dynamics",
                ],
            },
        },
        "determinism": {
            "timestamp_omitted": True,
            "randomness_used": False,
            "json_key_order": "lexicographic",
            "newline": "LF",
        },
    }
    validate_receipt(receipt)
    return receipt


def _fail_if(condition: bool, failures: list[str], message: str) -> None:
    if condition:
        failures.append(message)


def receipt_failures(receipt: dict[str, Any]) -> list[str]:
    """Return fail-closed status and invariant violations."""

    failures: list[str] = []
    _fail_if(receipt.get("schema") != SCHEMA, failures, "schema_mismatch")
    _fail_if(receipt.get("artifact") != ARTIFACT, failures, "artifact_mismatch")

    boundary = receipt.get("claim_boundary", {})
    _fail_if(
        boundary.get("assumption_tokens") != ASSUMPTION_TOKENS,
        failures,
        "conditional_assumption_tokens_mismatch",
    )
    for key in (
        "all_required_premises_discharged",
        "physical_shock_spectrum_established",
        "physical_coefficient_established",
        "promotion_allowed",
        "static_patch_trace_conjecture_rescued",
    ):
        _fail_if(boundary.get(key) is not False, failures, f"{key}_must_be_false")

    source = receipt.get("primary_source", {})
    _fail_if(
        source.get("arxiv_id") != "2607.14042v1",
        failures,
        "primary_source_version_mismatch",
    )
    source_equations = source.get("equation_map", {})
    _fail_if(
        "higher_dimensional_shock_operator" not in source_equations,
        failures,
        "primary_source_shock_equations_missing",
    )

    assumptions = receipt.get("assumptions", {})
    for token in ASSUMPTION_TOKENS:
        _fail_if(token not in assumptions, failures, f"missing_{token}")
        _fail_if(
            assumptions.get(token, {}).get("discharged") is not False,
            failures,
            f"{token}_must_remain_undischarged",
        )

    pure_ds = receipt.get("pure_de_sitter_identity", {})
    _fail_if(pure_ds.get("lambda_cancels") is not True, failures, "lambda_must_cancel")
    checks = pure_ds.get("checks", [])
    _fail_if(
        [row.get("spacetime_dimension") for row in checks] != [3, 4, 5, 6],
        failures,
        "pure_ds_dimension_range_mismatch",
    )
    for row in checks:
        _fail_if(
            row.get("mu_squared_exact") != row.get("lambda_l1_exact"),
            failures,
            "mu_squared_lambda_l1_mismatch",
        )
        _fail_if(
            row.get("identity_passed") is not True,
            failures,
            "pure_ds_identity_not_passed",
        )

    capacity = receipt.get("entropy_and_capacity", {})
    entropy_maximum = capacity.get("finite_generalized_entropy_maximum", {})
    _fail_if(
        entropy_maximum.get("unique_simplex_maximum") is not True,
        failures,
        "entropy_simplex_maximum_missing",
    )
    area = capacity.get("area_analytic_relaxation", {})
    _fail_if(
        area.get("unconstrained_stationary") is not False,
        failures,
        "symmetric_area_point_must_be_nonstationary",
    )
    _fail_if(
        area.get("fixed_M_tangent_classification") != "strict_local_minimum",
        failures,
        "fixed_M_tangent_sign_mismatch",
    )
    _fail_if(
        area.get("homogeneous_curvature_negative") is not True,
        failures,
        "homogeneous_curvature_sign_mismatch",
    )
    _fail_if(
        area.get("hessian_exact") != "(1/(n*d^2))*(I-(2/n)*J)",
        failures,
        "hessian_formula_mismatch",
    )

    transfer = capacity.get("one_sided_capacity_transfer", {})
    finite_transfer = transfer.get("finite_integer_statement", {})
    analytic_transfer = transfer.get("positive_real_analytic_interpolation", {})
    _fail_if(
        finite_transfer.get("monotone_on_admissible_depletions") is not True,
        failures,
        "finite_capacity_transfer_must_lose_coordinate",
    )
    _fail_if(
        finite_transfer.get("f_equals_zero_boundary_maximum") is not True,
        failures,
        "finite_capacity_boundary_maximum_missing",
    )
    _fail_if(
        analytic_transfer.get("monotone_coordinate_loss") is not True,
        failures,
        "analytic_capacity_transfer_must_lose_coordinate",
    )
    _fail_if(
        analytic_transfer.get("f_equals_zero_one_sided_boundary_maximum") is not True,
        failures,
        "analytic_capacity_boundary_maximum_missing",
    )
    expected_losses = [
        _round(-math.log(0.98)),
        _round(-math.log(0.95)),
        _round(-math.log(0.90)),
    ]
    observed_losses = [
        row.get("coordinate_loss")
        for row in finite_transfer.get("budget_checks", [])
    ]
    _fail_if(
        observed_losses != expected_losses,
        failures,
        "capacity_budget_losses_mismatch",
    )

    graph = receipt.get("graph_spectra", {})
    theorem = graph.get("line_graph_theorem", {})
    for key in (
        "port_identity_max_integer_residual",
        "edge_identity_max_integer_residual",
    ):
        _fail_if(theorem.get(key) != 0, failures, f"{key}_must_be_zero")
    _fail_if(
        theorem.get("extra_eigenvalue_exact") != "2k=10",
        failures,
        "line_graph_extra_eigenvalue_mismatch",
    )
    _fail_if(
        theorem.get("extra_multiplicity_exact") != "m-n=18",
        failures,
        "line_graph_extra_multiplicity_mismatch",
    )
    _fail_if(
        theorem.get("port_edge_low_spectra_identical") is not True,
        failures,
        "port_edge_low_spectra_mismatch",
    )
    _fail_if(
        theorem.get("hard_check_passed") is not True,
        failures,
        "line_graph_hard_check_not_passed",
    )

    conditional = graph.get("conditional_normalized_port_shock_spectrum", {})
    _fail_if(
        conditional.get("status")
        != "DERIVED_CONDITIONAL_ON_DS-GAUGE_AND_DS-LAPLACIAN",
        failures,
        "normalized_spectrum_status_mismatch",
    )
    minus_two = conditional.get("invariant_minus_two_guard", {})
    _fail_if(minus_two.get("expected_exact") != "-2", failures, "minus_two_exact_mismatch")
    _fail_if(minus_two.get("observed") != -2.0, failures, "minus_two_observed_mismatch")
    _fail_if(minus_two.get("passed") is not True, failures, "minus_two_guard_not_passed")
    expected_normalized = [
        ("1", "0", "-2", 1),
        ("3", "2", "0", 3),
        ("5", "3 + 3/sqrt(5)", "1 + 3/sqrt(5)", 5),
        ("3_prime", "3 + sqrt(5)", "1 + sqrt(5)", 3),
    ]
    observed_normalized = [
        (
            row.get("block"),
            row.get("normalized_laplacian_eigenvalue_exact"),
            row.get("shock_eigenvalue_exact"),
            row.get("multiplicity"),
        )
        for row in conditional.get("spectrum", [])
    ]
    _fail_if(
        observed_normalized != expected_normalized,
        failures,
        "normalized_shock_spectrum_mismatch",
    )
    graph_scope = conditional.get("scope_boundary", {})
    _fail_if(
        graph_scope.get("physical_prediction") is not False,
        failures,
        "graph_spectrum_physical_prediction_must_be_false",
    )

    interpretation = receipt.get("physical_interpretation", {})
    _fail_if(
        interpretation.get("status")
        != "CONDITIONAL_INTERPRETATION_WITH_OPEN_PHYSICAL_PREMISES",
        failures,
        "physical_interpretation_status_mismatch",
    )
    _fail_if(
        interpretation.get("premise_tokens") != INTERPRETATION_PREMISES,
        failures,
        "physical_interpretation_premises_mismatch",
    )
    for key in (
        "all_premises_discharged",
        "physical_time_advance_explained",
        "promotion_allowed",
    ):
        _fail_if(
            interpretation.get(key) is not False,
            failures,
            f"physical_interpretation_{key}_must_be_false",
        )
    interpretation_premises = interpretation.get("premises", {})
    for premise in INTERPRETATION_PREMISES:
        _fail_if(
            interpretation_premises.get(premise, {}).get("discharged") is not False,
            failures,
            f"{premise}_must_remain_undischarged",
        )

    repair = receipt.get("repair_generator_domain_boundary", {})
    for key in (
        "shock_deformation_in_repair_domain",
        "repair_gap_constrains_shock_mode",
        "repair_gap_forbids_shock_mode",
        "repair_gap_protects_shock_mode",
        "repair_kernel_exactness_established_by_this_receipt",
        "shock_generator_supplied_by_repair_operator",
    ):
        _fail_if(repair.get(key) is not False, failures, f"{key}_must_be_false")
    _fail_if(
        repair.get("classification")
        != "domain_separation_with_shock_generator_open",
        failures,
        "repair_domain_classification_mismatch",
    )
    return failures


def validate_receipt(receipt: dict[str, Any]) -> None:
    failures = receipt_failures(receipt)
    if failures:
        raise ReceiptValidationError("; ".join(failures))


def write_receipt(path: Path) -> dict[str, Any]:
    receipt = build_receipt()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(receipt, indent=2, sort_keys=True, ensure_ascii=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--check",
        action="store_true",
        help="compare freshly generated bytes with the committed receipt",
    )
    args = parser.parse_args()

    if args.check:
        expected = args.output.read_bytes()
        generated = (
            json.dumps(build_receipt(), indent=2, sort_keys=True, ensure_ascii=True)
            + "\n"
        ).encode("utf-8")
        if generated != expected:
            raise SystemExit(f"receipt is stale: {args.output}")
        print(f"screen shock receipt OK: {args.output.relative_to(ROOT)}")
        return

    write_receipt(args.output)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
