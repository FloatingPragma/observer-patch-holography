#!/usr/bin/env python3
"""Build the W5 residual-rotation stabiliser-spectrum receipt.

Here W5 is realised as the five-dimensional space of real symmetric traceless
3 by 3 matrices.  For a canonical rotation around the z axis the matrix splits
into a scalar planar block, a planar spin-2 anisotropy, and a transverse
spin-1 vector.  The fixed-vector determinants are

    det(R(theta) - I)  = 2(1 - cos(theta)),
    det(R(2 theta) - I) = 2(1 - cos(2 theta)).

Both determinants are nonzero for canonical rotations of order 3 and 5, so
only the scalar planar block remains.  Tracelessness then forces axial
symmetry and a double eigenvalue.  For order 2 the spin-2 determinant vanishes,
leaving a three-dimensional fixed locus and two projective parameters.

The receipt reclassifies W5_ORB as a potential-selection problem.  It does not
claim that simple-spectrum W5 points or predictive A5 constructions are
universally impossible.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "w5_stabiliser_spectrum_bound.json"

Matrix = list[list[float]]


def _matmul(left: Matrix, right: Matrix) -> Matrix:
    return [
        [
            sum(left[row][inner] * right[inner][column] for inner in range(3))
            for column in range(3)
        ]
        for row in range(3)
    ]


def _transpose(matrix: Matrix) -> Matrix:
    return [[matrix[column][row] for column in range(3)] for row in range(3)]


def _rotation_z(theta: float) -> Matrix:
    cosine = math.cos(theta)
    sine = math.sin(theta)
    return [
        [cosine, -sine, 0.0],
        [sine, cosine, 0.0],
        [0.0, 0.0, 1.0],
    ]


def _invariance_residual(rotation: Matrix, matrix: Matrix) -> Matrix:
    transformed = _matmul(_matmul(rotation, matrix), _transpose(rotation))
    return [
        [transformed[row][column] - matrix[row][column] for column in range(3)]
        for row in range(3)
    ]


def _max_abs(matrix: Matrix) -> float:
    return max(abs(entry) for row in matrix for entry in row)


def _w5_basis() -> list[Matrix]:
    """Basis for A=[[a,b,c],[b,d,e],[c,e,-a-d]]."""

    return [
        [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
        [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
        [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        [[0.0, 0.0, 1.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]],
    ]


def _constraint_matrix(rotation: Matrix) -> list[list[float]]:
    """Return six symmetric-entry equations in the five W5 coordinates."""

    positions = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
    basis_residuals = [_invariance_residual(rotation, basis) for basis in _w5_basis()]
    return [
        [residual[row][column] for residual in basis_residuals]
        for row, column in positions
    ]


def _matrix_rank(rows: Sequence[Sequence[float]], tolerance: float = 1e-10) -> int:
    matrix = [list(row) for row in rows]
    if not matrix:
        return 0
    row_count = len(matrix)
    column_count = len(matrix[0])
    rank = 0
    for column in range(column_count):
        pivot = max(range(rank, row_count), key=lambda row: abs(matrix[row][column]))
        if abs(matrix[pivot][column]) <= tolerance:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        pivot_value = matrix[rank][column]
        matrix[rank] = [entry / pivot_value for entry in matrix[rank]]
        for row in range(row_count):
            if row == rank:
                continue
            factor = matrix[row][column]
            if abs(factor) <= tolerance:
                continue
            matrix[row] = [
                entry - factor * pivot_entry
                for entry, pivot_entry in zip(matrix[row], matrix[rank], strict=True)
            ]
        rank += 1
        if rank == row_count:
            break
    return rank


def _rotation_record(order: int) -> dict[str, object]:
    theta = 2.0 * math.pi / order
    rotation = _rotation_z(theta)
    rank = _matrix_rank(_constraint_matrix(rotation))
    fixed_dimension = 5 - rank
    vector_determinant = 2.0 * (1.0 - math.cos(theta))
    spin2_determinant = 2.0 * (1.0 - math.cos(2.0 * theta))
    expected_dimension = {2: 3, 3: 1, 5: 1}[order]
    if fixed_dimension != expected_dimension:
        raise AssertionError(
            f"order-{order} fixed-locus dimension {fixed_dimension}, expected {expected_dimension}"
        )

    exact_determinants = {
        2: {
            "transverse_vector": "4",
            "planar_spin2_anisotropy": "0",
        },
        3: {
            "transverse_vector": "3",
            "planar_spin2_anisotropy": "3",
        },
        5: {
            "transverse_vector": "(5 - sqrt(5))/2",
            "planar_spin2_anisotropy": "(5 + sqrt(5))/2",
        },
    }
    if order in (3, 5) and (
        vector_determinant <= 1e-12 or spin2_determinant <= 1e-12
    ):
        raise AssertionError(f"order-{order} nondegeneracy determinant vanished")
    if order == 2 and abs(spin2_determinant) > 1e-12:
        raise AssertionError("order-2 planar spin-2 determinant must vanish")

    record: dict[str, object] = {
        "rotation_order": order,
        "canonical_axis": "(0,0,1)",
        "rotation_angle_degrees": round(math.degrees(theta), 12),
        "linear_constraint_rank_in_W5": rank,
        "fixed_locus_dimension": fixed_dimension,
        "nonzero_projective_dimension": fixed_dimension - 1,
        "fixed_block_determinants": {
            "transverse_vector_det_R_minus_I": {
                "formula": "2*(1-cos(theta))",
                "exact": exact_determinants[order]["transverse_vector"],
                "numeric": round(vector_determinant, 12),
            },
            "planar_spin2_det_R2_minus_I": {
                "formula": "2*(1-cos(2*theta))",
                "exact": exact_determinants[order]["planar_spin2_anisotropy"],
                "numeric": round(spin2_determinant, 12),
            },
        },
    }
    if order in (3, 5):
        record.update(
            {
                "general_invariant_form": "alpha*(n*n^T-I/3)",
                "canonical_matrix": "diag(-alpha/3,-alpha/3,2*alpha/3)",
                "spectrum": ["2*alpha/3", "-alpha/3", "-alpha/3"],
                "eigenvalue_multiplicities": [1, 2],
                "simple_spectrum_possible": False,
            }
        )
    else:
        record.update(
            {
                "general_invariant_form": "[[a,b,0],[b,d,0],[0,0,-a-d]]",
                "free_coordinates": ["a", "b", "d"],
                "simple_spectrum_possible": True,
                "simple_spectrum_witness": {
                    "matrix": [[1, 0, 0], [0, 2, 0], [0, 0, -3]],
                    "spectrum": [1, 2, -3],
                },
            }
        )
    return record


def build_artifact() -> dict[str, object]:
    records = {
        f"C{order}": _rotation_record(order)
        for order in (2, 3, 5)
    }

    alpha = 3.0
    axial_matrix = [
        [-alpha / 3.0, 0.0, 0.0],
        [0.0, -alpha / 3.0, 0.0],
        [0.0, 0.0, 2.0 * alpha / 3.0],
    ]
    for order in (3, 5):
        residual = _max_abs(
            _invariance_residual(_rotation_z(2.0 * math.pi / order), axial_matrix)
        )
        if residual > 1e-12:
            raise AssertionError(f"order-{order} axial witness is not invariant")

    c2_witness = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, -3.0]]
    c2_residual = _max_abs(_invariance_residual(_rotation_z(math.pi), c2_witness))
    if c2_residual > 1e-12:
        raise AssertionError("C2 simple-spectrum witness is not invariant")

    return {
        "artifact": "oph_w5_stabiliser_spectrum_bound",
        "schema_version": 1,
        "determinism": {
            "timestamp_omitted": True,
            "external_measured_inputs": [],
        },
        "representation": {
            "name": "W5",
            "realisation": "Symmetric traceless real 3x3 matrices",
            "ambient_dimension": 5,
            "group_action": "A -> R*A*R^T",
        },
        "symbolic_fixed_space_argument": {
            "canonical_axis_reduction": (
                "Every real three-dimensional rotation axis is conjugate to the z axis."
            ),
            "block_coordinates": {
                "planar_scalar": "t",
                "planar_spin2_anisotropy": ["p", "q"],
                "transverse_vector": ["u", "v"],
                "axial_entry": "-2*t by tracelessness",
            },
            "rotation_actions": {
                "planar_scalar": "fixed",
                "planar_spin2_anisotropy": "rotation by 2*theta",
                "transverse_vector": "rotation by theta",
            },
            "consequence_for_order_at_least_3": (
                "For canonical orders 3 and 5 both non-scalar fixed blocks vanish, "
                "leaving alpha*(n*n^T-I/3) with a double eigenvalue."
            ),
        },
        "canonical_rotation_checks": records,
        "hard_self_tests": {
            "orders_3_and_5_constraint_rank": 4,
            "orders_3_and_5_fixed_dimension": 1,
            "orders_3_and_5_axial_witness_max_residual": 0.0,
            "order_2_constraint_rank": 2,
            "order_2_fixed_dimension": 3,
            "order_2_projective_dimension": 2,
            "order_2_simple_spectrum_witness_max_residual": round(c2_residual, 12),
            "passed": True,
        },
        "w5_orb_reclassification": {
            "previous_gap_type_rejected": "symmetry_geometry_alone_selects_unique_simple_spectrum_orbit",
            "current_gap_type": "screen_derived_potential_selection_required",
            "reason": (
                "C3 and C5 invariance force a double eigenvalue. C2 invariance leaves "
                "a three-dimensional linear locus, hence two parameters after scaling, "
                "which is exactly enough freedom to carry two mass ratios without "
                "predicting them from symmetry alone."
            ),
            "symmetry_alone_sufficient": False,
            "screen_derived_potential_required": True,
        },
        "scope": {
            "universal_impossibility_claimed": False,
            "not_excluded": [
                "simple_spectrum_points_with_C2_or_trivial_stabiliser",
                "selection_by_a_specific_screen_derived_A5_invariant_potential",
                "predictive_models_with_additional_dynamical_input",
            ],
            "literature_claims_used": False,
        },
    }


def _write_json_lf(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build the W5 residual-rotation stabiliser-spectrum receipt."
    )
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()
    output = Path(args.output)
    _write_json_lf(output, build_artifact())
    print(f"saved: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
