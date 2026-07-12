#!/usr/bin/env python3
"""Verify the exact heat-spectrum identity for the transposition Cayley graph of S3.

This verifier establishes only finite-group linear algebra.  In particular, it
does not identify a physical OPH heat time.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import sympy as sp


CODE_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    CODE_ROOT
    / "particles"
    / "runs"
    / "flavor"
    / "s3_transposition_cayley_heat_shape_theorem.json"
)


Permutation = tuple[int, int, int]


def _compose(left: Permutation, right: Permutation) -> Permutation:
    """Return left after right, using image tuples on {0,1,2}."""

    return tuple(left[right[index]] for index in range(3))  # type: ignore[return-value]


def _permutations() -> list[Permutation]:
    from itertools import permutations

    return list(permutations(range(3)))


def _multiplicities(values: Iterable[sp.Expr]) -> dict[str, int]:
    result: dict[str, int] = {}
    for value in values:
        key = str(sp.simplify(value))
        result[key] = result.get(key, 0) + 1
    return dict(sorted(result.items(), key=lambda item: sp.sympify(item[0])))


def build_artifact() -> dict[str, object]:
    group = _permutations()
    generators: tuple[Permutation, ...] = (
        (1, 0, 2),
        (2, 1, 0),
        (0, 2, 1),
    )
    index = {element: position for position, element in enumerate(group)}

    adjacency = sp.zeros(len(group))
    for column, element in enumerate(group):
        for generator in generators:
            row = index[_compose(generator, element)]
            adjacency[row, column] = 1

    laplacian = 3 * sp.eye(len(group)) - adjacency
    adjacency_eigenvalues = adjacency.eigenvals()
    laplacian_eigenvalues = laplacian.eigenvals()

    tau = sp.symbols("tau", positive=True)
    gap_ratio = sp.simplify(
        (sp.exp(-3 * tau) - sp.exp(-6 * tau))
        / (1 - sp.exp(-3 * tau))
    )

    expected_adjacency = {-3: 1, 0: 4, 3: 1}
    expected_laplacian = {0: 1, 3: 4, 6: 1}
    checks = {
        "adjacency_spectrum": adjacency_eigenvalues == expected_adjacency,
        "laplacian_spectrum": laplacian_eigenvalues == expected_laplacian,
        "gap_ratio_identity": sp.simplify(gap_ratio - sp.exp(-3 * tau)) == 0,
        "three_regular": all(sum(adjacency[row, column] for row in range(6)) == 3 for column in range(6)),
    }

    return {
        "artifact": "oph_s3_transposition_cayley_heat_shape_theorem_v2",
        "claim_class": "exact_finite_group_linear_algebra_theorem",
        "proof_status": "exact_symbolic_verification_passed" if all(checks.values()) else "failed",
        "group_order": len(group),
        "transposition_count": len(generators),
        "adjacency_matrix": [[int(adjacency[i, j]) for j in range(6)] for i in range(6)],
        "adjacency_eigenvalues_with_multiplicity": _multiplicities(
            [value for value, count in adjacency_eigenvalues.items() for _ in range(count)]
        ),
        "laplacian_eigenvalues_with_multiplicity": _multiplicities(
            [value for value, count in laplacian_eigenvalues.items() for _ in range(count)]
        ),
        "distinct_heat_eigenvalues_sorted": ["exp(-6*tau)", "exp(-3*tau)", "1"],
        "adjacent_gap_ratio": "exp(-3*tau)",
        "checks": checks,
        "scope_boundary": (
            "The graph spectrum proves r=exp(-3*tau) after tau is supplied. "
            "It does not derive the proposed OPH identification "
            "tau_f=P/4-pi*alpha_U/5."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--print-json", action="store_true")
    args = parser.parse_args()

    artifact = build_artifact()
    text = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    if args.print_json:
        print(text, end="")
    else:
        print(f"saved: {args.output}")
    return 0 if artifact["proof_status"] == "exact_symbolic_verification_passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
