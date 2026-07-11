#!/usr/bin/env python3
"""Regression tests for GitHub issue #539's MaxEnt multiplier-family acceptance test."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from maxent_closure_acceptance import (  # noqa: E402
    cell_densities,
    decimate,
    duhamel_covariance,
    gibbs_state,
    global_sum_constraints,
    i_projection,
    independent_operator_count,
    relative_entropy,
    run_acceptance,
    run_lattice_pair,
    trace_norm,
)


def test_multiplier_count_is_cutoff_independent_for_global_sums() -> None:
    for n_sites in (3, 4, 6, 8):
        constraints = global_sum_constraints(n_sites)
        assert independent_operator_count(constraints) == 2
        assert sum(len(ops) for ops in cell_densities(n_sites)) == 2 * n_sites


def test_i_projection_moment_matches_and_is_unique() -> None:
    fine = global_sum_constraints(6)
    coarse = global_sum_constraints(3)
    omega, _ = gibbs_state(fine, np.array([0.7, 0.4]))
    sigma = decimate(omega, 6)
    lam_star, residual = i_projection(sigma, coarse)
    assert residual < 1e-9
    hessian_floor = np.linalg.eigvalsh(duhamel_covariance(coarse, lam_star)).min()
    assert hessian_floor > 1e-9


def test_pinsker_residual_bound() -> None:
    run = run_lattice_pair(6, np.array([0.7, 0.4]))
    assert run["trace_norm_residual"] <= run["pinsker_residual_bound"] + 1e-9
    assert run["closure_defect_nats"] > 1e-6


def test_transverse_field_subfamily_is_exactly_closed() -> None:
    run = run_lattice_pair(6, np.array([0.0, 0.4]))
    assert run["closure_defect_nats"] < 1e-10
    assert abs(run["induced_map_R_multipliers"][0]) < 1e-8
    assert abs(run["induced_map_R_multipliers"][1] - 0.4) < 1e-8


def test_relative_entropy_is_zero_only_on_the_diagonal() -> None:
    coarse = global_sum_constraints(3)
    rho, _ = gibbs_state(coarse, np.array([0.3, 0.5]))
    tau, _ = gibbs_state(coarse, np.array([0.4, 0.5]))
    assert relative_entropy(rho, rho) < 1e-12
    gap = relative_entropy(rho, tau)
    assert gap > 0
    assert trace_norm(rho - tau) <= (2 * gap) ** 0.5 + 1e-9


def test_full_acceptance_run_passes() -> None:
    results = run_acceptance()
    assert results["all_checks_pass"], results["checks"]


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"PASS  {name}")
    print("all tests passed")
