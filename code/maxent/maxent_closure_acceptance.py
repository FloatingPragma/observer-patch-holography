#!/usr/bin/env python3
"""Acceptance test for GitHub issue #539: MaxEnt multiplier counting and closure defect.

Issue #539 requires, on two successive regulator lattices:

1. The number of independent constraints and of Lagrange multipliers must match the
   parameter-space dimension displayed by the papers (N_con homogeneous global sums plus
   N_glob optional global charges), independent of the number of regulator cells. The
   rejected per-cell reading is counted alongside to show its dimension grows with cells.
2. The coarse-grained MaxEnt state must lie in the claimed family with a proved residual
   bound. This is the I-projection residual bound (spacetime and Einstein paper, Lemma
   ``lem:closure-residual``; synthesis paper, Lemma 2.6b): the moment-matching
   I-projection R(lambda) onto the coarse homogeneous family exists and is unique, the
   closure defect is eps = D(sigma || omega_L(R(lambda))), and
   || sigma - omega_L(R(lambda)) ||_1 <= sqrt(2 * eps).

The model is the transverse-field Ising family on a spin-1/2 ring: density labels
O_1(x) = Z_x Z_{x+1} and O_2(x) = X_x, constrained through their homogeneous global sums,
so N_con = 2 and N_glob = 0. The refinement channel is decimation (partial trace over odd
sites), a completely positive trace-preserving coarse-graining.

Both regimes relevant to the separate optimizer-pushforward interface are exhibited:

- generic multipliers: the closure defect is strictly positive, so closure under one fixed
  finite exponential family is a substantive renormalization condition, and the residual
  bound quantifies exactly how far the coarse-grained state sits from the family;
- the transverse-field product subfamily (lambda_ZZ = 0): decimation lands exactly in the
  coarse family, the defect vanishes, and R acts as the identity on the multipliers.

Run directly for a full receipt (written to runs/maxent_closure_acceptance_receipt.json):

    python3 maxent_closure_acceptance.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from quantum_information import (
    partial_trace,
    relative_entropy as _relative_entropy,
)

from maxent.information_projection import (
    constrained_hamiltonian, duhamel_covariance, gibbs_state,
    independent_operator_count, i_projection, project_information,
    projection_diagnostics,
)

HERE = Path(__file__).resolve().parent

# Roundoff allowance for a negative relative-entropy estimate; the state checks
# in quantum_information admit trace and positivity deviations up to 1e-12.
NEGATIVE_DEFECT_TOL = 1e-12

PAULI_I = np.eye(2, dtype=complex)
PAULI_X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
PAULI_Z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)


def site_operator(op: np.ndarray, site: int, n_sites: int) -> np.ndarray:
    """op acting on `site` of an n_sites spin-1/2 chain, identity elsewhere."""
    out = np.array([[1.0 + 0.0j]])
    for k in range(n_sites):
        out = np.kron(out, op if k == site else PAULI_I)
    return out


def cell_densities(n_sites: int) -> list[list[np.ndarray]]:
    """Per-cell density operators [O_1(x)]_x, [O_2(x)]_x on a ring of n_sites cells."""
    zz = [
        site_operator(PAULI_Z, x, n_sites) @ site_operator(PAULI_Z, (x + 1) % n_sites, n_sites)
        for x in range(n_sites)
    ]
    xs = [site_operator(PAULI_X, x, n_sites) for x in range(n_sites)]
    return [zz, xs]


def global_sum_constraints(n_sites: int) -> list[np.ndarray]:
    """The homogeneous global-sum constrained operators S_a = sum_x O_a(x)."""
    return [sum(ops) for ops in cell_densities(n_sites)]


def decimate(rho: np.ndarray, n_sites: int) -> np.ndarray:
    """Partial trace over the odd sites of a ring, keeping sites 0, 2, 4, ..."""
    if type(n_sites) is not int or n_sites <= 0 or n_sites % 2:
        raise ValueError("decimation channel needs an even number of fine sites")
    return partial_trace(rho, [2]*n_sites, list(range(0,n_sites,2)))


def relative_entropy(sigma: np.ndarray, rho: np.ndarray) -> float:
    """D(sigma || rho); preserve the historical keyword-argument convention."""
    return _relative_entropy(sigma,rho)


def trace_norm(delta: np.ndarray) -> float:
    return float(np.sum(np.abs(np.linalg.eigvalsh((delta + delta.conj().T) / 2))))


def closure_defect(sigma: np.ndarray, omega_coarse: np.ndarray) -> float:
    """Closure defect D(sigma || omega_coarse), fail-closed on negative drift.

    The Umegaki entropy is nonnegative for admissible states, so a value below
    roundoff means the pair handed in is not a valid comparison. Clamping it to
    zero would report exact closure for a failed projection, which is the one
    outcome this acceptance test must never fake.
    """
    raw = relative_entropy(sigma, omega_coarse)
    if not math.isfinite(raw):
        raise ValueError(f"non-finite closure defect nats: {raw}")
    if raw < -NEGATIVE_DEFECT_TOL:
        raise ValueError(
            "negative relative entropy in the closure defect: "
            f"{raw} nats is below -{NEGATIVE_DEFECT_TOL} nats")
    return abs(raw)


def run_lattice_pair(n_fine: int, lam_fine: np.ndarray) -> dict:
    """Run the full issue-#539 acceptance test on the lattice pair (n_fine, n_fine/2)."""
    n_coarse = n_fine // 2
    fine = global_sum_constraints(n_fine)
    coarse = global_sum_constraints(n_coarse)
    n_con = len(fine)

    fine_count = independent_operator_count(fine)
    coarse_count = independent_operator_count(coarse)
    per_cell_fine = sum(len(ops) for ops in cell_densities(n_fine))
    per_cell_coarse = sum(len(ops) for ops in cell_densities(n_coarse))

    omega_fine, _ = gibbs_state(fine, lam_fine)
    sigma = decimate(omega_fine, n_fine)
    projection = project_information(sigma, coarse)
    lam_star, omega_coarse = projection.multipliers, projection.state

    defect = closure_defect(sigma, omega_coarse)
    residual = trace_norm(sigma - omega_coarse)
    pinsker_bound = math.sqrt(2 * defect)
    hess_floor = float(np.linalg.eigvalsh(duhamel_covariance(coarse, lam_star)).min())

    return {
        "lattice_pair": [n_fine, n_coarse],
        "fine_multipliers": list(map(float, lam_fine)),
        "displayed_dimension_N_con_plus_N_glob": n_con,
        "independent_global_sum_constraints": {"fine": fine_count, "coarse": coarse_count},
        "rejected_per_cell_constraint_count": {"fine": per_cell_fine, "coarse": per_cell_coarse},
        "induced_map_R_multipliers": list(map(float, lam_star)),
        "moment_matching_residual": projection.raw_moment_residual,
        "normalized_moment_matching_residual": projection.normalized_residual,
        "projection_optimality_gap_bound_nats": projection.optimality_gap_bound,
        "projection_trace_distance_bound": projection.trace_distance_bound,
        "projection_iterations": projection.iterations,
        "duhamel_hessian_min_eigenvalue": hess_floor,
        "closure_defect_nats": defect,
        "trace_norm_residual": residual,
        "pinsker_residual_bound": pinsker_bound,
        "counts_match_displayed_dimension": fine_count == coarse_count == n_con,
        "projection_unique": coarse_count == n_con,
        "residual_bound_holds": residual <= pinsker_bound + 1e-9,
    }


def run_acceptance() -> dict:
    generic = np.array([0.7, 0.4])
    closed = np.array([0.0, 0.4])
    results = {
        "issue": 539,
        "model": "transverse-field Ising densities O_1(x)=Z_x Z_{x+1}, O_2(x)=X_x on a ring;"
        " homogeneous global-sum constraints; decimation refinement channel",
        "generic_branch": [run_lattice_pair(6, generic), run_lattice_pair(8, generic)],
        "closed_subfamily": run_lattice_pair(6, closed),
    }
    closed_run = results["closed_subfamily"]
    checks = {
        "multiplier_and_constraint_counts_cutoff_independent": all(
            r["counts_match_displayed_dimension"]
            for r in results["generic_branch"] + [closed_run]
        ),
        "i_projection_unique_everywhere": all(
            r["projection_unique"] for r in results["generic_branch"] + [closed_run]
        ),
        "pinsker_residual_bound_holds_everywhere": all(
            r["residual_bound_holds"] for r in results["generic_branch"] + [closed_run]
        ),
        "generic_closure_defect_strictly_positive": all(
            r["closure_defect_nats"] > 1e-6 for r in results["generic_branch"]
        ),
        "transverse_field_subfamily_exactly_closed": closed_run["closure_defect_nats"] < 1e-10
        and abs(closed_run["induced_map_R_multipliers"][0]) < 1e-8
        and abs(closed_run["induced_map_R_multipliers"][1] - 0.4) < 1e-8,
    }
    results["checks"] = checks
    results["all_checks_pass"] = all(checks.values())
    return results


def main() -> int:
    results = run_acceptance()
    receipt = HERE / "runs" / "maxent_closure_acceptance_receipt.json"
    receipt.parent.mkdir(exist_ok=True)
    receipt.write_text(json.dumps(results, indent=2) + "\n")
    for name, passed in results["checks"].items():
        print(f"{'PASS' if passed else 'FAIL'}  {name}")
    for run in results["generic_branch"]:
        print(
            f"lattices {run['lattice_pair']}: constraints/multipliers = "
            f"{run['independent_global_sum_constraints']} (displayed "
            f"{run['displayed_dimension_N_con_plus_N_glob']}; per-cell reading would be "
            f"{run['rejected_per_cell_constraint_count']}), closure defect = "
            f"{run['closure_defect_nats']:.6f} nats, trace-norm residual "
            f"{run['trace_norm_residual']:.6f} <= bound {run['pinsker_residual_bound']:.6f}"
        )
    print(f"receipt: {receipt}")
    return 0 if results["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
