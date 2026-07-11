#!/usr/bin/env python3
"""Acceptance test for GitHub issue #539: MaxEnt multiplier counting and closure defect.

Issue #539 requires, on two successive regulator lattices:

1. The number of independent constraints and of Lagrange multipliers must match the
   parameter-space dimension displayed by the papers (N_con homogeneous global sums plus
   N_glob optional global charges), independent of the number of regulator cells. The
   rejected per-cell reading is counted alongside to show its dimension grows with cells.
2. The coarse-grained MaxEnt state must lie in the claimed family with a proved residual
   bound. This is the I-projection residual bound (compact paper, Lemma
   ``lem:closure-residual``; synthesis paper, Lemma 2.6b): the moment-matching
   I-projection R(lambda) onto the coarse homogeneous family exists and is unique, the
   closure defect is eps = D(sigma || omega_L(R(lambda))), and
   || sigma - omega_L(R(lambda)) ||_1 <= sqrt(2 * eps).

The model is the transverse-field Ising family on a spin-1/2 ring: density labels
O_1(x) = Z_x Z_{x+1} and O_2(x) = X_x, constrained through their homogeneous global sums,
so N_con = 2 and N_glob = 0. The refinement channel is decimation (partial trace over odd
sites), a completely positive trace-preserving coarse-graining.

Both regimes stated by the refinement-closure clause are exhibited:

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
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

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


def independent_operator_count(operators: list[np.ndarray]) -> int:
    """Rank of the span of {operators, identity} minus one, via the HS Gram matrix."""
    dim = operators[0].shape[0]
    basis = [np.eye(dim, dtype=complex)] + list(operators)
    gram = np.array([[np.trace(a.conj().T @ b) for b in basis] for a in basis])
    return int(np.linalg.matrix_rank(gram, tol=1e-9)) - 1


def gibbs_state(constraints: list[np.ndarray], lam: np.ndarray) -> tuple[np.ndarray, float]:
    """omega(lambda) = exp(-sum_a lambda_a S_a)/Z and log Z, via eigendecomposition."""
    ham = sum(l * s for l, s in zip(lam, constraints))
    energies, vectors = np.linalg.eigh(ham)
    shifted = energies - energies.min()
    weights = np.exp(-shifted)
    log_z = math.log(weights.sum()) - energies.min()
    rho = (vectors * (weights / weights.sum())) @ vectors.conj().T
    return rho, log_z


def decimate(rho: np.ndarray, n_sites: int) -> np.ndarray:
    """Partial trace over the odd sites of a ring, keeping sites 0, 2, 4, ..."""
    if n_sites % 2:
        raise ValueError("decimation channel needs an even number of fine sites")
    tensor = rho.reshape([2] * (2 * n_sites))
    for site in reversed(range(1, n_sites, 2)):
        tensor = np.trace(tensor, axis1=site, axis2=site + tensor.ndim // 2)
    kept = n_sites // 2
    return tensor.reshape(2**kept, 2**kept)


def relative_entropy(sigma: np.ndarray, rho: np.ndarray) -> float:
    """D(sigma || rho) in nats; both states are full rank in this model."""
    sig_e, sig_v = np.linalg.eigh(sigma)
    rho_e, rho_v = np.linalg.eigh(rho)
    sig_e = np.clip(sig_e.real, 1e-300, None)
    rho_e = np.clip(rho_e.real, 1e-300, None)
    log_rho = (rho_v * np.log(rho_e)) @ rho_v.conj().T
    entropy_term = float(np.sum(sig_e * np.log(sig_e)))
    cross_term = float(np.real(np.trace(sigma @ log_rho)))
    return entropy_term - cross_term


def trace_norm(delta: np.ndarray) -> float:
    return float(np.sum(np.abs(np.linalg.eigvalsh((delta + delta.conj().T) / 2))))


def duhamel_covariance(constraints: list[np.ndarray], lam: np.ndarray) -> np.ndarray:
    """Kubo-Mori covariance matrix of the constraints: the Hessian of log Z(lambda).

    K_ab = d^2 log Z / (d lambda_a d lambda_b); positive definite iff the constrained
    operators together with the identity are linearly independent, which is the strict
    convexity input of the I-projection lemma.
    """
    ham = sum(l * s for l, s in zip(lam, constraints))
    energies, vectors = np.linalg.eigh(ham)
    shifted = energies - energies.min()
    probs = np.exp(-shifted)
    probs /= probs.sum()
    rho = (vectors * probs) @ vectors.conj().T
    centered = [
        vectors.conj().T @ (s - np.real(np.trace(rho @ s)) * np.eye(s.shape[0])) @ vectors
        for s in constraints
    ]
    logp = np.log(np.clip(probs, 1e-300, None))
    pi, pj = np.meshgrid(probs, probs, indexing="ij")
    li, lj = np.meshgrid(logp, logp, indexing="ij")
    with np.errstate(divide="ignore", invalid="ignore"):
        kernel = np.where(np.abs(li - lj) > 1e-12, (pi - pj) / (li - lj), pi)
    n_con = len(constraints)
    cov = np.empty((n_con, n_con))
    for a in range(n_con):
        for b in range(n_con):
            cov[a, b] = float(np.real(np.sum(kernel * centered[a] * centered[b].T)))
    return (cov + cov.T) / 2


def i_projection(
    sigma: np.ndarray, constraints: list[np.ndarray], tol: float = 1e-11, max_iter: int = 200
) -> tuple[np.ndarray, float]:
    """Unique minimizer lambda* of D(sigma || omega(lambda')) by damped Newton descent.

    The objective is log Z(lambda') + sum_a lambda'_a <S_a>_sigma (strictly convex); its
    gradient is the moment mismatch <S_a>_sigma - <S_a>_{omega(lambda')} and its Hessian
    is the Duhamel covariance, so the minimizer is the moment-matching multiplier vector.
    Returns (lambda*, final gradient norm).
    """
    targets = np.array([float(np.real(np.trace(sigma @ s))) for s in constraints])

    def objective(lam: np.ndarray) -> float:
        _, log_z = gibbs_state(constraints, lam)
        return log_z + float(lam @ targets)

    lam = np.zeros(len(constraints))
    for _ in range(max_iter):
        rho, _ = gibbs_state(constraints, lam)
        moments = np.array([float(np.real(np.trace(rho @ s))) for s in constraints])
        grad = targets - moments
        if np.linalg.norm(grad) < tol:
            break
        hess = duhamel_covariance(constraints, lam)
        step = np.linalg.solve(hess + 1e-14 * np.eye(len(lam)), -grad)
        scale, base = 1.0, objective(lam)
        while scale > 1e-8 and objective(lam + scale * step) > base + 1e-15:
            scale /= 2
        lam = lam + scale * step
    rho, _ = gibbs_state(constraints, lam)
    moments = np.array([float(np.real(np.trace(rho @ s))) for s in constraints])
    return lam, float(np.linalg.norm(targets - moments))


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
    lam_star, moment_residual = i_projection(sigma, coarse)
    omega_coarse, _ = gibbs_state(coarse, lam_star)

    defect = max(relative_entropy(sigma, omega_coarse), 0.0)
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
        "moment_matching_residual": moment_residual,
        "duhamel_hessian_min_eigenvalue": hess_floor,
        "closure_defect_nats": defect,
        "trace_norm_residual": residual,
        "pinsker_residual_bound": pinsker_bound,
        "counts_match_displayed_dimension": fine_count == coarse_count == n_con,
        "projection_unique": hess_floor > 1e-9,
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
