#!/usr/bin/env python3
"""Finite Z2 lattice-gauge diagnostic for the ground-state-transform receipt.

The Yang--Mills gap paper assumes a finite ground-state-transform and
cross-fiber receipt: a unitary ``U_r`` with ``U_r Omega_r = 1`` such that
``U_r H_r U_r^{-1} = sum_C D_C`` with each ``D_C = c_C (I - E_C)`` a heat-bath
collar generator whose rate ``c_C`` is independent of the repaired value.
This module evaluates that receipt numerically on two small gauge systems:
Z2 lattice gauge theory on an ``L x L`` periodic spatial torus in the
gauge-invariant (Gauss-law) sector.  The free control has an analytic exact
identity; interacting eigendecompositions and matrix logarithms use float64
arithmetic and the serialized tolerances.

Two transfer objects are tested.

* ``wilson``: the reflection-positive Wilson transfer matrix
  ``T = exp(beta_s P / 2) K exp(beta_s P / 2)`` with ``K(s, s') =
  prod_l exp(beta_t s_l s'_l)``, and ``H = -log(T / lambda_max)``.
* ``kogut_susskind``: the Hamiltonian ``H = -lam sum_l X_l - sum_p U_p``,
  shifted so that its ground energy is zero.

For each, ``Omega`` is the Perron ground state, ``pi = Omega^2`` the
stationary time-zero law, ``L = D_Omega^{-1} H D_Omega`` the ground-state
(Doob) transform, and ``E_l`` the ``pi``-preserving conditional expectation
on the single-link collar fiber ``{o, X_l o}``.

The script reports

* the relative Frobenius residual of the best constant-rate fit
  ``L ~ sum_l c_l (I - E_l)``;
* the exact fiber-dependent rates ``c_l(o)`` that reproduce ``L`` when ``L``
  is single-flip, and their spread ``max/min`` (the cross-fiber receipt
  requires spread 1);
* the Dobrushin influence ``eta_*`` of ``pi`` under single-link heat bath,
  the floor ``delta_* = c_*(1 - eta_*)``, and the exact spectral gaps of
  ``H`` and of the unit-rate heat-bath generator.

Everything here is a finite diagnostic on a toy gauge system.  It is not a
physical compact-simple-gauge receipt, and the output JSON says so.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

SCHEMA = "oph.yang_mills.z2_finite_transfer_receipt.v2"


# ----------------------------------------------------------------------------
# Lattice geometry
# ----------------------------------------------------------------------------


def link_index(L: int, x: int, y: int, direction: int) -> int:
    """Index of the link leaving site ``(x, y)`` in ``direction`` (0 = x, 1 = y)."""
    return 2 * ((x % L) * L + (y % L)) + direction


def plaquette_masks(L: int) -> list[int]:
    masks = []
    for x in range(L):
        for y in range(L):
            m = 0
            m |= 1 << link_index(L, x, y, 0)
            m |= 1 << link_index(L, x + 1, y, 1)
            m |= 1 << link_index(L, x, y + 1, 0)
            m |= 1 << link_index(L, x, y, 1)
            masks.append(m)
    return masks


def star_masks(L: int) -> list[int]:
    """Coboundary of each site: links with exactly that site as an endpoint."""
    masks = []
    for x in range(L):
        for y in range(L):
            m = 0
            m |= 1 << link_index(L, x, y, 0)
            m |= 1 << link_index(L, x, y, 1)
            m |= 1 << link_index(L, x - 1, y, 0)
            m |= 1 << link_index(L, x, y - 1, 1)
            masks.append(m)
    return masks


def gauge_group_masks(L: int) -> np.ndarray:
    """All gauge transformations as XOR masks (one per subset of sites)."""
    stars = star_masks(L)
    n_sites = len(stars)
    masks = np.zeros(1 << n_sites, dtype=np.int64)
    for subset in range(1 << n_sites):
        m = 0
        for site in range(n_sites):
            if subset >> site & 1:
                m ^= stars[site]
        masks[subset] = m
    return np.unique(masks)


def popcount(values: np.ndarray) -> np.ndarray:
    return np.bitwise_count(values.astype(np.uint64)).astype(np.int64)


# ----------------------------------------------------------------------------
# Orbit space
# ----------------------------------------------------------------------------


class Z2GaugeOrbits:
    """Gauge orbits of Z2 link configurations on the L x L periodic torus."""

    def __init__(self, L: int) -> None:
        self.L = L
        self.n_links = 2 * L * L
        self.n_configs = 1 << self.n_links
        self.plaquettes = plaquette_masks(L)
        self.gauge = gauge_group_masks(L)
        configs = np.arange(self.n_configs, dtype=np.int64)
        # canonical representative: minimum over the gauge orbit
        reps = configs.copy()
        for g in self.gauge:
            reps = np.minimum(reps, configs ^ int(g))
        self.rep_of_config = reps
        self.reps, inverse = np.unique(reps, return_inverse=True)
        self.orbit_of_config = inverse
        self.n_orbits = len(self.reps)
        self.orbit_size = np.bincount(inverse, minlength=self.n_orbits)
        if not np.all(self.orbit_size == self.orbit_size[0]):
            raise RuntimeError("gauge action is not free; weighted orbit basis needed")
        # members[o] = all configurations in orbit o
        self.members = (self.reps[:, None] ^ self.gauge[None, :]).astype(np.int64)
        # plaquette sum P(o) = sum_p U_p for the representative
        P = np.zeros(self.n_orbits, dtype=np.int64)
        for m in self.plaquettes:
            parity = popcount(self.reps & m) & 1
            P += 1 - 2 * parity
        self.plaquette_sum = P
        # single-link flip as an involution on orbits
        self.flip = np.zeros((self.n_links, self.n_orbits), dtype=np.int64)
        for l in range(self.n_links):
            self.flip[l] = self.orbit_of_config[self.reps ^ (1 << l)]
            if np.any(self.flip[l] == np.arange(self.n_orbits)):
                raise RuntimeError("single-link flip has a fixed orbit")

    # -- transfer objects ---------------------------------------------------

    def wilson_transfer(self, beta_s: float, beta_t: float) -> np.ndarray:
        """Symmetric gauge-invariant Wilson transfer matrix on orbit space."""
        n = self.n_orbits
        T = np.zeros((n, n))
        weight = np.exp(0.5 * beta_s * self.plaquette_sum)
        for o in range(n):
            d = popcount(self.members ^ self.reps[o])  # (n_orbits, |G|)
            kin = np.exp(beta_t * (self.n_links - 2 * d)).sum(axis=1)
            T[o] = weight[o] * kin * weight
        return T

    def kogut_susskind(self, lam: float) -> np.ndarray:
        n = self.n_orbits
        H = np.zeros((n, n))
        H[np.arange(n), np.arange(n)] = -self.plaquette_sum.astype(float)
        for l in range(self.n_links):
            H[np.arange(n), self.flip[l]] -= lam
        return H


# ----------------------------------------------------------------------------
# Receipt evaluation
# ----------------------------------------------------------------------------


def symmetric_log_hamiltonian(T: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    w, V = np.linalg.eigh(T)
    if w[0] <= 0:
        raise RuntimeError("transfer matrix is not positive definite")
    lam_max = w[-1]
    H = (V * (-np.log(w / lam_max))) @ V.T
    omega = V[:, -1]
    if omega.sum() < 0:
        omega = -omega
    if np.any(omega <= 0):
        raise RuntimeError("Perron vector is not strictly positive")
    return H, omega, float(lam_max)


def ground_state(H: np.ndarray) -> tuple[np.ndarray, float]:
    w, V = np.linalg.eigh(H)
    omega = V[:, 0]
    if omega.sum() < 0:
        omega = -omega
    if np.any(omega <= 0):
        raise RuntimeError("Perron vector is not strictly positive")
    return omega, float(w[0])


def doob_transform(H: np.ndarray, omega: np.ndarray, e0: float) -> np.ndarray:
    return (H - e0 * np.eye(len(omega))) * omega[None, :] / omega[:, None]


def heat_bath_projectors(orbits: Z2GaugeOrbits, pi: np.ndarray) -> list[np.ndarray]:
    n = orbits.n_orbits
    projectors = []
    for l in range(orbits.n_links):
        partner = orbits.flip[l]
        denom = pi + pi[partner]
        E = np.zeros((n, n))
        idx = np.arange(n)
        E[idx, idx] = pi / denom
        E[idx, partner] = pi[partner] / denom
        projectors.append(E)
    return projectors


def constant_rate_fit(Lgen: np.ndarray, projectors: list[np.ndarray]) -> dict[str, Any]:
    n = Lgen.shape[0]
    basis = np.stack([(np.eye(n) - E).ravel() for E in projectors], axis=1)
    target = Lgen.ravel()
    coeff, *_ = np.linalg.lstsq(basis, target, rcond=None)
    residual = target - basis @ coeff
    rel = float(np.linalg.norm(residual) / np.linalg.norm(target))
    return {
        "rates": [float(c) for c in coeff],
        "rate_min": float(coeff.min()),
        "rate_max": float(coeff.max()),
        "relative_frobenius_residual": rel,
    }


def fiber_dependent_rates(
    Lgen: np.ndarray, orbits: Z2GaugeOrbits, pi: np.ndarray
) -> dict[str, Any]:
    """Exact collar rates ``c_l(o)`` from the off-diagonal of ``L``.

    For ``o' = X_l o`` the heat-bath form gives
    ``L(o, o') = -c_l(o) pi(o') / (pi(o) + pi(o'))``.  The returned spread is
    ``max c / min c`` over all links and orbits; the cross-fiber receipt needs
    spread exactly one.  Off-diagonal mass of ``L`` outside single-link flips is
    returned separately; it must vanish for the single-flip form to be exact.
    """
    n = orbits.n_orbits
    idx = np.arange(n)
    rates = np.zeros((orbits.n_links, n))
    single_flip_mask = np.zeros((n, n), dtype=bool)
    for l in range(orbits.n_links):
        partner = orbits.flip[l]
        off = Lgen[idx, partner]
        rates[l] = -off * (pi + pi[partner]) / pi[partner]
        single_flip_mask[idx, partner] = True
    off_diag = Lgen - np.diag(np.diag(Lgen))
    outside = off_diag[~single_flip_mask]
    return {
        "rate_min": float(rates.min()),
        "rate_max": float(rates.max()),
        "spread_max_over_min": float(rates.max() / rates.min()) if rates.min() > 0 else math.inf,
        "all_rates_positive": bool(np.all(rates > 0)),
        "offdiagonal_mass_outside_single_flip": float(np.abs(outside).sum()),
        "offdiagonal_mass_single_flip": float(np.abs(off_diag[single_flip_mask]).sum()),
    }


def dobrushin_influence(orbits: Z2GaugeOrbits, pi: np.ndarray) -> dict[str, Any]:
    """Total-variation influence matrix of the single-link conditional kernels."""
    n_links = orbits.n_links
    idx = np.arange(orbits.n_orbits)
    # conditional probability that link l is in its representative state at o:
    # kernel(o) = pi(o) / (pi(o) + pi(X_l o)), as a function of the other links.
    kernel = np.zeros((n_links, orbits.n_orbits))
    for l in range(n_links):
        kernel[l] = pi / (pi + pi[orbits.flip[l]])
    influence = np.zeros((n_links, n_links))
    for l in range(n_links):
        for u in range(n_links):
            if u == l:
                continue
            # changing link u: compare kernel at o and at X_u o.  The kernel is
            # stated for the representative state of link l, which X_u does not
            # change, so the TV distance is the absolute difference.
            diff = np.abs(kernel[l] - kernel[l][orbits.flip[u]])
            influence[l, u] = float(diff.max())
    row_sums = influence.sum(axis=1)
    return {
        "eta_star": float(row_sums.max()),
        "row_sums": [float(r) for r in row_sums],
    }


def spectral_gap(M: np.ndarray, pi: np.ndarray | None = None) -> float:
    """Smallest nonzero eigenvalue of a generator symmetrisable by ``pi``."""
    if pi is not None:
        s = np.sqrt(pi)
        M = M * s[:, None] / s[None, :]
    w = np.linalg.eigvalsh(0.5 * (M + M.T))
    w = np.sort(w)
    return float(w[1])


def evaluate(orbits: Z2GaugeOrbits, transfer: str, **params: float) -> dict[str, Any]:
    if transfer == "wilson":
        T = orbits.wilson_transfer(params["beta_s"], params["beta_t"])
        H, omega, lam_max = symmetric_log_hamiltonian(T)
        e0 = 0.0
        extra = {"lambda_max": lam_max}
    elif transfer == "kogut_susskind":
        H = orbits.kogut_susskind(params["lam"])
        omega, e0 = ground_state(H)
        extra = {"ground_energy": e0}
    else:
        raise ValueError(transfer)
    pi = omega**2
    pi = pi / pi.sum()
    Lgen = doob_transform(H, omega, e0)
    projectors = heat_bath_projectors(orbits, pi)
    unit_heat_bath = sum(np.eye(orbits.n_orbits) - E for E in projectors)
    fit = constant_rate_fit(Lgen, projectors)
    fibre = fiber_dependent_rates(Lgen, orbits, pi)
    dob = dobrushin_influence(orbits, pi)
    eta = dob["eta_star"]
    result: dict[str, Any] = {
        "transfer": transfer,
        "parameters": params,
        "n_orbits": int(orbits.n_orbits),
        "n_links": int(orbits.n_links),
        "doob_generator_rows_sum_zero": bool(
            np.allclose(Lgen.sum(axis=1), 0, atol=1e-9)
        ),
        "doob_generator_offdiagonal_nonpositive": bool(
            np.all(Lgen - np.diag(np.diag(Lgen)) <= 1e-12)
        ),
        "constant_rate_fit": fit,
        "fiber_dependent_rates": fibre,
        "dobrushin": {
            "eta_star": eta,
            "dobrushin_condition_holds": bool(eta < 1),
            "unit_rate_floor_c_star_times_1_minus_eta": float(max(0.0, 1 - eta)),
        },
        "spectral": {
            "gap_H": spectral_gap(Lgen, pi),
            "gap_unit_rate_heat_bath": spectral_gap(unit_heat_bath, pi),
            "pi_min": float(pi.min()),
            "pi_max": float(pi.max()),
        },
    }
    if transfer == "kogut_susskind":
        analytic_floor = 2.0 * params["lam"]
        result["variable_rate_floor"] = {
            "identity": "c_l(o) = lambda * (r_l(o) + 1/r_l(o))",
            "analytic_lower_bound": "c_l(o) >= 2 * lambda by AM-GM",
            "lower_bound_value": analytic_floor,
            "numerical_min_respects_bound": bool(
                fibre["rate_min"] >= analytic_floor - 1e-10
            ),
            "scope": (
                "finite Kogut-Susskind Doob transform; quotient-space "
                "approximate tensorization and continuum transfer not proved"
            ),
        }
    result.update(extra)
    return result


def run(L_values: list[int], betas: list[float], lams: list[float]) -> dict[str, Any]:
    runs = []
    for L in L_values:
        orbits = Z2GaugeOrbits(L)
        runs.append({"L": L, "transfer": "wilson", "control": "beta_s_zero",
                     **evaluate(orbits, "wilson", beta_s=0.0, beta_t=0.5)})
        for beta in betas:
            runs.append({"L": L, **evaluate(orbits, "wilson", beta_s=beta, beta_t=beta)})
        for lam in lams:
            runs.append({"L": L, **evaluate(orbits, "kogut_susskind", lam=lam)})
    receipt = {
        "schema": SCHEMA,
        "scope": "finite_gauge_diagnostic",
        "physical_clay_receipt": False,
        "producer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "grid_scope": {
            "L": L_values,
            "wilson_diagonal_beta_s_eq_beta_t": betas,
            "kogut_susskind_lambda": lams,
            "universal_no_go": False,
        },
        "system": "Z2 lattice gauge theory, L x L periodic spatial torus, gauge-invariant sector",
        "receipt_under_test": (
            "finite ground-state-transform and cross-fiber receipt: "
            "U_r H_r U_r^{-1} = sum_C c_C (I - E_C) with c_C independent of the repaired value"
        ),
        "ground_state_transform": "Doob transform by the Perron vector, pi = Omega^2",
        "collars": "one collar per spatial link, fiber {o, X_l o}, pi-preserving heat bath",
        "runs": runs,
    }
    receipt["sha256_of_runs"] = hashlib.sha256(
        json.dumps(runs, sort_keys=True).encode("utf-8")
    ).hexdigest()
    return receipt


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--L", type=int, nargs="+", default=[2, 3])
    ap.add_argument("--beta", type=float, nargs="+", default=[0.1, 0.3, 0.5, 0.7, 1.0])
    ap.add_argument("--lam", type=float, nargs="+", default=[0.5, 1.0, 2.0])
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()
    receipt = run(args.L, args.beta, args.lam)
    text = json.dumps(receipt, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    for r in receipt["runs"]:
        tag = r["transfer"] + (" control" if "control" in r else "")
        print(
            f"L={r['L']} {tag:24s} params={r['parameters']} "
            f"fit_resid={r['constant_rate_fit']['relative_frobenius_residual']:.3e} "
            f"spread={r['fiber_dependent_rates']['spread_max_over_min']:.4f} "
            f"outside_single_flip={r['fiber_dependent_rates']['offdiagonal_mass_outside_single_flip']:.3e} "
            f"eta*={r['dobrushin']['eta_star']:.4f} "
            f"gapH={r['spectral']['gap_H']:.4f} gapHB={r['spectral']['gap_unit_rate_heat_bath']:.4f}"
        )


if __name__ == "__main__":
    main()
