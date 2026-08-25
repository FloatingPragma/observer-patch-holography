"""Numeric companion to ``Lean/Screen/GaugeOrbitQuotientGap.lean``.

Checks the worked two-fiber quotient instance against the committed
conventions of ``z2_finite_transfer_receipt.py``:

* the producer's ``star_masks`` at ``L = 1`` evaluates to one mask
  flipping both links, and ``gauge_group_masks`` to the order-two group
  ``{0, both}``, the generating data of the Lean instance;
* on the four-configuration product space with symmetric weights
  ``(1/2, 1/2)`` and one constant rate per fiber, the generator
  ``sum_l c_l (I - E_l)`` restricted to gauge-invariant mean-zero
  observables acts as the scalar ``c_0 + c_1``, so the exact quotient
  Dirichlet constant is ``c_0 + c_1``;
* the product-space Dirichlet constant on all mean-zero observables is
  ``min(c_0, c_1)``, so the minimum-rate lower bound holds on the
  quotient and is strict for positive rates.

The script asserts every claim and prints a JSON summary.  It writes no
receipt and registers nothing.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from z2_finite_transfer_receipt import gauge_group_masks, star_masks


def fiber_rate(lam: float, r: float) -> float:
    """Committed Kogut--Susskind collar rate ``lam * (r + 1/r)``."""
    return lam * (r + 1.0 / r)


def build_instance(c0: float, c1: float) -> dict[str, float]:
    n_links = 2
    n_cfg = 1 << n_links  # configurations as bitmasks 0..3

    stars = star_masks(1)
    assert stars == [0b11], f"star_masks(1) = {stars}, expected [0b11]"
    group = sorted(int(g) for g in gauge_group_masks(1))
    assert group == [0b00, 0b11], f"gauge_group_masks(1) = {group}"

    w = 0.5  # symmetric fiber weights, the descent condition on moved links
    pi = np.full(n_cfg, w * w)

    def cond_exp(link: int) -> np.ndarray:
        e = np.zeros((n_cfg, n_cfg))
        for x in range(n_cfg):
            for b in (0, 1):
                y = (x & ~(1 << link)) | (b << link)
                e[x, y] += w
        return e

    identity = np.eye(n_cfg)
    gen = c0 * (identity - cond_exp(0)) + c1 * (identity - cond_exp(1))

    # gauge translation: XOR with the star mask
    tau = np.zeros((n_cfg, n_cfg))
    for x in range(n_cfg):
        tau[x, x ^ 0b11] = 1.0

    # equivariance of the generator under the gauge translation
    assert np.allclose(tau @ gen, gen @ tau)

    # invariant subspace: eigenvectors of tau at eigenvalue +1
    # basis: indicator of orbit {00,11} and of orbit {01,10}
    v_diag = np.array([1.0, 0.0, 0.0, 1.0])
    v_off = np.array([0.0, 1.0, 1.0, 0.0])
    for v in (v_diag, v_off):
        assert np.allclose(tau @ v, v)

    # mean-zero invariant line: the Lean diagWitness
    witness = np.array([1.0, -1.0, -1.0, 1.0])
    assert abs(pi @ witness) < 1e-15
    assert np.allclose(tau @ witness, witness)

    # the generator acts on it as the scalar c0 + c1
    assert np.allclose(gen @ witness, (c0 + c1) * witness)

    def dirichlet(v: np.ndarray) -> float:
        return float((pi * v) @ (gen @ v)) / float((pi * v) @ v)

    quotient_gap = dirichlet(witness)
    assert abs(quotient_gap - (c0 + c1)) < 1e-12

    # product-space gap on all mean-zero observables: min fiber rate
    # (Walsh eigenvalues 0, c0, c1, c0+c1; drop the constant mode)
    d = np.diag(np.sqrt(pi))
    sym = d @ gen @ np.linalg.inv(d)
    eig = np.sort(np.linalg.eigvalsh((sym + sym.T) / 2.0))
    assert abs(eig[0]) < 1e-12
    product_gap = float(eig[1])
    assert abs(product_gap - min(c0, c1)) < 1e-12
    assert min(c0, c1) <= quotient_gap + 1e-12

    return {
        "c0": c0,
        "c1": c1,
        "product_gap": product_gap,
        "quotient_gap": quotient_gap,
        "min_rate_bound": min(c0, c1),
    }


def main() -> None:
    lam = 1.0
    r_moved = 1.0  # symmetric weights pin the moved-link ratio to one
    rate = fiber_rate(lam, r_moved)
    assert rate == 2.0 * lam

    committed = build_instance(rate, rate)
    assert committed["quotient_gap"] == 4.0
    generic = build_instance(2.5, 10.0 / 3.0)

    summary = {
        "schema": "oph.yang_mills.gauge_orbit_quotient_instances.v1",
        "lean_module": "Lean/Screen/GaugeOrbitQuotientGap.lean",
        "star_masks_L1": [3],
        "gauge_group_masks_L1": [0, 3],
        "committed_rates_lambda1": committed,
        "generic_rates": generic,
        "checks": "all assertions passed",
    }
    json.dump(summary, sys.stdout, indent=2)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
