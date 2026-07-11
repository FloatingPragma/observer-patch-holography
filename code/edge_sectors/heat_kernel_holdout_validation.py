#!/usr/bin/env python3
"""Held-out validation of the edge-sector heat-kernel law on finite gauge groups.

This script accompanies the "Numerical validation of the heat-kernel law"
section of the synthesis paper and implements the validation protocol
requested by paper-audit issue #540:

  * it separates symmetry-forced checks from overconstrained tests;
  * every substantive test fits the diffusion parameter t on a declared
    subset of spectral sectors and prints held-out residuals for the
    remaining sectors;
  * groups with degenerate nontrivial spectra (Z2, Z3) are labelled as
    implementation checks, because for them the multi-sector agreement is
    forced by charge conjugation plus eigenvalue degeneracy before any
    heat-kernel ansatz is imposed.

Models
------
Z_n  : 2x2 periodic lattice gauge theory (8 links) with Z_n link spaces and

           H = -K sum_p Re(B_p) - h sum_l Re(X_l) - Gamma sum_v Re(A_v),

       exactly the Hamiltonian displayed in the paper (K = 1, Gamma = 5).
       Region A consists of links whose tail has x = 0; the electric-center
       edge charge at a boundary vertex v is the restricted star
       Q_v = prod_{l in star(v) cap A} X_l^{+/-1}.

S3   : the exact single-plaquette reduction described in the paper.  With
       Gauss's law imposed at every vertex the physical space is the
       3-dimensional space of class functions of the plaquette holonomy,
       spanned by the normalized characters {chi_triv, chi_sign, chi_std}.
       The magnetic term is multiplication by Re chi_std(g); the electric
       term is diagonal in the character basis with the Cayley-graph
       Laplacian eigenvalues (transposition generating set)
       lambda_triv = 0, lambda_sign = 6, lambda_std = 3.

Extraction and held-out protocol
--------------------------------
The heat-kernel ansatz is p_R proportional to d_R exp(-t lambda_R).  For each
nontrivial sector, t_R = ln((p_0/d_0)/(p_R/d_R)) / lambda_R.  The fit sector
(declared below, always the lowest nonzero eigenvalue) determines t; every
other nontrivial sector's weight is then a parameter-free held-out
prediction, and the printed residual is

    residual(R) = ln(p_R_measured / p_R_predicted)  /  ln(p_R_predicted/p_0)

(a relative log-scale error), together with the eigenvalue-ratio diagnostic
log(p_R/p_0 d_R) / log(p_fit/p_0 d_fit) versus lambda_R / lambda_fit.

Only conventions internal to this script are used; overall normalizations of
the electric term rescale t but cancel in every ratio and residual.
"""

from __future__ import annotations

import argparse
import itertools
import math

import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh

K_PLAQ = 1.0
GAMMA = 5.0


# ----------------------------------------------------------------------
# Z_n on a 2x2 periodic lattice, full link basis, ground state by Lanczos.
# ----------------------------------------------------------------------

def _zn_links():
    """Link index table for the 2x2 torus: (x, y, dir) -> axis 0..7."""
    links = {}
    for x, y in itertools.product(range(2), range(2)):
        links[(x, y, "x")] = len(links)
        links[(x, y, "y")] = len(links)
    return links


def _zn_hamiltonian_apply(n, h):
    links = _zn_links()
    shape = (n,) * 8

    # Diagonal magnetic term: -K sum_p cos(2 pi (k1 + k2 - k3 - k4)/n).
    grids = np.meshgrid(*[np.arange(n)] * 8, indexing="ij", sparse=True)
    diag = np.zeros(shape)
    for x, y in itertools.product(range(2), range(2)):
        l1 = links[(x, y, "x")]
        l2 = links[((x + 1) % 2, y, "y")]
        l3 = links[(x, (y + 1) % 2, "x")]
        l4 = links[(x, y, "y")]
        phase = grids[l1] + grids[l2] - grids[l3] - grids[l4]
        diag = diag - K_PLAQ * np.cos(2.0 * np.pi * phase / n)

    # Gauss stars: A_v rolls outgoing links +1 and incoming links -1.
    stars = []
    for x, y in itertools.product(range(2), range(2)):
        out = [links[(x, y, "x")], links[(x, y, "y")]]
        inc = [links[((x - 1) % 2, y, "x")], links[(x, (y - 1) % 2, "y")]]
        stars.append((out, inc))

    def apply_h(vec):
        a = vec.reshape(shape)
        out = diag * a
        # Electric term: -h/2 (roll+1 + roll-1) per link.
        for ax in range(8):
            out = out - 0.5 * h * (np.roll(a, 1, axis=ax) + np.roll(a, -1, axis=ax))
        # Gauss term: -Gamma/2 (A_v + A_v^dagger).
        for out_axes, in_axes in stars:
            b = a
            for ax in out_axes:
                b = np.roll(b, 1, axis=ax)
            for ax in in_axes:
                b = np.roll(b, -1, axis=ax)
            c = a
            for ax in out_axes:
                c = np.roll(c, -1, axis=ax)
            for ax in in_axes:
                c = np.roll(c, 1, axis=ax)
            out = out - 0.5 * GAMMA * (b + c)
        return out.reshape(-1)

    return apply_h


def zn_edge_distribution(n, h):
    """Ground state of the 2x2 Z_n model and edge-charge distribution p_q.

    The distribution is measured at boundary vertex v = (0, 0) with the
    restricted star Q_v = X_(0,0,x) X_(0,0,y) X^dag_(0,1,y)  (outgoing X,
    incoming X^dagger, links in region A = {tail x = 0} only).
    """
    links = _zn_links()
    dim = n**8
    op = LinearOperator((dim, dim), matvec=_zn_hamiltonian_apply(n, h))
    # Ground state (smallest algebraic eigenvalue).
    _, vecs = eigsh(op, k=1, which="SA", maxiter=20000, tol=1e-12)
    psi = vecs[:, 0].reshape((n,) * 8)

    out_axes = [links[(0, 0, "x")], links[(0, 0, "y")]]
    in_axes = [links[(0, 1, "y")]]  # incoming y-link from (0, -1) = (0, 1)

    # <Q_v^k> for k = 0..n-1, then Fourier transform to p_q.
    expect = np.zeros(n, dtype=complex)
    for k in range(n):
        b = psi
        for ax in out_axes:
            b = np.roll(b, k, axis=ax)
        for ax in in_axes:
            b = np.roll(b, -k, axis=ax)
        expect[k] = np.vdot(psi, b)
    q = np.arange(n)
    omega = np.exp(-2j * np.pi * np.outer(q, q) / n)
    p = (omega @ expect).real / n
    p = np.clip(p, 0.0, None)
    return p / p.sum()


# ----------------------------------------------------------------------
# S3, exact single-plaquette class-function reduction.
# ----------------------------------------------------------------------

S3_IRREPS = ("triv", "sign", "std")
S3_DIMS = {"triv": 1, "sign": 1, "std": 2}
S3_LAMBDA = {"triv": 0.0, "sign": 6.0, "std": 3.0}
N_LINKS_S3 = 4  # single plaquette


def s3_edge_distribution(h):
    """Ground state of the reduced single-plaquette S3 model, p_R = |c_R|^2.

    Basis: normalized characters (chi_triv, chi_sign, chi_std).  Magnetic
    term = multiplication by Re chi_std(g); its matrix elements are the
    fusion multiplicities <chi_R, chi_std * chi_R'>:
    std x triv = std, std x sign = std, std x std = triv + sign + std.
    Electric term = h per link times the Cayley Laplacian, diagonal with
    eigenvalues (0, 6, 3).
    """
    magnetic = np.array(
        [
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
        ]
    )
    electric = np.diag([S3_LAMBDA[r] for r in S3_IRREPS])
    ham = -K_PLAQ * magnetic + h * N_LINKS_S3 * electric
    vals, vecs = np.linalg.eigh(ham)
    ground = vecs[:, np.argmin(vals)]
    p = ground**2
    return {r: p[i] for i, r in enumerate(S3_IRREPS)}


# ----------------------------------------------------------------------
# Held-out fitting.
# ----------------------------------------------------------------------

def fit_t(p0, p_fit, d_fit, lam_fit):
    return math.log((p0 / 1.0) / (p_fit / d_fit)) / lam_fit


def predict(p0, d, lam, t):
    return p0 * d * math.exp(-t * lam)


def report_zn(n, h_values):
    lam = [4.0 * math.sin(math.pi * q / n) ** 2 for q in range(n)]
    distinct = len({round(l, 12) for l in lam[1:]})
    if distinct < 2:
        kind = "SYMMETRY-FORCED IMPLEMENTATION CHECK"
    else:
        kind = "OVERCONSTRAINED HELD-OUT TEST"
    print(f"\n=== Z_{n} ({kind}) ===")
    print(f"eigenvalues: {[round(l, 6) for l in lam]}")
    if n == 3:
        print(
            "note: lambda_1 = lambda_2 and charge conjugation force "
            "p_1 = p_2 and t_{q=1} = t_{q=2} for ANY conjugation-invariant "
            "distribution; agreement below tests the implementation only."
        )
    header = "h      " + "  ".join(f"p_{q:<8d}" for q in range(n)) + "  t(fit q=1)  held-out residuals"
    print(header)
    for h in h_values:
        p = zn_edge_distribution(n, h)
        t = fit_t(p[0], p[1], 1.0, lam[1])
        residuals = []
        for q in range(2, (n // 2) + 1):
            pred = predict(p[0], 1.0, lam[q], t)
            res = math.log(p[q] / pred) / abs(math.log(pred / p[0]))
            ratio = math.log(p[q] / p[0]) / math.log(p[1] / p[0])
            residuals.append(
                f"q={q}: pred {pred:.3e} meas {p[q]:.3e} "
                f"res {res:+.2%} ratio {ratio:.4f} (target {lam[q]/lam[1]:.4f})"
            )
        row = f"{h:<6.2f} " + "  ".join(f"{p[q]:.4e}" for q in range(n)) + f"  {t:<10.4f}"
        print(row)
        for r in residuals:
            print(f"        {r}")


def report_s3(h_values):
    print("\n=== S_3 (OVERCONSTRAINED HELD-OUT TEST, nonabelian) ===")
    print("eigenvalues: triv 0, sign 6, std 3 (distinct nonzero pair)")
    print("fit sector: std; held-out sector: sign; target log-ratio = 2")
    print("h      p_triv      p_sign      p_std       t(std)   pred p_sign  res      log-ratio")
    for h in h_values:
        p = s3_edge_distribution(h)
        t = fit_t(p["triv"], p["std"], S3_DIMS["std"], S3_LAMBDA["std"])
        pred_sign = predict(p["triv"], S3_DIMS["sign"], S3_LAMBDA["sign"], t)
        res = math.log(p["sign"] / pred_sign) / abs(math.log(pred_sign / p["triv"]))
        log_ratio = math.log(p["sign"] / p["triv"]) / math.log(p["std"] / (2 * p["triv"]))
        print(
            f"{h:<6.2f} {p['triv']:.4e}  {p['sign']:.4e}  {p['std']:.4e}  "
            f"{t:<8.4f} {pred_sign:.3e}    {res:+.2%}  {log_ratio:.4f}"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--groups",
        nargs="+",
        default=["Z2", "Z3", "Z5", "S3"],
        choices=["Z2", "Z3", "Z5", "S3"],
    )
    parser.add_argument(
        "--h", nargs="+", type=float, default=None, help="electric couplings to scan"
    )
    args = parser.parse_args()

    for group in args.groups:
        if group == "S3":
            report_s3(args.h or [0.5, 1.0, 2.0, 5.0, 12.0, 100.0])
        else:
            n = int(group[1:])
            default_h = [0.05, 0.1, 0.2, 0.5] if n == 5 else [0.2, 0.5, 1.0, 2.0]
            report_zn(n, args.h or default_h)

    print(
        "\nSummary: Z2/Z3 rows are implementation checks (single or degenerate "
        "nontrivial eigenvalue).  Z5 and S3 fit t on one spectral sector and "
        "print held-out residuals for a second, distinct, nonzero eigenvalue "
        "sector.  The residuals converge in the declared regime of each model: "
        "h -> 0 for the Z5 torus model (weak field; at large h the ratio "
        "drifts toward the perturbative order-counting value 2), and large h "
        "for the single-plaquette S3 model."
    )


if __name__ == "__main__":
    main()
