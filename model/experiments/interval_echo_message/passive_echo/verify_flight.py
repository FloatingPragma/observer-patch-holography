"""Independent pair-array evolution and probability checks on the held-out station."""

import json
import math
from pathlib import Path

import numpy as np
from scipy.special import jv
from scipy.sparse import coo_matrix
from scipy.stats import binom

HERE = Path(__file__).parent


def propagate(edge_file, strength):
    with np.load(edge_file, allow_pickle=False) as data:
        a, b = data["edges"].T
    n = max(np.max(a), np.max(b)) + 1
    adj = coo_matrix(
        (np.ones(2 * len(a)), (np.r_[a, b], np.r_[b, a])), shape=(n, n)
    ).tocsr()
    contact = adj.toarray()
    deg = np.asarray(adj.sum(axis=0)).ravel()
    diag = strength * (deg[:, None] + deg[None, :] - contact) / 2 - 32 * contact

    def apply(psi):
        hop = adj @ psi + (adj @ psi.T).T + contact * psi.T
        out = ((diag + 10) * psi - strength * hop / 2) / 22
        np.fill_diagonal(out, 0)
        return out

    def step(psi):
        a, b = psi, apply(psi)
        out = jv(0, 22) * a - 2j * jv(1, 22) * b
        for k in range(2, 81):
            a, b = b, 2 * apply(b) - a
            out += (2 * (-1j) ** k * jv(k, 22)) * b
        return out

    # Gershgorin gives H in [-32,12]. The dropped scalar phase is common to
    # all arms of a given duration, including the two pieces of a no-mark arm.
    for source in (1, 182, 184):
        v = np.zeros(n, complex)
        v[source] = 1
        rows = [v]
        for k in range(8 if source == 1 else 4):
            psi = np.zeros((n, n), complex)
            psi[2, :] = rows[-1]
            whole = step(psi)
            if abs(np.vdot(whole, whole) - np.vdot(psi, psi)) > 2e-12:
                raise ValueError("one-interval norm")
            rows.append(whole[2, :].copy())
        yield source, rows


def check():
    rec = json.loads((HERE / "records.json").read_text())
    N = rec["shots"]
    edge = HERE.parents[2] / "executable_model/carriers/FEDERATION_LEVEL0.npz"
    tests = 0

    def compare(p, observed):
        nonlocal tests
        # Includes allowance for the producer's <=1e-15 per-bin rounding.
        lo = max(0.0, p - 1e-13)
        hi = min(1.0, p + 1e-13)
        assert (
            binom.ppf(0.01 / 20000, N, lo) <= observed <= binom.isf(0.01 / 20000, N, hi)
        )
        tests += 1

    for name, strength in (("native", 1), ("half_exchange", 0.5)):
        h = dict(propagate(edge, strength))
        for n in range(1, 5):
            b = rec["banks"][name]["182"][str(n)]
            for bit, source in enumerate((182, 184)):
                # Unconditioned first-leg probe, and both physically projected
                # second-leg arms using amplitudes independently reconstructed here.
                single = h[source][n]
                marked = h[1][n][source] * single
                missed = h[1][2 * n] - marked
                for vec, counts, tag in (
                    (single, b["one_way"][bit], ""),
                    (marked, b["echo"][bit], "marked/"),
                    (missed, b["echo"][bit], "unmarked/"),
                ):
                    mass = abs(vec) ** 2
                    for addr, p in enumerate(mass[:12]):
                        compare(p, counts[f"{tag}local/{addr}"])
                    compare(mass[12:].sum(), counts[tag + "quiet"])
    # |J_k(22)| <= 11^k/k! exp(121/(k+1)). Subsequent factorial ratios <=11/82.
    tail = 2 * 11**81 / math.factorial(81) * math.exp(121 / 82) / (1 - 11 / 82)
    print(
        f"{tests} held-out terminal counts agree with the pair-array calculation; "
        f"Chebyshev amplitude truncation <= {tail:.4g} per interval."
    )


if __name__ == "__main__":
    check()
