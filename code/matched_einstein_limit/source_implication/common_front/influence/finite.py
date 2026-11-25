"""Outward finite-frequency derivative at a unit-linewidth blank reference."""

import argparse
from math import isqrt
import json

import numpy as np
from mpmath import iv
from scipy import sparse
from scipy.sparse.linalg import splu
import sympy as sp

from . import measure, retarded

SCALE = 2**44


def kernels(h, zero, mark, scale=SCALE):
    # Complex shift has unit imaginary part. The exact inverse norm <= 1.
    n = h.shape[0]
    a = (zero * sparse.eye(n, dtype=np.int64) - h).tocsr()
    B = np.zeros((n, len(mark)), dtype=np.int64)
    B[mark] = np.eye(len(mark), dtype=np.int64)
    y = splu((a + 2j * sparse.eye(n)).tocsc()).solve(2 * B)
    u, w = [np.rint(x * scale).astype(np.int64) for x in (y.real, y.imag)]
    norm = int(abs(a).sum(axis=1).max())
    assert (norm + 2) * max(
        abs(u).max().item(), abs(w).max().item()
    ) + 2 * scale < 2**63
    ru = a @ u - 2 * w - 2 * scale * B
    rw = a @ w + 2 * u
    cap = max(abs(ru).max().item(), abs(rw).max().item())
    assert 2 * ru.size * cap**2 < 2**63
    total = int(np.square(ru).sum() + np.square(rw).sum())
    error = sp.Rational(isqrt(total) + 1, 2 * scale)
    assert error < sp.Rational("1e-8")
    return (sp.Matrix(u[mark]) + sp.I * sp.Matrix(w[mark])) / scale, error


def norm_bound(m):
    # Rational upper bound on the Frobenius norm avoids numerical eigenvalues.
    a = sp.re((m.H * m).trace())
    n = isqrt(int(sp.ceiling(a) * 10**8)) + 1
    return sp.Rational(n, 10**4)


def exterior_read(hs, zero, source, rho):
    # Direct sharp launch and return at the unchanged blank port. Centering
    # the common scalar is exact for trace-zero physical interventions.
    bmark = retarded.sites([5])
    rec = [kernels(h, zero, bmark, scale=2**52) for h in hs]
    err = max(x[1] for x in rec)
    scalar = sum(x[0].trace() for x in rec) / 8
    vals = measure.read([x[0] - scalar * sp.eye(4) for x in rec], source)
    rows = {}
    for k, cap in (("ordinary", 1000), ("pairing", 305), ("raw_o", 2), ("raw_p", 2)):
        w = iv.mpf([-1, 1]) * measure.action.front.full.ri(err * cap)
        rows[k] = [vals[k].real + w, vals[k].imag + w]
    out = {k: list(map(str, rows[k])) for k in ("ordinary", "pairing")}
    ratios = []
    for i in range(2):
        o, p, a, b = (rows[k][i] for k in ("ordinary", "pairing", "raw_o", "raw_p"))
        direct = p / o
        joint = 1 / rho + (b - a / rho) / o
        ratios.append(str(iv.mpf([max(direct.a, joint.a), min(direct.b, joint.b)])))
    out["collar_over_ordinary"] = ratios
    return out


def run():
    iv.dps = 32
    native = measure.action.native
    board = native.iso.pairing.formation.federation(1)
    route = native.iso.pairing.formation.region(board)
    mark = retarded.sites(route)
    source = measure.action.front.inputs()
    rho = measure.volume_ratio(source)
    out = {}
    ri = measure.action.front.full.ri
    for q in (-1, 1):
        hh = [native.generator(board, route, q, odd) for odd in (False, True)]
        zero = int(hh[0][20, 20])
        records = [kernels(h, zero, mark) for h in hh]
        error = max(v for _, v in records)
        # Qcenter is the exact reduced block at the rational head center and
        # dyadic pilot kernels. Every final solve reads the entire connected board.
        ks = [
            [
                [m[4 * i : 4 * i + 4, 4 * j : 4 * j + 4] for m, _ in records]
                for j in range(4)
            ]
            for i in range(4)
        ]
        Q = sp.Matrix(
            4, 4, lambda i, j: sum(sp.trace(e * x) for e, x in zip(ks[i][j], source[0]))
        ).applyfunc(sp.expand)
        inv = Q.inv().applyfunc(sp.expand)
        ell = norm_bound(inv)
        delta_q = sp.Rational("6e-7") + error
        delta_inv = ell * ell * delta_q / (1 - ell * delta_q)
        assert ell * delta_q < 1

        # Bilinear inverse weights are rounded outward to avoid carrying
        # unwieldy rational-decomposition coefficients into the source jets.
        def grid(z):
            return (int(sp.re(z) * SCALE) + sp.I * int(sp.im(z) * SCALE)) / SCALE

        inv0 = inv.applyfunc(grid)
        di = ri(delta_inv + sp.Rational(6, SCALE))
        L = ri(ell + sp.Rational(6, SCALE))
        D = {k: [] for k in ("ordinary", "pairing", "raw_o", "raw_p")}
        # Bound ||D Q|| by the entrywise, source-enclosed Frobenius norm.
        for row in ks:
            for elements in row:
                obs = measure.read(elements, source)
                for k in D:
                    D[k].append(abs(obs[k]).b)
        head_bound = {"ordinary": 1000, "pairing": 305, "raw_o": 2, "raw_p": 2}
        DK = {
            k: iv.sqrt(sum(u * u for u in samples)) + ri(4 * error * head_bound[k])
            for k, samples in D.items()
        }
        effect = [
            sum(
                (
                    inv0[3, i] * ks[i][j][p] * inv0[j, 0]
                    for i in range(4)
                    for j in range(4)
                ),
                sp.zeros(4),
            ).applyfunc(sp.expand)
            for p in (0, 1)
        ]
        nominal = measure.read(effect, source)
        for k in D:
            allowed = (2 * L + di) * di * DK[k] + L**2 * ri(error * head_bound[k])
            half = iv.mpf([-1, 1]) * allowed
            nominal[k] += iv.mpc(half, half)
        row = {}
        for comp in ("real", "imag"):
            v = {k: getattr(nominal[k], comp) for k in D}
            # The mean of this effect is not the self-energy; only its derivatives
            # enter the inverse derivative identity.
            ans = {k: str(v[k]) for k in ("ordinary", "pairing")}
            if 0 not in v["ordinary"]:
                direct = v["pairing"] / v["ordinary"]
                joint = 1 / rho + (v["raw_p"] - v["raw_o"] / rho) / v["ordinary"]
                ans["collar_over_ordinary"] = str(
                    iv.mpf([max(direct.a, joint.a), min(direct.b, joint.b)])
                )
            row[comp] = ans
        out[str(q)] = dict(
            exact_kernel_residual_error=str(error),
            dyadic_kernel_inverse_norm_upper=str(ell),
            true_inverse_radius=str(di),
            derivatives=row,
            exterior_5_return_derivative=exterior_read(hh, zero, source, rho),
        )
    return dict(blank_reference_port=5, omega=0, linewidth=1, charges=out)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    out = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (measure.HERE / "finite.json").write_text(out)
    else:
        print(out, end="")
