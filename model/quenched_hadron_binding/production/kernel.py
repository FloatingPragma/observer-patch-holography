"""Compiled kernels of the production (optional, need numba): APE smearing, Wilson loops with off-axis separations,
Wilson-Dirac operator, even-odd BiCGStab, Gaussian source smearing, meson and nucleon contractions.

Without numba the module imports and AVAILABLE is False. The check and the tests of the package use lattice.py.
Conventions (see PROOF.md):
- site index s = i0 + N0 (i1 + N1 (i2 + N2 i3)), axis 3 = time; links U[s, mu] from s to s + e_mu, the layout of
  model/colour_deconfinement/su3_numba.py;
- spinor field psi[s, spin, colour]; Dirac basis, gamma4 = diag(1, 1, -1, -1);
- M = 1 - kappa H, H psi(x) = sum_mu [(1 - g_mu) U_mu(x) psi(x+mu) + (1 + g_mu) U_mu(x-mu)^dag psi(x-mu)];
- propagator P[3 beta + b, x, alpha, a] = S(x, alpha, a; source, beta, b).
"""
import itertools

import numpy as np

try:
    from numba import njit
    AVAILABLE = True
except Exception:                                  # pragma: no cover
    AVAILABLE = False

    def njit(*args, **kwargs):
        if len(args) == 1 and callable(args[0]) and not kwargs:
            return args[0]
        return lambda f: f

# ----------------------------------------------------------------------------------------------------------------
# gamma matrices (numpy)
# ----------------------------------------------------------------------------------------------------------------


def gammas():
    s1 = np.array([[0, 1], [1, 0]], complex)
    s2 = np.array([[0, -1j], [1j, 0]], complex)
    s3 = np.array([[1, 0], [0, -1]], complex)
    Z = np.zeros((2, 2), complex)
    I2 = np.eye(2, dtype=complex)
    g = [np.block([[Z, -1j * s], [1j * s, Z]]) for s in (s1, s2, s3)]
    g.append(np.block([[I2, Z], [Z, -I2]]))
    g5 = g[0] @ g[1] @ g[2] @ g[3]
    return g, g5


def charge_conjugation():
    """C with C g_mu C^-1 = -g_mu^T. In this basis C = g4 g2 (checked in the gates)."""
    g, g5 = gammas()
    return g[3] @ g[1]


def parity_sites(dims):
    N0, N1, N2, N3 = dims
    i0, i1, i2, i3 = np.meshgrid(np.arange(N0), np.arange(N1), np.arange(N2), np.arange(N3), indexing='ij')
    s = (i0 + N0 * (i1 + N1 * (i2 + N2 * i3))).ravel()
    p = ((i0 + i1 + i2 + i3) % 2).ravel()
    par = np.empty(s.size, np.int64)
    par[s] = p
    even = np.nonzero(par == 0)[0].astype(np.int64)
    odd = np.nonzero(par == 1)[0].astype(np.int64)
    return even, odd


def fermion_links(U, dims):
    """Copy of the links with the antiperiodic time boundary: sign on the time-like links of the last time slice."""
    N0, N1, N2, N3 = dims
    V3 = N0 * N1 * N2
    Uf = U.copy()
    Uf[(N3 - 1) * V3:, 3] *= -1.0
    return Uf


# ----------------------------------------------------------------------------------------------------------------
# Wilson hopping term
# ----------------------------------------------------------------------------------------------------------------


@njit(cache=True, fastmath=True)
def hop(U, up, dn, psi, out, sites):
    """out(x) = sum_mu [(1 - g_mu) U_mu(x) psi(x+mu) + (1 + g_mu) U_mu(x-mu)^dag psi(x-mu)] for x in sites."""
    h0 = np.empty(3, np.complex128)
    h1 = np.empty(3, np.complex128)
    g0 = np.empty(3, np.complex128)
    g1 = np.empty(3, np.complex128)
    for n in range(sites.shape[0]):
        x = sites[n]
        for s in range(4):
            for c in range(3):
                out[x, s, c] = 0.0
        # ---------------- mu = 0 forward: (1 - g1): h0 = p0 + i p3, h1 = p1 + i p2; rows 2, 3 = -i h1, -i h0
        y = up[x, 0]
        for c in range(3):
            a = psi[y, 0, c]
            b = psi[y, 3, c]
            h0[c] = complex(a.real - b.imag, a.imag + b.real)
            a = psi[y, 1, c]
            b = psi[y, 2, c]
            h1[c] = complex(a.real - b.imag, a.imag + b.real)
        for i in range(3):
            g0[i] = U[x, 0, i, 0] * h0[0] + U[x, 0, i, 1] * h0[1] + U[x, 0, i, 2] * h0[2]
            g1[i] = U[x, 0, i, 0] * h1[0] + U[x, 0, i, 1] * h1[1] + U[x, 0, i, 2] * h1[2]
        for c in range(3):
            out[x, 0, c] += g0[c]
            out[x, 1, c] += g1[c]
            out[x, 2, c] += complex(g1[c].imag, -g1[c].real)
            out[x, 3, c] += complex(g0[c].imag, -g0[c].real)
        # ---------------- mu = 0 backward: (1 + g1): h0 = p0 - i p3, h1 = p1 - i p2; rows 2, 3 = i h1, i h0
        y = dn[x, 0]
        for c in range(3):
            a = psi[y, 0, c]
            b = psi[y, 3, c]
            h0[c] = complex(a.real + b.imag, a.imag - b.real)
            a = psi[y, 1, c]
            b = psi[y, 2, c]
            h1[c] = complex(a.real + b.imag, a.imag - b.real)
        for i in range(3):
            g0[i] = (np.conj(U[y, 0, 0, i]) * h0[0] + np.conj(U[y, 0, 1, i]) * h0[1]
                     + np.conj(U[y, 0, 2, i]) * h0[2])
            g1[i] = (np.conj(U[y, 0, 0, i]) * h1[0] + np.conj(U[y, 0, 1, i]) * h1[1]
                     + np.conj(U[y, 0, 2, i]) * h1[2])
        for c in range(3):
            out[x, 0, c] += g0[c]
            out[x, 1, c] += g1[c]
            out[x, 2, c] += complex(-g1[c].imag, g1[c].real)
            out[x, 3, c] += complex(-g0[c].imag, g0[c].real)
        # ---------------- mu = 1 forward: (1 - g2): h0 = p0 + p3, h1 = p1 - p2; rows 2, 3 = -h1, h0
        y = up[x, 1]
        for c in range(3):
            h0[c] = psi[y, 0, c] + psi[y, 3, c]
            h1[c] = psi[y, 1, c] - psi[y, 2, c]
        for i in range(3):
            g0[i] = U[x, 1, i, 0] * h0[0] + U[x, 1, i, 1] * h0[1] + U[x, 1, i, 2] * h0[2]
            g1[i] = U[x, 1, i, 0] * h1[0] + U[x, 1, i, 1] * h1[1] + U[x, 1, i, 2] * h1[2]
        for c in range(3):
            out[x, 0, c] += g0[c]
            out[x, 1, c] += g1[c]
            out[x, 2, c] -= g1[c]
            out[x, 3, c] += g0[c]
        # ---------------- mu = 1 backward: (1 + g2): h0 = p0 - p3, h1 = p1 + p2; rows 2, 3 = h1, -h0
        y = dn[x, 1]
        for c in range(3):
            h0[c] = psi[y, 0, c] - psi[y, 3, c]
            h1[c] = psi[y, 1, c] + psi[y, 2, c]
        for i in range(3):
            g0[i] = (np.conj(U[y, 1, 0, i]) * h0[0] + np.conj(U[y, 1, 1, i]) * h0[1]
                     + np.conj(U[y, 1, 2, i]) * h0[2])
            g1[i] = (np.conj(U[y, 1, 0, i]) * h1[0] + np.conj(U[y, 1, 1, i]) * h1[1]
                     + np.conj(U[y, 1, 2, i]) * h1[2])
        for c in range(3):
            out[x, 0, c] += g0[c]
            out[x, 1, c] += g1[c]
            out[x, 2, c] += g1[c]
            out[x, 3, c] -= g0[c]
        # ---------------- mu = 2 forward: (1 - g3): h0 = p0 + i p2, h1 = p1 - i p3; rows 2, 3 = -i h0, i h1
        y = up[x, 2]
        for c in range(3):
            a = psi[y, 0, c]
            b = psi[y, 2, c]
            h0[c] = complex(a.real - b.imag, a.imag + b.real)
            a = psi[y, 1, c]
            b = psi[y, 3, c]
            h1[c] = complex(a.real + b.imag, a.imag - b.real)
        for i in range(3):
            g0[i] = U[x, 2, i, 0] * h0[0] + U[x, 2, i, 1] * h0[1] + U[x, 2, i, 2] * h0[2]
            g1[i] = U[x, 2, i, 0] * h1[0] + U[x, 2, i, 1] * h1[1] + U[x, 2, i, 2] * h1[2]
        for c in range(3):
            out[x, 0, c] += g0[c]
            out[x, 1, c] += g1[c]
            out[x, 2, c] += complex(g0[c].imag, -g0[c].real)
            out[x, 3, c] += complex(-g1[c].imag, g1[c].real)
        # ---------------- mu = 2 backward: (1 + g3): h0 = p0 - i p2, h1 = p1 + i p3; rows 2, 3 = i h0, -i h1
        y = dn[x, 2]
        for c in range(3):
            a = psi[y, 0, c]
            b = psi[y, 2, c]
            h0[c] = complex(a.real + b.imag, a.imag - b.real)
            a = psi[y, 1, c]
            b = psi[y, 3, c]
            h1[c] = complex(a.real - b.imag, a.imag + b.real)
        for i in range(3):
            g0[i] = (np.conj(U[y, 2, 0, i]) * h0[0] + np.conj(U[y, 2, 1, i]) * h0[1]
                     + np.conj(U[y, 2, 2, i]) * h0[2])
            g1[i] = (np.conj(U[y, 2, 0, i]) * h1[0] + np.conj(U[y, 2, 1, i]) * h1[1]
                     + np.conj(U[y, 2, 2, i]) * h1[2])
        for c in range(3):
            out[x, 0, c] += g0[c]
            out[x, 1, c] += g1[c]
            out[x, 2, c] += complex(-g0[c].imag, g0[c].real)
            out[x, 3, c] += complex(g1[c].imag, -g1[c].real)
        # ---------------- mu = 3 forward: (1 - g4) = diag(0, 0, 2, 2)
        y = up[x, 3]
        for i in range(3):
            out[x, 2, i] += 2.0 * (U[x, 3, i, 0] * psi[y, 2, 0] + U[x, 3, i, 1] * psi[y, 2, 1]
                                   + U[x, 3, i, 2] * psi[y, 2, 2])
            out[x, 3, i] += 2.0 * (U[x, 3, i, 0] * psi[y, 3, 0] + U[x, 3, i, 1] * psi[y, 3, 1]
                                   + U[x, 3, i, 2] * psi[y, 3, 2])
        # ---------------- mu = 3 backward: (1 + g4) = diag(2, 2, 0, 0)
        y = dn[x, 3]
        for i in range(3):
            out[x, 0, i] += 2.0 * (np.conj(U[y, 3, 0, i]) * psi[y, 0, 0] + np.conj(U[y, 3, 1, i]) * psi[y, 0, 1]
                                   + np.conj(U[y, 3, 2, i]) * psi[y, 0, 2])
            out[x, 1, i] += 2.0 * (np.conj(U[y, 3, 0, i]) * psi[y, 1, 0] + np.conj(U[y, 3, 1, i]) * psi[y, 1, 1]
                                   + np.conj(U[y, 3, 2, i]) * psi[y, 1, 2])


def hop_reference(U, up, dn, psi):
    """Slow numpy reference of the hopping term with explicit gamma matrices."""
    g, g5 = gammas()
    one = np.eye(4)
    out = np.zeros_like(psi)
    for mu in range(4):
        fwd = np.einsum('xab,xsb->xsa', U[:, mu], psi[up[:, mu]])
        out += np.einsum('ts,xsa->xta', one - g[mu], fwd)
        Ub = np.conj(np.swapaxes(U[dn[:, mu], mu], 1, 2))
        bwd = np.einsum('xab,xsb->xsa', Ub, psi[dn[:, mu]])
        out += np.einsum('ts,xsa->xta', one + g[mu], bwd)
    return out


@njit(cache=True, fastmath=True)
def apply_M(U, up, dn, kappa, psi, out, sites):
    """out = (1 - kappa H) psi on all sites (sites = arange(V))."""
    hop(U, up, dn, psi, out, sites)
    for n in range(sites.shape[0]):
        x = sites[n]
        for s in range(4):
            for c in range(3):
                out[x, s, c] = psi[x, s, c] - kappa * out[x, s, c]


@njit(cache=True, fastmath=True)
def gamma5_inplace(psi):
    """psi <- gamma5 psi; gamma5 = [[0, -1], [-1, 0]] in 2 x 2 blocks."""
    for x in range(psi.shape[0]):
        for c in range(3):
            a = psi[x, 0, c]
            b = psi[x, 1, c]
            psi[x, 0, c] = -psi[x, 2, c]
            psi[x, 1, c] = -psi[x, 3, c]
            psi[x, 2, c] = -a
            psi[x, 3, c] = -b


@njit(cache=True, fastmath=True)
def apply_Mdag(U, up, dn, kappa, psi, out, sites, work):
    """out = M^dag psi = g5 M g5 psi."""
    for x in range(psi.shape[0]):
        for s in range(4):
            for c in range(3):
                work[x, s, c] = psi[x, s, c]
    gamma5_inplace(work)
    apply_M(U, up, dn, kappa, work, out, sites)
    gamma5_inplace(out)


@njit(cache=True, fastmath=True)
def vdot_s(a, b, sites):
    z = 0j
    for n in range(sites.shape[0]):
        x = sites[n]
        for s in range(4):
            for c in range(3):
                z += np.conj(a[x, s, c]) * b[x, s, c]
    return z


@njit(cache=True, fastmath=True)
def vnorm2_s(a, sites):
    z = 0.0
    for n in range(sites.shape[0]):
        x = sites[n]
        for s in range(4):
            for c in range(3):
                z += a[x, s, c].real ** 2 + a[x, s, c].imag ** 2
    return z


@njit(cache=True, fastmath=True)
def eo_op(U, up, dn, kappa, v, out, tmp, even, odd):
    """out_o = (1 - kappa^2 H_oe H_eo) v_o. Reads v on odd sites, writes out on odd sites, tmp on even sites."""
    hop(U, up, dn, v, tmp, even)
    hop(U, up, dn, tmp, out, odd)
    k2 = kappa * kappa
    for n in range(odd.shape[0]):
        x = odd[n]
        for s in range(4):
            for c in range(3):
                out[x, s, c] = v[x, s, c] - k2 * out[x, s, c]


@njit(cache=True, fastmath=True)
def bicgstab_eo(U, up, dn, kappa, b, x, even, odd, tol, maxit):
    """Solve M x = b by BiCGStab on the even-odd reduced system. Returns (iterations, reduced relative residual)."""
    bp = np.zeros_like(b)
    tmp = np.zeros_like(b)
    hop(U, up, dn, b, bp, odd)
    for n in range(odd.shape[0]):
        y = odd[n]
        for s in range(4):
            for c in range(3):
                bp[y, s, c] = b[y, s, c] + kappa * bp[y, s, c]
    r = np.zeros_like(b)
    rh = np.zeros_like(b)
    p = np.zeros_like(b)
    v = np.zeros_like(b)
    sv = np.zeros_like(b)
    t = np.zeros_like(b)
    xo = np.zeros_like(b)
    for n in range(odd.shape[0]):
        y = odd[n]
        for s in range(4):
            for c in range(3):
                r[y, s, c] = bp[y, s, c]
                rh[y, s, c] = bp[y, s, c]
    bn = np.sqrt(vnorm2_s(bp, odd))
    res = 1.0
    it = 0
    if bn > 0.0:
        rho = 1.0 + 0j
        alpha = 1.0 + 0j
        omega = 1.0 + 0j
        while it < maxit:
            it += 1
            rho_new = vdot_s(rh, r, odd)
            beta = (rho_new / rho) * (alpha / omega)
            for n in range(odd.shape[0]):
                y = odd[n]
                for s in range(4):
                    for c in range(3):
                        p[y, s, c] = r[y, s, c] + beta * (p[y, s, c] - omega * v[y, s, c])
            eo_op(U, up, dn, kappa, p, v, tmp, even, odd)
            alpha = rho_new / vdot_s(rh, v, odd)
            for n in range(odd.shape[0]):
                y = odd[n]
                for s in range(4):
                    for c in range(3):
                        sv[y, s, c] = r[y, s, c] - alpha * v[y, s, c]
            sn = np.sqrt(vnorm2_s(sv, odd))
            if sn < tol * bn:
                for n in range(odd.shape[0]):
                    y = odd[n]
                    for s in range(4):
                        for c in range(3):
                            xo[y, s, c] += alpha * p[y, s, c]
                res = sn / bn
                break
            eo_op(U, up, dn, kappa, sv, t, tmp, even, odd)
            omega = vdot_s(t, sv, odd) / vnorm2_s(t, odd)
            for n in range(odd.shape[0]):
                y = odd[n]
                for s in range(4):
                    for c in range(3):
                        xo[y, s, c] += alpha * p[y, s, c] + omega * sv[y, s, c]
                        r[y, s, c] = sv[y, s, c] - omega * t[y, s, c]
            rho = rho_new
            res = np.sqrt(vnorm2_s(r, odd)) / bn
            if res < tol or not np.isfinite(res):
                break
    # reconstruct: x_e = b_e + kappa H_eo x_o
    hop(U, up, dn, xo, x, even)
    for n in range(even.shape[0]):
        y = even[n]
        for s in range(4):
            for c in range(3):
                x[y, s, c] = b[y, s, c] + kappa * x[y, s, c]
    for n in range(odd.shape[0]):
        y = odd[n]
        for s in range(4):
            for c in range(3):
                x[y, s, c] = xo[y, s, c]
    return it, res


@njit(cache=True, fastmath=True)
def true_residual(U, up, dn, kappa, x, b, sites):
    """|b - M x| / |b| with the plain operator."""
    w = np.zeros_like(b)
    apply_M(U, up, dn, kappa, x, w, sites)
    num = 0.0
    den = 0.0
    for n in range(sites.shape[0]):
        y = sites[n]
        for s in range(4):
            for c in range(3):
                d = b[y, s, c] - w[y, s, c]
                num += d.real ** 2 + d.imag ** 2
                den += b[y, s, c].real ** 2 + b[y, s, c].imag ** 2
    return np.sqrt(num / den)


@njit(cache=True, fastmath=True)
def cgnr(U, up, dn, kappa, b, x, sites, tol, maxit):
    """Solve M x = b by conjugate gradient on M^dag M x = M^dag b, start x = 0. Returns (iterations, |b - M x|/|b|)."""
    res = b.copy()
    sv = np.zeros_like(b)
    q = np.zeros_like(b)
    work = np.zeros_like(b)
    for y in range(x.shape[0]):
        for s in range(4):
            for c in range(3):
                x[y, s, c] = 0.0
    apply_Mdag(U, up, dn, kappa, res, sv, sites, work)
    p = sv.copy()
    gam = vnorm2_s(sv, sites)
    bn = np.sqrt(vnorm2_s(b, sites))
    it = 0
    rel = 1.0
    while it < maxit:
        it += 1
        apply_M(U, up, dn, kappa, p, q, sites)
        alpha = gam / vnorm2_s(q, sites)
        for y in range(x.shape[0]):
            for s in range(4):
                for c in range(3):
                    x[y, s, c] += alpha * p[y, s, c]
                    res[y, s, c] -= alpha * q[y, s, c]
        rel = np.sqrt(vnorm2_s(res, sites)) / bn
        if rel < tol:
            break
        apply_Mdag(U, up, dn, kappa, res, sv, sites, work)
        gnew = vnorm2_s(sv, sites)
        bet = gnew / gam
        gam = gnew
        for y in range(x.shape[0]):
            for s in range(4):
                for c in range(3):
                    p[y, s, c] = sv[y, s, c] + bet * p[y, s, c]
    return it, rel


def solve(Uf, up, dn, kappa, b, even, odd, allsites, tol=1e-10, maxit=4000, tol_true=1e-8):
    """Even-odd BiCGStab with a check of the true residual and a CGNR fallback. Returns x and a dict of statistics."""
    x = np.zeros_like(b)
    it, red = bicgstab_eo(Uf, up, dn, kappa, b, x, even, odd, tol, maxit)
    tr = true_residual(Uf, up, dn, kappa, x, b, allsites)
    info = dict(iterations=int(it), reduced_residual=float(red), true_residual=float(tr), fallback=False)
    if not np.isfinite(tr) or tr > tol_true:
        x2 = np.zeros_like(b)
        it2, rel = cgnr(Uf, up, dn, kappa, b, x2, allsites, tol_true * 1e-2, 40000)
        tr2 = true_residual(Uf, up, dn, kappa, x2, b, allsites)
        info.update(fallback=True, cg_iterations=int(it2), true_residual=float(tr2))
        x = x2
    return x, info


# ----------------------------------------------------------------------------------------------------------------
# APE smearing of spatial links
# ----------------------------------------------------------------------------------------------------------------


@njit(cache=True, fastmath=True)
def ape_step_raw(U, up, dn, alpha, out):
    """out[x, i] = (1 - alpha) U_i(x) + (alpha / 4) sum of the four spatial staples, i = 0, 1, 2. out[x, 3] = U_3(x)."""
    V = U.shape[0]
    t1 = np.empty((3, 3), np.complex128)
    t2 = np.empty((3, 3), np.complex128)
    c1 = 1.0 - alpha
    c2 = alpha / 4.0
    for x in range(V):
        for a in range(3):
            for b in range(3):
                out[x, 3, a, b] = U[x, 3, a, b]
        for i in range(3):
            for a in range(3):
                for b in range(3):
                    out[x, i, a, b] = c1 * U[x, i, a, b]
            xpi = up[x, i]
            for j in range(3):
                if j == i:
                    continue
                xpj = up[x, j]
                xmj = dn[x, j]
                xmjpi = up[xmj, i]
                # forward staple U_j(x) U_i(x+j) U_j(x+i)^dag
                for a in range(3):
                    for b in range(3):
                        t1[a, b] = (U[x, j, a, 0] * U[xpj, i, 0, b] + U[x, j, a, 1] * U[xpj, i, 1, b]
                                    + U[x, j, a, 2] * U[xpj, i, 2, b])
                for a in range(3):
                    for b in range(3):
                        out[x, i, a, b] += c2 * (t1[a, 0] * np.conj(U[xpi, j, b, 0]) + t1[a, 1] * np.conj(U[xpi, j, b, 1])
                                                 + t1[a, 2] * np.conj(U[xpi, j, b, 2]))
                # backward staple U_j(x-j)^dag U_i(x-j) U_j(x-j+i)
                for a in range(3):
                    for b in range(3):
                        t1[a, b] = (np.conj(U[xmj, j, 0, a]) * U[xmj, i, 0, b] + np.conj(U[xmj, j, 1, a]) * U[xmj, i, 1, b]
                                    + np.conj(U[xmj, j, 2, a]) * U[xmj, i, 2, b])
                for a in range(3):
                    for b in range(3):
                        out[x, i, a, b] += c2 * (t1[a, 0] * U[xmjpi, j, 0, b] + t1[a, 1] * U[xmjpi, j, 1, b]
                                                 + t1[a, 2] * U[xmjpi, j, 2, b])


def project_su3(W):
    """Polar projection to U(3) and division by the cube root of the determinant. W has shape (..., 3, 3)."""
    u, s, vh = np.linalg.svd(W)
    R = u @ vh
    det = np.linalg.det(R)
    return R * np.exp(-1j * np.angle(det) / 3.0)[..., None, None]


def ape_smear(U, up, dn, alpha=0.5, nsteps=12):
    """APE smearing of the spatial links. Time-like links are returned unchanged."""
    cur = U.copy()
    raw = np.empty_like(U)
    for _ in range(nsteps):
        ape_step_raw(cur, up, dn, alpha, raw)
        cur[:, :3] = project_su3(raw[:, :3])
    return cur


# ----------------------------------------------------------------------------------------------------------------
# Wilson loops with on-axis and off-axis spatial separations
# ----------------------------------------------------------------------------------------------------------------


def displacement_classes(rmax_axis, diag2=(1, 2, 3), diag3=(1, 2), with_210=True):
    """List of classes. Each class: dict(name, r, vec, orientations=[paths array (npaths, nsteps)]).
    Step codes 0, 1, 2 = +x, +y, +z and 3, 4, 5 = -x, -y, -z."""
    cl = []
    for r in range(1, rmax_axis + 1):
        cl.append(dict(name=f'({r},0,0)', r=float(r), vec=(r, 0, 0),
                       orientations=[np.array([[i] * r], np.int64) for i in range(3)]))
    for n in diag2:
        ors = []
        for i, j in ((0, 1), (0, 2), (1, 2)):
            for sj in (0, 3):
                a, b = i, j + sj
                ors.append(np.array([[a, b] * n, [b, a] * n], np.int64))
        cl.append(dict(name=f'({n},{n},0)', r=float(np.sqrt(2.0) * n), vec=(n, n, 0), orientations=ors))
    for n in diag3:
        ors = []
        for s1 in (0, 3):
            for s2 in (0, 3):
                steps = (0, 1 + s1, 2 + s2)
                ors.append(np.array([list(p) * n for p in itertools.permutations(steps)], np.int64))
        cl.append(dict(name=f'({n},{n},{n})', r=float(np.sqrt(3.0) * n), vec=(n, n, n), orientations=ors))
    if with_210:
        ors = []
        for i in range(3):
            for j in range(3):
                if i == j:
                    continue
                for sj in (0, 3):
                    a, b = i, j + sj
                    ors.append(np.array([[a, a, b], [a, b, a], [b, a, a]], np.int64))
        cl.append(dict(name='(2,1,0)', r=float(np.sqrt(5.0)), vec=(2, 1, 0), orientations=ors))
    cl.sort(key=lambda c: c['r'])
    return cl


@njit(cache=True, fastmath=True)
def transporter(Vs, up, dn, paths, S_out, end_out):
    """S_out[x] = mean over the paths of the ordered link product from x; end_out[x] = end site of the paths."""
    V = Vs.shape[0]
    npth, nst = paths.shape
    M = np.empty((3, 3), np.complex128)
    Q = np.empty((3, 3), np.complex128)
    for x in range(V):
        for a in range(3):
            for b in range(3):
                S_out[x, a, b] = 0.0
        y = x
        for p in range(npth):
            for a in range(3):
                for b in range(3):
                    M[a, b] = 1.0 if a == b else 0.0
            y = x
            for k in range(nst):
                d = paths[p, k]
                if d < 3:
                    for a in range(3):
                        for b in range(3):
                            Q[a, b] = M[a, 0] * Vs[y, d, 0, b] + M[a, 1] * Vs[y, d, 1, b] + M[a, 2] * Vs[y, d, 2, b]
                    y = up[y, d]
                else:
                    y = dn[y, d - 3]
                    for a in range(3):
                        for b in range(3):
                            Q[a, b] = (M[a, 0] * np.conj(Vs[y, d - 3, b, 0]) + M[a, 1] * np.conj(Vs[y, d - 3, b, 1])
                                       + M[a, 2] * np.conj(Vs[y, d - 3, b, 2]))
                for a in range(3):
                    for b in range(3):
                        M[a, b] = Q[a, b]
            for a in range(3):
                for b in range(3):
                    S_out[x, a, b] += M[a, b] / npth
        end_out[x] = y


@njit(cache=True, fastmath=True)
def time_lines(U, up, Tmax):
    """TL[T-1, x] = product of T time-like links from x; sh[T, x] = x + T e_3."""
    V = U.shape[0]
    TL = np.empty((Tmax, V, 3, 3), np.complex128)
    sh = np.empty((Tmax + 1, V), np.int64)
    for x in range(V):
        sh[0, x] = x
        sh[1, x] = up[x, 3]
        for a in range(3):
            for b in range(3):
                TL[0, x, a, b] = U[x, 3, a, b]
    for T in range(2, Tmax + 1):
        for x in range(V):
            e = sh[T - 1, x]
            sh[T, x] = up[e, 3]
            for a in range(3):
                for b in range(3):
                    TL[T - 1, x, a, b] = (TL[T - 2, x, a, 0] * U[e, 3, 0, b] + TL[T - 2, x, a, 1] * U[e, 3, 1, b]
                                          + TL[T - 2, x, a, 2] * U[e, 3, 2, b])
    return TL, sh


@njit(cache=True, fastmath=True)
def loops_from_transporter(Sx, endx, TL, sh, out):
    """out[T-1] = average over x of Re tr[ S(x) TL_T(end x) (TL_T(x) S(x + T))^dag ] / 3."""
    V = Sx.shape[0]
    Tmax = TL.shape[0]
    A = np.empty((3, 3), np.complex128)
    B = np.empty((3, 3), np.complex128)
    for T in range(1, Tmax + 1):
        tot = 0.0
        for x in range(V):
            e = endx[x]
            xT = sh[T, x]
            for a in range(3):
                for b in range(3):
                    A[a, b] = (Sx[x, a, 0] * TL[T - 1, e, 0, b] + Sx[x, a, 1] * TL[T - 1, e, 1, b]
                               + Sx[x, a, 2] * TL[T - 1, e, 2, b])
                    B[a, b] = (TL[T - 1, x, a, 0] * Sx[xT, 0, b] + TL[T - 1, x, a, 1] * Sx[xT, 1, b]
                               + TL[T - 1, x, a, 2] * Sx[xT, 2, b])
            for a in range(3):
                for b in range(3):
                    tot += (A[a, b] * np.conj(B[a, b])).real
        out[T - 1] = tot / (3.0 * V)


def wilson_loops_classes(U, Vs, up, dn, classes, Tmax):
    """W[class, T-1] with smeared spatial transporters (links Vs) and thin time-like lines (links U)."""
    V = U.shape[0]
    TL, sh = time_lines(U, up, Tmax)
    W = np.zeros((len(classes), Tmax))
    Sx = np.empty((V, 3, 3), np.complex128)
    endx = np.empty(V, np.int64)
    tmp = np.empty(Tmax)
    for k, c in enumerate(classes):
        for paths in c['orientations']:
            transporter(Vs, up, dn, paths, Sx, endx)
            loops_from_transporter(Sx, endx, TL, sh, tmp)
            W[k] += tmp
        W[k] /= len(c['orientations'])
    return W


# ----------------------------------------------------------------------------------------------------------------
# Gaussian smearing of quark fields
# ----------------------------------------------------------------------------------------------------------------


@njit(cache=True, fastmath=True)
def smear_field(Vs, up, dn, psi, s, n):
    """psi <- (1 + s^2 / (4 n) Lap)^n psi with the spatial covariant Laplacian of the links Vs. psi[x, spin, colour]."""
    V = psi.shape[0]
    w = s * s / (4.0 * n)
    tmp = np.empty_like(psi)
    for _ in range(n):
        for x in range(V):
            for sp in range(4):
                for a in range(3):
                    z = -6.0 * psi[x, sp, a]
                    for i in range(3):
                        y = up[x, i]
                        z += Vs[x, i, a, 0] * psi[y, sp, 0] + Vs[x, i, a, 1] * psi[y, sp, 1] + Vs[x, i, a, 2] * psi[y, sp, 2]
                        y = dn[x, i]
                        z += (np.conj(Vs[y, i, 0, a]) * psi[y, sp, 0] + np.conj(Vs[y, i, 1, a]) * psi[y, sp, 1]
                              + np.conj(Vs[y, i, 2, a]) * psi[y, sp, 2])
                    tmp[x, sp, a] = psi[x, sp, a] + w * z
        for x in range(V):
            for sp in range(4):
                for a in range(3):
                    psi[x, sp, a] = tmp[x, sp, a]


# ----------------------------------------------------------------------------------------------------------------
# Contractions
# ----------------------------------------------------------------------------------------------------------------


def meson_correlators(P, dims):
    """Pion and rho correlators from P[12, V, 4, 3]. Returns dict of real arrays of length N3."""
    N0, N1, N2, N3 = dims
    V3 = N0 * N1 * N2
    g, g5 = gammas()
    # S[x, alpha, a, beta, b]
    S = np.transpose(P.reshape(4, 3, N3 * V3, 4, 3), (2, 3, 4, 0, 1))
    out = {}
    out['pion'] = (np.abs(S) ** 2).reshape(N3, -1).sum(1)
    rho = np.zeros(N3)
    Sc = np.conj(S)
    for i in range(3):
        A = g5 @ g[i]
        Ad = A.conj().T
        AS = np.einsum('ij,xjakb->xiakb', A, S)
        ASA = np.einsum('xiakb,kl->xialb', AS, Ad)
        rho += (ASA * Sc).real.reshape(N3, -1).sum(1)
    out['rho'] = rho / 3.0
    return out


def eps_terms():
    """Nonzero entries of eps_abc eps_def as arrays (a, b, c, d, e, f, sign)."""
    perms = [((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((0, 2, 1), -1), ((2, 1, 0), -1), ((1, 0, 2), -1)]
    rows = []
    for (p, s1) in perms:
        for (q, s2) in perms:
            rows.append(list(p) + list(q) + [s1 * s2])
    return np.array(rows, np.int64)


@njit(cache=True, fastmath=True)
def nucleon_kernel(P, G, Gb, eps, N3, V3, cp, cm):
    """Nucleon two-point function, operator eps_abc (u_a^T G d_b) u_c, degenerate quarks.
    C(x) = eps eps { S^{cc'} tr[Q^{bb'} (S^{aa'})^T] + S^{cc'} (Q^{bb'})^T S^{aa'} }, Q = G S Gb.
    cp[t], cm[t]: traces with (1 + g4)/2 (spin rows 0, 1) and (1 - g4)/2 (spin rows 2, 3)."""
    S = np.empty((3, 3, 4, 4), np.complex128)
    Q = np.empty((3, 3, 4, 4), np.complex128)
    T = np.empty((4, 4), np.complex128)
    for t in range(N3):
        tp = 0j
        tm = 0j
        for xs in range(V3):
            x = t * V3 + xs
            for a in range(3):
                for b in range(3):
                    for al in range(4):
                        for be in range(4):
                            S[a, b, al, be] = P[3 * be + b, x, al, a]
            for a in range(3):
                for b in range(3):
                    for i in range(4):
                        for j in range(4):
                            z = 0j
                            for k in range(4):
                                z += G[i, k] * S[a, b, k, j]
                            T[i, j] = z
                    for i in range(4):
                        for j in range(4):
                            z = 0j
                            for k in range(4):
                                z += T[i, k] * Gb[k, j]
                            Q[a, b, i, j] = z
            for n in range(eps.shape[0]):
                a = eps[n, 0]
                b = eps[n, 1]
                c = eps[n, 2]
                d = eps[n, 3]
                e = eps[n, 4]
                f = eps[n, 5]
                sg = eps[n, 6]
                # term 1: tr[Q^{be} (S^{ad})^T] times diagonal of S^{cf}
                tr = 0j
                for i in range(4):
                    for j in range(4):
                        tr += Q[b, e, i, j] * S[a, d, i, j]
                d_up = S[c, f, 0, 0] + S[c, f, 1, 1]
                d_dn = S[c, f, 2, 2] + S[c, f, 3, 3]
                # term 2: diagonal of S^{cf} (Q^{be})^T S^{ad}
                for i in range(4):
                    for j in range(4):
                        z = 0j
                        for k in range(4):
                            z += Q[b, e, k, i] * S[a, d, k, j]
                        T[i, j] = z
                u2 = 0j
                l2 = 0j
                for gsp in range(2):
                    for i in range(4):
                        u2 += S[c, f, gsp, i] * T[i, gsp]
                for gsp in range(2, 4):
                    for i in range(4):
                        l2 += S[c, f, gsp, i] * T[i, gsp]
                tp += sg * (tr * d_up + u2)
                tm += sg * (tr * d_dn + l2)
        cp[t] = tp
        cm[t] = tm


def nucleon_site_reference(S, G, Gb):
    """Brute-force Wick sum at one site. S[a, alpha, a', alpha'] (sink colour, sink spin, source colour, source spin).
    Returns the 4 x 4 matrix C[gamma, gamma'] = eps eps G Gb S_d (S_u S_u - S_u S_u)."""
    eps = np.zeros((3, 3, 3))
    eps[0, 1, 2] = eps[1, 2, 0] = eps[2, 0, 1] = 1.0
    eps[0, 2, 1] = eps[2, 1, 0] = eps[1, 0, 2] = -1.0
    t1 = np.einsum('abc,def,ij,kl,bjek,aidl,cgfh->gh', eps, eps, G, Gb, S, S, S, optimize=True)
    t2 = np.einsum('abc,def,ij,kl,bjek,aifh,cgdl->gh', eps, eps, G, Gb, S, S, S, optimize=True)
    return t1 - t2


def nucleon_correlators(P, dims):
    """Returns (C_plus[t], C_minus[t]) complex arrays."""
    N0, N1, N2, N3 = dims
    V3 = N0 * N1 * N2
    g, g5 = gammas()
    C = charge_conjugation()
    G = C @ g5
    Gb = g[3] @ G.conj().T @ g[3]
    cp = np.zeros(N3, np.complex128)
    cm = np.zeros(N3, np.complex128)
    nucleon_kernel(np.ascontiguousarray(P), np.ascontiguousarray(G), np.ascontiguousarray(Gb), eps_terms(), N3, V3, cp, cm)
    return cp, cm


# ----------------------------------------------------------------------------------------------------------------
# Propagators
# ----------------------------------------------------------------------------------------------------------------


def propagator(Uf, Vs, up, dn, dims, kappa, even, odd, allsites, smear=None, tol=1e-10, maxit=4000, src=0):
    """12 solves from the site src. smear = (s, n) for a Gaussian source built with the links Vs, or None for a point.
    Returns P[12, V, 4, 3] and the list of solver statistics."""
    V = Uf.shape[0]
    P = np.empty((12, V, 4, 3), np.complex128)
    stats = []
    shapes = []
    for b in range(3):
        f = np.zeros((V, 4, 3), np.complex128)
        f[src, 0, b] = 1.0
        if smear is not None:
            smear_field(Vs, up, dn, f, smear[0], smear[1])
        shapes.append(f[:, 0, :].copy())
    for beta in range(4):
        for b in range(3):
            rhs = np.zeros((V, 4, 3), np.complex128)
            rhs[:, beta, :] = shapes[b]
            x, info = solve(Uf, up, dn, kappa, rhs, even, odd, allsites, tol=tol, maxit=maxit)
            P[3 * beta + b] = x
            stats.append(info)
    return P, stats


def smear_sink(P, Vs, up, dn, s, n):
    """Apply the Gaussian smearing to the sink index of every column of P."""
    out = np.empty_like(P)
    for j in range(P.shape[0]):
        f = P[j].copy()
        smear_field(Vs, up, dn, f, s, n)
        out[j] = f
    return out


# ----------------------------------------------------------------------------------------------------------------
# Free field (analytic reference)
# ----------------------------------------------------------------------------------------------------------------


def free_propagator_momentum(dims, kappa):
    """M^-1(p) for unit links, antiperiodic in time. Returns Sx[i0, i1, i2, i3, 4, 4] = S(x; 0) (spin matrix)."""
    N0, N1, N2, N3 = dims
    g, g5 = gammas()
    p0 = 2 * np.pi * np.arange(N0) / N0
    p1 = 2 * np.pi * np.arange(N1) / N1
    p2 = 2 * np.pi * np.arange(N2) / N2
    p3 = 2 * np.pi * (np.arange(N3) + 0.5) / N3
    P0, P1, P2, P3 = np.meshgrid(p0, p1, p2, p3, indexing='ij')
    A = 1.0 - 2.0 * kappa * (np.cos(P0) + np.cos(P1) + np.cos(P2) + np.cos(P3))
    B = [2.0 * kappa * np.sin(q) for q in (P0, P1, P2, P3)]
    den = A ** 2 + sum(b ** 2 for b in B)
    Sp = np.zeros(P0.shape + (4, 4), complex)
    Sp += (A / den)[..., None, None] * np.eye(4)
    for mu in range(4):
        Sp += (-1j * B[mu] / den)[..., None, None] * g[mu]
    # S(x) = (1/V) sum_p exp(i p x) S(p); the half-integer time momenta give a phase exp(i pi t / N3)
    Sx = np.fft.ifftn(Sp, axes=(0, 1, 2, 3))
    ph = np.exp(1j * np.pi * np.arange(N3) / N3)
    Sx = Sx * ph[None, None, None, :, None, None]
    return Sx
