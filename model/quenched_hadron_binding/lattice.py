"""Quenched SU(3) lattice gauge theory with Wilson fermions in plain numpy.

Layout: links U[mu, x, y, z, t, 3, 3] from the site to its neighbour in direction mu; axis 3 is time.
Spinors psi[x, y, z, t, spin, colour]. Dirac basis with gamma4 = diag(1, 1, -1, -1).
Wilson-Dirac operator M = 1 - kappa H,
H psi(x) = sum_mu [(1 - g_mu) U_mu(x) psi(x + mu) + (1 + g_mu) U_mu(x - mu)^dag psi(x - mu)].
The production runs use a compiled kernel of the same operator (production/kernel.py, needs numba). This module is
the reference and runs the fast gates without numba.
"""
import itertools

import numpy as np


# ---------------------------------------------------------------------------------------------- gamma matrices
def gammas():
    s1 = np.array([[0, 1], [1, 0]], complex)
    s2 = np.array([[0, -1j], [1j, 0]], complex)
    s3 = np.array([[1, 0], [0, -1]], complex)
    zero = np.zeros((2, 2), complex)
    one = np.eye(2, dtype=complex)
    g = [np.block([[zero, -1j * s], [1j * s, zero]]) for s in (s1, s2, s3)]
    g.append(np.block([[one, zero], [zero, -one]]))
    return g, g[0] @ g[1] @ g[2] @ g[3]


def charge_conjugation():
    g, _ = gammas()
    return g[3] @ g[1]


_GAMMA, _GAMMA5 = gammas()
_PROJ_MINUS = [np.eye(4) - x for x in _GAMMA]
_PROJ_PLUS = [np.eye(4) + x for x in _GAMMA]


# ---------------------------------------------------------------------------------------------- SU(3) helpers
def dag(a):
    return np.conj(np.swapaxes(a, -1, -2))


def random_su3(shape, rng):
    z = rng.normal(size=tuple(shape) + (3, 3)) + 1j * rng.normal(size=tuple(shape) + (3, 3))
    q, r = np.linalg.qr(z)
    d = np.diagonal(r, axis1=-2, axis2=-1)
    q = q * (d / np.abs(d))[..., None, :]
    det = np.linalg.det(q)
    return q * np.exp(-1j * np.angle(det) / 3.0)[..., None, None]


def project_su3(w):
    u, _, vh = np.linalg.svd(w)
    r = u @ vh
    return r * np.exp(-1j * np.angle(np.linalg.det(r)) / 3.0)[..., None, None]


def cold_links(dims):
    u = np.zeros((4,) + tuple(dims) + (3, 3), complex)
    u[..., 0, 0] = u[..., 1, 1] = u[..., 2, 2] = 1.0
    return u


def hot_links(dims, rng):
    return random_su3((4,) + tuple(dims), rng)


def gauge_transform(u, g):
    """U_mu(x) -> g(x) U_mu(x) g(x + mu)^dag."""
    return np.stack([g @ u[mu] @ dag(np.roll(g, -1, axis=mu)) for mu in range(4)])


def staple_sum(u, mu):
    """A_mu(x) with local weight exp((beta / 3) Re tr(U_mu(x) A_mu(x)))."""
    a = np.zeros_like(u[mu])
    um = u[mu]
    for nu in range(4):
        if nu == mu:
            continue
        un = u[nu]
        a = a + np.roll(un, -1, axis=mu) @ dag(np.roll(um, -1, axis=nu)) @ dag(un)
        unm = np.roll(un, 1, axis=nu)
        a = a + dag(np.roll(unm, -1, axis=mu)) @ dag(np.roll(um, 1, axis=nu)) @ unm
    return a


def plaquette(u):
    """Average of Re tr U_P / 3 over all plaquettes."""
    tot = 0.0
    for mu in range(3):
        for nu in range(mu + 1, 4):
            p = u[mu] @ np.roll(u[nu], -1, axis=mu) @ dag(np.roll(u[mu], -1, axis=nu)) @ dag(u[nu])
            tot += np.trace(p, axis1=-2, axis2=-1).real.mean() / 3.0
    return tot / 6.0


def _x0_kennedy_pendleton(a, rng):
    """x0 in [-1, 1] with density sqrt(1 - x0^2) exp(a x0), elementwise."""
    x0 = np.empty(a.shape)
    todo = np.arange(a.size)
    af = a.ravel()
    xf = x0.ravel()
    while todo.size:
        aa = af[todo]
        r1 = 1.0 - rng.random(todo.size)
        r2 = rng.random(todo.size)
        r3 = 1.0 - rng.random(todo.size)
        delta = -(np.log(r1) + np.cos(2 * np.pi * r2) ** 2 * np.log(r3)) / aa
        ok = rng.random(todo.size) ** 2 <= 1.0 - 0.5 * delta
        xf[todo[ok]] = 1.0 - delta[ok]
        todo = todo[~ok]
    return x0


def _su2_left(m, i, j, a0, a1, a2, a3):
    """Rows i, j of m <- alpha rows, alpha = [[a0 + i a3, a2 + i a1], [-a2 + i a1, a0 - i a3]]."""
    ri = m[..., i, :].copy()
    rj = m[..., j, :].copy()
    m[..., i, :] = (a0 + 1j * a3)[..., None] * ri + (a2 + 1j * a1)[..., None] * rj
    m[..., j, :] = (-a2 + 1j * a1)[..., None] * ri + (a0 - 1j * a3)[..., None] * rj


def heatbath_sweep(u, beta, rng, overrelax=0):
    """One Cabibbo-Marinari heat-bath sweep over all links (checkerboard order), then `overrelax` microcanonical sweeps."""
    dims = u.shape[1:5]
    grid = np.indices(dims).sum(0) % 2
    for kind in ['heat'] + ['over'] * overrelax:
        for mu in range(4):
            for par in (0, 1):
                mask = grid == par
                w = (u[mu] @ staple_sum(u, mu))[mask]
                um = u[mu][mask]
                for i, j in ((0, 1), (0, 2), (1, 2)):
                    v0 = 0.5 * (w[:, i, i].real + w[:, j, j].real)
                    v3 = 0.5 * (w[:, i, i].imag - w[:, j, j].imag)
                    v1 = 0.5 * (w[:, i, j].imag + w[:, j, i].imag)
                    v2 = 0.5 * (w[:, i, j].real - w[:, j, i].real)
                    k = np.sqrt(v0 ** 2 + v1 ** 2 + v2 ** 2 + v3 ** 2)
                    v0, v1, v2, v3 = v0 / k, v1 / k, v2 / k, v3 / k
                    if kind == 'heat':
                        x0 = _x0_kennedy_pendleton(2.0 * beta * k / 3.0, rng)
                        r = np.sqrt(np.maximum(0.0, 1.0 - x0 ** 2))
                        ct = 2.0 * rng.random(k.size) - 1.0
                        st = np.sqrt(np.maximum(0.0, 1.0 - ct ** 2))
                        ph = 2 * np.pi * rng.random(k.size)
                        x1, x2, x3 = r * st * np.cos(ph), r * st * np.sin(ph), r * ct
                        a0 = x0 * v0 + x1 * v1 + x2 * v2 + x3 * v3
                        a1 = -x0 * v1 + v0 * x1 + (x2 * v3 - x3 * v2)
                        a2 = -x0 * v2 + v0 * x2 + (x3 * v1 - x1 * v3)
                        a3 = -x0 * v3 + v0 * x3 + (x1 * v2 - x2 * v1)
                    else:
                        a0 = v0 * v0 - v1 * v1 - v2 * v2 - v3 * v3
                        a1, a2, a3 = -2.0 * v0 * v1, -2.0 * v0 * v2, -2.0 * v0 * v3
                    _su2_left(um, i, j, a0, a1, a2, a3)
                    _su2_left(w, i, j, a0, a1, a2, a3)
                u[mu][mask] = um
    return u


# ---------------------------------------------------------------------------------------------- smearing and loops
def ape_smear(u, alpha=0.5, steps=12):
    """APE smearing of the spatial links; time-like links unchanged."""
    cur = u.copy()
    for _ in range(steps):
        new = cur.copy()
        for i in range(3):
            acc = (1.0 - alpha) * cur[i]
            for j in range(3):
                if j == i:
                    continue
                uj = cur[j]
                acc = acc + (alpha / 4.0) * (uj @ np.roll(cur[i], -1, axis=j) @ dag(np.roll(uj, -1, axis=i)))
                ujm = np.roll(uj, 1, axis=j)
                acc = acc + (alpha / 4.0) * (dag(ujm) @ np.roll(cur[i], 1, axis=j) @ np.roll(ujm, -1, axis=i))
            new[i] = project_su3(acc)
        cur = new
    return cur


def line(u, mu, n):
    """Product of n links in direction mu starting at every site."""
    out = u[mu].copy()
    for k in range(1, n):
        out = out @ np.roll(u[mu], -k, axis=mu)
    return out


def wilson_loop(u, v, r, t):
    """Average over sites and the three spatial directions of Re tr of the r x t loop / 3: spatial lines from the
    links v (smeared), time-like lines from the links u."""
    tot = 0.0
    tl = line(u, 3, t)
    for i in range(3):
        s = line(v, i, r)
        w = s @ np.roll(tl, -r, axis=i) @ dag(np.roll(s, -t, axis=3)) @ dag(tl)
        tot += np.trace(w, axis1=-2, axis2=-1).real.mean() / 3.0
    return tot / 3.0


# ---------------------------------------------------------------------------------------------- Wilson fermions
def fermion_links(u):
    """Antiperiodic time boundary: sign on the time-like links of the last time slice."""
    uf = u.copy()
    uf[3][:, :, :, -1] *= -1.0
    return uf


def hopping_blocks(uf):
    """The hopping term as site matrices on the 12 spin-colour components (index 3 spin + colour): for each direction
    the forward block (1 - g_mu) (x) U_mu(x) and the backward block (1 + g_mu) (x) U_mu(x - mu)^dag."""
    blocks = []
    for mu in range(4):
        shape = uf[mu].shape[:-2] + (12, 12)
        fwd = np.einsum('ts,...ab->...tasb', _PROJ_MINUS[mu], uf[mu]).reshape(shape)
        bwd = np.einsum('ts,...ab->...tasb', _PROJ_PLUS[mu], dag(np.roll(uf[mu], 1, axis=mu))).reshape(shape)
        blocks.append((fwd, bwd))
    return blocks


def hop(uf, psi, blocks=None):
    """H psi. A solver passes the site matrices of hopping_blocks(uf), built once for all its applications."""
    blocks = hopping_blocks(uf) if blocks is None else blocks
    v = psi.reshape(psi.shape[:-2] + (12, 1))
    out = np.zeros_like(v)
    for mu in range(4):
        out += blocks[mu][0] @ np.roll(v, -1, axis=mu)
        out += blocks[mu][1] @ np.roll(v, 1, axis=mu)
    return out.reshape(psi.shape)


def apply_m(uf, kappa, psi, blocks=None):
    return psi - kappa * hop(uf, psi, blocks)


def apply_mdag(uf, kappa, psi, blocks=None):
    return _GAMMA5 @ apply_m(uf, kappa, _GAMMA5 @ psi, blocks)


def solve_cg(uf, kappa, b, tol=1e-10, maxit=5000, blocks=None):
    """Solve M x = b by conjugate gradient on M^dag M x = M^dag b. Returns x, iterations, true relative residual."""
    blocks = hopping_blocks(uf) if blocks is None else blocks
    x = np.zeros_like(b)
    res = b.copy()
    s = apply_mdag(uf, kappa, res, blocks)
    p = s.copy()
    gam = np.vdot(s, s).real
    bn = np.linalg.norm(b)
    it = 0
    for it in range(1, maxit + 1):
        q = apply_m(uf, kappa, p, blocks)
        alpha = gam / np.vdot(q, q).real
        x = x + alpha * p
        res = res - alpha * q
        if np.linalg.norm(res) < tol * bn:
            break
        s = apply_mdag(uf, kappa, res, blocks)
        gnew = np.vdot(s, s).real
        p = s + (gnew / gam) * p
        gam = gnew
    return x, it, float(np.linalg.norm(b - apply_m(uf, kappa, x, blocks)) / bn)


def solve(uf, kappa, b, tol=1e-10, maxit=1000, blocks=None):
    """Solve M x = b by BiCGStab with the shadow vector M b, with conjugate gradient as the fallback.
    The shadow vector b itself breaks down for a point source: (1 - g_mu)(1 + g_mu) = 0 gives <b, H H b> = 0."""
    blocks = hopping_blocks(uf) if blocks is None else blocks
    x = np.zeros_like(b)
    r = b.copy()
    rh = apply_m(uf, kappa, b, blocks)
    p = np.zeros_like(b)
    v = np.zeros_like(b)
    rho = alpha = omega = 1.0 + 0j
    bn = np.linalg.norm(b)
    ok = False
    it = 0
    with np.errstate(all='ignore'):
        for it in range(1, maxit + 1):
            rho_new = np.vdot(rh, r)
            p = r + (rho_new / rho) * (alpha / omega) * (p - omega * v)
            v = apply_m(uf, kappa, p, blocks)
            alpha = rho_new / np.vdot(rh, v)
            s = r - alpha * v
            if not np.isfinite(alpha):
                break
            if np.linalg.norm(s) < tol * bn:
                x = x + alpha * p
                ok = True
                break
            t = apply_m(uf, kappa, s, blocks)
            omega = np.vdot(t, s) / np.vdot(t, t)
            x = x + alpha * p + omega * s
            r = s - omega * t
            rho = rho_new
            rn = np.linalg.norm(r)
            if not np.isfinite(rn):
                break
            if rn < tol * bn:
                ok = True
                break
    if ok:
        res = float(np.linalg.norm(b - apply_m(uf, kappa, x, blocks)) / bn)
        if np.isfinite(res) and res < 100 * tol:
            return x, it, res
    return solve_cg(uf, kappa, b, tol=tol, blocks=blocks)


def smear_field(v, psi, s=3.0, n=20):
    """(1 + s^2 / (4 n) Laplacian)^n with the spatial covariant Laplacian of the links v."""
    w = s * s / (4.0 * n)
    for _ in range(n):
        lap = -6.0 * psi
        for i in range(3):
            lap = lap + np.einsum('...ab,...sb->...sa', v[i], np.roll(psi, -1, axis=i))
            lap = lap + np.einsum('...ba,...sb->...sa', np.conj(np.roll(v[i], 1, axis=i)), np.roll(psi, 1, axis=i))
        psi = psi + w * lap
    return psi


def propagator(uf, kappa, v=None, smear=None, source=(0, 0, 0, 0), tol=1e-10):
    """S[x, y, z, t, alpha, a, beta, b] from the 12 sources at one site (Gaussian-smeared when smear = (s, n))."""
    dims = uf.shape[1:5]
    out = np.zeros(tuple(dims) + (4, 3, 4, 3), complex)
    blocks = hopping_blocks(uf)
    worst = 0.0
    for b in range(3):
        f = np.zeros(tuple(dims) + (4, 3), complex)
        f[source + (0, b)] = 1.0
        if smear is not None:
            f = smear_field(v, f, *smear)
        shape = f[..., 0, :]
        for beta in range(4):
            rhs = np.zeros(tuple(dims) + (4, 3), complex)
            rhs[..., beta, :] = shape
            x, _, res = solve(uf, kappa, rhs, tol=tol, blocks=blocks)
            worst = max(worst, res)
            out[..., beta, b] = x
    return out, worst


def free_propagator(dims, kappa):
    """S[x, y, z, t, alpha, beta] for unit links with the antiperiodic time boundary, from the momentum-space form
    M(p) = 1 - 2 kappa sum_mu cos p_mu + 2 i kappa sum_mu gamma_mu sin p_mu."""
    g, _ = gammas()
    ps = [2 * np.pi * np.arange(n) / n for n in dims[:3]] + [2 * np.pi * (np.arange(dims[3]) + 0.5) / dims[3]]
    grid = np.meshgrid(*ps, indexing='ij')
    a = 1.0 - 2.0 * kappa * sum(np.cos(q) for q in grid)
    b = [2.0 * kappa * np.sin(q) for q in grid]
    den = a ** 2 + sum(x ** 2 for x in b)
    sp = (a / den)[..., None, None] * np.eye(4) + sum((-1j * b[mu] / den)[..., None, None] * g[mu] for mu in range(4))
    sx = np.fft.ifftn(sp, axes=(0, 1, 2, 3))
    return sx * np.exp(1j * np.pi * np.arange(dims[3]) / dims[3])[None, None, None, :, None, None]


# ---------------------------------------------------------------------------------------------- correlators
def meson_correlators(s):
    """Pion and rho correlators at zero momentum from S[x, y, z, t, alpha, a, beta, b]."""
    g, g5 = gammas()
    pion = (np.abs(s) ** 2).sum(axis=(0, 1, 2, 4, 5, 6, 7))
    rho = 0.0
    for i in range(3):
        a = g5 @ g[i]
        asa = np.einsum('ij,...jakb,kl->...ialb', a, s, a.conj().T)
        rho = rho + (asa * np.conj(s)).real.sum(axis=(0, 1, 2, 4, 5, 6, 7)) / 3.0
    return pion, rho


def epsilon():
    e = np.zeros((3, 3, 3))
    for p in itertools.permutations(range(3)):
        e[p] = np.linalg.det(np.eye(3)[list(p)])
    return e


def nucleon_correlators(s):
    """Two-point function of eps_abc (u_a^T C g5 d_b) u_c with degenerate quarks, as the sum over the two Wick
    pairings of the up quarks. Returns the traces with (1 + g4)/2 and (1 - g4)/2 for every time slice."""
    g, g5 = gammas()
    big_g = charge_conjugation() @ g5
    big_gb = g[3] @ big_g.conj().T @ g[3]
    e = epsilon()
    # S[..., alpha, a, beta, b]: sink spin, sink colour, source spin, source colour
    q = np.einsum('ij,...jbke,kl->...bile', big_g, s, big_gb)       # (G S_d Gb)[b, i, e, l]
    # pairing 1: eps eps tr_spin[Q^{be} (S^{ad})^T] S^{cf}
    tr = np.einsum('...bile,...iald->...abde', q, s)
    k = np.einsum('abc,def,...abde->...cf', e, e, tr)
    t1 = np.einsum('...cf,...gchf->...gh', k, s)
    # pairing 2: eps eps S^{af}[i, h] Q^{be}[i, l] S^{cd}[g, l]
    w = np.einsum('abc,...bile->...acile', e, q)
    w = np.einsum('...acile,...iahf->...celhf', w, s)
    w = np.einsum('def,...celhf->...cdlh', e, w)
    t2 = np.einsum('...cdlh,...gcld->...gh', w, s)
    c = (t1 - t2).sum(axis=(0, 1, 2))
    return c[:, 0, 0] + c[:, 1, 1], c[:, 2, 2] + c[:, 3, 3]


def nucleon_site_reference(s):
    """The same two Wick pairings written as one sum over all indices, for one site: S[alpha, a, beta, b]."""
    g, g5 = gammas()
    big_g = charge_conjugation() @ g5
    big_gb = g[3] @ big_g.conj().T @ g[3]
    e = epsilon()
    t1 = np.einsum('abc,def,ij,kl,jbke,iald,gchf->gh', e, e, big_g, big_gb, s, s, s, optimize=True)
    t2 = np.einsum('abc,def,ij,kl,jbke,iahf,gcld->gh', e, e, big_g, big_gb, s, s, s, optimize=True)
    return t1 - t2
