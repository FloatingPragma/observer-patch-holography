"""Compiled update kernels for the three-dimensional SU(2) + Higgs lattice theory (optional, needs numba).

The production runs use these kernels. Without numba the package falls back to the numpy implementation in
`lattice.py`, which samples the same action. Storage here: U[3, N, N, N, 4] and phi[N, N, N, 4] with
q = (q0, q1, q2, q3) for the matrix q0 + i (q1 sigma_1 + q2 sigma_2 + q3 sigma_3).
"""
import math

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


# ----------------------------------------------------------------------------------------------
# quaternion helpers
# ----------------------------------------------------------------------------------------------
@njit(inline='always')
def qmul(a0, a1, a2, a3, b0, b1, b2, b3):
    return (a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3,
            a0 * b1 + b0 * a1 - (a2 * b3 - a3 * b2),
            a0 * b2 + b0 * a2 - (a3 * b1 - a1 * b3),
            a0 * b3 + b0 * a3 - (a1 * b2 - a2 * b1))


@njit(inline='always')
def nb(L, c0, c1, c2, mu, s):
    """Coordinates of the neighbour of (c0, c1, c2) in direction mu, step s = +1 or -1."""
    if mu == 0:
        return (c0 + s) % L, c1, c2
    elif mu == 1:
        return c0, (c1 + s) % L, c2
    else:
        return c0, c1, (c2 + s) % L


@njit(cache=True)
def staple(U, L, mu, c0, c1, c2):
    """Sum over the four staples of link (mu, c): the gauge action is -beta_G 1/2 Tr[U A] + const."""
    s0 = 0.0; s1 = 0.0; s2 = 0.0; s3 = 0.0
    f0, f1, f2 = nb(L, c0, c1, c2, mu, 1)
    for nu in range(3):
        if nu == mu:
            continue
        # upper: U_nu(x+mu) U_mu^dag(x+nu) U_nu^dag(x)
        g0, g1, g2 = nb(L, c0, c1, c2, nu, 1)
        a = U[nu, f0, f1, f2]
        b = U[mu, g0, g1, g2]
        c = U[nu, c0, c1, c2]
        t0, t1, t2, t3 = qmul(a[0], a[1], a[2], a[3], b[0], -b[1], -b[2], -b[3])
        t0, t1, t2, t3 = qmul(t0, t1, t2, t3, c[0], -c[1], -c[2], -c[3])
        s0 += t0; s1 += t1; s2 += t2; s3 += t3
        # lower: U_nu^dag(x+mu-nu) U_mu^dag(x-nu) U_nu(x-nu)
        h0, h1, h2 = nb(L, c0, c1, c2, nu, -1)
        k0, k1, k2 = nb(L, f0, f1, f2, nu, -1)
        a = U[nu, k0, k1, k2]
        b = U[mu, h0, h1, h2]
        c = U[nu, h0, h1, h2]
        t0, t1, t2, t3 = qmul(a[0], -a[1], -a[2], -a[3], b[0], -b[1], -b[2], -b[3])
        t0, t1, t2, t3 = qmul(t0, t1, t2, t3, c[0], c[1], c[2], c[3])
        s0 += t0; s1 += t1; s2 += t2; s3 += t3
    return s0, s1, s2, s3


@njit(cache=True)
def link_force(U, phi, L, mu, c0, c1, c2, bG, bH):
    """W with local weight exp(1/2 Tr[U W]) = exp(dot(U, W^dag))."""
    s0, s1, s2, s3 = staple(U, L, mu, c0, c1, c2)
    w0 = bG * s0; w1 = bG * s1; w2 = bG * s2; w3 = bG * s3
    if bH != 0.0:
        f0, f1, f2 = nb(L, c0, c1, c2, mu, 1)
        p = phi[f0, f1, f2]
        q = phi[c0, c1, c2]
        t0, t1, t2, t3 = qmul(p[0], p[1], p[2], p[3], q[0], -q[1], -q[2], -q[3])
        w0 += bH * t0; w1 += bH * t1; w2 += bH * t2; w3 += bH * t3
    return w0, w1, w2, w3


@njit(cache=True)
def sweep_links_hb(U, phi, L, bG, bH):
    """Kennedy-Pendleton heat bath on every link, checkerboard order."""
    for mu in range(3):
        for par in range(2):
            for c0 in range(L):
                for c1 in range(L):
                    for c2 in range((c0 + c1 + par) % 2, L, 2):
                        w0, w1, w2, w3 = link_force(U, phi, L, mu, c0, c1, c2, bG, bH)
                        k = math.sqrt(w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3)
                        if k < 1e-12:
                            # Haar random link
                            while True:
                                r0 = 2 * np.random.random() - 1; r1 = 2 * np.random.random() - 1
                                r2 = 2 * np.random.random() - 1; r3 = 2 * np.random.random() - 1
                                n = r0 * r0 + r1 * r1 + r2 * r2 + r3 * r3
                                if n <= 1.0 and n > 1e-12:
                                    break
                            n = math.sqrt(n)
                            U[mu, c0, c1, c2, 0] = r0 / n; U[mu, c0, c1, c2, 1] = r1 / n
                            U[mu, c0, c1, c2, 2] = r2 / n; U[mu, c0, c1, c2, 3] = r3 / n
                            continue
                        # X = U V0 has density exp(k x0) sqrt(1 - x0^2)
                        if k > 1.6:
                            while True:
                                r1 = 1.0 - np.random.random()
                                r2 = np.random.random()
                                r3 = 1.0 - np.random.random()
                                cc = math.cos(2.0 * math.pi * r2)
                                lam2 = -(math.log(r1) + cc * cc * math.log(r3)) / (2.0 * k)
                                r4 = np.random.random()
                                if r4 * r4 <= 1.0 - lam2:
                                    break
                            x0 = 1.0 - 2.0 * lam2
                        else:
                            # Creutz: x0 from exp(k x0) on [-1, 1], accept with sqrt(1 - x0^2)
                            em = math.exp(-2.0 * k)
                            while True:
                                r1 = np.random.random()
                                x0 = 1.0 + math.log(em + (1.0 - em) * r1) / k
                                r2 = np.random.random()
                                if r2 * r2 <= 1.0 - x0 * x0:
                                    break
                        rad = math.sqrt(max(0.0, 1.0 - x0 * x0))
                        ct = 2.0 * np.random.random() - 1.0
                        st = math.sqrt(max(0.0, 1.0 - ct * ct))
                        ph = 2.0 * math.pi * np.random.random()
                        x1 = rad * st * math.cos(ph); x2 = rad * st * math.sin(ph); x3 = rad * ct
                        # U = X V0^dag, V0 = W / k
                        u0, u1, u2, u3 = qmul(x0, x1, x2, x3, w0 / k, -w1 / k, -w2 / k, -w3 / k)
                        n = math.sqrt(u0 * u0 + u1 * u1 + u2 * u2 + u3 * u3)
                        U[mu, c0, c1, c2, 0] = u0 / n; U[mu, c0, c1, c2, 1] = u1 / n
                        U[mu, c0, c1, c2, 2] = u2 / n; U[mu, c0, c1, c2, 3] = u3 / n


@njit(cache=True)
def sweep_links_or(U, phi, L, bG, bH):
    """Overrelaxation U -> V0^dag U^dag V0^dag (action preserving), checkerboard order."""
    for mu in range(3):
        for par in range(2):
            for c0 in range(L):
                for c1 in range(L):
                    for c2 in range((c0 + c1 + par) % 2, L, 2):
                        w0, w1, w2, w3 = link_force(U, phi, L, mu, c0, c1, c2, bG, bH)
                        k = math.sqrt(w0 * w0 + w1 * w1 + w2 * w2 + w3 * w3)
                        if k < 1e-12:
                            continue
                        v0 = w0 / k; v1 = -w1 / k; v2 = -w2 / k; v3 = -w3 / k    # V0^dag
                        u = U[mu, c0, c1, c2]
                        t0, t1, t2, t3 = qmul(v0, v1, v2, v3, u[0], -u[1], -u[2], -u[3])
                        t0, t1, t2, t3 = qmul(t0, t1, t2, t3, v0, v1, v2, v3)
                        n = math.sqrt(t0 * t0 + t1 * t1 + t2 * t2 + t3 * t3)
                        U[mu, c0, c1, c2, 0] = t0 / n; U[mu, c0, c1, c2, 1] = t1 / n
                        U[mu, c0, c1, c2, 2] = t2 / n; U[mu, c0, c1, c2, 3] = t3 / n


@njit(cache=True)
def higgs_force(U, phi, L, c0, c1, c2):
    """b = sum_i [U_i(x) Phi(x+i) + U_i^dag(x-i) Phi(x-i)]; hopping part of the action is -beta_H dot(phi, b)."""
    s0 = 0.0; s1 = 0.0; s2 = 0.0; s3 = 0.0
    for mu in range(3):
        f0, f1, f2 = nb(L, c0, c1, c2, mu, 1)
        u = U[mu, c0, c1, c2]
        p = phi[f0, f1, f2]
        t0, t1, t2, t3 = qmul(u[0], u[1], u[2], u[3], p[0], p[1], p[2], p[3])
        s0 += t0; s1 += t1; s2 += t2; s3 += t3
        g0, g1, g2 = nb(L, c0, c1, c2, mu, -1)
        u = U[mu, g0, g1, g2]
        p = phi[g0, g1, g2]
        t0, t1, t2, t3 = qmul(u[0], -u[1], -u[2], -u[3], p[0], p[1], p[2], p[3])
        s0 += t0; s1 += t1; s2 += t2; s3 += t3
    return s0, s1, s2, s3


@njit(cache=True)
def sweep_higgs_hb(U, phi, L, bH, bR, ntry):
    """Heat bath for the Higgs field: Gaussian proposal with accept/reject for the quartic term.
    Returns the number of accepted sites."""
    nacc = 0
    for par in range(2):
        for c0 in range(L):
            for c1 in range(L):
                for c2 in range((c0 + c1 + par) % 2, L, 2):
                    b0, b1, b2, b3 = higgs_force(U, phi, L, c0, c1, c2)
                    F0 = bH * b0; F1 = bH * b1; F2 = bH * b2; F3 = bH * b3
                    FF = F0 * F0 + F1 * F1 + F2 * F2 + F3 * F3
                    if bR == 0.0:
                        al = 1.0
                        v = 0.0
                    else:
                        # alpha^3 - (1 - 2 bR) alpha^2 - 4 bR alpha - bR FF / 2 = 0, largest root
                        cc = 0.5 * bR * FF
                        al = 1.0 + 2.0 * math.sqrt(bR) + cc ** (1.0 / 3.0)
                        for it in range(8):
                            f = al * al * al - (1.0 - 2.0 * bR) * al * al - 4.0 * bR * al - cc
                            fp = 3.0 * al * al - 2.0 * (1.0 - 2.0 * bR) * al - 4.0 * bR
                            al -= f / fp
                        v = (al - 1.0 + 2.0 * bR) / (2.0 * bR)
                    sig = math.sqrt(0.5 / al)
                    m0 = 0.5 * F0 / al; m1 = 0.5 * F1 / al; m2 = 0.5 * F2 / al; m3 = 0.5 * F3 / al
                    for t in range(ntry):
                        p0 = m0 + sig * np.random.standard_normal()
                        p1 = m1 + sig * np.random.standard_normal()
                        p2 = m2 + sig * np.random.standard_normal()
                        p3 = m3 + sig * np.random.standard_normal()
                        ok = True
                        if bR != 0.0:
                            r2 = p0 * p0 + p1 * p1 + p2 * p2 + p3 * p3
                            d = r2 - v
                            ok = np.random.random() < math.exp(-bR * d * d)
                        if ok:
                            phi[c0, c1, c2, 0] = p0; phi[c0, c1, c2, 1] = p1
                            phi[c0, c1, c2, 2] = p2; phi[c0, c1, c2, 3] = p3
                            nacc += 1
                            break
    return nacc


@njit(cache=True)
def sweep_higgs_or(U, phi, L, bH, bR):
    """Overrelaxation of hep-lat/9510020 section 3: reflection of the components orthogonal to the
    force and the root X' of V(X') = V(X) along the force, accepted with |V'(X) / V'(X')|.
    Returns the number of accepted X moves."""
    nacc = 0
    for par in range(2):
        for c0 in range(L):
            for c1 in range(L):
                for c2 in range((c0 + c1 + par) % 2, L, 2):
                    b0, b1, b2, b3 = higgs_force(U, phi, L, c0, c1, c2)
                    F = bH * math.sqrt(b0 * b0 + b1 * b1 + b2 * b2 + b3 * b3)
                    if F < 1e-14:
                        continue
                    f0 = bH * b0 / F; f1 = bH * b1 / F; f2 = bH * b2 / F; f3 = bH * b3 / F
                    p = phi[c0, c1, c2]
                    X = p[0] * f0 + p[1] * f1 + p[2] * f2 + p[3] * f3
                    y0 = p[0] - X * f0; y1 = p[1] - X * f1; y2 = p[2] - X * f2; y3 = p[3] - X * f3
                    Y2 = y0 * y0 + y1 * y1 + y2 * y2 + y3 * y3
                    A = 1.0 + 2.0 * bR * (Y2 - 1.0)
                    # V(X) = -X F + A X^2 + bR X^4 ; other real root of V(X') = V(X)
                    Xn = F / A - X
                    if bR != 0.0:
                        for it in range(30):
                            g = bR * (Xn * Xn * Xn + Xn * Xn * X + Xn * X * X + X * X * X) + A * (Xn + X) - F
                            gp = bR * (3.0 * Xn * Xn + 2.0 * Xn * X + X * X) + A
                            dx = g / gp
                            Xn -= dx
                            if abs(dx) < 1e-14 * (1.0 + abs(Xn)):
                                break
                    dV = -F + 2.0 * A * X + 4.0 * bR * X * X * X
                    dVn = -F + 2.0 * A * Xn + 4.0 * bR * Xn * Xn * Xn
                    acc = True
                    if dVn == 0.0:
                        acc = False
                    else:
                        ratio = abs(dV / dVn)
                        if ratio < 1.0:
                            acc = np.random.random() < ratio
                    if acc:
                        nacc += 1
                        Xu = Xn
                    else:
                        Xu = X
                    phi[c0, c1, c2, 0] = Xu * f0 - y0; phi[c0, c1, c2, 1] = Xu * f1 - y1
                    phi[c0, c1, c2, 2] = Xu * f2 - y2; phi[c0, c1, c2, 3] = Xu * f3 - y3
    return nacc


@njit(cache=True)
def measure_higgs(U, phi, L):
    """Returns (H, S2, S4): H = sum_x sum_i 1/2 Tr Phi^dag(x) U_i(x) Phi(x+i), S2 = sum R^2, S4 = sum R^4."""
    H = 0.0; S2 = 0.0; S4 = 0.0
    for c0 in range(L):
        for c1 in range(L):
            for c2 in range(L):
                p = phi[c0, c1, c2]
                r2 = p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + p[3] * p[3]
                S2 += r2; S4 += r2 * r2
                for mu in range(3):
                    f0, f1, f2 = nb(L, c0, c1, c2, mu, 1)
                    u = U[mu, c0, c1, c2]
                    q = phi[f0, f1, f2]
                    t0, t1, t2, t3 = qmul(u[0], u[1], u[2], u[3], q[0], q[1], q[2], q[3])
                    H += p[0] * t0 + p[1] * t1 + p[2] * t2 + p[3] * t3
    return H, S2, S4


@njit(cache=True)
def measure_plaq(U, L):
    """Mean of 1/2 Tr P over all plaquettes."""
    s = 0.0
    for c0 in range(L):
        for c1 in range(L):
            for c2 in range(L):
                for mu in range(3):
                    for nu in range(mu + 1, 3):
                        f0, f1, f2 = nb(L, c0, c1, c2, mu, 1)
                        g0, g1, g2 = nb(L, c0, c1, c2, nu, 1)
                        a = U[mu, c0, c1, c2]; b = U[nu, f0, f1, f2]
                        c = U[mu, g0, g1, g2]; d = U[nu, c0, c1, c2]
                        t0, t1, t2, t3 = qmul(a[0], a[1], a[2], a[3], b[0], b[1], b[2], b[3])
                        t0, t1, t2, t3 = qmul(t0, t1, t2, t3, c[0], -c[1], -c[2], -c[3])
                        t0, t1, t2, t3 = qmul(t0, t1, t2, t3, d[0], -d[1], -d[2], -d[3])
                        s += t0
    return s / (3.0 * L * L * L)


@njit(cache=True)
def global_radial(phi, L, bH, bR, H, S2, S4, eps, nhit):
    """Global update R -> exp(xi) R of hep-lat/9510020 section 3. Returns updated sums and acceptances."""
    V = L * L * L
    nacc = 0
    for h in range(nhit):
        xi = eps * (2.0 * np.random.random() - 1.0)
        a = -bH * H + (1.0 - 2.0 * bR) * S2
        b = bR * S4
        e2 = math.exp(2.0 * xi)
        dS = a * (e2 - 1.0) + b * (e2 * e2 - 1.0)
        la = 4.0 * V * xi - dS
        if la >= 0.0 or np.random.random() < math.exp(la):
            e1 = math.exp(xi)
            for c0 in range(L):
                for c1 in range(L):
                    for c2 in range(L):
                        for k in range(4):
                            phi[c0, c1, c2, k] *= e1
            H *= e2; S2 *= e2; S4 *= e2 * e2
            nacc += 1
    return H, S2, S4, nacc


@njit(cache=True)
def seed_numba(s):
    np.random.seed(s)


@njit(cache=True)
def compound_sweep(U, phi, L, bG, bH, bR, n_or, eps, update_links, update_higgs):
    """One heat-bath sweep, n_or overrelaxation sweeps, global radial hits. Returns (H, S2, S4, acc_or, acc_gl)."""
    acc_or = 0
    if update_links:
        sweep_links_hb(U, phi, L, bG, bH)
    if update_higgs:
        sweep_higgs_hb(U, phi, L, bH, bR, 20)
    for i in range(n_or):
        if update_links:
            sweep_links_or(U, phi, L, bG, bH)
        if update_higgs:
            acc_or += sweep_higgs_or(U, phi, L, bH, bR)
    H, S2, S4 = measure_higgs(U, phi, L)
    acc_gl = 0
    if update_higgs and eps > 0.0:
        H, S2, S4, acc_gl = global_radial(phi, L, bH, bR, H, S2, S4, eps, 4)
    return H, S2, S4, acc_or, acc_gl


@njit(cache=True)
def run_series(U, phi, L, bG, bH, bR, n_or, eps, nsweep, plaq_every, update_links, update_higgs):
    out = np.zeros((nsweep, 4))
    acc = np.zeros(2)
    for n in range(nsweep):
        H, S2, S4, a1, a2 = compound_sweep(U, phi, L, bG, bH, bR, n_or, eps, update_links, update_higgs)
        out[n, 0] = H; out[n, 1] = S2; out[n, 2] = S4
        if plaq_every > 0 and n % plaq_every == 0:
            out[n, 3] = measure_plaq(U, L)
        else:
            out[n, 3] = np.nan
        acc[0] += a1; acc[1] += a2
    return out, acc


def cold_start(n_side, radius=1.0):
    U = np.zeros((3, n_side, n_side, n_side, 4)); U[..., 0] = 1.0
    phi = np.zeros((n_side, n_side, n_side, 4)); phi[..., 0] = radius
    return U, phi


def hot_start(n_side, rng, radius=1.0):
    U = rng.standard_normal((3, n_side, n_side, n_side, 4))
    U /= np.linalg.norm(U, axis=-1, keepdims=True)
    phi = rng.standard_normal((n_side, n_side, n_side, 4)) * (radius / 2.0)
    return U, phi
