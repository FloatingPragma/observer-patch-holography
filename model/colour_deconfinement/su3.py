"""Pure SU(3) lattice gauge theory with the Wilson plaquette action, numpy implementation.

Action S = beta sum_P (1 - Re tr U_P / 3) on a periodic lattice with even extents dims = (N0, N1, N2, N3).
Axis 3 is the time direction. Links U[mu] have shape dims + (3, 3) and point from x to x + e_mu.
Update: Cabibbo-Marinari heat bath over the SU(2) subgroups (0,1), (0,2), (1,2) and subgroup overrelaxation.
Links of one direction on sites of one parity share no plaquette. They are updated together.
The numba kernel in su3_numba.py runs the same algorithm in sequential order and is used when numba is present.
"""
import numpy as np

OMEGA = np.exp(2j * np.pi / 3.0)
SUBGROUPS = ((0, 1), (0, 2), (1, 2))


def dag(a):
    return np.conj(np.swapaxes(a, -1, -2))


def random_su3(n, rng):
    """n Haar-random SU(3) matrices."""
    z = rng.normal(size=(n, 3, 3)) + 1j * rng.normal(size=(n, 3, 3))
    q, r = np.linalg.qr(z)
    d = np.diagonal(r, axis1=1, axis2=2)
    q = q * (d / np.abs(d))[:, None, :]
    return q / (np.linalg.det(q) ** (1.0 / 3.0))[:, None, None]


def draw_x0(a, rng):
    """x0 in [-1, 1] with density sqrt(1 - x0^2) exp(a x0): Kennedy-Pendleton above a = 1.5, Creutz below."""
    a = np.asarray(a, float)
    x0 = np.empty(len(a))
    todo = np.arange(len(a))
    while len(todo):
        aa = a[todo]
        n = len(todo)
        r1 = 1.0 - rng.random(n)
        r2 = rng.random(n)
        r3 = 1.0 - rng.random(n)
        r4 = rng.random(n)
        delta = -(np.log(r1) + np.cos(2.0 * np.pi * r2) ** 2 * np.log(r3)) / np.maximum(aa, 1e-6)
        x_kp = 1.0 - delta
        acc_kp = r4 * r4 <= 1.0 - 0.5 * delta
        s = np.maximum(aa, 1e-6)
        ea = np.exp(-2.0 * s)
        x_cr = 1.0 + np.log(ea + r2 * (1.0 - ea)) / s
        acc_cr = r4 * r4 <= 1.0 - x_cr * x_cr
        big = aa > 1.5
        x = np.where(big, x_kp, x_cr)
        acc = np.where(big, acc_kp, acc_cr)
        x0[todo[acc]] = x[acc]
        todo = todo[~acc]
    return x0


class Lattice:
    def __init__(self, dims, beta, seed=0, start='cold'):
        if len(dims) != 4 or any(n % 2 for n in dims):
            raise ValueError('four even extents are required')
        self.dims = tuple(int(n) for n in dims)
        self.beta = float(beta)
        self.rng = np.random.default_rng(seed)
        self.volume = int(np.prod(self.dims))
        self.U = np.zeros((4,) + self.dims + (3, 3), complex)
        if start == 'cold':
            self.U[...] = np.eye(3)
        else:
            self.U[...] = random_su3(4 * self.volume, self.rng).reshape((4,) + self.dims + (3, 3))
        par = np.indices(self.dims).sum(0) % 2
        self.parity = (par == 0, par == 1)

    # ------------------------------------------------------------------------------------------------ update
    def staples(self, mu):
        """A(x) with local weight exp((beta/3) Re tr(U_mu(x) A(x)))."""
        U = self.U
        A = np.zeros_like(U[mu])
        for nu in range(4):
            if nu == mu:
                continue
            up = np.roll(U[nu], -1, axis=mu)
            A += up @ dag(np.roll(U[mu], -1, axis=nu)) @ dag(U[nu])
            A += np.roll(dag(up) @ dag(U[mu]) @ U[nu], 1, axis=nu)
        return A

    def _update(self, mu, par, heat):
        mask = self.parity[par]
        A = self.staples(mu)[mask]
        U = self.U[mu][mask]
        W = U @ A
        n = len(U)
        for i, j in SUBGROUPS:
            v0 = 0.5 * (W[:, i, i].real + W[:, j, j].real)
            v3 = 0.5 * (W[:, i, i].imag - W[:, j, j].imag)
            v1 = 0.5 * (W[:, i, j].imag + W[:, j, i].imag)
            v2 = 0.5 * (W[:, i, j].real - W[:, j, i].real)
            k = np.maximum(np.sqrt(v0 * v0 + v1 * v1 + v2 * v2 + v3 * v3), 1e-300)
            v0, v1, v2, v3 = v0 / k, v1 / k, v2 / k, v3 / k
            if heat:
                x0 = draw_x0(2.0 * self.beta * k / 3.0, self.rng)
                r = np.sqrt(np.maximum(0.0, 1.0 - x0 * x0))
                ct = self.rng.uniform(-1.0, 1.0, n)
                st = np.sqrt(np.maximum(0.0, 1.0 - ct * ct))
                ph = self.rng.uniform(0.0, 2.0 * np.pi, n)
                x1, x2, x3 = r * st * np.cos(ph), r * st * np.sin(ph), r * ct
                a0 = x0 * v0 + x1 * v1 + x2 * v2 + x3 * v3
                a1 = -x0 * v1 + v0 * x1 + (x2 * v3 - x3 * v2)
                a2 = -x0 * v2 + v0 * x2 + (x3 * v1 - x1 * v3)
                a3 = -x0 * v3 + v0 * x3 + (x1 * v2 - x2 * v1)
            else:
                a0 = v0 * v0 - v1 * v1 - v2 * v2 - v3 * v3
                a1, a2, a3 = -2.0 * v0 * v1, -2.0 * v0 * v2, -2.0 * v0 * v3
            c11, c12, c21, c22 = a0 + 1j * a3, a2 + 1j * a1, -a2 + 1j * a1, a0 - 1j * a3
            for M in (U, W):
                ri = M[:, i, :].copy()
                rj = M[:, j, :].copy()
                M[:, i, :] = c11[:, None] * ri + c12[:, None] * rj
                M[:, j, :] = c21[:, None] * ri + c22[:, None] * rj
        self.U[mu][mask] = U

    def sweep(self, heat=True):
        for mu in range(4):
            for par in (0, 1):
                self._update(mu, par, heat)

    def compound_sweeps(self, n, n_or=4):
        for _ in range(n):
            self.sweep(True)
            for _ in range(n_or):
                self.sweep(False)

    def reunitarize(self):
        M = self.U.reshape(-1, 3, 3)
        r0 = M[:, 0] / np.linalg.norm(M[:, 0], axis=1)[:, None]
        r1 = M[:, 1] - np.sum(np.conj(r0) * M[:, 1], axis=1)[:, None] * r0
        r1 = r1 / np.linalg.norm(r1, axis=1)[:, None]
        r2 = np.conj(np.cross(r0, r1))
        self.U = np.stack([r0, r1, r2], axis=1).reshape(self.U.shape)

    # ------------------------------------------------------------------------------------------ observables
    def plaquette_matrix(self, mu, nu):
        U = self.U
        return U[mu] @ np.roll(U[nu], -1, axis=mu) @ dag(np.roll(U[mu], -1, axis=nu)) @ dag(U[nu])

    def plaquette(self):
        """Average of Re tr U_P / 3 over sites and the six planes."""
        tot = 0.0
        for mu in range(3):
            for nu in range(mu + 1, 4):
                tot += np.trace(self.plaquette_matrix(mu, nu), axis1=-2, axis2=-1).real.mean()
        return tot / 18.0

    def polyakov_field(self):
        """tr prod_t U_3(x, t) / 3 on every spatial site."""
        P = self.U[3][:, :, :, 0]
        for t in range(1, self.dims[3]):
            P = P @ self.U[3][:, :, :, t]
        return np.trace(P, axis1=-2, axis2=-1) / 3.0

    def polyakov(self):
        return complex(self.polyakov_field().mean())

    def unitarity_defect(self):
        M = self.U.reshape(-1, 3, 3)
        return float(np.abs(M @ dag(M) - np.eye(3)).max()), float(np.abs(np.linalg.det(M) - 1.0).max())

    # -------------------------------------------------------------------------------------- transformations
    def gauge_transform(self, g):
        """U_mu(x) -> g(x) U_mu(x) g(x + mu)^dag, g of shape dims + (3, 3)."""
        for mu in range(4):
            self.U[mu] = g @ self.U[mu] @ dag(np.roll(g, -1, axis=mu))

    def centre_transform(self, t_slice=0, k=1, region=None):
        """Multiply the time-like links leaving the time slice by omega^k, on all spatial sites or on a region."""
        if region is None:
            self.U[3][:, :, :, t_slice] *= OMEGA ** k
        else:
            self.U[3][:, :, :, t_slice][region] *= OMEGA ** k

    def centre_flux(self, mu, nu):
        """Class of each plaquette holonomy in Z3 if it is a centre element: n with U_P = omega^n, else -1."""
        P = self.plaquette_matrix(mu, nu)
        out = np.full(self.dims, -1)
        for n in range(3):
            out[np.abs(P - OMEGA ** n * np.eye(3)).max(axis=(-1, -2)) < 1e-12] = n
        return out


def one_plaquette_average(beta, n=400):
    """<Re tr U / 3> with weight exp((beta/3) Re tr U) over SU(3): Weyl integration on the maximal torus."""
    t = 2.0 * np.pi * (np.arange(n) + 0.5) / n
    a, b = np.meshgrid(t, t, indexing='ij')
    c = -a - b
    e = (np.exp(1j * a), np.exp(1j * b), np.exp(1j * c))
    haar = np.abs(e[0] - e[1]) ** 2 * np.abs(e[0] - e[2]) ** 2 * np.abs(e[1] - e[2]) ** 2
    retr = np.cos(a) + np.cos(b) + np.cos(c)
    w = haar * np.exp(beta / 3.0 * retr)
    return float((w * retr / 3.0).sum() / w.sum())


def haar_moment_re_trace(power, n=240):
    """<(Re tr U)^power> over SU(3) by Weyl integration. Exact for trigonometric polynomials below degree n."""
    t = 2.0 * np.pi * np.arange(n) / n
    a, b = np.meshgrid(t, t, indexing='ij')
    c = -a - b
    e = (np.exp(1j * a), np.exp(1j * b), np.exp(1j * c))
    haar = np.abs(e[0] - e[1]) ** 2 * np.abs(e[0] - e[2]) ** 2 * np.abs(e[1] - e[2]) ** 2
    retr = np.cos(a) + np.cos(b) + np.cos(c)
    return float((haar * retr ** power).sum() / haar.sum())


def series(dims, beta, sweeps, n_or=4, seed=0, start='cold', kernel='auto'):
    """Plaquette and complex Polyakov loop after each compound sweep.

    kernel 'auto' uses the numba kernel when numba is installed and the numpy kernel otherwise.
    Returns (plaquette series, loop series, name of the kernel used).
    """
    K = None
    if kernel in ('auto', 'numba'):
        from . import su3_numba as K
        if not K.AVAILABLE:
            if kernel == 'numba':
                raise ImportError('the compiled kernel needs numba')
            K = None
    P = np.empty(sweeps)
    L = np.empty(sweeps, complex)
    if K is not None:
        up, dn = K.neighbours(*dims)
        K.seed_rng(seed)
        volume = int(np.prod(dims))
        U = K.cold_start(volume) if start == 'cold' else K.hot_start(volume)
        for n in range(sweeps):
            K.compound_sweeps(U, up, dn, beta, 1, n_or)
            if (n + 1) % 50 == 0:
                K.reunitarize(U)
            a, b = K.plaquette(U, up)
            P[n] = 0.5 * (a + b)
            L[n] = K.polyakov(U, up, dims[3])
        return P, L, 'numba'
    lat = Lattice(dims, beta, seed=seed, start=start)
    for n in range(sweeps):
        lat.compound_sweeps(1, n_or)
        if (n + 1) % 50 == 0:
            lat.reunitarize()
        P[n] = lat.plaquette()
        L[n] = lat.polyakov()
    return P, L, 'numpy'
