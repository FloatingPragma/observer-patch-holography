"""Three-dimensional SU(2) gauge theory with a fundamental Higgs doublet on a cubic lattice.

Self-contained numpy implementation with checkerboard updates. Action (hep-lat/9510020, section 2):

    S = beta_G sum_x sum_{i<j} [1 - 1/2 Tr P_ij(x)]
        - beta_H sum_x sum_i 1/2 Tr Phi^dag(x) U_i(x) Phi(x+i)
        + sum_x R^2(x) + beta_R sum_x [R^2(x) - 1]^2 ,      R^2 = 1/2 Tr Phi^dag Phi .

Links and the Higgs matrix are stored as four real numbers c = (c0, c1, c2, c3) for the matrix
[[p, q], [-conj(q), conj(p)]] with p = c0 + i c1, q = c2 + i c3; then 1/2 Tr A^dag B = dot(a, b).
"""
from __future__ import annotations

import math

import numpy as np

SIGMA = 3.175911535625          # hep-lat/9705003, appendix D
ZETA = 0.08849
DELTA = 1.942130
TWO_LOOP_G4 = 4.9941            # hep-lat/9705003, section 5.1
TWO_LOOP_LG2 = 5.2153
CONDENSATE_CONSTANT = ZETA + SIGMA * SIGMA / 4.0 - DELTA      # 0.66796


# --------------------------------------------------------------------------------------------
# relations between lattice and continuum parameters
# --------------------------------------------------------------------------------------------
def y_counterterm(x: float, beta_G: float) -> float:
    """Linear and logarithmic counterterms of the mass parameter in units of g3^4."""
    return (3.0 * SIGMA * beta_G / (32.0 * math.pi) * (1.0 + 4.0 * x)
            + ((51.0 / 16.0 + 9.0 * x - 12.0 * x * x) * (math.log(1.5 * beta_G) + ZETA)
               + TWO_LOOP_G4 + TWO_LOOP_LG2 * x) / (16.0 * math.pi ** 2))


def continuum_parameters(beta_G: float, beta_H: float, beta_R: float) -> tuple[float, float]:
    """(x, y) of the continuum theory for given lattice couplings."""
    x = beta_R * beta_G / (beta_H * beta_H)
    y = beta_G * beta_G / 8.0 * (1.0 / beta_H - 3.0 - 2.0 * x * beta_H / beta_G) + y_counterterm(x, beta_G)
    return x, y


def lattice_couplings(x: float, y: float, beta_G: float) -> tuple[float, float]:
    """(beta_H, beta_R) for given continuum parameters at lattice spacing a = 4 / (beta_G g3^2)."""
    bare = y - y_counterterm(x, beta_G)
    a = 3.0 + 8.0 * bare / (beta_G * beta_G)
    if x == 0.0:
        beta_H = 1.0 / a
    else:
        beta_H = (-a + math.sqrt(a * a + 8.0 * x / beta_G)) / (4.0 * x / beta_G)
    return beta_H, x * beta_H * beta_H / beta_G


def condensate(mean_R2: float, beta_G: float, beta_H: float) -> float:
    """<phi^dag phi> / g3^2 in the MS-bar scheme at the scale g3^2."""
    return (beta_G * beta_H / 8.0 * (mean_R2 - SIGMA / (math.pi * beta_H))
            - 3.0 / (16.0 * math.pi ** 2) * (math.log(1.5 * beta_G) + CONDENSATE_CONSTANT))


def susceptibility(var_R2bar: float, n_side: int, beta_G: float, beta_H: float) -> float:
    """chi_3 = g3^2 V <(phi^dag phi / g3^2 fluctuation of the volume average)^2>."""
    return n_side ** 3 * beta_H * beta_H / beta_G * var_R2bar


# --------------------------------------------------------------------------------------------
# quaternions
# --------------------------------------------------------------------------------------------
def qmul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Matrix product in the four-real storage: (p1, q1)(p2, q2) = (p1 p2 - q1 conj(q2), p1 q2 + q1 conj(p2))."""
    ac = np.ascontiguousarray(a, dtype=np.float64).view(np.complex128)
    bc = np.ascontiguousarray(b, dtype=np.float64).view(np.complex128)
    out = np.empty(np.broadcast_shapes(ac.shape, bc.shape), dtype=np.complex128)
    out[..., 0] = ac[..., 0] * bc[..., 0] - ac[..., 1] * np.conj(bc[..., 1])
    out[..., 1] = ac[..., 0] * bc[..., 1] + ac[..., 1] * np.conj(bc[..., 0])
    return out.view(np.float64)


def qdag(a: np.ndarray) -> np.ndarray:
    out = a.copy()
    out[..., 1:] *= -1.0
    return out


def random_su2(shape: tuple[int, ...], rng: np.random.Generator) -> np.ndarray:
    g = rng.standard_normal(shape + (4,))
    return g / np.linalg.norm(g, axis=-1, keepdims=True)


class Lattice:
    """State and updates of the lattice theory on N^3 sites with periodic boundaries (N even)."""

    def __init__(self, n_side: int, beta_G: float, beta_H: float, beta_R: float,
                 rng: np.random.Generator, start: str = "cold", radius: float = 1.0):
        if n_side % 2:
            raise ValueError("checkerboard updates need an even side length")
        self.N = n_side
        self.V = n_side ** 3
        self.beta_G, self.beta_H, self.beta_R = float(beta_G), float(beta_H), float(beta_R)
        self.rng = rng
        shape = (n_side,) * 3
        if start == "cold":
            self.U = np.zeros((3,) + shape + (4,))
            self.U[..., 0] = 1.0
            self.phi = np.zeros(shape + (4,))
            self.phi[..., 0] = radius
        else:
            self.U = random_su2((3,) + shape, rng)
            self.phi = rng.standard_normal(shape + (4,)) * (radius / 2.0)
        idx = np.indices(shape).sum(axis=0)
        self.parity = [(idx % 2) == 0, (idx % 2) == 1]
        self.update_links = True
        self.accept_overrelax = [0, 0]
        self.accept_global = [0, 0]

    # ---- forces ---------------------------------------------------------------------------
    def staple(self, mu: int) -> np.ndarray:
        """Sum of the four staples of every link in direction mu."""
        U = self.U
        total = np.zeros_like(U[mu])
        for nu in range(3):
            if nu == mu:
                continue
            u_nu_fwd = np.roll(U[nu], -1, axis=mu)                     # U_nu(x + mu)
            upper = qmul(qmul(u_nu_fwd, qdag(np.roll(U[mu], -1, axis=nu))), qdag(U[nu]))
            lower = qmul(qmul(qdag(u_nu_fwd), qdag(U[mu])), U[nu])     # evaluated at x - nu below
            total += upper + np.roll(lower, 1, axis=nu)
        return total

    def link_force(self, mu: int) -> np.ndarray:
        """W with local link weight exp(1/2 Tr[U W])."""
        w = self.beta_G * self.staple(mu)
        if self.beta_H != 0.0:
            w += self.beta_H * qmul(np.roll(self.phi, -1, axis=mu), qdag(self.phi))
        return w

    def higgs_force(self) -> np.ndarray:
        """b(x) = sum_i [U_i(x) Phi(x+i) + U_i^dag(x-i) Phi(x-i)]."""
        b = np.zeros_like(self.phi)
        for mu in range(3):
            b += qmul(self.U[mu], np.roll(self.phi, -1, axis=mu))
            b += np.roll(qmul(qdag(self.U[mu]), self.phi), 1, axis=mu)
        return b

    # ---- link updates ---------------------------------------------------------------------
    def _sample_x0(self, k: np.ndarray) -> np.ndarray:
        """x0 with density exp(k x0) sqrt(1 - x0^2) on [-1, 1] (Kennedy-Pendleton; Creutz for small k)."""
        rng = self.rng
        x0 = np.empty_like(k)
        todo = np.ones(k.shape, bool)
        while todo.any():
            kk = k[todo]
            n = kk.size
            big = kk > 1.6
            trial = np.empty(n)
            ok = np.zeros(n, bool)
            if big.any():
                kb = kk[big]
                r1 = 1.0 - rng.random(kb.size)
                r3 = 1.0 - rng.random(kb.size)
                c = np.cos(2.0 * np.pi * rng.random(kb.size))
                lam2 = -(np.log(r1) + c * c * np.log(r3)) / (2.0 * kb)
                trial[big] = 1.0 - 2.0 * lam2
                ok[big] = rng.random(kb.size) ** 2 <= 1.0 - lam2
            if (~big).any():
                ks = kk[~big]
                em = np.exp(-2.0 * ks)
                t = 1.0 + np.log(em + (1.0 - em) * rng.random(ks.size)) / ks
                trial[~big] = t
                ok[~big] = rng.random(ks.size) ** 2 <= 1.0 - t * t
            where = np.flatnonzero(todo)
            x0[where[ok]] = trial[ok]
            todo[where[ok]] = False
        return x0

    def heatbath_links(self) -> None:
        rng = self.rng
        for mu in range(3):
            for par in self.parity:
                w = self.link_force(mu)[par]
                k = np.linalg.norm(w, axis=-1)
                k = np.maximum(k, 1e-300)
                x0 = self._sample_x0(k)
                rad = np.sqrt(np.maximum(0.0, 1.0 - x0 * x0))
                ct = 2.0 * rng.random(k.size) - 1.0
                st = np.sqrt(np.maximum(0.0, 1.0 - ct * ct))
                ph = 2.0 * np.pi * rng.random(k.size)
                x = np.stack([x0, rad * st * np.cos(ph), rad * st * np.sin(ph), rad * ct], axis=-1)
                new = qmul(x, qdag(w / k[:, None]))
                new /= np.linalg.norm(new, axis=-1, keepdims=True)
                self.U[mu][par] = new

    def overrelax_links(self) -> None:
        for mu in range(3):
            for par in self.parity:
                w = self.link_force(mu)[par]
                k = np.maximum(np.linalg.norm(w, axis=-1, keepdims=True), 1e-300)
                v0d = qdag(w / k)
                new = qmul(qmul(v0d, qdag(self.U[mu][par])), v0d)
                new /= np.linalg.norm(new, axis=-1, keepdims=True)
                self.U[mu][par] = new

    # ---- Higgs updates --------------------------------------------------------------------
    def heatbath_higgs(self, max_tries: int = 20) -> None:
        """Gaussian proposal around F/(2 alpha), accepted with exp(-beta_R (R^2 - v)^2)."""
        rng = self.rng
        bR = self.beta_R
        for par in self.parity:
            F = self.beta_H * self.higgs_force()[par]
            FF = (F * F).sum(-1)
            if bR == 0.0:
                alpha = np.ones_like(FF)
                v = np.zeros_like(FF)
            else:
                c = 0.5 * bR * FF
                alpha = 1.0 + 2.0 * math.sqrt(bR) + np.cbrt(c)
                for _ in range(8):
                    f = alpha ** 3 - (1.0 - 2.0 * bR) * alpha ** 2 - 4.0 * bR * alpha - c
                    fp = 3.0 * alpha ** 2 - 2.0 * (1.0 - 2.0 * bR) * alpha - 4.0 * bR
                    alpha = alpha - f / fp
                v = (alpha - 1.0 + 2.0 * bR) / (2.0 * bR)
            centre = 0.5 * F / alpha[:, None]
            sigma = np.sqrt(0.5 / alpha)[:, None]
            new = self.phi[par].copy()
            todo = np.ones(FF.size, bool)
            for _ in range(max_tries):
                idx = np.flatnonzero(todo)
                if idx.size == 0:
                    break
                trial = centre[idx] + sigma[idx] * rng.standard_normal((idx.size, 4))
                if bR == 0.0:
                    ok = np.ones(idx.size, bool)
                else:
                    r2 = (trial * trial).sum(-1)
                    ok = rng.random(idx.size) < np.exp(-bR * (r2 - v[idx]) ** 2)
                new[idx[ok]] = trial[ok]
                todo[idx[ok]] = False
            self.phi[par] = new

    def overrelax_higgs(self) -> None:
        """Reflection of the components orthogonal to the force; along the force the other root of
        V(X') = V(X), accepted with |V'(X) / V'(X')| (hep-lat/9510020, section 3)."""
        rng = self.rng
        bR = self.beta_R
        for par in self.parity:
            b = self.beta_H * self.higgs_force()[par]
            F = np.linalg.norm(b, axis=-1)
            good = F > 1e-14
            f = np.zeros_like(b)
            f[good] = b[good] / F[good, None]
            p = self.phi[par]
            X = (p * f).sum(-1)
            Y = p - X[:, None] * f
            Y2 = (Y * Y).sum(-1)
            A = 1.0 + 2.0 * bR * (Y2 - 1.0)
            Xn = F / A - X
            if bR != 0.0:
                for _ in range(30):
                    g = bR * (Xn ** 3 + Xn * Xn * X + Xn * X * X + X ** 3) + A * (Xn + X) - F
                    gp = bR * (3.0 * Xn * Xn + 2.0 * Xn * X + X * X) + A
                    Xn = Xn - g / gp
            dV = -F + 2.0 * A * X + 4.0 * bR * X ** 3
            dVn = -F + 2.0 * A * Xn + 4.0 * bR * Xn ** 3
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = np.where(dVn != 0.0, np.abs(dV / dVn), 0.0)
            acc = (rng.random(X.size) < ratio) & good
            self.accept_overrelax[0] += int(acc.sum())
            self.accept_overrelax[1] += int(good.sum())
            Xu = np.where(acc, Xn, X)
            new = Xu[:, None] * f - Y
            new[~good] = p[~good]
            self.phi[par] = new

    def global_radial(self, width: float, hits: int = 4) -> None:
        """R -> exp(xi) R on all sites, accepted with exp(4 V xi - Delta S)."""
        H, S2, S4 = self.higgs_sums()
        for _ in range(hits):
            xi = width * (2.0 * self.rng.random() - 1.0)
            a = -self.beta_H * H + (1.0 - 2.0 * self.beta_R) * S2
            b = self.beta_R * S4
            e2 = math.exp(2.0 * xi)
            log_acc = 4.0 * self.V * xi - a * (e2 - 1.0) - b * (e2 * e2 - 1.0)
            self.accept_global[1] += 1
            if log_acc >= 0.0 or self.rng.random() < math.exp(log_acc):
                self.phi *= math.exp(xi)
                H *= e2
                S2 *= e2
                S4 *= e2 * e2
                self.accept_global[0] += 1

    def sweep(self, n_overrelax: int = 4, global_width: float | None = None) -> None:
        """One compound sweep: heat bath, n_overrelax overrelaxation sweeps, global radial hits."""
        if self.update_links:
            self.heatbath_links()
        self.heatbath_higgs()
        for _ in range(n_overrelax):
            if self.update_links:
                self.overrelax_links()
            self.overrelax_higgs()
        if global_width is None:
            global_width = 0.5 / math.sqrt(self.V)
        if global_width > 0.0:
            self.global_radial(global_width)

    # ---- observables ----------------------------------------------------------------------
    def higgs_sums(self) -> tuple[float, float, float]:
        """(H, S2, S4): hopping sum, sum of R^2, sum of R^4."""
        H = 0.0
        for mu in range(3):
            H += float((self.phi * qmul(self.U[mu], np.roll(self.phi, -1, axis=mu))).sum())
        r2 = (self.phi * self.phi).sum(-1)
        return H, float(r2.sum()), float((r2 * r2).sum())

    def plaquette(self) -> float:
        """Mean of 1/2 Tr P over all plaquettes."""
        U = self.U
        total = 0.0
        for mu in range(3):
            for nu in range(mu + 1, 3):
                p = qmul(qmul(U[mu], np.roll(U[nu], -1, axis=mu)),
                         qmul(qdag(np.roll(U[mu], -1, axis=nu)), qdag(U[nu])))
                total += float(p[..., 0].sum())
        return total / (3.0 * self.V)

    def action(self) -> float:
        H, S2, S4 = self.higgs_sums()
        return (self.beta_G * 3.0 * self.V * (1.0 - self.plaquette()) - self.beta_H * H + S2
                + self.beta_R * (S4 - 2.0 * S2 + self.V))

    def gauge_transform(self, G: np.ndarray) -> None:
        """Phi(x) -> G(x) Phi(x), U_i(x) -> G(x) U_i(x) G^dag(x+i)."""
        self.phi = qmul(G, self.phi)
        for mu in range(3):
            self.U[mu] = qmul(qmul(G, self.U[mu]), qdag(np.roll(G, -1, axis=mu)))

    def run(self, n_sweeps: int, n_overrelax: int = 4, plaquette: bool = False) -> np.ndarray:
        """Series of (H, S2, S4[, plaquette]) over n_sweeps compound sweeps."""
        rows = []
        for _ in range(n_sweeps):
            self.sweep(n_overrelax)
            row = list(self.higgs_sums())
            if plaquette:
                row.append(self.plaquette())
            rows.append(row)
        return np.array(rows)


def free_field_R2(n_side: int, beta_H: float) -> tuple[float, float]:
    """Exact <R^2> and <sum_i phi(x).phi(x+i)> per site for unit links and beta_R = 0."""
    k = 2.0 * np.pi * np.arange(n_side) / n_side
    c = np.cos(k)[:, None, None] + np.cos(k)[None, :, None] + np.cos(k)[None, None, :]
    return float(2.0 * np.mean(1.0 / (1.0 - beta_H * c))), float(2.0 * np.mean(c / (1.0 - beta_H * c)))


def strong_coupling_plaquette(beta_G: float) -> float:
    """u + 2 u^5 with u = I_2(beta)/I_1(beta): pure SU(2) gauge theory in three dimensions, error O(u^9)."""
    def bessel_i(n: int, z: float) -> float:
        term = (z / 2.0) ** n / math.factorial(n)
        total = term
        for m in range(1, 60):
            term *= (z / 2.0) ** 2 / (m * (m + n))
            total += term
        return total
    u = bessel_i(2, beta_G) / bessel_i(1, beta_G)
    return u + 2.0 * u ** 5
