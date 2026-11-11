"""Exact statements about the class balance of one seam link.

Conventions: weak factor SU(2), coupling g^2, alpha = g^2/(4 pi), S = 2 pi/(4 alpha) = 2 pi^2/g^2, unit heat time
tau = g^2/2. A sector law is p_j proportional to d_j^w exp(-t C_j), n = d_j = 2j + 1, C_j = (n^2 - 1)/4.
The class balance is <A> = sum_j p_j (-1)^(2j). Every result is written <A> = M exp(-c S).
"""
import math
from fractions import Fraction as F

import mpmath as mp
import numpy as np
from scipy import optimize, sparse
from scipy.sparse import linalg as splinalg

from . import complexes as C

PI = math.pi
X_U = 24.3166                      # 1/alpha_U of the chain
G2_U = 4 * PI / X_U
S_U = 2 * PI ** 2 / G2_U


# ----------------------------------------------------------------------------------------------- one register
def register_law(w, t, dps=40):
    """<A> for weights d^w exp(-t C) on the sectors of SU(2)"""
    with mp.workdps(dps):
        num = den = mp.mpf(0)
        for n in range(1, 20000):
            e = mp.mpf(n) ** w * mp.exp(-mp.mpf(t) * (n * n - 1) / 4)
            num += (-1) ** (n - 1) * e
            den += e
            if n > 10 and e < mp.mpf(10) ** (-dps - 10) * den:
                break
        return num / den


def cartan_law(t, dps=40):
    """<A> on the Cartan circle: charges m/2, weights exp(-t m^2/4), eigenvalue (-1)^m"""
    with mp.workdps(dps):
        num = den = mp.mpf(0)
        for m in range(0, 20000):
            e = mp.exp(-mp.mpf(t) * m * m / 4) * (1 if m == 0 else 2)
            num += (-1) ** m * e
            den += e
            if m > 10 and e < mp.mpf(10) ** (-dps - 10) * den:
                break
        return num / den


def heat_kernel_character_sum(theta, g2, dps=40):
    """K_tau(theta) = sum_j d_j chi_j(theta) exp(-tau C_j), tau = g2/2"""
    with mp.workdps(dps):
        tau = mp.mpf(g2) / 2
        th = mp.mpf(theta)
        tot = mp.mpf(0)
        for n in range(1, 4000):
            ch = mp.mpf(n) if th == 0 else mp.sin(n * th / 2) / mp.sin(th / 2)
            e = mp.exp(-tau * (n * n - 1) / 4)
            tot += n * ch * e
            if n > 10 and n * n * e < mp.mpf(10) ** (-dps - 10):
                break
        return +tot


def ln_heat_kernel(theta, g2):
    """dual (Poisson resummed) form with the images k = -1, 0, 1, in the normalisation of the character sum.
    Relative truncation error below exp(-16 pi^2/g2)."""
    tau = g2 / 2
    th = theta
    de = 2 * PI - th
    norm = 0.25 * tau + 0.5 * math.log(4 * PI / tau) - math.log(tau)
    if th < 1e-6:
        return -th * th / (2 * g2) + norm
    if de < 1e-7:
        return -th * th / (2 * g2) + math.log(2 * (4 * PI * PI / g2 - 1)) + norm
    u = 4 * PI * de / g2
    br = -de - (PI + 0.5 * de) * math.expm1(-u) + (2 * PI + 0.5 * th) * math.exp(-4 * PI * (2 * PI + th) / g2)
    s = math.sin(0.5 * de) if de < 1.0 else math.sin(0.5 * th)
    return -th * th / (2 * g2) + math.log(br / s) + norm


def migdal_sphere(ratio_of_n, faces):
    """Z[one twisted face]/Z on a closed sphere of `faces` plaquettes: sum_n n^2 (-1)^(n-1) r_n^faces / sum_n n^2 r_n^faces,
    r_n = f_n/n the plaquette coefficient per dimension (heat kernel: exp(-tau C); Wilson: 2 I_n(beta)/beta)"""
    num = den = mp.mpf(0)
    for n in range(1, 3000):
        e = mp.mpf(n) ** 2 * ratio_of_n(n) ** faces
        num += (-1) ** (n - 1) * e
        den += e
        if n > 10 and e < mp.mpf(10) ** (-45) * den:
            break
    return num / den


def migdal_hk(g2, faces=20):
    tau = mp.mpf(g2) / 2
    return migdal_sphere(lambda n: mp.exp(-tau * (n * n - 1) / 4), faces)


def migdal_wilson(beta, faces=20):
    i1 = mp.besseli(1, beta)
    return migdal_sphere(lambda n: mp.besseli(n, beta) / i1, faces)


# ----------------------------------------------------------------------------------------------- Cartan circle on a complex
def sp_cob(cx):
    rows, cols, vals = [], [], []
    for p in range(cx["nplaq"]):
        for k in range(cx["pn"][p]):
            rows.append(p); cols.append(int(cx["pl"][p, k])); vals.append(float(cx["ps"][p, k]))
    return sparse.csr_matrix((vals, (rows, cols)), shape=(cx["nplaq"], cx["nlink"]))


def project(D, v, cols=None):
    """orthogonal projection of v on the range of the coboundary D (restricted to the free links cols)"""
    Dm = sparse.csr_matrix(D) if cols is None else sparse.csr_matrix(D)[:, cols]
    sol = splinalg.lsqr(Dm, v, atol=1e-14, btol=1e-14, iter_lim=50000)[0]
    return Dm @ sol


def one_minus_pi(cx, P0=None, cols=None):
    """c = 1 - Pi_PP of the plaquette P0: its Gaussian flux variance is Pi_PP g^2"""
    P0 = cx["P0"] if P0 is None else P0
    e = np.zeros(cx["nplaq"]); e[P0] = 1.0
    return 1.0 - float(project(sp_cob(cx), e, cols)[P0])


def exact_one_minus_pi(D, P0=0):
    """the same in exact rational arithmetic for a small integer coboundary matrix D (list of rows)"""
    rows = [[F(x) for x in r] for r in D]
    nP, nL = len(rows), len(rows[0])
    # column basis by elimination on D^T D is avoided: solve min |D x - e| with a pivoted basis of columns
    basis = []
    red = []                                   # reduced copies of the chosen columns
    for j in range(nL):
        col = [rows[i][j] for i in range(nP)]
        v = col[:]
        for (piv, b) in red:
            if v[piv] != 0:
                f = v[piv] / b[piv]
                v = [a - f * c for a, c in zip(v, b)]
        nz = next((i for i, a in enumerate(v) if a != 0), None)
        if nz is not None:
            red.append((nz, v)); basis.append(col)
    r = len(basis)
    G = [[sum(a * b for a, b in zip(basis[i], basis[j])) for j in range(r)] for i in range(r)]
    rhs = [basis[i][P0] for i in range(r)]
    # solve G y = rhs
    M = [G[i][:] + [rhs[i]] for i in range(r)]
    for i in range(r):
        p = next(k for k in range(i, r) if M[k][i] != 0)
        M[i], M[p] = M[p], M[i]
        inv = 1 / M[i][i]
        M[i] = [a * inv for a in M[i]]
        for k in range(r):
            if k != i and M[k][i] != 0:
                f = M[k][i]
                M[k] = [a - f * b for a, b in zip(M[k], M[i])]
    y = [M[i][r] for i in range(r)]
    pi_pp = sum(y[i] * basis[i][P0] for i in range(r))
    return 1 - pi_pp


def theta_ratio(c, g2):
    """exact class sum when the classes are multiples of one plaquette (closed surfaces, a pinned face pair)"""
    num = den = mp.mpf(0)
    for N in range(-80, 81):
        num += mp.exp(-8 * mp.pi ** 2 * c * (N - mp.mpf(1) / 2) ** 2 / g2)
        den += mp.exp(-8 * mp.pi ** 2 * c * N ** 2 / g2)
    return num / den


def villain_weak(c, g2):
    """leading form on the Cartan circle: the two orientations of the half turn"""
    return 2.0 * math.exp(-c * 2 * PI ** 2 / g2)


def icosahedral_space_many_ticks():
    """c of the seam's time-like plaquette for the icosahedral surface as space and many ticks:
    (1/30) sum over the nonzero eigenvalues lam of the face Laplacian of sqrt(lam/(lam + 4))"""
    D = C.coboundary(C.ico2())
    lam = np.linalg.eigvalsh(D @ D.T)
    lam = lam[lam > 1e-9]
    return float(np.sqrt(lam / (lam + 4)).sum() / 30), float(np.sqrt(lam).sum() / 60), lam


def cubic_c(d, L=None):
    return 1 - F(2, d) if L is None else 1 - F(2, d) * (1 - F(1, L ** d))


def face_step(cx):
    """c for A read at the midpoint of the face step: time-like plaquette above the seam fully twisted, the faces of
    the seam in that time slice half twisted. Gaussian exponent in units of S."""
    D = sp_cob(cx)
    P0 = cx["P0"]
    l = int(cx["pl"][P0, 3]) if cx["name"].startswith("cub") else int(cx["pl"][P0, 0])
    u = np.zeros(cx["nplaq"])
    nface = 0
    for a in range(cx["lp_ptr"][l], cx["lp_ptr"][l + 1]):
        p = int(cx["lp_idx"][a]); s = float(cx["ps"][p, int(cx["lp_pos"][a])])
        if p == P0:
            u[p] = 2 * PI * s
        elif cx["ptype"][p] == 0:
            u[p] = PI * s; nface += 1
    quad = float(u @ (u - project(D, u.copy()))) / (4 * PI ** 2)
    return quad + nface * 0.25, nface, quad


def wilson_minimum(cx, starts=3, seed=5):
    """g^2 times the minimum of the Wilson action with one twisted plaquette, links on the Cartan circle"""
    Dm = sp_cob(cx)
    P0 = cx["P0"]
    sg = np.ones(cx["nplaq"]); sg[P0] = -1.0

    def fun(e):
        Fl = Dm @ e
        return 4.0 * (1.0 - sg * np.cos(Fl / 2)).sum(), Dm.T @ (2.0 * sg * np.sin(Fl / 2))
    rng = np.random.default_rng(seed)
    e0 = np.zeros(cx["nplaq"]); e0[P0] = 2 * PI
    gauss = splinalg.lsqr(Dm, e0, atol=1e-13, btol=1e-13, iter_lim=50000)[0]
    best = None
    for s0 in [gauss] + [gauss * rng.uniform(0.2, 1.2) + 0.05 * rng.normal(size=cx["nlink"]) for _ in range(starts)] + \
            [0.01 * rng.normal(size=cx["nlink"])]:
        r_ = optimize.minimize(fun, s0, jac=True, method="L-BFGS-B", options={"maxiter": 20000, "ftol": 1e-15, "gtol": 1e-10})
        if best is None or r_.fun < best[0]:
            best = (float(r_.fun), r_.x)
    return best[0], float(abs((Dm @ best[1])[P0]) / (2 * PI))


def massive_c(d, m2, L=16):
    """Gaussian estimate with a link mass m (the doublet on the Higgs side): c = 1 - (2/d) <khat^2/(khat^2 + m^2)>"""
    k = 2 * np.sin(PI * np.arange(L) / L)
    g = np.meshgrid(*([k] * d), indexing="ij")
    k2 = sum(x ** 2 for x in g)
    return 1.0 - (2 / d) * float((k2 / (k2 + m2)).mean())


def pinned_pair_quadrature(g2, n=20000):
    """two Villain faces on one seam, outer links pinned: Z[twist]/Z by direct quadrature over the one free angle"""
    def vil(Fl):
        return sum(np.exp(-(Fl - 4 * PI * k) ** 2 / (2 * g2)) for k in range(-3, 4))
    eps = np.linspace(0, 4 * PI, n, endpoint=False)
    return float((vil(eps + 2 * PI) * vil(-eps)).sum() / (vil(eps) * vil(-eps)).sum())


def clock_table(steps=(1, 2, 3, 4, 8)):
    """c = g^2 <E^2> of the seam link against the number of electric steps per tick (Cartan circle, exact per mode).
    A transverse mode of spatial frequency w contributes w/sqrt(a^2 w^2 + 4) at the midpoint of the electric step,
    w sqrt(a^2 w^2 + 4)/4 at the midpoint of the face step and w/2 in continuous time, a = 1/steps."""
    spaces = {}
    for ds, L in ((3, 32), (2, 128)):
        k = 2 * np.sin(PI * np.arange(L) / L)
        g = np.meshgrid(*([k] * ds), indexing="ij")
        w = np.sqrt(sum(x ** 2 for x in g)).ravel()
        spaces["cubic space, %d dimensions" % ds] = (w, np.full(w.size, (ds - 1) / ds / w.size))
    lam = icosahedral_space_many_ticks()[2]
    spaces["icosahedral space"] = (np.sqrt(lam), np.full(lam.size, 1 / 30))
    out = {}
    for name, (w, wt) in spaces.items():
        rows = {}
        for n in steps:
            a = 1.0 / n
            rows[str(n)] = {"electric_step": float((wt * w / np.sqrt(a * a * w * w + 4)).sum()),
                            "face_step": float((wt * w * np.sqrt(a * a * w * w + 4) / 4).sum())}
        rows["continuous"] = {"electric_step": float((wt * w / 2).sum()), "face_step": float((wt * w / 2).sum())}
        out[name] = rows
    return out


def step_coefficients(cx, steps):
    """exact Gaussian coefficients on a finite lattice with `steps` electric steps per tick: time-like plaquettes of
    variance g^2 a, faces of variance g^2/a, a = 1/steps. Returns (c electric step, c face step)."""
    a = 1.0 / steps
    sig = np.where(cx["ptype"] == 1, a, 1.0 / a)
    rows, cols, vals = [], [], []
    for p in range(cx["nplaq"]):
        for k in range(cx["pn"][p]):
            rows.append(p); cols.append(int(cx["pl"][p, k])); vals.append(float(cx["ps"][p, k]) / math.sqrt(sig[p]))
    Dw = sparse.csr_matrix((vals, (rows, cols)), shape=(cx["nplaq"], cx["nlink"]))

    def proj(v):
        return Dw @ splinalg.lsqr(Dw, v, atol=1e-14, btol=1e-14, iter_lim=100000)[0]
    P0 = cx["P0"]
    e = np.zeros(cx["nplaq"]); e[P0] = 1.0
    c_el = (1.0 - proj(e.copy())[P0]) / a
    l = int(cx["pl"][P0, 3]) if cx["name"].startswith("cub") else int(cx["pl"][P0, 0])
    u = np.zeros(cx["nplaq"])
    nface = 0
    for q in range(cx["lp_ptr"][l], cx["lp_ptr"][l + 1]):
        p = int(cx["lp_idx"][q]); s = float(cx["ps"][p, int(cx["lp_pos"][q])])
        if p == P0:
            u[p] = 2 * PI * s
        elif cx["ptype"][p] == 0:
            u[p] = PI * s; nface += 1
    uh = u / np.sqrt(sig)
    return float(c_el), float(uh @ (uh - proj(uh.copy()))) / (4 * PI ** 2) + nface * a / 4
