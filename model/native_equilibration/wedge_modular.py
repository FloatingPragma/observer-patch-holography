"""Zero-temperature wedge modular energy of the native free staggered field.

The architecture is the explicitly echoed, half-filled free CAR extension,
not the unechoed interacting exchange model or its product-empty vacuum.
In a scalar transverse channel H has hopping -c/(2a) and potential
m*(-1)**n.  The halfspace is n>=1, cut 0|1; its proper site and bond
coordinates are a*(n-1/2) and a*n.  Symmetric lapse coupling, specified
independently of entropy, gives B={X,H}/2 (bond coefficient -J*n).

Eisler arXiv:2410.16433, Eqs64,67,70--72 and87 prove for this precise
infinite, gapped, filled-negative-band problem K=gamma*B, with
gamma=4*kappa*K(kappa')/c and kappa=(1+(ma/c)**2)**(-1/2).
Its Eq72 gives signed levels delta*(2*l-1/2), l in Z, m>0.  The absolute
levels delta*(j+1/2), j>=0, each occur ONCE; no SSH/Ising degeneracy is
imported.  The sign-reversed mass has the same entropy.  This theorem,
not mere commutation or finite numerical agreement, supplies the infinite
modular identity. H, B and K_half here denote ONE-PARTICLE kernels. For a
finite regional covariance, k_A=log((I-C_A)/C_A) and the full modular
Hamiltonian is -log(rho_A)=dGamma(k_A)+log(Z_A)*I. The infinite-half-chain
theorem does not identify a finite interval's k_A with gamma*B_N.

The infinite signed kernel has infinitely many negative levels. Its full
modular Hamiltonian is represented relative to the filled negative MODULAR
sea, or equivalently by particle-hole excitation energies |epsilon_l| and
the finite excitation partition sum. Naive empty-Fock dGamma(K_half) and
log(Z) are not separately finite. A fixed reference normal-ordering changes
only the scalar. For trace-zero differentiable state tangents with finite
modular-energy/entropy derivative, this scalar cancels and the first law is
delta S=gamma*delta<dGamma(B)>. Finite interval kernels are computed separately.

The source 8-mode preparation is a literal finite Gaussian CAR word, with
all256 amplitudes retained.  Its GLOBAL fidelity uses an independent exact
many-body energy-gap certificate.  Large-volume controls are local covariance
calculations, not executed large-volume source preparations or global trace
certificates.  No temperature, measured masses, Newton constant, counterterm,
curvature or entropy-stationarity rule is fitted.  Bare planar entropy is
not renormalized generalized entropy; the wedge is not a small ball.

Analytic bounds (spectral calculus, including arbitrary reference systems):
For 0<m<=Lambda and m^2<=H^2<=Lambda^2, set
r=(Lambda-m)/(Lambda+m), z=(Lambda^2+m^2-2H^2)/(Lambda^2-m^2).
The Legendre generating function gives
(H^2)^(-1/2)=2/(Lambda+m)*sum(r^j P_j(z)).
Its degree-k tail is <=r^(k+1)/m.  P_-=(I-H/|H|)/2 therefore has a
degree-(2k+1) approximation with error <=Lambda*r^(k+1)/(2m).
If ALL monomial walks of that length from a window agree in two systems,
their restricted covariance difference is <=Lambda*r^(k+1)/m.  Matching
only the entries inside the window is insufficient.  This is operator
norm of the local one-body covariance, NOT global many-body trace distance.

AGM(1,kappa) lies between kappa and (1+kappa)/2.  Hence
0<=1-gamma/(2*pi/c)<=(1-kappa)/(1+kappa)<=mu^2/4.
For |p_y|,|p_z|<=P this is uniformly <=a^2*((m/c)^2+2P^2)/4.
It controls the low transverse band, not the full Brillouin zone.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.special import expit, xlogy
import sympy as s

from tools import strict_json
from ..common_source import repair_field_attachment as attachment
from ..common_source.chiral_candidate import gaussian_source as gaussian
from ..common_source.halfline import pi_interval
from ..observer.source_scalar import compiler

ROOT = Path(__file__).resolve().parents[2]
SOURCE_GATE = 1e-8
TRACE_GATE = F(1, 10**6)
NUMERIC_GATE = 2e-9
MIN_PURITY_GAP = 1e-13
_GRID = 10**36


def _positive(value, name):
    if type(value) is not F or value <= 0:
        raise ValueError(f"positive exact Fraction {name} required")
    return value


def _integer(value, low, high, name):
    if type(value) is not int or not low <= value <= high:
        raise ValueError(f"bounded integer {name} required")
    return value


def _numerical_scales(mass, speed, spacing):
    """Finite arithmetic domain; analytic rational bounds have no such gate."""
    for value, name in ((mass, "mass"), (speed, "speed"), (spacing, "spacing")):
        _positive(value, name)
    try:
        values = float(mass), float(speed/spacing), float(spacing)
    except (OverflowError, ValueError) as exc:
        raise ValueError("resolved finite numerical chain scales required") from exc
    # Products under square roots and weighted matrix entries must stay
    # representable. These are arithmetic admissions, not physical cutoffs.
    if any(not math.isfinite(value) or not 1e-150 <= value <= 1e150 for value in values):
        raise ValueError("resolved finite numerical chain scales required")


def _root(value):
    if type(value) is not F or value < 0:
        raise ValueError("nonnegative exact radicand required")
    lower = F(math.isqrt(value.numerator*_GRID*_GRID//value.denominator), _GRID)
    return lower, lower if lower*lower == value else lower+F(1, _GRID)


def _ceil(value, grid=_GRID):
    return F(-(-value.numerator*grid//value.denominator), grid)


def _box(value):
    return (F(value[0].numerator*_GRID//value[0].denominator, _GRID), _ceil(value[1]))


def _power_upper(base, exponent):
    """Directed rational binary powering, with bounded receipt size."""
    if type(base) is not F or not 0 <= base <= 1:
        raise ValueError("exact power base in [0,1] required")
    _integer(exponent, 0, 10**9, "power")
    value = F(1)
    while exponent:
        if exponent & 1:
            value = _ceil(value*base)
        base = _ceil(base*base)
        exponent //= 2
    return value


def parameters(sector=0):
    _integer(sector, 0, 2, "native spectator sector")
    native = attachment.logical_parameters()
    return F(native["positive_occupation_gaps"][sector]), 2*F(native["exchange_strength"])


def _agm(modulus):
    """Rational interval AGM; each root and arithmetic endpoint is directed."""
    a, b = (F(1), F(1)), modulus
    for _ in range(16):
        a, b = ((a[0]+b[0])/2, (a[1]+b[1])/2), (_root(a[0]*b[0])[0], _root(a[1]*b[1])[1])
    return b[0], a[1]


def elliptic_data(mu_squared, speed=F(1)):
    """Exact rational enclosures; elliptic K takes a modulus, not its square."""
    _positive(mu_squared, "mu squared")
    _positive(speed, "coordinate speed")
    if not F(1, 10**16) <= mu_squared <= 64:
        raise ValueError("resolved positive numerical mass domain required")
    kappa = _root(1/(1+mu_squared))
    complement = _root(mu_squared/(1+mu_squared))
    ag, agc, pi = _agm(kappa), _agm(complement), pi_interval()
    ratio = _box((kappa[0]/ag[1], kappa[1]/ag[0]))
    delta = _box((2*pi[0]*agc[0]/ag[1], 2*pi[1]*agc[1]/ag[0]))
    gamma = _box((2*pi[0]/speed*ratio[0], 2*pi[1]/speed*ratio[1]))
    return {"kappa_interval": kappa, "gamma_interval": gamma,
            "continuum_prefactor_ratio_interval": ratio,
            "relative_prefactor_error_upper": mu_squared/4,
            "entanglement_spacing_interval": delta}


def band_bound(mass, speed, spacing, transverse_cutoff):
    for value, name in ((mass, "mass"), (speed, "speed"), (spacing, "spacing"), (transverse_cutoff, "band cutoff")):
        _positive(value, name)
    if spacing*transverse_cutoff > F(1, 4):
        raise ValueError("declared low transverse momentum band required")
    return spacing**2*((mass/speed)**2+2*transverse_cutoff**2)/4


def transverse_algebra():
    """Exact native corner representation; the diagonalizer acts only in y,z.

    For a cut between x cells, Fourier/diagonalization of transverse labels
    is local with respect to the two halfspaces. B_perp has eigenvalues
    +/-M twice each, giving four scalar chains, or two3D Dirac tastes.
    The onsite i**(-n) gauge changes real hopping to the oriented imaginary
    hopping; its inverse is i**n. Both preserve every spatial region.
    """
    x, z, identity = s.Matrix([[0, 1], [1, 0]]), s.diag(1, -1), s.eye(2)
    mass, fy, fz, fx = s.symbols("m fy fz fx", real=True)
    b = mass*s.kronecker_product(x, x)+fy*s.kronecker_product(z, identity)+fz*s.kronecker_product(x, z)
    alphas = (s.kronecker_product(z, s.eye(4)), s.kronecker_product(x, z, identity), s.kronecker_product(x, x, z))
    beta = s.kronecker_product(x, x, x)
    h = fx*alphas[0]+fy*alphas[1]+fz*alphas[2]+mass*beta
    if h != fx*s.kronecker_product(z, s.eye(4))+s.kronecker_product(x, b):
        raise ArithmeticError("transverse halfspace-preserving decomposition failed")
    if s.simplify(b*b-(mass*mass+fy*fy+fz*fz)*s.eye(4)) != s.zeros(4) or s.trace(b) != 0:
        raise ArithmeticError("four-chain transverse spectral multiplicity failed")
    return {"transverse_mass_squared": "m^2+f_y^2+f_z^2",
            "positive_transverse_mass_multiplicity": 2,
            "negative_transverse_mass_multiplicity": 2,
            "scalar_chains_per_reduced_transverse_momentum": 4,
            "Dirac_tastes": 2, "halfspace_preserved_by_transverse_unitary": True,
            "imaginary_hopping_onsite_gauge": "G_n=i^(-n); G H_real G^*=H_imaginary"}


def projector_bound(mass, norm_upper, degree):
    """Two covariance windows, with common gap and degree-(2k+1) walk agreement.

    The caller must establish the walk/spectral hypotheses; this helper is
    an arithmetic theorem evaluation, not a validator of arbitrary matrices.
    Lambda=m is the exact scalar spectral case, handled without z's 0/0.
    """
    _positive(mass, "common gap")
    _positive(norm_upper, "common spectral norm upper bound")
    _integer(degree, 0, 10**8, "Legendre degree")
    if norm_upper < mass:
        raise ValueError("spectral norm upper bound must dominate gap")
    return norm_upper/mass*_power_upper((norm_upper-mass)/(norm_upper+mass), degree+1)


def chain_bound(half, width, mass, speed, spacing=F(1)):
    """Actual open/infinite chain supplies the spectral and walk premises.

    Its hopping A anticommutes with the alternating sign, so H^2=m^2+A^2;
    the row-sum bound ||A||<=c/a holds with either open or infinite ends.
    An even half length matches the mass parity at the central cut. Every
    walk of degree<=N-width stays away from either omitted outer boundary.
    """
    _integer(half, 4, 2**20, "half length")
    _integer(width, 1, 4, "window width")
    if half % 2 or half <= width:
        raise ValueError("even half length larger than window required")
    for value, name in ((mass, "mass"), (speed, "speed"), (spacing, "spacing")):
        _positive(value, name)
    limit = _root(mass*mass+(speed/spacing)**2)[1]
    degree = (half-width-1)//2
    return {"degree": degree, "maximum_walk_length": 2*degree+1,
            "distance_to_omitted_boundary": half-width,
            "local_covariance_operator_norm_upper": projector_bound(mass, limit, degree),
            "common_gap": mass, "common_norm_upper": limit}


def orbitals(half, mass, speed, spacing=F(1), *, rows=None, positive=False):
    """Analytic open2N-chain sine/mass frame; no eigensystem gauge selection."""
    _integer(half, 2, 2**20, "half length")
    if half % 2 or type(positive) is not bool:
        raise ValueError("even half length and Boolean band required")
    _numerical_scales(mass, speed, spacing)
    if rows is None:
        if half > 2048:
            raise ValueError("full orbital frame exceeds finite numerical domain")
        rows = np.arange(1, 2*half+1)
    elif (type(rows) not in (tuple, list) or not rows or len(rows) > 4096
          or any(type(n) is not int or not 1 <= n <= 2*half for n in rows)):
        raise ValueError("bounded one-based orbital rows required")
    q = math.pi*np.arange(1, half+1)/(2*half+1)
    energy = float(speed/spacing)*np.cos(q)
    omega = np.hypot(energy, float(mass))
    a = np.sqrt((omega+energy)/(2*omega))
    b = float(mass)/np.sqrt(2*omega*(omega+energy))
    n = np.asarray(rows, dtype=int)[:, None]
    sine = math.sqrt(2/(2*half+1))*np.sin(n*q)
    return sine*((-b+a*(-1.)**(n+1)) if positive else (a+b*(-1.)**(n+1)))


def chain_matrix(count, mass, speed, spacing=F(1)):
    _integer(count, 2, 4096, "chain length")
    _numerical_scales(mass, speed, spacing)
    return np.diag(float(mass)*(-1.)**np.arange(1, count+1))+np.diag(np.full(count-1, -float(speed/(2*spacing))), 1)+np.diag(np.full(count-1, -float(speed/(2*spacing))), -1)


def boost_matrix(count, mass, speed, spacing=F(1)):
    """Independent symmetric proper-lapse derivative {X,H}/2 on n>=1."""
    h = chain_matrix(count, mass, speed, spacing)
    x = float(spacing)*(np.arange(1, count+1)-.5)
    return (x[:, None]+x[None, :])*h/2


def infinite_window(mass, speed, spacing=F(1), *, width=4, samples=2**18):
    """Independent Eq87 Fourier kernel, not a Gibbs fit to the boost."""
    _integer(width, 1, 4, "window width")
    _integer(samples, 128, 2**20, "periodic quadrature samples")
    _numerical_scales(mass, speed, spacing)
    mu = mass*spacing/speed
    if not F(1, 10**5) <= mu <= 1:
        raise ValueError("resolved gapped Fourier domain required")
    q = -math.pi/2+math.pi*(np.arange(samples)+.5)/samples
    cq = np.cos(q)
    omega = np.hypot(cq, float(mu))
    entries = {}
    for r in range(width):
        entries[r] = float(np.mean(np.cos(r*q)*(float(mu) if r % 2 == 0 else cq)/omega)/2)
    value = np.empty((width, width))
    for i in range(width):
        for j in range(width):
            value[i, j] = ((.5 if i == j else 0)-(-1)**(i+1)*entries[abs(i-j)]
                           if (i-j) % 2 == 0 else entries[abs(i-j)])
    # Period-pi analytic trapezoid theorem; strip sigma=mu/4 contains no
    # branch zeros. Re(cos^2(z)+mu^2)>=3mu^2/4, cosh(sigma)<4/3.
    # The resulting uniform integrand bound2/(mu*(1-sigma*(width-1)))
    # includes the1/2 from dq/(2pi). exp lower uses only positive terms.
    sigma = mu/4
    sup = 2/(mu*(1-sigma*(width-1)))
    argument = 2*sigma*samples
    term, exponential_lower = F(1), F(1)
    for j in range(1, 129):
        term *= argument/j
        exponential_lower += term
    truncation = width*sup/(exponential_lower-1)
    return value, truncation


def entropy_series(mu_squared, terms=128):
    """Eq72 spectrum, numerical entropy and a separate exact analytic tail."""
    _integer(terms, 16, 512, "entropy terms")
    data = elliptic_data(mu_squared)
    low, high = data["entanglement_spacing_interval"]
    delta = float((low+high)/2)
    x = delta*(np.arange(terms)+.5)
    entropy = float(np.sum(np.logaddexp(0, -x)+x*expit(-x)))
    # s_F(x)<=(1+x)e^-x, decreasing for x>=0.  q=(1+d)^-1>=e^-d.
    d, q = low, 1/(1+low)
    tail = _ceil(_power_upper(q, terms)*((1+d*(terms+F(1, 2)))/(1-q)+d*q/(1-q)**2))
    return {"entropy_diagnostic": entropy, "analytic_series_tail_upper": tail,
            "terms": terms, "absolute_level_multiplicity": 1,
            "absolute_levels": "delta*(n+1/2), n=0,1,..."}


def finite_entropy(half, mass, speed, spacing=F(1)):
    """Pure global projector cross-block singular values; no eigenvalue clips."""
    if half > 1024:
        raise ValueError("bounded independent entropy control required")
    occupied = orbitals(half, mass, speed, spacing)
    singular = np.linalg.svd(occupied[:half]@occupied[half:].T, compute_uv=False)
    if np.max(singular) > .5 or not np.isfinite(singular).all():
        raise ArithmeticError("unresolved pure-state cross-block spectrum")
    small = 2*singular**2/(1+np.sqrt(1-4*singular**2))
    return float(np.sum(-xlogy(small, small)-xlogy(1-small, 1-small)))


def finite_modular(covariance):
    """Faithful finite-region log only; unresolved purity gaps are refused."""
    if np.ma.isMaskedArray(covariance):
        raise ValueError("masked covariance inadmissible")
    raw = np.asarray(covariance)
    if raw.ndim != 2 or raw.shape[0] != raw.shape[1] or not 1 <= len(raw) <= 128 or raw.dtype.kind not in "iufc" or not np.isfinite(raw).all():
        raise ValueError("finite Hermitian CAR covariance required")
    value = raw.astype(complex)
    if np.linalg.norm(value-value.conj().T, 2) > 1e-13:
        raise ValueError("Hermitian covariance required")
    vals, vecs = np.linalg.eigh(value)
    gap = min(float(vals[0]), float(1-vals[-1]))
    residual = float(np.linalg.norm(value@vecs-vecs*vals, 2))
    if gap < MIN_PURITY_GAP or residual > gap/20:
        raise ValueError("regional purity gap unresolved; no logarithm or clipping")
    return (vecs*np.log((1-vals)/vals))@vecs.conj().T


def _fock(h):
    modes = len(h)
    value = np.zeros((1 << modes,)*2, complex)
    for before in range(1 << modes):
        for j in range(modes):
            if not before >> j & 1:
                continue
            removed = before ^ (1 << j)
            sign = (-1)**((before & ((1 << j)-1)).bit_count())
            for i in range(modes):
                if not removed >> i & 1:
                    value[removed | (1 << i), before] += h[i, j]*sign*(-1)**((removed & ((1 << i)-1)).bit_count())
    return value


def _slater(frame):
    value = np.zeros(1 << len(frame), complex)
    for occupied in combinations(range(len(frame)), frame.shape[1]):
        value[sum(1 << j for j in occupied)] = np.linalg.det(frame[list(occupied)])
    return value


def _covariance(state):
    annihilated = []
    for j in range(8):
        value = np.zeros(256, complex)
        for mask in range(256):
            if mask >> j & 1:
                value[mask ^ (1 << j)] = (-1)**((mask & ((1 << j)-1)).bit_count())*state[mask]
        annihilated.append(value)
    return np.asarray([[np.vdot(annihilated[j], annihilated[i]) for j in range(8)] for i in range(8)])


def source_preparation():
    """Actual adjacent local Fock words on a fresh blank8CAR register."""
    mass, speed = parameters()
    negative = orbitals(4, mass, speed)
    frame = np.column_stack((negative, orbitals(4, mass, speed, positive=True)))
    work, calls = frame.copy(), []
    for column in range(8):
        for row in range(7, column, -1):
            a, b = work[row-1, column], work[row, column]
            radius = math.hypot(a, b)
            angle = math.atan2(b, a) if radius else 0.
            rotation = np.array([[math.cos(angle), math.sin(angle)], [-math.sin(angle), math.cos(angle)]])
            work[[row-1, row]] = rotation@work[[row-1, row]]
            calls.append((row-1, row, angle))
    if np.linalg.norm(work-np.diag(np.diag(work)), 2) > SOURCE_GATE:
        raise ArithmeticError("analytic orbital Givens elimination failed")
    state = np.zeros((256, 1), complex)
    state[0, 0] = 1
    counts = {key: 0 for key in gaussian._ADDITIVE_RESOURCES}
    def apply(gate, sites):
        nonlocal state
        state = attachment._little_apply(gate.fock_unitary, sites, 8, state)
        for key in counts:
            counts[key] += gate.resources[key]
    # Two pairing pulses prepare four charges. exp(-i*pi/2 pair) gives
    # (-i)^2=-1 on |1111>; retain that known global phase in the target.
    for sites in ((0, 1), (2, 3)):
        apply(gaussian.local_gate("pairing", 1, math.pi/2), sites)
    for mode, sign in enumerate(np.diag(work)):
        word = (compiler.pulse("I", (0. if sign >= 0 else math.pi)/2, [0]),
                compiler.pulse("Z", -(0. if sign >= 0 else math.pi)/2, [0]))
        actual = compiler.execute(word, 1)
        state = attachment._little_apply(actual, (mode,), 8, state)
        for key in counts:
            counts[key] += compiler.resources(word)[key]
    for left, right, angle in reversed(calls):
        apply(gaussian.local_gate("hopping", 1j, -angle), (left, right))
    target = -_slater(negative)
    output = state[:, 0]
    if np.linalg.norm(output-target) > SOURCE_GATE:
        raise ArithmeticError("literal native vacuum word exceeds unchanged source-state gate")
    certify_source(output)
    # QR angles are deterministic analytic inputs; integrated norm remains a
    # numerical cost diagnostic. Use a fixed pi envelope for portable receipts.
    costs = {key: value for key, value in counts.items() if key not in ("zero_angle_requests", "integrated_processor_norm_upper")}
    return output, target, costs


@lru_cache(maxsize=1)
def _ground_intervals():
    mass, speed = parameters()
    h = s.diag(*(s.Rational(mass)*(-1)**n for n in range(1, 9)))
    for n in range(7):
        h[n, n+1] = h[n+1, n] = -s.Rational(speed)/2
    intervals = s.Poly(h.charpoly().as_expr()).intervals(eps=s.Rational(1, 10**36))
    roots = tuple((F(a), F(b)) for (a, b), multiplicity in intervals for _ in range(multiplicity))
    if len(roots) != 8 or any(a <= 0 <= b for a, b in roots):
        raise ArithmeticError("exact finite-chain Fermi gap unresolved")
    ground = sum((a for a, b in roots if b < 0), F()), sum((b for a, b in roots if b < 0), F())
    gap = min(min(abs(a), abs(b)) for a, b in roots)
    return ground, gap


def certify_source(state):
    """Exact binary Rayleigh/gap proof for the retained raw source ray."""
    if np.ma.isMaskedArray(state):
        raise ValueError("masked source state inadmissible")
    raw = np.asarray(state)
    if raw.shape != (256,) or raw.dtype != np.dtype(complex) or not np.isfinite(raw).all():
        raise ValueError("original finite complex256 source amplitudes required")
    real, imag = tuple(F(float(x)) for x in raw.real), tuple(F(float(x)) for x in raw.imag)
    norm = sum((x*x+y*y for x, y in zip(real, imag, strict=True)), F())
    if norm == 0:
        raise ValueError("nonzero source ray required")
    mass, speed = parameters()
    expectation = F()
    for mask in range(256):
        expectation += mass*sum(((-1)**(j+1) for j in range(8) if mask >> j & 1), 0)*(real[mask]**2+imag[mask]**2)
        for j in range(7):
            if mask >> j & 1 and not mask >> (j+1) & 1:
                other = mask ^ (3 << j)
                expectation -= speed*(real[mask]*real[other]+imag[mask]*imag[other])
    ground, gap = _ground_intervals()
    infidelity = (expectation/norm-ground[0])/gap
    if infidelity < 0:
        raise ArithmeticError("exact ground-energy enclosure contradicted")
    trace = abs(norm-1)/2+norm*_root(infidelity)[1]
    if trace >= TRACE_GATE:
        raise ValueError("literal source global vacuum fidelity certificate failed")
    return {"raw_norm_squared_exact": norm, "Rayleigh_energy_exact": expectation/norm,
            "ground_energy_interval": ground, "gap_lower": gap,
            "infidelity_upper": infidelity, "global_raw_half_trace_norm_upper": trace}


def source_first_law(state):
    """Actual fullFock partial trace and non-Gaussian mixture entropy secants."""
    certify_source(state)
    tensor = state.reshape(16, 16)
    # The physical comparison uses the ray's explicit norm. Original binary
    # source arrays are retained unchanged; this is not a source correction.
    norm = float(np.vdot(state, state).real)
    rho = (tensor@tensor.conj().T)/norm
    vals, vecs = np.linalg.eigh(rho)
    if vals[0] < MIN_PURITY_GAP:
        raise ValueError("source regional fullFock purity gap unresolved")
    log_rho = (vecs*np.log(vals))@vecs.conj().T
    phase = np.exp(1j*.25*(np.arange(16) & 1))
    moved = phase[:, None]*rho*phase.conj()[None, :]
    tangent = moved-rho
    derivative = float(-np.trace(log_rho@tangent).real)
    entropy = float(-np.dot(vals, np.log(vals)))
    secants = []
    for epsilon in (1e-2, 1e-3, 1e-4):
        spectrum = np.linalg.eigvalsh(rho+epsilon*tangent)
        if spectrum[0] <= 0:
            raise ArithmeticError("physical mixture lost faithful support")
        secants.append(float((-np.dot(spectrum, np.log(spectrum))-entropy)/epsilon))
    mass, speed = parameters()
    c = _covariance(state)[4:, 4:]/norm
    modular = finite_modular(c)
    u = np.diag([np.exp(.25j), 1, 1, 1])
    dc = u@c@u.conj().T-c
    quadratic_derivative = float(np.trace(modular@dc).real)
    gamma = float(sum(elliptic_data((mass/speed)**2, speed)["gamma_interval"])/2)
    bw_response = gamma*float(np.trace(boost_matrix(4, mass, speed)@dc).real)
    if abs(quadratic_derivative-derivative) > NUMERIC_GATE:
        raise ArithmeticError("independent fullFock and CAR first laws disagree")
    return {"entropy": entropy, "full_Fock_first_law_derivative": derivative,
            "quadratic_first_law_derivative": quadratic_derivative,
            "mixture_secants": secants, "mixture_parameters": [.01, .001, .0001],
            "infinite_wedge_weight_finite_chain_response": bw_response,
            "finite_interval_BW_equality_claimed": False,
            "normalized_ray_used_only_for_state_function_comparison": True,
            "original_source_amplitudes_modified": False}


def planar_entropy(sector=0, transverse_cells=8, spacing=F(1)):
    """Finite transverse torus, infinite normal halfchain; bare entropy/area."""
    _integer(transverse_cells, 4, 32, "even transverse cells")
    if transverse_cells % 4:
        raise ValueError("transverse cell count divisible by four required for the declared reduced-zone grid")
    _positive(spacing, "proper lattice spacing")
    mass, speed = parameters(sector)
    entropy, max_tail = 0., F()
    for ny in range(-transverse_cells//4, transverse_cells//4):
        for nz in range(-transverse_cells//4, transverse_cells//4):
            # p is reduced sublattice momentum: range[-pi/(2a),pi/(2a)).
            sine_squared = math.sin(2*math.pi*ny/transverse_cells)**2+math.sin(2*math.pi*nz/transverse_cells)**2
            mu2 = (mass*spacing/speed)**2+F(sine_squared)
            entry = entropy_series(mu2)
            entropy += 4*entry["entropy_diagnostic"]
            max_tail += 4*entry["analytic_series_tail_upper"]
    area = transverse_cells**2*spacing**2
    return {"transverse_cells_per_axis": transverse_cells,
            "proper_area_exact": area, "scalar_chains_per_reduced_transverse_momentum": 4,
            "Dirac_tastes": 2, "entropy_diagnostic": entropy,
            "bare_entropy_per_proper_area_diagnostic": entropy/float(area),
            "analytic_level_series_tail_upper": max_tail,
            "transverse_momentum_trigonometry_is_numerical": True,
            "renormalized_area_coefficient_or_Newton_constant": False}


def finite_planar_control():
    """Independent128-site3D projector versus sixteen scalar-chain cuts.

    Open normal direction8, periodic transverse4x4, cut4|4. The native
    real-space operator is used directly; no transverse spectrum is used
    to construct its projector. Only this finite numerical check uses eigh.
    """
    from ..common_source import native_staggered_field as field
    coordinates = field._coordinates((8, 4, 4))
    h = np.asarray(field.lattice_operator((8, 4, 4)), complex)
    # Remove only the x wrap bonds. Transverse boundaries stay periodic.
    for i, x in enumerate(coordinates):
        if x[0] != 0:
            continue
        j = coordinates.index((7, x[1], x[2]))
        h[i, j] = h[j, i] = 0
    eigenvalues, eigenvectors = np.linalg.eigh(h)
    occupied = eigenvectors[:, eigenvalues < 0]
    if occupied.shape != (128, 64):
        raise ArithmeticError("native3D half-filled spectral count changed")
    cross = occupied[:64]@occupied[64:].conj().T
    singular = np.linalg.svd(cross, compute_uv=False)
    if np.max(singular) > .5:
        raise ArithmeticError("3D finite projector purity calculation unresolved")
    small = 2*singular**2/(1+np.sqrt(1-4*singular**2))
    direct = float(np.sum(-xlogy(small, small)-xlogy(1-small, 1-small)))
    mass, speed = parameters()
    reduced = 0.
    for ny in (-1, 0):
        for nz in (-1, 0):
            transverse_mass = F(math.hypot(float(mass), float(speed)*math.sin(ny*math.pi/2), float(speed)*math.sin(nz*math.pi/2)))
            reduced += 4*finite_entropy(4, transverse_mass, speed)
    if abs(direct-reduced) > NUMERIC_GATE:
        raise ArithmeticError("native3D finite cut disagrees with four scalar-chain transverse reduction")
    return {"lattice_shape": [8, 4, 4], "normal_boundary": "open", "transverse_boundary": "periodic",
            "normal_cut": "4|4", "full_one_body_modes": 128,
            "occupied_modes": 64, "scalar_chains": 16,
            "direct_native_projector_entropy_diagnostic": direct,
            "four_chain_reduction_entropy_diagnostic": reduced,
            "large_Fock_source_preparation_executed": False}


def bindings():
    names = set(attachment.source_bindings())|{
        "model/native_equilibration/wedge_modular.py", "model/native_equilibration/test_wedge_modular.py",
        "model/common_source/native_staggered_field.py", "model/common_source/test_native_staggered_field.py",
        "tools/strict_json.py"}
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def run():
    pins = bindings()
    state, target, resources = source_preparation()
    arrays = {"source_state": state, "source_target": target, "source_covariance": _covariance(state)}
    controls = []
    for sector in range(3):
        mass, speed = parameters(sector)
        infinite, quadrature = infinite_window(mass, speed)
        arrays[f"infinite_window_sector{sector}"] = infinite
        local = []
        for half in (32, 512, 8192, 32768):
            frame = orbitals(half, mass, speed, rows=list(range(half+1, half+5)))
            covariance = frame@frame.T
            error = float(np.linalg.norm(covariance-infinite, 2))
            bound = chain_bound(half, 4, mass, speed)
            if error > float(bound["local_covariance_operator_norm_upper"]+quadrature)+NUMERIC_GATE:
                raise ArithmeticError("finite/thermodynamic local covariance exceeds analytic IR bound")
            local.append({"half_length": half, "local_covariance_error_diagnostic": error, **bound})
        entropy = entropy_series((mass/speed)**2)
        entropies = [{"half_length": half, "entropy_diagnostic": finite_entropy(half, mass, speed)} for half in (8, 32, 128, 512)]
        # sigma_e=(1-e)sigma+e U sigma U*, U=exp(i theta n_1), is a
        # legitimate trace-zero tangent. Its B difference uses only bond1|2:
        # Delta B=2J(1-cos theta) C_12>0. The mixture need not be Gaussian.
        # This is a theorem-derived infinite first law; independent entropy
        # secants are computed on the actually prepared finite state above.
        theta = .25
        delta_b = float(float(speed)*(1-math.cos(theta))*infinite[0, 1])
        gamma = float(sum(elliptic_data((mass/speed)**2, speed)["gamma_interval"])/2)
        phase_control = {"theta_exact": "1/4", "unitary": "exp(i*theta*n_1)",
            "physical_tangent": "sigma_e=(1-e)*sigma+e*U*sigma*U^*, e>=0",
            "independent_energy_moment_change_diagnostic": delta_b,
            "infinite_modular_first_law_derivative_diagnostic": gamma*delta_b,
            "positive_energy_moment_identity": "Delta B=2J*(1-cos(theta))*C_12>0",
            "continuum_relative_error_upper_exact": (mass/speed)**2/4,
            "numerical_infinite_entropy_secants_claimed": False}
        controls.append({"sector": sector, "mass_exact": mass, "speed_exact": speed,
            "mu_squared_exact": (mass/speed)**2, **elliptic_data((mass/speed)**2, speed),
            "infinite_entropy": entropy, "independent_finite_entropy_controls": entropies,
            "Fourier_quadrature_analytic_operator_tail_upper": quadrature,
            "local_covariance_controls": local,
            "infinite_local_phase_mixture_first_law": phase_control,
            "planar_entropy": planar_entropy(sector),
            "low_transverse_band_relative_error_upper_a1_16_P1": band_bound(mass, speed, F(1, 16), F(1))})
    # Independent gapped control makes the exact level offset visible well
    # before native low-mass correlation lengths fit in small finite chains.
    mass, speed = parameters()
    gapped = {"mu_exact": "1/4", **entropy_series(F(1, 16)),
              "independent_half32_entropy_diagnostic": finite_entropy(32, speed/4, speed)}
    finite_controls = []
    for half in (4, 8, 16, 32, 64, 128, 256):
        occupied = orbitals(half, mass, speed)
        covariance = occupied[half:]@occupied[half:].T
        gamma = float(sum(elliptic_data((mass/speed)**2, speed)["gamma_interval"])/2)
        b = boost_matrix(half, mass, speed)
        eigen, frame = eigh_tridiagonal(np.diag(b), np.diag(b, 1))
        prediction = (frame*expit(-gamma*eigen))@frame.T
        try:
            finite_modular(covariance)
            status = "resolved finite regional modular logarithm"
        except ValueError:
            status = "purity gap unresolved; finite modular logarithm refused"
        finite_controls.append({"half_length": half, "finite_log_status": status,
            "first_four_covariance_BW_error_diagnostic": float(np.linalg.norm(covariance[:4, :4]-prediction[:4, :4], 2))})
    report = {"schema": "native-free-wedge-modular-v1", "source_bindings": pins,
        "primary_theorem": {"url": "https://arxiv.org/abs/2410.16433", "equations": [64, 67, 70, 71, 72, 87],
            "domain": "infinite gapped free scalar staggered chain, filled negative band, halfspace n>=1, cut0|1"},
        "independent_geometry": {"site_position": "a*(n-1/2)", "bond_position": "a*n", "lapse_coupling": "H[N]={N,H}/2", "energy_first_moment": "B={X,H}/2"},
        "exact_result": "K_half=gamma*B; 0<=1-gamma/(2*pi/c)<=mu^2/4",
        "native_transverse_reduction": transverse_algebra(),
        "independent_finite_planar_control": finite_planar_control(),
        "arithmetic": {"analytic_inequalities_evaluated_with_directed_rationals": True,
            "Fourier_and_entropy_partial_sums_outward_certified": False,
            "floating_diagnostic_absolute_admission": NUMERIC_GATE,
            "source_global_energy_certificate_exact_for_original_binary_arrays": True},
        "first_law_domain": "trace-zero differentiable state tangents with finite modular-energy/entropy derivative; finite-support phase mixtures are explicit controls",
        "source_preparation": {"CAR_modes": 8, "retained_Fock_amplitudes": 256, "filled_negative_modes": 4,
            "representation": "real-hopping scalar transverse chain; onsite G_n=i^(-n) identifies the oriented imaginary-hopping representation and preserves the cut",
            "source_arrays_are_in_real_hopping_representation": True,
            "blank_and_pair_seed_services_supplied": True, "native_CE_cooling_or_vacuum_selection": False,
            "source_state_error_admission": SOURCE_GATE, "global_raw_half_trace_norm_strict_upper": str(TRACE_GATE),
            "global_certificate": "exact binary Rayleigh quotient, exact algebraic ground-energy enclosure and many-body Fermi gap",
            "own_binary_certificate": certify_source(state),
            "resources": resources, "per_request_phase_angle_envelope": "pi", "service_error_accumulation": "at most sum of admitted per-request channel errors; finite replay accuracy is not a physical service theorem"},
        "source_full_Fock_first_law": source_first_law(state),
        "native_sector_controls": controls, "independent_gapped_level_control": gapped,
        "finite_interval_controls": finite_controls,
        "scope": {"zero_temperature_reference": True, "supplied_inverse_temperature": False,
            "free_echo_and_parity_architecture_required": True,
            "four_scalar_chains_and_two_Dirac_tastes": True,
            "native_spectator_sector_values_are_three_independent_species": False,
            "source_executed_scalar_channel_sector0_at_zero_transverse_momentum": True,
            "source_executed_large_volume_vacuum": False,
            "local_covariance_bound_implies_global_many_body_trace_bound": False,
            "full_Brillouin_zone_continuum_BW_prefactor": False,
            "small_ball_modular_or_curved_geometry_stationarity": False,
            "universal_area_density_or_Newton_constant": False,
            "Einstein_premise_closed": False,
            "measured_parameter_inputs": False}}
    if bindings() != pins:
        raise RuntimeError("wedge source changed during execution")
    return _json(report), arrays


def _json(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, dict):
        return {key: _json(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json(item) for item in value]
    return value


def compare(retained, fresh, location="root"):
    if type(retained) is not type(fresh):
        raise ValueError("wedge replay type changed: "+location)
    if isinstance(fresh, dict):
        if retained.keys() != fresh.keys():
            raise ValueError("wedge receipt fields changed: "+location)
        for key in fresh:
            compare(retained[key], fresh[key], location+"."+key)
    elif isinstance(fresh, list):
        if len(retained) != len(fresh):
            raise ValueError("wedge list length changed: "+location)
        for index, (left, right) in enumerate(zip(retained, fresh, strict=True)):
            compare(left, right, location+f"[{index}]")
    elif isinstance(fresh, float):
        if not math.isfinite(retained) or not math.isfinite(fresh):
            raise ValueError("nonfinite wedge diagnostic")
        # Fixed gates are exact fields. Numerical backend diagnostics are
        # portable but never substitute for the retained binary source gates.
        if "admission" in location:
            good = retained == fresh
        else:
            good = math.isclose(retained, fresh, abs_tol=NUMERIC_GATE, rel_tol=NUMERIC_GATE)
        if not good:
            raise ValueError("wedge numerical replay changed: "+location)
    elif retained != fresh:
        raise ValueError("wedge exact replay changed: "+location)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--output", type=Path)
    group.add_argument("--check", type=Path)
    args = parser.parse_args()
    report, arrays = run()
    path = args.output or args.check
    artifact = path.with_suffix(".npz")
    if args.output:
        if path.exists() or artifact.exists():
            parser.error("fresh output/array paths required")
        with artifact.open("xb") as stream:
            np.savez_compressed(stream, **arrays)
        report["artifact"] = {"filename": artifact.name, "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest()}
        path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+"\n")
    else:
        old = strict_json.load(path)
        if old.pop("artifact", None) != {"filename": artifact.name, "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest()}:
            raise ValueError("wedge artifact custody changed")
        with np.load(artifact, allow_pickle=False) as saved:
            if len(saved.files) != len(set(saved.files)) or set(saved.files) != set(arrays):
                raise ValueError("wedge array menu changed")
            for key, value in arrays.items():
                archived = saved[key]
                if archived.shape != value.shape or archived.dtype != value.dtype or not np.isfinite(archived).all() or np.linalg.norm(archived-value) > SOURCE_GATE:
                    raise ValueError("wedge retained binary replay failed: "+key)
            retained_certificate = _json(certify_source(saved["source_state"]))
            compare(old["source_preparation"]["own_binary_certificate"], retained_certificate)
            if np.linalg.norm(saved["source_state"]-arrays["source_target"]) > SOURCE_GATE or np.linalg.norm(saved["source_target"]-arrays["source_target"]) > SOURCE_GATE:
                raise ValueError("retained source exceeds independent target gate")
            if np.linalg.norm(saved["source_covariance"]-_covariance(saved["source_state"]), 2) > NUMERIC_GATE:
                raise ValueError("retained covariance does not belong to retained source state")
            retained_law = _json(source_first_law(saved["source_state"]))
            compare(old["source_full_Fock_first_law"], retained_law)
        # Exact binary certificates describe their OWN retained amplitudes.
        # Fresh and archived arrays independently pass the same fixed gates;
        # a platform's different binary rounding is not forged exact science.
        report["source_preparation"]["own_binary_certificate"] = retained_certificate
        compare(old, report)
        print("native free wedge modular/stress, entropy levels and source vacuum replay: PASS")


if __name__ == "__main__":
    main()
