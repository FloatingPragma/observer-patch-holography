"""Compressed canonical scalar Hamiltonian on a supplied open cubic cell.

The eight sites are the lexicographically ordered coordinates in {0, 1}^3,
with site zero at the origin.  Occupations are little endian: the flat index
is sum(n_j * product(dimensions[:j])).  The canonical infinite oscillator
operators are compressed *after* forming q^2, p^2 and q^4.  Squaring a
truncated q instead would change the Hamiltonian at its cutoff boundary.

The reference action, spatial lattice and model time are supplied.  Finite
matrices do not obey the exact canonical commutator, and this finite cell
does not establish a Lorentz invariant or interacting continuum limit.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import product
from math import isfinite, prod, sqrt
from numbers import Integral, Real

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply


CUBE_SITES = tuple(product(range(2), repeat=3))
CUBE_EDGES = tuple(
    (i, j) for i, x in enumerate(CUBE_SITES)
    for j, y in enumerate(CUBE_SITES) if i < j
    and sum(abs(a-b) for a, b in zip(x, y)) == 1
)


def _real(value, name, *, positive=False, nonnegative=False):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Real)
            or not isfinite(value) or (positive and value <= 0)
            or (nonnegative and value < 0)):
        qualification = "positive " if positive else "nonnegative " if nonnegative else ""
        raise ValueError(f"{name} must be a finite {qualification}real")
    return float(value)


def _dimension(value):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral)
            or value < 1):
        raise ValueError("each occupation dimension must be a positive integer")
    return int(value)


def _numeric(value, name):
    if np.ma.isMaskedArray(value):
        raise ValueError(f"{name} cannot be masked")
    result = np.asarray(value)
    if result.dtype.kind not in "fiuc" or not np.isfinite(result).all():
        raise ValueError(f"{name} must contain finite numeric entries")
    return result


@dataclass(frozen=True)
class ScalarParameters:
    """Supplied canonical action and oscillator basis, with unit hbar.

    ``coupling`` denotes lambda, so the coefficient of q^4 is
    g=lambda/(24*a^3), while kappa=c^2/a^2.  The basis frequency ``omega``
    is not the interacting excitation gap.  Its default value two makes
    the quadratic onsite term diagonal for the other default parameters.
    """

    mass: float = 1.0
    speed: float = 1.0
    spacing: float = 1.0
    coupling: float = 1.0
    omega: float = 2.0

    def __post_init__(self):
        for name in ("mass", "speed", "spacing", "omega"):
            object.__setattr__(self, name, _real(getattr(self, name), name, positive=True))
        object.__setattr__(self, "coupling", _real(self.coupling, "coupling", nonnegative=True))
        # Reject finite inputs whose derived action cannot be represented.
        try:
            scales = (self.mass**2, self.kappa, self.g, self.omega**2)
        except (OverflowError, ZeroDivisionError) as exc:
            raise ValueError("derived scalar action coefficients must be finite") from exc
        if not all(isfinite(x) for x in scales) or self.kappa <= 0:
            raise ValueError("derived scalar action coefficients must be finite and kappa positive")

    @property
    def kappa(self):
        return (self.speed/self.spacing)**2

    @property
    def g(self):
        return self.coupling/(24*self.spacing**3)


@lru_cache(maxsize=32)
def _local_operators(dimension, omega):
    # Four ladder steps suffice for every matrix element of q^4 with both
    # endpoints below dimension.  Extra levels are only assembly workspace.
    extended = dimension+4
    annihilation = np.diag(np.sqrt(np.arange(1, extended, dtype=float)), 1)
    q = (annihilation+annihilation.T)/sqrt(2*omega)
    p = 1j*sqrt(omega/2)*(annihilation.T-annihilation)
    q2, p2 = q@q, (p@p).real
    operators = {
        "a": annihilation[:dimension, :dimension].copy(),
        "q": q[:dimension, :dimension].copy(),
        "p": p[:dimension, :dimension].copy(),
        "q2": q2[:dimension, :dimension].copy(),
        "p2": p2[:dimension, :dimension].copy(),
        "q4": (q2@q2)[:dimension, :dimension].copy(),
        "number": np.diag(np.arange(dimension, dtype=float)),
    }
    for matrix in operators.values():
        matrix.flags.writeable = False
    return operators


def local_operators(dimension, omega=2.0):
    """Return independent arrays for compressed oscillator polynomials.

    Powers refer to the infinite oscillator before compression, not powers
    of the returned finite q or p.  The annihilation matrix itself is the
    usual compressed lowering operator.
    """
    dimension = _dimension(dimension)
    omega = _real(omega, "oscillator frequency", positive=True)
    return {name: matrix.copy() for name, matrix in _local_operators(dimension, omega).items()}


def free_dispersion(wavevector, mass=1.0, speed=1.0, spacing=1.0):
    """Infinite/periodic free spatial lattice frequency for three-momenta.

    Returns sqrt(m^2 + 4*c^2/a^2 sum_i sin^2(a*k_i/2)).  Open-cell modes
    instead diagonalize ``ScalarLattice.stiffness``.  This formula concerns
    the untruncated quadratic action, not a finite occupation spectrum.
    """
    parameters = ScalarParameters(mass=mass, speed=speed, spacing=spacing, coupling=0)
    k = _numeric(wavevector, "three-momentum")
    if k.ndim < 1 or k.shape[-1] != 3 or np.iscomplexobj(k):
        raise ValueError("three-momenta must have final dimension three and real entries")
    with np.errstate(over="ignore", invalid="ignore"):
        phase = parameters.spacing*k/2
    if not np.isfinite(phase).all():
        raise ValueError("lattice momentum phases must be finite")
    result = np.sqrt(parameters.mass**2+4*parameters.kappa*np.sum(np.sin(phase)**2, axis=-1))
    return float(result) if result.ndim == 0 else result


class ScalarLattice:
    """Sparse compression P H P of the canonical repulsive scalar action.

    H = 1/2 sum p_x^2 + m^2/2 sum q_x^2
        + kappa/2 sum_edges (q_x-q_y)^2 + g sum q_x^4.

    Per-site cutoffs may differ.  They declare a rectangular tensor-product
    compression, not an invariant sector of the infinite interacting field.
    No cutoff-error estimate or normalization is inserted by this class.
    """

    def __init__(self, dimensions=(2,)*8, parameters=None):
        if not isinstance(dimensions, (tuple, list)) or len(dimensions) != 8:
            raise ValueError("the open cubic cell needs exactly eight occupation dimensions")
        self.dimensions = tuple(_dimension(d) for d in dimensions)
        self.dimension = prod(self.dimensions)
        self.parameters = ScalarParameters() if parameters is None else parameters
        if not isinstance(self.parameters, ScalarParameters):
            raise ValueError("parameters must be a ScalarParameters instance")
        self.sites, self.edges = CUBE_SITES, CUBE_EDGES
        self.volume, self.origin = 8, 0
        self.site_index = {site: i for i, site in enumerate(self.sites)}
        self.degree = np.bincount(np.asarray(self.edges).ravel(), minlength=self.volume)
        self.strides = tuple(prod(self.dimensions[:i]) for i in range(self.volume))
        self.operators = tuple(local_operators(d, self.parameters.omega) for d in self.dimensions)
        adjacency = np.zeros((self.volume, self.volume))
        for x, y in self.edges:
            adjacency[x, y] = adjacency[y, x] = 1
        self.stiffness = (self.parameters.mass**2*np.eye(self.volume)
                          + self.parameters.kappa*(np.diag(self.degree)-adjacency))
        self.hamiltonian = self._hamiltonian()

    def _site(self, site):
        if isinstance(site, Integral) and not isinstance(site, (bool, np.bool_)):
            if 0 <= site < self.volume:
                return int(site)
        elif (isinstance(site, (tuple, list)) and len(site) == 3
              and all(isinstance(x, Integral) and not isinstance(x, (bool, np.bool_)) for x in site)):
            if tuple(site) in self.site_index:
                return self.site_index[tuple(site)]
        raise ValueError("site must be an index or a coordinate in the cubic cell")

    def _sites(self, sites):
        # A three-element sequence here denotes three site indices; use
        # _site / onsite for a coordinate or pass a coordinate inside a list.
        if isinstance(sites, Integral) and not isinstance(sites, (bool, np.bool_)):
            sites = (sites,)
        if not isinstance(sites, (tuple, list)) or not sites:
            raise ValueError("a nonempty sequence of distinct sites is required")
        selected = tuple(self._site(site) for site in sites)
        if len(set(selected)) != len(selected):
            raise ValueError("local sites must be distinct")
        return selected

    def _state(self, state):
        state = _numeric(state, "state")
        if state.ndim not in (1, 2) or state.shape[0] != self.dimension:
            raise ValueError("state must have one amplitude row per tensor occupation")
        return np.asarray(state, dtype=complex)

    def _matrix(self, matrix, sites):
        if sparse.issparse(matrix):
            if matrix.dtype.kind not in "fiuc" or not np.isfinite(matrix.data).all():
                raise ValueError("local matrix must have finite numeric entries")
            matrix = matrix.tocsr()
        else:
            matrix = _numeric(matrix, "local matrix")
        size = prod(self.dimensions[site] for site in sites)
        if matrix.shape != (size, size):
            raise ValueError("local matrix shape must match the selected occupation dimensions")
        return matrix

    def _tensor_product(self, operators):
        result = sparse.csr_matrix([[1.]])
        for site in reversed(range(self.volume)):
            factor = operators.get(site)
            if factor is None:
                factor = sparse.eye(self.dimensions[site], format="csr")
            result = sparse.kron(result, sparse.csr_matrix(factor), format="csr")
        result.eliminate_zeros()
        return result

    def onsite(self, site):
        """Compressed onsite summand, including incident gradient diagonals."""
        site = self._site(site)
        p, op = self.parameters, self.operators[site]
        return .5*op["p2"]+.5*(p.mass**2+p.kappa*self.degree[site])*op["q2"]+p.g*op["q4"]

    def _hamiltonian(self):
        result = sparse.csr_matrix((self.dimension, self.dimension), dtype=float)
        for site in range(self.volume):
            result += self._tensor_product({site: self.onsite(site)})
        for x, y in self.edges:
            result -= self.parameters.kappa*self._tensor_product(
                {x: self.operators[x]["q"], y: self.operators[y]["q"]})
        result.eliminate_zeros()
        result.sort_indices()
        if not np.isfinite(result.data).all():
            raise ValueError("compressed action entries must be finite")
        return result

    def embed_local(self, matrix, sites):
        """Embed a local matrix; the first supplied site is least significant.

        This also accepts nonadjacent and reversed site lists.  The sparse
        construction retains every spectator configuration without projecting
        or normalizing it.
        """
        sites = self._sites(sites)
        matrix = sparse.coo_matrix(self._matrix(matrix, sites))
        if len(sites) == 1:
            return self._tensor_product({sites[0]: matrix})
        spectators = tuple(site for site in range(self.volume) if site not in sites)
        offsets = np.zeros(prod(self.dimensions[site] for site in spectators), dtype=np.int64)
        remainder = np.arange(len(offsets), dtype=np.int64)
        for site in spectators:
            offsets += (remainder % self.dimensions[site])*self.strides[site]
            remainder //= self.dimensions[site]

        def expand(indices):
            remainder = indices.copy()
            result = np.zeros(len(indices), dtype=np.int64)
            for site in sites:
                result += (remainder % self.dimensions[site])*self.strides[site]
                remainder //= self.dimensions[site]
            return result

        rows = (expand(matrix.row)[:, None]+offsets).ravel()
        cols = (expand(matrix.col)[:, None]+offsets).ravel()
        values = np.repeat(matrix.data, len(offsets))
        result = sparse.coo_matrix((values, (rows, cols)),
                                  shape=(self.dimension, self.dimension)).tocsr()
        result.eliminate_zeros()
        return result

    def apply_local(self, state, sites, matrix):
        """Apply a local matrix to amplitudes or columns, with no renormalization."""
        state = self._state(state)
        sites = self._sites(sites)
        matrix = self._matrix(matrix, sites)
        shape = self.dimensions+state.shape[1:]
        tensor = state.reshape(shape, order="F")
        moved = np.moveaxis(tensor, sites, tuple(range(len(sites))))
        block = moved.reshape((matrix.shape[1], -1), order="F")
        transformed = (matrix@block).reshape(moved.shape, order="F")
        restored = np.moveaxis(transformed, tuple(range(len(sites))), sites)
        return restored.reshape(state.shape, order="F")

    def basis_state(self, occupations):
        if not isinstance(occupations, (tuple, list)) or len(occupations) != self.volume:
            raise ValueError("one occupation is required for each of the eight sites")
        index = 0
        for n, d, stride in zip(occupations, self.dimensions, self.strides):
            if (isinstance(n, (bool, np.bool_)) or not isinstance(n, Integral)
                    or not 0 <= n < d):
                raise ValueError("occupations must be integers within the declared site cutoffs")
            index += int(n)*stride
        result = np.zeros(self.dimension, dtype=complex)
        result[index] = 1
        return result

    def vacuum_state(self):
        """Product oscillator reference blank, not the interacting vacuum."""
        return self.basis_state((0,)*self.volume)

    def free_frequencies(self):
        """Normal frequencies of the untruncated free open-cell action."""
        return np.sqrt(np.linalg.eigvalsh(self.stiffness))

    def evolve(self, state, model_time):
        """Independent sparse evolution in the declared model-time coordinate."""
        state = self._state(state)
        model_time = _real(model_time, "model time")
        if model_time == 0:
            return state.copy()
        scale = -1j*model_time
        return expm_multiply(scale*self.hamiltonian, state,
                             traceA=scale*self.hamiltonian.diagonal().sum())
