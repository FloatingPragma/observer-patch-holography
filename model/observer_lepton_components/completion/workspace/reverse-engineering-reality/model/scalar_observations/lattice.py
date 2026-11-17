"""Seven-mode Galerkin scalar regulator on the fixed torus [0, 2*pi)^3.

The real orthonormal mode functions are 1/sqrt(V) and
sqrt(2/V) cos(x_j), sqrt(2/V) sin(x_j), with V=(2*pi)^3.  The momentum
band stays fixed when the spatial quadrature grid changes.  Grids with at
least five vertices per direction integrate the quartic exactly.  Only the
finite-difference gradient coefficient changes under this refinement.

Occupation matrices are compressions of the full oscillator polynomials.
In particular q2 and q4 are not powers of the already truncated q matrix.
The finite Ramsey read is instead a declared function of that finite q;
its relation to the unrestricted instrument needs a separate estimate.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import isfinite, pi, prod, sin, sqrt
from numbers import Integral, Real

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from ..relativistic_scalar.lattice import local_operators


MODE_NAMES = ("zero", "cos_x", "sin_x", "cos_y", "sin_y", "cos_z", "sin_z")
PAIRS = tuple(combinations(range(7), 2))
PHYSICAL_VOLUME = (2*pi)**3


def _real(value, name, *, positive=False, nonnegative=False):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Real)
            or not isfinite(value) or (positive and value <= 0)
            or (nonnegative and value < 0)):
        qualification = "positive " if positive else "nonnegative " if nonnegative else ""
        raise ValueError(f"{name} must be a finite {qualification}real")
    return float(value)


def _grid(value):
    if value is None:
        return None
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value < 5:
        raise ValueError("spatial grid must be None or an integer of at least five")
    return int(value)


def _mode(value):
    if isinstance(value, Integral) and not isinstance(value, (bool, np.bool_)) and 0 <= value < 7:
        return int(value)
    if isinstance(value, str) and value in MODE_NAMES:
        return MODE_NAMES.index(value)
    raise ValueError("mode must be a seven-mode index or a declared mode name")


def _dimension(value):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value < 1:
        raise ValueError("occupation dimensions must be positive integers")
    return int(value)


def _numeric(value, name):
    if np.ma.isMaskedArray(value):
        raise ValueError(f"{name} cannot be masked")
    result = np.asarray(value)
    if result.dtype.kind not in "fiuc" or not np.isfinite(result).all():
        raise ValueError(f"{name} must contain finite numeric entries")
    return result


def stiffness(grid=None):
    """First-harmonic finite-difference stiffness; None uses the exact gradient.

    For spacing a=2*pi/n, this is 4*sin(a/2)^2/a^2.  Both sine and cosine
    modes on all three axes have the same coefficient.  No mass or coupling
    is fitted when the grid changes.
    """
    grid = _grid(grid)
    if grid is None:
        return 1.0
    # This equivalent form avoids overflow in n**2 on very fine grids.
    phase = pi/grid
    return (sin(phase)/phase)**2


@dataclass(frozen=True)
class GalerkinParameters:
    """Supplied action coefficients with unit hbar and speed c=1.

    ``coupling`` is g=lambda/V, not the dimensionless physical lambda.
    The mass is the positive bare quadratic coefficient, not a measured
    interacting pole mass.  ``omega`` selects only the oscillator code.
    """

    grid: int | None = None
    mass: float = 1.0
    coupling: float = 1.0
    omega: float = 2.0

    def __post_init__(self):
        object.__setattr__(self, "grid", _grid(self.grid))
        for name in ("mass", "omega"):
            object.__setattr__(self, name, _real(getattr(self, name), name, positive=True))
        object.__setattr__(self, "coupling", _real(self.coupling, "coupling", nonnegative=True))
        try:
            scales = (self.mass**2, self.omega**2, self.physical_coupling)
        except OverflowError as exc:
            raise ValueError("derived action coefficients must be finite") from exc
        if not all(isfinite(value) for value in scales):
            raise ValueError("derived action coefficients must be finite")

    @property
    def kappa(self):
        return stiffness(self.grid)

    @property
    def volume(self):
        return PHYSICAL_VOLUME

    @property
    def physical_coupling(self):
        return self.coupling*PHYSICAL_VOLUME

    def frequency_squared(self, mode):
        return self.mass**2+(self.kappa if _mode(mode) else 0.)

    def onsite_quartic(self, mode):
        return self.coupling/(16 if _mode(mode) else 24)

    def pair_coefficient(self, left, right):
        left, right = _mode(left), _mode(right)
        if left == right:
            raise ValueError("an interaction pair requires two distinct modes")
        same_axis = left > 0 and right > 0 and (left-1)//2 == (right-1)//2
        return self.coupling/(8 if same_axis else 4)


def mode_functions(points):
    """Evaluate the seven real, physically normalized torus mode functions."""
    points = _numeric(points, "torus points")
    if points.ndim < 1 or points.shape[-1] != 3 or np.iscomplexobj(points):
        raise ValueError("torus points need final dimension three and real coordinates")
    values = [np.ones(points.shape[:-1])]
    for axis in range(3):
        values.extend((sqrt(2)*np.cos(points[..., axis]), sqrt(2)*np.sin(points[..., axis])))
    return np.stack(values, axis=-1)/sqrt(PHYSICAL_VOLUME)


def quartic_polynomial(coordinates):
    """Return W=V*integral(phi^4 dx) in the seven real mode coordinates."""
    coordinates = _numeric(coordinates, "mode coordinates")
    if coordinates.ndim < 1 or coordinates.shape[-1] != 7 or np.iscomplexobj(coordinates):
        raise ValueError("seven real mode coordinates are required")
    zero = coordinates[..., 0]
    radii = [coordinates[..., 2*axis+1]**2+coordinates[..., 2*axis+2]**2
             for axis in range(3)]
    return (zero**4+6*zero**2*sum(radii)+1.5*sum(radius**2 for radius in radii)
            +6*sum(radii[left]*radii[right] for left, right in combinations(range(3), 2)))


class SevenModeLattice:
    """Sparse P H P compression on a common seven-mode tensor-product code.

    The first mode is the least significant tensor index.  All 21 pairs
    have q2*q2 interactions.  These pair generators commute because every
    appearance of a given mode uses the same compressed full q2 matrix.
    """

    def __init__(self, dimensions=(2,)*7, parameters=None):
        if not isinstance(dimensions, (tuple, list)) or len(dimensions) != 7:
            raise ValueError("exactly seven occupation dimensions are required")
        self.dimensions = tuple(_dimension(value) for value in dimensions)
        self.dimension = prod(self.dimensions)
        self.parameters = GalerkinParameters() if parameters is None else parameters
        if not isinstance(self.parameters, GalerkinParameters):
            raise ValueError("parameters must be a GalerkinParameters instance")
        self.modes, self.pairs = MODE_NAMES, PAIRS
        self.strides = tuple(prod(self.dimensions[:index]) for index in range(7))
        self.operators = tuple(local_operators(dimension, self.parameters.omega)
                               for dimension in self.dimensions)
        self.hamiltonian = self._hamiltonian()

    def onsite(self, mode):
        mode = _mode(mode)
        op = self.operators[mode]
        return (op["p2"]/2+self.parameters.frequency_squared(mode)*op["q2"]/2
                +self.parameters.onsite_quartic(mode)*op["q4"])

    def pair_coefficient(self, left, right):
        return self.parameters.pair_coefficient(left, right)

    def _tensor_product(self, operators):
        result = sparse.csr_matrix([[1.]])
        for mode in reversed(range(7)):
            factor = operators.get(mode)
            if factor is None:
                factor = sparse.eye(self.dimensions[mode], format="csr")
            result = sparse.kron(result, sparse.csr_matrix(factor), format="csr")
        result.eliminate_zeros()
        return result

    def _hamiltonian(self):
        result = sparse.csr_matrix((self.dimension, self.dimension), dtype=float)
        for mode in range(7):
            result += self._tensor_product({mode: self.onsite(mode)})
        for left, right in self.pairs:
            coefficient = self.pair_coefficient(left, right)
            if coefficient:
                result += coefficient*self._tensor_product(
                    {left: self.operators[left]["q2"], right: self.operators[right]["q2"]})
        result.eliminate_zeros()
        result.sort_indices()
        if not np.isfinite(result.data).all():
            raise ValueError("compressed action entries must be finite")
        return result

    def _sites(self, sites):
        if isinstance(sites, (Integral, str)) and not isinstance(sites, (bool, np.bool_)):
            sites = (sites,)
        if not isinstance(sites, (tuple, list)) or not sites:
            raise ValueError("a nonempty sequence of distinct modes is required")
        selected = tuple(_mode(site) for site in sites)
        if len(set(selected)) != len(selected):
            raise ValueError("local mode indices must be distinct")
        return selected

    def _state(self, state):
        state = _numeric(state, "state")
        if state.ndim not in (1, 2) or state.shape[0] != self.dimension:
            raise ValueError("state must have one amplitude row per tensor occupation")
        return np.asarray(state, dtype=complex)

    def _matrix(self, matrix, sites):
        if sparse.issparse(matrix):
            if matrix.dtype.kind not in "fiuc" or not np.isfinite(matrix.data).all():
                raise ValueError("local matrix must contain finite numeric entries")
            matrix = matrix.tocsr()
        else:
            matrix = _numeric(matrix, "local matrix")
        size = prod(self.dimensions[site] for site in sites)
        if matrix.shape != (size, size):
            raise ValueError("local matrix shape must match the selected occupation dimensions")
        return matrix

    def embed_local(self, matrix, sites):
        """Embed with the first listed mode least significant in the local block."""
        sites = self._sites(sites)
        matrix = sparse.coo_matrix(self._matrix(matrix, sites))
        if len(sites) == 1:
            return self._tensor_product({sites[0]: matrix})
        spectators = tuple(mode for mode in range(7) if mode not in sites)
        offsets = np.zeros(prod(self.dimensions[mode] for mode in spectators), dtype=np.int64)
        remainder = np.arange(len(offsets), dtype=np.int64)
        for mode in spectators:
            offsets += (remainder % self.dimensions[mode])*self.strides[mode]
            remainder //= self.dimensions[mode]

        def expand(indices):
            remainder = indices.copy()
            result = np.zeros(len(indices), dtype=np.int64)
            for mode in sites:
                result += (remainder % self.dimensions[mode])*self.strides[mode]
                remainder //= self.dimensions[mode]
            return result

        rows = (expand(matrix.row)[:, None]+offsets).ravel()
        cols = (expand(matrix.col)[:, None]+offsets).ravel()
        values = np.repeat(matrix.data, len(offsets))
        result = sparse.coo_matrix((values, (rows, cols)),
                                  shape=(self.dimension, self.dimension)).tocsr()
        result.eliminate_zeros()
        return result

    def apply_local(self, state, sites, matrix):
        """Apply the local matrix to amplitudes or columns without normalization."""
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
        if not isinstance(occupations, (tuple, list)) or len(occupations) != 7:
            raise ValueError("one occupation is required for each of the seven modes")
        index = 0
        for value, dimension, stride in zip(occupations, self.dimensions, self.strides):
            if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral)
                    or not 0 <= value < dimension):
                raise ValueError("occupations must be integers within the declared cutoffs")
            index += int(value)*stride
        state = np.zeros(self.dimension, dtype=complex)
        state[index] = 1
        return state

    def vacuum_state(self):
        """Return the supplied product basis blank, not the interacting vacuum."""
        return self.basis_state((0,)*7)

    def evolve(self, state, model_time):
        """Independent sparse Hamiltonian evolution in supplied model time."""
        state = self._state(state)
        model_time = _real(model_time, "model time")
        if model_time == 0:
            return state.copy()
        scale = -1j*model_time
        return expm_multiply(scale*self.hamiltonian, state,
                             traceA=scale*self.hamiltonian.diagonal().sum())

    def metadata(self):
        return {"physical_torus_side": 2*pi, "physical_volume": PHYSICAL_VOLUME,
                "modes": list(self.modes), "spatial_grid_per_axis": self.parameters.grid,
                "spatial_spacing": None if self.parameters.grid is None else 2*pi/self.parameters.grid,
                "gradient_stiffness": self.parameters.kappa,
                "bare_mass": self.parameters.mass, "coupling_lambda_over_volume": self.parameters.coupling,
                "physical_lambda": self.parameters.physical_coupling,
                "basis_frequency": self.parameters.omega, "local_dimensions": list(self.dimensions),
                "state_dimension": self.dimension, "momentum_band_changes_with_grid": False,
                "interacting_pole_mass_calibrated": False,
                "spatial_ultraviolet_cutoff_removed": False}
