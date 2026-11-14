"""Exact Gaussian-integer energy currents for the retained Weyl/contact field.

Energy entries have denominator 32; current entries have denominator 1024.
Integer sparse arithmetic, with checked bounds before products, certifies the
operator identities. Hamiltonian model time is the only time parameter here.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy import sparse
from ..interacting_field.lattice import WeylLattice

SCALE = 32


def _clean(value):
    value = sparse.csr_matrix(value, dtype=np.int64)
    value.eliminate_zeros()
    return value


@dataclass(frozen=True)
class GaussianInteger:
    real: sparse.csr_matrix
    imag: sparse.csr_matrix

    @classmethod
    def zero(cls, dimension):
        z = sparse.csr_matrix((dimension, dimension), dtype=np.int64)
        return cls(z, z.copy())

    @classmethod
    def from_exact(cls, matrix):
        matrix = sparse.csr_matrix(matrix, dtype=complex)
        if not np.isfinite(matrix.data).all():
            raise ValueError('finite Gaussian-integer entries required')
        for part in (matrix.data.real, matrix.data.imag):
            if not np.array_equal(part, np.rint(part)) or np.max(abs(part), initial=0) > 2**30:
                raise ValueError('exact bounded Gaussian-integer entries required')
        return cls(_clean(matrix.real), _clean(matrix.imag))

    def __add__(self, other):
        if self.maximum()+other.maximum() >= 2**62:
            raise OverflowError('Gaussian-integer addition bound exceeded')
        return GaussianInteger(_clean(self.real+other.real), _clean(self.imag+other.imag))

    def __sub__(self, other):
        if self.maximum()+other.maximum() >= 2**62:
            raise OverflowError('Gaussian-integer subtraction bound exceeded')
        return GaussianInteger(_clean(self.real-other.real), _clean(self.imag-other.imag))

    def __matmul__(self, other):
        # Each output component is a difference/sum of two integer dot products.
        # This deliberately conservative bound covers every intermediate sum.
        size = self.real.shape[0]
        if 2*size*self.maximum()*other.maximum() >= 2**62:
            raise OverflowError('Gaussian-integer multiplication bound exceeded')
        return GaussianInteger(_clean(self.real@other.real-self.imag@other.imag),
                               _clean(self.real@other.imag+self.imag@other.real))

    def times_i(self):
        return GaussianInteger(_clean(-self.imag), self.real)

    def adjoint(self):
        return GaussianInteger(_clean(self.real.T), _clean(-self.imag.T))

    def maximum(self):
        return max(int(np.max(abs(self.real.data), initial=0)),
                   int(np.max(abs(self.imag.data), initial=0)))

    def is_zero(self):
        return self.real.nnz == self.imag.nnz == 0

    def matrix(self, denominator=1):
        return (self.real.astype(complex)+1j*self.imag)/denominator


def commutator(first, second):
    return first@second-second@first


def local_energies(lattice):
    """h_x=U D_x plus half of each incident hopping bond, on the exact N=2 sector."""
    if lattice.speed != 1/8 or lattice.interaction not in (0., 1.):
        raise ValueError('this integer certificate covers c=1/8 and U=0 or 1')
    single = lattice.single_h.tocoo()
    result = []
    for site in range(lattice.volume):
        keep = (single.row//2 == site) | (single.col//2 == site)
        one = sparse.coo_matrix((single.data[keep]/2,
                                (single.row[keep], single.col[keep])),
                               shape=single.shape)
        matrix = lattice.lift_one_body(one)
        idx = int(lattice.contact_indices[site])
        contact = sparse.csr_matrix(([lattice.interaction], ([idx], [idx])), shape=matrix.shape)
        result.append(GaussianInteger.from_exact(SCALE*(matrix+contact)))
    return result


def sum_exact(operators, dimension):
    total = GaussianInteger.zero(dimension)
    for operator in operators:
        total = total+operator
    return total


def certify(lattice=None):
    lattice = WeylLattice(1) if lattice is None else lattice
    local = local_energies(lattice)
    total = sum_exact(local, lattice.dimension)
    independent = GaussianInteger.from_exact(SCALE*lattice.hamiltonian)
    if not (total-independent).is_zero():
        raise ArithmeticError('local energies do not partition the independently assembled action')
    divergences = [GaussianInteger.zero(lattice.dimension) for _ in local]
    currents = {}
    maximum_range = 0
    for x in range(lattice.volume):
        if not (local[x]-local[x].adjoint()).is_zero():
            raise ArithmeticError('non-Hermitian local energy')
        for y in range(x+1, lattice.volume):
            current = commutator(local[x], local[y]).times_i()
            if not (current-current.adjoint()).is_zero():
                raise ArithmeticError('non-Hermitian energy current')
            if not (current+commutator(local[y], local[x]).times_i()).is_zero():
                raise ArithmeticError('energy-current antisymmetry failed')
            if current.is_zero():
                continue
            distance = sum(abs(a-b) for a,b in zip(lattice.sites[x], lattice.sites[y], strict=True))
            if distance > 2:
                raise ArithmeticError('current exceeds the support-overlap range')
            maximum_range = max(maximum_range, distance)
            currents[x,y] = current
            divergences[x] = divergences[x]+current
            divergences[y] = divergences[y]-current
    for energy, divergence in zip(local, divergences, strict=True):
        if not (commutator(total, energy).times_i()+divergence).is_zero():
            raise ArithmeticError('exact local energy continuity failed')
    if not sum_exact(divergences, lattice.dimension).is_zero():
        raise ArithmeticError('total energy flux failed to cancel')
    # A complete half-box region includes every crossing current, including
    # range-two currents generated by our half-bond energy assignment.
    region = {i for i,x in enumerate(lattice.sites) if x[0] <= 0}
    regional = sum_exact([local[x] for x in region], lattice.dimension)
    boundary = GaussianInteger.zero(lattice.dimension)
    crossing = 0
    for (x,y), current in currents.items():
        if (x in region) != (y in region):
            boundary = boundary+current if x in region else boundary-current
            crossing += 1
    if not (commutator(total, regional).times_i()+boundary).is_zero():
        raise ArithmeticError('regional energy balance failed')
    return {'sites': lattice.volume, 'dimension': lattice.dimension,
            'local_energy_denominator': SCALE, 'current_denominator': SCALE**2,
            'operator_identities_exact': True,
            'integer_arithmetic': 'signed int64 real and imaginary sparse matrices; products bounded below 2^62',
            'local_continuity_equations_checked': lattice.volume,
            'all_unordered_site_pairs_checked': lattice.volume*(lattice.volume-1)//2,
            'nonzero_unordered_currents': len(currents),
            'maximum_current_range_in_lattice_steps': maximum_range,
            'region_sites': len(region), 'boundary_current_pairs': crossing,
            'nonzero_regional_flux_operator': not boundary.is_zero(),
            'largest_current_numerator_component': max(c.maximum() for c in currents.values()),
            'time_parameter': 'Hamiltonian model time; hbar=1 is a supplied action convention'}, local, currents
