"""A full finite source specified by sparse CAR words and scalar columns.

The exact scalar energy is the unique positive square root of the integer
stiffness matrix. Its displayed large-volume columns use a numerical square
root with a retained residual; the one-cell identity is exact.
"""
from __future__ import annotations

from collections import defaultdict
from functools import lru_cache
from itertools import product

import numpy as np
from scipy import linalg as la
import sympy as s

from ..compact_defect_pair.topology import direction
from ..coherent_defect_wave.chiral_joint.higgs_composites import model as car
from ..coherent_defect_wave.chiral_joint.model import sequence
from ..coherent_defect_wave.chiral_higgs.model import higgs_memory, PAULI
from ..coherent_defect_wave.physical_attachment.chiral import representation

CHANNELS = ('up', 'down', 'charged_lepton')
LEPTONS = (9, 10, 11, 24, 25, 26)


def sparse(matrix):
    return s.ImmutableSparseMatrix(matrix)


def clean(poly):
    return {key: s.expand(value) for key, value in poly.items() if s.expand(value) != 0}


@lru_cache(maxsize=1)
def local_interactions():
    """All three actual source polynomials at primitive coefficient one."""
    parts, coefficients = car.quantum_yukawa()
    return {name: car.clean({key: value.subs({y: 1 for y in coefficients.values()})
                            for key, value in parts[name].items()}) for name in CHANNELS}


def compact_points(radius=3):
    if type(radius) is not int or radius < 3:
        raise ValueError('the compact source needs its complete exterior layer')
    return tuple(product(range(-radius, radius+1), repeat=3))


def spatial_operators(points, mass=64, texture=True):
    """Literal compact hopping and its scalar stiffness on the same graph."""
    points = tuple(tuple(x) for x in points)
    if (not points or len(set(points)) != len(points)
            or any(len(p) != 3 or any(type(x) is not int for x in p) for p in points)
            or type(mass) is not int or mass < 8):
        raise ValueError('distinct integer spatial sites and integer texture amplitude >=8 required')
    n, address = len(points), {p: i for i, p in enumerate(points)}
    h = s.MutableSparseMatrix(2*n, 2*n, {})
    stiffness = s.MutableSparseMatrix.eye(n)
    edges, occupied_texture = [], []
    for x, p in enumerate(points):
        raw = direction(p) if texture else (0, 0, 0)
        length = s.sqrt(sum(v*v for v in raw))
        active = int(length != 0)
        occupied_texture.append(active)
        stiffness[x, x] += mass*mass*active
        if active:
            block = sum((mass*raw[a]*PAULI[a]/length for a in range(3)), s.zeros(2))
            for (i, j), value in sparse(block).todok().items():
                h[2*x+i, 2*x+j] = value
        for axis in range(3):
            following = tuple(v+int(a == axis) for a, v in enumerate(p))
            if following not in address:
                continue
            y = address[following]
            edges.append((x, y, axis))
            stiffness[x, x] += 1
            stiffness[y, y] += 1
            stiffness[x, y] = stiffness[y, x] = -1
            block = s.I*PAULI[axis]/2
            for (i, j), value in sparse(block).todok().items():
                h[2*x+i, 2*y+j] = value
                h[2*y+j, 2*x+i] = s.conjugate(value)
    if h != h.H or stiffness != stiffness.T:
        raise ArithmeticError('the actual spatial coefficients are not Hermitian')
    return {'points': points, 'h': sparse(h), 'K': sparse(stiffness),
            'edges': tuple(edges), 'texture_active': tuple(occupied_texture)}


class Source:
    """Complete finite Hamiltonian, with all occupation and scalar states.

    mode=30*site+15*spin+internal. Scalar labels are0, then all particle
    doublets, then all antiparticle doublets. The only scalar cutoff is
    total scalar occupation<=1 across this declared finite member.
    """

    def __init__(self, points, *, mass=64, texture=True, channels=CHANNELS):
        channels = tuple(channels)
        if len(set(channels)) != len(channels) or not set(channels) <= set(CHANNELS):
            raise ValueError('an explicit subset of the three unit native channels required')
        self.spatial = spatial_operators(points, mass, texture)
        self.points = self.spatial['points']
        self.sites = len(self.points)
        self.modes, self.scalar_dimension = 30*self.sites, 1+4*self.sites
        self.channels, self.mass, self.texture = channels, mass, texture
        self.stagger = tuple(1 if sum(p) % 2 == 0 else -1 for p in self.points)
        self.h = np.asarray(self.spatial['h'], complex)
        self.K = np.asarray(self.spatial['K'], float)
        eigenvalues, vectors = la.eigh(self.K)
        if eigenvalues[0] < 1-1e-10:
            raise ArithmeticError('the positive stiffness bound failed')
        self.frequencies = np.sqrt(eigenvalues)
        self.scalar_modes = vectors
        self.Omega = (vectors*self.frequencies)@vectors.T
        self.scalar_root_residual = float(la.norm(self.Omega@self.Omega-self.K, 2))
        self.h_columns = defaultdict(list)
        self.exact_h_columns = defaultdict(list)
        for (row, col), value in self.spatial['h'].todok().items():
            self.h_columns[col].append((row, complex(value)))
            self.exact_h_columns[col].append((row, value))
        self.local_words = car.add(*(local_interactions()[name] for name in channels))
        self.pair_columns = defaultdict(list)
        self.exact_pair_columns = defaultdict(list)
        for (word, row, col), value in self.local_words.items():
            self.pair_columns[col].append((word, row, complex(value)))
            self.exact_pair_columns[col].append((word, row, value))

    def scalar_label(self, site, local):
        if local == 0:
            return 0
        if not 0 <= site < self.sites or local not in (1, 2, 3, 4):
            raise ValueError('valid site and local scalar state required')
        return 1+2*site+(local-1)%2+(2*self.sites if local > 2 else 0)

    def decode_scalar(self, label):
        if label == 0:
            return None, 0
        if type(label) is not int or not 0 < label < self.scalar_dimension:
            raise ValueError('full scalar coordinate required')
        orientation, remainder = divmod(label-1, 2*self.sites)
        site, component = divmod(remainder, 2)
        return site, 1+component+2*orientation

    def column_terms(self, mask, scalar, *, exact_omega=None):
        """Every literal term of a full column, before coincident rows combine.

        An optional exact scalar square root supplies algebraic coefficients.
        Its caller must verify its positivity and square-root identity once.
        """
        if type(mask) is not int or not 0 <= mask < 1 << self.modes:
            raise ValueError('a full source CAR occupation mask required')
        scalar_site, local_scalar = self.decode_scalar(scalar)
        h_columns = self.h_columns if exact_omega is None else self.exact_h_columns
        pair_columns = self.pair_columns if exact_omega is None else self.exact_pair_columns
        omega = self.Omega if exact_omega is None else exact_omega
        occupied = mask
        while occupied:
            bit = occupied & -occupied
            old = bit.bit_length()-1
            spatial, internal = divmod(old, 15)
            for following, coefficient in h_columns[spatial]:
                step = sequence(mask, ((15*following+internal, True), (old, False)))
                if step:
                    yield (step[0], scalar), coefficient*int(step[1])
            occupied ^= bit
        if scalar:
            for following in range(self.sites):
                coefficient = omega[following, scalar_site]
                if coefficient:
                    yield (mask, self.scalar_label(following, local_scalar)), coefficient
        sites = range(self.sites) if scalar == 0 else (scalar_site,)
        for site in sites:
            for word, following, coefficient in pair_columns[local_scalar]:
                shifted = tuple((mode+30*site, creation) for mode, creation in word)
                step = sequence(mask, shifted)
                if step:
                    yield ((step[0], self.scalar_label(site, following)),
                           self.stagger[site]*coefficient*int(step[1]))

    def column(self, mask, scalar):
        """One complete numerical column, with no particle-number truncation."""
        out = defaultdict(complex)
        for target, coefficient in self.column_terms(mask, scalar):
            out[target] += coefficient
        return {key: value for key, value in out.items() if value != 0}

    def apply(self, state):
        out = defaultdict(complex)
        for key, amplitude in state.items():
            for target, coefficient in self.column(*key).items():
                out[target] += amplitude*coefficient
        return {key: value for key, value in out.items() if value != 0}

    def norm_bound(self):
        """A conservative bound on the exact entire finite Hamiltonian."""
        degree = max(sum(x in (a, b) for a, b, _ in self.spatial['edges']) for x in range(self.sites))
        kinetic = self.modes*((self.mass if self.texture else 0)+s.Rational(degree, 2))
        scalar = s.sqrt(1+self.mass**2*max(self.spatial['texture_active'])+2*degree)
        interaction = self.sites*sum(abs(value) for value in self.local_words.values())
        return s.simplify(kinetic+scalar+interaction)

    def descriptor(self):
        return {'points': [list(p) for p in self.points], 'sites': self.sites,
            'CAR_modes': self.modes, 'scalar_dimension': self.scalar_dimension,
            'complete_Hilbert_dimension': f'2^{self.modes}*{self.scalar_dimension}',
            'channels_at_unit_primitive_coefficient': list(self.channels),
            'texture_amplitude': self.mass, 'texture_enabled': self.texture,
            'scalar_stiffness': 'K=I+L_graph+M^2 diag(|n_x|^2)',
            'scalar_energy': 'Omega is the unique positive square root of K',
            'Yukawa_site_phase': '(-1)^(x1+x2+x3)',
            'Hamiltonian': 'dGamma(h tensor I15)+diag(0,Omega tensor I2,conjugate(Omega) tensor I2)+sum_x D_x H_Y,x',
            'global_current_symmetry': 'the same twelve native R15/Higgs currents',
            'local_quantum_gauge_links': False,
            'norm_bound': str(self.norm_bound()),
            'numeric_scalar_root_residual': self.scalar_root_residual,
            'numeric_columns_are_outward_certified': False,
            'full_occupation_carrier_is_materialized': False}


@lru_cache(maxsize=2)
def compact_source():
    return Source(compact_points())


@lru_cache(maxsize=2)
def one_cell(channels=CHANNELS):
    return Source(((0, 0, 0),), texture=False, channels=channels)


def one_cell_polynomial(channels=CHANNELS):
    """The exact normal-ordered member, independent of numerical columns."""
    return car.add(*(local_interactions()[name] for name in channels),
                   car.monomial(matrix=higgs_memory()[2]))


def lepton_mask(mask):
    return sum(((mask >> i) & 1) << m for i, m in enumerate(LEPTONS))


def map_lepton_column_back(column):
    out = {}
    for (mask, scalar), value in column.items():
        if mask & ~sum(1 << i for i in LEPTONS):
            raise ArithmeticError('a quark occupation prevents the lepton-only inclusion')
        reduced = sum(((mask >> m) & 1) << i for i, m in enumerate(LEPTONS))
        out[5*reduced+scalar] = value
    return out
