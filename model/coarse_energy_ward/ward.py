"""Exact rational CAR currents, nearest-face routing and block aggregation.

Clifford monomials form a faithful basis of the complete finite Fock algebra;
no number sector or state is selected.  Floating matrices are independent
small-volume diagnostics only.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from itertools import permutations
from math import prod

from ..full_fock_field.field import FieldGeometry
from ..full_fock_field.local_simulation import local_tail_boundary_bound, split_bound


def _add(a, b):
    return a[0]+b[0], a[1]+b[1]


def _mul(a, b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def _coefficient(value):
    if type(value) in (int, F):
        return F(value), F()
    if type(value) is tuple and len(value) == 2 and all(type(x) in (int, F) for x in value):
        return F(value[0]), F(value[1])
    raise ValueError('exact Gaussian rational coefficient required')


def multiply_masks(left, right):
    """gamma_i^2=1; gamma_i gamma_j=-gamma_j gamma_i for i != j."""
    inversions = 0
    bits = left
    while bits:
        bit = bits & -bits
        inversions += (right & (bit-1)).bit_count()
        bits ^= bit
    return left ^ right, -1 if inversions % 2 else 1


@dataclass(frozen=True)
class Clifford:
    terms: dict

    def __post_init__(self):
        if any(type(mask) is not int or mask < 0 for mask in self.terms):
            raise ValueError('nonnegative integer Clifford masks required')
        object.__setattr__(self, 'terms', {m: _coefficient(c) for m, c in self.terms.items()
                                         if _coefficient(c) != (0, 0)})

    @classmethod
    def scalar(cls, value):
        return cls({0: _coefficient(value)})

    def __add__(self, other):
        terms = self.terms.copy()
        for mask, value in other.terms.items():
            terms[mask] = _add(terms.get(mask, (F(), F())), value)
        return Clifford(terms)

    def scale(self, value):
        value = _coefficient(value)
        return Clifford({m: _mul(c, value) for m, c in self.terms.items()})

    def __neg__(self):
        return self.scale(-1)

    def __sub__(self, other):
        return self+-other

    def __matmul__(self, other):
        terms = {}
        for left, a in self.terms.items():
            for right, b in other.terms.items():
                mask, sign = multiply_masks(left, right)
                value = _mul(a, b)
                terms[mask] = _add(terms.get(mask, (F(), F())), (sign*value[0], sign*value[1]))
        return Clifford(terms)

    def adjoint(self):
        return Clifford({m: ((-1)**(m.bit_count()*(m.bit_count()-1)//2)*c[0],
                             -(-1)**(m.bit_count()*(m.bit_count()-1)//2)*c[1])
                         for m, c in self.terms.items()})

    def is_zero(self):
        return not self.terms

    def norm_bound(self):
        """Rational triangle bound, since each Majorana monomial is unitary."""
        return sum((abs(a)+abs(b) for a, b in self.terms.values()), F())

    def support(self):
        mask = 0
        for term in self.terms:
            mask |= term
        return frozenset(i//4 for i in range(mask.bit_length()) if mask >> i & 1)

    def matrix(self, modes):
        """Independent JW representation; only small full-Fock diagnostics."""
        import numpy as np
        from scipy import sparse
        from ..full_fock_field.field import FullFock
        if type(modes) is not int or not 2 <= modes <= 8:
            raise ValueError('matrix diagnostic limited to two through eight modes')
        if any(m.bit_length() > 2*modes for m in self.terms):
            raise ValueError('Clifford mask exceeds the stated carrier')
        space = FullFock(modes)
        identity = sparse.eye(space.dimension, dtype=complex, format='csr')
        gammas = []
        for p in range(modes):
            a = space.annihilator(p)
            gammas.extend((a+a.conj().T, -1j*(a-a.conj().T)))
        result = sparse.csr_matrix(identity.shape, dtype=complex)
        for mask, (real, imag) in self.terms.items():
            monomial = identity
            for i, gamma in enumerate(gammas):
                if mask >> i & 1:
                    monomial = monomial@gamma
            result += complex(float(real), float(imag))*monomial
        return np.asarray(result.toarray())


def total(values):
    result = Clifford({})
    for value in values:
        result = result+value
    return result


def annihilator(mode):
    return Clifford({1 << (2*mode): (F(1, 2), 0),
                     1 << (2*mode+1): (0, F(1, 2))})


def commutator(a, b):
    return a@b-b@a


def energies(geometry):
    if not isinstance(geometry, FieldGeometry) or geometry.volume > 64:
        raise ValueError('exact symbolic certificate admits one through 64 sites')
    # The source exposes binary floating coefficients. Interpret those exact
    # values, never a rounded decimal approximation to another action.
    c, u = F(geometry.speed), F(geometry.interaction)
    if max(c.denominator, u.denominator) > 2**16:
        raise ValueError('binary coefficient denominator exceeds certificate domain')
    ladders = [annihilator(p) for p in range(geometry.modes)]
    local = [Clifford({}) for _ in geometry.sites]
    for group in geometry.groups:
        for p, q, z in group:
            term = (ladders[p].adjoint()@ladders[q]).scale((F(z.real), F(z.imag)))
            half = (term+term.adjoint()).scale(F(1, 2))
            local[p//2] = local[p//2]+half
            local[q//2] = local[q//2]+half
    for x in range(geometry.volume):
        n0 = ladders[2*x].adjoint()@ladders[2*x]
        n1 = ladders[2*x+1].adjoint()@ladders[2*x+1]
        local[x] = local[x]+(n0@n1).scale(u)
    return tuple(local)


def shortest_paths(start, end):
    """Every shortest rectangular-lattice path, without an axis-priority choice."""
    if (len(start) != 3 or len(end) != 3 or
            any(type(x) is not int for x in (*start, *end))):
        raise ValueError('three integer coordinates required')
    steps = []
    for axis in range(3):
        steps.extend([(axis, 1 if end[axis] > start[axis] else -1)]*abs(end[axis]-start[axis]))
    if not 1 <= len(steps) <= 2:
        raise ValueError('routing certificate requires separation one or two')
    paths = []
    for order in sorted(set(permutations(steps))):
        point, path = start, [start]
        for axis, sign in order:
            next_point = list(point)
            next_point[axis] += sign
            point = tuple(next_point)
            path.append(point)
        paths.append(tuple(path))
    return tuple(paths)


def route_pair(sites, x, y):
    indices = {site: i for i, site in enumerate(sites)}
    paths = shortest_paths(sites[x], sites[y])
    weights = {}
    for path in paths:
        for left, right in zip(path, path[1:]):
            a, b = indices[left], indices[right]
            edge, sign = ((a, b), 1) if a < b else ((b, a), -1)
            weights[edge] = weights.get(edge, F())+F(sign, len(paths))
    return weights


def divergence(flows, vertices):
    result = [Clifford({}) for _ in range(vertices)]
    for (x, y), current in flows.items():
        result[x] = result[x]+current
        result[y] = result[y]-current
    return tuple(result)


def aggregate(values, labels):
    if len(values) != len(labels) or any(type(a) is not int or a < 0 for a in labels):
        raise ValueError('one nonnegative integer block label per value required')
    if not labels:
        raise ValueError('nonempty finite block carrier required')
    result = [Clifford({}) for _ in range(max(labels)+1)]
    for label, value in zip(labels, values, strict=True):
        result[label] = result[label]+value
    return tuple(result)


def aggregate_flows(flows, labels):
    if not labels or any(type(a) is not int or a < 0 for a in labels):
        raise ValueError('nonempty nonnegative block labels required')
    result = {}
    for (x, y), current in flows.items():
        if not 0 <= x < y < len(labels):
            raise ValueError('canonical oriented edge outside the carrier')
        a, b = labels[x], labels[y]
        if a == b:
            continue
        edge, sign = ((a, b), 1) if a < b else ((b, a), -1)
        result[edge] = result.get(edge, Clifford({}))+current.scale(sign)
    return {edge: value for edge, value in result.items() if not value.is_zero()}


def block_labels(sites, widths):
    if len(widths) != 3 or any(type(w) is not int or w < 1 for w in widths):
        raise ValueError('three positive integer block widths required')
    blocks = [tuple(x//w for x, w in zip(site, widths, strict=True)) for site in sites]
    index = {block: i for i, block in enumerate(sorted(set(blocks)))}
    return tuple(index[block] for block in blocks)


def certify(geometry=FieldGeometry((3, 3, 3))):
    local = energies(geometry)
    h = total(local)
    currents, routed = {}, {}
    for x in range(geometry.volume):
        for y in range(x+1, geometry.volume):
            current = commutator(local[x], local[y]).scale((0, 1))
            if current.is_zero():
                continue
            if current != current.adjoint():
                raise ArithmeticError('energy current lost Hermiticity')
            currents[x, y] = current
            for edge, weight in route_pair(geometry.sites, x, y).items():
                routed[edge] = routed.get(edge, Clifford({}))+current.scale(weight)
    original, nearest = divergence(currents, geometry.volume), divergence(routed, geometry.volume)
    for x, hx in enumerate(local):
        if hx != hx.adjoint() or original[x] != nearest[x]:
            raise ArithmeticError('exact local energy or routing identity failed')
        if not (commutator(h, hx).scale((0, 1))+nearest[x]).is_zero():
            raise ArithmeticError('same-action nearest-face Ward identity failed')
    # All rectangular widths, including ragged edge blocks; no state required.
    from itertools import product
    partitions = 0
    for widths in product(*(range(1, n+1) for n in geometry.sides)):
        labels = block_labels(geometry.sites, widths)
        lhs = divergence(aggregate_flows(routed, labels), max(labels)+1)
        if lhs != aggregate(nearest, labels):
            raise ArithmeticError('block divergence failed to commute')
        partitions += 1
    return {'local': local, 'H': h, 'currents': currents, 'routed': routed,
            'rectangular_partitions': partitions}


def readout_error(radius, duration, steps, speed=F(1, 8), interaction=F(1),
                  *, samples=None, endpoint_readout=F(), face_readout=F()):
    """Uniform operator error for two energy endpoints and integrated face flux.

    Source Strang bound and local tail are inherited from full_fock_field.
    Arbitrary input states/references are allowed, but their patch restriction
    and the observable readout service are inputs, not prepared here.
    """
    if (type(radius) is not int or not 3 < radius <= 4096 or
            type(duration) not in (int, F) or duration < 0 or
            type(steps) is not int or steps < 1 or
            any(type(x) not in (int, F) or x < 0 for x in (speed, interaction))):
        raise ValueError('radius > 3, nonnegative rational data and positive steps required')
    if (samples is not None and (type(samples) is not int or samples < 1)):
        raise ValueError('positive integer midpoint sample count required')
    if any(type(x) not in (int, F) or x < 0 for x in (endpoint_readout, face_readout)):
        raise ValueError('exact nonnegative observable readout error bounds required')
    duration, speed, interaction = F(duration), F(speed), F(interaction)
    a = interaction+3*speed
    j = 14*a*a
    sides = (2*radius+1,)*3
    split = min(F(2), split_bound(sides, duration, steps, speed, interaction))
    energy_tail = local_tail_boundary_bound(radius, 1, 7, duration, speed)
    flux_tail = local_tail_boundary_bound(radius, 3, 343, duration, speed)
    energy = a*min(F(2), energy_tail+2*split)
    integrated_flux = duration*j*min(F(2), flux_tail+2*split)
    # The face support lies in a radius-three cube, with at most 343 sites.
    # ||dJ/dt|| <= 2 ||J|| * 343*(U+6c); midpoint integration uses only
    # this Lipschitz bound, not numerical smoothness or a sampled fit.
    quadrature = (duration**2*j*343*(interaction+6*speed)/(2*samples)
                  if samples is not None else F())
    energy += F(endpoint_readout)
    integrated_flux += quadrature+duration*F(face_readout)
    # One four-dimensional cell has two energy endpoints and six spatial faces.
    return {'local_energy_norm': a, 'routed_face_norm': j,
            'source_unitary_error': split, 'midpoint_quadrature_error': quadrature,
            'energy_endpoint_error': energy,
            'integrated_face_error': integrated_flux,
            'cell_Ward_residual_bound': 2*energy+6*integrated_flux}
