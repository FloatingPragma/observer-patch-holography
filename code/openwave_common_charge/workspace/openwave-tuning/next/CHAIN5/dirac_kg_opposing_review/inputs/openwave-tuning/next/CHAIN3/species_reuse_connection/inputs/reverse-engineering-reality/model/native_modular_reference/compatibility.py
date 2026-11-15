"""All faithful references compatible with every bare two-port algebra.

Let H=-log(sigma) on the binary total-two shell over N>=3 ports. A faithful sigma
preserves each number-preserving bare A_ij under its modular group iff

    H = sum_i a_i N_i, equivalently sigma_{ij}=lambda*r_i*r_j > 0.

This includes initially nondiagonal references. Indeed N_i belongs to every
A_ij with j!=i, so invariance implies [H,N_i] belongs to their intersection.
Two distinct partners have disjoint off-diagonal supports. The intersection
is diagonal, whereas every diagonal entry of [H,N_i] is zero. Thus H commutes
with every N_i, whose joint eigenspaces are the individual configurations.
For N=3 the intersection is the FULL three-dimensional diagonal algebra,
not span{I,N_i}; the same zero-diagonal argument still applies.

For diagonal H, the matrix unit sum |i,k><j,k| over outside spectators k
belongs to A_ij. Its commutator belongs to A_ij precisely when H_ik-H_jk
does not depend on k. Exact row reduction below verifies these linear
constraints have kernel H_ij=a_i+a_j. Their integer matrix has the same rank
over Q, R and C, so this verifies all Hermitian generators, not only rational
examples. The explicit left inverse of the
pair-incidence matrix proves this kernel has dimension N. Exponentiating
gives the product-weight criterion already used by common_caps.py.
Conversely additive H preserves every bare regional algebra: outside
energies cancel in each local matrix unit. This is an algebraic theorem,
not a construction of controls or a classification of nonfaithful states.

The finite verifier supports 3..16 ports. The proof above is uniform in
N>=3. Every native off-diagonal coordinate receives an explicit forbidden
commutator-entry witness; diagonal ranks use Fractions, never eigenvalues
or a numerical rank threshold. Exact complex entries are represented by
(rational real part, rational imaginary part), with Hermiticity checked.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import combinations
import json
from typing import NamedTuple


ZERO = (F(), F())
MAX_PORTS = 16


def _vertices(vertices):
    if type(vertices) is not int or not 3 <= vertices <= MAX_PORTS:
        raise ValueError('an original integer from 3 to 16 ports required')
    return vertices


@lru_cache(maxsize=14)
def _basis(vertices):
    return tuple(combinations(range(vertices), 2))


def basis(vertices=12):
    return _basis(_vertices(vertices))


def _ports(vertices, ports):
    if (type(ports) is not tuple or len(ports) != 2
            or any(type(p) is not int or not 0 <= p < vertices for p in ports)
            or ports[0] == ports[1]):
        raise ValueError('two distinct original integer ports required')
    return tuple(sorted(ports))


@lru_cache(maxsize=1024)
def _bare_basis(vertices, ports):
    states = _basis(vertices)
    index = {state: i for i, state in enumerate(states)}
    outside = tuple(p for p in range(vertices) if p not in ports)
    supports = []
    for charge in range(3):
        local = tuple(combinations(ports, charge))
        spectators = tuple(combinations(outside, 2-charge))
        if not spectators:
            continue
        for x in local:
            for y in local:
                supports.append(tuple((index[tuple(sorted(x+z))],
                                       index[tuple(sorted(y+z))]) for z in spectators))
    return tuple(supports)


def bare_basis(vertices, ports):
    """Disjoint exact supports of the bare algebra's matrix-unit basis."""
    vertices = _vertices(vertices)
    return _bare_basis(vertices, _ports(vertices, ports))


def _component(value):
    if type(value) in (int, F):
        return F(value), F()
    if (type(value) is tuple and len(value) == 2
            and all(type(part) in (int, F) for part in value)):
        return tuple(map(F, value))
    raise ValueError('original rational or rational real/imaginary pair required')


def _operator(vertices, entries, *, hermitian=False):
    d = len(_basis(vertices))
    if type(entries) is not dict:
        raise ValueError('a plain exact sparse operator dictionary required')
    result = {}
    for address, value in entries.items():
        if (type(address) is not tuple or len(address) != 2
                or any(type(i) is not int or not 0 <= i < d for i in address)):
            raise ValueError('original integer matrix addresses in the shell required')
        value = _component(value)
        if value != ZERO:
            result[address] = value
    if hermitian:
        for (i, j), (real, imag) in result.items():
            if result.get((j, i), ZERO) != (real, -imag):
                raise ValueError('an exact Hermitian generator required')
    return result


def _member(entries, supports):
    allowed = set()
    for support in supports:
        coefficient = entries.get(support[0], ZERO)
        if any(entries.get(address, ZERO) != coefficient for address in support):
            return False
        allowed.update(support)
    return all(address in allowed for address in entries)


def in_bare_algebra(vertices, ports, entries):
    vertices = _vertices(vertices)
    return _member(_operator(vertices, entries), bare_basis(vertices, ports))


def _commutator(entries, support):
    out = {}
    for (i, j), value in entries.items():
        for a, b in support:
            if j == a:
                old = out.get((i, b), ZERO)
                out[i, b] = old[0]+value[0], old[1]+value[1]
            if b == i:
                old = out.get((a, j), ZERO)
                out[a, j] = old[0]-value[0], old[1]-value[1]
    return {address: value for address, value in out.items() if value != ZERO}


def directly_normalizes(vertices, entries):
    """Independent exact commutator test on every bare matrix-unit support."""
    vertices = _vertices(vertices)
    operator = _operator(vertices, entries, hermitian=True)
    for ports in combinations(range(vertices), 2):
        supports = bare_basis(vertices, ports)
        if any(not _member(_commutator(operator, support), supports) for support in supports):
            return False
    return True


class OffDiagonalWitness(NamedTuple):
    row: int
    column: int
    number_port: int
    algebra_partner: int
    commutator_coefficient: int


@lru_cache(maxsize=14)
def _off_diagonal(vertices):
    states = tuple(frozenset(x) for x in _basis(vertices))
    witnesses = []
    for row, x in enumerate(states):
        for column, y in enumerate(states):
            if row == column:
                continue
            port = min(x-y)
            partner = next(j for j in range(vertices)
                           if j != port and not (x ^ y) <= {port, j})
            coefficient = int(port in y)-int(port in x)
            if not coefficient or (x-{port, partner}) == (y-{port, partner}):
                raise ArithmeticError('off-diagonal entry has no forbidden commutator witness')
            witnesses.append(OffDiagonalWitness(row, column, port, partner, coefficient))
    return tuple(witnesses)


def off_diagonal_witnesses(vertices=12):
    return _off_diagonal(_vertices(vertices))


@lru_cache(maxsize=14)
def _diagonal_constraints(vertices):
    states = _basis(vertices)
    index = {pair: i for i, pair in enumerate(states)}
    rows = []
    for i, j in combinations(range(vertices), 2):
        outside = [k for k in range(vertices) if k not in (i, j)]
        k = outside[0]
        for l in outside[1:]:
            row = [0]*len(states)
            for sign, pair in ((1, (i, k)), (-1, (j, k)), (-1, (i, l)), (1, (j, l))):
                row[index[tuple(sorted(pair))]] += sign
            rows.append(tuple(row))
    return tuple(rows)


def diagonal_constraints(vertices=12):
    """Integer equations enforcing a common coefficient for all spectators."""
    return _diagonal_constraints(_vertices(vertices))


def _row_echelon(rows):
    pivots = {}
    for values in rows:
        row = {i: F(value) for i, value in enumerate(values) if value}
        while row:
            pivot = min(row)
            if pivot not in pivots:
                divisor = row[pivot]
                pivots[pivot] = {i: value/divisor for i, value in row.items()}
                break
            scale = row[pivot]
            for i, value in pivots[pivot].items():
                row[i] = row.get(i, F())-scale*value
                if not row[i]:
                    del row[i]
    return pivots


def _onsite_rows(vertices):
    states = _basis(vertices)
    index = {pair: i for i, pair in enumerate(states)}
    first = [0]*len(states)
    for pair, coefficient in (((0, 1), 1), ((0, 2), 1), ((1, 2), -1)):
        first[index[pair]] = coefficient
    rows = [first]
    for i in range(1, vertices):
        row = [-v for v in first]
        row[index[0, i]] += 2
        rows.append(row)
    return tuple(tuple(row) for row in rows)


def _intersection_dimension(vertices, port):
    """Exact diagonal intersection, using equality classes of local records."""
    states = _basis(vertices)
    parent = list(range(len(states)))

    def find(i):
        while i != parent[i]:
            i = parent[i]
        return i

    for partner in range(vertices):
        if partner == port:
            continue
        groups = {}
        for i, state in enumerate(states):
            local = tuple(p for p in state if p in (port, partner))
            if local in groups:
                parent[find(i)] = find(groups[local])
            else:
                groups[local] = i
    return len({find(i) for i in range(len(states))})


def exact_certificate(vertices=12):
    vertices = _vertices(vertices)
    states = basis(vertices)
    rows = diagonal_constraints(vertices)
    witnesses = off_diagonal_witnesses(vertices)
    incidence = tuple(tuple(int(i in pair) for i in range(vertices)) for pair in states)
    inverse = _onsite_rows(vertices)
    # inverse/2 is a rational left inverse: rank(incidence)=N exactly.
    for i, row in enumerate(inverse):
        if any(sum(value*incidence[k][j] for k, value in enumerate(row)) != 2*int(i == j)
               for j in range(vertices)):
            raise ArithmeticError('pair-incidence left inverse failed')
    if any(sum(value*incidence[k][j] for k, value in enumerate(row))
           for row in rows for j in range(vertices)):
        raise ArithmeticError('additive energies fail a spectator equation')
    pivots = _row_echelon(rows)
    d = len(states)
    if len(pivots) != d-vertices or len(witnesses) != d*(d-1):
        raise ArithmeticError('exact compatibility constraints have the wrong dimension')
    encoded = [(w.row, w.column, w.number_port, w.algebra_partner, w.commutator_coefficient)
               for w in witnesses]
    digest = lambda value: hashlib.sha256(json.dumps(value, separators=(',', ':')).encode()).hexdigest()
    return {'ports': vertices, 'shell_dimension': d,
            'off_diagonal_ordered_entries_forced_zero': len(witnesses),
            'off_diagonal_witness_sha256': digest(encoded),
            'spectator_equations_checked': len(rows),
            'diagonal_constraint_sha256': digest(rows),
            'diagonal_exact_rank': len(pivots),
            'diagonal_pivot_columns': sorted(pivots),
            'number_operator_intersection_dimensions': [_intersection_dimension(vertices, i)
                                                       for i in range(vertices)],
            'hermitian_generator_real_dimension': vertices,
            'hermitian_real_constraints': d*d-vertices,
            'all_faithful_references_covered': True,
            'nonfaithful_references_classified': False,
            'numerical_rank_used': False}


def analyze_generator(vertices, entries):
    """Classify H exactly; exp(-H)/Tr(exp(-H)) is a faithful reference."""
    vertices = _vertices(vertices)
    operator = _operator(vertices, entries, hermitian=True)
    for witness in off_diagonal_witnesses(vertices):
        value = operator.get((witness.row, witness.column), ZERO)
        if value != ZERO:
            return {'compatible': False, 'obstruction': 'forbidden off-diagonal commutator entry',
                    'entry': [witness.row, witness.column],
                    'number_port': witness.number_port, 'algebra_partner': witness.algebra_partner,
                    'commutator_entry': [str(witness.commutator_coefficient*x) for x in value]}
    diagonal = tuple(operator.get((i, i), ZERO)[0] for i in range(len(basis(vertices))))
    for row in diagonal_constraints(vertices):
        value = sum((coefficient*entry for coefficient, entry in zip(row, diagonal, strict=True)), F())
        if value:
            return {'compatible': False, 'obstruction': 'spectator-dependent modular difference',
                    'integer_equation': list(row), 'residual': str(value)}
    onsite = tuple(sum((F(coefficient, 2)*entry for coefficient, entry in zip(row, diagonal, strict=True)), F())
                   for row in _onsite_rows(vertices))
    if any(diagonal[k] != onsite[i]+onsite[j] for k, (i, j) in enumerate(basis(vertices))):
        raise ArithmeticError('spectator equations did not reconstruct additive energies')
    return {'compatible': True, 'onsite_energies': list(map(str, onsite)),
            'generator': 'sum_i a_i N_i', 'reference_weights': 'lambda*exp(-a_i)*exp(-a_j)'}
