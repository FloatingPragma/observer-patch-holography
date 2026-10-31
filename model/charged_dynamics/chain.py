"""The complete joined four-site Weyl/contact chain and its charged band.

The link parameter multiplies only the middle bond. At one this is exactly
the existing uniform open four-site action. The product dimer state is an
initial state, not the coupled ground reference. No number states are removed.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
import math

import numpy as np

from ..full_fock_field.field import FieldGeometry, FullFock
from ..interacting_field.lattice import SIGMA
from ..particle_response import dimer


def frozen(value):
    value = np.asarray(value)
    return np.frombuffer(value.tobytes(), dtype=value.dtype).reshape(value.shape)


def _axis(axis):
    if type(axis) is not int or axis not in range(3):
        raise ValueError('integer axis in 0..2 required')
    literal = (np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]))
    if len(SIGMA) != 3 or any(not np.array_equal(value, expected) for value, expected in zip(SIGMA, literal, strict=True)):
        raise ArithmeticError('the committed Weyl spin chart changed')
    return axis


def _link(value):
    if type(value) not in (int, F):
        raise ValueError('exact middle-link strength in [0,1] required')
    numerator, denominator = (value, 1) if type(value) is int else (value.numerator, value.denominator)
    if type(numerator) is not int or type(denominator) is not int or denominator <= 0:
        raise ValueError('canonical integer rational components required')
    value = F(numerator, denominator)
    if not 0 <= value <= 1:
        raise ValueError('exact middle-link strength in [0,1] required')
    return value


def _kind(kind):
    if type(kind) is not str or kind not in ('hole', 'particle'):
        raise ValueError('the existing hole or particle doublet required')
    return kind


def geometry(axis=0):
    sides = [1, 1, 1]; sides[_axis(axis)] = 4
    return FieldGeometry(tuple(sides), speed=1/8, interaction=1.)


def _local(axis):
    sides = [1, 1, 1]; sides[_axis(axis)] = 2
    return dimer.setup(FieldGeometry(tuple(sides), speed=1/8, interaction=1.))


def ground_column(axis=0):
    """Analytic dimer ground with positive separated-singlet coefficient."""
    axis = _axis(axis)
    space = FullFock(4)
    c = [space.annihilator(j).toarray() for j in range(4)]
    turn = 1j*SIGMA[axis]
    d = c[:2]+[sum(turn[s, j]*c[2+j] for j in range(2)) for s in range(2)]
    create = [value.conj().T for value in d]
    blank = np.eye(16, dtype=complex)[:, 0]
    singlet = (create[0]@create[3]-create[1]@create[2])@blank/math.sqrt(2)
    doublon = (create[0]@create[1]+create[2]@create[3])@blank/math.sqrt(2)
    radical = math.sqrt(17)/4
    a = math.sqrt((1+1/radical)/2)
    b = -1/(8*radical*a)
    return frozen(a*singlet+b*doublon)


def local_columns(kind='hole', axis=0):
    """Same phase-fixed lower doublets as the existing dimer action."""
    _kind(kind); _axis(axis)
    result = np.zeros((16, 2), complex)
    result[[1, 2]] = np.eye(2)/math.sqrt(2)
    result[[4, 8]] = 1j*SIGMA[axis]/math.sqrt(2)
    if kind == 'particle':
        result = dimer.charge_conjugation(_local(axis))@result
    return frozen(result)


candidate_columns = local_columns


def bare_columns(kind='hole', axis=0):
    """A spin doublet at A or B, with the other block in its dimer ground."""
    candidate, ground = local_columns(kind, axis), ground_column(axis)
    return frozen(np.column_stack((np.kron(ground[:, None], candidate),
                                   np.kron(candidate, ground[:, None]))))


def gauge(axis=0):
    """Exact monomial Fock lift from d_j=(i sigma_axis)^j c_j coordinates."""
    axis = _axis(axis)
    images, phases = [], []
    for site in range(4):
        turn = np.linalg.matrix_power((1j*SIGMA[axis]).conj().T, site)
        for spin in range(2):
            occupied = np.flatnonzero(turn[:, spin])
            if len(occupied) != 1 or abs(turn[occupied[0], spin]) != 1:
                raise ArithmeticError('the fixed onsite gauge must be monomial')
            images.append(2*site+int(occupied[0]))
            phases.append(complex(turn[occupied[0], spin]))
    result = np.zeros((256, 256), complex)
    for address in range(256):
        occupied = [j for j in range(8) if address >> j & 1]
        moved = [images[j] for j in occupied]
        crossings = sum(a > b for i, a in enumerate(moved) for b in moved[i+1:])
        result[sum(1 << j for j in moved), address] = (-1)**crossings*math.prod(phases[j] for j in occupied)
    return frozen(result)


def hubbard_integer_parts(up, down):
    """Independent exact occupation matrices: 16 K(lambda)=A+lambda B.

    Fixed N_up/N_down is a symmetry block of the full carrier, not a physical
    support deletion. Entries and every hopping sign are Python integers.
    """
    if any(type(n) is not int or not 0 <= n <= 4 for n in (up, down)):
        raise ValueError('integer component populations in 0..4 required')
    addresses = tuple(mask for mask in range(256)
        if sum(mask >> (2*j) & 1 for j in range(4)) == up
        and sum(mask >> (2*j+1) & 1 for j in range(4)) == down)
    index = {address: j for j, address in enumerate(addresses)}
    base = [[0]*len(addresses) for _ in addresses]
    link = [[0]*len(addresses) for _ in addresses]
    for column, mask in enumerate(addresses):
        base[column][column] = 16*sum((mask >> (2*j) & 3) == 3 for j in range(4))-8*(up+down)
        for site in range(3):
            for spin in (0, 1):
                a, b = 2*site+spin, 2*(site+1)+spin
                if ((mask >> a) & 1) != ((mask >> b) & 1):
                    target = mask ^ (1 << a) ^ (1 << b)
                    sign = (-1)**sum(mask >> j & 1 for j in range(a+1, b))
                    (link if site == 1 else base)[index[target]][column] += sign
    return addresses, tuple(map(tuple, base)), tuple(map(tuple, link))


@dataclass(frozen=True)
class Chain:
    axis: int
    link_strength: F
    K: np.ndarray
    K0: np.ndarray
    V: np.ndarray
    N: np.ndarray
    NA: np.ndarray
    NB: np.ndarray
    ground_energy: float
    ground: np.ndarray
    ground_projector: np.ndarray
    energies: np.ndarray
    vectors: np.ndarray
    sectors: np.ndarray


def setup(axis=0, link_strength=F(1)):
    axis, link_strength = _axis(axis), _link(link_strength)
    space = FullFock(8)
    local = _local(axis)
    K0 = np.kron(np.eye(16), local.generator)+np.kron(local.generator, np.eye(16))
    ladders = [space.annihilator(j).toarray() for j in range(8)]
    transfer = sum(1j*SIGMA[axis][a, b]/16*ladders[2+a].conj().T@ladders[4+b]
                   for a in range(2) for b in range(2))
    V = transfer+transfer.conj().T
    K = K0+float(link_strength)*V
    if any(not np.isfinite(value).all() for value in (K0, V, K)):
        raise ArithmeticError('finite complete chain generators required')
    numbers = np.array([j.bit_count() for j in range(256)])
    nA = np.array([(j & 15).bit_count() for j in range(256)])
    rows = []
    for number in range(9):
        addresses = np.flatnonzero(numbers == number)
        values, vectors = np.linalg.eigh(K[np.ix_(addresses, addresses)])
        for energy, vector in zip(values, vectors.T, strict=True):
            full = np.zeros(256, complex); full[addresses] = vector
            rows.append((float(energy), number, full))
    rows.sort(key=lambda row: (row[0], row[1]))
    energies = np.array([row[0] for row in rows])
    vectors = np.column_stack([row[2] for row in rows])
    sectors = np.array([row[1] for row in rows])
    if (not np.isfinite(energies).all() or not np.isfinite(vectors).all() or sectors[0] != 4):
        raise ArithmeticError('complete chain spectral admission failed')
    with np.errstate(over='ignore', invalid='ignore'):
        residuals = (float(np.linalg.norm(K@vectors-vectors*energies)),
                     float(np.linalg.norm(vectors.conj().T@vectors-np.eye(256))))
    if any(not math.isfinite(value) or value > 2e-11 for value in residuals):
        raise ArithmeticError('complete chain spectral admission failed')
    ground = vectors[:, 0]
    return Chain(axis, F(link_strength.numerator, link_strength.denominator),
        *(frozen(value) for value in (K, K0, V, np.diag(numbers), np.diag(nA), np.diag(numbers-nA))),
        float(energies[0]), frozen(ground), frozen(np.outer(ground, ground.conj())),
        frozen(energies), frozen(vectors), frozen(sectors))


@dataclass(frozen=True)
class Band:
    energies: np.ndarray
    columns: np.ndarray
    projector: np.ndarray
    gap_to_complement: float


def band(kind='hole', axis=0, link_strength=F(1)):
    model = setup(axis, link_strength)
    selected = np.flatnonzero(model.sectors == (3 if _kind(kind) == 'hole' else 5))
    columns = model.vectors[:, selected[:4]]
    return Band(frozen(model.energies[selected[:4]]), frozen(columns),
                frozen(columns@columns.conj().T),
                float(model.energies[selected[4]]-model.energies[selected[3]]))


def dressed_columns(kind='hole', axis=0, link_strength=F(1)):
    projector = band(kind, axis, link_strength).projector
    bare = bare_columns(kind, axis)
    gram = bare.conj().T@projector@bare
    values, vectors = np.linalg.eigh(gram)
    if (not np.isfinite(values).all() or not np.isfinite(vectors).all()
            or np.any(values <= 1e-12)):
        raise ArithmeticError('bare-to-band polar chart is singular or unresolved')
    result = projector@bare@((vectors*(1/np.sqrt(values)))@vectors.conj().T)
    if not np.isfinite(result).all():
        raise ArithmeticError('dressed band chart failed its isometry admission')
    with np.errstate(over='ignore', invalid='ignore'):
        residual = float(np.linalg.norm(result.conj().T@result-np.eye(4)))
    if not math.isfinite(residual) or residual > 2e-11:
        raise ArithmeticError('dressed band chart failed its isometry admission')
    return frozen(result)
