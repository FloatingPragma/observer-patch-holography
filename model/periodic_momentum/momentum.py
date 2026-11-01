"""Exact periodic Weyl/contact translation charges on the complete CAR algebra.

The polynomial logarithms use Arg in (-pi,pi], including +pi at -1.
No finite matrix logarithm or eigenvalue tolerance defines a branch.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
from math import prod

import sympy as s

from ..coarse_energy_ward import ward as w
from ..lattice_momentum_balance.momentum import lift


def cleaned(matrix):
    return s.ImmutableSparseMatrix(matrix.rows, matrix.cols,
        {key: s.expand(value) for key, value in matrix.todok().items()})


def is_zero(matrix):
    return not cleaned(matrix).todok()


@dataclass(frozen=True)
class PeriodicGeometry:
    """Declared periodic boundary; lengths two are excluded (duplicate bonds)."""
    sides: tuple[int, int, int] = (3, 3, 3)
    speed: F = F(1, 8)
    interaction: F = F(1)

    def __post_init__(self):
        if (type(self.sides) is not tuple or len(self.sides) != 3 or
                any(type(n) is not int or n not in (1, *range(3, 65)) for n in self.sides) or
                prod(self.sides) > 64):
            raise ValueError('three periodic lengths 1 or 3..64, at most 64 sites required')
        if any(type(v) not in (int, F) or v < 0 for v in (self.speed, self.interaction)):
            raise ValueError('nonnegative exact rational source coefficients required')
        object.__setattr__(self, 'speed', F(self.speed))
        object.__setattr__(self, 'interaction', F(self.interaction))

    @property
    def sites(self):
        return tuple(product(*(range(n) for n in self.sides)))

    @property
    def modes(self):
        return 2*prod(self.sides)

    def shift(self, axis):
        if type(axis) is not int or not 0 <= axis < 3:
            raise ValueError('axis 0, 1 or 2 required')
        indices = {site: x for x, site in enumerate(self.sites)}
        result = []
        for site in self.sites:
            target = list(site)
            target[axis] = (target[axis]+1) % self.sides[axis]
            x = indices[tuple(target)]
            result.extend((2*x, 2*x+1))
        return tuple(result)


def one_particle(geometry):
    if not isinstance(geometry, PeriodicGeometry):
        raise ValueError('periodic source geometry required')
    pauli = (s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1))
    indices = {site: x for x, site in enumerate(geometry.sites)}
    h = s.MutableSparseMatrix(geometry.modes, geometry.modes, {})
    for x, site in enumerate(geometry.sites):
        for axis, length in enumerate(geometry.sides):
            if length == 1:
                continue
            neighbor = list(site); neighbor[axis] = (neighbor[axis]+1) % length
            y = indices[tuple(neighbor)]
            for a, b in product(range(2), repeat=2):
                z = s.I*s.Rational(geometry.speed)*pauli[axis][a, b]/2
                h[2*x+a, 2*y+b] += z
                h[2*y+b, 2*x+a] += s.conjugate(z)
    return s.ImmutableSparseMatrix(h)


def action(geometry):
    h0 = lift(one_particle(geometry))
    v = w.total((w.annihilator(2*x).adjoint()@w.annihilator(2*x)@
                 w.annihilator(2*x+1).adjoint()@w.annihilator(2*x+1)).scale(geometry.interaction)
                for x in range(len(geometry.sites)))
    return h0+v, h0, v


def valid_permutation(permutation):
    if (type(permutation) is not tuple or not permutation or
            any(type(x) is not int for x in permutation) or
            sorted(permutation) != list(range(len(permutation)))):
        raise ValueError('nonempty exact mode permutation required')


def permute_clifford(operator, permutation):
    """CAR automorphism gamma_(2p+a) -> gamma_(2 permutation[p]+a)."""
    valid_permutation(permutation)
    if any(mask.bit_length() > 2*len(permutation) for mask in operator.terms):
        raise ValueError('operator exceeds permutation carrier')
    terms = {}
    for mask, value in operator.terms.items():
        moved = [2*permutation[j//2]+j % 2 for j in range(mask.bit_length()) if mask >> j & 1]
        inversions = sum(a > b for i, a in enumerate(moved) for b in moved[i+1:])
        new_mask = sum(1 << j for j in moved)
        sign = (-1)**inversions
        terms[new_mask] = (sign*value[0], sign*value[1])
    return w.Clifford(terms)


def graded_columns(permutation, *, grading=True):
    """Signed permutation columns on all Fock states, including the vacuum."""
    valid_permutation(permutation)
    if len(permutation) > 12 or type(grading) is not bool:
        raise ValueError('at most 12 modes and a Boolean grading choice required')
    result = []
    for mask in range(1 << len(permutation)):
        moved = [permutation[j] for j in range(len(permutation)) if mask >> j & 1]
        inversions = sum(a > b for i, a in enumerate(moved) for b in moved[i+1:])
        result.append((sum(1 << j for j in moved), (-1)**inversions if grading else 1))
    return tuple(result)


def graded_matrix(permutation, *, grading=True):
    columns = graded_columns(permutation, grading=grading)
    return s.ImmutableSparseMatrix(len(columns), len(columns),
                                   {(row, col): sign for col, (row, sign) in enumerate(columns)})


def adjacent_word(permutation):
    """Exact adjacent FSWAP word for a mode permutation; no new primitive."""
    valid_permutation(permutation)
    inverse = [permutation.index(j) for j in range(len(permutation))]
    current = list(range(len(permutation)))
    word = []
    for destination, label in enumerate(inverse):
        origin = current.index(label)
        for left in range(origin-1, destination-1, -1):
            current[left], current[left+1] = current[left+1], current[left]
            word.append(left)
    if current != inverse:
        raise ArithmeticError('adjacent FSWAP synthesis failed')
    return tuple(word)


def word_columns(modes, word):
    if type(modes) is not int or not 1 <= modes <= 12 or any(
            type(j) is not int or not 0 <= j < modes-1 for j in word):
        raise ValueError('bounded adjacent FSWAP word required')
    result = []
    for mask in range(1 << modes):
        state, sign = mask, 1
        for left in word:
            a, b = (state >> left) & 1, (state >> (left+1)) & 1
            sign *= -1 if a and b else 1
            if a != b:
                state ^= 3 << left
        result.append((state, sign))
    return tuple(result)


def branch_labels(length):
    """Integer labels for Arg in (-pi,pi]; even Nyquist label is +L/2."""
    if type(length) is not int or not 3 <= length <= 127:
        raise ValueError('integer period 3..127 required')
    return tuple(j if 2*j <= length else j-length for j in range(length))


def centered_integer(value, length):
    if type(value) is not int:
        raise ValueError('integer total momentum label required')
    return branch_labels(length)[value % length]


def spectral_log(shift, length):
    """Exact P/pi and every spectral projector, for periods 3 or 4.

    General periods use the same proved finite Fourier formula; this exact
    full-matrix evaluator is intentionally limited to simple cyclotomic fields.
    """
    if length not in (3, 4) or shift.rows != shift.cols or shift.rows > 1024:
        raise ValueError('square bounded exact shift with period 3 or 4 required')
    identity = s.eye(shift.rows)
    if shift**length != identity:
        raise ValueError('shift does not obey its declared finite order')
    root = s.I if length == 4 else (-1+s.sqrt(3)*s.I)/2
    powers = [s.ImmutableSparseMatrix(identity)]
    for _ in range(1, length):
        powers.append(s.ImmutableSparseMatrix(powers[-1]*shift))
    projectors = []
    for j in range(length):
        matrix = s.MutableSparseMatrix(shift.rows, shift.cols, {})
        for r in range(length):
            coefficient = s.expand(root**((-j*r) % length))/length
            matrix += coefficient*powers[r]
        projectors.append(cleaned(matrix))
    p = sum((s.Rational(2*m, length)*projectors[j]
             for j, m in enumerate(branch_labels(length))),
            s.MutableSparseMatrix(shift.rows, shift.cols, {}))
    return cleaned(p), tuple(projectors)


def clifford_matrix(operator, modes):
    """Exact occupation matrix directly from Majorana bit actions."""
    if type(modes) is not int or not 1 <= modes <= 10:
        raise ValueError('one through ten modes required')
    result = {}
    for column in range(1 << modes):
        for mask, (a, b) in operator.terms.items():
            row, coefficient = column, s.Rational(a)+s.I*s.Rational(b)
            for gamma in reversed(range(mask.bit_length())):
                if not mask >> gamma & 1:
                    continue
                mode, kind = divmod(gamma, 2)
                occupied = (row >> mode) & 1
                coefficient *= (-1)**((row & ((1 << mode)-1)).bit_count())
                if kind:
                    coefficient *= -s.I if occupied else s.I
                row ^= 1 << mode
            result[row, column] = result.get((row, column), 0)+coefficient
    return cleaned(s.ImmutableSparseMatrix(1 << modes, 1 << modes, result))


def lift_matrix(matrix):
    """Independent exact occupation lift for a bounded general algebraic matrix."""
    modes = matrix.rows
    if matrix.cols != modes or not 1 <= modes <= 10:
        raise ValueError('bounded square one-particle matrix required')
    result = {}
    for column in range(1 << modes):
        for (p, q), z in matrix.todok().items():
            if not column >> q & 1:
                continue
            row = column ^ (1 << q)
            if row >> p & 1:
                continue
            sign = (-1)**((column & ((1 << q)-1)).bit_count()+(row & ((1 << p)-1)).bit_count())
            row ^= 1 << p
            result[row, column] = result.get((row, column), 0)+sign*z
    return cleaned(s.ImmutableSparseMatrix(1 << modes, 1 << modes, result))


def fourier_pair(length, first, second):
    """Normalized up/down pair in the source's exp(-ikx) Fourier convention."""
    if length != 4 or any(type(j) is not int for j in (first, second)):
        raise ValueError('period-four exact Fourier pair required')
    result = {}
    for x, y in product(range(length), repeat=2):
        a, b = 2*x, 2*y+1
        coefficient = s.I**((-first*x-second*y) % 4)/4
        result[(1 << a)+(1 << b), 0] = coefficient*(1 if a < b else -1)
    return s.ImmutableSparseMatrix(1 << (2*length), 1, result)


def principal_kernel(length, separation):
    """Closed one-particle kernel divided by pi, including the even diagonal."""
    branch_labels(length)
    if type(separation) is not int:
        raise ValueError('integer separation required')
    r = separation % length
    if not r:
        return s.Rational(1, length) if length % 2 == 0 else s.S.Zero
    theta = s.pi*s.Rational(r, length)
    if length % 2 == 0:
        return s.Rational((-1)**r, length)*(1+s.I*s.cot(theta))
    return s.I*(-1)**r/(length*s.sin(theta))


def prime_branch_certificate(length):
    """Minimal-polynomial obstruction for every spin-two logarithm branch.

    If an off-diagonal spin trace vanished, all integer coefficients
    a_j=2*j+L*n_j would be equal. Their residues 2*j mod L are distinct.
    The quotient calculation retains symbolic, unrestricted integers n_j.
    """
    if type(length) is not int or not 3 <= length <= 97 or not s.isprime(length):
        raise ValueError('odd prime period 3..97 required')
    x = s.Symbol('x')
    integers = s.symbols(f'n0:{length}', integer=True)
    polynomial = sum((2*j+length*integers[j])*x**j for j in range(length))
    phi = s.Poly(s.cyclotomic_poly(length, x), x)
    expected_phi = s.Poly(sum(x**j for j in range(length)), x)
    remainder = s.rem(polynomial, phi.as_expr(), x).expand()
    expected = sum((2*(j-length+1)+length*(integers[j]-integers[-1]))*x**j
                   for j in range(length-1))
    if phi != expected_phi or s.expand(remainder-expected) != 0:
        raise ArithmeticError('prime cyclotomic coefficient reduction failed')
    residues = tuple((2*j) % length for j in range(length))
    if len(set(residues)) != length:
        raise ArithmeticError('spin-trace branch residues failed to separate')
    return {'period': length, 'cyclotomic_coefficients': phi.all_coeffs(),
            'spin_trace_residues': residues,
            'constant_remainder_mod_period': (-2*(length-1)) % length,
            'every_nonzero_separation_is_primitive': True}


def smooth_sector_budget(maximum_number, cutoff, duration, interaction=F(1)):
    """Volume-independent fixed-number local-sine approximation bound.

    Omega restricts N<=maximum_number and every occupied k_j to |k_j|<=cutoff.
    The conservative rational domain N*cutoff<3 guarantees no principal-branch
    crossing since pi>3. The exact theorem permits N*cutoff<pi. This is an
    initial-state sector, not an invariant interacting truncation.
    """
    if (type(maximum_number) is not int or not 1 <= maximum_number <= 128 or
            any(type(v) not in (int, F) or v < 0 for v in (cutoff, duration, interaction)) or
            maximum_number*cutoff >= 3):
        raise ValueError('positive bounded number, nonnegative exact data and N*cutoff<3 required')
    n, kappa, t, u = maximum_number, F(cutoff), F(duration), F(interaction)
    initial = n*kappa**3/6
    leakage = min(F(1), t*u*(n//2))
    full_norm = F(22, 7)+n
    return {'initial_sine_error': initial, 'interacting_sector_leakage_upper': leakage,
            'evolved_sine_error_upper': min(full_norm, initial+full_norm*leakage)}
