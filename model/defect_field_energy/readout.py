"""Charge-preserving grouped energy pointers and a Pauli negative control."""
from __future__ import annotations

from functools import lru_cache
from itertools import product
import math

import numpy as np
import sympy as s
from scipy.linalg import expm

from ..common_source.chiral_candidate import gaussian_source
from ..defect_field_bridge import joint
from ..defect_field_interactions import interaction
from ..observer.source_scalar import compiler
from .energy import small_energy_oracle

SOURCE_GUARD = 2e-7


def pauli_matrix(label):
    matrices = {'I': np.eye(2), 'X': np.array([[0, 1], [1, 0]]),
                'Y': np.array([[0, -1j], [1j, 0]]), 'Z': np.diag([1, -1])}
    result = np.ones((1, 1), complex)
    for letter in label[::-1]:
        result = np.kron(result, matrices[letter])
    return result


def pointer_word(label):
    """Supplementary Pauli pointer; individual words need not preserve charge."""
    if not isinstance(label, str) or len(label) != 4 or any(x not in 'IXYZ' for x in label):
        raise ValueError('one complete four-mode Pauli label required')
    basis, parity = [], []
    for bit, letter in enumerate(label):
        if letter == 'X':
            basis.append(compiler.pulse('H', math.pi/2, [bit]))
        elif letter == 'Y':
            basis += [compiler.pulse('Z', math.pi/4, [bit]),
                      compiler.pulse('H', math.pi/2, [bit])]
        if letter != 'I':
            parity += compiler.cnot(bit, 4)
    return basis+parity+compiler.inverse(basis)


@lru_cache(maxsize=256)
def pointer_gate(label):
    word = pointer_word(label)
    actual = compiler.execute(word, 5)
    pauli = pauli_matrix(label)
    plus, minus = (np.eye(16)+pauli)/2, (np.eye(16)-pauli)/2
    expected = np.block([[plus, minus], [minus, plus]])
    if np.linalg.norm(actual-expected, 2) > SOURCE_GUARD:
        raise ArithmeticError('complete Pauli pointer differs from its full32-column source law')
    return actual, word


def observe(label, state):
    state = np.asarray(state, complex)
    if state.ndim != 2 or state.shape[0] != 16:
        raise ValueError('the full16-row field state and all reference columns required')
    initial = np.vstack((state, np.zeros_like(state)))
    gate, word = pointer_gate(label)
    output = (gate@initial).reshape(2, 16, state.shape[1])
    probabilities = np.sum(abs(output)**2, axis=(1, 2))
    return {'amplitudes': output, 'probabilities': probabilities,
            'expectation': float(probabilities[0]-probabilities[1]),
            'resources': compiler.resources(word)}


@lru_cache(maxsize=1)
def exact_pauli_terms():
    """Exact support from the fixed rational-complex chart, with symbolic mass.

    No numerical coefficient threshold decides which energy term is read.
    All omitted Pauli coefficients are proved exactly zero in this oracle.
    """
    m = s.Symbol('m', real=True)
    half = s.Rational(1, 2)
    R = s.Matrix([[half-s.I*half, 0, 0, half+s.I*half],
                  [0, 1, 0, 0], [0, 0, 1, 0],
                  [-half+s.I*half, 0, 0, half+s.I*half]])
    bigR = s.kronecker_product(R, R)
    ladders = []
    for bit in range(4):
        a = s.zeros(16)
        for mask in range(16):
            if mask >> bit & 1:
                a[mask ^ (1 << bit), mask] = (-1)**(mask & ((1 << bit)-1)).bit_count()
        ladders.append(a)
    hopping = s.zeros(16)
    for a, b, value in ((0, 2, s.I/2), (1, 3, -s.I/2)):
        hopping += value*ladders[a].H*ladders[b]+s.conjugate(value)*ladders[b].H*ladders[a]
    hopping = bigR.H*hopping*bigR
    pair0 = m*s.diag(-1, 0, 0, 1)
    derivative = s.zeros(4); derivative[0, 3] = derivative[3, 0] = m
    h0 = s.kronecker_product(s.eye(4), pair0)+hopping/2
    h1 = s.kronecker_product(pair0, s.eye(4))+hopping/2
    matrices = {'H': h0+h1, 'h0': h0, 'h1': h1,
                'J0to1': s.I*(h0*h1-h1*h0),
                'drive_power_per_phase_rate': s.kronecker_product(s.eye(4), derivative)
                                                +s.kronecker_product(derivative, s.eye(4))}
    ring_current = s.zeros(16)
    for a, b in ((0, 2), (1, 3)):
        ring_current += -s.I/16*ladders[a].H*ladders[b]+s.I/16*ladders[b].H*ladders[a]
    matrices['ring_J0to2'] = bigR.H*ring_current*bigR
    pauli = {'I': s.eye(2), 'X': s.Matrix([[0, 1], [1, 0]]),
             'Y': s.Matrix([[0, -s.I], [s.I, 0]]), 'Z': s.diag(1, -1)}
    result = {}
    mass = small_energy_oracle()['energy']
    for name, matrix in matrices.items():
        terms = []
        reconstructed = s.zeros(16)
        for letters in product('IXYZ', repeat=4):
            label = ''.join(letters)
            P = s.kronecker_product(*(pauli[x] for x in label[::-1]))
            coefficient = s.simplify(s.trace(P*matrix)/16)
            if coefficient != 0:
                reconstructed += coefficient*P
                terms.append({'label': label, 'coefficient_symbolic': str(coefficient),
                              'coefficient': float(coefficient.subs(m, mass))})
        if any(s.simplify(x) != 0 for x in matrix-reconstructed):
            raise ArithmeticError('exact Pauli expansion lost an observable term')
        result[name] = tuple(terms)
    return result


@lru_cache(maxsize=1)
def exact_quadratic_terms():
    """Exact grouping into individually charge-invariant CAR observables."""
    m = s.Symbol('m', real=True)
    mass = small_energy_oracle()['energy']
    annihilators = []
    for mode in range(4):
        a = s.zeros(16)
        for mask in range(16):
            if mask >> mode & 1:
                a[mask ^ (1 << mode), mask] = (-1)**(mask & ((1 << mode)-1)).bit_count()
        annihilators.append(a)
    Q = s.diag(*[sum((1 if i % 2 == 0 else -1)*((mask >> i) & 1) for i in range(4)) for mask in range(16)])
    exact_pauli = {'I': s.eye(2), 'X': s.Matrix([[0, 1], [1, 0]]),
                   'Y': s.Matrix([[0, -s.I], [s.I, 0]]), 'Z': s.diag(1, -1)}
    result = {}
    for name, pauli_terms in exact_pauli_terms().items():
        matrix = s.zeros(16)
        for row in pauli_terms:
            coefficient = s.sympify(row['coefficient_symbolic'], locals={'m': m})
            P = s.kronecker_product(*(exact_pauli[x] for x in row['label'][::-1]))
            matrix += coefficient*P
        constant = s.simplify(matrix[0, 0])
        reconstruction = constant*s.eye(16)
        groups = []
        for i, a in enumerate(annihilators):
            value = s.simplify(matrix[1 << i, 1 << i]-constant)
            if value:
                operator = value*a.H*a
                groups.append(('number', (i,), value, operator))
                reconstruction += operator
            for j in range(i+1, 4):
                b = annihilators[j]
                for kind, coefficient, term in (
                    ('hopping', s.simplify(matrix[1 << i, 1 << j]), a.H*b),
                    ('pairing', s.simplify(matrix[(1 << i)+(1 << j), 0]), a.H*b.H)):
                    if coefficient:
                        operator = coefficient*term+s.conjugate(coefficient)*term.H
                        groups.append((kind, (i, j), coefficient, operator))
                        reconstruction += operator
        if any(s.simplify(x) != 0 for x in reconstruction-matrix):
            raise ArithmeticError('grouped CAR readout does not reconstruct the exact energy observable')
        terms = []
        for kind, modes, coefficient, operator in groups:
            if Q*operator != operator*Q:
                raise ArithmeticError('an energy readout group does not conserve field charge')
            z = complex(coefficient.subs(m, mass))
            terms.append({'kind': kind, 'modes': modes, 'coefficient': (z.real, z.imag),
                          'coefficient_symbolic': str(coefficient)})
        result[name] = {'constant': float(constant.subs(m, mass)),
                        'constant_symbolic': str(constant), 'terms': tuple(terms)}
    return result


def group_operator(kind, modes, coefficient):
    a = joint.ladders()
    z = complex(*coefficient)
    if kind == 'number':
        return z.real*a[modes[0]].conj().T@a[modes[0]]
    left, right = modes
    term = a[left].conj().T@(a[right] if kind == 'hopping' else a[right].conj().T)
    return z*term+z.conjugate()*term.conj().T


@lru_cache(maxsize=64)
def grouped_pointer(kind, modes, coefficient):
    """All64 source columns of a charge-preserving spectral-group pointer."""
    operator = group_operator(kind, modes, coefficient)
    z = complex(*coefficient)
    if kind == 'number':
        basis = []
        U = np.eye(16, dtype=complex)
        copied = modes
    elif kind in ('hopping', 'pairing'):
        left, right = modes
        route = []
        for neighbor in range(right-1, left, -1):
            route += compiler.remap(interaction.fermionic_swap_word(), [neighbor, neighbor+1])
        perpendicular = 1j*z/abs(z)
        diagonalizer = gaussian_source.local_gate(kind, perpendicular, math.pi/4)
        basis = route+compiler.remap(diagonalizer.word, [left, left+1])
        generator = group_operator(kind, modes, (perpendicular.real, perpendicular.imag))
        U = expm(-1j*math.pi/4*generator)
        copied = (left, left+1)
    else:
        raise ValueError('a number, hopping or neutral-pair energy group required')
    copy = []
    for bit, mode in enumerate(copied):
        copy += compiler.cnot(mode, 4+bit)
    word = basis+copy+compiler.inverse(basis)
    actual = compiler.execute(word, 6)
    diagonal = U@operator@U.conj().T
    if np.linalg.norm(diagonal-np.diag(np.diag(diagonal))) > 2e-10:
        raise ArithmeticError('charge-invariant CAR group was not diagonalized')
    projectors, eigenvalues = [], []
    for outcome in range(4):
        support = [mask for mask in range(16)
                   if sum(((mask >> mode) & 1) << bit for bit, mode in enumerate(modes)) == outcome]
        P = np.zeros((16, 16), complex)
        P[support, support] = 1
        projectors.append(U.conj().T@P@U)
        values = diagonal.diagonal()[support].real
        if len(values) and np.max(abs(values-values[0])) > 2e-10:
            raise ArithmeticError('pointer outcome is not a spectral subspace of the readout group')
        eigenvalues.append(float(values[0]) if len(values) else 0.)
    target = np.block([[projectors[a ^ b] for b in range(4)] for a in range(4)])
    if np.linalg.norm(actual-target, 2) > SOURCE_GUARD:
        raise ArithmeticError('complete grouped source pointer differs from its spectral instrument')
    Q = joint.charge()
    K = tuple(actual[16*i:16*(i+1), :16] for i in range(4))
    if any(np.linalg.norm(k@Q-Q@k, 2) > SOURCE_GUARD for k in K):
        raise ArithmeticError('grouped neutral-pointer instrument changed field charge')
    return {'unitary': actual, 'K': K, 'eigenvalues': tuple(eigenvalues),
            'word': word, 'operator': operator}


def observe_group(term, state):
    pointer = grouped_pointer(term['kind'], tuple(term['modes']), tuple(term['coefficient']))
    output = np.stack([K@state for K in pointer['K']])
    probabilities = np.sum(abs(output)**2, axis=(1, 2))
    return {'amplitudes': output, 'probabilities': probabilities,
            'expectation': float(probabilities@pointer['eigenvalues'])}
