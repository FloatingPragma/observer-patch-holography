"""Exact sparse preparation identity and all-line dark spectral response."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
import math

import numpy as np
import sympy as s

from ..common_spin_protection import spin
from ..repair_spectral_sectors import algebra
from . import circuits


@lru_cache(maxsize=1)
def target_numerator():
    # The four empty sites have eta-pair signs +,-,+,-.  The base contains
    # two spin-up particles on sites0,1; pair creation crosses even parity.
    hubbard = np.zeros(4096, dtype=np.int64)
    for site in range(2, 6):
        hubbard[5 | (3 << (2*site))] = (-1)**site
    old, pairs = algebra.dark_columns()
    if pairs[0] != (0, 1) or not np.array_equal(old[:, 0].toarray().ravel(), hubbard):
        raise ArithmeticError('sparse source target is not the exact old dark column')
    frame = spin.operators().frame.adjoint()
    real = frame.real@hubbard
    imag = frame.imag@hubbard
    return hubbard, real, imag


def target():
    _, real, imag = target_numerator()
    return (real+1j*imag)/2


@lru_cache(maxsize=1)
def local_gate_certificate():
    identity, x, y, z = s.eye(2), s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1)
    paulis = {'I': identity, 'X': x, 'Y': y, 'Z': z}
    def rotation(label, multiple):
        p = s.kronecker_product(*(paulis[letter] for letter in reversed(label)))
        return s.cos(s.pi*multiple)*s.eye(p.rows)+s.I*s.sin(s.pi*multiple)*p
    # Pauli factors commute, but multiply in chronological order anyway.
    split = rotation('XY', -s.Rational(1, 8))*rotation('YX', s.Rational(1, 8))
    root = s.sqrt(2)/2
    expected = s.Matrix([[1, 0, 0, 0], [0, root, -root, 0], [0, root, root, 0], [0, 0, 0, 1]])
    if (split-expected).applyfunc(s.simplify) != s.zeros(4):
        raise ArithmeticError('two commuting Pauli rotations do not implement the claimed split')
    cnot = s.eye(4)
    for label, angle in [('II', -s.Rational(1, 4)), ('ZI', s.Rational(1, 4)),
                         ('IX', s.Rational(1, 4)), ('ZX', -s.Rational(1, 4))]:
        cnot = rotation(label, angle)*cnot
    expected_cnot = s.Matrix([[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]])
    if (cnot-expected_cnot).applyfunc(s.simplify) != s.zeros(4):
        raise ArithmeticError('phase-exact copy gate identity failed')
    return {'split_identity_exact': True, 'copy_identity_exact': True,
            'controlled_rotation_identity': 'C exp(i theta P)=exp(i theta P/2) exp(-i theta P Z_pointer/2)',
            'scalar_phases_paid_before_controlling': True}


def preparation_certificate():
    local = local_gate_certificate()
    # Exact sparse execution of the three proved split gates.
    amplitudes = {21: s.Integer(1)}  # up occupations0,1,2
    for first, second in ((4, 8), (4, 6), (8, 10)):
        out = {}
        for mask, amplitude in amplitudes.items():
            bits = ((mask >> first) & 1, (mask >> second) & 1)
            if bits in ((0, 0), (1, 1)):
                out[mask] = out.get(mask, 0)+amplitude
            else:
                exchanged = mask ^ (1 << first) ^ (1 << second)
                out[mask] = out.get(mask, 0)+amplitude/s.sqrt(2)
                out[exchanged] = out.get(exchanged, 0)+amplitude/s.sqrt(2)*(1 if bits == (1, 0) else -1)
        amplitudes = out
    copied = {}
    for mask, amplitude in amplitudes.items():
        for bit in (4, 6, 8, 10):
            if mask >> bit & 1: mask ^= 1 << (bit+1)
        sign = (-1)**(((mask >> 6) & 1)+((mask >> 10) & 1))
        copied[mask] = s.simplify(2*sign*amplitude)
    hubbard, real, imag = target_numerator()
    expected = {int(j): s.Integer(hubbard[j]) for j in np.flatnonzero(hubbard)}
    if copied != expected or int(np.dot(real, real)+np.dot(imag, imag)) != 4:
        raise ArithmeticError('the sparse deterministic preparation differs from D|0,1>')
    plan = circuits.preparation_plan()
    controlled = circuits.controlled(plan)
    for (p, theta), ((a, left), (b, right)) in zip(plan, zip(controlled[::2], controlled[1::2]), strict=True):
        if a != p+'I' or b != p+'Z' or left+right != 0 or left-right != theta:
            raise ArithmeticError('pointer control changes the zero branch or drops a relative phase')
    return {**local, 'field_modes': 12, 'pointer_modes': 1,
            'input': 'initialized field vacuum and pointer zero',
            'exact_target': 'frame^dagger eta_plus|0up,1up>/2',
            'Hubbard_occupation_numerators': [[j, int(v)] for j, v in sorted(copied.items())],
            'common_amplitude_denominator': 2,
            'split_gates': 3, 'copied_pair_occupations': 4,
            'abstract_Pauli_factors': len(plan),
            'generic_amplitude_tree_used': False,
            'charge_postselection_used': False,
            'target_number': 4, 'target_spin': 1, 'target_eta': 2,
            'physical_charge_phase_reference_derived': False}


def spectral_lines():
    rows = []
    for p, q in combinations(range(1, 7), 2):
        a, b = p*math.pi/7, q*math.pi/7
        weight = 16/49*math.sin(a)**2*math.sin(b)**2*(math.cos(b)-math.cos(a))**2
        rows.append({'p': p, 'q': q, 'H_energy': 1+(math.cos(a)+math.cos(b))/8,
                     'K_energy': -1+(math.cos(a)+math.cos(b))/8,
                     'weight': weight, 'exact_energy_group': 'central' if p+q == 7 else f'{p},{q}'})
    if min(row['weight'] for row in rows) <= 0 or abs(math.fsum(row['weight'] for row in rows)-1) > 2e-15:
        raise ArithmeticError('all-line normalized spectral response check failed')
    return rows


def amplitude(time):
    return sum(row['weight']*np.exp(-1j*time*row['H_energy']) for row in spectral_lines())


def spectral_certificate():
    algebra.dark_certificate()  # Retain the exact distinct-root proof.
    _, real, imag = target_numerator()
    h = spin.operators().H
    image_real, image_imag = h.real@real-h.imag@imag, h.real@imag+h.imag@real
    mean = F(int(real@image_real+imag@image_imag), 4*spin.operators().scale)
    second = F(int(image_real@image_real+image_imag@image_imag), 4*spin.operators().scale**2)
    if mean != 1 or second-mean**2 != F(1, 256):
        raise ArithmeticError('the native target lost its exact nonzero Hamiltonian variance')
    rows = spectral_lines()
    groups = {}
    for row in rows:
        groups.setdefault(row['exact_energy_group'], []).append(row)
    grouped = [{'group': key, 'mode_pairs': [[r['p'], r['q']] for r in members],
                'H_energy': members[0]['H_energy'],
                'weight': math.fsum(r['weight'] for r in members)} for key, members in groups.items()]
    return {'energy_reference': 'actual H=K+N/2; vacuum H=0, dark N4 energies are K+2',
            'exact_mean_H_energy': str(mean), 'exact_H_energy_variance': str(second-mean**2),
            'H_energy_formula': '1+(cos(p*pi/7)+cos(q*pi/7))/8',
            'weight_formula': '(16/49)*sin(p*pi/7)^2*sin(q*pi/7)^2*(cos(q*pi/7)-cos(p*pi/7))^2',
            'exact_positive_weight_reason': 'sine strictly positive and cosine strictly decreasing on (0,pi)',
            'exact_normalization_reason': 'orthonormal sine orbitals and exterior-power Parseval identity',
            'all15_mode_pairs': rows, 'all13_distinct_energy_groups': grouped,
            'Ramsey_formula': 'P0(phi,t)=(1+Re(exp(i phi)*C_H(t)))/2; P1=1-P0',
            'characteristic_function': 'C_H(t)=sum_(p<q) w_pq exp(-i t H_energy_pq)',
            'independent_determinant': 'C_H(t)=exp(-it)*det([exp(-it A_chain/16)]_{0,1;0,1})',
            'spectral_lines_removed_by_threshold': 0,
            'absolute_physical_mass_identified': False}
