"""Small executed generalization: finite Hermitian matrices to source words.

Exponential cost is explicit. This removes a target-unitary execution oracle
for the tested finite blocks; it does not execute the enormous spatial model.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import product
import json
import math
import numpy as np
import sympy as s
from scipy.linalg import expm

from .construction import HERE, require, source, yukawa_code


PAULI = {'I': s.eye(2), 'X': s.Matrix([[0, 1], [1, 0]]),
         'Y': s.Matrix([[0, -s.I], [s.I, 0]]), 'Z': s.diag(1, -1)}


def matrix(label):
    """Little endian exactly matches the source compiler's memory addresses."""
    out = s.Matrix([[1]])
    for letter in reversed(label):
        out = s.kronecker_product(out, PAULI[letter])
    return out


def decompose(hamiltonian):
    h = s.ImmutableMatrix(hamiltonian)
    dimension = h.rows
    if h.cols != dimension or dimension < 2 or dimension & (dimension - 1) or dimension > 16:
        raise ValueError('an exact 2, 4, 8 or 16 dimensional Hamiltonian required for this executor')
    if h != h.H or any(v.has(s.Float) for v in h):
        raise ValueError('an exact Hermitian Hamiltonian required')
    bits = dimension.bit_length() - 1
    terms = []
    rebuilt = s.zeros(dimension)
    for letters in product('IXYZ', repeat=bits):
        label = ''.join(letters)
        p = matrix(label)
        coefficient = s.simplify(s.trace(p * h) / dimension)
        if coefficient:
            if coefficient.is_real is not True:
                raise ValueError('provably real Pauli coefficients required')
            terms.append((label, coefficient))
            rebuilt += coefficient * p
    require(rebuilt == h, 'Pauli reconstruction failed')
    return terms


def anticommutes(a, b):
    return sum(x != y and x != 'I' and y != 'I' for x, y in zip(a, b, strict=True)) % 2 == 1


def product_word(hamiltonian, time=F(1, 3), steps=16):
    if type(time) not in (int, F) or type(steps) is not int or steps < 1:
        raise ValueError('rational time and positive integer steps required')
    terms = decompose(hamiltonian)
    one_step = [event for label, coefficient in terms
                for event in source.pauli_rotation(label, -float(time) * float(coefficient) / steps)]
    commutator_sum = s.simplify(sum(abs(a * b) for i, (la, a) in enumerate(terms)
                                  for lb, b in terms[i + 1:] if anticommutes(la, lb)))
    # ||U - product|| <= t^2/(2r) sum_{j<k} ||[H_j,H_k]||;
    # anticommuting Paulis give 2 |a_j a_k| exactly, others zero.
    bound = s.Rational(time.numerator, time.denominator) ** 2 * commutator_sum / steps
    return one_step * steps, {
        'time': str(time), 'steps': steps,
        'terms': [{'Pauli_little_endian': label, 'coefficient': str(c)} for label, c in terms],
        'exact_anticommutator_coefficient_sum': str(commutator_sum),
        'operator_error_upper': str(bound),
        'complete_instrument_half_diamond_upper': str(s.Min(1, bound)),
        'exact_decomposition': True,
        'source_pulses': len(one_step) * steps,
        'proper_code_service_events': 2 * len(one_step) * steps,
        'basis_changes_are_emitted_source_pulses': True,
        'physical_pulse_or_transport_errors_included': False,
        'intermediate_Pauli_controls_need_not_preserve_native_Gauss': True,
        'exact_Gauss_during_every_control_claimed': False,
    }


def prepare(vector):
    """Explicit source state synthesis, including every relative phase.

    The caller supplies the target coefficients. It is not installed as a
    state: the emitted unitary is applied to a separately reset blank.
    """
    vector = np.asarray(vector, dtype=complex)
    if vector.ndim != 1 or len(vector) < 2 or len(vector) & (len(vector) - 1) or len(vector) > 256:
        raise ValueError('power-of-two finite target through eight qubits required')
    if not np.isfinite(vector).all() or abs(np.linalg.norm(vector) - 1) > 1e-12:
        raise ValueError('normalized finite target required')
    return source.compile_real_state(abs(vector)) + source.compile_diagonal(np.angle(vector))


def run():
    _, _, _ = yukawa_code()
    from model.coherent_defect_wave.chiral_higgs.model import quantum_lepton_block
    _, _, support, block, _ = quantum_lepton_block()
    padded = s.zeros(8)
    padded[:5, :5] = block
    word, receipt = product_word(padded)
    actual = source.execute(word, 3)
    target = expm(-1j * np.asarray(padded.tolist(), complex) / 3)
    error = float(np.linalg.norm(actual - target, 2))
    require(error <= float(s.sympify(receipt['operator_error_upper'])) + 1e-12, 'product formula bound failed')
    state = np.array([1, 1j, -1, -1j, 2j, 2, -2j, -2], complex) / np.sqrt(20)
    prep = prepare(state)
    raw = source.execute(prep, 3, np.eye(8, dtype=complex)[:, 0])
    prep_error = float(np.linalg.norm(raw - state))
    require(prep_error < 2e-12, 'emitted complex state preparation failed')
    return {
        'schema': 'oph.join.finite_matrix_compiler.v1',
        'native_source_support': support, 'padding_addresses': [5, 6, 7],
        'padding_retained_in_every_source_execution': True,
        'product_compiler': receipt,
        'executed_full_input_columns': 8,
        'source_vs_target_operator_error_diagnostic': error,
        'raw_unitarity_defect': float(np.linalg.norm(actual.conj().T @ actual - np.eye(8), 2)),
        'source_word': word,
        'complex_preparation': {'source_pulses': len(prep), 'error_norm_diagnostic': prep_error,
                                'state_substituted_for_source_word': False},
        'general_finite_theorem': {
            'dimension': 'D <= 2^n; embed H as H direct_sum 0, retain padding outputs',
            'Pauli_terms_max': '4^n', 'source_pulses_per_Pauli_rotation_max': 'max(1,8n-5)',
            'steps_for_operator_error_epsilon': 'r >= t^2 sum_{anticommuting j<k}|h_j h_k|/epsilon',
            'total_source_pulses_max': 'r 4^n max(1,8n-5)',
            'memory_work_qubits': 'n; two active auxiliary source bits for native pointer, with its full native carrier retained',
            'matrix_entries_and_coherent_routing_are_required_inputs': True,
            'this_software_decomposition_cap': 4,
            'arbitrary_primary_spatial_H_constructed_or_executed': False,
            'efficient_scaling_claimed': False,
        },
    }
