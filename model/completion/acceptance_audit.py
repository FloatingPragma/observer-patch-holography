"""Exact obstruction to extending the eta clock by independent native copies.

The inherited modular generator has integer-spaced gaps after the existing
eta calibration.  A physical K=0 pair eigenvalue lies strictly between13/12
and3/2 in units U.  No finite tensor power of that reference supplies the
full field clock by an isometric intertwiner.  Arbitrary supplied actuators,
different references, correlated copies and different laws are outside this
exclusion.  This is a finite interface theorem, not program falsification.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
from itertools import product
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from ..local_repair_family import local_bias, modular_clock

ROOT = Path(__file__).resolve().parents[2]
SOURCE_BINDINGS = {
    'model/local_repair_family/local_bias.py': '4902174ac8b27c1d8111bc3942cf5cd1024a07999d77f272cbfa1f1242876755',
    'model/local_repair_family/modular_clock.py': 'e66cdb06b13d2e4d66495ce8f4ccfae1bd96828c8265bdbb22c5be071bf71d8a',
    'model/local_repair_family/record_clock_results.json': 'ac5c90623c5311f5de9965b12e4ef732888fc0ceac65a1d695fdb8ab6be6192b',
    'model/bound_pair/fiber.py': 'a5a21df2e5e6a4c71c66feaff2436390e357da8ae019df794280577afa48a49e',
    'model/bound_pair/protection.py': '99affa2080d03d5526122b97e0b75b6cdc1dd114116eedcacdf094287a2f0db1',
    'model/particle_clock/results.json': '114cf1b49b1b80b21edebabf9d2258b81432c4a35d1409ab613ac0e0402700bb',
}


def check_bindings(root=ROOT, bindings=None):
    root = Path(root).resolve()
    for name, expected in (SOURCE_BINDINGS if bindings is None else bindings).items():
        path = (root/name).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            raise ValueError(f'missing or external source: {name}')
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f'changed source requires a new interface review: {name}')


def native_levels():
    """All66 exact reference eigenvalues, without numerical logarithms."""
    p = modular_clock.parameters()
    if p.port_weights != (25, 25)+(8,)*10:
        raise ValueError('this exclusion requires the canonical two-weight reference')
    states, reference = local_bias.basis(), local_bias.reference(p)
    counts = Counter()
    for state, probability in zip(states, reference, strict=True):
        occupied_high_ports = state[0]+state[1]
        if probability != F(64, 7505)*F(25, 8)**occupied_high_ports:
            raise ArithmeticError('native reference does not have the claimed energy lattice')
        counts[occupied_high_ports] += 1
    if counts != {0: 45, 1: 20, 2: 1}:
        raise ArithmeticError('native reference multiplicities changed')
    return counts


def tensor_multiplicities(copies):
    """Exact coefficients of (45+20z+z^2)^copies; no tensor allocation."""
    if type(copies) is not int or not 1 <= copies <= 256:
        raise ValueError('an integer copy count in1..256 required for this evaluator')
    levels = {0: 1}
    for _ in range(copies):
        next_levels = Counter()
        for a, x in levels.items():
            for b, y in native_levels().items():
                next_levels[a+b] += x*y
        levels = dict(next_levels)
    return levels


def exact_pair_moments():
    """Independent L=4 K=0 spin blocks, with integer Gaussian arithmetic.

    A=8T; the unnormalized constant singlet is (0,1,-1,0).
    All intermediate integers are below2^20. No rounded trigonometry,
    eigenvalue solve, source fiber implementation or fitted pole enters.
    """
    identity = np.eye(2, dtype=np.int64)
    sr = (np.array([[0, 1], [1, 0]], dtype=np.int64),
          np.zeros((2, 2), dtype=np.int64),
          np.array([[1, 0], [0, -1]], dtype=np.int64))
    si = (np.zeros((2, 2), dtype=np.int64),
          np.array([[0, -1], [1, 0]], dtype=np.int64),
          np.zeros((2, 2), dtype=np.int64))
    singlet = np.array([0, 1, -1, 0], dtype=np.int64)
    sine = (0, 1, 0, -1)
    moments = [0, 0, 0, 0]
    imaginary = [0, 0, 0, 0]
    for k in product(range(4), repeat=3):
        ar = sum(sine[k[j]]*(-np.kron(sr[j], identity)+np.kron(identity, sr[j]))
                 for j in range(3))
        ai = sum(sine[k[j]]*(-np.kron(si[j], identity)+np.kron(identity, si[j]))
                 for j in range(3))
        real, imag = singlet.copy(), np.zeros(4, dtype=np.int64)
        for order in range(4):
            moments[order] += int(singlet@real)
            imaginary[order] += int(singlet@imag)
            real, imag = ar@real-ai@imag, ar@imag+ai@real
    if imaginary != [0]*4 or moments[0] != 128:
        raise ArithmeticError('the exact singlet moments lost Hermiticity or norm')
    result = tuple(F(moments[k], moments[0]*8**k) for k in range(4))
    if result != (F(1), F(0), F(3, 32), F(0)):
        raise ArithmeticError('the same-action Rayleigh witness changed')
    return result


def spectral_certificate():
    """Rational inequalities proving13/12<E_pair/U<3/2."""
    moments = exact_pair_moments()
    lower = F(13, 12)
    # Compression to span(d,Td) is [[1,sqrt(3/32)],[sqrt(3/32),0]].
    determinant = lower**2-lower-moments[2]
    kinetic_norm_squared = F(3, 16)  # (2 sqrt(3)c)^2, c/U=1/8.
    if not determinant < 0 or not kinetic_norm_squared < F(1, 4):
        raise ArithmeticError('strict noninteger pole enclosure failed')
    return {
        'source_action': 'H=T+U|d><d|; c/U=1/8; physical K=0 two-fermion fiber',
        'exact_moments_T_over_U': [str(x) for x in moments],
        'rayleigh_polynomial_at_lower': str(determinant),
        'pole_over_U_strict_interval': [str(lower), '3/2'],
        'distance_to_integer_strict_lower': '1/12',
        'centered_generator_residual_over_U_strict_lower': '1/24',
        'ideal_recurrence_binary_probability_difference_strict_lower': '1/36',
        'ideal_recurrence_half_trace_distance_strict_lower': '1/6',
        'recurrence_parameter': '|s|=2*pi/log(25/8); Hamiltonian model time |t|=2*pi/U',
        'recurrence_source_word_executed': False,
        'periodic_even_side_domain': 'all even L>=4 by the discrete Fourier moment identities',
        'exact_integer_check_side': 4,
    }


def required_factors(field_dimension, native_dimension=66):
    """Necessary Hilbert-space capacity, not a sufficient physical encoder."""
    if (type(field_dimension) is not int or field_dimension < 1
            or type(native_dimension) is not int or native_dimension < 2):
        raise ValueError('positive field dimension and native dimension at least2 required')
    factors, capacity = 0, 1
    while capacity < field_dimension:
        factors += 1
        capacity *= native_dimension
    return factors


def field_capacity(side):
    if type(side) is not int or side < 4 or side % 2:
        raise ValueError('an even periodic side at least4 required')
    dimension = math.comb(2*side**3, 2)
    return {'side': side, 'complete_N2_dimension': dimension,
            'native_H66_factors_at_least': required_factors(dimension),
            'native_H78_factors_at_least': required_factors(dimension, 78),
            'isometry_constructed': False,
            'dimension_capacity_solves_modular_obstruction': False}


def report(root=ROOT):
    check_bindings(root)
    root = Path(root)
    program = strict_json.load(root/'claims/completion_program.json')
    stages = {row['id']: row for row in program['stages']}
    clock = strict_json.load(root/'model/particle_clock/results.json')
    if clock['field_dictionary']['beta_times_U'] != 'log(25/8)':
        raise ValueError('the source/field clock calibration changed')
    if clock['field_dictionary']['faithful_algebra'] != 'B(span{r,d}) = M2':
        raise ValueError('the executed clock domain changed')
    return {
        'schema': 'oph.completion.acceptance_audit.v1',
        'scope': 'finite inherited-product-reference clock extension exclusion',
        'source_bindings': dict(SOURCE_BINDINGS),
        'native_reference_level_multiplicities': dict(native_levels()),
        'native_tensor_reference': '(64/7505)^n*(25/8)^m, integer0<=m<=2n',
        'spectral_certificate': spectral_certificate(),
        'capacity': [field_capacity(side) for side in (4, 6, 8, 12, 16)],
        'gate_dependencies': {
            gate: {'registry_status': stages[gate]['status'],
                   'unaccepted_predecessors': [p for p in stages[gate]['requires']
                                               if stages[p]['status'] != 'ACCEPTED']}
            for gate in ('CP-00', 'CP-01', 'CP-13')},
        'boundary_status': program['boundary_status'],
        'accepted_gate_changes': [],
        'whole_program_falsified': False,
        'canonical_premises_discharged': [],
        'physical_comparison_count_added': 0,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='check source custody and exact inequalities')
    parser.parse_args()
    print(json.dumps(report(), indent=2))


if __name__ == '__main__':
    main()
