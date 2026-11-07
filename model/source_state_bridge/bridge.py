"""Exact state, modular-frequency and quantum-carrier boundary checks.

These concern the named biased native member, not every OPH source family.
The mathematical proofs are in PROOF.md; finite controls do not replace them.
"""
from collections import Counter
from fractions import Fraction as F
import json
from math import prod
from pathlib import Path

from ..local_repair_family import local_bias as law
from ..native_program import interpreter

ROOT = Path(__file__).resolve().parents[2]
SOURCE_FILES = (
    'model/source_state_bridge/__init__.py',
    'model/source_state_bridge/bridge.py',
    'model/source_state_bridge/test_bridge.py',
    'model/source_state_bridge/PROOF.md',
    'model/local_repair_family/local_bias.py',
    'model/local_repair_family/mixing.py',
    'model/local_repair_family/learned_flow.py',
    'model/local_repair_family/record_clock.py',
    'model/local_repair_family/README.md',
    'model/native_program/interpreter.py',
    'model/native_program/members.py',
    'model/interacting_extension/higgs_join/joined_basis.json',
)


def positive_integer(value):
    if type(value) is not int or value < 1:
        raise ValueError('a positive integer dimension is required')
    return value


def native_spectrum():
    """Derive the actual member's weights from its program and native basis."""
    value = interpreter.program()
    weights = tuple(value['parameters']['port_weights'])
    if weights != (25, 25)+(8,)*10 or value['carrier'] != {'kind': 'local', 'shell': 'H66'}:
        raise ArithmeticError('the named native source changed')
    parameters = law.Parameters(port_weights=weights)
    states = law.basis()
    numerators = tuple(prod(w**n for w, n in zip(weights, x, strict=True)) for x in states)
    normalizer = sum(numerators)
    probabilities = tuple(F(x, normalizer) for x in numerators)
    if len(states) != 66 or probabilities != law.reference(parameters):
        raise ArithmeticError('the native reference differs from its product law')
    # Remove the scalar -log(625/Z). The ideal target is delta*levels,
    # where delta=log(25/8)>0; no floating logarithm is used.
    levels = tuple(2-x[0]-x[1] for x in states)
    if any(F(625, n) != F(25, 8)**level for n, level in zip(numerators, levels, strict=True)):
        raise ArithmeticError('the exact modular level identity failed')
    census = Counter(numerators)
    if census != {625: 1, 200: 20, 64: 45} or normalizer != 7505:
        raise ArithmeticError('native reference spectrum changed')
    gaps = Counter(a-b for a in levels for b in levels)
    if gaps != {-2: 45, -1: 920, 0: 2426, 1: 920, 2: 45}:
        raise ArithmeticError('complete modular frequency census changed')
    return {'weights': weights, 'normalizer': normalizer, 'numerators': numerators,
            'probabilities': probabilities, 'levels': levels, 'gaps': gaps}


def support_state_distance_bound(rank):
    """Lower bound against rank-at-most-r targets for the full isometric image.

    The input is the exact stationary state, not its finite approximation.
    Extra residual quantum systems are part of the compared full state.
    """
    positive_integer(rank)
    probabilities = sorted(native_spectrum()['probabilities'], reverse=True)
    return F(1)-sum(probabilities[:rank], F())


def carrier_channel_distance_bound(target_dimension, quantum_dimension=66):
    """Half-diamond lower bound for transmission through a quantum carrier.

    Arbitrary retained classical flags are allowed. Additional transmitted
    quantum systems and entanglement-assisted side channels are excluded.
    The comparison is all-input identity transmission, not one preparation.
    """
    target_dimension = positive_integer(target_dimension)
    quantum_dimension = positive_integer(quantum_dimension)
    return max(F(), F(1)-F(quantum_dimension, target_dimension))


def exact_flagged_code(target_dimension, quantum_dimension):
    """Small explicit encode/decode control, retaining each classical flag."""
    import sympy as s
    positive_integer(target_dimension); positive_integer(quantum_dimension)
    if target_dimension > 16 or quantum_dimension > 16:
        raise ValueError('the exact finite control admits dimensions through 16')
    rows = []
    for start in range(0, target_dimension, quantum_dimension):
        width = min(quantum_dimension, target_dimension-start)
        encode = s.zeros(quantum_dimension, target_dimension)
        for j in range(width):
            encode[j, start+j] = 1
        # The main decoder is isometric on the reached code. Unreached
        # carrier directions are reset to the same output basis state.
        decoder = [encode.T]
        for j in range(width, quantum_dimension):
            extra = s.zeros(target_dimension, quantum_dimension)
            extra[start, j] = 1
            decoder.append(extra)
        rows.append((encode, tuple(decoder)))
    return tuple(rows)


def certificate():
    source = native_spectrum()
    basis = json.loads((ROOT/'model/interacting_extension/higgs_join/joined_basis.json').read_text())
    dimensions = []
    for charge in ('0', '1'):
        row = basis[charge]
        d = row['dimension']
        if type(d) is not int or any(len(row[key]) != d for key in ('labels', 'raw_norms', 'basis')):
            raise ArithmeticError('the complete Higgs basis inventory is inconsistent')
        dimensions.append(d)
    if dimensions != [438, 1596]:
        raise ArithmeticError('the declared Higgs sectors changed')
    spectrum = Counter(source['numerators'])
    return {
        'schema': 'oph.source_state_bridge.direct_bounds.v1',
        'native_port_weights': list(source['weights']),
        'native_dimension': 66,
        'stationary_normalizer': source['normalizer'],
        'stationary_spectrum': [{'numerator': n, 'multiplicity': spectrum[n],
                                'eigenvalue': str(F(n, source['normalizer']))}
                               for n in sorted(spectrum, reverse=True)],
        'ideal_modular_gap': 'log(25/8)',
        'ideal_modular_levels_modulo_scalar': [0, 1, 2],
        'ideal_Bohr_multiplicities': [{'gap_multiple': gap, 'multiplicity': count}
                                    for gap, count in sorted(source['gaps'].items())],
        'complete_operator_directions': sum(source['gaps'].values()),
        'full_isometric_state_to_pure_target_trace_distance_lower': str(support_state_distance_bound(1)),
        'full_isometric_state_to_rank20_target_trace_distance_lower': str(support_state_distance_bound(20)),
        'all_input_channel_half_diamond_lower': {
            str(d): str(carrier_channel_distance_bound(d)) for d in dimensions},
        'classical_records_allowed_in_channel_bound': True,
        'scope': {
            'port_weights_are_declared_source_parameters': True,
            'pure_target_compares_full_state_not_only_a_reduced_subsystem': True,
            'additional_quantum_memory_or_active_preparation_excluded_by_state_bound': False,
            'channel_bound_covers_one_selected_preparation_only': False,
            'clock_frequencies_are_ideal_modular_target_not_exact_compiled_pulses': True,
            'Higgs_ground_uniqueness_interval_certified': False,
            'general_multicarrier_or_encoded_field_identification_excluded': False,
            'source_selected_interacting_Higgs_action_or_state_derived': False,
            'physical_particle_mass_or_clock_derived': False,
        },
    }
