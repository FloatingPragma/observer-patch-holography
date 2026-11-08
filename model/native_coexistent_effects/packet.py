"""Reproduce complete bounded native weak-measurement menus."""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np

from tools import strict_json
from model.native_effect_synthesis import packet as parent
from . import menus as m

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_FILES = tuple('model/native_coexistent_effects/'+name for name in (
    '__init__.py','menus.py','packet.py','test_menus.py','PROOF.md','README.md'))+(
    'Lean/EventAlgebra/OperationalBornBridge.lean','tools/strict_json.py','model/requirements.txt')


def source_bindings():
    paths = set(SOURCE_FILES)|set(parent.source_bindings())
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def example():
    first = tuple((F(i == 0), F()) for i in range(66))
    second = tuple((F(i in (0, 65)), F()) for i in range(66))
    return ((first, F(3, 5)), (second, F(4, 5)))


def matrix_rows(matrix):
    return [[i, j, str(x), str(y)] for (i, j), (x, y) in sorted(matrix.items())]


def bounds_record(stages):
    value = m.menu_bounds(stages)
    # Stable rational envelopes are checked against the full exact expression.
    if value['complete_halfdiamond_upper'] >= F(1, 50000):
        raise ArithmeticError('the complete bounded menu exceeds its error gate')
    if value['complete_gain_upper'] >= F(10001, 10000):
        raise ArithmeticError('the complete bounded menu exceeds its gain gate')
    return {'stages': stages, 'words': value['native_words'], 'pulses': value['native_pulses'],
            'complete_halfdiamond_strict_upper': '1/50000',
            'complete_gain_strict_upper': '10001/10000',
            'per_word_target_grid_units': [24, 23],
            'per_primitive_operator_gate': str(m.s.PULSE_GATE),
            'conditional_on_all_emitted_primitive_checks': True}


def ideal_output(stages, columns):
    data = m.ideal_menu(stages); blocks = []
    for (vector, survival), prefix in zip(stages, data['prefixes'], strict=True):
        norm = sum((x*x+y*y for x, y in vector), F())
        row = np.array([complex(float(x), -float(y)) for x, y in vector])/np.sqrt(float(norm))
        target = row@m.dense(prefix)@columns
        blocks.append(float(1-survival**2)*np.outer(target, target.conj()))
        blocks.append(np.zeros_like(blocks[-1]))
    return tuple(blocks), m.dense(data['residual'])@columns


def produce():
    before = source_bindings(); stages = example()
    exact = m.ideal_menu(stages); words = m.compile_menu(stages); checks = []
    m.s.certify_primitive.cache_clear()
    for index, pair in enumerate(words):
        for component, word in enumerate(pair):
            distinct = m.s.certify_word(word)
            admitted = m.s.interpreter.source_word(word, 9)
            stopped = m.s.interpreter.execute_word(admitted, np.array([[1.0]]))
            if stopped['status'] != 'RESOURCE_STOP' or stopped['pulses_executed'] != 0:
                raise ArithmeticError('the generic evaluator admission boundary changed')
            checks.append({'stage': index, 'component': component,
                           'word_sha256': parent.word_digest(word),
                           'distinct_emitted_primitives_checked': distinct,
                           'generic_evaluator_status': stopped['status']})
    commutator = m.matrix_sum(m.product(exact['effects'][0], exact['effects'][1]),
                             {key: m.scale(x, -1) for key, x in
                              m.product(exact['effects'][1], exact['effects'][0]).items()})
    if commutator[(0, 65)] != (F(216, 3125), F()):
        raise ArithmeticError('the coexistent example lost its noncommuting effects')
    rng = np.random.default_rng(271828)
    columns = rng.normal(size=(66, 2))+1j*rng.normal(size=(66, 2))
    columns /= np.linalg.norm(columns)
    aborts, residual = m.diagnostic_output(words, columns)
    ideal_aborts, ideal_residual = ideal_output(stages, columns)
    # These floating checks are diagnostics against proved fixed envelopes.
    if (np.linalg.norm(residual-ideal_residual) >= 1e-8
            or any(np.linalg.norm(a-b, 'nuc') >= 1e-8
                   for a, b in zip(aborts, ideal_aborts, strict=True))):
        raise ArithmeticError('the complete weak menu diagnostic exceeds its gate')
    value = {'schema': 'oph.native_coexistent_effects.v1',
             'construction': {'source_dimension': 66, 'work_qubits': 9, 'work_dimension': 512,
                              'stages_max': m.MAX_STAGES, 'source_words_per_stage': 2,
                              'pulses_per_source_word': m.WORD_PULSES,
                              'admission_pulse_ceiling': m.s.interpreter.MAX_WORD,
                              'input_rational_bits_max': m.s.MAX_RATIONAL_BITS,
                              'angle_grid_bits': m.iv.ANGLE_BITS,
                              'trigonometric_enclosure_bits': m.iv.BITS},
             'noncommuting_example': {
                 'stages': [{'vector': [[str(x), str(y)] for x, y in v], 'survival': str(s)}
                            for v, s in stages],
                 'effects': [matrix_rows(effect) for effect in exact['effects']],
                 'residual_operator': matrix_rows(exact['residual']),
                 'commutator_0_65': ['216/3125', '0'],
                 'exact_effect_sum_is_identity': True,
                 'coherent_source_and_reference_diagnostic_gate': '1/100000000',
                 'actual_word_checks': checks, 'outcomes': m.outcomes(words)},
             'uniform_bounded_menu_certificate': bounds_record(m.MAX_STAGES),
             'scope': {'complete_noncommuting_coexistent_family': True,
                       'successful_Luders_state_update_derived': True,
                       'all_processor_and_padding_failures_retained': True,
                       'final_quantum_residual_retained': True,
                       'arbitrary_untouched_reference_in_proof': True,
                       'stored_primitive_acceptance_required': True,
                       'rounded_maps_asserted_exactly_trace_preserving': False,
                       'floating_propagation_arrays_claimed_certified': False,
                       'free_classical_randomizer_or_additional_register': False,
                       'generic_host_limits_raised': False,
                       'arbitrary_POVM_or_universal_coexistence_compiler': False,
                       'arbitrarily_small_error_at_fixed_budget': False,
                       'Born_semantics_derived_from_bare_repair': False,
                       'growing_capacity_quantum_refinement': False,
                       'physical_QFT_or_particle_identification': False},
             'source_bindings': before}
    if source_bindings() != before:
        raise ArithmeticError('source inputs changed during complete menu replay')
    return value


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def verify(path=OUTPUT):
    if canonical(strict_json.load(Path(path))) != canonical(produce()):
        raise ValueError('native coexistent-effect receipt differs from source replay')
    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(produce()))
        print('Native coexistent-effect receipt written.')
    else:
        verify()
        print('Native coexistent-effect receipt verified.')
