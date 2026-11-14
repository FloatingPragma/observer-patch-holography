"""Outward full-source error certificate for the frozen LOCAL12 word.

Uses exact binary64 rational controls and the Q(sqrt(5)) response certificate.
No floating matrix residual is used as an error bound.
"""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PROVIDER = HERE.parent / 'native_walk_bridge'
sys.path.insert(0, str(ROOT / 'reverse-engineering-reality'))
sys.path.insert(0, str(PROVIDER))
import compile_walk as compiled
from model.certified_reference import certificate as cert, interval as iv, source
from model.observer.source_scalar import compiler

CONTROL_TABLE = {}


def remap(template, mapping):
    return [(k, tuple(mapping[j] for j in targets), a) for k, targets, a in template]


def exact_template_bound(word, template):
    end, fixed, variables = cert.consume(word, 0, template)
    assert end == len(word) and not variables
    return fixed + cert.controls_bound(word, CONTROL_TABLE)


def diagonal_template(pi_multiples):
    bits = (len(pi_multiples) - 1).bit_length()
    result = []
    for mask in range(len(pi_multiples)):
        coefficient = sum((F(x) * (-1)**((i & mask).bit_count())
                           for i, x in enumerate(pi_multiples)), F()) / len(pi_multiples)
        label = ''.join('Z' if mask >> j & 1 else 'I' for j in range(bits))
        result += cert.pauli(label, cert.PiMultiple(coefficient))
    return result


def uniform_bound(word):
    # Gray-order elimination of [1,1,1,0] gives these three two-plane
    # rotations. Parse their literal commuting factors independently.
    offset, bound = 0, F()
    rows = []
    for a, b, left_weight, right_weight in ((0, 1, 1, 2), (1, 3, 1, 1), (3, 2, 0, 1)):
        difference = a ^ b
        target = difference.bit_length() - 1
        orientation = (-1)**((a >> target) & 1)
        coefficients = []
        for mask in range(4):
            if mask & difference:
                continue
            label = ''.join('Y' if j == target else 'Z' if mask >> j & 1 else 'I'
                            for j in range(2))
            template = remap(cert.pauli(label), {0: 2, 1: 3})
            offset, fixed, variables = cert.consume(word, offset, template)
            assert len(variables) == 1
            bound += fixed
            coefficients.append((mask, variables[0]))
        theta = 2 * orientation * coefficients[0][1]
        assert all(value == orientation * theta * (-1)**((a & mask).bit_count()) / 2
                   for mask, value in coefficients)
        error = cert.two_vector_error(theta, left_weight, right_weight, False)
        bound += error
        rows.append({'plane': [a, b], 'normalized_squared_weights': [left_weight, right_weight],
                     'literal_ordered_plane_angle': str(theta), 'trigonometric_error_upper_exact': str(error)})
    for bit in (0, 1):
        template = [('I', (bit,), cert._pi(-1, 2)), ('H', (bit,), cert._pi(1, 2))]
        offset, fixed, variables = cert.consume(word, offset, template)
        assert not variables
        bound += fixed
    assert offset == len(word)
    return bound + cert.controls_bound(word, CONTROL_TABLE), rows


def run():
    custody_path = PROVIDER / 'COMMON_SOURCE_CUSTODY.json'
    custody = json.loads(custody_path.read_text())
    paths = [Path(__file__), Path(cert.__file__), Path(iv.__file__), Path(source.__file__),
             custody_path, PROVIDER / 'compile_walk.py', PROVIDER / 'common_walk.py',
             PROVIDER / 'COMMON_WALK_RECEIPT.json', PROVIDER / 'COMMON_WALK_PROGRAM.json']
    pins = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    for relative, expected in custody['files_sha256'].items():
        p = ROOT / relative
        assert hashlib.sha256(p.read_bytes()).hexdigest() == expected
    built = compiled.controls()
    template = [('H', (2,), cert._pi(1, 2))] + diagonal_template([0] * 7 + [1])
    template += [('H', (2,), cert._pi(-1, 2))]
    primitive = {
        'ccx': exact_template_bound(built['primitive_words']['ccx'], template),
        'cx': exact_template_bound(built['primitive_words']['cx'], cert.cnot(0, 1)),
        'flip': exact_template_bound(built['primitive_words']['flip'],
                    [('I', (0,), cert._pi(-1, 2)), ('X', (0,), cert._pi(1, 2))]),
        'hadamard': exact_template_bound(built['primitive_words']['hadamard'],
                    [('I', (0,), cert._pi(-1, 2)), ('H', (0,), cert._pi(1, 2))])}
    uniform, rotations = uniform_bound(built['preparation_word'])
    reflection_word = built['reflection_word']
    length = len(built['preparation_word'])
    assert reflection_word[:length] == compiler.inverse(built['preparation_word'])
    assert reflection_word[-length:] == built['preparation_word']
    diagonal = exact_template_bound(reflection_word[length:-length], diagonal_template([0] + [1] * 15))
    reflection = 2 * uniform + diagonal
    counts = Counter(kind for kind, targets in built['source_program'])
    swap = sum((number * primitive[kind] for kind, number in counts.items()), F())
    walk = swap + reflection
    preparation = 2 * uniform + primitive['hadamard']
    assert 0 < walk < F(1, 10**7)
    # A dropped scalar sign is rejected by the independent diagonal template.
    mutated = [dict(p, angle=-p['angle']) for p in reflection_word[length:-length]]
    # Negating phases is not a valid sign mutation (pi and -pi agree).
    # Instead replace diag(1,-1,...,-1) by its scalar negative explicitly.
    wrong = compiler.compile_diagonal([math.pi] + [0.] * 15)
    mutation = exact_template_bound(wrong, diagonal_template([0] + [1] * 15))
    assert mutation > 1
    assert pins == {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    result = {'schema': 'oph.local12.full_source_error.v1',
              'candidate_id': 'native-local12-current-v0-proposed',
              'checks': {'frozen_transitive_source_custody_verified': True,
                         'all_reversible_primitive_templates_consumed': True,
                         'uniform_three_component_rotations_certified': True,
                         'full_reflection_scalar_phase_certified': True,
                         'all_actual_swap_routine_uses_counted': True,
                         'integer_outward_Qsqrt5_response_bounds': True,
                         'global_sign_mutation_rejected': True, 'source_bytes_unchanged': True},
              'primitive_operator_errors_exact': {k: str(v) for k, v in primitive.items()},
              'uniform_preparation_operator_error_exact': str(uniform),
              'reflection_operator_error_exact': str(reflection),
              'swap_operator_error_exact': str(swap),
              'whole_W_operator_error_exact': str(walk),
              'whole_W_operator_error_upper': math.nextafter(float(walk), math.inf),
              'whole_source_preparation_error_exact': str(preparation),
              'whole_source_preparation_error_upper': math.nextafter(float(preparation), math.inf),
              'source_W_pulses': built['source_swap_native_pulses'] + len(reflection_word),
              'counts': dict(counts), 'uniform_rotations': rotations,
              'complete_channel_bound': 'Retained source/export Stinespring norm error <= delta; diamond norm <=2delta and total variation <=delta for arbitrary reference inputs in the declared clean work sector. For n steps use n*whole_W_error; add source preparation, apparatus and explicit reset timeout separately.',
              'physical_scope': 'Exact retained Q(sqrt(5)) response evaluated at literal binary64 integrated controls, with declared ideal coherent load/export service. No uncontrolled laboratory noise or physical time calibration is bounded.',
              'source_sha256': pins, 'official_rows_changed': []}
    (HERE / 'SOURCE_ERROR.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('checks', 'whole_W_operator_error_upper',
                                           'whole_source_preparation_error_upper', 'source_W_pulses')}))


if __name__ == '__main__':
    run()
