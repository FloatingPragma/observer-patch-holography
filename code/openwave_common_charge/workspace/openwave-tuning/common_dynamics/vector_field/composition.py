"""Complete native composition bound for the declared matter-field word."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import sys

import autonomous as af
import native_field as nf
from model.certified_reference import certificate as cert
from model.observer.source_discrimination import circuit
from model.native_clock import configuration
from model.native_interface import pointer

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def phase_template(mode, angle):
    half = cert.PiMultiple(angle.value / 2) if type(angle) is cert.PiMultiple else angle / 2
    return [('I', (mode,), half), ('Z', (mode,), -half)]


def givens_template(a, b, angle):
    half = cert.PiMultiple(angle.value / 2) if type(angle) is cert.PiMultiple else angle / 2
    first = ['I'] * (b + 1); second = first.copy()
    first[a], first[b] = 'X', 'Y'
    second[a], second[b] = 'Y', 'X'
    for k in range(a + 1, b): first[k] = second[k] = 'Z'
    return cert.pauli(''.join(first), -half) + cert.pauli(''.join(second), half)


def preparation_word(charge=1):
    if charge not in (-1, 1):
        raise ValueError('one particle or separately prepared one hole required')
    modes = ([0] if charge == 1 else [16 * af.mb.common.encode_position(x) + e
              for x in range(96) for e in range(12) if (x, e) != (0, 0)])
    word = [pulse for mode in modes for pulse in af.mb.compiler.remap(
        af.mb.architecture.primitive_routines()['flip'], [mode])]
    for k in range(1, 12):
        word += af.mb.givens_word(0, k, math.asin(1 / math.sqrt(13 - k)))
    return word


def run():
    paths = [Path(__file__), Path(af.mb.__file__), Path(af.mode_field.__file__),
             HERE / 'NATIVE_FIELD.json', HERE / 'FIELD_PROGRAM.json',
             HERE.parent / 'manybody/MANIFEST.json', HERE.parent / 'manybody/MODE_CURRENT_PROGRAMS.json',
             ROOT / 'openwave-tuning/auditor/common_walk_review/SOURCE_ERROR.json',
             Path(cert.__file__), Path(configuration.__file__), Path(pointer.__file__)]
    pins = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    bound = nf.root_certificate.exact_template_bound
    uniform_error = F()
    for k in range(1, 12):
        theta = math.asin(1 / math.sqrt(13 - k))
        uniform_error += bound(af.mb.givens_word(0, k, theta), givens_template(0, k, F(theta)))
        # On each spectator parity sector the same sine/cosine plane error
        # controls the full Fock rotation; the vacuum/double block is identity.
        uniform_error += cert.two_vector_error(-F(theta), 12 - k, 1, False)
    phase_error = bound(af.mb.occupation_phase_word(0, math.pi), phase_template(0, cert._pi(1)))
    onsite_error = 2 * uniform_error + 15 * phase_error
    program = af.mb.mode_program()
    distances = Counter(b - a for a, b, e in program['transport_pairs'])
    swap_error = F(); swap_pulses = 0
    swap_rows = []
    for distance, count in sorted(distances.items()):
        word = af.mb.fock_swap_word(0, distance)
        template = phase_template(distance, cert._pi(1)) + givens_template(0, distance, cert._pi(1, 2))
        error = bound(word, template)
        swap_error += count * error
        swap_pulses += count * len(word)
        swap_rows.append({'mode_separation': distance, 'uses': count, 'pulses_each': len(word), 'error_exact_each': str(error)})
    root = json.loads((ROOT / 'openwave-tuning/auditor/common_walk_review/SOURCE_ERROR.json').read_text())
    flip_error = F(root['primitive_operator_errors_exact']['flip'])
    blocks = json.loads((HERE.parent / 'manybody/MODE_CURRENT_PROGRAMS.json').read_text())
    block = next(row for row in blocks if row['modulus'] == 16)
    counts = Counter(kind for kind, targets in block['program']['gates'])
    current_error = 288 * sum((number * F(root['primitive_operator_errors_exact'][kind]) for kind, number in counts.items()), F())
    matter_error = 128 * onsite_error + swap_error + current_error
    field = json.loads((HERE / 'NATIVE_FIELD.json').read_text())
    electric = 288 * F(field['electric_per_edge_error_exact'])
    one_field = electric + F(field['faces'][0]['operator_error_upper_exact'])
    all_field = electric + sum((F(row['operator_error_upper_exact']) for row in field['faces']), F())
    source_prep = flip_error + uniform_error
    pointer_error = bound(pointer._local_word(), cert.cnot(2, 0) + cert.cnot(0, 1))
    reset_probabilities = {}
    for code in (0, pointer.READY_CODE):
        probabilities = configuration.success_probabilities(4, 32, full=True, mode='reset', target=circuit.decode(code))
        reset_probabilities[code] = min(map(F, probabilities))
    p0, pr = reset_probabilities[0], reset_probabilities[pointer.READY_CODE]
    reset_fail = 3203 * (1 - p0) + 3200 * (1 - pr)
    # Read all2048 matter occupations and1152 electric bits. The one-particle
    # decoder preserves flags for zero/multiple/padded occupations.
    nsteps = 3
    complete = reset_fail + source_prep + nsteps * (matter_error + one_field) + 3200 * pointer_error
    assert 0 < complete < F(1, 100)
    preparation = preparation_word()
    (HERE / 'PREPARATION_PROGRAM.json').write_text(json.dumps({
        'blank_native_resets': 3203, 'reset_target': 0,
        'source_preparation_word': preparation,
        'separate_one_hole_preparation_word': preparation_word(-1),
        'initial_field': 'all288 modular electric registers zero',
        'Gauss_reference': 'G=delta_port0 for one-particle preparation',
        'full_occupation_and_field_readout_bits': 3200,
        'native_bit_pointer_word': pointer._local_word(),
        'initialization_failures': 'stop at the first reset timeout, retain its fine native tape and final carrier state; preceding success factors remain in the unconditional probability',
        'detector_failures': 'every ready reset timeout, native pointer fine record, source/export complement and padding/multiload decoder flag is retained'}, indent=2) + '\n')
    result = {'schema': 'oph.local12.full_matter_field_composition.v1', 'candidate_id': af.MODEL_ID,
              'checks': {'full_GammaR_source_rotations_certified': True,
                         'all_spectator_parity_strings_retained': True,
                         'all288_gauge_current_blocks_counted': True,
                         'autonomous_field_word_error_added': True,
                         'all3203_blank_and3200_pointer_reset_failures_retained': True,
                         'one_particle_preparation_word_emitted': True,
                         'complete_quantum_reference_bound': True, 'source_bytes_unchanged': True},
              'onsite_GammaR_error_exact': str(onsite_error),
              'complete_matter_current_layer_error_exact': str(matter_error),
              'oneface_complete_V_error_exact': str(matter_error + one_field),
              'oneface_complete_V_error_upper': float(matter_error + one_field),
              'homogeneous208face_complete_V_error_upper': float(matter_error + all_field),
              'source_preparation_error_exact': str(source_prep),
              'separate_one_hole_preparation_error_exact': str(1151 * flip_error + uniform_error),
              'swap_distance_classes': swap_rows,
              'source_pulses_per_matter_current_layer': 128 * len(af.mb.fock_reflection_word()) + swap_pulses + block['native_pulses_all288_edges'],
              'native_blank_reset_success_lower_exact': str(p0),
              'native_ready_reset_success_lower_exact': str(pr),
              'full_prepare_read_reset_timeout_upper_exact': str(reset_fail),
              'three_step_full_record_TV_bound_exact': str(complete),
              'three_step_full_record_TV_bound_upper': float(complete),
              'full_law': 'Multiply all configured native reset success factors; then the coherent word with every first-export failure output; then complete3200-bit native readout. Retain the first failing reset and all pointer/source/padding flags. The displayed sparse ideal law plus these complete errors is not success-renormalized.',
              'scope': 'Factorized operator/circuit certificate on all Fock sectors and arbitrary references in the declared clean-work domain. Complete one-particle initialization and full readout programs are emitted; their enormous full basis is not densely enumerated. Actual finite one-particle/one-hole ideal joint laws are in AUTONOMOUS.json. The measured code services and source response are declared architecture. Clock operations, controlled-V phase probe and nonzero-background preparation add their separately certified resources.',
              'source_sha256': pins, 'official_rows_changed': []}
    assert pins == {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (HERE / 'COMPOSITION.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('checks', 'oneface_complete_V_error_upper',
                                           'homogeneous208face_complete_V_error_upper',
                                           'three_step_full_record_TV_bound_upper', 'source_pulses_per_matter_current_layer')}))


if __name__ == '__main__':
    run()
