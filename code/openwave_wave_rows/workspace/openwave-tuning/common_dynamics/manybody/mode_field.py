"""Shared modular current transport for the explicit fermion-mode amendment."""
from collections import Counter
from itertools import product
from pathlib import Path
import hashlib
import json

import numpy as np

import manybody as mb

HERE = Path(__file__).resolve().parent


def port_loads(occupied):
    if len(set(occupied)) != len(occupied) or any(type(x) is not int or not 0 <= x < 2048 for x in occupied):
        raise ValueError('distinct padded matter modes required')
    n = np.zeros(96, dtype=np.int64)
    for mode in occupied:
        address, event = divmod(mode, 16)
        vertex = mb.common.decode_position(address)
        if vertex is not None and event < 12:
            n[vertex] += 1
    return n


def mode_transport(occupied, F, modulus):
    """Full padded matter permutation, CAR phase and one shared field update."""
    if type(modulus) is not int or modulus < 2:
        raise ValueError('finite declared modulus required')
    before = tuple(sorted(occupied))
    port_loads(before)
    F = np.asarray(F, dtype=np.int64)
    if F.shape != (288,) or np.any(F < 0) or np.any(F >= modulus):
        raise ValueError('288 admitted modular field labels required')
    index = mb.common.expected_swap_permutation()
    image = [int(index[x]) for x in before]
    sign = (-1) ** sum(image[i] > image[j] for i in range(len(image)) for j in range(i + 1, len(image)))
    occupied_set = set(before)
    current = np.zeros(288, dtype=np.int64)
    for d, matching in enumerate(mb.matchings()):
        for a, b, edge in matching:
            left = 16 * mb.common.encode_position(a) + 2 * d + 1
            right = 16 * mb.common.encode_position(b) + 2 * d + 1
            current[edge] = int(left in occupied_set) - int(right in occupied_set)
    return tuple(sorted(image)), (F + current) % modulus, sign, current


def field_block(modulus):
    """Finite current-addition native circuit; pair occupations are controls.

    Field padding is an explicit fixed complement. The matter swap and its
    CAR phase follow this block in the full layer.
    """
    bits = (modulus - 1).bit_length()
    size = 1 << (bits + 2)
    permutation = []
    mask = (1 << bits) - 1
    for index in range(size):
        f = index & mask
        a, b = (index >> bits) & 1, (index >> (bits + 1)) & 1
        output = (f + a - b) % modulus if f < modulus else f
        permutation.append((index & ~mask) | output)
    if modulus == 1 << bits:
        # Descending targets leave every lower control bit at its old value.
        # Increment for10 and decrement for01; equal occupations do nothing.
        scratch = tuple(range(bits + 2, bits + 2 + max(0, bits - 1)))
        gates = []
        for a, b, lower in ((1, 0, 1), (0, 1, 0)):
            for target in reversed(range(bits)):
                controls = {bits: a, bits + 1: b}
                controls.update({q: lower for q in range(target)})
                gates += mb.compiled.multi_control(controls, target, scratch)
        program = {'address_bits': bits + 2, 'scratch_bits': len(scratch), 'gates': gates}
        strategy = 'controlled binary increment/decrement with descending targets'
    else:
        program = mb.architecture.compile_permutation(tuple(permutation))
        strategy = 'complete finite permutation compiler with fixed field padding'
    actual = mb.architecture.execute_basis(program)
    assert np.array_equal(actual[:size], permutation)
    assert len(np.unique(actual)) == len(actual)
    gates = Counter(kind for kind, targets in program['gates'])
    words = mb.architecture.primitive_routines()
    pulses = sum(number * len(words[kind]) for kind, number in gates.items())
    return {'modulus': modulus, 'field_bits': bits, 'pair_control_bits': [bits, bits + 1],
            'program': program, 'primitive_words': words, 'compiler_strategy': strategy,
            'clean_columns': size, 'all_dirty_columns': len(actual),
            'native_pulses_per_edge': pulses, 'native_pulses_all288_edges': 288 * pulses,
            'padding_labels': list(range(modulus, 1 << bits))}


def run():
    paths = [Path(__file__), Path(mb.__file__), Path(mb.field.__file__), Path(mb.common.__file__),
             Path(mb.architecture.__file__), Path(mb.compiler.__file__)]
    pins = {str(p.relative_to(mb.ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    g = mb.field.geometry()
    B, P, C = g['B'], g['P'], g['C']
    rng = np.random.default_rng(3106)
    cases = []
    # Every physical seam's empty, either-single and double mode occupations.
    for d, matching in enumerate(mb.matchings()):
        for a, b, edge in matching:
            left = 16 * mb.common.encode_position(a) + 2 * d + 1
            right = 16 * mb.common.encode_position(b) + 2 * d + 1
            cases.extend(((), (left,), (right,), (left, right)))
    cases.extend(tuple(sorted(map(int, rng.choice(2048, count, replace=False))))
                 for count in (2, 3, 6, 12) for _ in range(16))
    checked = 0
    for modulus in (6, 13, 16):
        F = rng.integers(0, modulus, size=288, dtype=np.int64)
        for occupied in cases:
            n = port_loads(occupied)
            after, F_after, sign, j = mode_transport(occupied, F, modulus)
            returned, F_returned, reverse_sign, reverse_j = mode_transport(after, F_after, modulus)
            n_after = port_loads(after)
            assert returned == occupied and np.array_equal(F_returned, F)
            assert sign * reverse_sign == 1 and np.array_equal(reverse_j, -j)
            assert np.array_equal(n_after + B @ j, n)
            assert np.array_equal((n_after + B @ F_after) % modulus, (n + B @ F) % modulus)
            particular = P @ n + j - P @ n_after
            assert np.array_equal(C @ particular[g['chords']], particular)
            checked += 1
    # A lost native current breaks the same-source Gauss identity.
    before = (1,)
    after, _, _, j = mode_transport(before, np.zeros(288, dtype=np.int64), 13)
    changed = np.flatnonzero(j)
    assert len(changed) == 1
    wrong = j.copy()
    wrong[changed[0]] = 0
    assert np.any(port_loads(after) + B @ wrong != port_loads(before))
    blocks = [field_block(modulus) for modulus in (6, 13, 16)]
    (HERE / 'MODE_CURRENT_PROGRAMS.json').write_text(json.dumps(blocks, indent=2) + '\n')
    result = {'schema': 'oph.local12.fermion_mode_field.v1', 'candidate_id': mb.MODE_ID,
              'checks': {'complete_mode_transport_involution': True,
                         'fermion_phase_and_all_spectator_signs_retained': True,
                         'shared_native_current_and_Gauss': True,
                         'integer_affine_cycle_identity': True,
                         'modular_field_block_full_dirty_extension': True,
                         'dropped_current_mutation_rejected': True,
                         'source_pins_unchanged': True},
              'full_matter_field_cases': checked,
              'moduli': [6, 13, 16], 'known_reference_minimum_local_decoder_modulus': 13,
              'native_field_blocks': [{k: v for k, v in block.items() if k not in ('program', 'primitive_words')}
                                      for block in blocks],
              'scope': 'Finite current word generators and exact basis checks. Fermionic statistics and mode carrier are declared amendments. Full matter-field preparation/reset, aggregate pulse-error bound and autonomous field energy are not executed here.',
              'shared_state_api': 'mode_transport(occupied_modes,field288,modulus) -> modes,field,CAR_sign,current288',
              'source_sha256': pins, 'official_rows_changed': []}
    assert pins == {str(p.relative_to(mb.ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (HERE / 'MODE_FIELD.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('checks', 'full_matter_field_cases', 'native_field_blocks')}))


if __name__ == '__main__':
    run()
