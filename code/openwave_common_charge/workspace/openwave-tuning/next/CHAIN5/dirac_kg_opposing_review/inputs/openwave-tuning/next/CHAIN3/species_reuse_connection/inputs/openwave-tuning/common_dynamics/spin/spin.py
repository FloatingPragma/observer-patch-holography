"""Geometric rotations and native exchange records of the shared local walk."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path[:0] = [str(ROOT/'reverse-engineering-reality'),
               str(ROOT/'openwave-tuning/auditor/native_walk_bridge'),
               str(ROOT/'openwave-tuning/auditor/common_walk_particle'),
               str(HERE.parent/'manybody'), str(HERE.parent/'continuum')]
import common_walk as common
import compile_walk as compiled
import manybody
import continuum
from native_phase_probe import control_word, actual_native_pointer
from model.openwave_coherent_native_architecture import source as architecture
from model.observer.source_scalar import compiler
from model.reference_recovery import recovery


def rotations():
    geometry = common.spatial
    carrier = geometry.slow_geometry.source()
    frame = geometry.frame()
    zero = carrier.q5(0)
    rows = []
    for shift in range(3):
        for signs in product((-1, 1), repeat=3):
            if math.prod(signs) != 1:
                continue
            O = np.zeros((3, 3), dtype=np.int64)
            for a in range(3):
                O[a, (a+shift) % 3] = signs[a]
            perm = []
            for v in frame:
                image = tuple(v[(a+shift) % 3] if signs[a] == 1
                              else carrier.q5_sub(zero, v[(a+shift) % 3]) for a in range(3))
                perm.append(frame.index(image))
            rows.append((O, tuple(perm)))
    assert len(rows) == 12
    return rows


def geometric_permutation(O, port_permutation):
    g = common.graph()
    coordinates = {c: i for i, c in enumerate(g.cell_coordinates)}
    cells = [coordinates[tuple(map(int, O @ np.array(c) % 2))] for c in g.cell_coordinates]
    vertices = [12*cells[x//12] + port_permutation[x % 12] for x in range(96)]
    neighbours, _ = common.geometry_tables()
    result = np.empty(2048, dtype=np.int64)
    for position in range(128):
        cell, port = divmod(position, 16)
        destination = 16*cells[cell] + (port_permutation[port] if port < 12 else port)
        for event in range(16):
            target_event = event
            if port < 12 and event < 12:
                x = 12*cell+port
                direction, sigma = divmod(event, 2)
                image_neighbour = vertices[int(neighbours[x, direction])]
                image_direction = list(neighbours[vertices[x]]).index(image_neighbour)
                target_event = 2*image_direction+sigma
            result[16*position+event] = 16*destination+target_event
    assert sorted(result) == list(range(2048))
    return result


def exterior_action(occupied, permutation):
    image = [int(permutation[x]) for x in sorted(occupied)]
    inversions = sum(image[a] > image[b] for a in range(len(image)) for b in range(a+1, len(image)))
    return tuple(sorted(image)), (-1)**inversions


def native_controlled_records(word, system, system_bits):
    control = system_bits
    H = compiled.controls()['primitive_words']['hadamard']
    program = (compiler.remap(H, [control]) + control_word(word, control)
               + compiler.remap(H, [control]))
    initial = np.zeros((2**(system_bits+1), 1), complex)
    initial[:2**system_bits, 0] = system
    run = recovery.execute_complete(program, system_bits+1, initial)
    state = run.success[:, 0].reshape(2, 2**system_bits)
    rho = state @ state.conj().T
    outputs, pointer = actual_native_pointer()
    laws = {}
    for law, outcomes in outputs.items():
        fine, bits = [], [0., 0.]
        for meta, kraus in outcomes:
            p = sum(float(np.trace(K @ rho @ K.conj().T).real) for K in kraus)
            fine.append({'record': meta, 'probability': p})
            bits[meta['bit']] += p
        laws[law] = {'fine_records': fine, 'bits': bits}
    return {'native_pulses_before_pointer': len(program), 'pointer_native_pulses': 6,
            'laws': laws, 'export_failure_trace': float(np.trace(run.failure_gram).real),
            'pointer_export_failure_trace': float(np.trace(pointer.failure_gram).real),
            'reset_scope': 'fresh source/probe/pointer initialization needs the common native reset contract; none is conditioned away in a complete implementation'}


def run():
    paths = [Path(__file__), Path(common.__file__), Path(compiled.__file__),
             Path(manybody.__file__), Path(continuum.__file__),
             HERE.parent/'manybody/MANIFEST.json', HERE.parent/'continuum/MANIFEST.json']
    pins = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    group = rotations()
    permutations = [geometric_permutation(*row) for row in group]
    lookup = {tuple(O.reshape(-1)): i for i, (O, _) in enumerate(group)}
    swap = common.expected_swap_permutation()
    for i, (O, _) in enumerate(group):
        p = permutations[i]
        assert np.array_equal(p[swap], swap[p])
        # The uniform-event reflection is invariant under the complete
        # within-source permutation of admitted events; padding stays typed.
        assert all(sorted(int(p[16*x+e] % 16) for e in range(12)) == list(range(12)) for x in range(128))
        assert all(p[16*x+e]//16 == p[16*x]//16 for x in range(128) for e in range(16))
        for j, (Q, _) in enumerate(group):
            k = lookup[tuple((O @ Q).reshape(-1))]
            assert np.array_equal(permutations[i][permutations[j]], permutations[k])

    data = continuum.rational_bloch_data()
    e = data['scalar']
    B = np.column_stack([e]+[np.sqrt(48)*v for v in data['vectors']])
    assert np.max(abs(B.T @ B - np.eye(4))) < 1e-13
    covariance = []
    for (O, _), p in zip(group, permutations):
        small = np.array([12*(p[16*a+event]//16)+p[16*a+event] % 16 for a in range(12) for event in range(12)])
        assert sorted(small) == list(range(144))
        GB = np.empty_like(B)
        GB[small] = B
        target = np.eye(4)
        target[1:, 1:] = O
        error = float(np.max(abs(GB-B @ target)))
        assert error < 1e-13
        covariance.append(error)

    # Native compiler for two actual geometric generators. The local8-bit
    # program retains every dirty scratch input; cell-bit permutations are
    # separately exact reversible wiring operations.
    targets = [np.diag([-1, -1, 1]), np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])]
    programs, cycles = [], []
    primitives = compiled.controls()['primitive_words']
    for O in targets:
        index = lookup[tuple(O.reshape(-1))]
        p = permutations[index]
        local = list(map(int, p[:256]))
        assert max(local) < 256
        program = architecture.compile_permutation(local)
        execution = architecture.execute_basis(program)
        assert np.array_equal(execution[:256], local)
        assert len(np.unique(execution)) == 8192
        order = 2 if np.array_equal(O, targets[0]) else 3
        composed = np.arange(2048)
        for _ in range(order):
            composed = p[composed]
        assert np.array_equal(composed, np.arange(2048))
        phase_rows = []
        for occupied in ((0,), (0, 17), (0, 17, 83), (1, 7, 31, 143)):
            state, phase = occupied, 1
            for _ in range(order):
                state, sign = exterior_action(state, p)
                phase *= sign
            assert state == tuple(sorted(occupied)) and phase == 1
            phase_rows.append({'occupied': occupied, 'return_phase': phase})
        counts = Counter(kind for kind, _ in program['gates'])
        programs.append({'rotation': O.tolist(), 'order': order, 'local_program': program,
                         'local_native_pulses': sum(n*len(primitives[k]) for k,n in counts.items()),
                         'source_data_bits': list(range(11)), 'fresh_scratch_bits': list(range(23,28)),
                         'local_scratch_remap': {str(8+j): 23+j for j in range(5)},
                         'cell_action': 'the stated3-by3 signed permutation modulo2; cyclic coordinates use two swaps, six native CNOT words; sign flips are identity atsize2',
                         'native_cell_CNOT_count': 0 if order == 2 else 6})
        cycles.append({'rotation_order': order, 'all2048_one_particle_return_columns': True, 'exterior_controls': phase_rows})

    # Continuous extension on the actual scalar/vector acoustic coordinates.
    J = []
    for axis in range(3):
        generator = np.zeros((4,4), complex)
        for a in range(3):
            for b in range(3):
                epsilon = 0 if len({axis,a,b}) < 3 else (1 if (axis,a,b) in ((0,1,2),(1,2,0),(2,0,1)) else -1)
                generator[a+1,b+1] = -1j*epsilon
        J.append(generator)
    assert max(np.max(abs(J[a] @ J[b]-J[b] @ J[a]-1j*J[c])) for a,b,c in ((0,1,2),(1,2,0),(2,0,1))) < 1e-14
    casimir = sum(x @ x for x in J)
    assert np.array_equal(casimir, np.diag([0,2,2,2]))
    assert np.max(abs(expm(-2j*np.pi*J[2])-np.eye(4))) < 1e-14
    alpha = [np.sqrt(48)*B.T @ H @ B for H in data['H']]
    assert max(np.max(abs(a @ a - np.eye(4))) for a in alpha) > .99
    mass_R = B.T @ data['R'] @ B
    assert np.max(abs(mass_R-np.diag([1,-1,-1,-1]))) < 1e-13
    assert max(np.max(abs(mass_R @ j-j @ mass_R)) for j in J) < 1e-13
    for sign in (-1,1):
        longitudinal = np.array([1,0,0,sign],complex)/np.sqrt(2)
        assert np.linalg.norm(J[2] @ longitudinal) < 1e-14

    exchange = native_controlled_records(manybody.fock_swap_word(0,1), np.array([0,0,0,1],complex), 2)
    double_exchange = native_controlled_records(manybody.fock_swap_word(0,1)*2, np.array([0,0,0,1],complex), 2)
    for law in ('C','E'):
        assert abs(exchange['laws'][law]['bits'][1]-1) < 2e-11
        assert abs(double_exchange['laws'][law]['bits'][0]-1) < 2e-11
    result = {'schema':'oph.local12.rotation_exchange.v1', 'base_candidate':common.MODEL_ID,
              'checks':{'actual_tetrahedral_port_cell_rotations':True,'complete_padded_group_representation':True,
                        'same_walk_and_mass_coin_covariance':True,'native_geometric_generator_circuits':True,
                        'scalar_vector_acoustic_representation':True,'integer_Casimir_and_helicity':True,
                        'Clifford_mutation_rejected':True,'native_fermion_exchange_readout':True,
                        'exchange_is_not_geometric_half_spin':True,'source_bytes_unchanged':True},
              'group_order':12,'group_product_checks':144,'complete_source_dimension':2048,
              'native_local_clean_columns':512,'native_local_dirty_columns':16384,
              'maximum_acoustic_rotation_error':max(covariance),'rotation_cycles':cycles,
              'continuum_rotation':'scalar direct sum ordinaryvector; physicalhelicity0 on longitudinal propagating branch',
              'Casimir_eigenvalues':[0,2,2,2], 'two_pi_geometric_action':'identity',
              'native_controlled_exchange':exchange,'native_controlled_double_exchange':double_exchange,
              'scope':'The geometric action is single-valued and its exterior lift returns+1 on these complete rotation cycles. The CAR minus sign is supplied exchange statistics. Nohalfspin, intrinsic magneticmoment or g2 follows from this source. Full rotation apparatus resets and accumulated pulse errors use the common source contracts and are not inferred from numerical residuals.',
              'source_sha256':pins,'official_rows_changed':[]}
    assert pins == {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (HERE/'ROTATION_PROGRAM.json').write_text(json.dumps({'primitive_words':primitives,'generators':programs},indent=2)+'\n')
    (HERE/'ROTATION_EXCHANGE.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('checks','group_order','maximum_acoustic_rotation_error','Casimir_eigenvalues')}))


if __name__ == '__main__':
    run()
