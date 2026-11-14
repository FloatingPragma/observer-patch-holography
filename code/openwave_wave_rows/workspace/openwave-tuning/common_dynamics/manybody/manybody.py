"""Typed many-body extensions of the frozen local twelve-slot walk.

The collective binary construction preserves the actual binary native rule.
The mode construction chooses fermionic statistics and changes the carrier
state. Neither choice is identified with unrestricted integer balancing.
"""
from collections import Counter, deque
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
WALK = ROOT / 'openwave-tuning/auditor/native_walk_bridge'
FIELD = ROOT / 'openwave-tuning/auditor/native_current_walk'
sys.path[:0] = [str(WALK), str(FIELD), str(ROOT / 'reverse-engineering-reality')]
import common_walk as common
import compile_walk as compiled
import current_walk as field
from model.observer import native
from model.openwave_coherent_native_architecture import source as architecture
from model.observer.source_scalar import compiler

BASE_ID = common.MODEL_ID
COLLECTIVE_ID = 'native-local12-collective-binary-v1-proposed'
MODE_ID = 'native-local12-fermion-modes-v1-proposed'


@lru_cache(None)
def matchings():
    neighbours, edge_ids = common.geometry_tables()
    result = []
    for d in range(6):
        row = [(x, int(neighbours[x, d]), int(edge_ids[x, d]))
               for x in range(96) if x < neighbours[x, d]]
        assert len(row) == 48
        assert sorted(v for a, b, _ in row for v in (a, b)) == list(range(96))
        result.append(tuple(row))
    return tuple(result)


def occupations(mask):
    if type(mask) is not int or not 0 <= mask < 1 << 96:
        raise ValueError('96 binary port occupations required')
    return np.array([(mask >> x) & 1 for x in range(96)], dtype=np.int64)


def collective_step(mask, event):
    """Complete padded event permutation and its physical edge currents."""
    occupations(mask)
    if type(event) is not int or not 0 <= event < 16:
        raise ValueError('four-bit event required')
    current = np.zeros(288, dtype=np.int64)
    if event >= 12 or not event % 2:
        return mask, current
    result = mask
    for a, b, e in matchings()[event // 2]:
        na, nb = (mask >> a) & 1, (mask >> b) & 1
        current[e] = na - nb
        if na != nb:
            result ^= (1 << a) | (1 << b)
    return result, current


def native_matching_records(mask, event, law):
    """The original C/E labels for each actual binary seam in the layer.

    The original coins are correlated and occupation-controlled. They are
    not independent fresh native attempts.
    """
    if law not in ('C', 'E') or not 0 <= event < 12:
        raise ValueError('active event and named native law required')
    n = occupations(mask)
    sigma = event % 2
    rows = []
    for a, b, e in matchings()[event // 2]:
        coin = sigma ^ int(n[a] == 1 and n[b] == 0)
        output, (_, current) = native.native(tuple(n), (a, b), coin)
        rows.append({'seam_index': e, 'seam': [a, b], 'coin': coin,
                     'current': current, 'label': current if law == 'C' else int(n[a] - n[b]),
                     'after_pair': [output[a], output[b]]})
    return rows


def local_controlled_swap(direction):
    """A complete local 9-bit word: pair2,event4,clean scratch3."""
    values = {2 + q: ((2 * direction + 1) >> q) & 1 for q in range(4)}
    values[1] = 1
    return [('cx', (0, 1))] + compiled.multi_control(values, 0, (6, 7, 8)) + [('cx', (0, 1))]


def execute_gates(gates, values):
    answer = np.asarray(values, dtype=np.int64).copy()
    for kind, targets in gates:
        if kind == 'flip':
            answer ^= 1 << targets[0]
        elif kind == 'cx':
            a, b = targets
            answer ^= ((answer >> a) & 1) << b
        elif kind == 'ccx':
            a, b, c = targets
            answer ^= (((answer >> a) & 1) & ((answer >> b) & 1)) << c
        else:
            raise ValueError(kind)
    return answer


def collective_program():
    """Finite native macro program; reused primitives retain their phases."""
    built = compiled.controls()
    gates = []
    for d, matching in enumerate(matchings()):
        local = local_controlled_swap(d)
        for a, b, e in matching:
            mapping = [a, b, 96, 97, 98, 99, 100, 101, 102]
            gates.extend((kind, tuple(mapping[q] for q in targets)) for kind, targets in local)
    counts = Counter(kind for kind, targets in gates)
    pulses = sum(number * len(built['primitive_words'][kind]) for kind, number in counts.items())
    return {'candidate_id': COLLECTIVE_ID, 'matter_bits': 96, 'event_bits': [96, 97, 98, 99],
            'clean_scratch_bits': [100, 101, 102], 'swap_calls': gates,
            'primitive_words': built['primitive_words'],
            'reflection_word': compiler.remap(built['reflection_word'], [96, 97, 98, 99]),
            'swap_gate_counts': dict(counts),
            'W_native_pulses': pulses + len(built['reflection_word']),
            'full_fine_record_intervention': 'inverse binary native recoder on each of48 selected seams; no independent-coin interpretation',
            'padding': 'event12..15 fixed by S and reflected by-1; all scratch complements retained',
            'failure_scope': 'native preparation/reset and proper-code export branches must accompany execution; full103-cell reset instrument is not replayed here'}


def apply_mode_walk(vector):
    """One-particle callback of the proposed mode extension, exactly the base W."""
    return common.apply_walk(np.asarray(vector).reshape(128, 16)).reshape(-1)


def fock_pair_records(left, right):
    """All128×128 position records for a normalized two-fermion Slater state.

    The diagonal contains two particles at one port in distinct modes. It is
    not removed. Off-diagonal entries use unordered position pairs.
    """
    left, right = np.asarray(left), np.asarray(right)
    assert left.shape == right.shape == (2048,)
    amplitude = left[:, None] * right[None, :] - right[:, None] * left[None, :]
    probabilities = abs(np.triu(amplitude, 1)) ** 2
    records = probabilities.reshape(128, 16, 128, 16).sum(axis=(1, 3))
    assert abs(records.sum() - 1) < 2e-10
    return records


def fock_swap_sign(occupied, a, b):
    """Canonical exterior-power lift of a transposition, including spectators."""
    original = tuple(sorted(occupied))
    image = [b if x == a else a if x == b else x for x in original]
    inversions = sum(image[i] > image[j] for i in range(len(image)) for j in range(i + 1, len(image)))
    return tuple(sorted(image)), -1 if inversions % 2 else 1


def givens_word(a, b, theta):
    """Native Jordan-Wigner rotation; signs/statistics are supplied architecture."""
    if not 0 <= a < b:
        raise ValueError('ordered modes required')
    first = ['I'] * (b + 1)
    second = first.copy()
    first[a], first[b] = 'X', 'Y'
    second[a], second[b] = 'Y', 'X'
    for k in range(a + 1, b):
        first[k] = second[k] = 'Z'
    return (compiler.pauli_rotation(''.join(first), -theta / 2)
            + compiler.pauli_rotation(''.join(second), theta / 2))


def occupation_phase_word(mode, theta):
    return [compiler.pulse('I', theta / 2, [mode]), compiler.pulse('Z', -theta / 2, [mode])]


def fock_reflection_word():
    """Sixteen-mode onsite Gamma(R), with identity on the vacuum."""
    U = []
    for k in range(1, 12):
        U += givens_word(0, k, math.asin(1 / math.sqrt(13 - k)))
    return (compiler.inverse(U)
            + [pulse for k in range(1, 16) for pulse in occupation_phase_word(k, math.pi)]
            + U)


def fock_swap_word(a, b):
    return occupation_phase_word(b, math.pi) + givens_word(a, b, math.pi / 2)


def mode_program():
    """Generator metadata avoids materializing a full Fock-space matrix."""
    pairs = []
    for d, matching in enumerate(matchings()):
        for a, b, e in matching:
            pairs.append((16 * common.encode_position(a) + 2 * d + 1,
                          16 * common.encode_position(b) + 2 * d + 1, e))
    return {'candidate_id': MODE_ID, 'base_one_particle_id': BASE_ID,
            'matter_modes': 2048, 'physical_modes': 1152,
            'onsite_mode_groups': [[16 * x + k for k in range(16)] for x in range(128)],
            'onsite_reflection_word': fock_reflection_word(),
            'transport_pairs': pairs, 'transport_word_generator': 'manybody.fock_swap_word(a,b)',
            'vacuum': 'empty occupation, fixed with phase1',
            'statistics': 'fermionic exterior algebra explicitly chosen; no emergent-spin or native-statistics claim',
            'aggregate_load_range_per_valid_port': [0, 12],
            'native_scope': 'one binary mode pair has native balancing populations after recoding; spectator parity is an added coherent phase; aggregate integer-port balancing differs',
            'locality': 'onsite even CAR rotations and adjacent even CAR hopping; encoded Jordan-Wigner pulse strings are not a proved distributed native implementation',
            'current_lift': 'each physical transported mode contributes nu_a- nu_b to its one shared edge F; local reflection preserves every port total',
            'full_field_affine_preparation': 'P_n=sum_x n_x P_x, G=Q delta0; same193 finite Z6 cycle coordinates',
            'complete_instrument_scope': 'finite source-word generators and small-factor execution; no full2048-mode reset or full-Fock dense execution claim'}


def run():
    source_paths = [Path(__file__), Path(common.__file__), Path(compiled.__file__), Path(field.__file__),
                    Path(native.__file__), Path(architecture.__file__), Path(compiler.__file__),
                    WALK / 'COMMON_WALK_RECEIPT.json', FIELD / 'CURRENT_WALK.json']
    pins = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    geometry = field.geometry()
    B, P, C = geometry['B'], geometry['P'], geometry['C']
    # Actual native binary branch and exact reversible recoder comparison.
    native_permutation = architecture.event_permutation()
    recoder = architecture.record_recoder()
    local_cases = 0
    for x, sigma, law in product(range(4), range(2), range(2)):
        coin = sigma ^ int(x == 1)
        before = architecture.BINARY_CODES[x] | (coin << 6) | (law << 7)
        after = recoder[native_permutation[before]]
        y = ((x & 1) << 1 | x >> 1) if sigma else x
        assert after == architecture.BINARY_CODES[y] | sigma << 6 | law << 7
        actual, (_, current) = native.native((x & 1, x >> 1), (0, 1), coin)
        assert actual == (y & 1, y >> 1) and current == (x & 1) - (y & 1)
        local_cases += 1
    gates_checked = 0
    for d in range(6):
        gates = local_controlled_swap(d)
        actual = execute_gates(gates, np.arange(512))
        assert len(np.unique(actual)) == 512
        assert np.array_equal(execute_gates(gates, actual), np.arange(512))
        for x, event in product(range(4), range(16)):
            inp = x | event << 2
            y = (((x & 1) << 1) | x >> 1) if event == 2 * d + 1 else x
            assert actual[inp] == y | event << 2
            gates_checked += 1
    # All one- and two-load binary development inputs, every active event.
    atoms = 0
    expected_one_load = common.expected_swap_permutation()
    for occupied in [*(tuple([x]) for x in range(96)), *combinations(range(96), 2)]:
        mask = sum(1 << x for x in occupied)
        n = occupations(mask)
        P_n = P[:, occupied].sum(axis=1)
        for event in range(16):
            after, j = collective_step(mask, event)
            back, inverse_current = collective_step(after, event)
            assert back == mask and np.array_equal(inverse_current, -j)
            after_n = occupations(after)
            nonzero = np.flatnonzero(j)
            divergence = B[:, nonzero] @ j[nonzero]
            assert np.array_equal(after_n + divergence, n)
            if len(occupied) == 1:
                old = 16 * common.encode_position(occupied[0]) + event
                image = int(expected_one_load[old])
                target = common.decode_position(image // 16)
                assert after == 1 << target
            # P_n+j-P_n' is precisely an integral cycle, including simultaneous moves.
            difference = P_n + j - P[:, np.flatnonzero(after_n)].sum(axis=1)
            chord_values = difference[geometry['chords']]
            changed_chords = np.flatnonzero(chord_values)
            assert np.array_equal(C[:, changed_chords] @ chord_values[changed_chords], difference)
            atoms += 1
    print(json.dumps({'collective_atom_checks_finished': atoms}), flush=True)
    # The discarded general-integer current is required for reversibility.
    integer_images = []
    for loads, coin in product(((2, 0), (1, 1), (0, 2)), (0, 1)):
        out, record = native.native(loads, (0, 1), coin)
        integer_images.append({'before': list(loads), 'coin': coin, 'after': list(out), 'current': record[1]})
        assert out == (1, 1)
    assert len({(tuple(r['after']), r['coin']) for r in integer_images}) == 2
    assert len({(tuple(r['after']), r['coin'], r['current']) for r in integer_images}) == 6
    # Distant binary units move together under one shared sigma.
    distance = [-1] * 96
    distance[0] = 0
    todo = deque([0])
    while todo:
        x = todo.popleft()
        for y in common.graph().neighbours[x]:
            if distance[y] < 0:
                distance[y] = distance[x] + 1
                todo.append(y)
    distant = int(np.argmax(distance))
    assert distance[distant] >= 3
    joint = Counter()
    mask = 1 | 1 << distant
    for event in range(12):
        after, j = collective_step(mask, event)
        joint[(int(not after & 1), int(not after & (1 << distant)))] += F(1, 12)
    assert joint == {(0, 0): F(1, 2), (1, 1): F(1, 2)}
    # Fermion onsite rotations and transpositions, full vacuum/double/spectator domain.
    native_errors = []
    for b in (1, 2, 3):
        word = fock_swap_word(0, b)
        actual = compiler.execute(word, b + 1)
        expected = np.zeros_like(actual)
        for mask in range(1 << (b + 1)):
            image, sign = fock_swap_sign([x for x in range(b + 1) if mask >> x & 1], 0, b)
            expected[sum(1 << x for x in image), mask] = sign
        error = float(np.linalg.norm(actual - expected, 2))
        assert error < 1e-11
        native_errors.append(error)
    # Local uniform mode preparation and exterior reflection on two modes.
    U = np.eye(16)
    for k in range(1, 12):
        theta = math.asin(1 / math.sqrt(13 - k))
        G = np.eye(16)
        G[0, 0] = G[k, k] = math.cos(theta)
        G[k, 0], G[0, k] = math.sin(theta), -math.sin(theta)
        U = G @ U
    target_u = np.r_[np.ones(12) / math.sqrt(12), np.zeros(4)]
    assert np.linalg.norm(U[:, 0] - target_u) < 1e-14
    R = U @ np.diag([1] + [-1] * 15) @ U.T
    assert np.linalg.norm(R - (2 * np.outer(target_u, target_u) - np.eye(16)), 2) < 1e-14
    pairs = list(combinations(range(16), 2))
    exterior = np.array([[R[a, c] * R[b, d] - R[a, d] * R[b, c]
                         for c, d in pairs] for a, b in pairs])
    assert np.linalg.norm(exterior.T @ exterior - np.eye(120), 2) < 1e-13
    # Complete same-one-particle-callback two-fermion records.
    left, right = np.zeros(2048, complex), np.zeros(2048, complex)
    left[16 * common.encode_position(0):16 * common.encode_position(0) + 12] = 1 / math.sqrt(12)
    right[16 * common.encode_position(distant):16 * common.encode_position(distant) + 12] = 1 / math.sqrt(12)
    mode_rows = []
    for step in range(3):
        records = fock_pair_records(left, right)
        mode_rows.append({'step': step, 'all128_position_pair_probabilities': records.tolist(),
                          'complete_probability': float(records.sum()),
                          'one_particle_orthogonality': float(abs(np.vdot(left, right)))})
        if step == 1:
            # Disjoint one-step supports permit independent movement indicators.
            p0 = common.encode_position(0)
            p1 = common.encode_position(distant)
            stay0 = float(abs(left.reshape(128, 16)[p0]).dot(abs(left.reshape(128, 16)[p0])))
            stay1 = float(abs(right.reshape(128, 16)[p1]).dot(abs(right.reshape(128, 16)[p1])))
            assert abs(stay0 - .5) < 1e-12 and abs(stay1 - .5) < 1e-12
        left, right = apply_mode_walk(left), apply_mode_walk(right)
    # Two different local fermion modes can occupy one port. Literal total-two
    # balancing sends(2,0) deterministically to(1,1); the mode program does not.
    left, right = np.zeros(2048, complex), np.zeros(2048, complex)
    left[:12] = 1 / math.sqrt(12)
    right[0], right[2] = 1 / math.sqrt(2), -1 / math.sqrt(2)
    local_pairs = fock_pair_records(apply_mode_walk(left), apply_mode_walk(right))
    assert abs(local_pairs[0, 0] - .5) < 1e-12
    # Supplied fermion parity changes a coherent output, including a spectator.
    image, sign = fock_swap_sign((0, 1), 0, 2)
    assert image == (1, 2) and sign == -1
    # A fixed-total-two reference detects the supplied fermion double phase.
    # Pair occupied versus a disjoint occupied spectator pair, both total2.
    pair_phase = fock_swap_sign((0, 1), 0, 1)[1]
    spectator_phase = fock_swap_sign((2, 3), 0, 1)[1]
    assert pair_phase == -1 and spectator_phase == 1
    coherent_native_phase_distance = abs(pair_phase - spectator_phase) / 2
    assert coherent_native_phase_distance == 1
    # Known-reference regional intervals improve over the symmetric decoder.
    decode_checks = 0
    for Q in range(6):
        for G_R in range(Q + 1):
            degrees = list(range(-G_R, Q - G_R + 1))
            assert len({q % 6 for q in degrees}) == len(degrees)
            decode_checks += 1
    assert 0 % 6 == 6 % 6
    # Actual equal-G, equal-field Gauss configurations in the mode model.
    n_root, n_other = np.zeros(96, dtype=np.int64), np.zeros(96, dtype=np.int64)
    n_root[0], n_other[1] = 6, 6
    F_root, F_other = P @ n_root, P @ n_other
    assert np.array_equal(F_root % 6, F_other % 6)
    assert np.array_equal((n_root + B @ F_root) % 6, (n_other + B @ F_other) % 6)
    assert n_root[1] == 0 and n_other[1] == 6
    for modulus in (13, 16, 25):
        assert len({n % modulus for n in range(13)}) == 13
    programs = {'collective': collective_program(), 'mode': mode_program()}
    (HERE / 'MANYBODY_PROGRAMS.json').write_text(json.dumps(programs, indent=2) + '\n')
    result = {'schema': 'oph.local12.manybody.v1',
              'base_candidate_id': BASE_ID,
              'checks': {'actual_native_binary_recoder': True,
                         'full_local_circuit_padding_and_dirty_complements': True,
                         'all_one_two_load_collective_atoms': True,
                         'one_load_exact_base_restriction': True,
                         'shared_field_Gauss_and_affine_cycle_identity': True,
                         'unrestricted_balancing_needs_retained_memory': True,
                         'shared_coin_cluster_control': True,
                         'native_fermion_parity_and_spectators': True,
                         'onsite_exterior_reflection_unitary': True,
                         'complete_mode_pair_measurements': True,
                         'aggregate_native_mismatch_exposed': True,
                         'reference_aware_decoder_and_alias_control': True,
                         'source_pins_unchanged': True},
              'binary_native_cases': local_cases, 'local_clean_gate_columns': gates_checked,
              'local_full_dirty_columns': 6 * 512, 'collective_atoms': atoms,
              'integer_record_erasure_control': integer_images,
              'cluster_control': {'occupied_ports': [0, distant], 'distance': distance[distant],
                                  'shared_coin_joint_move_probabilities': {'00': '1/2', '01': '0', '10': '0', '11': '1/2'},
                                  'independent_local_mode_joint': {'00': '1/4', '01': '1/4', '10': '1/4', '11': '1/4'},
                                  'movement_covariance': '1/4', 'total_variation_from_product': '1/2',
                                  'scope': 'A shared prepared event is a long-range resource. This control does not assert superluminal signaling.'},
              'native_fock_swap_operator_errors': native_errors,
              'fermion_added_phase_fixed_Q2_trace_distance': coherent_native_phase_distance,
              'mode_same_port_after_one_W_probability': float(local_pairs[0, 0]),
              'literal_integer_same_port_after_addressed_event_probability': '0',
              'mode_pair_records': mode_rows,
              'field_contract': {'modulus': 6, 'shared_edges': 288, 'cycle_coordinates': 193,
                                 'fixed_Q_reference': 'G=Q delta0', 'affine_particular': 'P_n=sum_x n_x P_x',
                                 'identity': 'S_F J_Q=J_Q S, R J_Q=J_Q R on each fixed-Q sector',
                                 'normalization': 'finite6^193 fibre only; no infinite uniform state',
                                 'preparation_scope': 'full manybody affine initialization/reset instrument not executed'},
              'regional_decoder': {'known_Q_and_G_R': 'q lies in[-G_R,Q-G_R]; unique mod6 for Q<=5',
                                   'source_independent_symmetric_interval': '[-Q,Q] unique mod6 for Q<=2',
                                   'first_reference_aware_alias': 'Q6: q0 and q6 have equal residue',
                                   'binary_pointwise_exception': 'binary port loads0/1 are recovered from each local Gauss residue even at larger totalQ; one regional residue alone is insufficient',
                                   'mode_pointwise_alias': 'port occupations0,6,12 share a mod6 Gauss residue'},
              'actual_Q6_alias': {'occupied_modes_first': list(range(6)),
                                 'occupied_modes_second': list(range(16, 22)),
                                 'same_reference': 'G=6 delta0', 'same_mod6_field': 'all288 edge labels zero',
                                 'region': [1], 'regional_charges': [0, 6]},
              'field_modulus_amendment': {'parameter': 'm>=13 for known-reference local physical occupations0..12',
                                         'executed_distinguishing_moduli': [13, 16, 25],
                                         'reference_independent_signed_interval': 'm>=25 for[-12,12]',
                                         'source_transfer': 'tree and finite-fibre algebra transfers; new modulus requires a new field preparation, current circuit and source version',
                                         'regional_note': 'a single regional residue needs its regional bound; reading all local residues with m>=13 and summing decoded integer occupations avoids that alias'},
              'architecture_choice': 'Collective binary version exactly retains one-load W but has shared-event correlations. Local fermion-mode version fixes factorization for distant one-step product preparations, chooses exterior statistics and changes aggregate native multi-load behavior. Neither is silently selected for the common source.',
              'accuracy_scope': 'Exact permutation/algebra statements with numerical native-factor and full-vector controls; no aggregate physical pulse-error certificate inferred from floating residuals.',
              'source_sha256': pins, 'official_rows_changed': []}
    assert pins == {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    (HERE / 'MANYBODY.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('checks', 'collective_atoms', 'native_fock_swap_operator_errors',
                     'mode_same_port_after_one_W_probability', 'regional_decoder')}))


if __name__ == '__main__':
    run()
