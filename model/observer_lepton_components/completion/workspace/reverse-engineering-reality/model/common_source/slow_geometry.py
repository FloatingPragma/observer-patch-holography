"""Attach the existing rank-three native response carrier to joint records.

Uniform all-seam calibration and the actual two-seam feedback experiment are
distinct protocols on the SAME native primitive and port chart. Only the
former has mean I-L/60. All geometric coordinates below are dimensionless
response readouts; no physical position, momentum or clock is inferred.
"""
from __future__ import annotations

from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import importlib
import json
from pathlib import Path
import re
import sys

import numpy as np

from ..observer import native
from . import native_join


ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = ROOT/'evidence/observer_dynamics_20260925/oph-physics-sim'
SPEC = ROOT/'evidence/native_geometric_source_20260925/spec.json'
SOURCE_TO_RER = (8, 10, 1, 3, 7, 11, 0, 4, 6, 9, 2, 5)
THIRD_SEAM = (0, 7)


@lru_cache(maxsize=1)
def source():
    """Reuse the existing exact carrier, checking its retained public bytes."""
    spec = json.loads(SPEC.read_text(encoding='utf-8'))
    for name in ('oph_exact/carrier.py',
                 'oph_fpe/dynamics/self_readback_repair_closure.py'):
        if hashlib.sha256((ARCHIVE/name).read_bytes()).hexdigest() != spec['source_sha256'][name]:
            raise ValueError('the pinned native geometry source has changed')
    location = str(ARCHIVE)
    sys.path.insert(0, location)
    try:
        module = importlib.import_module('oph_exact.carrier')
    finally:
        sys.path.remove(location)
    if Path(module.__file__).resolve() != (ARCHIVE/'oph_exact/carrier.py').resolve():
        raise ValueError('a foreign native carrier module is already loaded')
    closure = sys.modules['oph_fpe.dynamics.self_readback_repair_closure']
    if Path(closure.__file__).resolve() != (ARCHIVE/'oph_fpe/dynamics/self_readback_repair_closure.py').resolve():
        raise ValueError('a foreign native incidence module is already loaded')
    return module


def raw_frame():
    """The existing source chart, exactly over Q(sqrt5), without its radius.

    This coordinate witness is the one in PrimitivePortFrameQuotient and the
    pinned source icosahedron. The invariant geometry is its checked Gram.
    """
    c = source()
    zero, one, minus, phi = c.q5(0), c.q5(1), c.q5(-1), c.GOLDEN
    negative_phi = c.q5_sub(zero, phi)
    return ((minus, phi, zero), (one, phi, zero),
            (minus, negative_phi, zero), (one, negative_phi, zero),
            (zero, minus, phi), (zero, one, phi),
            (zero, minus, negative_phi), (zero, one, negative_phi),
            (phi, zero, minus), (phi, zero, one),
            (negative_phi, zero, minus), (negative_phi, zero, one))


def _add(x, y):
    return tuple(source().q5_add(a, b) for a, b in zip(x, y, strict=True))


def _subtract(x, y):
    return tuple(source().q5_sub(a, b) for a, b in zip(x, y, strict=True))


def _scale(value, x):
    return tuple(source().q5_mul(source().q5(value), a) for a in x)


def _dot(x, y):
    c = source()
    result = c.q5(0)
    for a, b in zip(x, y, strict=True):
        result = c.q5_add(result, c.q5_mul(a, b))
    return result


def _encode_vector(vector):
    return [source().q5_str(value) for value in vector]


def _determinant(rows):
    c = source()
    a, b, d = rows
    cross = tuple(c.q5_sub(c.q5_mul(b[(i+1) % 3], d[(i+2) % 3]),
                           c.q5_mul(b[(i+2) % 3], d[(i+1) % 3])) for i in range(3))
    return _dot(a, cross)


def load_readout(loads):
    """External working-state audit; the controller does not read these loads."""
    if (not isinstance(loads, (tuple, list)) or len(loads) != 12
            or any(type(n) is not int or n < 0 for n in loads) or sum(loads) != 2):
        raise ValueError('a native total-two integer load is required')
    result = (source().q5(0),)*3
    for count, direction in zip(loads, raw_frame(), strict=True):
        result = _add(result, _scale(count, direction))
    return result


def current_displacement(seam, current):
    """Readable from the admitted seam and committed current alone."""
    if (not isinstance(seam, (tuple, list)) or len(seam) != 2
            or any(type(i) is not int for i in seam)
            or tuple(sorted(seam)) not in source().seams()):
        raise ValueError('a seam in the native source chart is required')
    if type(current) is not int or current not in (-1, 0, 1):
        raise ValueError('a total-two native current is required')
    a, b = seam
    return _scale(current, _subtract(raw_frame()[b], raw_frame()[a]))


def registered_binding():
    """Check all 144 Gram entries against the actual registered Lean chart."""
    c = source()
    quotient = (ROOT/'Lean/Screen/PrimitivePortFrameQuotient.lean').read_text(encoding='utf-8')
    match = re.search(r'def sourceToRERPort[^=]*:=\s*!\[([^]]+)\]', quotient)
    if match is None or tuple(map(int, match.group(1).split(','))) != SOURCE_TO_RER:
        raise ValueError('registered source-to-RER relabel changed')
    lean = (ROOT/'Lean/Screen/PortFrameGram.lean').read_text(encoding='utf-8')
    block = lean.split('def neighbors :', 1)[1].split('/-- Adjacency relation', 1)[0]
    neighbors = {int(i): tuple(map(int, values.split(',')))
                 for i, values in re.findall(r'\|\s*(\d+)\s*=>\s*\[([^]]+)\]', block)}
    if set(neighbors) != set(range(12)):
        raise ValueError('registered incidence table is incomplete')
    gram = c.intrinsic_gram_exact()
    for i, registered_i in enumerate(SOURCE_TO_RER):
        for j, registered_j in enumerate(SOURCE_TO_RER):
            expected = (c.q5(1) if i == j else c.q5(0, F(1, 5))
                        if registered_j in neighbors[registered_i] else c.q5(-1)
                        if registered_j == 11-registered_i else c.q5(0, F(-1, 5)))
            if gram[i][j] != expected:
                raise ValueError('native source Gram differs from the registered Lean chart')
    return {'source_to_registered_lean_port': list(SOURCE_TO_RER),
            'registered_gram_entries_checked_exactly': 144,
            'coordinate_witness_theorem': 'OPH.PrimitivePortFrameQuotient.rawAxis_dot_eq_sourceGram',
            'all_pass': True}


@lru_cache(maxsize=1)
def _calibration():
    c = source()
    seams, states = c.seams(), native.states()
    lap = c.laplacian()
    mean = tuple(tuple(F(int(i == j))-F(int(lap[i, j]), 60)
                       for j in range(12)) for i in range(12))
    coordinates = {state: load_readout(state) for state in states}
    count = 0
    for state in states:
        total = [0]*12
        for seam in seams:
            for coin in (0, 1):
                final, record = native.native(state, seam, coin)
                if _subtract(coordinates[final], coordinates[state]) != current_displacement(seam, record[1]):
                    raise AssertionError('native record and actual load displacement disagree')
                total = [x+y for x, y in zip(total, final, strict=True)]
                count += 1
        expected = tuple(sum(mean[i][j]*state[j] for j in range(12)) for i in range(12))
        if tuple(F(n, 60) for n in total) != expected:
            raise AssertionError('uniform native primitive does not implement I-L/60')
    projector = c.slow_band_projector_exact()
    if c.q5_matmul(projector, projector) != [list(row) for row in projector]:
        raise AssertionError('existing slow projector failed exact idempotence')
    trace = c.q5(0)
    for i in range(12):
        trace = c.q5_add(trace, projector[i][i])
    if trace != c.q5(3):
        raise AssertionError('existing slow projector does not have rank three')
    radius_squared = c.q5_add(c.GOLDEN, c.q5(2))
    gram = c.intrinsic_gram_exact()
    inverse_radius_squared = c.q5_inv(radius_squared)
    for i in range(12):
        for j in range(12):
            if c.q5_mul(_dot(raw_frame()[i], raw_frame()[j]), inverse_radius_squared) != gram[i][j]:
                raise AssertionError('coordinate witness and exact source response Gram disagree')
    lap_exact = [[c.q5(int(value)) for value in row] for row in lap]
    expected_eigen = [[c.q5_mul(c.q5(5, -1), value) for value in row] for row in projector]
    if c.q5_matmul(lap_exact, projector) != expected_eigen:
        raise AssertionError('Gram band does not have the source slow eigenvalue')
    third = THIRD_SEAM
    witness = (current_displacement(native.E, 1), current_displacement(native.F, 1),
               current_displacement(third, 1))
    determinant = _determinant(witness)
    if determinant != c.q5(4):
        raise AssertionError('the three admitted calibration seams do not span the carrier')
    # E has only nonzero x; F has a nonzero z. Production uses only E and F.
    if witness[0][1:] != (c.q5(0), c.q5(0)) or witness[0][0] == c.q5(0) or witness[1][2] == c.q5(0):
        raise AssertionError('actual controlled displacement span is not rank two')
    diagnostics = [{'uniform_attempt_count': n,
                    'finite_nonconstant_rank_exact': 11,
                    'normalized_kernel_distance_to_gram_frobenius': float(np.linalg.norm(
                        c.normalized_response_kernel(n)-c.intrinsic_gram()))}
                   for n in (0, 16, 64, 256)]
    return {'protocol': 'uniform independent choice of all30 native seams and both fair placement coins',
            'primitive_branches_checked_exactly': count,
            'native_state_count': len(states), 'seams': [list(x) for x in seams],
            'mean_operator_exact': [[str(x) for x in row] for row in mean],
            'laplacian': lap.tolist(), 'registered_binding': registered_binding(),
            'slow_eigenvalue_Qsqrt5': c.q5_str(c.q5(5, -1)),
            'slow_response_Qsqrt5': c.q5_str(c.q5(F(55, 60), F(1, 60))),
            'slow_projector_Qsqrt5': [[c.q5_str(x) for x in row] for row in projector],
            'gram_Qsqrt5': [[c.q5_str(x) for x in row] for row in gram],
            'raw_coordinate_frame_Qsqrt5': [_encode_vector(x) for x in raw_frame()],
            'squared_raw_radius_Qsqrt5': c.q5_str(radius_squared),
            'calibration_three_seam_witness': [list(native.E), list(native.F), list(third)],
            'three_seam_raw_determinant_Qsqrt5': c.q5_str(determinant),
            'response_carrier_rank_exact': 3, 'original_production_displacement_span_exact': 2,
            'finite_kernel_diagnostics': diagnostics,
            'checks': {'native_mean_operator_exact': True, 'all_native_current_displacements_exact': True,
                       'projector_idempotent_trace_three_exact': True,
                       'frame_and_response_gram_identical_exact': True,
                       'slow_eigenprojector_exact': True}, 'all_pass': True}


def calibration():
    return deepcopy(_calibration())


def attach_joint(joint):
    """Decode relative carrier displacement from each actual retained history.

    The quantum state references and their complete vectors are retained.
    Working-load coordinates only audit the readout; its increments depend
    on the two already committed currents and their actual seams.
    """
    if joint.get('schema') != 'oph.common_source.native_quantum_join.v1':
        raise ValueError('the actual native/quantum joint history is required')
    if len(joint['quantum_histories']) != 3:
        raise ValueError('exactly three current-conditioned histories are required')
    histories = {row['native_current']: row for row in joint['quantum_histories']}
    _, indexed, _ = native_join._validated_histories(histories)
    rows, seen = [], set()
    c = source()
    for path in joint['paths']:
        if type(path['initial_source']) is not int or not 0 <= path['initial_source'] < len(native_join.ALL_STATES):
            raise ValueError('a native source index is required')
        initial = native_join.ALL_STATES[path['initial_source']]
        first_coin, first_current = path['first_native_record']
        a, b = path['quantum_record']
        second_coin, second_current = path['second_native_record']
        if any(type(bit) is not int or bit not in (0, 1) for bit in (first_coin, a, b, second_coin)):
            raise ValueError('committed binary records are required')
        probability_weight = native_join._number(path['history_probability'], 'history probability')
        key = (initial, first_coin, a, b, second_coin)
        if initial not in native_join.BINARY_STATES or key in seen:
            raise ValueError('complete joint paths cannot be duplicated or change preparation')
        seen.add(key)
        after_first, first = native.native(initial, native.E, first_coin)
        final, second = native_join.native_continuation(after_first, b, second_coin)
        branch_index, probability = indexed[first[1]][(a, b)]
        reference = {'current': first[1], 'branch_index': branch_index, 'record': [a, b]}
        if (first != tuple(path['first_native_record']) or second != tuple(path['second_native_record'])
                or list(final) != path['final_native_loads']
                or path['after_first_source'] != native_join.STATE_INDEX[after_first]
                or path['final_source'] != native_join.STATE_INDEX[final]
                or path['quantum_branch_reference'] != reference
                or path['second_seam'] != list(native_join.actuator_seam(b))
                or path['terminal_density_prefactor'] != '1/264'
                or abs(probability_weight-probability/264) > native_join.TOLERANCE):
            raise ValueError('a geometry path is detached from the executed joint history')
        first_shift = current_displacement(native.E, first_current)
        second_shift = current_displacement(path['second_seam'], second_current)
        displacement = _add(first_shift, second_shift)
        difference = _subtract(load_readout(final), load_readout(initial))
        if displacement != difference:
            raise AssertionError('readable current sum differs from final load audit')
        norm_squared = c.q5_mul(_dot(displacement, displacement), c.q5_inv(c.q5_add(c.GOLDEN, c.q5(2))))
        rows.append({**deepcopy(path),
                     'carrier_readout': {
                         'initial_raw_load_coordinate_Qsqrt5': _encode_vector(load_readout(initial)),
                         'final_raw_load_coordinate_Qsqrt5': _encode_vector(load_readout(final)),
                         'first_raw_record_displacement_Qsqrt5': _encode_vector(first_shift),
                         'second_raw_record_displacement_Qsqrt5': _encode_vector(second_shift),
                         'total_raw_record_displacement_Qsqrt5': _encode_vector(displacement),
                         'total_raw_centroid_displacement_Qsqrt5': _encode_vector(_scale(F(1, 2), displacement)),
                         'normalized_load_displacement_squared_Qsqrt5': c.q5_str(norm_squared)}})
    if len(rows) != 1056 or abs(sum(row['history_probability'] for row in rows)-1) > native_join.TOLERANCE:
        raise ValueError('all1056 normalized complete joint histories are required')
    return {'paths': rows, 'quantum_histories': deepcopy(joint['quantum_histories']),
            'path_count': len(rows), 'history_probability_sum': sum(row['history_probability'] for row in rows),
            'production_seams': [list(native.E), list(native.F)],
            'actual_displacement_span_rank': 2,
            'uniform_calibration_schedule_used_by_production': False,
            'all_record_load_displacement_identities_exact': True,
            'reference_to_full_quantum_terminal_states_preserved': True,
            'all_pass': True}


def three_direction_continuation(joint):
    """A second, explicit protocol: first E, then bit 0 -> F or bit 1 -> G.

    Every new second event is executed on its actual after-first work state.
    The quantum instrument and all earlier records are unchanged. The returned
    references address the complete quantum histories in ``attach_joint``.
    """
    original = attach_joint(joint)
    c = source()
    rows = []
    active = set()
    for old in original['paths']:
        row = deepcopy(old)
        bit = row['quantum_record'][1]
        seam = native.F if bit == 0 else THIRD_SEAM
        after = native_join.ALL_STATES[row['after_first_source']]
        final, record = native.native(after, seam, row['second_native_record'][0])
        initial = native_join.ALL_STATES[row['initial_source']]
        first_shift = current_displacement(native.E, row['first_native_record'][1])
        second_shift = current_displacement(seam, record[1])
        displacement = _add(first_shift, second_shift)
        if displacement != _subtract(load_readout(final), load_readout(initial)):
            raise AssertionError('three-direction record displacement failed its load audit')
        row.update(second_seam=list(seam), second_native_record=list(record),
                   final_native_loads=list(final), final_source=native_join.STATE_INDEX[final])
        readout = row['carrier_readout']
        readout.update(final_raw_load_coordinate_Qsqrt5=_encode_vector(load_readout(final)),
                       second_raw_record_displacement_Qsqrt5=_encode_vector(second_shift),
                       total_raw_record_displacement_Qsqrt5=_encode_vector(displacement),
                       total_raw_centroid_displacement_Qsqrt5=_encode_vector(_scale(F(1, 2), displacement)),
                       normalized_load_displacement_squared_Qsqrt5=c.q5_str(c.q5_mul(
                           _dot(displacement, displacement), c.q5_inv(c.q5_add(c.GOLDEN, c.q5(2))))))
        if row['history_probability'] > 0:
            if row['first_native_record'][1] != 0:
                active.add(native.E)
            if record[1] != 0:
                active.add(seam)
        rows.append(row)
    directions = (native.E, native.F, THIRD_SEAM)
    determinant = _determinant(tuple(current_displacement(seam, 1) for seam in directions))
    if active != set(directions) or determinant != c.q5(4):
        raise ValueError('the supplied quantum preparation does not activate all three directions')
    return {'protocol': 'first E; committed scalar bit 0 selects F, bit 1 selects G',
            'paths': rows, 'path_count': len(rows),
            'quantum_histories_reference': 'joint_readouts.quantum_histories',
            'history_probability_sum': original['history_probability_sum'],
            'production_seams': [list(seam) for seam in directions],
            'positive_probability_active_seams': [list(seam) for seam in directions],
            'actual_displacement_span_rank': 3,
            'three_direction_raw_determinant_Qsqrt5': c.q5_str(determinant),
            'uniform_calibration_schedule_used_by_production': False,
            'earlier_records_and_quantum_terminal_states_unchanged': True,
            'all_record_load_displacement_identities_exact': True, 'all_pass': True}


def report(joint=None):
    if joint is None:
        joint = native_join.run()
    return {'schema': 'oph.common_source.slow_geometry.v1',
            'calibration': calibration(), 'joint_readouts': attach_joint(joint),
            'three_direction_continuation': three_direction_continuation(joint),
            'source': {'native_carrier': 'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
                       'registered_binding': 'Lean/Screen/PrimitivePortFrameQuotient.lean',
                       'coordinate_field': 'Q(sqrt5); each pair encodes a+b*sqrt(5)',
                       'unit_direction': 'raw_direction/sqrt(phi+2)',
                       'readout': 'sum n_p*u_p; total-two centroid is half this vector'},
            'scope': {'physical_position_identified': False, 'physical_momentum_or_ccr_constructed': False,
                      'physical_length_or_time_scale_calibrated': False,
                      'rank_three_two_oscillator_interaction_graph_claimed': False,
                      'exact_finite_step_rank_three_claimed': False,
                      'autonomous_three_dimensional_motion_established': False,
                      'meaning': 'existing native response quotient and actual record displacements on the same joint histories'},
            'all_pass': True}
