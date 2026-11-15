"""A specified twelve-current observer interface on the persistent native pair.

Native incidence labels its geometric triplets. This architecture admits
complex unitary currents on the first triplet and real rotations on the
second, including the observable relative central phase. It is a chosen
control grammar, not a claim that recorded C/E repair selects this grammar
from the separately tested seventeen-dimensional exchange response.

The same protected pair carries the current, source preparations, retained
native repair records and detector interface. The existing coherent compiler
and its load/export services remain apparatus. Public current requests are
typed separately from those apparatus operations. This does not make the
controller physically inaccessible or establish a local gauge redundancy.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm
import sympy as s

from tools import strict_json
from . import native_gauge_object as objects
from . import native_code_transport as transport
from .mass_response import history
from ..observer.source_scalar import compiler

ROOT = Path(__file__).resolve().parents[2]
SOURCE_ADMISSION = 4e-8
REPLAY_ADMISSION = 3e-9
STEP = F(1, 64)


def _clean(matrix):
    return objects._clean(matrix)


def _zero(matrix):
    expanded = s.Matrix(matrix).applyfunc(s.expand)
    return not any(expanded) or objects._zero(expanded)


def _entries(matrix):
    return objects._entries(matrix)


def _hat(vector):
    x, y, z = vector
    return s.Matrix([[0, -z, y], [z, 0, -x], [-y, x, 0]])


def _real_flat(matrix):
    return s.Matrix([s.re(x).expand() for x in matrix]+[s.im(x).expand() for x in matrix])


def _exact_blocks(blocks):
    result = s.zeros(6)
    for block, rows in enumerate(blocks):
        for i, row in enumerate(rows):
            for j, entry in enumerate(row):
                a, b, c, d = map(s.Rational, entry)
                result[3*block+i, 3*block+j] = a+b*s.sqrt(5)+s.I*(c+d*s.sqrt(5))
    return _clean(result)


@lru_cache(maxsize=1)
def exact_architecture():
    """Construct currents from native coordinates before comparing the fixture.

    The first triplet is the higher eigenvalue of the near-fivefold class
    sum, defined by actual port incidence. The axes, handed coordinate frame,
    unit response scale and complex-versus-real control choice are the stated
    architecture. Neither a measured mass nor a particle label is an input.
    """
    data = objects.geometry()
    frame = data['frame']
    r2 = s.expand((frame[0, :]*frame[0, :].T)[0])
    generators = []
    for port in range(12):
        v = frame[port, :].T
        sigma = v.xreplace({s.sqrt(5): -s.sqrt(5)})
        Q = v*v.T/2-r2*s.eye(3)/6+s.eye(3)/12
        generators.append(_clean(s.diag(_hat(v/2)+s.I*Q, _hat(sigma/2))))
    frames = objects.geometric_intertwiners()
    T = _clean(s.Matrix.hstack(*(frames[name]['frame']/s.sqrt(frames[name]['metric'])
                                 for name in ('three', 'three_prime'))))
    if not _zero(T.T*T-s.eye(6)) or not _zero(T*T.T-objects.geometric_projectors()['six']):
        raise ArithmeticError('the protected object does not match its actual geometric coordinates')
    packet = strict_json.load(ROOT/'code/source_selection_model/response.json')
    fixture = tuple(_exact_blocks(blocks) for blocks in packet['generators'])
    fixture_vectors = tuple(s.Matrix([2*s.re(g[2, 1]), 2*s.re(g[0, 2]), 2*s.re(g[1, 0])])
                            for g in fixture)
    port_map = []
    for port, generator in enumerate(generators):
        matches = [q for q, vector in enumerate(fixture_vectors) if _zero(vector-frame[port, :].T)]
        if len(matches) != 1 or not _zero(generator-fixture[matches[0]]):
            raise ArithmeticError('native geometric current differs from the finite response fixture')
        port_map.append(matches[0])
    if sorted(port_map) != list(range(12)):
        raise ArithmeticError('native-to-fixture port chart is not a bijection')
    basis = s.Matrix.hstack(*[_real_flat(g) for g in generators])
    pivot_rows = tuple(basis.T.to_DM(extension=True).rref()[1])
    if len(pivot_rows) != 12:
        raise ArithmeticError('the specified twelve-port current is not faithful')
    inverse = basis.extract(pivot_rows, range(12)).to_DM(extension=True).inv().to_Matrix()
    return {'generators': tuple(generators), 'geometric_frame': T, 'port_frame': frame,
            'native_to_fixture_port': tuple(port_map), 'fixture_intertwiner': s.eye(6),
            'central': s.diag(s.I, s.I, s.I, 0, 0, 0),
            'real_basis': basis, 'pivot_rows': pivot_rows, 'coefficient_inverse': inverse}


def admit_generator(matrix):
    """Accept precisely the public current tangent, excluding apparatus gates."""
    if not isinstance(matrix, s.MatrixBase) or matrix.shape != (6, 6):
        raise ValueError('an exact six-dimensional public current is required')
    if any(x.has(s.Float) or x.free_symbols or x.is_finite is not True for x in matrix):
        raise ValueError('finite exact numerical current entries are required')
    matrix = _clean(matrix)
    if (not _zero(matrix+matrix.conjugate().T) or not _zero(matrix[:3, 3:])
            or not _zero(matrix[3:, :3]) or not _zero(s.im(matrix[3:, 3:]))):
        raise ValueError('current lies outside the declared complex-triplet/real-triplet grammar')
    data = exact_architecture()
    flat = _real_flat(matrix)
    coefficients = _clean(data['coefficient_inverse']*flat.extract(data['pivot_rows'], [0]))
    if not _zero(data['real_basis']*coefficients-flat):
        raise ValueError('current is not a real combination of the admitted port controls')
    return tuple(coefficients)


def _controls(values):
    if type(values) is not tuple or len(values) != 12:
        raise ValueError('a tuple of twelve exact rational port controls is required')
    if any(type(value) not in (int, F) for value in values):
        raise ValueError('port controls must be original integers or Fractions')
    return tuple(F(value) for value in values)


def public_generator(controls):
    controls = _controls(controls)
    result = _clean(sum((s.Rational(v.numerator, v.denominator)*g
                         for v, g in zip(controls, exact_architecture()['generators'], strict=True)), s.zeros(6)))
    admit_generator(result)
    return result


@lru_cache(maxsize=1)
def _exact_certificate():
    data = exact_architecture()
    generators, T = data['generators'], data['geometric_frame']
    brackets = []
    for p in range(12):
        admit_generator(generators[p])
        for q in range(p):
            bracket = _clean(generators[p]*generators[q]-generators[q]*generators[p])
            coefficients = admit_generator(bracket)
            brackets.append({'ports': [p, q], 'coefficients': [str(x) for x in coefficients]})
    rotations = objects.geometric_projectors()['rotation_rows']
    near_sum = s.zeros(6); order3 = []; cayley = []
    seams = set(objects.geometry()['seams'])
    for row in rotations:
        # geometric_intertwiners already checks all60 exact identities
        # (axis(g) tensor axis(g))T=T diag(R_g,sigma(R_g)). Reuse that
        # certified coordinate action rather than expanding sqrt metrics
        # through the36-dimensional multiplication again.
        record = row['rotation3']
        R = s.Matrix(*record['shape'], lambda i, j: 0)
        for i, j, value in record['entries']:
            R[i, j] = s.sympify(value)
        U = s.diag(R, R.xreplace({s.sqrt(5): -s.sqrt(5)}))
        g = tuple(row['port_action'])
        for p in range(12):
            if not _zero(U*generators[p]*U.T-generators[g[p]]):
                raise ArithmeticError('the chosen current does not covary with actual native port rotations')
        if row['order'] == 5 and sum(tuple(sorted((p, g[p]))) in seams for p in range(12)) == 10:
            near_sum += U
        if row['order'] == 3:
            K = _clean((U-s.eye(6))*(U+s.eye(6)).inv())
            coefficients = admit_generator(K)
            if not _zero((s.eye(6)+K)*(s.eye(6)-K).inv()-U):
                raise ArithmeticError('same-current native Cayley implementation differs')
            order3.append(g); cayley.append({'port_action': list(g), 'coefficients': [str(x) for x in coefficients]})
    near_sum = _clean(near_sum)
    expected = s.diag(*([2+2*s.sqrt(5)]*3+[2-2*s.sqrt(5)]*3))
    if not _zero(near_sum-expected):
        raise ArithmeticError('native incidence no longer labels the selected complex triplet')
    paths = {tuple(range(12))}|set(order3)|{tuple(a[b[i]] for i in range(12)) for a in order3 for b in order3}
    if paths != set(objects.geometry()['rotations']):
        raise ArithmeticError('two native order-three factors do not cover the proper rotation group')
    # The kernel of conjugation on End(C6) is exactly the scalar line;
    # on End(C3)+End(C3) each block scalar is invisible independently.
    full_probes = [_real_flat(g-s.trace(g)*s.eye(6)/6) for g in generators]
    block_probes = [_real_flat(g-s.diag(s.trace(g[:3, :3])*s.eye(3)/3,
                                      s.trace(g[3:, 3:])*s.eye(3)/3)) for g in generators]
    visible = s.Matrix.hstack(*full_probes).to_DM(extension=True).rank()
    block_visible = s.Matrix.hstack(*block_probes).to_DM(extension=True).rank()
    if (visible, block_visible) != (12, 11):
        raise ArithmeticError('full-matrix current visibility or relative-phase deletion control failed')
    protection = objects.recorded_protection()
    return {'real_current_dimension': 12, 'commutator_closure': brackets,
            'complete_EndM6_tomography_rank': int(visible), 'block_diagonal_tomography_rank': int(block_visible),
            'native_port_covariance_identities': 720,
            'near_fivefold_class_sum': _entries(near_sum),
            'complex_triplet_is_higher_near_class_eigenvalue': True,
            'actual_order_three_cayley_factors': cayley,
            'all_60_proper_actions_use_at_most_two_cayley_factors': True,
            'complete_recorded_C_E_scalar_Gram_branches_per_coloring': protection['complete_raw_branches'],
            'C_E_distinct_branch_weights': protection['probability_counts_by_law'],
            'ideal_corrected_C_C_C_E_E_E_maps_all_apply_the_same_chosen_current': True,
            'grammar': 'u(3) on the selected geometric triplet; so(3) on its Galois-conjugate triplet',
            'grammar_is_selected_architecture_not_repair_selection_of_native17': True,
            'apparatus_compilation_gates_are_public_current_requests': False,
            'interface_restriction_proves_physical_inaccessibility_or_gauge_redundancy': False}


def _parameter(value):
    if type(value) not in (int, F):
        raise ValueError('an original exact rational model-coordinate parameter is required')
    return F(value)


def current_request(kind, parameter, *, port=None):
    """The fixed repair law accepts a port pulse or the uniform port field.

    These are its control inputs. Apparatus compilation tapes and arbitrary
    matrix generators are not accepted as repair requests. The old weighted
    exchange interaction is not part of this revised repair law.
    """
    parameter = _parameter(parameter)
    if kind == 'port' and type(port) is int and 0 <= port < 12:
        return {'kind': kind, 'port': port, 'parameter': str(parameter)}
    if kind == 'uniform' and port is None:
        return {'kind': kind, 'port': None, 'parameter': str(parameter)}
    raise ValueError('a fixed-law native port or uniform-field request is required')


def _request(request):
    if type(request) is not dict or set(request) != {'kind', 'port', 'parameter'}:
        raise ValueError('the fixed repair request menu differs')
    if type(request['parameter']) is not str:
        raise ValueError('the model-coordinate parameter must retain its exact rational string')
    try:
        parameter = F(request['parameter'])
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError('finite exact model-coordinate parameter required') from error
    normalized = current_request(request['kind'], parameter, port=request['port'])
    if normalized != request:
        raise ValueError('the repair request is not canonically typed')
    return normalized


def request_generator(request):
    request = _request(request)
    data = exact_architecture()
    return data['central'] if request['kind'] == 'uniform' else data['generators'][request['port']]


def _logical_word(request, targets=(3, 4, 5)):
    """Exact commuting real-rotation and symmetric-phase factorization.

    For a single port Q_p is a polynomial in v_p v_p^T, so it commutes with
    hat(v_p). The uniform field is the relative central phase. Every admitted
    primitive therefore factors exactly; no finite Trotter error is hidden.
    The numerical compiler is tested against this exact factorization.
    """
    request = _request(request)
    generator = request_generator(request)
    parameter = float(F(request['parameter']))
    K = s.re(generator); S = s.im(generator)
    if not _zero(K*S-S*K):
        raise ArithmeticError('an admitted primitive lost its exact commuting factorization')
    if request['kind'] == 'uniform':
        phases = np.zeros(8); phases[:3] = parameter
        word = compiler.compile_diagonal(phases)
    else:
        rotation = np.eye(8)
        rotation[:6, :6] = expm(parameter*np.asarray(K, float))
        # Q_p has a double perpendicular eigenvalue. Its basis must not be
        # chosen by a platform-dependent eigensolver in that degenerate plane.
        v = exact_architecture()['port_frame'][request['port'], :].T
        squared = s.expand((v.T*v)[0])
        direction = v/s.sqrt(squared)
        axis = next(i for i in range(3) if v[i] == 0)
        perpendicular = s.eye(3)[:, axis]
        second = direction.cross(perpendicular)
        exact_basis = s.Matrix.hstack(perpendicular, second, direction)
        if not _zero(exact_basis.T*exact_basis-s.eye(3)):
            raise ArithmeticError('the deterministic port eigenframe is not orthonormal')
        eigenvalues = (s.Rational(1, 12)-squared/6,)*2+(s.Rational(1, 12)+squared/3,)
        if not _zero(S[:3, :3]*exact_basis-exact_basis*s.diag(*eigenvalues)):
            raise ArithmeticError('the analytic native port phase frame differs')
        values, vectors = np.asarray(eigenvalues, float), np.asarray(exact_basis, float)
        basis = np.eye(8); basis[:3, :3] = vectors
        phases = np.zeros(8); phases[:3] = parameter*values
        word = (compiler.compile_real_orthogonal(basis.T)+compiler.compile_diagonal(phases)
                +compiler.compile_real_orthogonal(basis)+compiler.compile_real_orthogonal(rotation))
    return compiler.remap(word, targets)


def current_word(requests):
    """Source word on all64 decoded-pair rows for the declared repair control.

    The apparatus conjugates the three-bit current by the already compiled
    full64-dimensional geometric preparation extension. This retains padded
    rows and every complementary code state. Its public derivative on the
    six geometric input columns is exactly the stated current map.
    """
    if type(requests) is not tuple or not requests:
        raise ValueError('a nonempty fixed sequence of repair requests is required')
    preparation = [event for stage in transport._preparation_stages() for event in stage]
    word = compiler.inverse(preparation)
    for request in requests:
        word += _logical_word(request)
    return word+preparation


def _ideal_current(requests):
    result = np.eye(6, dtype=complex)
    for request in requests:
        result = expm(float(F(request['parameter']))*np.asarray(request_generator(request), complex))@result
    return result


def _native_factors():
    """Both complete native instruments, then their actual local decoders."""
    arrays, resources = {}, {}
    for law in ('C', 'E'):
        code = transport.PAIR_CODES[law]
        initial = code._encode(np.eye(512, 8, dtype=complex))
        rows = code.repair_branches(law, transport.SEAMS[law], initial)
        factors = np.asarray([row['amplitude'] for row in rows])
        decoded = code.decoding(np.concatenate(list(factors), axis=1))
        decoded = np.asarray([decoded[:, 8*i:8*(i+1)] for i in range(len(rows))])
        arrays['native_'+law] = factors
        arrays['decoded_'+law] = decoded
        arrays['labels_'+law] = np.asarray([(row['coin'], row['record_code']) for row in rows], np.int64)
        resources[law] = {'encoding': code.encoding_resources(), 'decoding': code.decoding_resources(),
            'repair': code.repair_resources(law, transport.SEAMS[law]),
            'recovery_max_native_pulses': max(row['recovery_resources']['native_pulses']
                for row in rows if not row['abort'])}
    return arrays, resources


def _low_channel(factor):
    # Trace all64 high-register syndromes, keeping each in the retained
    # full512-row factors. No syndrome is conditioned on or projected away.
    blocks = factor.reshape(64, 8, 8)
    return np.einsum('hia,hjc->ijac', blocks.conj(), blocks, optimize=True)


def record_effects(arrays):
    """Pull the actual pointer effects through both complete local channels."""
    left = [_low_channel(factor) for factor in arrays['decoded_C']]
    right = [_low_channel(factor) for factor in arrays['decoded_E']]
    effects = [matrix.conj().T@matrix for matrix in arrays['detector']]
    output = np.empty((len(left), len(right), len(effects), 64, 64), complex)
    for a, A in enumerate(left):
        for bit, effect in enumerate(effects):
            first = np.einsum('ijac,ikjl->ackl', A, effect.reshape(8, 8, 8, 8), optimize=True)
            for b, B in enumerate(right):
                output[a, b, bit] = np.einsum('ackl,klbd->abcd', first, B, optimize=True).reshape(64, 64)
    return output


def execute():
    """One fixed current-controlled C/E repair law with an actual pointer."""
    arrays, native_resources = _native_factors()
    primitives = tuple(current_request('port', STEP, port=p) for p in range(12))+(current_request('uniform', F(1, 2)),)
    primitive_words = tuple(_logical_word(request, targets=(0, 1, 2)) for request in primitives)
    arrays['primitive_source'] = np.asarray([compiler.execute(word, 3) for word in primitive_words])
    ideals = []
    for request in primitives:
        value = np.eye(8, dtype=complex)
        value[:6, :6] = _ideal_current((request,))
        ideals.append(value)
    arrays['primitive_ideal'] = np.asarray(ideals)
    schedules = (
        (current_request('uniform', 0),),
        (current_request('uniform', F(1, 2)),),
        (current_request('port', STEP, port=0), current_request('port', -STEP, port=1)),
        (current_request('port', -STEP, port=1), current_request('port', STEP, port=0)))
    preparation = [event for stage in transport._preparation_stages() for event in stage]
    logical_words = tuple([event for request in schedule for event in _logical_word(request)]
                          for schedule in schedules)
    words = tuple(compiler.inverse(preparation)+word+preparation for word in logical_words)
    initial = np.eye(64, dtype=complex)[:, 8*np.arange(6)]
    source = compiler.execute(preparation, 6, initial)
    decoded_source = compiler.execute(compiler.inverse(preparation), 6, source)
    moved = np.concatenate([compiler.execute(word, 6, decoded_source) for word in logical_words], axis=1)
    encoded = compiler.execute(preparation, 6, moved)
    arrays['controlled'] = np.asarray([encoded[:, 6*i:6*(i+1)] for i in range(len(schedules))])
    arrays['prepared'] = source
    arrays['decoded_preparation'] = decoded_source
    T = np.zeros((64, 6), complex); T[transport.product_addresses()] = transport.geometric_frame()
    arrays['ideal_controlled'] = np.asarray([T@_ideal_current(schedule) for schedule in schedules])
    marker = np.zeros(6); marker[0] = marker[3] = 1/math.sqrt(2)
    vector = (T@marker).real
    detector = history.projector_read(np.eye(64, dtype=complex), vector)
    if detector.failure_outputs:
        raise ArithmeticError('six-bit pointer unexpectedly added unregistered scratch records')
    arrays['detector'] = np.asarray(detector.outputs)
    pulled = record_effects(arrays)
    arrays['record_effects'] = np.asarray([
        np.einsum('ax,ijkab,by->ijkxy', A.conj(), pulled, A, optimize=True)
        for A in arrays['controlled']])
    arrays['raw_labels'] = np.asarray([(cc, cr, ec, er, bit)
        for cc, cr in arrays['labels_C'].tolist() for ec, er in arrays['labels_E'].tolist()
        for bit in (0, 1)], np.int64)
    report = {'schema': 'oph.native_sm_current.v1', 'exact': _exact_certificate(),
        'architecture': {'native_ports': 12, 'matching_choices': [0, 1],
            'clock_parameter': 'declared coherent response parameter; no laboratory time or coupling calibration',
            'repair_law': 'chosen current pulse on the geometric pair, then both complete fixed native C/E attempts and recorded recovery',
            'executed_load_schedule': 'geometric source preparation and actual inverse, current on decoded unknown inputs, geometric reencoding, both native encoders and fixed C/E recovery, both native decoders, source pointer',
            'geometric_encoder_and_inverse_executed_on_all_six_unknown_inputs': True,
            'current_schedule_includes_a_prior_native_repair_record': False,
            'terminal_policy': 'both attempts and terminal pointer are scheduled on every raw label; invalid labels remain failure outcomes and select no recovery',
            'old_seventeen_current_exchange_is_part_of_this_revised_law': False,
            'preparation_current_and_decoder_compiler': 'existing coherent twelve-port apparatus with declared load/export services',
            'passive_CE_derives_the_apparatus_current_primitives': False,
            'completeness_domain': 'the specified fixed repair law and its admitted port-control compositions',
            'all_physically_possible_controller_interventions_declared_inaccessible': False},
        'native_to_fixture_port': list(exact_architecture()['native_to_fixture_port']),
        'generators': [_entries(g) for g in exact_architecture()['generators']],
        'fixture_intertwiner': _entries(exact_architecture()['fixture_intertwiner']),
        'schedules': [list(schedule) for schedule in schedules],
        'complete_primitive_menu_execution': list(primitives),
        'resources': {'initial_geometric_preparation': compiler.resources(preparation),
            'complete_primitive_menu': [compiler.resources(word) for word in primitive_words],
            'source_current_words': [compiler.resources(word) for word in words],
            'native_factors': native_resources, 'pointer': detector.resources},
        'admissions': {'source': SOURCE_ADMISSION, 'replay': REPLAY_ADMISSION},
        'scope': {'fixed_current_repair_architecture_executed': True,
            'full_model_A1_complete_accessible_response_discharged': False,
            'complete_current_grammar_is_a_declared_operational_architecture': True,
            'physical_gauge_redundancy_or_particle_masses_derived': False,
            'spin_half_lift_or_exterior_matter_selected_by_CE': False,
            'measured_particle_parameters_used': False}}
    report['execution'] = metrics(arrays)
    return report, arrays


def metrics(arrays):
    shapes = {'native_C': (8, 512, 8), 'decoded_C': (8, 512, 8), 'labels_C': (8, 2),
              'native_E': (16, 512, 8), 'decoded_E': (16, 512, 8), 'labels_E': (16, 2),
              'prepared': (64, 6), 'decoded_preparation': (64, 6),
              'primitive_source': (13, 8, 8), 'primitive_ideal': (13, 8, 8),
              'controlled': (4, 64, 6), 'ideal_controlled': (4, 64, 6),
              'detector': (2, 64, 64), 'record_effects': (4, 8, 16, 2, 6, 6), 'raw_labels': (256, 5)}
    if type(arrays) is not dict or set(arrays) != set(shapes):
        raise ValueError('complete current-repair source array menu differs')
    for key, shape in shapes.items():
        value = arrays[key]
        dtype = np.dtype('int64') if key.startswith('labels_') or key == 'raw_labels' else np.dtype(complex)
        if type(value) is not np.ndarray or value.shape != shape or value.dtype != dtype or not np.isfinite(value).all():
            raise ValueError('original current-repair source array schema differs: '+key)
    grams, errors, totalweights = [], [], []
    for law in ('C', 'E'):
        code = transport.PAIR_CODES[law]
        record_bits = 2 if law == 'C' else 3
        labels = np.asarray([(coin, raw) for coin in (0, 1) for raw in range(1 << record_bits)], np.int64)
        if not np.array_equal(arrays['labels_'+law], labels):
            raise ValueError('a complete native coin/raw record label is missing')
        target = np.zeros((len(labels), 512, 6), complex)
        decoded_target = np.zeros_like(target)
        weights = []
        for j, (coin, raw) in enumerate(labels.tolist()):
            weight = float(code.branch_weight(law, transport.SEAMS[law], coin, raw))
            target[j] = math.sqrt(weight)*code.ideal_code_matrix()
            decoded_target[j, :6] = math.sqrt(weight)*np.eye(6)
            weights.append(weight)
        factors = arrays['native_'+law]; decoded = arrays['decoded_'+law]
        errors += [float(np.linalg.norm(factors[:, :, :6]-target)),
                   float(np.linalg.norm(decoded[:, :, :6]-decoded_target))]
        grams += [np.einsum('ria,rib->ab', factors.conj(), factors, optimize=True),
                  np.einsum('ria,rib->ab', decoded.conj(), decoded, optimize=True)]
        totalweights.append(weights)
    gram_defect = max(float(np.linalg.norm(G-np.eye(8), 2)) for G in grams)
    T = np.zeros((64, 6), complex); T[transport.product_addresses()] = transport.geometric_frame()
    preparation_error = float(np.linalg.norm(arrays['prepared']-T))
    inverse_error = float(np.linalg.norm(arrays['decoded_preparation']-np.eye(64)[:, 8*np.arange(6)]))
    current_error = max(float(np.linalg.norm(a-b)) for a, b in
                        zip(arrays['controlled'], arrays['ideal_controlled'], strict=True))
    primitive_errors = [float(np.linalg.norm(a-b, 2)) for a, b in
                        zip(arrays['primitive_source'], arrays['primitive_ideal'], strict=True)]
    marker = np.zeros(6); marker[0] = marker[3] = 1/math.sqrt(2)
    vector = T@marker
    projection = np.outer(vector, vector.conj())
    detector_target = np.asarray([projection, np.eye(64)-projection])
    detector_error = float(np.linalg.norm(arrays['detector']-detector_target))
    if max(errors+primitive_errors+[gram_defect, preparation_error, inverse_error, current_error, detector_error]) > SOURCE_ADMISSION:
        raise ValueError('the complete native current-repair history exceeds source admission')
    pulled = record_effects(arrays)
    recomputed = np.asarray([np.einsum('ax,ijkab,by->ijkxy', A.conj(), pulled, A, optimize=True)
                             for A in arrays['controlled']])
    if np.linalg.norm(recomputed-arrays['record_effects']) > 5e-12:
        raise ValueError('record effects differ from their own full native factors and actual source pointer')
    labels = np.asarray([(cc, cr, ec, er, bit) for cc, cr in arrays['labels_C'].tolist()
                         for ec, er in arrays['labels_E'].tolist() for bit in (0, 1)], np.int64)
    if not np.array_equal(arrays['raw_labels'], labels):
        raise ValueError('complete joint raw record and pointer menu differs')
    ideal = np.empty_like(recomputed)
    weights = np.outer(*totalweights)
    pointer = np.outer(marker, marker)
    for j, incoming in enumerate(arrays['ideal_controlled']):
        U = T.conj().T@incoming
        for bit, effect in enumerate((pointer, np.eye(6)-pointer)):
            ideal[j, :, :, bit] = weights[:, :, None, None]*(U.conj().T@effect@U)[None, None]
    effect_error = float(np.linalg.norm(recomputed-ideal))
    normalization = max(float(np.linalg.norm(row.sum(axis=(0, 1, 2))-np.eye(6), 2)) for row in recomputed)
    if effect_error > 12*SOURCE_ADMISSION or normalization > 12*SOURCE_ADMISSION:
        raise ValueError('complete unknown-input record response exceeds its channel admission')
    probabilities = [float(np.real(marker@row[:, :, 0].sum(axis=(0, 1))@marker)) for row in recomputed]
    expected_central = math.cos(.25)**2
    if abs(probabilities[0]-1) > 12*SOURCE_ADMISSION or abs(probabilities[1]-expected_central) > 12*SOURCE_ADMISSION:
        raise ValueError('native source records do not resolve the relative central current')
    ordered_difference = float(np.linalg.norm(arrays['controlled'][2]-arrays['controlled'][3]))
    if ordered_difference < 1e-5:
        raise ValueError('the ordered response control is not noncommuting at the fixed parameter')
    return {'unknown_geometric_input_columns': 6, 'source_current_output_rows': 64,
        'local_native_rows': 512, 'local_decoded_input_columns': 8,
        'retained_high_decoder_syndromes_per_carrier': 64,
        'joint_raw_native_records': 128, 'pointer_outcomes_per_record': 2,
        'complete_record_pointer_outcomes': 256,
        'local_native_and_decoder_dilation_errors_diagnostic': errors,
        'complete_primitive_menu_dilation_errors_diagnostic': primitive_errors,
        'source_executed_port_primitives': 12, 'source_primitive_input_and_output_dimension': 8,
        'local_complete_Gram_defect_diagnostic': gram_defect,
        'preparation_dilation_error_diagnostic': preparation_error,
        'source_preparation_inverse_error_diagnostic': inverse_error,
        'current_dilation_error_diagnostic': current_error,
        'pointer_dilation_error_diagnostic': detector_error,
        'record_effect_error_diagnostic': effect_error,
        'record_normalization_defect_diagnostic': normalization,
        'cross_triplet_marker_probabilities_diagnostic': probabilities,
        'uniform_half_parameter_marker_probability_exact': 'cos(1/4)^2',
        'ordered_response_difference_diagnostic': ordered_difference,
        'pointer_effect_uses_ideal_target_instead_of_source_word': False,
        'any_native_raw_record_or_decoder_syndrome_postselected': False,
        'full64_preparation_input_basis_source_executed': False,
        'complete_reference_linear_domain': 'six arbitrary geometric inputs, with all64 current outputs and full local8-input factors',
        'ideal_instrument_derivative': 'dK_rs/df_p at0=sqrt(w_r*w_s)(V0 tensor V1) T D_p',
        'C_C_C_E_E_E_corrected_ideal_current_response_is_same': True,
        'different_C_E_raw_records_have_been_identified': False,
        'physical_device_precision_or_finite_shot_error_certified': False}


def source_bindings():
    paths = set(objects.source_bindings())|set(transport.source_bindings())|{
        'model/common_source/native_sm_current.py', 'model/common_source/test_native_sm_current.py',
        'model/observer/source_scalar/compiler.py', 'model/common_source/mass_response/history.py',
        'code/source_selection_model/response.json', 'code/m1_source_realization/receipt.json', 'tools/strict_json.py'}
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def _compare(expected, actual, path='receipt'):
    if type(expected) is not type(actual):
        raise ValueError('current receipt original type differs: '+path)
    if type(expected) is dict:
        if expected.keys() != actual.keys():
            raise ValueError('current receipt complete keys differ: '+path)
        for key in expected:
            _compare(expected[key], actual[key], path+'.'+key)
    elif type(expected) is list:
        if len(expected) != len(actual):
            raise ValueError('current receipt complete list differs: '+path)
        for i, (a, b) in enumerate(zip(expected, actual, strict=True)):
            _compare(a, b, path+'['+str(i)+']')
    elif type(expected) is float:
        exact = path.startswith('receipt.admissions.')
        if not math.isfinite(actual) or (expected != actual if exact else not math.isclose(expected, actual, rel_tol=2e-9, abs_tol=1e-10)):
            raise ValueError('current receipt numerical value differs: '+path)
    elif expected != actual:
        raise ValueError('current receipt exact value differs: '+path)


def compare_execution(report, arrays, fresh, fresh_arrays):
    if type(arrays) is not dict or type(fresh_arrays) is not dict or arrays.keys() != fresh_arrays.keys():
        raise ValueError('retained and fresh complete current arrays differ')
    own = metrics(arrays); metrics(fresh_arrays)
    for key in arrays:
        if arrays[key].dtype.kind in 'iu':
            if not np.array_equal(arrays[key], fresh_arrays[key]):
                raise ValueError('current raw record labels differ')
        elif np.linalg.norm(arrays[key]-fresh_arrays[key]) > REPLAY_ADMISSION:
            raise ValueError('complete source current arrays differ: '+key)
    normalized = deepcopy(report)
    for key, value in own.items():
        if key.endswith('_diagnostic'):
            old = report['execution'].get(key)
            if type(value) is list:
                if type(old) is not list or len(value) != len(old) or any(type(x) is not float or
                    not math.isclose(a, x, rel_tol=2e-10, abs_tol=64*np.finfo(float).eps)
                    for a, x in zip(value, old, strict=True)):
                    raise ValueError('current diagnostics differ from their own source arrays')
            elif type(old) is not float or not math.isclose(value, old, rel_tol=2e-10, abs_tol=64*np.finfo(float).eps):
                raise ValueError('current diagnostic differs from its own source arrays')
            normalized['execution'][key] = fresh['execution'][key]
    _compare(fresh, normalized)


def check(path):
    path = Path(path); pins = source_bindings(); report = strict_json.load(path)
    if type(report) is not dict or report.get('source_bindings') != pins:
        raise ValueError('native current source bindings differ')
    artifact = path.with_suffix('.npz')
    if report.get('artifact') != {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}:
        raise ValueError('native current binary custody differs')
    with np.load(artifact, allow_pickle=False) as stored:
        if len(stored.files) != len(set(stored.files)):
            raise ValueError('duplicate native current artifact keys')
        arrays = {key: stored[key] for key in stored.files}
    fresh, output = execute()
    compare_execution({key: value for key, value in report.items() if key not in ('artifact', 'source_bindings')}, arrays, fresh, output)
    if source_bindings() != pins:
        raise ValueError('native current sources changed during replay')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--output', type=Path); mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    if args.check:
        check(args.check); print('Specified native twelve-current repair instrument: PASS'); return
    artifact = args.output.with_suffix('.npz')
    if args.output.exists() or artifact.exists():
        parser.error('fresh native current receipt and artifact paths required')
    pins = source_bindings(); report, arrays = execute()
    if source_bindings() != pins:
        raise ValueError('native current sources changed during generation')
    report['source_bindings'] = pins
    with artifact.open('xb') as output:
        np.savez_compressed(output, **arrays)
    report['artifact'] = {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}
    with args.output.open('x') as output:
        json.dump(report, output, sort_keys=True, indent=2, allow_nan=False); output.write('\n')


if __name__ == '__main__':
    main()
