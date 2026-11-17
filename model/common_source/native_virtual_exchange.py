"""Virtual interactions from twelve corresponding native-port exchanges.

Two matching-code carriers have 24 binary native ports and total load four.
The complete load-four space has dimension10626; the work code has dimension36.
V=sum_a(I-SWAP_(a,a+12)) acts on this complete space. Its first code compression
is10/3 times identity, but B=P V Q V P is nonscalar. The intermediate vectors
are retained, including changes of each carrier's native load.

The declared protection law is H0=U[(I-P_A)+(I-P_B)], with the actual matching
projectors. It is NOT the recorded C/E repair cycle, which is logically identity.
Second-order elimination gives (10t/3)I-t^2 B/(2U), with an explicit finite
domain and spectral error. A separate finite port-rotation/code-phase word
produces the same B without assuming a continuous penalty. Its controls are
declared; neither the code phases nor a physical time/energy scale are selected
by C/E repair. The port word is a coherent 24-port calculation, not a new
compilation into the existing load-two nine-qubit native service.

The existing logical occupation is tested, not renamed electric charge. The
actual response M6 gauge generators are also tested under an explicitly
declared coordinate attachment to code6. Raw B preserves neither that gauge
action nor logical occupation. Its conditional gauge projection is calculated
exactly using the already established u3+so3 source representation. No exterior
matter lift, chiral field, Yukawa law or measured particle mass is inferred.

Nearest framework: Bravyi, DiVincenzo and Loss, arXiv1105.0675, degenerate
Schrieffer-Wolff theory. Here the model mapping, second-order matrix, elementary
Feshbach bound and finite commutator word are calculated explicitly.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import eigh, expm
from scipy.sparse.linalg import eigsh
import sympy as s

from tools import strict_json
from . import record_code as code, parameter_selection as selection, slow_geometry
from . import repair_field_attachment as field
from ..observer.source_scalar import compiler

ROOT = Path(__file__).resolve().parents[2]
ANGLE = F(1, 256)
PENALTY = F(1)
EXCHANGE_DIVISOR = 64
NUMERICAL_GATE = 2e-10


def _positive(value):
    if type(value) not in (int, F) or value <= 0:
        raise ValueError('an original positive exact rational is required')
    return F(value)


@lru_cache(maxsize=1)
def geometry():
    """All native load-four configurations and the actual matching supports."""
    masks = tuple(sum(1 << p for p in occupied)
                  for occupied in itertools.combinations(range(24), 4))
    index = {mask: i for i, mask in enumerate(masks)}
    supports = tuple(tuple(sum(bit << p for p, bit in enumerate(state))
                           for state in row) for row in code.SUPPORTS)
    if len(masks) != 10626 or len(supports) != 6 or any(len(row) != 6 for row in supports):
        raise ArithmeticError('native load or matching-code carrier changed')
    if any(mask.bit_count() != 2 for row in supports for mask in row):
        raise ArithmeticError('matching work states must have native load two')
    work = np.array([[index[a | b << 12] for a in left for b in right]
                     for left in supports for right in supports], dtype=np.int64)
    if len(set(work.ravel())) != 1296:
        raise ArithmeticError('work-code supports are not disjoint')
    permutations = []
    for port in range(12):
        row = []
        for mask in masks:
            other = mask
            if ((mask >> port) ^ (mask >> (port+12))) & 1:
                other ^= (1 << port) | (1 << (port+12))
            row.append(index[other])
        permutations.append(np.array(row, dtype=np.int64))
    native_two = tuple(sum(1 << p for p in pair)
                       for pair in itertools.combinations(range(12), 2))
    groups_a = np.array([[index[a | b << 12] for a in row]
                         for row in supports for b in native_two], dtype=np.int64)
    groups_b = np.array([[index[a | b << 12] for b in row]
                         for a in native_two for row in supports], dtype=np.int64)
    return masks, work, tuple(permutations), groups_a, groups_b


def _group_matrix(groups, dimension):
    size = groups.shape[1]
    return sparse.csr_matrix((np.ones(len(groups)*size*size, dtype=np.int64),
                             (np.repeat(groups, size, axis=1).ravel(),
                              np.tile(groups, (1, size)).ravel())),
                            shape=(dimension, dimension))


@lru_cache(maxsize=1)
def source_matrices():
    masks, work, permutations, groups_a, groups_b = geometry()
    dimension = len(masks)
    eye = sparse.eye(dimension, format='csr', dtype=np.int64)
    v = sparse.csr_matrix((dimension, dimension), dtype=np.int64)
    for permutation in permutations:
        swap = sparse.csr_matrix((np.ones(dimension, dtype=np.int64),
                                  (permutation, np.arange(dimension))),
                                 shape=(dimension, dimension))
        v = v+eye-swap
    frame_integer = sparse.csr_matrix((np.ones(work.size, dtype=np.int64),
                                      (work.ravel(), np.repeat(np.arange(36), 36))),
                                     shape=(dimension, 36))
    pa_integer = _group_matrix(groups_a, dimension)
    pb_integer = _group_matrix(groups_b, dimension)
    return v, frame_integer, pa_integer, pb_integer


@lru_cache(maxsize=1)
def exact_interaction():
    """Integer source calculation; no eigenvalue fitting or discarded shell."""
    v, u, pa, pb = source_matrices()
    first = (u.T@v@u).toarray()
    if not np.array_equal(first, 120*np.eye(36, dtype=np.int64)):
        raise ArithmeticError('first-order native compression changed')
    # U=u/6, Q V U=(3Vu-10u)/18. Native matrices are exact integers.
    a_integer = 3*(v@u)-10*u
    if (pa@a_integer).nnz or (pb@a_integer).nnz:
        raise ArithmeticError('virtual transitions are not in both code complements')
    gram = (a_integer.T@a_integer).toarray()
    b = s.ImmutableMatrix(36, 36, lambda i, j: s.Rational(int(gram[i, j]), 324))
    if b != b.T or any(value < 0 for value in b):
        raise ArithmeticError('exact virtual interaction lost its Gram structure')
    return b, a_integer


def protection_bound(coupling, penalty):
    """Ordered36 low eigenvalues, with a uniform Feshbach remainder.

    H0 has gap U; V>=0 and ||V||<=8 on total native load four. Min-max gives
    0<=E_j<=a*t, a=10/3. QVU lies in H0's2U eigenspace. The resolvent identity
    bounds the error from replacing (QHQ-E)^-1 by (QH0Q)^-1. Its sandwich norm
    is at most beta*(8t+at)/(2U*(U-at)), beta>=||B||. Ordered Feshbach eigenvalue
    branches decrease with E, so the same bound applies to all36 ordered roots.
    """
    t, U = _positive(coupling), _positive(penalty)
    b, _ = exact_interaction()
    beta = max(sum(F(x) for x in b.row(i)) for i in range(36))
    upper = F(10, 3)*t
    admitted = upper < U
    bound = t*t*beta*(8*t+upper)/(2*U*(U-upper)) if admitted else None
    return {'admitted': admitted, 'coupling': str(t), 'code_penalty': str(U),
            'low_band_count': 36, 'low_band_energy_upper': str(upper),
            'native_exchange_operator_norm_upper': '8',
            'virtual_Gram_operator_norm_upper': str(beta),
            'ordered_second_order_energy_error_upper': str(bound) if admitted else None,
            'continuous_penalty_is_the_recorded_CE_cycle': False}


def finite_penalty_spectrum():
    """All36 low eigenvalues of the complete10626-dimensional native law."""
    t = F(field.logical_parameters()['exchange_strength'])/EXCHANGE_DIVISOR
    bound = protection_bound(t, PENALTY)
    v, u, pa, pb = source_matrices()
    h0 = 2*sparse.eye(v.shape[0], format='csr')-(pa+pb)/6
    h = float(PENALTY)*h0+float(t)*v
    values, vectors = eigsh(h, k=36, which='SA', tol=1e-12, ncv=100,
                            v0=np.linspace(1, 2, h.shape[0]))
    order = np.argsort(values)
    values, vectors = values[order], vectors[:, order]
    residual = np.linalg.norm(h@vectors-vectors*values, 2)
    if residual > NUMERICAL_GATE:
        raise ArithmeticError('complete native low-band eigenpair admission failed')
    b, _ = exact_interaction()
    approx = eigh(float(F(10, 3)*t)*np.eye(36)
                  -float(t*t/(2*PENALTY))*np.asarray(b, float), eigvals_only=True)
    distance = float(np.max(np.abs(values-approx)))
    if distance > float(F(bound['ordered_second_order_energy_error_upper']))+NUMERICAL_GATE:
        raise ArithmeticError('complete native spectrum exceeds the Feshbach bound')
    return {'protection_certificate': bound, 'complete_native_dimension': h.shape[0],
            'complete_low_energies': list(map(float, values)),
            'second_order_low_energies': list(map(float, approx)),
            'spectra_are_physical_particle_masses': False}, {
                'complete_low_energies': values, 'second_order_low_energies': approx}


def _code_phase(columns, angle):
    work = geometry()[1]
    result = np.array(columns, dtype=np.complex128, copy=True)
    averages = columns[work].mean(axis=1)
    result[work] += (np.exp(1j*angle)-1)*averages[:, None, :]
    return result


def _port_rotation(columns, parameter):
    result = np.array(columns, dtype=np.complex128, copy=True)
    cosine, sine = math.cos(parameter), math.sin(parameter)
    phase = np.exp(-1j*parameter)
    for permutation in geometry()[2]:
        result = phase*(cosine*result+1j*sine*result[permutation])
    return result


def _odd_word(columns, parameter):
    # Symmetric splitting for K=(V-RVR)/2, R=I-2P.
    result = _port_rotation(columns, parameter/4)
    result = _code_phase(result, math.pi)
    result = _port_rotation(result, -parameter/2)
    result = _code_phase(result, math.pi)
    return _port_rotation(result, parameter/4)


def commutator_word(angle):
    """All36 work inputs, all10626 output amplitudes, without code projection."""
    theta = _positive(angle)
    if theta > F(1, 8):
        raise ValueError('the finite-word comparison domain requires angle<=1/8')
    _, u, _, _ = source_matrices()
    output = np.asarray(u.toarray(), dtype=np.complex128)/6
    for axis, signed in (('Y', -theta), ('K', -theta), ('Y', theta), ('K', theta)):
        if axis == 'Y':
            output = _code_phase(output, -math.pi/2)
        output = _odd_word(output, float(signed))
        if axis == 'Y':
            output = _code_phase(output, math.pi/2)
    b, _ = exact_interaction()
    target = (u/6)@expm(2j*float(theta*theta)*np.asarray(b, float))
    error = float(np.linalg.norm(output-target, 2))
    # B has nonnegative entries and row sums<=80/9<9, so ||K||<3.
    # The exact four-factor word has third derivative norm<=12^3.
    # One Strang remainder is<=(8^3+3^3)theta^3/6. The target
    # exponential's second remainder is<=2*9^2 theta^4.
    allowance = F(1942, 3)*theta**3+162*theta**4
    gram_error = float(np.linalg.norm(output.conj().T@output-np.eye(36), 2))
    if error > float(allowance)+NUMERICAL_GATE or gram_error > NUMERICAL_GATE:
        raise ArithmeticError('complete port-level commutator word failed its domain')
    return {'angle_exact': str(theta), 'effective_interaction_parameter': str(2*theta*theta),
            'effective_parameter_is_physical_time': False,
            'complete_operator_error_upper': str(allowance),
            'port_exchange_rotations': 144, 'joint_code_projector_phases': 12,
            'input_columns': 36, 'retained_native_output_rows': 10626,
            'code_projection_or_postselection_applied': False,
            'continuous_penalty_required_for_this_word': False,
            'nine_qubit_load_two_native_service_compilation_executed': False}, {
                'port_word_output': output, 'port_word_target': np.asarray(target)}


def port_macro_controls():
    """All four columns of the existing source compiler's local port macro.

    The full24-port word uses its ideal comparison law. The numerical local
    macro admission can be propagated by a144-term telescoping inequality;
    joint matching-code phase services and physical service accuracy remain
    separate premises. No leakage is projected out of these local4x4 controls.
    """
    basis = np.eye(4)
    basis[1:3, 1:3] = np.array([[1, 1], [-1, 1]])/math.sqrt(2)
    encoder = compiler.compile_real_orthogonal(basis)
    swap = np.eye(4)
    swap[[1, 2]] = swap[[2, 1]]
    rows, arrays = [], {}
    for index, angle in enumerate((-ANGLE/2, -ANGLE/4, ANGLE/4, ANGLE/2)):
        phases = [0., -2*float(angle), 0., 0.]
        tape = compiler.inverse(encoder)+compiler.compile_diagonal(phases)+encoder
        actual = compiler.execute(tape, 2)
        target = expm(-1j*float(angle)*(np.eye(4)-swap))
        if np.linalg.norm(actual-target, 2) > NUMERICAL_GATE:
            raise ArithmeticError('source-executed native-port macro exceeds admission')
        rows.append({'angle_exact': str(angle), 'source_pulses': len(tape),
                     'input_columns': 4, 'output_rows': 4})
        arrays[f'port_macro_{index}'] = actual
    return {'local_source_macros': rows, 'numerical_operator_admission': NUMERICAL_GATE,
            'whole_word_local_macro_telescoping_allowance': 144*NUMERICAL_GATE,
            'word_port_macro_source_pulses': 144*rows[0]['source_pulses'],
            'whole_word_joint_code_phases_compiled_here': False,
            'local_macro_admission_is_a_physical_service_certificate': False}, arrays


def source_generators():
    path = ROOT/'code/m1_source_realization/model.py'
    spec = importlib.util.spec_from_file_location('native_virtual_source_current', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return tuple(s.ImmutableMatrix(x) for x in module.source_generators())


def gauge_projection(matrix):
    """Exact conditional Haar projection for the established response M6.

    M6=color3 plus weak spin1. The two-carrier decomposition has color6 and
    color3bar, weak spins0/1/2, and two equivalent mixed (3,3) copies. Thus its
    commutant has dimension9. The coordinate map from code6 to this response
    carrier is an explicit extra attachment; equal dimensions do not derive it.
    """
    if not isinstance(matrix, s.MatrixBase) or matrix.shape != (36, 36):
        raise ValueError('an exact matrix on the full two-carrier space is required')
    eye = s.eye(36)
    c = s.diag(1, 1, 1, 0, 0, 0)
    w = s.eye(6)-c
    swap = s.zeros(36)
    for i in range(6):
        for j in range(6):
            swap[6*j+i, 6*i+j] = 1
    cc, ww = s.kronecker_product(c, c), s.kronecker_product(w, w)
    cw, wc = s.kronecker_product(c, w), s.kronecker_product(w, c)
    singlet = s.Matrix([int(i == j and i >= 3) for i in range(6) for j in range(6)])
    singlet = singlet*singlet.T/3
    projectors = (cc*(eye+swap)/2, cc*(eye-swap)/2, singlet,
                  ww*(eye-swap)/2, ww*(eye+swap)/2-singlet)
    names = ('color_symmetric6', 'color_antisymmetric3bar', 'weak_spin0', 'weak_spin1', 'weak_spin2')
    units = (cw, wc, cw*swap*wc, wc*swap*cw)
    projected, rows = s.zeros(36), []
    for name, p in zip(names, projectors, strict=True):
        value = s.trace(p*matrix)/s.trace(p)
        projected += value*p
        rows.append({'response_channel': name, 'dimension': int(s.trace(p)),
                     'interaction_coefficient_exact': str(value)})
    mixed = []
    for unit in units:
        value = s.trace(unit.T*matrix)/9
        projected += value*unit
        mixed.append(str(value))
    return s.ImmutableMatrix(projected), rows, mixed


def symmetry_contract():
    b, _ = exact_interaction()
    p = field.logical_parameters()
    occupations = tuple(int(label in p['occupied_labels']) for label in range(6))
    n = s.diag(*(occupations[a]+occupations[c] for a in range(6) for c in range(6)))
    parity = s.diag(*((-1)**n[i, i] for i in range(36)))
    averaged = s.ImmutableMatrix(36, 36, lambda i, j: b[i, j] if n[i, i] == n[j, j] else 0)
    bg, rows, mixed = gauge_projection(b)
    raw_counts, projected_counts = [], []
    for d in source_generators():
        current = s.kronecker_product(d, s.eye(6))+s.kronecker_product(s.eye(6), d)
        raw = (b*current-current*b).applyfunc(s.expand)
        repaired = (bg*current-current*bg).applyfunc(s.expand)
        raw_counts.append(sum(x != 0 for x in raw))
        projected_counts.append(sum(x != 0 for x in repaired))
    if any(projected_counts) or not all(raw_counts):
        raise ArithmeticError('source-current covariance diagnostic changed')
    spectrum = [{'value_exact': str(value), 'multiplicity': int(multiplicity)}
                for value, multiplicity in sorted(bg.eigenvals().items(), key=lambda x: x[0])]
    return {'logical_occupation_per_code_label': list(occupations),
            'raw_occupation_commutator_nonzeros': sum(x != 0 for x in b*n-n*b),
            'raw_fermion_parity_commutator_nonzeros': sum(x != 0 for x in b*parity-parity*b),
            'three_phase_occupation_average_is_a_separate_control': True,
            'occupation_average_is_an_SM_gauge_attachment': False,
            'actual_response_generator_count': 12,
            'code6_to_response_M6_coordinate_attachment_is_declared': True,
            'raw_source_current_commutator_nonzeros': raw_counts,
            'projected_source_current_commutator_nonzeros': projected_counts,
            'conditional_gauge_commutant_dimension': 9,
            'gauge_projection_is_an_additional_control': True,
            'gauge_projection_channels': rows,
            'mixed_copy_coefficients_exact': mixed,
            'conditional_gauge_interaction_spectrum': spectrum,
            'gauge_projected_occupation_commutator_nonzeros': sum(x != 0 for x in bg*n-n*bg),
            'M15_exterior_Higgs_or_chiral_attachment_derived': False,
            'logical_occupation_identified_with_electric_charge': False}, {
                'virtual_Gram': np.asarray(b, float),
                'occupation_averaged_Gram': np.asarray(averaged, float),
                'conditional_gauge_Gram': np.asarray(bg, float)}


@lru_cache(maxsize=1)
def persistent_port_symmetry():
    """Actual code stabilizer and failure of its symmetry-closed enlargement.

    Each orbit vector is a normalized uniform perfect matching of nonnative
    pairs. Its span is calculated exactly. For a native seam (0,b), the C
    branch with coin0/current+1 has effect n0/2 on this nonedge-supported
    span. The original code obeys V*E V=I/12. The enlarged span fails that
    same condition, with a displayed exact off-diagonal witness.
    """
    pairs = tuple(itertools.combinations(range(12), 2))
    supports = tuple(frozenset(row) for row in code.MATCHINGS)
    transformations = slow_geometry.source().incidence_automorphisms()
    orbit, stabilizers, port_stabilizers = [], [], []
    for permutation in transformations:
        images = tuple(frozenset(tuple(sorted((permutation[a], permutation[b])))
                                 for a, b in support) for support in supports)
        columns = s.Matrix([[int(pair in support) for support in images] for pair in pairs])
        orbit.append(columns)
        if all(support in supports for support in images):
            stabilizers.append(tuple(supports.index(support) for support in images))
            port_stabilizers.append(tuple(permutation))
    complete = s.Matrix.hstack(*orbit)
    pivots = complete.rref()[1]
    basis = complete[:, list(pivots)]
    gram = basis.T*basis
    projector = basis*gram.inv()*basis.T
    if projector*projector != projector or projector.T != projector:
        raise ArithmeticError('exact native symmetry-orbit projector failed')
    logical = tuple(sorted(set(stabilizers)))
    b, _ = exact_interaction()
    for permutation in logical:
        d = s.zeros(6)
        for before, after in enumerate(permutation):
            d[after, before] = 1
        action = s.kronecker_product(d, d)
        if action*b != b*action:
            raise ArithmeticError('virtual interaction lost its actual port stabilizer')
    seam = next(edge for edge in code.SEAMS if edge[0] == 0)
    # Explicit C branch: raw code1 means current+1; fair coin contributes1/2.
    effect = s.diag(*(s.Rational(int(seam[0] in pair and seam[1] not in pair), 2)
                     for pair in pairs))
    # The first automorphism need not be identity, so construct the actual code.
    original = s.Matrix([[int(pair in support) for support in supports] for pair in pairs])
    if original.T*effect*original != original.T*original/12:
        raise ArithmeticError('original matching code lost recorded-branch correctability')
    failure = basis.T*effect*basis-gram/12
    witness = next(((i, j) for i in range(basis.cols) for j in range(basis.cols)
                    if failure[i, j]), None)
    if witness is None:
        raise ArithmeticError('the symmetry enlargement unexpectedly retained correctability')
    i, j = witness
    defect = failure[i, j]/6
    first_pairs = [list(pair) for row, pair in enumerate(pairs) if basis[row, i]]
    second_pairs = [list(pair) for row, pair in enumerate(pairs) if basis[row, j]]
    report = {'native_incidence_automorphisms': len(transformations),
              'actual_code_stabilizer_order': len(stabilizers),
              'distinct_logical_stabilizer_actions': [list(row) for row in logical],
              'logical_stabilizer_group': 'Klein four',
              'virtual_interaction_preserves_actual_stabilizer': True,
              'distinct_matching_orbit_vectors': len(set(tuple(complete[:, j]) for j in range(complete.cols))),
              'exact_symmetry_closed_span_dimension': basis.cols,
              'symmetry_closed_span_is_correctable_for_complete_recorded_menu': False,
              'correctability_witness': {'law': 'C', 'seam': list(seam), 'coin': 0, 'raw_current': 1,
                  'first_matching_pairs': first_pairs, 'second_matching_pairs': second_pairs,
                  'normalized_vector_overlap': str(gram[i, j]/6),
                  'normalized_branch_off_diagonal': str((basis.T*effect*basis)[i, j]/6),
                  'required_scalar_branch_probability': '1/12',
                  'Knill_Laflamme_defect_exact': str(defect)},
              'full_response_M6_action_on_this_object_derived': False,
              'continuous_gauge_action_inferred_from_finite_stabilizer': False}
    return report, {'symmetry_orbit_basis': np.array(basis, dtype=np.int64),
                    'symmetry_orbit_projector': np.array(projector, dtype=np.float64)}


def build():
    b, a = exact_interaction()
    masks = geometry()[0]
    changed_load = np.array([(mask & 4095).bit_count() != 2 for mask in masks])
    if a[changed_load].nnz == 0 or a[~changed_load].nnz == 0:
        raise ArithmeticError('both inter-load and same-shell syndrome transitions must remain')
    spectrum, arrays = finite_penalty_spectrum()
    word, word_arrays = commutator_word(ANGLE)
    macros, macro_arrays = port_macro_controls()
    symmetries, symmetry_arrays = symmetry_contract()
    persistent, persistent_arrays = persistent_port_symmetry()
    arrays.update(word_arrays)
    arrays.update(macro_arrays)
    arrays.update(symmetry_arrays)
    arrays.update(persistent_arrays)
    report = {'schema': 'oph.native_virtual_exchange.v1',
              'source_action': 'V=sum_12_ports(I-SWAP_(a,a+12))',
              'native_total_load': 4, 'native_port_count': 24,
              'complete_native_dimension': 10626, 'work_dimension': 36,
              'first_order_code_compression_exact': '10/3*I36',
              'second_order_code_interaction_exact': '-t^2/(2U)*B',
              'virtual_Gram_exact': [[str(x) for x in b.row(i)] for i in range(36)],
              'load_changes_and_same_shell_syndromes_retained': True,
              'virtual_range_is_in_both_code_complements': True,
              'load_penalty_alone_selects_the_matching_work_code': False,
              'complete_penalty_spectrum': spectrum,
              'finite_port_word': word, 'symmetry_contract': symmetries,
              'source_executed_local_port_macros': macros,
              'persistent_object_and_actual_port_symmetry': persistent,
              'empirical_inputs_to_prediction': [],
              'numerical_admission': NUMERICAL_GATE,
              'nearest_framework': 'Bravyi-DiVincenzo-Loss, Schrieffer-Wolff degenerate perturbation theory',
              'primary_source': 'https://arxiv.org/abs/1105.0675',
              'scope': {'recorded_CE_repair_selects_this_coherent_law': False,
                        'recorded_corrected_CE_cycle_is_logically_identity': True,
                        'static_protection_or_pulse_controls_are_declared': True,
                        'complete_24_port_and_joint_code_phase_compilation_executed_here': False,
                        'observed_particle_spectrum_or_physical_mass_scale_derived': False}}
    return report, arrays


def source_bindings():
    paths = set(selection.source_bindings()) | {
        'model/common_source/native_virtual_exchange.py',
        'model/common_source/test_native_virtual_exchange.py',
        'model/common_source/repair_field_attachment.py',
        'code/m1_source_realization/model.py', 'code/source_selection_model/response.json',
        'model/observer_particles/gauge_particle_identification.py',
        'model/observer_particles/representation_identification.py'}
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def _compare(expected, actual, path='receipt'):
    if type(expected) is not type(actual):
        raise ValueError(path+' original type differs')
    if type(expected) is dict:
        if set(expected) != set(actual):
            raise ValueError(path+' keys differ')
        for key in expected:
            _compare(expected[key], actual[key], path+'.'+key)
    elif type(expected) is list:
        if len(expected) != len(actual):
            raise ValueError(path+' length differs')
        for i, (a, b) in enumerate(zip(expected, actual, strict=True)):
            _compare(a, b, path+f'[{i}]')
    elif type(expected) is float:
        exact = path.endswith(('numerical_admission', 'numerical_operator_admission',
                               'whole_word_local_macro_telescoping_allowance'))
        if not math.isfinite(actual) or (expected != actual if exact else abs(expected-actual) > 2e-11):
            raise ValueError(path+' numerical comparison failed')
    elif expected != actual:
        raise ValueError(path+' differs')


def check(path):
    path = Path(path)
    report = strict_json.load(path)
    if type(report) is not dict or report.get('source_bindings') != source_bindings():
        raise ValueError('virtual-exchange source bindings differ')
    artifact = path.with_suffix('.npz')
    if report.get('artifact') != {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}:
        raise ValueError('virtual-exchange binary custody differs')
    fresh, arrays = build()
    _compare(fresh, {key: value for key, value in report.items() if key not in ('source_bindings', 'artifact')})
    with np.load(artifact, allow_pickle=False) as stored:
        if len(stored.files) != len(set(stored.files)) or set(stored.files) != set(arrays):
            raise ValueError('virtual-exchange complete array keys differ')
        for key, value in arrays.items():
            old = stored[key]
            if old.shape != value.shape or old.dtype != value.dtype or not np.isfinite(old).all():
                raise ValueError('virtual-exchange array schema differs: '+key)
            if np.linalg.norm((old-value).ravel()) > NUMERICAL_GATE:
                raise ValueError('virtual-exchange array comparison failed: '+key)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--check', type=Path)
    mode.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.check:
        check(args.check)
        print('Native virtual exchange, source-current controls and complete port word: PASS')
        return
    artifact = args.output.with_suffix('.npz')
    if args.output.exists() or artifact.exists():
        parser.error('choose fresh receipt and artifact paths')
    report, arrays = build()
    report['source_bindings'] = source_bindings()
    with artifact.open('xb') as output:
        np.savez_compressed(output, **arrays)
    report['artifact'] = {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}
    with args.output.open('x') as output:
        json.dump(report, output, indent=2, sort_keys=True, allow_nan=False)
        output.write('\n')


if __name__ == '__main__':
    main()
