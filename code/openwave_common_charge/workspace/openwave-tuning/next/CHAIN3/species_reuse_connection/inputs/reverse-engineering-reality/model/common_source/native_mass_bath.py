"""A fixed Hermitian incidence bath from a complete native selector attempt.

The stochastic source determines D=I-T=F*F. Choosing its Hermitian Dirac
dilation as a coherent bath is an additional repair/service choice, not an
identification of dissipation with energy. All port, seam, timeout, padding,
current-dark and light-response states remain in the declared operator.
No measured mass, spectral tuning, source-word execution or physical pole
identification is used or claimed here.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
import math

import numpy as np
from scipy import sparse
import sympy as s

from . import native_gauge_object as objects, native_mass_dressing as dressing
from . import native_repair_mass as feedback
from ..native_modular_reference import local_preparation as preparation
from ..observer import native

ROOT = objects.ROOT
ADMISSION = 2e-10
BATH_DIMENSION = 64
SCHEMA = {
    'bath_incidence': ((21, 12), 'float64'),
    'bath_dirichlet': ((12, 12), 'float64'),
    'bath_hamiltonian': ((64, 64), 'float64'),
    'bath_family_embedding': ((64, 3), 'float64'),
    'bath_current_frame': ((12, 15, 15), 'complex128'),
    'bath_bare_response': ((90, 90), 'float64'),
    'bath_poles': ((64,), 'float64'),
    'bath_charged_hamiltonians': ((3, 134, 134), 'float64'),
    'bath_charged_poles': ((3, 134), 'float64'),
    'bath_charged_residues': ((3, 47, 6, 6), 'float64'),
    'bath_neutral_hamiltonian': ((67, 67), 'float64'),
    'bath_neutral_poles': ((67,), 'float64'),
    'bath_neutral_residues': ((23, 3, 3), 'float64')}


def _zero(value):
    return all(s.expand(x) == 0 or s.simplify(x) == 0 for x in value)


@lru_cache(maxsize=1)
def source():
    """One unconditional selector attempt; no continuation after timeout."""
    packet = preparation.selectors(1)
    weights = tuple(sorted(packet['completion_accepted_seam_weights'].items()))
    timeout = s.Rational(packet['completion_timeout_probability'])
    if len(weights) != 20 or timeout != s.Rational(5, 32) or sum(w for _, w in weights)+timeout != 1:
        raise ArithmeticError('the frozen complete selector law changed')
    F = s.zeros(21, 12); T = timeout*s.eye(12)
    for row, ((a, b), weight) in enumerate(weights):
        F[row, a], F[row, b] = s.sqrt(s.Rational(weight)/2), -s.sqrt(s.Rational(weight)/2)
        endpoint = s.eye(12)
        endpoint[a, a] = endpoint[a, b] = endpoint[b, a] = endpoint[b, b] = s.Rational(1, 2)
        T += s.Rational(weight)*endpoint
    # Independently enumerate every actual native input, coin and timeout.
    states = native.states(binary=True); index = {x: i for i, x in enumerate(states)}
    transition = s.zeros(66)
    for column, state in enumerate(states):
        for seam, weight in weights:
            for coin in (0, 1):
                output, _ = native.native(state, seam, coin)
                transition[index[output], column] += s.Rational(weight)/2
        transition[column, column] += timeout
    loads = s.Matrix(12, 66, lambda i, j: states[j][i])
    D = s.eye(12)-T
    if (s.ones(1, 66)*transition != s.ones(1, 66) or loads*transition != T*loads
            or D != D.T or F.T*F != D or D.rank() != 9):
        raise ArithmeticError('complete source first-moment/positive Gram identity failed')
    B = objects.geometry()['frame']; metric = s.expand((B.T*B)[0, 0]); V = B/s.sqrt(metric)
    if not _zero(V.T*V-s.eye(3)):
        raise ArithmeticError('the actual family source map is not an isometry')
    core = list(preparation.CORE); outer = [i for i in range(12) if i not in core]
    zero = s.zeros(12)
    for i in core: zero[i, i] = 1
    for i in outer:
        for j in outer: zero[i, j] = s.Rational(1, len(outer))
    if zero*zero != zero or D*zero != s.zeros(12) or s.trace(zero) != 3:
        raise ArithmeticError('all source zero directions were not retained')
    overlap = (B.T*zero*B).applyfunc(s.simplify)
    null = overlap.nullspace()
    if overlap.rank() != 2 or len(null) != 1:
        raise ArithmeticError('the complete bath-to-family zero-mode rank changed')
    Y = feedback.architecture()['Y']
    q = s.simplify((null[0].T*Y*null[0])[0])
    if q.is_positive is not True:
        raise ArithmeticError('charged zero-mode exclusion lost positivity')
    G = (V.T*D*V).applyfunc(s.simplify)
    commutator = (G*Y-Y*G).applyfunc(s.simplify)
    defect = s.simplify(s.trace(commutator.T*commutator))
    if defect.is_positive is not True:
        raise ArithmeticError('the complete source bath is no longer noncommuting')
    H = s.zeros(64); H[:12, 12:33] = F.T; H[12:33, :12] = F
    W = s.zeros(64, 3); W[:12, :] = V
    return {'weights': weights, 'timeout': timeout, 'F': s.ImmutableMatrix(F),
            'D': s.ImmutableMatrix(D), 'T': s.ImmutableMatrix(T),
            'H': s.ImmutableMatrix(H), 'W': s.ImmutableMatrix(W),
            'zero_projector': s.ImmutableMatrix(zero), 'commutator_squared': defect,
            'zero_overlap_rank': 2, 'charged_zero_direction_quadratic': q}


def charged_block(left, right):
    a = source(); Y = np.asarray(feedback.architecture()['Y'], float)
    H = np.zeros((134, 134)); W = np.asarray(a['W'], float)
    H[:3, 3:6] = Y; H[3:6, :3] = Y
    H[6:70, 6:70] = np.asarray(a['H'], float); H[70:, 70:] = H[6:70, 6:70]
    H[6:70, :3] = math.sqrt(float(left))*W; H[:3, 6:70] = H[6:70, :3].T
    H[70:, 3:6] = math.sqrt(float(right))*W; H[3:6, 70:] = H[70:, 3:6].T
    return H


def neutral_block():
    a = source(); H = np.zeros((67, 67)); H[3:, 3:] = np.asarray(a['H'], float)
    H[3:, :3] = math.sqrt(float(dressing.architecture()['neutral_casimir']))*np.asarray(a['W'], float)
    H[:3, 3:] = H[3:, :3].T
    return H


def groups(positive, zero):
    """Exact zero multiplicity is grouped; nonzero poles are simple here."""
    return tuple((i,) for i in range(positive))+ (tuple(range(positive, positive+zero)),) + tuple(
        (i,) for i in range(positive+zero, 2*positive+zero))


def _spectrum(H, light, positive, zero):
    poles, vectors = np.linalg.eigh(H)
    if (np.linalg.norm(H@vectors-vectors*poles, 2) > ADMISSION
            or np.linalg.norm(vectors.T@vectors-np.eye(len(H)), 2) > ADMISSION
            or np.count_nonzero(abs(poles) <= ADMISSION) != zero
            or np.count_nonzero(poles > ADMISSION) != positive):
        raise ArithmeticError('complete bath eigenpair/zero-multiplicity admission failed')
    projected = vectors[:light]
    residues = np.asarray([projected[:, g]@projected[:, g].T for g in groups(positive, zero)])
    if np.linalg.norm(residues.sum(axis=0)-np.eye(light), 2) > ADMISSION:
        raise ArithmeticError('grouped complete bath residues lost the light identity')
    return poles, residues


def full_sparse_operator(arrays):
    """All23130 columns, including current-dark, timeout and padding states."""
    W = sparse.csr_matrix(arrays['bath_family_embedding'])
    K = []
    for rho in arrays['bath_current_frame']:
        nambu = sparse.block_diag((sparse.csr_matrix(rho), sparse.csr_matrix(rho.conj())), format='csr')
        K.append(sparse.kron(nambu, W, format='csr'))
    K = sparse.vstack(K, format='csr')
    bath = sparse.kron(sparse.eye(360, format='csr'), sparse.csr_matrix(arrays['bath_hamiltonian']), format='csr')
    return sparse.bmat([[sparse.csr_matrix(arrays['bath_bare_response']), K.conj().T], [K, bath]], format='csr')


def execute():
    a = source(); channels = dressing.architecture()['channels']
    hs, poles, residues = [], [], []
    for row in channels:
        H = charged_block(row['left'], row['right']); values, R = _spectrum(H, 6, 23, 88)
        hs.append(H); poles.append(values); residues.append(R)
    neutral = neutral_block(); npoles, nresidues = _spectrum(neutral, 3, 11, 45)
    arrays = {'bath_incidence': np.asarray(a['F'], float), 'bath_dirichlet': np.asarray(a['D'], float),
        'bath_hamiltonian': np.asarray(a['H'], float), 'bath_family_embedding': np.asarray(a['W'], float),
        'bath_current_frame': np.asarray(dressing.current_data()['matter_frame'], complex),
        'bath_bare_response': dressing.full_factors()[0].copy(),
        'bath_poles': np.linalg.eigvalsh(np.asarray(a['H'], float)),
        'bath_charged_hamiltonians': np.asarray(hs), 'bath_charged_poles': np.asarray(poles),
        'bath_charged_residues': np.asarray(residues), 'bath_neutral_hamiltonian': neutral,
        'bath_neutral_poles': npoles, 'bath_neutral_residues': nresidues}
    validate_arrays(arrays)
    return _report(arrays), arrays


def validate_arrays(arrays):
    if set(arrays) != set(SCHEMA): raise ValueError('complete bath source array menu required')
    for key, (shape, dtype) in SCHEMA.items():
        value = arrays[key]
        if (type(value) is not np.ndarray or value.shape != shape or value.dtype != np.dtype(dtype)
                or not np.isfinite(value).all()): raise ValueError('invalid complete bath array '+key)
    a = source()
    for key, expected in (('bath_incidence', a['F']), ('bath_dirichlet', a['D']),
            ('bath_hamiltonian', a['H']), ('bath_family_embedding', a['W']),
            ('bath_current_frame', dressing.current_data()['matter_frame']),
            ('bath_bare_response', dressing.full_factors()[0]),
            ('bath_poles', np.linalg.eigvalsh(np.asarray(a['H'], float)))):
        if np.max(abs(arrays[key]-np.asarray(expected, dtype=arrays[key].dtype))) > ADMISSION:
            raise ValueError('primitive full bath factor differs: '+key)
    for i, row in enumerate(dressing.architecture()['channels']):
        H = charged_block(row['left'], row['right']); poles, residues = _spectrum(H, 6, 23, 88)
        for key, expected in (('bath_charged_hamiltonians', H), ('bath_charged_poles', poles),
                              ('bath_charged_residues', residues)):
            if np.max(abs(arrays[key][i]-expected)) > ADMISSION:
                raise ValueError('complete charged bath spectrum differs: '+key)
    H = neutral_block(); poles, residues = _spectrum(H, 3, 11, 45)
    for key, expected in (('bath_neutral_hamiltonian', H), ('bath_neutral_poles', poles), ('bath_neutral_residues', residues)):
        if np.max(abs(arrays[key]-expected)) > ADMISSION:
            raise ValueError('complete neutral bath spectrum differs: '+key)


def _report(arrays):
    a = source(); bath_poles = arrays['bath_poles']
    threshold = float(bath_poles[bath_poles > ADMISSION][0]); rows = []
    for i, row in enumerate(dressing.architecture()['channels']):
        poles = arrays['bath_charged_poles'][i]; positive = poles[-23:]
        R = arrays['bath_charged_residues'][i][-23:]
        low = positive < threshold
        rows.append({'Higgs_channel': row['name'], 'left_Casimir_exact': str(row['left']),
            'right_Casimir_exact': str(row['right']), 'complete_Nambu_color_multiplicity': row['multiplicity'],
            'all_positive_poles': positive.tolist(), 'all_positive_light_residue_traces': np.trace(R, axis1=1, axis2=2).tolist(),
            'frozen_window_positive_poles': positive[low].tolist(),
            'frozen_window_light_residue_traces': np.trace(R[low], axis1=1, axis2=2).tolist(),
            'frozen_window_pole_count': int(np.count_nonzero(low)),
            'exact_zero_multiplicity': 88, 'grouped_zero_light_residue_diagnostic': float(np.linalg.norm(arrays['bath_charged_residues'][i][23], 2))})
    return {'classification': 'fixed coherent incidence bath with complete finite charged response',
        'source': {'selector_rounds': 1, 'timeout_probability_exact': str(a['timeout']),
            'accepted_seams': [[*edge, str(w)] for edge, w in a['weights']],
            'all_native_binary_load_two_inputs': 66, 'complete_reduced_one_attempt_records': 66*(2*20+1),
            'reduced_record_labels': 'input,seam,repair_coin,output or terminal timeout',
            'fine_selector_cap_address_histories_reduced_by_existing_contract': True,
            'timeout_is_terminal_in_original_native_history': True,
            'timeout_outcomes_conditioned_away': False,
            'D12': 'I-T12=F^T F; real symmetric positive semidefinite',
            'F_seam_row': 'sqrt(w_e/2)*(e_a-e_b)^T', 'F_timeout_row': 'zero; retained',
            'rank_D12': 9, 'zero_port_dimension': 3,
            'family_compressed_D_commutator_with_Y_squared_exact': str(a['commutator_squared'])},
        'fixed_rule': {'Hamiltonian': '[[H_Higgs(Y),K*],[K,I_adjoint tensor I_NambuMatter tensor H_b]]',
            'bath': 'H_b=[[0,F^T],[F,0]] plus31 zero padding rows',
            'coupling': 'K_a=rho(E_a) tensor W; W=(V_native,0); g=1',
            'family_domain_compression_used_as_closed_bath': False,
            'coherent_Hermitian_lift_is_an_added_repair_service_choice': True,
            'common_coherent_model_unit_is_physical_observer_clock': False},
        'complete_domain': {'light_Nambu_matter_family': 90, 'adjoint_current_indices': 12,
            'Nambu_matter_dimension': 30, 'bath_dimension': 64, 'unpadded_bath_dimension': 33,
            'port_rows': 12, 'seam_rows': 20, 'timeout_row': 1, 'zero_padding_rows': 31,
            'total_dimension': 23130, 'current_image_reduction_dimension': 2010,
            'current_orthogonal_dimension': 21120,
            'current_dark_copy_rule': '330 complete copies of the64-state bath',
            'current_dark_each_nonzero_bath_pole_multiplicity': 330,
            'current_dark_zero_multiplicity': 15180, 'full_zero_multiplicity': 16502,
            'full_positive_pole_multiplicity': 3314, 'full_negative_pole_multiplicity': 3314,
            'all_full_operator_columns_determined_by_retained_factors': True,
            'all_bath_poles': bath_poles.tolist(),
            'internal_domain': 'selected M15; other17 exterior states outside this declared response domain'},
        'exact_reduction': {'current_intertwiner': '(Ad(g) tensor rho(g))K_current=K_current rho(g)',
            'same_actual_current_identities_reused': 144,
            'self_energy': 'C_R V^T z*(z^2 I12-D12)^(-1) V',
            'Schur_domain': 'z outside bath spectrum; complete operator defines every bath pole too',
            'charged_zero_count_reason': 'rank(F)=9; rank(V^T PkerD V)=2; remaining scalar v^T Y v>0',
            'neutral_zero_count_reason': 'rank([V,F^T])=11, so rank(H_neutral)=22',
            'auxiliary_modes_eliminated_from_retained_spectrum': False},
        'comparison_window': {'rule_fixed_before_computation': '0<E<sqrt(lambda_min_positive(D12))',
            'upper_endpoint': threshold, 'three_charged_modes_selected_by_this_window': all(r['frozen_window_pole_count'] == 3 for r in rows),
            'extra_charged_modes_discarded': False, 'measured_mass_values_used': False},
        'charged_channels': rows,
        'neutral_channel': {'full_Nambu_multiplicity': 2, 'positive_poles': arrays['bath_neutral_poles'][-11:].tolist(),
            'zero_multiplicity_per_copy': 45, 'zero_light_residue_trace': float(np.trace(arrays['bath_neutral_residues'][11]))},
        'numerical_operator_admission': ADMISSION,
        'scope': {'stationary_averaged_Y_is_the_bare_Higgs_input': True,
            'simultaneous_finite_feedback_and_bath_executed': False,
            'coherent_bath_words_or_detector_records_executed': False,
            'auxiliary_Weyl_kinetics_physical_spin_and_particle_hole_history_attached': False,
            'physical_interacting_particle_spectrum_derived': False,
            'numerical_intervals_are_outward_physical_certificates': False,
            'rule_or_window_tuned_after_spectrum': False}}


def report_from_arrays(arrays):
    validate_arrays(arrays)
    return _report(arrays)


def report(): return execute()[0]


def source_bindings():
    paths = set(dressing.source_bindings()) | set(preparation.SOURCE_PATHS) | {
        'model/common_source/native_mass_bath.py', 'model/common_source/test_native_mass_bath.py'}
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}
