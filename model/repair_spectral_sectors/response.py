"""Actual complete repair energy transfer and active/dark particle responses."""
from __future__ import annotations

from functools import lru_cache
from itertools import combinations
import math

import numpy as np
from scipy import linalg, sparse
from scipy.sparse.linalg import expm_multiply

from ..charged_defect import dictionary
from ..common_spin_protection import spin, spectrum
from ..full_fock_field.field import FieldGeometry, FullFock
from . import algebra


ADMISSION = 3e-10


@lru_cache(maxsize=1)
def active_chart():
    """Numerical chart on the EXACT ker(S+) intersect ker(eta-) domain."""
    h, raising, lowering, _, _ = algebra.exact_data()
    columns = linalg.null_space(np.vstack((raising, lowering)), rcond=1e-13)
    if columns.shape != (120, 90):
        raise ArithmeticError('numerical active chart differs from the exact multiplicity')
    matrix = columns.T@h@columns/16
    energies, vectors = np.linalg.eigh(matrix)
    full = np.zeros((4096, 90), complex)
    full[algebra.WHERE] = columns@vectors
    return spin.operators().frame.adjoint().matrix()@full, energies


def complete_work_identity():
    """Unital completeness binds W and its exact finite transfer formulas."""
    generator = spin.operators().K.matrix(16)
    outcomes = spin.repair_outcomes('hole', 'E', (0, 1))
    work = -generator.copy()
    symmetrized = sparse.csr_matrix(generator.shape, dtype=complex)
    quadratic = 0.
    for _, operator in outcomes:
        work += operator.conj().T@generator@operator
        comm = generator@operator-operator@generator
        symmetrized += (operator.conj().T@comm+comm.conj().T@operator)/2
        quadratic += float(sparse.linalg.norm(comm)**2)
    identity_error = float(sparse.linalg.norm(work-symmetrized))
    hilbert_pair = float((generator.multiply(work.T)).sum().real)
    dirichlet_error = abs(hilbert_pair+quadratic/2)
    if max(identity_error, dirichlet_error) > ADMISSION:
        raise ArithmeticError('complete Kraus energy transfer identity failed')
    # Exact two-dimer instance: every native branch holds on P.  Thus the
    # actual full instrument is precisely pinching {P,Q}, not an arbitrary
    # thermodynamic ansatz.  P=W W*/4 and16K are Gaussian dyadic integers.
    local = np.zeros((16, 2), complex)
    local[1, 0] = local[2, 1] = 1
    local[4, 1] = local[8, 0] = 1j
    integer_v = np.column_stack([np.kron(local[:, j], local[:, i]) for i in range(2) for j in range(2)])
    p4 = integer_v@integer_v.conj().T
    field = FullFock(8)
    number = sparse.diags([x.bit_count() for x in range(256)], dtype=float)
    k16 = 16*(field.hamiltonian(FieldGeometry((4, 1, 1)))-number/2).toarray()
    if (not np.array_equal(k16.real, np.rint(k16.real))
            or not np.array_equal(k16.imag, np.rint(k16.imag))):
        raise ArithmeticError('exact two-dimer transfer coefficients changed')
    # 256W = 2 P4 K16 P4 - 4(P4 K16 + K16 P4).
    transfer256 = 2*p4@k16@p4-4*(p4@k16+k16@p4)
    squared = int(round(np.vdot(transfer256, transfer256).real))
    if squared != 1536 or np.trace(transfer256) != 0:
        raise ArithmeticError('exact two-dimer work norm or zero trace changed')
    q4 = 4*np.eye(256)-p4
    leakage128 = q4@k16@integer_v
    singlet = np.array([0, 1, -1, 0])
    if not np.array_equal(leakage128.conj().T@leakage128,
                          32*(np.eye(4)+np.outer(singlet, singlet))):
        raise ArithmeticError('exact signed work spectrum lost its leakage attachment')
    actual = sum((op.conj().T@(k16/16)@op for _, op in dictionary.full_extension(2)),
                 start=np.zeros((256, 256), complex))-k16/16
    if np.linalg.norm(actual-transfer256/256) > 1e-12:
        raise ArithmeticError('actual complete two-dimer repair is not the proved pinching')
    return {'work_definition': 'W=sum_r L_r^dagger K L_r-K',
            'exact_operator_identity': 'W=1/2 sum_r (L_r^dagger[K,L_r]+[L_r^dagger,K]L_r)',
            'unital_trace_identity': 'Tr(K W)=-1/2 sum_r ||[K,L_r]||_HS^2',
            'six_site_complete_outcomes': len(outcomes),
            'six_site_work_Frobenius_norm': float(sparse.linalg.norm(work)),
            'six_site_symmetrized_identity_residual': identity_error,
            'six_site_dirichlet_identity_residual': dirichlet_error,
            'two_dimer_work_Frobenius_squared': '3/128',
            'two_dimer_nonzero_work_eigenvalues': ['+sqrt(3/512)', '-sqrt(3/512)',
                                                '+sqrt(1/512)', '-sqrt(1/512)'],
            'two_dimer_corresponding_multiplicities': [1, 1, 3, 3],
            'two_dimer_zero_work_eigenvalue_multiplicity': 248,
            'two_dimer_exact_max_absolute_energy_transfer': 'sqrt(6)/32',
            'input_can_gain_or_lose_energy': True,
            'transfer_identified_as_heat_or_controller_work': False,
            'controller_energy_model_supplied': False}


def active_response():
    columns, energies = active_chart()
    complete = np.zeros((90, 90))
    outside = np.zeros_like(complete)
    operator_residual = 0.
    for label, operator in spin.repair_outcomes('hole', 'E', (0, 1)):
        matrix = columns.conj().T@(operator@columns)
        operator_residual = max(operator_residual, float(np.linalg.norm(operator@columns-columns@matrix)))
        weights = abs(matrix)**2
        complete += weights
        if label == ('outside',):
            outside += weights
    if (operator_residual > ADMISSION or np.max(abs(complete.sum(axis=0)-1)) > ADMISSION
            or np.max(abs(complete.sum(axis=1)-1)) > ADMISSION):
        raise ArithmeticError('complete active energy-response law lost a state or outcome')
    energy_transfer = energies@complete-energies
    variance = (energies[:, None]-energies[None, :])**2
    spreads = np.sum(complete*variance, axis=0)
    if energy_transfer[0] <= 1e-8 or np.min(spreads) <= 1e-12:
        raise ArithmeticError('active spectral-disturbance control disappeared')
    return {'domain': 'N=4,S=1,eta=1 highest-weight multiplicity90',
            'all90_input_and90_output_energies': energies.tolist(),
            'complete_energy_transition_matrix': complete.tolist(),
            'outside_energy_transition_matrix': outside.tolist(),
            'mean_K_transfer_per_energy_input': energy_transfer.tolist(),
            'mean_squared_energy_change_per_input': spreads.tolist(),
            'lowest_energy_input_mean_transfer': float(energy_transfer[0]),
            'all90_eigenvectors_have_resolved_energy_spread': True,
            'operator_sector_residual': operator_residual,
            'numerical_admission': ADMISSION,
            'one_repair_instrument': 'complete E on local blocks(0,1)',
            'spectral_grouping_or_threshold_discard_used': False}


def dark_response():
    numerator, pairs = algebra.dark_columns()
    # Two exact spinless-chain Slater eigenfunctions; their eta-raised
    # spin-one highest columns retain nontrivial coherent relative phase.
    eigen_columns = []
    selected = ((1, 2), (1, 3))
    energies = []
    for p, q in selected:
        u = [math.sqrt(2/7)*math.sin(p*(x+1)*math.pi/7) for x in range(6)]
        v = [math.sqrt(2/7)*math.sin(q*(x+1)*math.pi/7) for x in range(6)]
        wedge = np.array([u[x]*v[y]-u[y]*v[x] for x, y in pairs])
        eigen_columns.append(numerator@wedge/2)
        energies.append((math.cos(p*math.pi/7)+math.cos(q*math.pi/7))/8-1)
    frame_back = spin.operators().frame.adjoint().matrix()
    modes = frame_back@np.column_stack(eigen_columns)
    initial = modes/math.sqrt(2)  # external reference labels preserve full coherence
    generator = spin.operators().K.matrix(16)
    actual = expm_multiply(-1.25j*generator, initial)
    total_native_probability = 0.
    retained = []
    # Complete stopping instrument: terminate on the first native record,
    # otherwise continue the unnormalized outside branch.  Its40 native
    # records and final outside record are all retained, with no conditioning.
    for attempt, (law, seam) in enumerate((('E', (0, 2)), ('C', (0, 1)), ('E', (1, 2)))):
        for label, operator in spin.repair_outcomes('hole', law, seam):
            output = operator@actual
            if label == ('outside',):
                continuation = output
            else:
                probability = float(np.vdot(output, output).real)
                total_native_probability += probability
                retained.append({'attempt': attempt, 'law': law, 'seam': list(seam),
                                 'record': list(label), 'probability': probability})
        actual = continuation
    retained.append({'attempt': 3, 'law': 'stop', 'seam': [], 'record': ['outside'],
                     'probability': float(np.vdot(actual, actual).real)})
    target = initial*np.exp(-1.25j*np.array(energies))
    residual = float(np.linalg.norm(actual-target))
    complete_error = math.sqrt(residual**2+total_native_probability)
    if (complete_error > ADMISSION or total_native_probability > 1e-24
            or abs(math.fsum(row['probability'] for row in retained)-1) > ADMISSION):
        raise ArithmeticError('the exact dark spectral history is not repair transparent')
    # All simple two-annihilator probes commute with eta- and cannot take an
    # eta-singlet ground to eta2 at N4.  Check the complete66-probe menu.
    _, ground = spectrum.ground()
    ops = spin.operators()
    eta_minus = (ops.frame@ops.eta_minus@ops.frame.adjoint()).matrix()
    reference_eta_residual = float(np.linalg.norm(eta_minus@ground))
    field = FullFock(12)
    dark = numerator/2
    largest = 0.
    for i, j in combinations(range(12), 2):
        response = field.annihilator(j)@(field.annihilator(i)@ground)
        largest = max(largest, float(np.linalg.norm(dark.conj().T@response)))
    if max(largest, reference_eta_residual) > ADMISSION:
        raise ArithmeticError('simple ground pair probe reached the eta2 dark sector')
    return {'selected_slater_mode_pairs': [list(x) for x in selected],
            'selected_exact_energies': ['(cos(pi/7)+cos(2*pi/7))/8-1',
                                      '(cos(pi/7)+cos(3*pi/7))/8-1'],
            'numerical_energies': energies,
            'nonzero_beat_frequency': abs(energies[0]-energies[1]),
            'supplied_Hamiltonian_duration': 1.25,
            'complete_repair_attempts': 3,
            'complete_stopping_history_outputs': retained,
            'native_record_stops_without_postselection': True,
            'native_outcome_probability': total_native_probability,
            'final_outside_outcome_probability': retained[-1]['probability'],
            'full_state_reference_error': complete_error,
            'simple_two_annihilator_probes_checked': 66,
            'largest_ground_probe_dark_overlap': largest,
            'supplied_ground_eta_singlet_residual': reference_eta_residual,
            'exact_selection_rule': '[eta_minus,c_i c_j]=0 and eta_minus|ground>=0; N4 lowest eta weight implies eta=1',
            'dark_state_preparation_and_detector_source_executed': False,
            'localized_reference_like_particle_identified': False,
            'absolute_physical_masses_predicted': 0}
