"""Exact stationary-reference premises and fixed-sector clock classification.

The scalable construction uses total CAR number, not an unproved thermal
fixed point. It applies to every finite original Weyl/contact geometry and
every unchanged complete number-preserving repair placed on its field modes.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
import math

import numpy as np
from scipy import sparse

from ..charged_defect import branches, dictionary
from ..common_spin_protection import spin
from ..field_stress.energy import GaussianInteger as GI, commutator
from ..field_symmetry import symmetry
from ..full_fock_field.field import FieldGeometry
from ..observer import native
from ..repair_spectral_sectors import algebra as spectral


def exact_number_premises():
    """Exact local number conservation plus the full original finite action."""
    number = np.diag([j.bit_count() for j in range(16)]).astype(np.int64)
    ground = np.zeros((16, 2), dtype=complex)
    ground[5, 0], ground[10, 0] = -1j, 1j
    ground[3, 1] = ground[12, 1] = 1
    if np.any(number@ground-2*ground):
        raise ArithmeticError('the analytic ground-generating span lost number two')
    checked = []
    for kind, charge, entries in (
        ('hole', 1, ((1, 0, 1), (2, 1, 1), (4, 1, 1j), (8, 0, 1j))),
        ('particle', 3, ((7, 0, -1j), (11, 1, 1j), (13, 1, -1), (14, 0, 1)))):
        v = np.zeros((16, 2), dtype=complex)
        for row, column, value in entries:
            v[row, column] = value
        if np.any(number@v-charge*v):
            raise ArithmeticError('the actual analytic charged doublet lost its fixed number')
        checked.append({'kind': kind, 'number': charge})
    # Every monomial c_i^* c_j and n_up n_down preserves number. This
    # concrete full-Fock check attaches that factorwise proof to the source.
    field = spin.operators()
    n = GI.from_exact(sparse.diags([j.bit_count() for j in range(4096)], dtype=np.int64))
    if not commutator(field.H, n).is_zero() or not commutator(field.K, n).is_zero():
        raise ArithmeticError('the unchanged field does not conserve total number')
    laws = [branches.decomposition_certificate(law) for law in ('C', 'E')]
    return {
        'analytic_ground_generating_columns': 2,
        'analytic_charged_doublets': checked,
        'full_Fock_action_dimension_checked': 4096,
        'exact_H_and_K_number_commutators_zero': True,
        'arbitrary_finite_geometry_proof': '[N,c_i^*c_j]=0 and [N,n_i n_j]=0 term by term',
        'arbitrary_repair_placement_proof': 'fixed-number local projectors and number-preserving graded block permutations; outside projector I-P retained',
        'literal_native_ports': 12,
        'complete_native_laws': laws,
        'full_extension_unitality': 'sum native L L*=P and outside Q Q*=I-P',
        'additional_ports_or_modified_repair_law_assumed': False,
    }


def fixed_reference_classification():
    """Exact active moments supplement the existing90/15 joint algebra proof."""
    h, raising, lowering, _, _ = spectral.exact_data()
    p4 = 4*np.eye(len(h), dtype=np.int64)-raising.T@raising-lowering.T@lowering
    if np.any(p4@p4-4*p4) or np.any(h@p4-p4@h) or np.trace(p4) != 360:
        raise ArithmeticError('exact active projector or its action invariance failed')
    if np.any(raising@p4) or np.any(lowering@p4):
        raise ArithmeticError('the active projector violates its spin/eta constraints')
    trace = F(int(np.trace(p4@h)), 64)
    square_trace = F(int(np.trace(p4@h@h)), 1024)
    variance = square_trace/90-(trace/90)**2
    if (trace, square_trace, variance) != (F(-135), F(1815, 8), F(13, 48)):
        raise ArithmeticError('the exact active energy variance changed')
    cyclic = spectral.cyclic_certificate()
    dark = spectral.dark_certificate()
    return {
        'domain': 'same six-site N=4,S=1 sector, dimension315',
        'fixed_point_identity': 'sum_r ||[L_r,X]||_HS^2 = 2(||X||_HS^2-Re Tr(X* Phi(X))) for Hermitian X and a complete unital channel',
        'common_fixed_states_equal_joint_commutant': True,
        'requires_all_Hamiltonian_times_and_each_complete_repair': True,
        'single_stroboscopic_composite_fixed_points_not_classified': True,
        'active_algebra': cyclic['active_joint_star_algebra'],
        'active_commutant': cyclic['active_joint_commutant'],
        'dark_distinct_energies': dark['distinct_energies'],
        'dark_energy_block_dimensions': [9]+[3]*12,
        'general_faithful_common_reference': 'I90 tensor A direct_sum B0 direct_sum B1 ... direct_sum B12, A>0 in M3, B0>0 in M9, Bj>0 in M3',
        'normalization': '90 Tr(A)+sum_j Tr(Bj)=1',
        'commutant_complex_dimension': 9+81+12*9,
        'active_exact_projector': '(4I-Splus* Splus-eta_minus* eta_minus)/4 on the120-dimensional N4,Mspin1 space',
        'active_K_trace': trace, 'active_K_squared_trace': square_trace,
        'active_K_variance': variance,
        'active_energy_modular_distance_lower': '|beta| sqrt(13/48) > |beta|/2 for beta !=0, even after arbitrary scalar and protected-spin terms',
        'nonzero_beta_original_Gibbs_is_not_common_stationary': True,
        'positive_nontracial_sector_reference': 'choose A diagonal(1,3,9) and strictly positive dark energy weights, then divide by total trace',
        'arbitrary_sector_reference_source_preparation_claimed': False,
        'full_Fock_number_reference_is_scalar_on_this_fixed_N_sector': True,
    }


def _triangle(law='E', seam=(0, 1)):
    pairs = tuple(combinations(range(3), 2)); ports = spin.BLOCK_PORTS
    matrices = []
    for coin in (0, 1):
        for raw in range(4 if law == 'C' else 8):
            matrix = np.zeros((3, 3), dtype=np.int64)
            for column, pair in enumerate(pairs):
                state = tuple(int(p in (ports[pair[0]], ports[pair[1]])) for p in range(12))
                output, (_, current) = native.native(state, tuple(sorted(ports[j] for j in seam)), coin)
                record = {0: 0, 1: 1, -1: 2}[current]
                if law == 'E':
                    record += 4*((state[ports[seam[0]]]+state[ports[seam[1]]]) % 2)
                if record != raw:
                    continue
                target = tuple(j for j, p in enumerate(ports) if output[p])
                labels = ([seam[1] if j == seam[0] else seam[0] if j == seam[1] else j
                           for j in pair] if current else list(pair))
                matrix[pairs.index(target), column] = -1 if labels[0] > labels[1] else 1
            matrices.append((('native', coin, raw), matrix))
    matrices.append((('outside',), np.zeros((3, 3), dtype=np.int64)))
    if (np.any(sum(m.T@m for _, m in matrices)-2*np.eye(3))
            or np.any(sum(m@m.T for _, m in matrices)-2*np.eye(3))):
        raise ArithmeticError('literal native triangle ceased to be complete and unital')
    return matrices


def feedback_counterexample():
    """A complete retained-record feedback schedule need not be unital."""
    first, second = _triangle(seam=(0, 1)), _triangle(seam=(0, 2))
    gram4 = np.zeros((3, 3), dtype=np.int64)
    output4 = gram4.copy(); records = 0
    for label, a in first:
        if label == ('native', 0, 5):
            for _, b in second:
                c = b@a
                gram4 += c.T@c; output4 += c@c.T; records += 1
        else:
            gram4 += 2*a.T@a; output4 += 2*a@a.T; records += 1
    if np.any(gram4-4*np.eye(3)) or np.any(output4-np.diag([5, 4, 3])):
        raise ArithmeticError('complete adaptive counterexample no longer has its exact bias')
    return {'domain': 'highest-spin three-position repair corner',
            'first_complete_E_seam': [0, 1],
            'feedback': 'after native(coin0,raw5) apply complete E seam(0,2); otherwise stop',
            'complete_retained_record_count': records,
            'trace_preserving_Gram_numerator_over4': gram4.tolist(),
            'unitality_numerator_over4': output4.tolist(),
            'initial_corner_state': ['1/3']*3,
            'final_corner_state': ['5/12', '1/3', '1/4'],
            'reference_trace_distance': F(1, 12),
            'all_zero_and_outside_outcomes_counted': True,
            'record_independent_scheduler_premise_is_necessary': True}


def numerical_stationarity():
    """Independent full-space diagnostic, not the exact stationarity proof."""
    rows = []
    for blocks in (2, 3):
        count = 4*blocks; dimension = 1 << count
        # Exact rational ideal number reference; actual source changes only
        # the common ratio and scalar, to which the same proof applies.
        weights = np.array([3**j.bit_count()/4**count for j in range(dimension)])
        rho = sparse.diags(weights, format='csr')
        for kind in ('hole', 'particle'):
            for law in ('C', 'E'):
                output = sparse.csr_matrix(rho.shape, dtype=complex)
                for _, op in dictionary.full_extension(blocks, kind, 0, (0, 1), law):
                    output += op@rho@op.conj().T
                error = float(sparse.linalg.norm(output-rho))
                if not math.isfinite(error) or error > 1e-12:
                    raise ArithmeticError('actual complete repair stationarity diagnostic failed')
                rows.append({'blocks': blocks, 'law': law, 'kind': kind,
                             'full_Fock_dimension': dimension, 'Frobenius_residual': error})
    return rows


def aligned_spin_reference():
    """A separate exact ideal reference with nontrivial fixed-N modular action."""
    local = spin.local_certificate()
    return {'finite_aligned_chain_target': 'rho=3^(2S_z_hidden)/(16/3)^L',
            'local_mode_probabilities': '3/4 and1/4, exchanged on alternating sites',
            'faithful_local_site_eigenvalues': ['3/16', '9/16', '1/16', '3/16'],
            'local_modular_generator_up_to_scalar': '-2 log(3) S_z_hidden',
            'exact_target_fixed_by_field_and_each_aligned_complete_repair': True,
            'local_analytic_proof': local,
            'same_literal_pulse_exact_preparation_not_claimed_for_this_spin_target': True,
            'nonaligned_3D_swaps_not_included': True}


def history_theorem():
    # For two dimers the charged-hole projector has rank4 at N=2. Each
    # nonzero native hold outcome is P/sqrt2; all other native slots are zero.
    # These exact probabilities expose the distinction between field and
    # conditioned/reference-record stationarity without any eigensolver.
    mass = F(4*3**2, 4**8)
    return {
        'reference_success_block': 'sigma_M=diag(a,b) tensor M = a^M exp(gamma N), gamma=log(b/a)',
        'reference_success_trace': '(a+b)^M; never replaced by one in the complete output',
        'faithful_on_full_proper_field_algebra': True,
        'field_stationarity': '[H,N]=0 implies exp(-itH) sigma_M exp(itH)=sigma_M for all real t',
        'repair_stationarity': '[L_r,N]=0 and sum L_r L_r*=I imply sum L_r sigma_M L_r*=sigma_M',
        'record_retaining_covariance': '(sigma_M^(it) tensor I_record) V = V sigma_M^(it), V=sum_r L_r tensor |r>',
        'schedule': 'any finite record-independent sequence of original Hamiltonian durations and complete unital repairs',
        'abort_policy': 'prepare every mode; retain every local export record; freeze every register sector containing an abort flag and apply field/repair only to the all-proper sector',
        'complete_field_plus_abort_marginal_stationary': True,
        'marginal_traces_all_successful_occupation_readout_and_later_repair_records': True,
        'growing_joint_record_state_stationary': False,
        'conditioned_native_hold_probability_for_ideal_two_dimer_reference': mass/2,
        'conditioned_native_hold_state': 'P/4, supported at N=2',
        'conditioned_hold_distance_from_unconditioned_reference': 1-mass,
        'ideal_history_does_not_accumulate_preparation_error': True,
        'implemented_clock_or_field_word_errors_must_be_added_per_use': True,
        'physical_duration_or_Hamiltonian_temperature_identified': False,
    }
