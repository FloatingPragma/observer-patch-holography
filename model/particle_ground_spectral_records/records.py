"""Same-source ground intervals, compensated probes and complete energy records.

Exact certificates use rational arithmetic. Floating spectral lines are
diagnostics of the 208-state source and carry no particle identification.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from pathlib import Path

import sympy as s

from model.particle_candidate_atlas.interacting_chiral import creator
from model.particle_field_state_response import state_response as source


ROOT = Path(__file__).resolve().parents[2]
SOURCE_INPUTS = (
    'model/particle_field_state_response/state_response.py',
    'model/particle_observer_outputs/field_route.py',
    'model/particle_candidate_atlas/interacting_chiral.py',
    'model/particle_spatial_ground_state/ground.py',
    'Lean/QFT/GroundRecordTail.lean',
)


def _rational(value, name):
    if isinstance(value, bool) or isinstance(value, float):
        raise ValueError(f'{name} must be exact rational data')
    try:
        return F(value)
    except (TypeError, ValueError, ZeroDivisionError) as error:
        raise ValueError(f'{name} must be exact rational data') from error


def _bits(bits):
    if type(bits) is not int or not 4 <= bits <= 128:
        raise ValueError('interval bits must be an integer between4 and128')
    return bits


@lru_cache(maxsize=1)
def rational_blocks():
    """All208 states, partitioned only by exact reducing components."""
    _, h = source.real_computation_chart()
    result = []
    for indices in source.connected_blocks():
        occupations = {(i // 13).bit_count() for i in indices}
        if len(occupations) != 1:
            raise ArithmeticError('a source reducing block changes occupation')
        matrix = tuple(tuple(F(h[i, j]) for j in indices) for i in indices)
        result.append((indices, occupations.pop(), matrix))
    if sorted(i for indices, _, _ in result for i in indices) != list(range(208)):
        raise ArithmeticError('the source partition omits a basis state')
    return tuple(result)


@lru_cache(maxsize=8192)
def _shift_is_positive(matrix, shift):
    """Exact Sylvester/LDL test for H-shift I>0; no numerical proposal."""
    n = len(matrix)
    diagonal = []
    lower = [[F(0) for _ in range(n)] for _ in range(n)]
    for j in range(n):
        pivot = matrix[j][j]-shift-sum(lower[j][k]**2*diagonal[k] for k in range(j))
        if pivot <= 0:
            return False
        diagonal.append(pivot)
        lower[j][j] = F(1)
        for i in range(j+1, n):
            lower[i][j] = (matrix[i][j]-sum(
                lower[i][k]*lower[j][k]*diagonal[k] for k in range(j)))/pivot
    return True


@lru_cache(maxsize=16)
def sector_minima(bits=32):
    """Certified absolute minima of every exact occupation sector."""
    _bits(bits)
    result = {}
    for occupation in range(5):
        matrices = tuple(set(matrix for _, n, matrix in rational_blocks() if n == occupation))
        lower, upper = F(-4), F(22)
        if not all(_shift_is_positive(matrix, lower) for matrix in matrices):
            raise ArithmeticError('the strict source lower bracket fails')
        if all(_shift_is_positive(matrix, upper) for matrix in matrices):
            raise ArithmeticError('the source upper bracket fails')
        while upper-lower > F(1, 1 << bits):
            middle = (lower+upper)/2
            if all(_shift_is_positive(matrix, middle) for matrix in matrices):
                lower = middle
            else:
                upper = middle
        result[occupation] = (lower, upper)
    return result


@lru_cache(maxsize=16)
def global_ground_interval(bits=32):
    """The global minimum is in occupation2, without a sector projection."""
    sectors = sector_minima(bits)
    lower, upper = sectors[2]
    if any(lo <= upper for n, (lo, _) in sectors.items() if n != 2):
        raise ArithmeticError('the ground charge sector is not separated')
    return lower, upper


@lru_cache(maxsize=16)
def charged_thresholds(bits=32):
    """Adjacent-charge sector minima relative to the full source ground."""
    e0 = global_ground_interval(bits)
    sectors = sector_minima(bits)
    return {n: (sectors[n][0]-e0[1], sectors[n][1]-e0[0]) for n in (1, 3)}


@lru_cache(maxsize=1)
def instrument():
    """Four-mode addition/removal, retaining every mode and both signs."""
    rows = []
    for sign in (1, -1):
        for mode in range(4):
            create = creator(4, mode)
            op = create if sign == 1 else create.H
            rows.append({'sign': sign, 'mode': mode,
                         'source_charge_change': 6*sign,
                         'reference_charge': -6*sign,
                         'reference_parity': 1,
                         'K': s.ImmutableSparseMatrix(s.kronecker_product(op/2, s.eye(13)))})
    effect = sum((row['K'].H*row['K'] for row in rows), s.zeros(208))
    gauss = source.exact_control()['Gauss']
    parity = s.diag(*[(-1)**((i//13).bit_count()) for i in range(208)])
    if effect != s.eye(208):
        raise ArithmeticError('the probe instrument is incomplete')
    for row in rows:
        k = row['K']
        if gauss*k-k*gauss != row['source_charge_change']*k or parity*k != -k*parity:
            raise ArithmeticError('charge or parity compensation does not close')
    return tuple(rows)


@lru_cache(maxsize=1)
def charge_conjugation():
    """Exact antiunitary C=matrix times complex conjugation on native208."""
    gamma = source.field_route.rotor_coefficients()['gamma']
    creators = tuple(creator(4, j) for j in range(4))
    particle_hole = s.eye(16)
    for op in creators:
        particle_hole = particle_hole*(op+op.H)
    exterior = s.zeros(16)
    for mask in range(16):
        vector = s.eye(16)[:, 0]
        for j in reversed(range(4)):
            if mask >> j & 1:
                transformed = sum((gamma[i, j]*creators[i] for i in range(4)), s.zeros(16))
                vector = transformed*vector
        exterior[:, mask] = vector
    flux_reversal = s.zeros(13)
    for j in range(13):
        flux_reversal[12-j, j] = 1
    matrix = s.ImmutableSparseMatrix(s.kronecker_product(exterior*particle_hole, flux_reversal))
    control = source.exact_control()
    identities = {
        'antiunitary_square_identity': matrix*matrix.conjugate() == s.eye(208),
        'Hamiltonian_invariant': matrix*control['H'].conjugate()*matrix.H == control['H'],
        'Gauss_charge_reversed': matrix*control['Gauss'].conjugate()*matrix.H == -control['Gauss'],
        'electric_flux_reversed': matrix*control['E'].conjugate()*matrix.H == -control['E'],
    }
    if not all(identities.values()):
        raise ArithmeticError('native charge conjugation failed')
    return {'matrix': matrix, 'identities': identities}


def negative_transfer_bound(mean_excess, threshold, readout_error=0):
    """Complete original instrument weight, with optional channel error."""
    epsilon = _rational(mean_excess, 'mean ground excess')
    eta = _rational(threshold, 'negative-transfer threshold')
    error = _rational(readout_error, 'complete readout error')
    if epsilon < 0 or eta <= 0 or not 0 <= error <= 1:
        raise ValueError('nonnegative excess, positive threshold and error in[0,1] are required')
    return min(F(1), epsilon/eta+error)


def dyadic_interval(value, bits=40):
    """Outward interval for an exact rational, independent of floating output."""
    value = F(value)
    scale = 1 << _bits(bits)
    numerator = (value.numerator*scale)//value.denominator
    lower = F(numerator, scale)
    return lower, lower if lower == value else F(numerator+1, scale)


def _real_probe_blocks(k):
    phases, _ = source.real_computation_chart()
    blocks = source.connected_blocks()
    owner = {index: b for b, indices in enumerate(blocks) for index in indices}
    pairs = {(owner[i], owner[j]) for i, j in k.todok()}
    result = []
    for out_block, in_block in sorted(pairs):
        output, input_ = blocks[out_block], blocks[in_block]
        part = s.ImmutableMatrix(len(output), len(input_), lambda i, j:
            s.conjugate(phases[output[i]])*k[output[i], input_[j]]*phases[input_[j]])
        result.append((output, input_, part))
    return result


def _norm_squared(matrix):
    return s.expand(sum(s.conjugate(value)*value for value in matrix))


@lru_cache(maxsize=2048)
def _real_trace_pair(q2, kernel):
    """Reuse repeated exact source blocks and discard zero summands only."""
    return sum(q2[i, j]*value for (j, i), value in kernel.todok().items())


def complete_filter_records(prepared=None, duration='1/64', steps=4):
    """Two-bin Fourier apparatus on source and its retained clock reference.

    The reference undergoes the conjugate source unitary. This declared
    counter-evolution reads energy transfer on vec(KQ), using the same H.
    External untouched references retain the complete-channel error bound.
    """
    if prepared is None:
        from model.particle_spatial_ground_state import ground
        prepared = ground.prepared_control()
    amplitudes = dict(prepared['amplitude_blocks'])
    expected = set(source.connected_blocks())
    if set(amplitudes) != expected:
        raise ValueError('all exact source blocks must be retained')
    if tuple(prepared['computational_basis_phases']) != source.real_computation_chart()[0]:
        raise ValueError('the amplitude chart differs from the native source')
    h = source.real_computation_chart()[1]
    for indices, q in amplitudes.items():
        hb = h.extract(indices, indices)
        if q.shape != hb.shape or q != q.H or q*hb != hb*q:
            raise ValueError('the source amplitude must be Hermitian and commute with the same H')
        if any(value.is_Rational is not True for value in q):
            raise ValueError('the source polynomial requires rational amplitudes in its real chart')
    squares = {indices: s.ImmutableMatrix(q*q) for indices, q in amplitudes.items()}
    norm_squared = sum(s.trace(q2) for q2 in squares.values())
    if norm_squared <= 0 or norm_squared != prepared['norm_squared']:
        raise ValueError('the complete preparation normalization differs')
    unitary = {}
    for indices in expected:
        real, imag = source.cayley_parts(indices, duration, steps)
        unitary[indices] = real, imag
    result = []
    for row in instrument():
        norm = s.S.Zero
        overlap = s.S.Zero
        for output, input_, k in _real_probe_blocks(row['K']):
            # Q commutes with the source clock. Moving both Q factors
            # together avoids multiplying its large rational entries for
            # every detector branch; this leaves the original trace intact.
            q2 = squares[input_]
            # In this native chart each CAR component has one phase, either
            # real or imaginary. Its phase cancels from K* U K U*. Rational
            # real/imaginary products avoid symbolic complex expansion.
            if all(value.is_Rational for value in k):
                kr = k
            elif all((value/s.I).is_Rational for value in k):
                kr = k/s.I
            else:
                raise ArithmeticError('a native CAR component has mixed computation-chart phases')
            norm += _real_trace_pair(q2, s.ImmutableMatrix(kr.T*kr))
            ro, io = unitary[output]
            ri, ii = unitary[input_]
            real_transfer = s.ImmutableMatrix(kr.T*ro*kr*ri.T+kr.T*io*kr*ii.T)
            overlap += _real_trace_pair(q2, real_transfer)
        for phase_bit, numerator in enumerate((norm+overlap, norm-overlap)):
            weight = F(numerator/(2*norm_squared))
            if not 0 <= weight <= 1:
                raise ArithmeticError('the complete Fourier record is not a probability')
            result.append({'sign': row['sign'], 'mode': row['mode'],
                           'phase_bit': phase_bit, 'original_probability': weight})
    if sum(row['original_probability'] for row in result) != 1:
        raise ArithmeticError('the complete original Fourier law is not normalized')
    delta = source.cayley_bound(duration, steps)
    epsilon = F(prepared['preparation_trace_error_upper'])
    return {'records': result, 'source_dimension': 208,
            'duration': F(duration), 'Cayley_steps': steps, 'Fourier_bins': 2,
            'probability_sum': F(1), 'complete_original_branch_count': len(result),
            'Gibbs_exact_clock_output_half_trace_error_upper': min(F(1), epsilon+2*delta),
            'prepared_exact_clock_output_half_trace_error_upper': min(F(1), 2*delta),
            'same_H_clock': True, 'reference_counter_evolution': True,
            'all_quantum_references_retained': True,
            'mode_and_sign_are_apparatus_records': True,
            'negative_bin_is_negative_energy_measurement': False}


def resolved_negative_record_bound(mean_excess, threshold='1/4', resolution='1/8',
                                   count=16384, duration='1/16', steps=4096):
    """Constructive resolution plan, distinct from the executed two-bin law."""
    from model.spectral_readout.kernels import real_energy_leakage
    epsilon = _rational(mean_excess, 'mean ground excess')
    eta = _rational(threshold, 'negative-transfer threshold')
    delta = _rational(resolution, 'energy resolution')
    tau = _rational(duration, 'clock step')
    if epsilon < 0 or not 0 < delta < eta:
        raise ValueError('positive resolution below threshold and nonnegative excess are required')
    leakage = real_energy_leakage(count, tau, delta, F(26**2), support_upper=26)
    clock_error = min(F(1), 2*(count-1)*source.cayley_bound(tau, steps))
    return {'Fourier_bins': count, 'clock_step': tau, 'energy_resolution': delta,
            'negative_transfer_threshold': eta,
            'entire_source_transfer_support_bound': F(26),
            'alias_probability_upper': leakage['alias'],
            'filter_leakage_upper': leakage['total'],
            'complete_Cayley_clock_error_upper': clock_error,
            'original_negative_record_probability_upper': min(F(1),
                negative_transfer_bound(epsilon, eta-delta)+leakage['total']+clock_error),
            'maximum_source_time': (count-1)*tau,
            'Cayley_substeps_per_clock_step': steps,
            'executed_full_distribution': False}


@lru_cache(maxsize=1)
def cooled_control_certificate():
    """Exact cooled-source excess and a complete finite clock apparatus."""
    from model.particle_spatial_ground_state import ground
    prepared = ground.prepared_control()
    mean = F(prepared['mean_energy'])
    e0 = global_ground_interval()
    excess = mean-e0[0]
    if excess < 0:
        raise ArithmeticError('a state lies below the certified source minimum')
    law = complete_filter_records(prepared)
    for row in law['records']:
        row['original_probability_interval'] = dyadic_interval(row.pop('original_probability'))
    return {
        'beta': prepared['bound']['beta'],
        'Taylor_order': prepared['bound']['Taylor_order'],
        'dyadic_scale': prepared['bound']['dyadic_scale'],
        'mean_source_energy_interval': dyadic_interval(mean),
        'mean_ground_excess_interval': (dyadic_interval(mean-e0[1])[0], dyadic_interval(excess)[1]),
        'negative_transfer_threshold': F(1, 4),
        'complete_negative_transfer_probability_upper': dyadic_interval(
            negative_transfer_bound(excess, F(1, 4)))[1],
        'resolved_readout_plan_from_executed_control': resolved_negative_record_bound(
            dyadic_interval(excess)[1], resolution=F(1, 16)),
        'complete_clock_apparatus': law,
    }


def spectral_diagnostics(betas=('1/64', '1', '4', '16')):
    """Floating Lehmann weights, with complete sum and exact bound comparison.

    Each input state is the Gibbs state of the identical full source. The
    reported decimal lines and residues are numerical diagnostics. Exact
    finite-state preparation errors are supplied independently by ground.py.
    """
    import numpy as np
    from scipy.linalg import eigh
    from model.particle_spatial_ground_state import ground
    phases, h = source.real_computation_chart()
    energies, occupations, columns = [], [], []
    for indices, occupation, _ in rational_blocks():
        values, vectors = eigh(np.asarray(h.extract(indices, indices).tolist(), dtype=float))
        for j, energy in enumerate(values):
            column = np.zeros(208, dtype=complex)
            column[list(indices)] = np.asarray([complex(phases[i]) for i in indices])*vectors[:, j]
            energies.append(energy); occupations.append(occupation); columns.append(column)
    order = np.argsort(energies)
    energy = np.asarray(energies)[order]
    occupation = np.asarray(occupations)[order]
    frame = np.stack(columns, axis=1)[:, order]
    transitions = []
    for row in instrument():
        k = np.asarray(row['K'].tolist(), dtype=complex)
        weights = np.abs(frame.conj().T@k@frame)**2
        transitions.append((row['sign'], row['mode'], weights))
    total = sum(weights for _, _, weights in transitions)
    if np.max(np.abs(total.sum(axis=0)-1)) > 1e-12:
        raise ArithmeticError('numerical complete instrument normalization failed')
    omega = energy[:, None]-energy[None, :]
    rows = []
    for value in betas:
        beta = _rational(value, 'diagnostic beta')
        if beta <= 0:
            raise ValueError('diagnostic beta must be positive')
        probability = np.exp(-float(beta)*(energy-energy[0])); probability /= probability.sum()
        joint = total*probability[None, :]
        excess = float(probability@(energy-energy[0]))
        rows.append({'beta': str(beta), 'mean_ground_excess_decimal': excess,
                     'negative_transfer_at_least_quarter_weight_decimal': float(joint[omega <= -.25].sum()),
                     'mean_excess_over_threshold_decimal': min(1., excess/.25),
                     'all_transition_weight_sum_decimal': float(joint.sum()),
                     'polynomial_canonical_trace_error_bound': str(ground.polynomial_bound(beta, -4, 22)['canonical_state_half_trace_error_upper'])})
    lines = {}
    for sign, _, weights in transitions:
        for m in range(208):
            key = (sign, int(occupation[m]), round(float(omega[m, 0]), 10))
            lines[key] = lines.get(key, 0.)+float(weights[m, 0])
    ranked = sorted(({'sign': key[0], 'target_occupation': key[1],
                      'source_energy_transfer_decimal': key[2],
                      'original_ground_probe_weight_decimal': value}
                     for key, value in lines.items() if value > 1e-12),
                    key=lambda row: (-row['original_ground_probe_weight_decimal'], row['sign']))
    return {'arithmetic': 'floating diagnostic, not interval-admitted residues',
            'ground_energy_decimal': float(energy[0]), 'ground_occupation': int(occupation[0]),
            'cooling': rows, 'ground_lines_above_display_threshold': ranked,
            'display_threshold': 1e-12,
            'omitted_display_weight_decimal': max(0., 1-sum(row['original_ground_probe_weight_decimal'] for row in ranked)),
            'complete_208_by_208_transitions_used': True,
            'negative_frequency_weights_discarded': False,
            'physical_particle_masses': []}


def serial(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(key): serial(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serial(item) for item in value]
    return value


def certificate(include_control=True, include_diagnostics=False):
    """Exact source facts; empirical masses and couplings are absent."""
    minima = sector_minima()
    ground = global_ground_interval()
    conjugation = charge_conjugation()
    rows = instrument()
    result = {
        'schema': 'oph.native_ground_spectral_records.v1',
        'source_files': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCE_INPUTS},
        'source_dimension': 208,
        'all_fock_occupation_states': 16,
        'all_quantum_flux_labels': list(range(-6, 7)),
        'reducing_block_sizes': dict(sorted(Counter(len(indices) for indices, _, _ in rational_blocks()).items())),
        'sector_minimum_intervals': minima,
        'global_ground_interval': ground,
        'global_ground_occupation': 2,
        'adjacent_charge_minimum_transfer_intervals': charged_thresholds(),
        'charge_conjugation': conjugation['identities'],
        'instrument': [{key: value for key, value in row.items() if key != 'K'} for row in rows],
        'instrument_effect_sum_identity': True,
        'complete_original_negative_transfer_bound': 'min(1,mean_ground_excess/eta)',
        'resolved_readout_plan_for_mean_excess_one_hundredth': resolved_negative_record_bound(F(1, 100)),
        'same_source_H_and_clock': True,
        'accepted_branch_normalization': False,
        'three_dimensional_spatial_source': False,
        'relativistic_mass_shell_derived': False,
        'observed_particle_identification': False,
    }
    if include_control:
        result['cooled_source_records'] = cooled_control_certificate()
    if include_diagnostics:
        result['spectral_diagnostics'] = spectral_diagnostics()
    return serial(result)
