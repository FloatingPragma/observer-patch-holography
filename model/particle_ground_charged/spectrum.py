"""Integer-admitted charged spectra and weights of the complete-source ground.

The numerical diagonalization only proposes a complete eigenframe. Exact
integer residual and Gram bounds certify the spectrum and full projectors.
The actual prepared-ground error is propagated without normalizing a branch.
"""
from __future__ import annotations

from collections import deque
from fractions import Fraction as F
from functools import lru_cache
from math import lcm
import json

import numpy as np
from scipy import sparse
from scipy.linalg import eigh
import sympy as s

from model.particle_scalar_ground import ground


BITS = 27
H_NORM_UPPER = 9
SOURCE_INPUTS = (
    'model/particle_scalar_ground/ground.py',
    'model/particle_scalar_ground/PROOF.md',
    'model/particle_calibration_selection/phase_attenuation.py',
    'model/particle_observer_outputs/CHARGED_HIGGS_RECORDS.md',
    'model/particle_observer_outputs/READOUT_CONTRACT.md',
    'model/standard_model_source_map/FERMION_ATTACHMENT.md',
    'model/particle_candidate_atlas/INTERACTING_CHIRAL_JOIN.md',
)


@lru_cache(maxsize=None)
def literal_vertex(scalar, mask):
    return tuple(ground.source.vertex_action(s.diag(1, 2, 3), scalar, mask).items())


@lru_cache(maxsize=1)
def ground_input():
    trial = ground.singlet_trial()
    denominator = lcm(*(value.denominator for value in trial['source_state'].values()))
    u = {(0, mask): int(value*denominator) for mask, value in trial['source_state'].items()}
    image = {(scalar, mask): int(value*denominator) for (scalar, mask), value in trial['source_image'].items()}
    winning = ground.all_active_blocks()['winning']
    exact_empty_support = {(0, ground.original_column(ground.WINNING_TYPES, winning['states'][index])[0])
                           for index in winning['indices']}
    if set(u) != exact_empty_support:
        raise ArithmeticError('the charged trial support differs from the exact Perron ground support')
    possible_image = {destination for state in u for destination, _ in literal_vertex(*state)}
    if set(image) != possible_image:
        raise ArithmeticError('a retained source image column vanished from the charged preparation support')
    n = sum(value*value for value in u.values())
    m = sum(value*value for value in image.values())
    r = F(m, n)
    if r != trial['Gram_Rayleigh_quotient']:
        raise ArithmeticError('the charged program changed its certified preparation')
    delta = ground.iv.sqrt_rational(1+4*r, 28)
    alpha = (1+delta.reciprocal())/2
    return {'u': u, 'Au': image, 'norm': n, 'image_norm': m,
            'r': r, 'delta': delta, 'alpha': alpha, 'beta': 1-alpha,
            'ground_trace_error': ground.complete_ground_records()['pure_preparation_trace_distance_upper']}


def probe(state, mode, addition):
    if type(mode) is not int or not 0 <= mode < 18 or type(addition) is not bool:
        raise ValueError('one original CAR mode and an addition/removal flag are required')
    result = {}
    for (scalar, mask), value in state.items():
        after = (ground.source._create if addition else ground.source._annihilate)(mask, mode)
        if after is not None:
            target, sign = after
            result[scalar, target] = value*sign
    return result


@lru_cache(maxsize=36)
def source_block(mode, addition):
    """Exact reducing closure of the complete probe vector, with its complement."""
    initial = ground_input()
    x = probe(initial['u'], mode, addition)
    y = probe(initial['Au'], mode, addition)
    seen = set(x) | set(y); pending = deque(seen)
    while pending:
        state = pending.popleft()
        for destination, _ in literal_vertex(*state):
            if destination not in seen:
                seen.add(destination); pending.append(destination)
    basis = tuple(sorted(seen)); indices = {state: i for i, state in enumerate(basis)}
    rows, columns, values = [], [], []
    for column, state in enumerate(basis):
        if state[0]:
            rows.append(column); columns.append(column); values.append(1)
        for destination, coefficient in literal_vertex(*state):
            if destination not in indices:
                raise ArithmeticError('the charged source subspace is not invariant')
            if not coefficient.is_Integer:
                raise ArithmeticError('the original finite action is not the declared integer fixture')
            rows.append(indices[destination]); columns.append(column); values.append(int(coefficient))
    h = sparse.csr_matrix((values, (rows, columns)), shape=(len(basis), len(basis)), dtype=np.int64)
    if (h-h.T).nnz or {mask.bit_count()-2*(scalar != 0) for scalar, mask in basis} != {6+(1 if addition else -1)}:
        raise ArithmeticError('the actual charged block changed source Hermiticity or charge')
    return {'mode': mode, 'addition': addition, 'basis': basis, 'H': h,
            'x': tuple(x.get(state, 0) for state in basis),
            'y': tuple(y.get(state, 0) for state in basis),
            'complete_source_dimension': 786432,
            'orthogonal_complement_dimension': 786432-len(basis)}


def exact_proposal(h, bits=BITS):
    """Untrusted numerical proposal, quantized before any admission check."""
    values, vectors = eigh(h.toarray().astype(float))
    scale = 1 << bits
    if not np.isfinite(values).all() or not np.isfinite(vectors).all():
        raise ArithmeticError('the numerical proposal is not finite')
    return np.rint(vectors*scale).astype(np.int64), np.rint(values*scale).astype(np.int64), scale


def _int_bound(values):
    return max((abs(int(value)) for value in np.asarray(values).flat), default=0)


def admit_frame(h, q, energies, scale):
    """Weyl and polar bounds from complete, exact integer matrix identities."""
    n = h.shape[0]
    if (type(h) is not sparse.csr_matrix or h.dtype != np.dtype(np.int64) or
            n == 0 or h.shape != (n, n) or type(q) is not np.ndarray or q.dtype != np.dtype(np.int64) or
            q.shape != (n, n) or type(energies) is not np.ndarray or
            energies.dtype != np.dtype(np.int64) or energies.shape != (n,) or
            type(scale) is not int or scale <= 0):
        raise ValueError('the complete integer source and proposal domain is required')
    pointers = tuple(map(int, h.indptr))
    if (len(pointers) != n+1 or pointers[0] != 0 or pointers[-1] != len(h.data) or
            len(h.indices) != len(h.data) or any(a > b for a, b in zip(pointers, pointers[1:]))):
        raise ValueError('complete canonical CSR storage is required')
    entries = {}
    row_sums = []
    for row in range(n):
        columns = tuple(int(h.indices[k]) for k in range(pointers[row], pointers[row+1]))
        if any(not 0 <= column < n for column in columns) or any(a >= b for a, b in zip(columns, columns[1:])):
            raise ValueError('canonical source rows without duplicate addresses are required')
        row_sums.append(sum(abs(int(h.data[k])) for k in range(pointers[row], pointers[row+1])))
        entries.update({(row, column): int(h.data[k]) for column, k in
                        zip(columns, range(pointers[row], pointers[row+1]), strict=True)})
    if (any(entries.get((column, row), 0) != value for (row, column), value in entries.items()) or
            any(int(a) > int(b) for a, b in zip(energies, energies[1:]))):
        raise ValueError('a Hermitian source and sorted proposal energies are required')
    maxq = _int_bound(q); maxhrow = max(row_sums)
    maximum = np.iinfo(np.int64).max
    if maxhrow*maxq > maximum or n*maxq*maxq+scale*scale > maximum:
        raise OverflowError('exact integer source/Gram multiplication could overflow')
    hq = h@q
    if _int_bound(hq)*scale+maxq*_int_bound(energies) > maximum:
        raise OverflowError('exact source residual subtraction could overflow')
    gram_error = q.T@q-scale*scale*np.eye(n, dtype=np.int64)
    residual = hq*scale-q*energies[np.newaxis, :]
    eps = F(max(sum(abs(int(x)) for x in row) for row in gram_error), scale*scale)
    if eps >= F(1, 100):
        raise ArithmeticError('the proposed complete frame is too far from orthogonal')
    distance = eps/(2-eps)
    row_bound = max(sum(abs(int(x)) for x in row) for row in residual)
    col_bound = max(sum(abs(int(x)) for x in row) for row in residual.T)
    residual_norm = ground.iv.sqrt_rational(F(row_bound*col_bound, scale**4), 28).hi
    emax = F(_int_bound(energies), scale)
    # This verifier is valid for any admitted integer matrix: its own exact
    # absolute-row bound supplies B. The sharper full-source bound9 is used
    # separately for source-moment errors, after the source binding.
    matrix_norm_upper = maxhrow
    radius = residual_norm+(matrix_norm_upper+emax)*distance
    if radius >= F(1, 100):
        raise ArithmeticError('the exact frame residual does not resolve the stated spectral domain')
    groups = []
    for index, energy in enumerate(energies):
        center = F(int(energy), scale)
        if groups and center-radius <= groups[-1]['upper']:
            groups[-1]['indices'].append(index)
            groups[-1]['upper'] = center+radius
        else:
            groups.append({'indices': [index], 'lower': center-radius, 'upper': center+radius})
    for group in groups:
        indices = group['indices']; inside = set(indices)
        other = [i for i in range(n) if i not in inside]
        if other:
            separation = min(abs(F(int(energies[i])-int(energies[j]), scale))
                             for i in indices for j in other)-radius
            if separation <= 0:
                raise ArithmeticError('the retained complete spectral groups overlap their complement')
            frobenius2 = sum(int(residual[row, column])**2 for row in range(n) for column in indices)
            r_frob = ground.iv.sqrt_rational(F(frobenius2, scale**4), 28).hi
            root_rank = ground.iv.sqrt_rational(len(indices), 28).hi
            projector_error = (r_frob+(matrix_norm_upper+emax)*root_rank*distance)/separation
            projector_error += distance*(2+eps)
        else:
            projector_error = distance*(2+eps)
        group['projector_operator_error'] = min(1+eps, projector_error)
    return {'q': q, 'energy_numerators': energies, 'scale': scale,
            'gram_error_norm_upper': eps, 'polar_frame_distance_upper': distance,
            'eigenvalue_radius': radius, 'groups': groups,
            'exact_matrix_norm_upper': matrix_norm_upper,
            'source_matrix_dimension': n, 'exact_integer_admission': True}


def quadratic_interval(a, b, c):
    """Ground-trial quadratic form from exact rational empty/occupied terms."""
    state = ground_input()
    return (state['alpha']*F(a, state['norm'])
            +state['beta']*F(b, state['image_norm'])
            -2*state['delta'].reciprocal()*F(c, state['norm']))


def original_probe_weight(block):
    return quadratic_interval(sum(v*v for v in block['x']), sum(v*v for v in block['y']), 0)


def matrix_moments(block, maximum=4):
    """Unnormalized, same-H moments; never a normalized accepted branch."""
    h = block['H']; x, y = block['x'], block['y']
    hx, hy = list(x), list(y)
    result = []
    for power in range(maximum+1):
        a = sum(v*w for v, w in zip(x, hx, strict=True))
        b = sum(v*w for v, w in zip(y, hy, strict=True))
        c = sum(v*w for v, w in zip(x, hy, strict=True))
        value = quadratic_interval(a, b, c)
        error = ground_input()['ground_trace_error']*(1 if power == 0 else 2*H_NORM_UPPER**power)
        result.append({'power': power, 'trial_interval': value,
                       'true_ground_interval': value+ground.iv.I(-error, error)})
        hx = [sum(int(h.data[k])*hx[int(h.indices[k])]
                  for k in range(int(h.indptr[i]), int(h.indptr[i+1]))) for i in range(h.shape[0])]
        hy = [sum(int(h.data[k])*hy[int(h.indices[k])]
                  for k in range(int(h.indptr[i]), int(h.indptr[i+1]))) for i in range(h.shape[0])]
    return result


@lru_cache(maxsize=12)
def certified_channel(family, kind, addition):
    if type(family) is not int or not 0 <= family < 3 or kind not in ('left', 'singlet'):
        raise ValueError('one original family and gauge multiplet are required')
    mode = 6*family+(0 if kind == 'left' else 4)
    multiplicity = 4 if kind == 'left' else 2
    block = source_block(mode, addition)
    frame = admit_frame(block['H'], *exact_proposal(block['H']))
    q, scale = frame['q'], frame['scale']; n = q.shape[0]
    xq = [sum(int(q[row, column])*block['x'][row] for row in range(n)) for column in range(n)]
    yq = [sum(int(q[row, column])*block['y'][row] for row in range(n)) for column in range(n)]
    ground_energy = ground.exact_root_and_gap()['global_ground_interval']
    rows = []
    for group in frame['groups']:
        indices = group['indices']
        a = F(sum(xq[i]*xq[i] for i in indices), scale*scale)
        b = F(sum(yq[i]*yq[i] for i in indices), scale*scale)
        c = F(sum(xq[i]*yq[i] for i in indices), scale*scale)
        weight = quadratic_interval(a, b, c)
        error = group['projector_operator_error']+ground_input()['ground_trace_error']
        lo, hi = max(F(0), weight.lo-error), min(F(1), weight.hi+error)
        if lo > hi:
            raise ArithmeticError('the complete spectral probability enclosure is empty')
        rows.append({'rank_in_reducing_source_block': len(indices),
                     'energy_interval': (group['lower'], group['upper']),
                     'excitation_interval': (group['lower']-ground_energy[1], group['upper']-ground_energy[0]),
                     'unscaled_single_mode_weight': (lo, hi),
                     'original_complete_instrument_weight': (lo*multiplicity/18, hi*multiplicity/18),
                     'projector_error_upper': group['projector_operator_error'],
                     'zero_or_small_weights_retained': True})
    weight = original_probe_weight(block)
    error = ground_input()['ground_trace_error']
    return {'family': family, 'source_gauge_multiplet': kind,
            'direction': 'addition' if addition else 'removal',
            'representative_original_mode': mode, 'multiplicity': multiplicity,
            'same_H_Qstar': 7 if addition else 5,
            'complete_source_dimension': 786432,
            'actual_probe_reducing_block_dimension': n,
            'retained_outside_block_dimension': block['orthogonal_complement_dimension'],
            'outside_block_weight_for_this_preparation': 0,
            'trial_original_probe_weight_interval': (weight.lo*multiplicity/18, weight.hi*multiplicity/18),
            'complete_original_probe_weight_interval': (max(F(0), weight.lo-error)*multiplicity/18,
                                                        min(F(1), weight.hi+error)*multiplicity/18),
            'frame_eigenvalue_radius': frame['eigenvalue_radius'],
            'frame_gram_error_norm_upper': frame['gram_error_norm_upper'],
            'all_spectral_groups': rows,
            'original_same_H_moments': matrix_moments(block),
            'dark_or_zero_groups_discarded': False}


def offset_obstruction():
    """Translation-invariant controls for the same prepared Q*=6 reference."""
    rows = []
    for family in range(3):
        for kind in ('left', 'singlet'):
            distributions = {}
            for addition in (False, True):
                channel = certified_channel(family, kind, addition)
                moments = [row['true_ground_interval'] for row in channel['original_same_H_moments']]
                if moments[0].lo <= 0:
                    raise ArithmeticError('a normalized-shape diagnostic has no admitted denominator')
                mean = moments[1]/moments[0]
                variance = moments[2]/moments[0]-mean**2
                distributions[addition] = {'weight': moments[0], 'mean_H': mean, 'variance_H': variance}
            minus, plus = distributions[False], distributions[True]
            weight_gap = plus['weight']-minus['weight']
            variance_gap = plus['variance_H']-minus['variance_H']
            rows.append({'family': family, 'source_gauge_multiplet': kind,
                         'addition_minus_removal_raw_weight': weight_gap,
                         'equal_full_raw_laws_excluded': weight_gap.lo > 0 or weight_gap.hi < 0,
                         'addition_minus_removal_conditional_variance': variance_gap,
                         'equal_translated_conditional_shapes_excluded': variance_gap.lo > 0 or variance_gap.hi < 0,
                         'offset_required_to_equalize_conditional_means': (minus['mean_H']-plus['mean_H'])/2})
    if not any(row['equal_translated_conditional_shapes_excluded'] for row in rows):
        raise ArithmeticError('the source moments do not certify the stated spectral-shape obstruction')
    allowed = ground.exact_root_and_gap()['full_ground_gap_lower']/12
    return {'declared_offset': 'H_epsilon=H+epsilon Qstar, with unchanged pair cost and Y',
            'same_reference_addition_excitation_shift': '+epsilon',
            'same_reference_removal_excitation_shift': '-epsilon',
            'original_weights_and_conditional_variances_unchanged': True,
            'conditional_shape_diagnostics_replace_original_histograms': False,
            'Qstar6_stays_unique_global_ground_if_abs_epsilon_less_than': allowed,
            'outside_that_window_ground_must_be_recomputed': True,
            'all_offset_values_with_the_same_reference_match_complete_conjugate_laws': False,
            'physical_CPT_violation_claimed': False,
            'scope': 'these finite charged probes about this fixed singlet reference; no Lorentz vacuum or physical charge-conjugation identification',
            'channel_controls': rows}


def certificate():
    channels = [certified_channel(family, kind, addition) for family in range(3)
                for kind in ('left', 'singlet') for addition in (False, True)]
    state = ground_input()
    sum_weights = sum((ground.iv.I(*channel['complete_original_probe_weight_interval']) for channel in channels), ground.iv.I(0))
    if not sum_weights.lo <= 1 <= sum_weights.hi:
        raise ArithmeticError('the original complete CAR instrument lost probability')
    return serial({'schema': 'oph.global_ground_charged_readout.v1',
        'source_inputs': list(SOURCE_INPUTS),
        'same_source_H': 'W(diag(1,2,3))+W^dagger+NH; epsilonF=0',
        'prepared_reference': 'certified full 786432-state Qstar6 gauge/spin-singlet ground',
        'complete_probe': 'Kplus=sum_i c_i^dagger tensor|i* > /sqrt18; Kminus=sum_i c_i tensor|i> /sqrt18',
        'complete_probe_effects': ['I-NF/18', 'NF/18'],
        'probe_reference_compensates_gauge_representation': True,
        'retained_pointer_labels': 'addition/removal, family, full gauge multiplet, every source spectral group and all other energy outcomes',
        'reducing_block_projection_is_a_computational_reduction_not_an_extra_gauge_pointer': True,
        'individual_gauge_component_projectors_assumed_physical': False,
        'common_clock': 'exp(-itH) from the same full source',
        'initial_state_trace_error_upper': state['ground_trace_error'],
        'floating_eigenvalues_accepted_as_proof': False,
        'full_original_instrument_probability_interval': (sum_weights.lo, sum_weights.hi),
        'channels': channels,
        'addition_removal_offset_obstruction': offset_obstruction(),
        'conditional_branch_histograms_substituted': False,
        'all_quantum_references_retained': True,
        'physical_particle_assignments': [], 'physical_mass_predictions': [],
        'measured_particle_inputs': []})


def serial(value):
    if isinstance(value, dict):
        return {str(k): serial(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [serial(v) for v in value]
    if isinstance(value, ground.iv.I):
        return value.strings()
    if isinstance(value, (F, s.Basic)):
        return str(value)
    return value


if __name__ == '__main__':
    print(json.dumps(certificate(), indent=2, sort_keys=True))
