"""Complete correlated-memory weak feedback with the existing flavor means.

The controller and its coherent model clock are explicit supplied wiring.
All error bounds use unhalved trace/diamond distance and original records.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import product
from math import isqrt

import numpy as np
from scipy.linalg import expm
import sympy as s

from model.common_source import native_family_record_memory as memory
from model.common_source import native_gauge_family_records as gates
from model.native_flavor_orbit import algebra as flavor

POLICIES = ('always', 'C', 'E', 'CE')
FLAGS = ((False, False), (False, True), (True, False), (True, True))


def _rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, str, F)):
        raise ValueError('an exact rational is required')
    try:
        return F(value.numerator, value.denominator) if isinstance(value, F) else F(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError('an exact rational is required') from exc


def _factors(p, q):
    p, q = _rational(p), _rational(q)
    if not 0 <= p < 1 or not 0 <= q <= 1:
        raise ValueError('0<=p<1 and 0<=q<=1 are required')
    return p, q


def enabled(policy, c, e):
    if policy not in POLICIES or type(c) is not bool or type(e) is not bool:
        raise ValueError('one declared duty and two boolean events are required')
    return policy == 'always' or policy == 'C' and c or policy == 'E' and e or policy == 'CE' and c and e


def sqrt_upper(value, bits=40):
    """Rational outward square root, without a floating admission boundary."""
    value = _rational(value)
    if value < 0 or type(bits) is not int or not 0 <= bits <= 256:
        raise ValueError('nonnegative value and integer precision0..256 required')
    scale = 1 << bits
    target = (value.numerator*scale*scale+value.denominator-1)//value.denominator
    root = isqrt(target)
    if root*root < target:
        root += 1
    return F(root, scale)


def averaging_bound(steps, time, p, *, feedback_norm=1, stationary_product=False):
    """Reference-safe reduced and full-archive bounds, with any common H0.

    The smaller reduced bound for a stationary product is not used for the
    archive: the archive reveals the initial memory, so its proof uses the
    uniform fixed-memory estimate. Arbitrary initial CQ correlations are kept
    in the comparison instrument. No stationary-product archive is asserted.
    """
    if type(steps) is not int or steps <= 0 or type(stationary_product) is not bool:
        raise ValueError('positive integer steps and explicit product flag required')
    time, norm, p = _rational(time), _rational(feedback_norm), _rational(p)
    if time < 0 or norm < 0 or not 0 <= p < 1:
        raise ValueError('nonnegative time/norm and p<1 required')
    u = norm*time
    transient = 4*p*u/(steps*(1-p))
    quadratic = 4*(1+p)*u*u/(steps*(1-p))
    arbitrary = transient+quadratic
    reduced = quadratic if stationary_product else arbitrary
    return {'reduced_distance_upper': min(F(2), reduced),
            'reduced_uncapped_bound': reduced,
            'uniform_fixed_memory_bound': arbitrary,
            'archive_distance_squared_upper': min(F(4), 2*arbitrary),
            'archive_distance_upper': min(F(2), sqrt_upper(2*arbitrary)),
            'distance_is_unhalved': True,
            'common_H0_norm_enters_bound': False,
            'archive_uses_same_record_law_and_initial_correlations': True}


def steps_for_archive(error, time, p, *, feedback_norm=1):
    error = _rational(error)
    if not 0 < error <= 2:
        raise ValueError('a positive unhalved error at most2 is required')
    base = averaging_bound(1, time, p, feedback_norm=feedback_norm)['uniform_fixed_memory_bound']
    required = 2*base/(error*error)
    return max(1, (required.numerator+required.denominator-1)//required.denominator)


def record_model(p, q, policy='always'):
    """Exact48 coarse rows, retaining every old memory and all four events.

    These are exact sufficient-statistic groups of the declared fine source
    records. The fine labels are restored by a conditional CPTP refinement.
    """
    p, q = _factors(p, q)
    enabled(policy, False, False)
    ps, qs = s.Rational(p.numerator, p.denominator), s.Rational(q.numerator, q.denominator)
    data = memory.architecture_from_factors(ps, qs)
    weights = ((1-p)*(1-q), (1-p)*q, p*(1-q), p*q)
    rows = []
    mean = s.zeros(3)
    transition = s.zeros(12)
    for old in range(12):
        for (c, e), weight in zip(FLAGS, weights, strict=True):
            dest = memory._destination(old, c, e, 2, data['gC'], data['gE'])
            on = enabled(policy, c, e)
            ws = s.Rational(weight.numerator, weight.denominator)
            transition[dest, old] += ws
            if on:
                mean += data['stationary'][old]*ws*data['projectors'][dest]
            rows.append({'old': old, 'C': c, 'E': e, 'destination': dest,
                         'probability': F(weight.numerator, weight.denominator), 'enabled': on})
    mean = mean.applyfunc(s.simplify)
    expected = flavor.coefficient(ps, qs, policy)
    if transition != data['T'] or any(s.simplify(x) != 0 for x in mean-expected):
        raise ArithmeticError('complete memory or native flavor mean mismatch')
    return {'rows': tuple(rows), 'stationary': data['stationary'],
            'projectors': data['projectors'], 'mean': s.ImmutableMatrix(mean),
            'transition': s.ImmutableMatrix(transition), 'policy': policy}


def path_instrument(p, q, policy, *, steps, time, initial_memory, common_hamiltonian=None):
    """Small exhaustive control; every event path, including zero paths.

    Returns work-space Kraus matrices and ideal matrices with the SAME path
    probabilities. It does not allocate a native fine-record tape.
    """
    if type(steps) is not int or not 1 <= steps <= 5:
        raise ValueError('exhaustive control supports one through five steps')
    if type(initial_memory) is not int or not 0 <= initial_memory < 12:
        raise ValueError('all twelve memory inputs are indexed0..11')
    time = _rational(time)
    if time < 0:
        raise ValueError('nonnegative exact time required')
    data = record_model(p, q, policy)
    h0 = np.zeros((6, 6), complex) if common_hamiltonian is None else np.asarray(common_hamiltonian, complex)
    if h0.shape != (6, 6) or not np.isfinite(h0).all() or np.max(abs(h0-h0.conj().T)) > 1e-13:
        raise ValueError('a finite Hermitian6-by6 common Hamiltonian is required')
    X = np.array([[0, 1], [1, 0]], complex)
    generators = [np.kron(X, np.asarray(P, complex)) for P in data['projectors']]
    rows = {(r['old'], r['C'], r['E']): r for r in data['rows']}
    pulse = {(dest, on): expm(-1j*float(time/steps)*(h0+(generators[dest] if on else 0)))
             for dest in range(12) for on in (False, True)}
    ideal = expm(-1j*float(time)*(h0+np.kron(X, np.asarray(data['mean'], complex))))
    result = []
    for path in product(FLAGS, repeat=steps):
        state, probability, U = initial_memory, F(1), np.eye(6, dtype=complex)
        for c, e in path:
            row = rows[state, c, e]
            state = row['destination']
            probability *= row['probability']
            U = pulse[state, row['enabled']]@U
        scale = np.sqrt(float(probability))
        result.append({'path': path, 'destination': state, 'probability': probability,
                       'kraus': scale*U, 'ideal_kraus': scale*ideal})
    return tuple(result)


def pairing_lift(kernel):
    """Finite physical CAR pairing symbol; no dense2^90 Fock allocation.

    H=1/2 sum Delta_ij a_i^dagger a_j^dagger+h.c. has norm
    ||Delta||_*/2=||K||_* for Delta=epsilon tensor conjugate(K).
    The displayed numeric norm is diagnostic; analytic budgets use rank.
    """
    K = np.asarray(kernel, complex)
    if K.shape != (45, 45) or not np.isfinite(K).all() or np.max(abs(K-K.T)) > 1e-13:
        raise ValueError('a finite symmetric45-by45 action kernel is required')
    epsilon = np.array([[0, 1], [-1, 0]], complex)
    return {'creation_kernel': np.kron(epsilon, K.conj()),
            'Fock_operator_norm_diagnostic': float(np.linalg.svd(K, compute_uv=False).sum()),
            'physical_odd_CAR_modes': 90,
            'record_conditioned_controller_is_configured': True}


def physical_budget(steps, time, factors, *, cells=1, higgs_norm=1):
    """Conservative finite-volume budget for three independent sector memories.

    The unit Higgs directed ranks3,3,1 give pairing norms6,6,2 per
    cell. A simultaneous global reset supplies the exact contraction used
    here. Each factor may be a certified upper bound on its C probability,
    since the bound increases with p. No input dominance is inferred by this
    function. This estimate is explicitly not uniform in growing volume.
    """
    if type(cells) is not int or cells <= 0 or not isinstance(factors, (tuple, list)) or len(factors) != 3:
        raise ValueError('positive cells and three sector C probabilities required')
    probabilities = [_factors(p, 0)[0] for p in factors]
    amplitude = _rational(higgs_norm)
    if amplitude < 0:
        raise ValueError('nonnegative configured Higgs amplitude required')
    reset = F(1)
    for p in probabilities:
        reset *= 1-p
    p_joint = 1-reset**cells
    norm = 14*cells*amplitude
    result = averaging_bound(steps, time, p_joint, feedback_norm=norm)
    result.update(joint_C_probability=p_joint, pairing_norm_upper=norm,
                  physical_CAR_modes=90*cells, source_volume_uniform=False,
                  fixed_neutral_Higgs_is_full_unbroken_gauge_symmetry=False)
    return result


def physical_lift_certificate():
    """Replay the source partial-isometry normalization used by L=14."""
    tensors = gates.higgs_tensor_certificate()
    norms = {row['action_channel']: row['neutral_rank'] for row in tensors}
    if norms != {'up': 6, 'down': 6, 'charged_lepton': 2}:
        raise ArithmeticError('source neutral-Higgs pairing rank changed')
    if any(row['directed_nonzero_singular_value_exact'] != '1'
           or row['neutral_nonzero_singular_value_exact'] != '1'
           or not row['dictionary_scale_is_sign_only'] for row in tensors):
        raise ArithmeticError('a sector-dependent source modulus appeared')
    return {'rank_one_family_Fock_pairing_norms_exact': norms,
            'one_cell_pairing_norm_upper_exact': 14,
            'source_Higgs_tensor_normalization_replayed': True,
            'generic_complex_creation_convention': 'epsilon_spin tensor conjugate(K45)',
            'fixed_neutral_Higgs_preserves_electric_charge': True,
            'all_twelve_currents_require_simultaneous_Higgs_transformation': True}


def finite_history_control():
    """Replay complete short histories on a reference-entangled logical input.

    Binary64 controls have a fixed rejection tolerance. The convergence
    theorem and exact probability strings do not depend on their last digits.
    """
    p, q, time = F(1, 12), F(1, 6), F(1, 16)
    h0 = np.diag([20, -7, 3, 5, -11, 9]).astype(complex)
    h0[0, 4], h0[4, 0] = 2j, -2j
    psi = np.arange(1, 13)+1j*np.arange(12, 0, -1)
    psi = psi/np.linalg.norm(psi)
    bounds = averaging_bound(3, time, p)
    probabilities = None
    for policy in POLICIES:
        paths = path_instrument(p, q, policy, steps=3, time=time,
                                initial_memory=0, common_hamiltonian=h0)
        gram = np.zeros((6, 6), complex)
        actual, ideal = np.zeros((12, 12), complex), np.zeros((12, 12), complex)
        archive = 0.
        for row in paths:
            gram += row['kraus'].conj().T@row['kraus']
            v = np.kron(row['kraus'], np.eye(2))@psi
            w = np.kron(row['ideal_kraus'], np.eye(2))@psi
            a, b = np.outer(v, v.conj()), np.outer(w, w.conj())
            actual += a
            ideal += b
            archive += np.abs(np.linalg.eigvalsh(a-b)).sum()
        reduced = np.abs(np.linalg.eigvalsh(actual-ideal)).sum()
        if (np.linalg.norm(gram-np.eye(6), 2) > 1e-12
                or reduced > float(bounds['reduced_distance_upper'])+1e-12
                or archive > float(bounds['archive_distance_upper'])+1e-12
                or archive <= 1e-8):
            raise ArithmeticError('complete noncommuting finite-history control failed')
        row_probabilities = tuple(str(row['probability']) for row in paths)
        if probabilities is not None and probabilities != row_probabilities:
            raise ArithmeticError('duty feedback changed the original archive probabilities')
        probabilities = row_probabilities
    return {'histories_per_policy': 64, 'policies': list(POLICIES),
            'total_complete_coarse_histories': 256,
            'initial_memory': 0, 'time_exact': str(time), 'steps': 3,
            'original_path_probabilities_exact': list(probabilities),
            'all_input_completeness_rejection_threshold': '1/1000000000000',
            'reference_dimension': 2, 'common_H0_noncommuting': True,
            'finite_archive_differs_from_mean_unitary_control': True,
            'native_fine_tape_executed': False}


def certificate():
    p, q = F(1, 12), F(1, 6)
    rows = []
    for policy in POLICIES:
        data = record_model(p, q, policy)
        rows.append({'policy': policy, 'coarse_rows': len(data['rows']),
                     'mean_exact': [[str(x) for x in data['mean'].row(i)] for i in range(3)],
                     'unconditional_trace': str(s.trace(data['mean'])),
                     'equals_existing_native_flavor_coefficient': True})
    def strings(value):
        if isinstance(value, F):
            return str(value)
        if isinstance(value, dict):
            return {k: strings(v) for k, v in value.items()}
        return value
    return {'all_four_duties': rows,
            'exact_baseline_event_factors': [str(p), str(q)],
            'gated_event_factor_domain': '0<=p<=1/12,0<=q<=1/6 from the existing current instrument',
            'native_raw_labels_per_old_memory': 128,
            'gated_fine_labels_per_old_memory': 21632,
            'gated_binary_labels_per_old_memory': 512,
            'all_old_memory_inputs': 12,
            'complete_finite_history_control': finite_history_control(),
            'physical_CAR_lift': physical_lift_certificate(),
            'logical_example_4096_steps': strings(averaging_bound(4096, 1, p)),
            'physical_example_65536_steps': strings(physical_budget(65536, 1, (p, p, p))),
            'physical_example_uses_source_C_probability_upper_1_over_12': True,
            'logical_archive_error_1_over_16_steps': steps_for_archive(F(1, 16), 1, p),
            'finite_step_channel_is_exact_mean_unitary': False,
            'complete_archive_is_discarded_or_replaced_by_stationary_product': False,
            'physical_record_feedback_wiring_is_configured': True,
            'physical_feedback_tape_executed': False,
            'physical_clock_or_vacuum_selected': False,
            'source_growth_uniform_limit_proved': False,
            'flavor_singular_values_modified_by_limit': False,
            'Koide_relation_or_nonzero_CP_generated': False,
            'interacting_particle_poles_identified': False}
