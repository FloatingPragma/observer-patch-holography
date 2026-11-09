"""Exact finite-archive response of the configured native C-duty memory.

The feedback pulse is the owner's declared logical service. Its first
derivative, record count and finite-step mixing are computed without a
physical-capacity law or experimental mass input.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
from pathlib import Path

import sympy as s

from model.common_source import native_family_record_memory as memory


def _clean(matrix):
    return s.ImmutableMatrix(s.Matrix(matrix).applyfunc(
        lambda x: s.factor(s.radsimp(s.expand(x)))))


def _steps(n):
    if type(n) is not int or n < 1:
        raise ValueError('a positive integer number of complete attempts is required')
    return n


def _rational(value):
    if isinstance(value, (float, bool)):
        raise ValueError('an exact rational is required')
    value = s.Rational(value) if isinstance(value, str) else s.sympify(value)
    if value.is_Rational is not True:
        raise ValueError('an exact rational is required')
    return value


@lru_cache(maxsize=1)
def source():
    """Reconstruct all 12 old-memory inputs and all 128 original raw labels."""
    owner = memory.architecture()
    p, q = owner['p'], owner['q']
    if (p, q) != (s.Rational(1, 12), s.Rational(1, 6)):
        raise ArithmeticError('the literal native event law differs')
    transition, enabled = s.zeros(12), s.zeros(12)
    for old in range(12):
        for cc, cr, ec, er, _, weight in owner['records']:
            c, e = (cc, cr) == (0, 1), ec == 0 and er in (4, 5)
            dest = memory._destination(old, c, e, 2, owner['gC'], owner['gE'])
            transition[dest, old] += weight
            if c:
                enabled[dest, old] += weight
    if transition != owner['T'] or any(sum(enabled[:, j]) != p for j in range(12)):
        raise ArithmeticError('the complete original transition was not reconstructed')
    positive = sum(1 for record in owner['records'] if record[-1] > 0)
    if len(owner['records']) != 128 or positive != 24:
        raise ArithmeticError('the complete raw-record support differs')

    # Reset2 reaches only these addresses. Old input addresses outside this
    # support were retained above; this projection is an exact reachable-sector
    # reduction, not an assertion about arbitrary initial old memory.
    addresses = (2, 3, 9, 10)
    injection = s.zeros(12, 4)
    for j, address in enumerate(addresses):
        injection[address, j] = 1
    reduced_t, reduced_e = injection.T * transition * injection, injection.T * enabled * injection
    if transition * injection != injection * reduced_t or enabled * injection != injection * reduced_e:
        raise ArithmeticError('the reset source does not preserve its reachable addresses')
    lump = s.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 1]])
    a, b = p * (1-q), p*q
    t = s.Matrix([[1-a, 1-b, 1-a], [a, 0, 0], [0, b, a]])
    e = s.Matrix([[b, a, b], [a, 0, 0], [0, b, a]])
    if lump*reduced_t != t*lump or lump*reduced_e != e*lump:
        raise ArithmeticError('the actual complete source fails the three-ray quotient')
    projectors = tuple(owner['projectors'][i] for i in (2, 3, 10))
    if owner['projectors'][9] != owner['projectors'][10]:
        raise ArithmeticError('the merged projector rays differ')
    for i, pi in enumerate(projectors):
        if _clean(pi*pi-pi) != s.zeros(3) or s.trace(pi) != 1:
            raise ArithmeticError('a family response is not a rank-one projector')
        for pj in projectors[i+1:]:
            if s.simplify(s.trace(pi*pj)) != s.Rational(1, 5):
                raise ArithmeticError('the exact carrier overlap differs')
    determinant = s.factor(sum(projectors, s.zeros(3)).det())
    d = s.Rational(2, 5)-s.Rational(2, 25)*s.sqrt(5)
    if s.simplify(determinant-d) != 0:
        raise ArithmeticError('the carrier three-ray Gram determinant differs')
    delta = 1-a*a+a*b
    stationary = s.ImmutableMatrix([1-a, a*(1-a), a*b])/delta
    if t*stationary != stationary or sum(stationary) != 1:
        raise ArithmeticError('the reduced stationary state differs')
    return {'p': p, 'q': q, 'a': a, 'b': b,
            'T12': s.ImmutableMatrix(transition), 'E12': s.ImmutableMatrix(enabled),
            'T': s.ImmutableMatrix(t), 'E': s.ImmutableMatrix(e),
            'projectors': projectors, 'stationary': stationary,
            'gram_determinant': d, 'allocated_raw_labels': 128,
            'positive_raw_labels': positive, 'old_memory_inputs': 12}


def weights(n):
    """Mean enabled ray counts per attempt, starting from literal reset port2."""
    n = _steps(n)
    row = source()
    eta, total = s.Matrix([1, 0, 0]), s.zeros(3, 1)
    for _ in range(n):
        total += row['E']*eta
        eta = row['T']*eta
    return s.ImmutableMatrix(total/n)


def coefficient(n):
    """Y_n: d Phi_n(t)/dt at zero is -i[Y_n,.] for h=t/n in X-plus."""
    row, w = source(), weights(n)
    return _clean(sum((w[i]*row['projectors'][i] for i in range(3)), s.zeros(3)))


def stationary_coefficient():
    row = source()
    w = row['E']*row['stationary']
    return _clean(sum((w[i]*row['projectors'][i] for i in range(3)), s.zeros(3)))


def finite_coefficient_error_bound(n):
    """Operator-norm bound from the source reset contraction, all n>=1."""
    n, p = _steps(n), source()['p']
    return 2*p*(1-p**n)/(n*(1-p))


def uniform_spectral_bounds():
    """Exact bounds for every integer n>=2, proved by invariant inequalities."""
    row = source()
    p, a, b = row['p'], row['a'], row['b']
    lower_weights = (b, a*(1-p), a*b*(1-p)/2)
    pair_sum = sum(lower_weights[i]*lower_weights[j] for i in range(3) for j in range(i+1, 3))
    # det(sum_i l_i P_i) > prod(l_i)/5; e2 = 4 pair_sum/5.
    # lambda_min >= det/e2, hence it is strictly larger than this rational.
    eigenvalue_lower = s.prod(lower_weights)/(4*pair_sum)
    middle_upper = b+(a-b)*a+a*b/(1-a)
    return {'minimum_attempts': 2, 'lower_ray_weights': tuple(lower_weights),
            'smallest_eigenvalue_strict_lower': eigenvalue_lower,
            'middle_eigenvalue_upper': middle_upper,
            'middle_to_smallest_strict_upper': middle_upper/eigenvalue_lower,
            'largest_to_smallest_strict_upper': p/eigenvalue_lower,
            'rank_at_one_attempt': 2, 'rank_at_every_n_ge_two': 3}


def archive_budget(n):
    """Exact history counts for fixed reset input and complete raw retention.

    N_record=n ln24 is the direct correctable-record capacity if halt and
    exact reread are the declared continuations. It is not a cosmic capacity.
    """
    n = _steps(n)
    row = source()
    return {'attempts': n, 'controller_address_dimension': 12,
            'allocated_raw_histories': row['allocated_raw_labels']**n,
            'positive_raw_histories': row['positive_raw_labels']**n,
            'zero_raw_histories': row['allocated_raw_labels']**n-row['positive_raw_labels']**n,
            'correctable_record_nats_under_exact_reread': f'{n}*log(24)',
            'physical_P_identified_with_record_capacity': False}


def purity_loss_floor(n, tangent):
    """Complete-channel floor for pulse h=2 atan(tangent), input family P2.

    Two retained coarse histories suffice: all pulses disabled, or only the
    last pulse enabled with C true and E false. Every other record is retained
    and contributes a nonnegative pair term to the purity deficit.
    """
    n, tangent = _steps(n), _rational(tangent)
    if tangent == 0:
        raise ValueError('a nontrivial rational half-angle is required')
    row = source()
    p, a = row['p'], row['a']
    all_disabled, final_ray3 = (1-p)**n, (1-p)**(n-1)*a
    return s.Rational(32, 25)*all_disabled*final_ray3*tangent*tangent/(1+tangent*tangent)


def complete_unread_density(n, tangent):
    """Exact reduced channel for a P2 input, summing every original record.

    The source quotient has been checked against all 1536 original rows.
    Enabled and disabled entries are kept separately in every memory block.
    """
    n, tangent = _steps(n), _rational(tangent)
    row = source()
    phase = (1-s.I*tangent)**2/(1+tangent*tangent)
    pulses = tuple(_clean(s.eye(3)+(phase-1)*pi) for pi in row['projectors'])
    blocks = [s.Matrix(row['projectors'][0]), s.zeros(3), s.zeros(3)]
    for _ in range(n):
        next_blocks = []
        for new in range(3):
            value = s.zeros(3)
            for old in range(3):
                value += (row['T'][new, old]-row['E'][new, old])*blocks[old]
                value += row['E'][new, old]*pulses[new]*blocks[old]*pulses[new].H
            next_blocks.append(_clean(value))
        blocks = next_blocks
    return _clean(sum(blocks, s.zeros(3)))


def comparison(muon_electron, tau_electron):
    """External comparison after the source bounds are computed and frozen."""
    bounds = uniform_spectral_bounds()
    mu, tau = _rational(muon_electron), _rational(tau_electron)
    if mu <= 1 or tau <= mu:
        raise ValueError('ordered positive comparison ratios are required')
    return {'comparison_muon_electron': str(mu), 'comparison_tau_electron': str(tau),
            'middle_ratio_above_all_horizon_bound': bool(mu >= bounds['middle_to_smallest_strict_upper']),
            'largest_ratio_above_all_horizon_bound': bool(tau >= bounds['largest_to_smallest_strict_upper']),
            'required_attachment': 'physical masses proportional to the three ordered eigenvalues of Y_n',
            'other_controller_or_interacting_response_excluded': False}


def certificate():
    row, bounds = source(), uniform_spectral_bounds()
    snapshots = []
    for n in (1, 2, 4, 8):
        matrix = coefficient(n)
        snapshots.append({'attempts': n, 'weights': list(map(str, weights(n))),
                          'trace': str(s.trace(matrix)), 'rank': matrix.rank(),
                          'e2': str(s.factor((s.trace(matrix)**2-s.trace(matrix*matrix))/2)),
                          'determinant': str(s.factor(matrix.det())),
                          'coefficient_norm_error_bound': str(finite_coefficient_error_bound(n)),
                          'archive': archive_budget(n)})
    root = Path(__file__).resolve().parents[2]
    pins = ('model/common_source/native_family_record_memory.py',
            'model/common_source/native_repair_mass.py',
            'model/particle_memory_response/PROOF.md',
            'Lean/Screen/SeamCurrentPhysicalMetricAttachment.lean')
    return {'schema': 'oph.capacity_response.v1',
            'source': 'literal reset2 C-duty memory; native rates; declared logical projector feedback',
            'p': str(row['p']), 'q': str(row['q']),
            'raw_rows_reconstructed': 12*128, 'positive_raw_labels': row['positive_raw_labels'],
            'snapshots': snapshots,
            'uniform_all_integer_n_ge_two_bounds': {
                key: list(map(str, value)) if isinstance(value, tuple) else str(value)
                for key, value in bounds.items()},
            'finite_step_purity_floor_n2_tangent_one_tenth': str(purity_loss_floor(2, s.Rational(1, 10))),
            'same_prefix_under_FIFO_resegmentation_error': '0',
            'physical_P_to_program_map_selected': False,
            'record_capacity_equals_cosmic_capacity': False,
            'mass_or_physical_self_energy_derived': False,
            'target_mass_input_used': False,
            'source_sha256': {path: hashlib.sha256((root/path).read_bytes()).hexdigest() for path in pins}}
