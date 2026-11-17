"""Exact joint fields/clocks under native feedback on overlapping regions.

The same diagonal C process is evaluated directly by endpoint recursion and
by a source-derived, weighted effective transfer. Every timeout and terminal
population is retained. The finite-cap effective cone is transported: its
matrices are not asserted to be entrywise nonnegative.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import product
from pathlib import Path

import numpy as np
import sympy as sp

from ..observer import native
from ..observer.scalar_response.experiment import canonical, compare_receipt
from ..observer_time_fields import core
from ..native_action import history_action
from ..observer_reconstruction import native_geometry as geometry
from tools import strict_json
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
REGIONS = ((0, 1), (0, 5))
UNION = (0, 1, 5)
SOURCE_PATHS = tuple(sorted(set((*core.SOURCE_PATHS, *history_action.SOURCE_PATHS,
    'model/observer_correspondence/__init__.py',
    'model/observer_correspondence/__main__.py',
    'model/observer_correspondence/construction.py',
    'model/observer_reconstruction/native_geometry.py',
    'Lean/Screen/RegionalContinuity.lean'))))


def _zero():
    return np.zeros(66, dtype=object)


def _mark(seam, current):
    return (tuple(int(bool(set(seam).intersection(region))) for region in REGIONS),
            tuple(current*(int(v == seam[1])-int(v == seam[0])) for v in UNION))


def update_field(field, mark):
    """Only disclosed initial local fields and signed current increments."""
    return tuple(a+b for a, b in zip(field, mark[1], strict=True))


def observer_view(initial, word, observer):
    """Restrict to local initial bits and incident increments only.

    Common history cuts are supplied by the comparison service. Neither the
    hidden source state nor the other observer's nonincident records enter.
    """
    if type(observer) is not int or observer not in (0, 1):
        raise ValueError('observer index zero or one required')
    positions = tuple(UNION.index(port) for port in REGIONS[observer])
    initial_local = tuple(initial[index] for index in positions)
    records = []
    for kind, mark, _ in word:
        delta = tuple(mark[1][index] for index in positions)
        if mark[0][observer]:
            records.append((kind, delta))
        elif any(delta):
            raise ArithmeticError('a supposedly nonincident event changes a local field')
    return initial_local, tuple(records)


def decode_local(view):
    """Incident ticks count zero currents; entry is read from the local field."""
    initial, records = view
    field, entry = initial, None
    for tick, (_, delta) in enumerate(records, 1):
        field = tuple(a+b for a, b in zip(field, delta, strict=True))
        if not initial[1] and field[1] and entry is None:
            entry = tick
    return field, len(records), entry


@lru_cache(maxsize=1)
def source():
    base = core.source()
    states, seams = base.states, base.seams
    index = {state: i for i, state in enumerate(states)}
    remote = np.zeros((66, 66), dtype=object)
    incident = {}
    transitions = []
    feedback = {}
    local_groups = dict(history_action.kraus_numerators('C'))
    checks = 0
    for column, state in enumerate(states):
        outgoing = []
        for seam in seams:
            for coin in (0, 1):
                output, (_, current) = native.native(state, seam, coin)
                row = index[output]
                local_before = state[seam[0]]+2*state[seam[1]]
                local_after = output[seam[0]]+2*output[seam[1]]
                w = local_groups[coin, current, current]
                if w[local_after, local_before] != 1:
                    raise ArithmeticError('native event and C action Kraus branch disagree')
                checks += 1
                is_incident = bool(set(seam).intersection(UNION))
                mark = _mark(seam, current)
                if is_incident:
                    incident.setdefault(mark, np.zeros((66, 66), dtype=object))[row, column] += 1
                    if update_field(tuple(state[v] for v in UNION), mark) != tuple(output[v] for v in UNION):
                        raise ArithmeticError('local current reconstruction failed')
                else:
                    remote[row, column] += 1
                outgoing.append((row, mark if is_incident else None))
                if seam in ((0, 1), (0, 5)):
                    feedback.setdefault(seam, {}).setdefault(mark,
                        np.zeros((66, 66), dtype=object))[row, column] += 1
        transitions.append(tuple(outgoing))
    identity = np.eye(66, dtype=object)
    d = 60*identity-remote
    inverse = sp.Matrix(d.tolist()).inv(method='DM')
    denominator = int(sp.ilcm(*(entry.q for entry in inverse)))
    numerator = np.array([[int(entry*denominator) for entry in row]
                          for row in inverse.tolist()], dtype=object)
    if not np.array_equal(d@numerator, denominator*identity):
        raise ArithmeticError('exact remote resolvent inverse failed')
    unread = sum(incident.values(), np.zeros((66, 66), dtype=object))
    if (not np.array_equal(unread, unread.T) or not np.array_equal(remote, remote.T)
            or set(unread.sum(axis=0)) != {24} or set(remote.sum(axis=0)) != {36}):
        raise ArithmeticError('the uniform overlapping-observer source changed')
    for maps in feedback.values():
        for matrix in maps.values():
            if not np.array_equal(d@matrix, matrix@d):
                raise ArithmeticError('internal record-selected feedback must commute with remote repair')
    return {'states': states, 'seams': seams, 'remote': remote, 'incident': incident,
            'feedback': feedback, 'transitions': tuple(transitions), 'D': d,
            'inverse_numerator': numerator, 'inverse_denominator': denominator,
            'unread_incident': unread, 'action_branches_checked': checks}


def _cap(cap):
    if type(cap) is not int or not 1 <= cap <= 3:
        raise ValueError('an exact per-interval attempt cap from one through three is required')


@lru_cache(maxsize=3)
def effective_maps(cap=2):
    _cap(cap)
    data = source()
    a, d, e, q = (data[key] for key in
                  ('remote', 'D', 'inverse_numerator', 'inverse_denominator'))
    power = np.eye(66, dtype=object)
    series = np.zeros((66, 66), dtype=object)
    for waiting in range(cap):
        series += 60**(cap-1-waiting)*power
        power = a@power
    result = {}
    for mark, b in data['incident'].items():
        numerator = e@(b@series)@d
        result[mark] = numerator
    # D commutes with the timeout A^cap. Its conjugate is therefore itself.
    return result, q*60**cap, power, 60**cap


def _apply(matrix, state):
    return matrix@state if any(state) else _zero()


def direct_interval(state, cap):
    """Independent endpoint recursion through all remote waits, not K matrices."""
    _cap(cap)
    data = source()
    outputs = {mark: _zero() for mark in data['incident']}
    waiting = state
    for step in range(cap):
        next_waiting = _zero()
        multiplier = 60**(cap-step-1)
        for column, count in enumerate(waiting):
            if not count:
                continue
            for row, mark in data['transitions'][column]:
                if mark is None:
                    next_waiting[row] += count
                else:
                    outputs[mark][row] += count*multiplier
        waiting = next_waiting
    return outputs, 60**cap, waiting, 60**cap


def _interval(state, cap, effective):
    if not effective:
        return direct_interval(state, cap)
    maps, denominator, timeout, timeout_denominator = effective_maps(cap)
    return ({mark: _apply(matrix, state) for mark, matrix in maps.items()}, denominator,
            _apply(timeout, state), timeout_denominator)


def _feedback(state, seam, effective):
    data = source()
    if effective:
        return {mark: _apply(matrix, state) for mark, matrix in data['feedback'][seam].items()}
    result = {mark: _zero() for mark in data['feedback'][seam]}
    index = {value: i for i, value in enumerate(data['states'])}
    for column, count in enumerate(state):
        if count:
            for coin in (0, 1):
                output, (_, current) = native.native(data['states'][column], seam, coin)
                result[_mark(seam, current)][index[output]] += count
    return result


def joint_law(*, effective=False, cap=2, feedback=True):
    """All declared output words, zero outcomes, timeouts and terminal states.

    Uniform sector preparation discloses only the three union bits. Two
    scheduled union incidents are requested. After the first, feedback acts
    once on (0,1) if n0==1, otherwise (0,5); its current and both ticks remain.
    Invalid reconstructed bits occur only in impossible zero-weight menus.
    """
    _cap(cap)
    if type(effective) is not bool or type(feedback) is not bool:
        raise ValueError('explicit Boolean representation and feedback choices required')
    data = source()
    result = {}
    for initial in product((0, 1), repeat=3):
        prepared = np.array([int(tuple(state[v] for v in UNION) == initial)
                             for state in data['states']], dtype=object)
        # G=D/24 fixes each initial disclosed-label vector: no reset or lift.
        if not np.array_equal(data['D']@prepared, 24*prepared):
            raise ArithmeticError('preparation does not have the stated stationary boundary')
        first, d1, timeout, dt1 = _interval(prepared, cap, effective)
        result[initial, (), 'timeout_first'] = timeout, 66*dt1
        for mark, state in first.items():
            field = update_field(initial, mark)
            event = ('scheduled', mark, field)
            seam = (0, 1) if field[0] == 1 else (0, 5)
            following = (_feedback(state, seam, effective) if feedback
                         else {None: state})
            for control_mark, controlled in following.items():
                current_field = update_field(field, control_mark) if feedback else field
                prefix = (event, ('feedback', control_mark, current_field)) if feedback else (event,)
                control_denominator = 2 if feedback else 1
                second, d2, timeout, dt2 = _interval(controlled, cap, effective)
                result[initial, prefix, 'timeout_second'] = (
                    timeout, 66*d1*control_denominator*dt2)
                for last_mark, terminal in second.items():
                    last_field = update_field(current_field, last_mark)
                    word = prefix+(('scheduled', last_mark, last_field),)
                    result[initial, word, 'complete'] = (
                        terminal, 66*d1*control_denominator*d2)
    return result


def compare_laws(native_law, effective_law):
    if native_law.keys() != effective_law.keys():
        raise ArithmeticError('the two representations have different complete record menus')
    data = source()
    total, timeout, positive, zeros = F(0), F(0), 0, 0
    digest_rows = []
    for record, (left, ld) in native_law.items():
        right, rd = effective_law[record]
        decoded = data['D']@right
        if not np.array_equal(24*rd*left, ld*decoded):
            raise ArithmeticError('decoded terminal population differs on a complete record')
        probability = F(sum(left), ld)
        if probability != F(sum(right), rd) or probability < 0:
            raise ArithmeticError('complete joint-record probability mismatch')
        total += probability
        timeout += probability if record[2].startswith('timeout') else 0
        positive += probability > 0
        zeros += probability == 0
        if probability:
            for cursor, (_, _, field) in enumerate(record[1], 1):
                if any(bit not in (0, 1) for bit in field):
                    raise ArithmeticError('a positive-weight current word has invalid local fields')
                views = [decode_local(observer_view(record[0], record[1][:cursor], observer))
                         for observer in (0, 1)]
                if (views[0][0] != (field[0], field[1])
                        or views[1][0] != (field[0], field[2])
                        or views[0][0][0] != views[1][0][0]):
                    raise ArithmeticError('observer restrictions fail common-prefix field agreement')
            if record[1]:
                final_field = record[1][-1][2]
                for count, state in zip(left, data['states'], strict=True):
                    if count and tuple(state[v] for v in UNION) != final_field:
                        raise ArithmeticError('readback and continuing source disagree')
        digest_rows.append([record, str(probability),
            hashlib.sha256(canonical([str(F(value, ld)) for value in left]).encode()).hexdigest()])
    if total != 1:
        raise ArithmeticError('complete experiment lost probability mass')
    return {'record_menu_size': len(native_law), 'positive_records': int(positive),
            'zero_records_retained': int(zeros), 'total_probability': str(total),
            'timeout_probability': str(timeout), 'terminal_population_equality_exact': True,
            'local_restrictions_and_overlap_equal_at_every_positive_prefix': True,
            'complete_law_sha256': hashlib.sha256(canonical(digest_rows).encode()).hexdigest()}


def geometry_and_action():
    """Same source's spatial dictionary and local action, with fixed scope."""
    data = source()
    laplacian = [[0]*12 for _ in range(12)]
    for a, b in data['seams']:
        laplacian[a][a] += 1; laplacian[b][b] += 1
        laplacian[a][b] -= 1; laplacian[b][a] -= 1
    checks, projector = geometry._band_certificate(laplacian)
    if not all(checks.values()):
        raise ArithmeticError('same-source spatial band certificate failed')
    groups, w, d = history_action.completion_blocks('C')
    if (not np.array_equal(w.T@w, 2*np.eye(4, dtype=object))
            or not np.array_equal(d@d, d) or any((d@w).flat)):
        raise ArithmeticError('source action dilation identities failed')
    # The union is a native triangle; its sum of three seam exchanges has
    # eigenvalues 3 or 6 in integer coin-count units, hence B >= 3I/60.
    tri = np.zeros((66, 66), dtype=object)
    for maps in [data['feedback'][(0, 1)], data['feedback'][(0, 5)]]:
        tri += sum(maps.values(), np.zeros((66, 66), dtype=object))
    index = {state: i for i, state in enumerate(data['states'])}
    for column, state in enumerate(data['states']):
        for coin in (0, 1):
            output, _ = native.native(state, (1, 5), coin)
            tri[index[output], column] += 1
    identity = np.eye(66, dtype=object)
    if any(((tri-3*identity)@(tri-6*identity)).flat):
        raise ArithmeticError('source triangle lower bound failed')
    # Every remaining seam is (I+SWAP)/2, so PSD; its complement is PSD too.
    for a, b in data['seams']:
        swap = np.zeros((66, 66), dtype=object)
        actual = np.zeros((66, 66), dtype=object)
        for column, state in enumerate(data['states']):
            target = list(state); target[a], target[b] = target[b], target[a]
            swap[index[tuple(target)], column] = 1
            for coin in (0, 1):
                output, _ = native.native(state, (a, b), coin)
                actual[index[output], column] += 1
        if (not np.array_equal(swap@swap, identity) or not np.array_equal(swap, swap.T)
                or not np.array_equal(actual, identity+swap)):
            raise ArithmeticError('native exchange positivity certificate failed')
    return {'spatial_readout': 'selected rank-three slow response Gram, restricted to ports 0,1,5',
            'local_gram_qsqrt5': geometry._encode_q5([
                [geometry._qscale(projector[a][b], 4) for b in UNION] for a in UNION]),
            'spatial_coordinates_reconstructed_from_this_finite_word_sample': False,
            'spatial_dictionary_uses_same_supplied_source_incidence': True,
            'source_C_action_branches_checked': data['action_branches_checked'],
            'local_Kraus_count': len(groups), 'dilation_dimension': int(d.shape[0]+4),
            'all_actual_seam_count_matrices_equal_I_plus_swap': True,
            'action': 'A0=pi(I-U)/2; exp(-i A0)=U at native event boundaries',
            'unread_completed_transfer': 'F=D^-1 B, metric G=D/(2/5)',
            'positive_self_adjoint_in_G': True,
            'effective_spectral_lower_bound': '1/20',
            'effective_spectral_upper_bound': '1',
            'Euclidean_action': 'H=-log(F) in the G inner product, per completed union incident tick',
            'native_repairs_equal_unitary_continuation_of_H': False}


def local_entry_records(law):
    """Entry of occupation at port 5, derived from the same local clock word.

    This is current transport in a stationary preparation, not a separate
    prepared signal or a physical propagation-speed assertion.
    """
    masses = defaultdict(F)
    for (initial, word, status), (state, denominator) in law.items():
        probability = F(sum(state), denominator)
        if not probability:
            continue
        _, _, entry = decode_local(observer_view(initial, word, 1))
        if initial[2]:
            key = 'initially_occupied'
        elif entry is not None:
            key = f'first_entry_at_B_incident_tick_{entry}'
        else:
            key = 'no_entry_timeout' if status.startswith('timeout') else 'no_entry_complete_word'
        masses[key] += probability
    if sum(masses.values()) != 1:
        raise ArithmeticError('local entry readback dropped a censored or initial outcome')
    return {key: str(value) for key, value in sorted(masses.items())}


def clock_control():
    """Same single-clock marginals do not specify the joint clock law."""
    actual = {(0, 0): F(0), (1, 0): F(1, 4), (0, 1): F(1, 4), (1, 1): F(1, 2)}
    factorized = {(a, b): (F(3, 4) if a else F(1, 4))*(F(3, 4) if b else F(1, 4))
                  for a, b in product((0, 1), repeat=2)}
    data = source()
    measured = defaultdict(F)
    for mark, matrix in data['incident'].items():
        measured[mark[0]] += F(sum(matrix.flat), 66*24)
    if any(measured[key] != value for key, value in actual.items()):
        raise ArithmeticError('joint local clock marks changed')
    variation = sum(abs(actual[key]-factorized[key]) for key in actual)/2
    return {'conditional_on_first_union_incident': [[list(key), str(value)] for key, value in actual.items()],
            'both_local_marginals': '3/4', 'joint_total_variation_from_independent_clocks': str(variation),
            'first_interval_clock_output_TV_including_timeout': str(F(16, 25)*variation),
            'independent_clocks_preserve_each_marginal_but_change_joint_law': True}


def run():
    cap = 2
    native_law = joint_law(cap=cap)
    effective_law = joint_law(cap=cap, effective=True)
    comparison = compare_laws(native_law, effective_law)
    if comparison['timeout_probability'] != '369/625':
        raise ArithmeticError('all timeout prefixes must retain their exact mass')
    maps, _, _, _ = effective_maps(cap)
    return {'schema': 'oph.observer_correspondence.v1',
            'source_sha256': bound_pins({path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                              for path in SOURCE_PATHS}),
            'contract': {'law': 'C, diagonal binary total two', 'states': 66,
                'regions': [list(region) for region in REGIONS], 'union_region': list(UNION),
                'scheduler': 'uniform independent native seam and fair coin',
                'preparation': 'uniform sector with three local initial bits disclosed',
                'scheduled_union_ticks': 2, 'attempt_cap_per_interval': cap,
                'feedback': 'one native repair after first event; (0,1) if reconstructed n0=1, else (0,5)',
                'feedback_is_an_explicit_tick_for_both_observers': True,
                'timeout_service_and_joint_record_comparison_supplied': True},
            'joint_correspondence': comparison,
            'finite_cap_representation': {
                'A': 'remote count matrix / 60', 'D': 'I-A (normalized, not the integer code matrix)',
                'G': 'D/(2/5)', 'native_state': 'p=G e',
                'effective_capped_mark': 'G^-1 B_m sum_(r=0)^(cap-1) A^r G',
                'effective_timeout': 'A^cap', 'cone': 'G^-1 R_+^66',
                'negative_effective_matrix_entries': sum(int(value < 0) for matrix in maps.values() for value in matrix.flat),
                'entrywise_positive_capped_maps_claimed': False,
                'each_marked_map_is_symmetric_or_PSD_claimed': False},
            'geometry_and_action': geometry_and_action(), 'clock_negative_control': clock_control(),
            'local_current_entry_records': local_entry_records(native_law),
            'acceptance': {'exact_same_protocol_joint_diagonal_observer_law': True,
                'arbitrary_coherent_state_or_reference_certificate': False,
                'same_history_observer_relative_QFT_established': False,
                'Lorentzian_spacetime_or_physical_units_established': False,
                'observed_Universe_identification': False, 'federation_required': False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    destination = parser.add_mutually_exclusive_group()
    destination.add_argument('--check', type=Path)
    destination.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.output and args.output.exists():
        parser.error('--output requires a new destination')
    receipt = run()
    if args.check:
        expected = strict_json.loads(args.check.read_text())
        compare_receipt(expected, receipt)
        print('observer correspondence receipt matches')
    elif args.output:
        with args.output.open('x') as stream:
            stream.write(canonical(receipt))
    else:
        print(canonical(receipt))
