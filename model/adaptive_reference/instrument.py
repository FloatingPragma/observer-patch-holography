"""Exact native attachment, all-input balance, and finite-corner mixing."""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product

import numpy as np

from ..charged_defect import branches
from ..observer import native
from . import controller as c


PAIRS = tuple(combinations(range(3), 2))


@lru_cache(maxsize=2)
def corner(selector):
    if type(selector) is not int or selector not in (0, 1):
        raise ValueError('selector bit required')
    law, seam = c.selector((selector, 0))
    result = []
    for record in c.LABELS:
        a = np.zeros((3, 3), dtype=np.int64)
        if record != ('outside',):
            _, coin, raw = record
            for col, pair in enumerate(PAIRS):
                loads = tuple(int(j in tuple(c.PORTS[k] for k in pair)) for j in range(12))
                updated, (_, current) = native.native(loads, seam, coin)
                emitted = {0: 0, 1: 1, -1: 2}[current]
                if law == 'E':
                    emitted += 4*((loads[seam[0]]+loads[seam[1]]) % 2)
                if raw != emitted:
                    continue
                destination = tuple(j for j, p in enumerate(c.PORTS) if updated[p])
                labels = [c.PORTS[k] for k in pair]
                if current:
                    labels = [seam[1] if x == seam[0] else seam[0] if x == seam[1] else x for x in labels]
                a[PAIRS.index(destination), col] = -1 if labels[0] > labels[1] else 1
        a.flags.writeable = False
        result.append(a)
    return tuple(result)


def attachment_certificate():
    indices = [4*branches.INDEX[tuple(c.PORTS[k] for k in pair)] for pair in PAIRS]
    checked = []
    for selector in (0, 1):
        law, seam = c.selector((selector, 0))
        columns = np.zeros(264, dtype=np.int64)
        rows = columns.copy()
        positive = 0
        for record, small in zip(c.LABELS, corner(selector), strict=True):
            if record == ('outside',):
                if np.any(small):
                    raise ArithmeticError('the invariant corner must have zero outside outcome')
                continue
            _, coin, raw = record
            if law == 'C' and raw >= 4:
                if np.any(small):
                    raise ArithmeticError('padded C record should be zero')
                continue
            full = branches.branch(law, seam, coin, raw).matrix()
            if not np.array_equal(full[indices][:, indices].toarray(), small):
                raise ArithmeticError('independent native endpoint reconstruction differs from the actual charged lift')
            if full[[j for j in range(264) if j not in indices]][:, indices].nnz:
                raise ArithmeticError('selected repair corner leaks under its native branch')
            columns += np.asarray(full.multiply(full).sum(axis=0)).ravel()
            rows += np.asarray(full.multiply(full).sum(axis=1)).ravel()
            positive += bool(full.nnz)
        if np.any(columns != 2) or np.any(rows != 2):
            raise ArithmeticError('full264 charged native instrument is not complete and unital')
        if selector == 1 and any(np.count_nonzero(a) > 1 for a in corner(selector)):
            raise ArithmeticError('E corner outcome ceased to be rank one')
        checked.append({'selector': selector, 'law': law, 'native_seam': list(seam),
                        'full_charged_dimension': 264, 'positive_raw_slots': positive,
                        'input_and_output_column_counts': 2,
                        'independent_signed_corner_columns_bound': True})
    return {'selected_instruments': checked,
            'full_Fock_extension': 'L_native on P plus Q=I-P; Q retained at unit amplitude',
            'complete_alphabet': [list(r) for r in c.LABELS],
            'zero_C_padding_is_not_an_extra_process': True,
            'corner': 'three positions, two identical highest-spin charged defects on ports0,1,5',
            'all_E_corner_records_rank_at_most_one': True,
            'joined_H_evolution_preserves_this_corner': False}


def exact_channels():
    """Twice the CQ channel on4 M3 blocks and its population restriction."""
    channel = np.zeros((36, 36), dtype=np.int64)
    population = np.zeros((12, 12), dtype=np.int64)
    for first, second in product(range(2), repeat=2):
        before = 2*first+second
        for record, op in zip(c.LABELS, corner(first), strict=True):
            after_bits = c.advance((first, second), record)
            after = 2*after_bits[0]+after_bits[1]
            channel[9*after:9*after+9, 9*before:9*before+9] += np.kron(op, op)
            population[3*after:3*after+3, 3*before:3*before+3] += op*op
    return channel, population


def balance_certificate():
    channel, population = exact_channels()
    identity = np.tile(np.eye(3, dtype=np.int64).ravel(), 4)
    if np.any(channel@identity != 2*identity) or np.any(identity@channel != 2*identity):
        raise ArithmeticError('CQ channel lost all-input completeness or unitality')
    if np.any(population.sum(axis=0) != 2) or np.any(population.sum(axis=1) != 2):
        raise ArithmeticError('native joint populations lost double stochasticity')
    return {'CQ_algebra': 'direct_sum_4 M3 for the finite mixing check; full field uses direct_sum_4 M_(2^M)',
            'complete_and_unital_exact_integer_check': True,
            'joint_population_numerator_over2': population.tolist(),
            'general_incoming_row_identity': 'sum_r L_r^(v) L_r^(v)* = I at every output(u,v)',
            'outgoing_column_identity': 'sum_r L_r^(c)* L_r^(c) = I at every input(c,d)',
            'arbitrary_finite_selector_family': '(c,d)->(d+f(c,r) mod q,c) on q^2 classical working states',
            'arbitrary_quantum_field_inputs_and_untouched_references': True,
            'coherent_quantum_controller_source_execution_claimed': False}


def mixing_certificate():
    channel, pop = exact_channels()
    no_e = channel.copy(); no_e[:, 18:] = 0
    no_e_pop = pop.copy(); no_e_pop[:, 6:] = 0
    eighth = np.linalg.matrix_power(channel, 8)
    ranked = eighth-np.linalg.matrix_power(no_e, 8)
    counts = np.linalg.matrix_power(pop, 8)-np.linalg.matrix_power(no_e_pop, 8)
    diagonal = [9*m+4*j for m in range(4) for j in range(3)]
    expected = np.zeros_like(ranked)
    expected[np.ix_(diagonal, diagonal)] = counts
    if not np.array_equal(ranked, expected) or int(counts.min()) != 1:
        raise ArithmeticError('rank-one paths do not minorize the complete depolarizing channel')
    if np.any(counts < 1):
        raise ArithmeticError('some input/output population has no E-measured path')
    # Each measured history has rank-one Kraus operator. Subtract one such
    # |j><i|/16 from each of144 pairs; the remaining channel is still CP.
    return {'domain': 'repair-only invariant three-position highest-spin corner with two classical controller bits',
            'quantum_algebra_dimension': 36, 'Hilbert_space_population_dimension': 12,
            'ticks': 8, 'path_denominator': 256,
            'rank_one_path_counts': counts.tolist(),
            'minimum_rank_one_path_count': 1,
            'minorization': 'Phi^8-(3/64)D is completely positive, D(X)=Tr(X) I12/12 on the CQ algebra',
            'depolarizing_weight': F(3, 64), 'remaining_CPTP_weight': F(61, 64),
            'unique_faithful_stationary_state': 'I12/12',
            'nontrivial_spectrum_radius_upper': '(61/64)^(1/8)',
            'arbitrary_reference_trace_error': '||Phi^n(rho)-I12/12 tensor rho_R||1 <= (61/64)^floor(n/8) ||rho-I12/12 tensor rho_R||1',
            'all_off_diagonal_field_inputs_included': True,
            'full_Fock_or_interleaved_H_mixing_claimed': False}


def negative_controls():
    one = np.zeros((6, 6), dtype=np.int64)
    for first in (0, 1):
        for record, op in zip(c.LABELS, corner(first), strict=True):
            after = first ^ c.feedback(record)
            one[3*after:3*after+3, 3*first:3*first+3] += op*op
    output = one.sum(axis=1)
    if list(output) != [3, 1, 2, 1, 3, 2] or np.any(one.sum(axis=0) != 2):
        raise ArithmeticError('one-bit reversible-feedback bias changed')
    distance = sum(abs(F(int(x), 12)-F(1, 6)) for x in output)/2
    if distance != F(1, 6):
        raise ArithmeticError('one-bit bias distance changed')
    return {'one_bit_rule': 'c -> c XOR f(r), same selector c and same complete raw outcomes',
            'one_bit_map_reversible_for_each_fixed_record': True,
            'complete_one_bit_output_row_numerators_over2': output.tolist(),
            'uniform_input_trace_distance_after_one_tick': distance,
            'mere_recordwise_reversibility_does_not_imply_joint_unitality': True,
            'all_raw_records_or_outside_removed_is_not_the_declared_channel': True,
            'fixed_zero_working_memory_is_not_uniform_stationary_input': True}
