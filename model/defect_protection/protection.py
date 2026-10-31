"""Exact native ordered-cover sector; the added memory is part of the domain.

This is a structural composition of the existing ordered-pair lift and the
complete C/E instrument, not a new spectral fit.  The source has 66 states.
One persistent bit records the ordering of its two occupied ports.  Every
native move transports that ordering; its parity is derived from the move.
"""
from __future__ import annotations

from collections import deque
from fractions import Fraction as F
import hashlib
from itertools import combinations
import json

from ..exchange_sector import statistics as exchange


DIMENSION = 66
COVER_DIMENSION = 132


def branches(law):
    """All raw slots, with integer partial maps of amplitude 1/sqrt(2)."""
    return exchange.branches(law, -1)


def lifted_pairs(row, chart=None):
    """The ordered-cover map in a possibly recharted orientation basis.

    A chart is one bit at each unordered source pair.  Changing it conjugates
    the complete map by G= sum_x |x><x| tensor X**chart[x].
    """
    if chart is None:
        chart = (0,)*DIMENSION
    if (type(chart) is not tuple or len(chart) != DIMENSION
            or any(type(x) is not int or x not in (0, 1) for x in chart)):
        raise ValueError('one exact orientation bit per source pair required')
    return tuple((r*DIMENSION+i,
                  (r ^ int(sign == -1) ^ chart[i] ^ chart[j])*DIMENSION+j)
                 for i, j, sign in zip(row.inputs, row.outputs, row.signs, strict=True)
                 for r in (0, 1))


def exchange_loop(law):
    """An actual three-move closed native history with odd ordered lift."""
    pair = (0, 11)
    index = exchange.INDEX[pair]
    rows, parity = [], 0
    for seam, target in zip(((0, 5), (0, 11), (5, 11)),
                            ((5, 11), (0, 5), (0, 11)), strict=True):
        target_index = exchange.INDEX[target]
        matches = [row for row in branches(law) if row.seam == seam
                   and any(i == index and j == target_index
                           for i, j in zip(row.inputs, row.outputs, strict=True))]
        if len(matches) != 1:
            raise ArithmeticError('the closed exchange has no unique native moving flag')
        row = matches[0]
        parity ^= int(row.signs[row.inputs.index(index)] == -1)
        rows.append((row.law, row.seam, row.coin, row.raw))
        index = target_index
    if index != exchange.INDEX[pair] or parity != 1:
        raise ArithmeticError('native closed exchange did not flip the ordering')
    history = exchange.history_map(tuple(rows), -1)
    source_index = exchange.INDEX[pair]
    if (history.inputs != (source_index,) or history.outputs != (source_index,)
            or history.signs != (-1,)):
        raise ArithmeticError('complete closed-history source support changed')
    return {'initial_pair': list(pair), 'history': [list(row) for row in rows],
            'source_history_rank': 1, 'orientation_action': 'X',
            'probability_on_initial_pair': str(history.squared_amplitude),
            'character_on_X_plus': 1, 'character_on_X_minus': -1}


def algebra_certificate(law):
    """Construct M66 tensor span{I,X}, without numerical rank decisions.

    Effects separate all 66 source addresses, hence generate P_x tensor I.
    Tree edges generate transport matrix units.  One odd loop gives
    P_x tensor X; connected transport gives it at every vertex.  Multiplying
    yields every E_xy tensor I and E_xy tensor X.  All original generators
    commute with I tensor X, proving equality, not only a lower bound.
    """
    rows = branches(law)
    supports = [set(row.inputs) for row in rows]
    unresolved, selected = list(combinations(range(DIMENSION), 2)), []
    while unresolved:
        scores = [sum((a in support) != (b in support) for a, b in unresolved)
                  for support in supports]
        chosen = max(range(len(rows)), key=scores.__getitem__)
        if scores[chosen] == 0:
            raise ArithmeticError('native effects no longer separate source addresses')
        selected.append(chosen)
        unresolved = [(a, b) for a, b in unresolved
                      if (a in supports[chosen]) == (b in supports[chosen])]
    signatures = [sum(int(i in supports[k]) << n for n, k in enumerate(selected))
                  for i in range(DIMENSION)]
    if len(set(signatures)) != DIMENSION:
        raise ArithmeticError('source projector construction failed')

    adjacency = {i: {} for i in range(DIMENSION)}
    cover = {i: set() for i in range(COVER_DIMENSION)}
    incoming, outgoing = [0]*COVER_DIMENSION, [0]*COVER_DIMENSION
    for k, row in enumerate(rows):
        pairs = lifted_pairs(row)
        if len({a for a, _ in pairs}) != len(pairs) or len({b for _, b in pairs}) != len(pairs):
            raise ArithmeticError('lift is not a partial injection')
        for a, b in pairs:
            incoming[b] += 1
            outgoing[a] += 1
            cover[a].add(b)
        # Commutation with X is an exact equality of all columns, not only
        # equality on the two selected states of a witness experiment.
        lookup = dict(pairs)
        flip = lambda i: (i+DIMENSION) % COVER_DIMENSION
        if any(lookup.get(flip(a)) != flip(b) for a, b in pairs):
            raise ArithmeticError('a complete native branch fails memory-charge conservation')
        for a, b, sign in zip(row.inputs, row.outputs, row.signs, strict=True):
            if a != b:
                adjacency[a].setdefault(b, (k, int(sign == -1)))
                adjacency[b].setdefault(a, (k, int(sign == -1)))
    # Each source has two coin outputs per seam.  This is a separate uniform
    # all-seam channel; arbitrary scheduling needs only branch commutation.
    if incoming != [60]*COVER_DIMENSION or outgoing != [60]*COVER_DIMENSION:
        raise ArithmeticError('uniform lifted native channel is not bistochastic')
    if any(a not in cover[b] for a, neighbors in cover.items() for b in neighbors):
        raise ArithmeticError('ordered native movement graph lost reversal support')
    parents, gauge, queue = {}, {0: 0}, deque([0])
    while queue:
        parent = queue.popleft()
        for child, (k, parity) in sorted(adjacency[parent].items()):
            if child in gauge:
                continue
            gauge[child] = gauge[parent] ^ parity
            parents[child] = [parent, k, parity]
            queue.append(child)
    if len(gauge) != DIMENSION:
        raise ArithmeticError('native source transport graph disconnected')
    reached, queue = {0}, deque([0])
    while queue:
        a = queue.popleft()
        for b in cover[a]:
            if b not in reached:
                reached.add(b)
                queue.append(b)
    if len(reached) != COVER_DIMENSION:
        raise ArithmeticError('ordered cover ceased to be connected')
    loop = exchange_loop(law)
    return {'raw_slots_including_zero_and_abort': len(rows),
            'zero_slots': sum(not row.inputs for row in rows),
            'abort_slots': sum(row.abort for row in rows),
            'source_dimension': DIMENSION, 'joint_dimension': COVER_DIMENSION,
            'lifted_partial_injection_columns': sum(len(lifted_pairs(row)) for row in rows),
            'complete_branch_X_commutation': True,
            'uniform_channel_denominator': 60,
            'uniform_channel_bistochastic': True,
            'separating_effect_indices': selected, 'source_effect_signatures': signatures,
            'source_spanning_tree': [[child, *parents[child]] for child in sorted(parents)],
            'ordered_cover_reachable_states': len(reached),
            'classical_configuration_invariants': 'constants only on the connected 132-state cover',
            'odd_closed_native_history': loop,
            'generated_star_algebra': 'M66(C) tensor span{I,X} = M66(C) direct_sum M66(C)',
            'generated_star_algebra_complex_dimension': 2*DIMENSION**2,
            'commutant': 'span{I132, I66 tensor X}', 'commutant_complex_dimension': 2,
            'uniform_channel_fixed_observables': 'span{I132, I66 tensor X}',
            'raw_maps_sha256': hashlib.sha256(json.dumps([
                [row.law, row.seam, row.coin, row.raw, lifted_pairs(row)] for row in rows],
                separators=(',', ':')).encode()).hexdigest()}


def memory_and_localization_certificate():
    """Exact reduced-state and memory-error statements on this same domain."""
    # rho_+/- = I66/66 tensor |+/- ><+/-|.  Every branch sends its source
    # factor to the same diagonal partial-injection image, at any history.
    checked = 0
    for law in ('C', 'E'):
        for row in branches(law):
            plus = {j: F(1, 132) for j in row.outputs}
            minus = {j: F(sign*sign, 132)
                     for j, sign in zip(row.outputs, row.signs, strict=True)}
            if plus != minus:
                raise ArithmeticError('the two reference sectors have unequal visible native output')
            checked += 1
    return {'reference': 'I66/66 on the binary native source',
            'sharp_states': 'rho_plus/minus = I66/66 tensor (I +/- X)/2',
            'source_reduced_state_both_sectors': 'I66/66',
            'joint_state_trace_distance': '1', 'source_trace_distance': '0',
            'source_plus_complete_classical_native_tape_distance_any_finite_history': '0',
            'diagonal_source_branch_checks': checked,
            'memory_dephasing_Z_basis': 'both sharp states become I66/66 tensor I2/2',
            'memory_reset_to_zero': 'both sharp states become I66/66 tensor |0><0|',
            'random_Z_error_probability_p': 'X expectation after n cycles is (1-2p)^n times initial expectation',
            'exact_preservation_condition': 'each admitted joint memory/service Kraus operator commutes with I66 tensor X',
            'sector_coherence_protected': False,
            'source_localized_defect': False,
            'remaining_join': 'embed the orientation/charge algebra in a specified local core and prove the complete memory/repair law preserves it; an observer bit with reference-like source exterior is insufficient'}
