"""Faithful native references and exact retained-record acceptance filters.

The sixty mixing records stay in the actual archive. The stated reference
observation keeps the tree and flag and marginalizes those sixty values.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import prod

from model.native_fair_local_preparation import preparation as p

MIXING = p.SEAMS * 2
DENOMINATOR = 1 << 60


def tree_weights():
    """A declared observer protocol: one plus the actual native tree depth."""
    depth = {0: 0}
    for a, b in p.TREE:
        if a in depth and b not in depth:
            depth[b] = depth[a] + 1
        elif b in depth and a not in depth:
            depth[a] = depth[b] + 1
        else:
            raise ArithmeticError('the source tree lost its rooted traversal')
    return tuple(1 + depth[i] for i in range(p.PORTS))


def occupation_weights(total):
    p.integer(total, 1, 11)
    weights = tree_weights()
    return tuple(prod(weights[i] for i in range(p.PORTS) if mask >> i & 1) for mask in p.masks(total))


def direct_reference(terminal):
    p.integer(terminal, 1, p.FULL - 1)
    total = terminal.bit_count()
    states = p.masks(total)
    column = states.index(terminal)
    return tuple(F(int(x), DENOMINATOR) for x in p.block_numerators(total)[:, column])


def filtered_reference(total):
    weights = occupation_weights(total)
    return tuple(F(w, sum(weights)) for w in weights)


@lru_cache(maxsize=11, typed=True)
def _filter_budget(total):
    weights = occupation_weights(total)
    matrix = p.block_numerators(total)
    multiplicity = min(int(matrix[i, :].min()) // weight for i, weight in enumerate(weights))
    acceptance = F(multiplicity * sum(weights), DENOMINATOR)
    if multiplicity <= 0 or not F(1, 12) <= acceptance <= 1:
        raise ArithmeticError('the all-input exact acceptance lower bound failed')
    return {'total': total, 'dimension': len(weights), 'path_multiplicity': multiplicity,
            'partition_weight': sum(weights), 'acceptance': acceptance,
            'rejection': 1 - acceptance, 'source_eigenvalue_floor': F(min(weights), sum(weights))}


def filter_budget(total):
    return dict(_filter_budget(total))


@lru_cache(maxsize=32, typed=True)
def suffix_counts(total, final):
    """Exact counts of suffix coin words to one final occupation.

    Sixty-one rows of at most924 integers, each in[0,2^60], suffice.
    The computation is finite classical record postprocessing, not a
    native repair opcode or a simulated sixty-coin enumeration.
    """
    p.integer(total, 1, 11)
    states = p.masks(total)
    if type(final) is not int or final not in states:
        raise ValueError('a final occupation in the stated sector is required')
    lookup = {state: i for i, state in enumerate(states)}
    rows = [None] * (len(MIXING) + 1)
    rows[-1] = tuple(int(state == final) for state in states)
    for t in range(len(MIXING)-1, -1, -1):
        edge = MIXING[t]
        following = rows[t+1]
        rows[t] = tuple(sum(following[lookup[p.readout.native.image(state, edge, coin, p.PORTS)[0]]]
                            for coin in (0, 1)) for state in states)
    return tuple(rows)


def _coins(coins):
    if (type(coins) is not tuple or len(coins) != len(MIXING)
            or any(type(coin) is not int or coin not in (0, 1) for coin in coins)):
        raise ValueError('sixty literal fair-coin bits are required')


def path_rank(initial, coins):
    """Lexicographic rank among all coin words with this exact final state."""
    p.integer(initial, 1, p.FULL - 1)
    _coins(coins)
    state = initial
    for edge, coin in zip(MIXING, coins, strict=True):
        state = p.readout.native.image(state, edge, coin, p.PORTS)[0]
    final = state
    states = p.masks(initial.bit_count())
    lookup = {state: i for i, state in enumerate(states)}
    counts = suffix_counts(initial.bit_count(), final)
    rank = 0
    state = initial
    for t, (edge, coin) in enumerate(zip(MIXING, coins, strict=True)):
        if coin:
            alternate = p.readout.native.image(state, edge, 0, p.PORTS)[0]
            rank += counts[t+1][lookup[alternate]]
        state = p.readout.native.image(state, edge, coin, p.PORTS)[0]
    count = counts[0][lookup[initial]]
    if state != final or not 0 <= rank < count:
        raise ArithmeticError('native path ranking lost its endpoint or multiplicity')
    return final, rank, count


def unrank(initial, final, rank):
    """Inverse path construction, used to falsify endpoint and threshold errors."""
    p.integer(initial, 1, p.FULL - 1)
    p.integer(rank)
    total = initial.bit_count()
    states = p.masks(total)
    lookup = {state: i for i, state in enumerate(states)}
    counts = suffix_counts(total, final)
    if rank >= counts[0][lookup[initial]]:
        raise ValueError('rank exceeds the exact native endpoint path count')
    state = initial
    coins = []
    for t, edge in enumerate(MIXING):
        zero = p.readout.native.image(state, edge, 0, p.PORTS)[0]
        count_zero = counts[t+1][lookup[zero]]
        if rank < count_zero:
            coin = 0
            state = zero
        else:
            coin = 1
            rank -= count_zero
            state = p.readout.native.image(state, edge, 1, p.PORTS)[0]
        coins.append(coin)
    if state != final or rank:
        raise ArithmeticError('native path unranking failed')
    return tuple(coins)


def classify(tree_slots, mixing_slots):
    """A complete classical flag from original native records alone.

    Constants retain their two-dimensional coherent code. No local total
    is inferred there. Invalid formal words stay explicit zero outcomes.
    """
    p.readout._slots(mixing_slots, len(MIXING))
    weight, pairs, tree_records = p.readout.branch(p.FAMILY, tree_slots)
    if not pairs:
        return {'status': 'ZERO_RAW_WORD', 'accepted': False}
    rows = [p.readout.alphabet(p.FAMILY, edge)[slot] for edge, slot in zip(MIXING, mixing_slots, strict=True)]
    if any(not factor for _, factor in rows):
        return {'status': 'ZERO_RAW_WORD', 'accepted': False}
    records = tuple(record for record, _ in rows)
    survivors = []
    for initial, terminal in pairs:
        state = terminal
        for edge, record in zip(MIXING, records, strict=True):
            following, current, contrast = p.readout.native.image(state, edge, record.coin, p.PORTS)
            if (current, contrast) != (record.current, record.label):
                break
            state = following
        else:
            survivors.append((initial, terminal, state))
    if not survivors:
        return {'status': 'ZERO_RAW_WORD', 'accepted': False}
    if len(pairs) == 2:
        if len(survivors) != 2:
            raise ArithmeticError('a fair E record resolved the coherent constant code')
        return {'status': 'CONSTANT_CODE', 'accepted': False,
                'constant_code_preserved': True, 'local_total_measured': False}
    initial, terminal, final = survivors[0]
    total = terminal.bit_count()
    computed_final, rank, count = path_rank(terminal, tuple(record.coin for record in records))
    if computed_final != final:
        raise ArithmeticError('path rank and authenticated native records disagree')
    budget = filter_budget(total)
    states = p.masks(total)
    threshold = budget['path_multiplicity'] * occupation_weights(total)[states.index(final)]
    return {'status': 'ACCEPT' if rank < threshold else 'REJECT', 'accepted': rank < threshold,
            'initial_mask': initial, 'tree_terminal_mask': terminal, 'terminal_mask': final,
            'local_total': total, 'path_rank': rank, 'endpoint_paths': count,
            'acceptance_threshold': threshold, 'fine_word_weight': weight * F(1, DENOMINATOR)}


def marginal(total, probabilities, keep):
    """Actual occupation partial trace, retaining every regional zero outcome."""
    p.integer(total, 1, 11)
    if (type(keep) is not tuple or tuple(sorted(set(keep))) != keep
            or any(type(bit) is not int or not 0 <= bit < p.PORTS for bit in keep)):
        raise ValueError('sorted distinct local port indices required')
    if (type(probabilities) is not tuple or len(probabilities) != len(p.masks(total))
            or any(type(x) is not F or x < 0 for x in probabilities) or sum(probabilities) != 1):
        raise ValueError('a complete normalized exact sector probability vector required')
    result = [F()] * (1 << len(keep))
    for mask, value in zip(p.masks(total), probabilities, strict=True):
        result[sum(((mask >> bit) & 1) << i for i, bit in enumerate(keep))] += value
    return tuple(result)


def modular_ratio(total, left, right, outside):
    """Compare a same-local-number coherence in one exterior context."""
    states = p.masks(total)
    if (type(left) is not int or type(right) is not int or type(outside) is not int
            or min(left, right, outside) < 0 or outside & (left | right)
            or left.bit_count() != right.bit_count()
            or (left | outside) not in states or (right | outside) not in states):
        raise ValueError('two equal local occupations and one disjoint feasible exterior required')
    weights = occupation_weights(total)
    return F(weights[states.index(left | outside)], weights[states.index(right | outside)])
