"""Exact recovery after classical commitment of complete native C/E records.

No quantum register is protected between repair attempts. A decoder may use
the final load register and the entire classical tape. It reverses the known
translations, recovering every initial population and the coherences within
each history fiber. It cannot undo the committed history measurement.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

from . import native
from ..common_source import slow_geometry


STATES = native.states()
INDEX = {x: i for i, x in enumerate(STATES)}


def candidate(value):
    if type(value) is not str or value not in ("C", "E"):
        raise ValueError("candidate C or E required")
    return value


def carrier(binary=False):
    if type(binary) is not bool:
        raise ValueError("binary must be a boolean")
    return native.states(binary=binary)


def seams():
    return tuple(slow_geometry.source().seams())


def seam(value):
    if (not isinstance(value, (list, tuple)) or len(value) != 2
            or any(type(v) is not int for v in value)
            or tuple(sorted(value)) not in seams()):
        raise ValueError("a native oriented seam required")
    return tuple(value)


def attempts(value):
    if type(value) is not int or value < 0:
        raise ValueError("a nonnegative integer attempt count required")
    return value


def load(value):
    if (not isinstance(value, (tuple, list)) or len(value) != 12
            or any(type(v) is not int or v < 0 for v in value)
            or sum(value) != 2):
        raise ValueError("a nonnegative integer total-two load required")
    return tuple(value)


def recorded_step(law, state, edge, coin):
    law, state, edge = candidate(law), load(state), seam(edge)
    if type(coin) is not int or coin not in (0, 1):
        raise ValueError("a binary integer coin required")
    output, emitted = native.native(state, edge, coin)
    output = load(output)
    total = state[edge[0]] + state[edge[1]]
    expected = list(state)
    expected[edge[0]], expected[edge[1]] = ((total + coin) // 2,
                                          (total + 1 - coin) // 2)
    current = state[edge[0]] - output[edge[0]]
    if (not isinstance(emitted, (tuple, list)) or len(emitted) != 2
            or any(type(v) is not int for v in emitted)
            or tuple(emitted) != (coin, current) or output != tuple(expected)):
        raise ArithmeticError("native endpoint output or emitted record changed")
    value = current if law == "C" else state[edge[0]] - state[edge[1]]
    return output, (*edge, coin, value)


def record_current(law, record):
    law = candidate(law)
    if (not isinstance(record, (list, tuple)) or len(record) != 4
            or any(type(v) is not int for v in record)):
        raise ValueError("a complete integer seam/coin/value record required")
    a, b, coin, value = record
    seam((a, b))
    if coin not in (0, 1) or value not in (range(-1, 2) if law == "C" else range(-2, 3)):
        raise ValueError("invalid native record value")
    return value if law == "C" else (value - (2 * coin - 1) * (abs(value) % 2)) // 2


def restore_load(law, terminal, history):
    """Invert an admitted complete history, rejecting forged local records.

    This is a basis-state decoder, not a claim that a hidden load is a public
    readout or that phase coherence between committed outcomes is recoverable.
    """
    law, state = candidate(law), load(terminal)
    if not isinstance(history, (tuple, list)):
        raise ValueError("a complete history sequence required")
    for record in reversed(history):
        current = record_current(law, record)
        a, b, coin, _ = record
        previous = list(state)
        previous[a] += current
        previous[b] -= current
        previous = load(previous)
        output, actual = recorded_step(law, previous, (a, b), coin)
        if output != state or actual != tuple(record):
            raise ValueError("record is incompatible with its terminal load")
        state = previous
    return state


@dataclass(frozen=True)
class Branch:
    history: tuple[tuple[int, int, int, int], ...]
    inputs: tuple[int, ...]
    outputs: tuple[int, ...]
    weight: F


def observer_instrument(law, *, binary=False):
    """The existing two-attempt current-controlled observer, fully committed.

    E uses its decoded CURRENT in the same selector, not the contrast sign.
    Indices in the result use the declared input carrier, which is invariant.
    """
    law, states = candidate(law), carrier(binary)
    index = {x: i for i, x in enumerate(states)}
    groups = defaultdict(list)
    for i, state in enumerate(states):
        for coins in product((0, 1), repeat=2):
            first, r1 = recorded_step(law, state, native.E, coins[0])
            edge = native.F if record_current(law, r1) > 0 else native.E
            final, r2 = recorded_step(law, first, edge, coins[1])
            history = (r1, r2)
            if restore_load(law, final, history) != state:
                raise ArithmeticError("native history does not invert")
            groups[history].append((i, index[final]))
    result = []
    for history, pairs in sorted(groups.items()):
        inputs, outputs = zip(*pairs, strict=True)
        if len(set(inputs)) != len(inputs) or len(set(outputs)) != len(outputs):
            raise ArithmeticError("record branch is not a partial permutation")
        result.append(Branch(history, inputs, outputs, F(1, 4)))
    if any(sum((b.weight for b in result if i in b.inputs), F(0)) != 1
           for i in range(len(states))):
        raise ArithmeticError("incomplete retained instrument")
    return tuple(result)


def observer_kernel(law, *, binary=False):
    states = carrier(binary)
    matrix = [[F(0)] * len(states) for _ in states]
    for branch in observer_instrument(law, binary=binary):
        for i in branch.inputs:
            for j in branch.inputs:
                matrix[i][j] += branch.weight
    return tuple(map(tuple, matrix))


def observer_report(law, *, binary=False):
    branches = observer_instrument(law, binary=binary)
    kernel = observer_kernel(law, binary=binary)
    dimension = len(kernel)
    blocks = {tuple(j for j in range(dimension) if kernel[i][j] == 1)
              for i in range(dimension)}
    if sum(len(block) for block in blocks) != dimension:
        raise ArithmeticError("full-input correctable algebra is not a partition")
    return {
        "dimension": dimension, "histories": len(branches),
        "branch_rank_counts": {str(k): v for k, v in sorted(Counter(len(b.inputs) for b in branches).items())},
        "full_input_correctable_block_sizes": sorted(map(len, blocks), reverse=True),
        "full_input_correctable_algebra_dimension": sum(len(b)**2 for b in blocks),
        "phase_survival_counts": {str(k): v for k, v in sorted(Counter(x for row in kernel for x in row).items())},
        "optimal_maximally_mixed_entanglement_fidelity": str(
            sum((b.weight * len(b.inputs)**2 for b in branches), F(0)) / dimension**2),
        "all_initial_populations_recovered": True,
    }


@lru_cache(maxsize=None)
def _table(law, edge, binary):
    states = carrier(binary)
    index = {x: i for i, x in enumerate(states)}
    return tuple(tuple((index[out], record[-1]) for state in states
                       for out, record in (recorded_step(law, state, edge, coin),))
                 for coin in (0, 1))


def fixed_word_kernel(law, word, *, binary=False):
    """All pair-history collision probabilities without materializing tapes."""
    law, states = candidate(law), carrier(binary)
    if not isinstance(word, (tuple, list)):
        raise ValueError("a finite seam word required")
    word = tuple(seam(edge) for edge in word)
    matrix = [[F(1)] * len(states) for _ in states]
    for edge in reversed(word):
        table = _table(law, edge, binary)
        following = [[F(0)] * len(states) for _ in states]
        for coin in (0, 1):
            for i, (x, rx) in enumerate(table[coin]):
                for j, (y, ry) in enumerate(table[coin]):
                    if rx == ry:
                        following[i][j] += matrix[x][y] / 2
        matrix = following
    return tuple(map(tuple, matrix))


def iid_rates(law, *, binary=False):
    """Exact same-record probabilities under uniform independent seam coins.

    Matching records translate both states identically, so their initial
    difference stays fixed. The t-attempt collision probability is rate**t.
    """
    law, states = candidate(law), carrier(binary)
    edges = seams()
    rates = []
    for left in states:
        row = []
        for right in states:
            delta = tuple(x-y for x, y in zip(left, right, strict=True))
            loss = sum(min(abs(delta[a]-delta[b]), 2) if law == "C"
                       else 2 * (delta[a] != delta[b]) for a, b in edges)
            row.append(1 - F(loss, 2 * len(edges)))
        rates.append(tuple(row))
    return tuple(rates)


def iid_report(law, *, binary=False):
    rates = iid_rates(law, binary=binary)
    dimension = len(rates)
    counts = Counter(x for row in rates for x in row)
    if counts[F(1)] != dimension or any(not 0 <= q <= 1 for q in counts):
        raise ArithmeticError("iid records do not distinguish all input pairs")
    maximum = max(q for q in counts if q < 1)
    return {
        "dimension": dimension,
        "rate_multiplicities": {str(q): n for q, n in sorted(counts.items())},
        "largest_offdiagonal_rate": str(maximum),
        "optimal_entanglement_fidelity": {
            str(t): str(sum((n*q**t for q, n in counts.items()), F(0)) / dimension**2)
            for t in (0, 1, 2, 4, 8, 16)},
        "limiting_optimal_entanglement_fidelity": str(F(1, dimension)),
        "reference_trace_distance_bound": f"min(1,({dimension-1}/2)*({maximum})^t)",
        "finite_time_entanglement_breaking": False,
    }


def one_event_code(law):
    """Exact KL Gram test for a logical qubit, not an idle protected factor.

    The two codewords have disjoint two-configuration support. Every nonzero
    recorded branch maps them to orthogonal vectors of equal norm. A final
    branch-dependent isometry therefore restores an arbitrary logical qubit.
    Its compilation is not claimed by this algebraic certificate.
    """
    law = candidate(law)
    pairs = (((0, 2), (1, 3)), ((0, 3), (1, 2)))
    supports = tuple(tuple(tuple(int(v in pair) for v in range(12))
                           for pair in logical) for logical in pairs)
    if any(tuple(sorted(pair)) in seams() for logical in pairs for pair in logical):
        raise ArithmeticError("the logical nonedge trade is not in this carrier")
    total = F(0)
    nonzero = 0
    probability_counts = Counter()
    for edge in seams():
        for coin in (0, 1):
            groups = defaultdict(lambda: [[], []])
            for logical, states in enumerate(supports):
                for state in states:
                    output, record = recorded_step(law, state, edge, coin)
                    groups[record][logical].append(output)
            for images in groups.values():
                left, right = images
                if (len(left) != len(right) or set(left).intersection(right)
                        or len(set(left)) != len(left) or len(set(right)) != len(right)):
                    raise ArithmeticError("native branch violates logical KL Gram identity")
                probability = F(len(left), 4 * len(seams()))
                probability_counts[probability] += 1
                total += probability
                nonzero += 1
    if total != 1:
        raise ArithmeticError("logical recovery certificate omitted outcomes")
    return {
        "logical_basis_pairs": [[list(pair) for pair in logical] for logical in pairs],
        "amplitudes": "1/sqrt(2) on each of the two disjoint supports",
        "logical_dimension": 2, "nonzero_record_branches": nonzero,
        "branch_probability_counts": {str(p): n for p, n in sorted(probability_counts.items())},
        "exact_correctability_after_one_event": True,
        "recovery": "branch-dependent isometry on the same load carrier",
        "compiled_logical_decoder": False,
    }
