"""Independent rational return enclosures on the enrolled physical tree."""

from collections import deque
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, pairwise
import json
from math import comb, factorial
from pathlib import Path

from model.executable_model import federation


HERE = Path(__file__).resolve().parent
ROOT = 1045
BRANCHES = ((1010,), (1044, 1028), (1052, 1047, 1032))
EDGES = tuple(edge for branch in BRANCHES for edge in pairwise((ROOT, *branch)))
INITIAL_HEADS = (1045, 1044, 1028, 1045, 1052, 1047)
INITIAL_PAIR = (1028, 1045)
HOP = 4
EVIDENCE_ALLOWANCE = Fraction(1, 10**12)


def head_at_occupation(pair, sign):
    """Solve Gauss on each leaf-side subtree from the enrolled word."""
    word = 0
    edge_number = 0
    for branch in BRANCHES:
        for level, (parent, child) in enumerate(pairwise((ROOT, *branch))):
            subtree = branch[level:]
            before = sum(v in INITIAL_PAIR for v in subtree)
            after = sum(v in pair for v in subtree)
            head = INITIAL_HEADS[edge_number]
            assert head in (parent, child)
            inward = int(head == child)
            if sign == -1:
                inward = 1 - inward
            inward += sign * (after - before)
            if inward not in (0, 1):
                return None
            word |= inward << edge_number
            edge_number += 1
    return word


def occupation_graph(degrees, sign, owner_guard, cut):
    """First solve all tree sectors, then connect admissible occupation pairs."""
    words = {}
    for pair in combinations(sorted(degrees), 2):
        if owner_guard and pair[0] // 12 == pair[1] // 12:
            continue
        word = head_at_occupation(pair, sign)
        if word is not None:
            words[pair] = word
    physical_edges = EDGES[int(cut) :]
    physical_bits = {frozenset(edge): EDGES.index(edge) for edge in physical_edges}
    neighbors = {pair: {} for pair in words}
    for first, last in combinations(words, 2):
        moved = frozenset(first) ^ frozenset(last)
        if moved not in physical_bits:
            continue
        bit = physical_bits[moved]
        assert words[first] ^ words[last] == 1 << bit
        (old,) = (p for p in first if p not in last)
        (new,) = (p for p in last if p not in first)
        current_head = EDGES[bit][int(bool(words[first] & (1 << bit)))]
        assert (current_head == old) == (sign == 1)
        # Follow two individually labelled particles for this one move.
        # Reordering the outgoing labels to the ordered wedge basis gives
        # the sign. No simulator's fermionic generator enters this check.
        labelled = tuple(new if v == old else v for v in first)
        reversal = labelled != last
        neighbors[first][last] = reversal
        neighbors[last][first] = reversal

    initial = (INITIAL_PAIR, False)
    distance = {initial: 0}
    queue = deque([initial])
    while queue:
        pair, ordering = queue.popleft()
        for other, reversal in neighbors[pair].items():
            target = (other, ordering ^ reversal)
            if target not in distance:
                distance[target] = distance[pair, ordering] + 1
                queue.append(target)
    states = sorted({state for state, _ in distance})
    indices = {state: i for i, state in enumerate(states)}
    native = [[0] * len(states) for _ in states]
    polarized = [[0] * len(states) for _ in states]
    for i, state in enumerate(states):
        diagonal = HOP * len(neighbors[state]) + sum(
            8 * (6 - degrees[p]) for p in state
        )
        native[i][i] = polarized[i][i] = diagonal
        for other, reversal in neighbors[state].items():
            j = indices[other]
            native[i][j] = -HOP
            polarized[i][j] = HOP if reversal else -HOP
    lag = distance.get((INITIAL_PAIR, True))
    if lag is None:
        # A cut removes exchange. Every reachable pair has exactly one
        # label parity, which rephases the polarized matrix to the native one.
        parity = [1 if (state, False) in distance else -1 for state in states]
        assert all(
            native[i][j] == parity[i] * polarized[i][j] * parity[j]
            for i in range(len(states))
            for j in range(len(states))
        )
    return (
        [[list(pair), words[pair]] for pair in states],
        indices[INITIAL_PAIR],
        lag,
        native,
        polarized,
    )


def arctan_interval(x):
    terms = [(-1) ** k * x ** (2 * k + 1) / (2 * k + 1) for k in range(50)]
    lower = sum(terms)
    return lower, lower + x**101 / 101


def waiting_interval():
    # Machin's identity: tan(4 arctan(1/5) - arctan(1/239)) = 1.
    # Both small positive arcs place that difference in the first quadrant.
    def tan_add(a, b):
        return (a + b) / (1 - a * b)

    fifth, last = Fraction(1, 5), Fraction(1, 239)
    twice = tan_add(fifth, fifth)
    assert tan_add(tan_add(twice, twice), -last) == 1
    a, b = arctan_interval(fifth)
    c, d = arctan_interval(last)
    lo, hi = 6 * a - 3 * d / 2, 6 * b - 3 * c / 2
    scale = 10**24
    start = Fraction(lo * scale // 1, scale)
    assert start <= lo <= hi <= start + Fraction(1, scale)
    return start, start + Fraction(1, scale)


def moments(matrix, source, order):
    vector = [int(i == source) for i in range(len(matrix))]
    result = [1]
    for _ in range(order):
        vector = [
            sum(a * b for a, b in zip(row, vector, strict=True)) for row in matrix
        ]
        result.append(vector[source])
    return result


def doubled_moment(values, k, coefficients):
    """Diagonal moment of H⊗I−I⊗H; no doubled matrix is needed."""
    return sum(
        (-1) ** j * coefficient * values[j] * values[k - j]
        for j, coefficient in enumerate(coefficients)
    )


def return_enclosures(native, polarized, source, lag):
    """Use a real cosine polynomial with operator-norm Taylor remainder."""
    time, upper_time = waiting_interval()
    diagonal = [row[i] for i, row in enumerate(native)]
    midpoint = Fraction(min(diagonal) + max(diagonal), 2)
    assert midpoint.denominator == 1
    shifted = [
        [
            [a - int(i == j) * int(midpoint) for j, a in enumerate(row)]
            for i, row in enumerate(matrix)
        ]
        for matrix in (native, polarized)
    ]
    diameter = 2 * max(sum(abs(a) for a in row) for matrix in shifted for row in matrix)
    scaled = diameter * time
    order = 16
    limit = Fraction(1, 10**20)
    while scaled ** (order + 2) / factorial(order + 2) > limit:
        order += 2
    # The scalar cosine remainder gives the same bound on a selfadjoint
    # matrix, and the unknown fraction of the time unit costs at most
    # ||H⊗I−I⊗H|| times the time uncertainty in the expectation.
    error = scaled ** (order + 2) / factorial(order + 2) + diameter * (
        upper_time - time
    )
    bases = [moments(matrix, source, order) for matrix in shifted]
    integer_coefficients = [
        tuple(comb(k, j) for j in range(k + 1)) for k in range(0, order + 1, 2)
    ]
    doubled = [
        [
            doubled_moment(values, 2 * k, binomial)
            for k, binomial in enumerate(integer_coefficients)
        ]
        for values in bases
    ]
    difference = [b - a for a, b in zip(*doubled, strict=True)]
    first = next((2 * k for k, p in enumerate(difference) if p), None)
    assert first == lag
    coefficient = None
    if lag is not None:
        coefficient = Fraction(
            (-1) ** (lag // 2) * difference[lag // 2], factorial(lag)
        )
        assert bases[0][:lag] == bases[1][:lag]
    polynomials = [
        sum(
            Fraction((-1) ** k * value, factorial(2 * k)) * time ** (2 * k)
            for k, value in enumerate(sequence)
        )
        for sequence in doubled
    ]
    assert error < EVIDENCE_ALLOWANCE
    return [(value - error, value + error) for value in polynomials], coefficient


def decimals(interval):
    scale = 10**15
    lo, hi = interval
    endpoints = (int(lo * scale // 1), -int(-hi * scale // 1))
    return "[" + ", ".join(f"{x // scale}.{x % scale:015d}" for x in endpoints) + "]"


def main():
    data = json.loads((HERE / "records.json").read_text(), parse_float=Fraction)
    assert (
        sha256((HERE / "DESIGN.md").read_bytes()).hexdigest() == data["design_sha256"]
    )
    assert data["edges_and_starting_excess_head"] == [
        [*edge, head] for edge, head in zip(EDGES, INITIAL_HEADS, strict=True)
    ]
    ports = sorted({p for edge in EDGES for p in edge})
    board = federation(2)
    assert all(v in board.neighbours[u] for u, v in EDGES)
    degrees = {p: int(board.degree[p]) for p in ports}
    assert list(map(list, degrees.items())) == data["ports_degrees"]

    for case in data["cases"]:
        guard, cut = case["distinct_owner_guard"], case["spare_arm_cut"]
        states, source, lag, native, polarized = occupation_graph(
            degrees, 1, guard, cut
        )
        conjugate = occupation_graph(degrees, -1, guard, cut)
        assert conjugate[0] == [[pair, word ^ 63] for pair, word in states]
        assert conjugate[1:] == (source, lag, native, polarized)
        assert states == case["address_pair_and_excess_head_word"]
        assert source == case["physical_return_bin_index"]
        assert lag == case["shortest_odd_exchange_walk"]
        bounds, coefficient = return_enclosures(native, polarized, source, lag)
        assert (str(coefficient) if coefficient is not None else None) == case[
            "exact_first_polarized_minus_native_return_coefficient"
        ]
        (report,) = (
            row for row in case["instruments"] if row["wait_in_pi_over_16"] == 6
        )
        probability_banks = report["native_polarized_and_labelled"]
        for interval, arm in zip(bounds, probability_banks[:2], strict=True):
            saved = arm["address_probabilities_given_both_inputs"][source]
            assert max(abs(saved - end) for end in interval) < EVIDENCE_ALLOWANCE
            assert (
                abs(
                    data["first_arrival_per_original"] * saved
                    - arm["return_per_original"]
                )
                < EVIDENCE_ALLOWANCE
            )
        delta = (bounds[1][0] - bounds[0][1], bounds[1][1] - bounds[0][0])
        assert (
            max(
                abs(report["polarized_minus_native_return_given_both"] - end)
                for end in delta
            )
            < 2 * EVIDENCE_ALLOWANCE
        )
        if lag is None:
            assert bounds[1] == bounds[0]
        else:
            assert delta[1] < 0
        print(
            f"guard={guard}, cut={cut}, states={len(states)}, exchange order={lag}: "
            f"native {decimals(bounds[0])}, polarized {decimals(bounds[1])}"
        )


if __name__ == "__main__":
    main()
