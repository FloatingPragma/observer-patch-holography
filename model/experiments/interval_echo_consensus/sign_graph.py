"""Unit-rate local permutations with shared or privately addressable sign memories."""
import json
from collections import Counter
from itertools import combinations, product

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.executable_model import Carrier, cell, federation
from .record_action import addresses, changes, signed_table
from .record_release import source


FREQUENCY_LEAKAGE_LIMIT = 1e-23


def common_width(board):
    """Largest number of distinct bystander sign addresses injectable at every contact."""
    return min(len(addresses(board, e))-1 for e in board.seams)


def questions(board, aperture="assembled", reverse=False):
    """A physical witness address for every comparison and single-load reference unit."""
    width = common_width(board)
    slots = {e: i for i, e in enumerate(board.seams)}
    offers, onsite = [], np.zeros(board.n, dtype=int)
    for e in board.seams:
        nearby = [f for f in addresses(board, e) if f != e]
        if aperture == "primitive":
            # Any injection works in the all-odd sector, including at an intercell contact.
            nearby = sorted(nearby, reverse=reverse)[:width]
            if len(nearby) != width:
                raise ValueError("Primitive contact chart cannot fit this carrier")
        elif aperture != "assembled":
            raise ValueError("Unknown contact aperture")
        offers.extend((e[0], e[1], slots[f]) for f in nearby)
        onsite[list(e)] += len(nearby)
    for v, count in enumerate(onsite.max()-onsite):
        local = sorted((slots[tuple(sorted((u, v)))] for u in board.neighbours[v]), reverse=reverse)
        offers.extend((v, None, local[j % len(local)]) for j in range(count))
    return offers, int(onsite.max())


class PrivateRecordWalk:
    """Lazy microscopic simple graph. An integer stores every privately addressed sign bit."""

    def __init__(self, board, aperture="primitive", reverse=False):
        self.board = board
        self.offers, self.rest = questions(board, aperture, reverse)
        self.port_masks = [0]*board.n
        for j, (a, b, _) in enumerate(self.offers):
            self.port_masks[a] |= 1 << j
            if b is not None:
                self.port_masks[b] |= 1 << j

    def active(self, pair):
        a, b = pair
        if not 0 <= a < b < self.board.n:
            raise ValueError("Use two distinct ordered port indices")
        return self.port_masks[a] ^ self.port_masks[b]

    def neighbors(self, pair, bank):
        mask = self.active(pair)
        if not isinstance(bank, int) or not 0 <= bank < 1 << len(self.offers):
            raise ValueError("Sign bank outside the finite private register")
        a, b = pair
        for moving, fixed in ((a, b), (b, a)):
            for dest in self.board.neighbours[moving]:
                if dest != fixed:
                    yield tuple(sorted((dest, fixed))), bank
        while mask:
            flip = mask & -mask
            yield pair, bank ^ flip
            mask ^= flip


def frequency_read(walsh, moves):
    """Eight controlled exposures and a Fourier read of the clock register."""
    size = len(walsh)
    if len(set(moves)) != len(moves) or any(not 0 < s < size for s in moves):
        raise ValueError("This phase reader requires distinct nonidentity unit flips")
    times = np.arange(size)*2*np.pi/size
    unit = np.exp(-1j*times)[:, None]
    # For each physical flip expose its even/odd subspaces, coherently
    # controlled by the digital clock. No spectral value is supplied to the reader.
    flight = np.ones((size, size), complex)
    for move in moves:
        parity = walsh[move]
        flight *= (1+parity)/2 + unit*(1-parity)/2
    probabilities = abs(np.fft.ifft(flight, axis=0))**2
    np.testing.assert_allclose(probabilities.sum(axis=0), 1, rtol=0, atol=2e-13)
    records = probabilities.argmax(axis=0)
    residual = probabilities.copy()
    residual[records, np.arange(size)] = 0.
    leakage = float(residual.sum(axis=0).max())
    if not leakage < FREQUENCY_LEAKAGE_LIMIT or np.any(records > size//2):
        raise ArithmeticError("A Fourier read is ambiguous on the finite unit-rate graph")
    return records, leakage


def witness_masks():
    """Exact character responses; single and joint masks challenge duplicated or coupled reads."""
    inventory = {}
    for board in (cell(), *(federation(i) for i in range(3))):
        counts = Counter(len(addresses(board, e))-1 for e in board.seams)
        width = common_width(board)
        if width != min(counts) or width != common_width(cell()):
            raise ArithmeticError("The portable address count differs from the primitive contact")
        inventory[board.cells] = dict(sorted(counts.items()))
    width = common_width(cell())
    # Distinct sign markers commute. The complete collection of independently set
    # characters separates all diagonal interactions, including nonlinear responses.
    board = cell()
    edge = board.seams[0]
    readers = [f for f in addresses(board, edge) if f != edge]
    state = np.zeros(board.n, dtype=int)
    state[edge[0]] = 1
    delta = changes(readers, edge, state, signed_table())
    weights = np.array(delta) != 0
    masks = np.array(list(product((0, 1), repeat=width)), dtype=int)
    response = masks @ weights
    if not np.array_equal(response, masks.sum(axis=1)):
        raise ArithmeticError("Masked word tests disagree with the registered addresses")
    # Literal sign permutations, before selecting any Fourier character of the bank.
    bank = np.arange(2**width)
    walsh = 1-2*((masks @ masks.T) % 2)
    twice_action = np.zeros_like(walsh)
    for j, changed in enumerate(delta):
        if changed:
            twice_action += walsh - walsh[bank ^ (1 << (width-1-j))]
    if not np.array_equal(twice_action, 2*walsh*response):
        raise ArithmeticError("Microscopic sign flips give the wrong joint phase response")
    # Recover every rate for a simultaneous flip of any nonempty subset from
    # the *executed* single-load character energies. No sparse-rate ansatz.
    energies = twice_action[0]//2
    coefficients = -2*(walsh @ energies)
    if np.any(coefficients % len(bank)):
        raise ArithmeticError("Nonintegral inverse character response")
    rates = (coefficients//len(bank))[1:]
    expected = np.array([int(s.bit_count() == 1) for s in range(1,len(bank))])
    if not np.array_equal(rates, expected):
        raise ArithmeticError("A correlated clock flip or a missing local rate escaped")
    # Test all four binary occupations under a single held contact. The two
    # equal-load inputs would expose a pair-only backreaction invisible at one load.
    occupations, wrong_outcome = {}, 0.
    for a, b in product((0, 1), repeat=2):
        state[list(edge)] = (a, b)
        word_changes = changes(readers, edge, state, signed_table())
        character_action = np.zeros_like(walsh)
        for j, changed in enumerate(word_changes):
            if changed:
                character_action += walsh - walsh[bank ^ (1 << (width-1-j))]
        costs = character_action[0]//2
        if not np.array_equal(costs, response*int(a != b)):
            raise ArithmeticError("Binary contact action depends on equal loads or load orientation")
        actual_moves = [1 << (width-1-j) for j, changed in enumerate(word_changes) if changed]
        records, leakage = frequency_read(walsh, actual_moves)
        wrong_outcome = max(wrong_outcome, leakage)
        if not np.array_equal(records, costs):
            raise ArithmeticError("The phase-register observations differ from the generator costs")
        occupations[f"{a}{b}"] = int(records[-1])
    # Competing simple graph: sixteen distinguishable flips enabled only on
    # two occupied endpoints. It agrees on the entire one-load experiment.
    joint_moves = [1 << j for j in range(width)]
    joint_moves += [sum(1 << ((j+k) % width) for k in range(3)) for j in range(width)]
    if len(set(joint_moves)) != len(joint_moves):
        raise ArithmeticError("Pair-only challenge reused a microscopic destination")
    extra = sum((1-walsh[s])//2 for s in joint_moves)
    competing_records, leakage = frequency_read(walsh, joint_moves)
    extreme, check_leakage = frequency_read(walsh, list(range(1, len(bank))))
    wrong_outcome = max(wrong_outcome, leakage, check_leakage)
    if (not np.array_equal(extra, competing_records) or not np.all(extreme[1:] == len(bank)//2)
            or extreme[0] or extra[-1] != 2*width or extra[0] != 0 or occupations["11"] != 0):
        raise ArithmeticError("The equal-load phase reference failed to distinguish the pair-only walk")
    duplicate = response + masks[:, 0]-masks[:, 1]
    nonlinear = response + masks[:, 0]*masks[:, 1]
    differing = [int(np.count_nonzero(other != response)) for other in (duplicate, nonlinear)]
    if not all(differing) or duplicate[-1] != response[-1] or np.any((nonlinear != response) &
                                                                      (masks.sum(axis=1) == 1)):
        raise ArithmeticError("The masking controls were not discriminating")
    return dict(width=width, inventories=inventory, sign_patterns=len(masks), rejected_masks=differing,
                recovered_single_rates=rates[expected == 1].tolist(),
                recovered_zero_joint_rates=int(np.count_nonzero(rates == 0)),
                equal_and_unequal_occupations=occupations, rival_equal_load_phase_rate=int(extra[-1]),
                frequency_register_bits=width, maximum_forbidden_read_probability=wrong_outcome)


def defects(pair, offers):
    return [(a in pair) != (b in pair) for a, b, _ in offers]


def microscopic_path(private=False):
    board = Carrier("open path", 3, [(0, 1), (1, 2)])
    pairs = list(combinations(range(board.n), 2))
    offers, _ = questions(board)
    if private:
        offers = [(a, b, j) for j, (a, b, _) in enumerate(offers)]
    signs = 2**(len(offers) if private else len(board.seams))
    dim = signs*len(pairs)
    h2 = np.zeros((dim, dim), dtype=int)      # Twice the Hamiltonian, in unit parity-gap convention.
    for j, pair in enumerate(pairs):
        for bank in range(signs):
            start = signs*j+bank
            for a, b in board.seams:
                move = {a: b, b: a}
                end = pairs.index(tuple(sorted(move.get(v, v) for v in pair)))
                h2[start, start] += 1
                h2[signs*end+bank, start] -= 1
            for active, (_, _, sign) in zip(defects(pair, offers), offers, strict=True):
                end = bank ^ (int(active) << sign)
                h2[start, start] += 1
                h2[signs*j+end, start] -= 1
    if not np.array_equal(h2, h2.T) or np.any(h2.sum(axis=0)) or np.any((h2-np.diag(h2.diagonal())) > 0):
        raise ArithmeticError("Invalid unweighted microscopic walk")
    max_multiplicity = int(-np.min(h2))
    if private and max_multiplicity != 1:
        raise ArithmeticError("Two distinct private-record transitions have the same endpoints")
    kinetic2 = np.array([[1, -1, 0], [-1, 2, -1], [0, -1, 1]])
    evolved = []
    for mask in range(signs):
        character = np.array([(-1)**(s & mask).bit_count() for s in range(signs)])
        V = np.kron(np.eye(len(pairs), dtype=int), character[:, None])
        penalties = [sum(bit*((mask >> sign) & 1) for bit, (_, _, sign) in zip(defects(pair, offers), offers,
                                                                               strict=True)) for pair in pairs]
        reduced2 = kinetic2+2*np.diag(penalties)
        if not np.array_equal(h2 @ V, V @ reduced2):
            raise ArithmeticError("A microscopic character block differs from its material walk")
        evolved.append(expm_multiply(-1j*reduced2/2, np.eye(len(pairs))[:, 0]))
    # The first sign is sharp and all others are odd. Retain every final sign outcome.
    prepared = np.zeros(dim)
    prepared[:signs:2] = [(-1)**s.bit_count()/np.sqrt(signs/2) for s in range(0, signs, 2)]
    final = expm_multiply(-1j*h2/2, prepared).reshape(len(pairs), signs)
    terminal = (abs(final)**2).sum(axis=1)
    np.testing.assert_allclose(terminal, (abs(evolved[-2])**2+abs(evolved[-1])**2)/2, atol=2e-14, rtol=0)
    if np.max(abs(abs(evolved[-2])**2-abs(evolved[-1])**2)) <= 1e-3:
        raise ArithmeticError("The wrong-character control has no resolved effect")
    return dict(states=dim, character_blocks=signs, max_multiplicity=max_multiplicity, sharp_marker_material=terminal.tolist())


def full_projection(arm, reverse=False):
    board = federation(0)
    walk = PrivateRecordWalk(board, arm, reverse)
    offers, rest = walk.offers, walk.rest
    positions = np.array(list(combinations(range(board.n), 2)))
    i, j = positions.T
    diagonal = np.zeros(len(i), dtype=int)
    for a, b, _ in offers:
        occupied_a = (i == a) | (j == a)
        occupied_b = False if b is None else (i == b) | (j == b)
        diagonal += occupied_a ^ occupied_b
    private_cost = np.fromiter((walk.active((int(a), int(b))).bit_count() for a, b in positions), int)
    if not np.array_equal(private_cost, diagonal):
        raise ArithmeticError("Private simple-graph degree and compact named-clock cost differ")
    # PairSector supplies the native exchange and its upper-triangle ordering.
    from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
    sector = PairSector(board)
    if not np.array_equal(i, sector.aa) or not np.array_equal(j, sector.bb):
        raise ArithmeticError("Independent pair order differs")
    projected = sector.exchange()+sparse.diags(diagonal-2*rest, dtype=np.int64, format="csr")
    reference, _, _, _ = source(arm)
    if (projected != reference).nnz:
        raise ArithmeticError("Microscopic graph disagrees on a material matrix element")
    return dict(private_signs=len(offers), shared_sign_memories=len(board.seams), load_exchanges=len(board.seams), sign_checks=len(offers),
                material_states=len(i), one_load_reference=rest)


if __name__ == "__main__":
    print(json.dumps(dict(path=microscopic_path(), private_path=microscopic_path(True), masking=witness_masks(), primitive=full_projection("primitive"),
                          other_chart=full_projection("primitive", True), assembled=full_projection("assembled")),
                     allow_nan=False))
