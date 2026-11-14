"""Phase-exact source words for the complete 78-state C/E instruments.

The seam occupies bits 0..2 and a six-bit spectator label occupies bits
3..8. Spectator labels can be reused between different occupation shells:
the seam total determines which shell they belong to. Proper-code load and
export, reset, computational pointer readout and classical control remain
the apparatus services declared by the existing source compiler.
"""
from __future__ import annotations

import math

from model.observer import native
from model.observer.source_scalar import compiler as source


SOURCE_QUBITS = 9
POINTER = 9
SEAM_STATES = ((0, 0), (0, 1), (1, 0), (0, 2), (1, 1), (2, 0))
CURRENT_CODES = {0: 0, 1: 1, -1: 2}
RECORD_BITS = {"C": 2, "E": 3}
SETTINGS = ("phase_0", "phase_pi", "phase_pi_over_2", "diagonal_control")

# All endpoints are distinct, and the destination of a moving admitted
# column has a nonblank record. Thus these disjoint swaps complete the
# current-record partial injection to a unitary on the entire local register.
# Contrast recording uses this same repair followed by one parity copy.
REPAIR_SWAPS = {
    ("C", 0): ((2, 9), (3, 20), (5, 12)),
    ("C", 1): ((1, 18), (3, 20), (5, 12)),
}


def record_values(candidate, coin):
    """Decode admitted raw records; absent codes are retained abort outcomes.

    C stores the current code. E stores current_code + 4*parity, with
    contrast d = 2*current + (2*coin-1)*parity. The reserved current code
    and the impossible contrasts +/-3 are excluded from the native carrier.
    """
    if candidate not in RECORD_BITS or type(coin) is not int or coin not in (0, 1):
        raise ValueError("candidate C/E and a binary integer coin required")
    if candidate == "C":
        return {0: 0, 1: 1, 2: -1}
    if coin == 0:
        return {0: 0, 1: 2, 2: -2, 4: -1, 5: 1}
    return {0: 0, 1: 2, 2: -2, 4: 1, 6: -1}


def _spectator_shells():
    shells = []
    preferred = {1: (1,) + (0,) * 9, 2: (1, 1) + (0,) * 8}
    for total in range(3):
        entries = sorted({state[2:] for state in native.states()
                          if sum(state[2:]) == total})
        if total in preferred:
            entries.remove(preferred[total])
            entries.insert(0, preferred[total])
        shells.append(tuple(entries))
    return tuple(shells)


SPECTATOR_SHELLS = _spectator_shells()


def encode(state):
    """Embed a native total-two state; the prediction's X and Z map to 0,1."""
    if (not isinstance(state, (list, tuple)) or len(state) != 12
            or any(type(value) is not int or value < 0 for value in state)
            or sum(state) != 2):
        raise ValueError("a twelve-port nonnegative integer total-two state required")
    seam = SEAM_STATES.index(tuple(state[:2]))
    spectators = tuple(state[2:])
    label = SPECTATOR_SHELLS[sum(spectators)].index(spectators)
    return seam + 8 * label


def decode(code):
    """Reject padding; labels are interpreted in the seam's spectator shell."""
    if type(code) is not int or not 0 <= code < 1 << SOURCE_QUBITS:
        raise ValueError("a native code integer in the nine-qubit register required")
    seam, label = code % 8, code // 8
    if seam >= len(SEAM_STATES):
        raise ValueError("source code is outside the 78-state carrier")
    endpoints = SEAM_STATES[seam]
    shell = SPECTATOR_SHELLS[2 - sum(endpoints)]
    if label >= len(shell):
        raise ValueError("source code is outside the 78-state carrier")
    return endpoints + shell[label]


SOURCE_CODES = tuple(encode(state) for state in native.states())


def _integer(value, low, high, name):
    if type(value) is not int or not low <= value < high:
        raise ValueError(f"{name} must be an integer in [{low}, {high})")
    return value


def controlled_flip_word(bits, target, pattern):
    """Flip target exactly when all other bits match pattern.

    The target bit of pattern must be zero. For the matching-controls
    projector P, the exact word is exp(-i*pi*P/2) exp(i*pi*P*X/2).
    Its commuting Walsh expansion retains the scalar term for exact
    unitary equality. Omitting only that term changes a global phase;
    the standalone instrument is insensitive to a common global phase.
    """
    bits = _integer(bits, 1, 11, "register size")
    target = _integer(target, 0, bits, "target")
    pattern = _integer(pattern, 0, 1 << bits, "control pattern")
    if pattern >> target & 1:
        raise ValueError("the control pattern's target bit must be zero")
    controls = [bit for bit in range(bits) if bit != target]
    word = []
    for mask in range(1 << len(controls)):
        label = ["I"] * bits
        sign = 1
        for position, bit in enumerate(controls):
            if mask >> position & 1:
                label[bit] = "Z"
                sign *= (-1) ** (pattern >> bit & 1)
        angle = sign * math.pi / (1 << bits)
        word.extend(source.pauli_rotation(label, -angle))
        label[target] = "X"
        word.extend(source.pauli_rotation(label, angle))
    return word


def basis_swap_word(left, right, bits):
    """Swap two basis states with phase +1, fixing every other basis state."""
    bits = _integer(bits, 1, 11, "register size")
    left = _integer(left, 0, 1 << bits, "left endpoint")
    right = _integer(right, 0, 1 << bits, "right endpoint")
    if left == right:
        raise ValueError("distinct basis endpoints required")
    path = [left]
    for bit in range(bits):
        if (left ^ right) >> bit & 1:
            path.append(path[-1] ^ (1 << bit))
    edges = list(zip(path[:-1], path[1:], strict=True))
    word = []
    for a, b in edges + edges[-2::-1]:
        target = (a ^ b).bit_length() - 1
        word.extend(controlled_flip_word(bits, target, a & ~(1 << target)))
    return word


def repair_word(candidate, coin):
    """Local dilation: seam bits 0..2 and blank record bits starting at 3.

    The external spectator register is untouched. Coin is already committed
    to a classical record and selects one of these fixed source words.
    E performs C's repair and copies seam parity to one further record bit.
    For the repaired total-two carrier, source bit `coin` equals that parity;
    measuring the joint current/parity record is exactly contrast recording.
    """
    if candidate not in RECORD_BITS or type(coin) is not int or coin not in (0, 1):
        raise ValueError("candidate C/E and a binary integer coin required")
    if candidate == "E":
        return repair_word("C", coin) + source.cnot(coin, 5)
    return [pulse for a, b in REPAIR_SWAPS["C", coin]
            for pulse in basis_swap_word(a, b, 5)]


def preparation_word(setting):
    """Prepare one specified pure phase from the nine-qubit blank.

    A unitary word cannot prepare diagonal_control from a pure blank.
    The experiment mixes phase_0 and phase_pi with a separately charged
    source randomizer and retains that classical preparation record.
    Pure settings retain three requests, including zero-angle phase pulses.
    """
    if setting not in SETTINGS[:3]:
        raise ValueError("a pure phase setting required; diagonal_control needs a retained mixture")
    theta = {"phase_0": 0., "phase_pi": math.pi,
             "phase_pi_over_2": math.pi / 2}[setting]
    return [source.pulse("Y", -math.pi / 4, [0]),
            source.pulse("I", theta / 2, [0]),
            source.pulse("Z", -theta / 2, [0])]


def readout_word():
    """Copy the rank-one plus projector to blank pointer bit 9.

    Measuring only this pointer gives the complete Lueders instrument,
    retaining coherence throughout the minus subspace and preserving old
    records. The source basis changes are undone before pointer readout.
    """
    return ([source.pulse("Y", math.pi / 4, [0])]
            + controlled_flip_word(SOURCE_QUBITS + 1, POINTER, 0)
            + [source.pulse("Y", -math.pi / 4, [0])])
