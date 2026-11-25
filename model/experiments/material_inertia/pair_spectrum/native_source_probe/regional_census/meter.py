"""Sharp parity and sign acquisition with a reset conditioned on the first public word."""

from fractions import Fraction
from itertools import product

import numpy as np

from .. import circuit, physical
from ..measure import check

pointer = circuit.pointer
BRIGHT, DARK = 6, 2  # load words 20, 02
OTHER = [i for i in range(9) if i not in (BRIGHT, DARK)]
PARITY, SIGN = Fraction(2, 65), Fraction(1, 65)


def read_word(targets, ratio, sign=False, omit_y=False):
    basis = pointer.neutral_basis(3, 4)
    word = list(basis)
    for target in targets:
        word += pointer.diagonal_word(target, 4, ratio)
    if sign and not omit_y:
        word.append(pointer.pulse("reference", 4, 4, Fraction(1, 64)))
    return word + pointer.inverse(basis)


def two_target_calibration(scale=1, omit_y=False):
    """Complete five-qutrit replay, with an untouched compensator."""
    inputs = []
    charges = []
    for a, b in product(range(3), repeat=2):
        q = a + b - 2
        if abs(q) <= 1:
            inputs.append([a, b, 1 - q, 2, 0])
            charges.append(q)
    states = np.zeros((243, len(inputs)), complex)
    for col, loads in enumerate(inputs):
        states[int(np.array(loads) @ [81, 27, 9, 3, 1]), col] = 1
    parity_final = physical.primitive(read_word((0, 1), scale * PARITY), states)
    words = 3 * physical.WORDS[:, 3] + physical.WORDS[:, 4]
    first = np.array(
        [np.sum(abs(parity_final[words == w]) ** 2, axis=0) for w in range(9)]
    )
    odd = parity_final.copy()
    odd[words != DARK] = 0
    reset = [pointer.pulse("contact", 3, 4, 1)]
    signed = physical.primitive(
        reset + read_word((0, 1), scale * SIGN, True, omit_y), odd
    )
    second = np.array([np.sum(abs(signed[words == w]) ** 2, axis=0) for w in range(9)])
    # Every first read other than DARK ends the program. DARK has nine possible second reads.
    hist = np.vstack([first[[w for w in range(9) if w != DARK]], second])
    check(
        np.max(abs(hist.sum(axis=0) - 1)) < 2e-12,
        "native parity/sign instrument is incomplete",
    )
    # Decoder: parity BRIGHT is zero. Sign BRIGHT is positive with this Y reference.
    decoder = [None if w != BRIGHT else 0 for w in range(9) if w != DARK] + [
        -1 if w == DARK else 1 if w == BRIGHT else None for w in range(9)
    ]
    success = [float(hist[decoder.index(q), col]) for col, q in enumerate(charges)]
    if scale == 1 and not omit_y:
        check(
            min(success) > 1 - 3e-12, "physical sign decoder failed a compensated input"
        )
    return dict(
        inputs=inputs,
        charges=charges,
        correct_probabilities=success,
        complete_pointer_probabilities=hist.T.tolist(),
        decoder=decoder,
    )


def local_blocks():
    """Inspect the native three-qutrit refocusing on every endpoint occupation."""
    rows = []
    for ratio in (PARITY, SIGN):
        word = pointer.diagonal_word(0, 2, ratio)
        block = circuit.apply(word, np.eye(27))
        loads = circuit.STATES
        diag = np.diag(circuit.law.seam_term(8))[3 * loads[:, 0] + loads[:, 2]]
        err = np.linalg.norm(
            block - np.diag(np.exp(-1j * np.pi * float(ratio) * diag)), 2
        )
        expected_ratio = []
        for n, idx in enumerate(circuit.CODE):
            relative = block[idx[1], idx[1]] / block[idx[0], idx[0]]
            want = np.exp(1j * np.pi * float(ratio) * 65 / 2 * (n - 1))
            expected_ratio.append(float(abs(relative - want)))
        check(
            max(err, *expected_ratio) < 3e-12,
            "refocused native contact or arm phase disagrees",
        )
        rows.append(
            dict(
                dwell_over_pi=str(ratio),
                full_operator_error=float(err),
                relative_arm_phase_errors=expected_ratio,
            )
        )
    return rows
