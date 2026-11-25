"""Complete ternary native source circuit, with every off-code row retained."""

from fractions import Fraction
from itertools import product

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.experiments.interval_echo_consensus.signed_effects import fixed_pointer as pointer

STATES = np.array(list(product(range(3), repeat=3)))
# Local labels: target0, source park1, source cross2; contacts0-2 and1-2.
CODE = [[n * 9 + 6, n * 9 + 2] for n in range(3)]


def check(value, message):
    if not value:
        raise ArithmeticError(message)


def apply(word, columns):
    out = np.asarray(columns, complex).copy()
    for kind, a, b, ratio in word:
        time = float(ratio) * np.pi
        if kind == 'reference':
            out *= np.exp(-8j * time * STATES[:, a]**2)[:, None]
        else:
            u = expm(-1j * time * law.seam_term(8))
            result = np.zeros_like(out)
            for src, row in enumerate(STATES):
                col = law.SRC.index(tuple(row[[a, b]]))
                for dest in range(9):
                    if abs(u[dest, col]) > 0:
                        target = row.copy()
                        target[[a, b]] = law.SRC[dest]
                        index = int(target @ [9, 3, 1])
                        result[index] += u[dest, col] * out[src]
            out = result
    return out


def audit():
    H = apply(pointer.neutral_basis(1, 2), np.eye(27))
    h = H[np.ix_(CODE[1], CODE[1])]
    check(np.linalg.norm(h.conj().T @ h - np.eye(2)) < 2e-13, 'source basis leaks')
    inverse = apply(pointer.inverse(pointer.neutral_basis(1, 2)), np.eye(27))
    check(np.linalg.norm(inverse - H.conj().T, 2) < 3e-12, 'native source inverse differs')
    yword = [pointer.pulse('reference', 2, 2, Fraction(1, 64))]
    physical_y = apply(yword, np.eye(27))
    quadrature = physical_y[np.ix_(CODE[1], CODE[1])]
    check(np.linalg.norm(quadrature - np.diag([1, -1j]), 2) < 3e-13,
          'native degree-five Y reference differs')
    meters, diagnostics = {}, []
    for multiplier in (1, 2):
        word = pointer.diagonal_word(0, 2, Fraction(8 * multiplier, 31))
        diagonal = apply(word, np.eye(27))
        expected = np.diag(np.exp(-1j * float(Fraction(8 * multiplier, 31)) * np.pi
                                  * np.diag(law.seam_term(8))[3 * STATES[:, 0] + STATES[:, 2]]))
        check(np.linalg.norm(diagonal - expected, 2) < 3e-12, 'full ternary refocusing differs')
        phases = [diagonal[np.ix_(indices, indices)] for indices in CODE]
        relative = [d / phases[1][0, 0] for d in phases]
        check(np.linalg.norm(relative[1] - np.eye(2)) < 3e-12, 'host source reference differs')
        check(np.linalg.norm(relative[0] @ relative[2] - np.eye(2)) < 3e-12,
              'common charge-square phase is not cancelled')
        meters[multiplier] = {}
        for charge in (-1, 0, 1):
            phase = np.linalg.matrix_power(relative[2], charge)
            meters[multiplier][charge] = {}
            for name, extra in (('X', np.eye(2)), ('Y', quadrature)):
                state = h.conj().T @ extra @ phase @ h[:, 0]
                # Complete pointer words0..8, no discarded off-code outcome.
                p = np.zeros(9)
                p[6], p[2] = abs(state)**2
                check(abs(sum(p) - 1) < 3e-12, 'source meter loses a word')
                meters[multiplier][charge][name] = p
        bad = apply(pointer.diagonal_word(0, 2, Fraction(8 * multiplier, 31), False), H[:, [v[0] for v in CODE]])
        leakage = [float(1 - np.sum(abs(bad[indices, i])**2)) for i, indices in enumerate(CODE)]
        diagnostics.append(dict(multiplier=multiplier, native_pulses=len(word),
                                full_operator_error=float(np.linalg.norm(diagonal - expected, 2)),
                                extreme_refocus_omission_leakage=leakage,
                                sharp_calibration={str(q): {axis: p.tolist() for axis, p in setting.items()}
                                                   for q, setting in meters[multiplier].items()}))
    return meters, diagnostics
