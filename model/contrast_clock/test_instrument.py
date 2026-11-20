"""Independent dense operators for the complete sequential reader instrument."""

from functools import lru_cache
from itertools import combinations_with_replacement, product
import unittest

import numpy as np
from scipy.linalg import expm

from model.contrast_memory_spherical_v1.carrier import CELL_SEAMS

from . import controller


def literal_states():
    states = []
    for a, b in combinations_with_replacement(range(12), 2):
        state = [0] * 12
        state[a] += 1
        state[b] += 1
        states.append(state)
    return np.asarray(states, dtype=int)


@lru_cache(maxsize=3)
def dense_sequential_kraus(kappa):
    """Literal qutrit generators and matrix exponentials, with no gate helpers."""
    states = literal_states()
    index = {tuple(state): i for i, state in enumerate(states)}
    terms = []
    for a, b in CELL_SEAMS:
        generator = np.zeros((78, 78), dtype=complex)
        for column, state in enumerate(states):
            generator[column, column] += .5 + kappa * (int(state[a]) - int(state[b])) ** 2
            swapped = list(state)
            swapped[a], swapped[b] = swapped[b], swapped[a]
            generator[index[tuple(swapped)], column] -= .5
        terms.append(generator)
    terms.extend(np.diag(kappa * states[:, port] ** 2) for port in range(12))

    half_steps = [expm(-1j * generator / 16) for generator in terms[:-1]]
    full_step = expm(-1j * terms[-1] / 8)
    period = np.eye(78, dtype=complex)
    for gate in half_steps + [full_step] + half_steps[::-1]:
        period = gate @ period
    between_readings = np.linalg.matrix_power(period, 8)

    first_flag = (states[:, 0] == 1) & (states[:, 1] == 1)
    second_flag = (states[:, 0] == 1) & (states[:, 5] == 1)
    return {
        (first, second): np.diag(second_flag == second) @ between_readings
        @ np.diag(first_flag == first) @ between_readings
        for first, second in product((0, 1), repeat=2)
    }


def fresh_reader_word(flag, control):
    """Independent complete-word encoding, including private failure records."""
    value = flag ^ int(control == 'wrong_agreement')
    observer_a = 2 | value
    observer_b = 0 if control == 'unavailable' else 2 | value
    courier = observer_a ^ int(control == 'mismatch')
    decision = {'nominal': 1, 'mismatch': 2, 'unavailable': 3, 'wrong_agreement': 1}[control]
    public = 7 * value if decision == 1 else 0
    return observer_a | (observer_b << 2) | (courier << 4) | (decision << 6) | (public << 8)


class SequentialInstrumentTests(unittest.TestCase):
    def test_all_source_columns_match_the_independent_two_reading_instrument(self):
        np.testing.assert_array_equal(controller.source_sector(), literal_states())
        cases = [(8., control) for control in ('nominal', 'mismatch', 'unavailable', 'wrong_agreement')]
        cases += [(0., 'nominal'), (-8., 'nominal')]
        for kappa, control in cases:
            with self.subTest(kappa=kappa, control=control):
                expected = dense_sequential_kraus(kappa)
                actual, _ = controller.run_history(controller.program(), kappa, control, full_operator=True)
                expected_words = set()
                for flags, operator in expected.items():
                    words = tuple(fresh_reader_word(flag, control) for flag in flags)
                    expected_words.add(words)
                    self.assertIn(words, actual)
                    # Full Kraus-operator equality covers coherences and every
                    # untouched quantum reference, rather than basis probabilities.
                    self.assertLess(np.linalg.norm(actual[words] - operator, 2), 3e-12)
                for words, operator in actual.items():
                    if words not in expected_words:
                        self.assertEqual(float(np.linalg.norm(operator, 2)), 0.)
                completeness = sum(operator.conj().T @ operator for operator in expected.values())
                self.assertLess(np.linalg.norm(completeness - np.eye(78), 2), 3e-12)


if __name__ == '__main__':
    unittest.main()
