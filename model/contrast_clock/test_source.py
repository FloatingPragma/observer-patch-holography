"""Independent source-domain checks for the finite contrast controller."""

from itertools import combinations_with_replacement, product
import math
import unittest
from unittest.mock import patch
import warnings

import numpy as np
from scipy.linalg import expm

from model.contrast_memory_spherical_v1 import carrier, law, reader
from model.contrast_memory_spherical_v1.sectors import MultisetSector, PairSector, SeamEvents

from . import controller


def raw_memory_law(kappa, exchange=True):
    """Literal 99-state write/phase/unwrite from the declared labels."""
    loads = list(product(range(3), repeat=2))
    labels = list(product((0, 1), range(-2, 3)))
    write = np.zeros((90, 9))
    for col, (a, b) in enumerate(loads):
        for coin in (0, 1):
            low, odd = divmod(a + b, 2)
            output = (low + coin * odd, low + (1 - coin) * odd)
            row = 9 * labels.index((coin, a - b)) + loads.index(output)
            write[row, col] = 1 / math.sqrt(2)
    echo = np.zeros((99, 99))
    echo[:9, 9:] = write.T
    echo[9:, :9] = write
    echo[9:, 9:] = np.eye(90) - write @ write.T
    memory = np.zeros((99, 99))
    for label, (coin, contrast) in enumerate(labels):
        for load in range(9):
            col = 9 + 9 * label + load
            memory[col, col] = kappa * contrast ** 2 + .5 * exchange
            row = 9 + 9 * labels.index((coin, -contrast)) + load
            memory[row, col] -= .5 * exchange
    return echo @ memory @ echo.T, echo


def literal_embedding(local, seam, states):
    """Embed by literal qutrit labels without the controller's index tables."""
    lookup = {tuple(state): index for index, state in enumerate(states)}
    result = np.zeros((len(states), len(states)), dtype=complex)
    a, b = seam
    for col, state in enumerate(states):
        local_col = 3 * int(state[a]) + int(state[b])
        for left, right in product(range(3), repeat=2):
            amplitude = local[3 * left + right, local_col]
            if abs(amplitude) < 1e-14:
                continue
            output = list(state)
            output[a], output[b] = left, right
            if tuple(output) not in lookup:
                raise ValueError('local operator leaves the declared full two-load sector')
            result[lookup[tuple(output)], col] += amplitude
    return result


def literal_seam_generator(kappa, shifted=False):
    result = np.zeros((9, 9))
    for a, b in product(range(3), repeat=2):
        col = 3 * a + b
        result[col, col] += .5 + kappa * (-2 * a * b if shifted else (a - b) ** 2)
        result[3 * b + a, col] -= .5
    return result


class SourceTests(unittest.TestCase):
    def test_native_gate_rejects_unrepresentable_phases_without_warnings(self):
        for kappa, angle in ((1e308, 1.), (-1e308, 1.), (1e308, 1e308)):
            with self.subTest(kappa=kappa, angle=angle), warnings.catch_warnings():
                warnings.simplefilter('error', RuntimeWarning)
                with self.assertRaisesRegex(ValueError, 'phase'):
                    controller.native_gate(kappa, angle)
        # Large finite parameters remain valid when every required phase fits.
        for kappa, angle in ((np.finfo(float).max / 8, 1.), (0., 1e308)):
            gate = controller.native_gate(kappa, angle)
            self.assertTrue(np.isfinite(gate).all())
            np.testing.assert_allclose(gate.conj().T @ gate, np.eye(9), atol=5e-16, rtol=0)

    def test_raw_99_gate_reduces_on_all_nine_columns_and_is_unitary_when_dirty(self):
        blank = np.eye(99)[:, :9]
        source_parts = law.local_term_parts()
        for kappa, angle in ((0., .125), (8., .125), (-8., .2), (8., -.2)):
            with self.subTest(kappa=kappa, angle=angle):
                generator, echo = raw_memory_law(kappa)
                np.testing.assert_allclose(generator, source_parts[0] + kappa * source_parts[1],
                                           atol=2e-13, rtol=0)
                np.testing.assert_allclose(echo, law.echo(), atol=1e-15, rtol=0)
                gate = expm(-1j * angle * generator)
                reduced = controller.native_gate(kappa, angle)
                np.testing.assert_allclose(gate @ blank, blank @ reduced, atol=3e-13, rtol=0)
                np.testing.assert_allclose(gate.conj().T @ gate, np.eye(99), atol=3e-13, rtol=0)
                np.testing.assert_allclose(echo @ echo, np.eye(99), atol=2e-15, rtol=0)
                # Every dirty input column is retained, including sectors which
                # cannot be represented by a blank-memory load-only simulation.
                self.assertEqual(gate.shape, (99, 99))
                self.assertGreater(np.linalg.norm(gate[9:, 9:] - np.eye(90)), .01)

    def test_doubled_load_is_swapped_as_one_qutrit_and_reference_has_no_exchange(self):
        kappa, angle = 8., .125
        gate = controller.native_gate(kappa, angle)
        phase = np.exp(-1j * angle)
        expected = np.zeros(9, complex)
        expected[6] = np.exp(-4j * kappa * angle) * (1 + phase) / 2
        expected[2] = np.exp(-4j * kappa * angle) * (1 - phase) / 2
        np.testing.assert_allclose(gate[:, 6], expected, atol=2e-14, rtol=0)
        self.assertEqual(float(abs(gate[4, 6])), 0.)

        reference, _ = raw_memory_law(kappa, exchange=False)
        zero_reference = np.eye(99)[:, [0, 3, 6]]
        actual = expm(-1j * angle * reference) @ zero_reference
        expected = zero_reference * np.exp(-1j * angle * kappa * np.array([0, 1, 4]))
        np.testing.assert_allclose(actual, expected, atol=2e-14, rtol=0)
        # Including the exchange in a degree/reference pair leaks load into
        # the reference. The primary source explicitly excludes this term.
        self.assertGreater(abs(gate[1, 3]) ** 2, .001)

    def test_sector_and_literal_qutrit_embedding_keep_all_78_states(self):
        states = np.asarray(controller.source_sector())
        expected = []
        for a, b in combinations_with_replacement(range(12), 2):
            state = [0] * 12
            state[a] += 1
            state[b] += 1
            expected.append(tuple(state))
        self.assertEqual(states.shape, (78, 12))
        self.assertEqual(set(map(tuple, states)), set(expected))
        self.assertEqual(sum(2 in state for state in states), 12)
        for seam, kappa, angle in (((0, 1), 8., .125), ((1, 5), -8., .2),
                                   ((10, 11), 0., -.25)):
            with self.subTest(seam=seam, kappa=kappa, angle=angle):
                local = expm(-1j * angle * literal_seam_generator(kappa))
                literal = literal_embedding(local, seam, states)
                np.testing.assert_allclose(controller.embed_local(local, seam, states), literal,
                                           atol=3e-14, rtol=0)
                np.testing.assert_allclose(controller.event_matrix(seam, kappa, angle), literal,
                                           atol=3e-13, rtol=0)
                np.testing.assert_allclose(literal.conj().T @ literal, np.eye(78),
                                           atol=3e-14, rtol=0)

    def test_embedding_rejects_incomplete_duplicate_and_malformed_source_bases(self):
        states = controller.source_sector()
        duplicate = states.copy()
        duplicate[-1] = duplicate[0]
        wrong_total = states.copy()
        wrong_total[0] = 0
        wrong_load = states.copy()
        wrong_load[0] = 0
        wrong_load[0, :2] = (-1, 3)
        noninteger = states.astype(float)
        noninteger[0, 0] = .5
        invalid = (states[:-1], states[:, :-1], np.vstack((states, states[0])),
                   states[np.all(states <= 1, axis=1)], duplicate, wrong_total,
                   wrong_load, noninteger, states.astype(bool), states.astype(complex),
                   np.full((78, 12), np.nan))
        for basis in invalid:
            with self.subTest(shape=basis.shape, dtype=str(basis.dtype)):
                with self.assertRaisesRegex(ValueError, 'source basis'):
                    controller.embed_local(np.eye(9), (0, 1), basis)
        # A reordered complete basis is valid and carries the same operator.
        order = np.random.default_rng(1209).permutation(78)
        local = expm(-.17j * literal_seam_generator(8.))
        canonical = literal_embedding(local, (0, 1), states)
        reordered = controller.embed_local(local, (0, 1), states[order])
        np.testing.assert_array_equal(reordered, canonical[np.ix_(order, order)])

    def test_degree_completion_and_doubled_energy_reference(self):
        states = np.asarray(controller.source_sector())
        car, kappa = carrier.cell(), 8.
        native = sum(literal_embedding(literal_seam_generator(kappa), seam, states)
                     for seam in car.seams)
        degree = np.array([sum((6 - int(car.degree[v])) * state[v] ** 2 for v in range(12))
                           for state in states])
        source = MultisetSector(car, 2)
        order = [source.index[tuple(state)] for state in states]
        complete = source.hamiltonian(kappa)[np.ix_(order, order)]
        np.testing.assert_array_equal(native + kappa * np.diag(degree), complete)
        shifted = sum(literal_embedding(literal_seam_generator(kappa, shifted=True), seam, states)
                      for seam in car.seams)
        square_sum = np.sum(states ** 2, axis=1)
        np.testing.assert_array_equal(shifted, complete - 6 * kappa * np.diag(square_sum))
        wrong = complete - 12 * kappa * np.eye(78)
        self.assertAlmostEqual(np.linalg.norm(shifted - wrong, 2), 12 * kappa)
        doubled = np.flatnonzero(np.any(states == 2, axis=1))
        ports = [int(np.flatnonzero(states[index] == 2)[0]) for index in doubled]
        expected = car.laplacian().toarray()[np.ix_(ports, ports)] / 2 + 24 * kappa * np.eye(12)
        np.testing.assert_array_equal(complete[np.ix_(doubled, doubled)], expected)

    def test_shifted_binary_event_is_not_the_native_event_on_overlapping_seams(self):
        car, kappa, angle = carrier.cell(), 8., .2
        sector = PairSector(car)
        events = SeamEvents(sector, kappa)
        full = np.asarray(controller.source_sector())
        rows = []
        lookup = {tuple(state): index for index, state in enumerate(full)}
        for a, b in zip(sector.aa, sector.bb):
            state = [0] * 12
            state[a] = state[b] = 1
            rows.append(lookup[tuple(state)])
        native = np.zeros(66, complex)
        native[sector.index(0, 5)] = 1
        shifted = native.copy()
        for seam in ((0, 1), (1, 5), (0, 1), (1, 5)):
            native_gate = controller.event_matrix(seam, kappa, angle)[np.ix_(rows, rows)]
            shifted_gate = events.matrix(car.seams.index(seam), angle)
            local_number = full[rows, seam[0]] + full[rows, seam[1]]
            correction = np.diag(np.exp(1j * kappa * angle * local_number))
            np.testing.assert_allclose(correction @ native_gate, shifted_gate, atol=3e-14, rtol=0)
            native = native_gate @ native
            shifted = shifted_gate @ shifted
        target = sector.index(0, 1)
        self.assertGreater(abs(native[target]) ** 2, .016)
        self.assertLess(abs(shifted[target]) ** 2, .0003)
        self.assertGreater(np.abs(abs(native) ** 2 - abs(shifted) ** 2).sum() / 2, .03)

    def test_program_actions_preserve_dirty_archives_and_reference_columns(self):
        states = np.asarray(controller.source_sector())
        rng = np.random.default_rng(417)
        blocks = {words: rng.normal(size=(78, 3)) + 1j * rng.normal(size=(78, 3))
                  for words in ((17, 2047), (2047, 39))}
        seam = (0, 1)
        local = expm(-.125j * literal_seam_generator(8.))
        literal = literal_embedding(local, seam, states)
        steps = ((controller.Step('seam', seam, .125), literal),
                 (controller.Step('reference', (5,), .125),
                  np.diag(np.exp(-1j * states[:, 5] ** 2))))
        for step, gate in steps:
            actual = controller.gate_action(step, blocks, 8., states, None)
            self.assertEqual(set(actual), set(blocks))
            for words, initial in blocks.items():
                np.testing.assert_allclose(actual[words], gate @ initial, atol=2e-14, rtol=0)
            inverse = controller.Step(step.kind, step.ports, -step.angle)
            restored = controller.gate_action(inverse, actual, 8., states, None)
            for words, initial in blocks.items():
                np.testing.assert_allclose(restored[words], initial, atol=2e-14, rtol=0)
        permutations = controller.reader_permutations('nominal')
        for stage in range(4):
            step = controller.Step('reader', seam, reader=1, stage=stage)
            actual = controller.gate_action(step, blocks, 8., states, permutations)
            restored = controller.gate_action(step, actual, 8., states, permutations)
            for words in set(restored) | set(blocks):
                np.testing.assert_allclose(restored.get(words, np.zeros((78, 3))),
                                           blocks.get(words, np.zeros((78, 3))), atol=0, rtol=0)

    def test_controller_rechecks_primary_reader_permutations_under_optimization(self):
        broken = [np.arange(2048) for _ in range(4)]
        broken[2][0] = 1
        with patch.object(controller, 'stages', return_value=broken):
            with self.assertRaisesRegex(ValueError, 'reversible involution'):
                controller.reader_permutations('nominal')

    def test_primary_reader_all_dirty_words_and_controls_without_assert_statements(self):
        # Test all combinations of supplied switches, beyond the four named
        # controls. These checks execute under python -O as unittest methods.
        words = np.arange(2048)
        for values in product((False, True), repeat=5):
            controls = dict(zip(('read_a', 'read_b', 'communicate', 'flip_message', 'common_flip'), values))
            for flag in (0, 1):
                stages = reader.stages(flag, **controls)
                actual = words.copy()
                for stage in stages:
                    np.testing.assert_array_equal(stage[stage], words)
                    actual = stage[actual]
                expected = []
                for word in range(2048):
                    a, b, courier = word & 3, (word >> 2) & 3, (word >> 4) & 3
                    decision, public = (word >> 6) & 3, (word >> 8) & 7
                    published = 2 | (flag ^ controls['common_flip'])
                    a ^= published if controls['read_a'] else 0
                    b ^= published if controls['read_b'] else 0
                    courier ^= a if controls['communicate'] else 0
                    courier ^= int(controls['flip_message'])
                    valid = bool(b & 2) and bool(courier & 2)
                    agree = (b & 1) == (courier & 1)
                    decision ^= 1 if valid and agree else (2 if valid else 3)
                    public ^= 7 * (b & 1) if valid and agree else 0
                    expected.append(a | (b << 2) | (courier << 4) | (decision << 6) | (public << 8))
                np.testing.assert_array_equal(actual, expected)
                np.testing.assert_array_equal([reader.gate(flag, int(word), **controls) for word in words], expected)
                np.testing.assert_array_equal([reader.gate(flag, int(word), inverse=True, **controls)
                                               for word in actual], words)
        # Agreement alone does not identify the source predicate: both readers
        # can share the same wrong value while every validity bit is set.
        wrong = reader.gate(0, 0, common_flip=True)
        self.assertEqual((wrong >> 6) & 3, 1)
        self.assertEqual((wrong >> 8) & 7, 7)
        unavailable = reader.gate(1, 0, communicate=False)
        self.assertEqual((unavailable >> 6) & 3, 3)
        self.assertEqual((unavailable >> 8) & 7, 0)


if __name__ == '__main__':
    unittest.main()
