"""Literal Fock/pointer checks independent of covariance propagation."""

from itertools import combinations, product
import math
import unittest

import numpy as np
from scipy.linalg import expm

from . import compiler as c
from . import instruments as i


def fermions(n):
    result = []
    for mode in range(n):
        operator = np.zeros((2 ** n, 2 ** n), complex)
        for word in range(2 ** n):
            if word & (1 << mode):
                operator[word ^ (1 << mode), word] = (-1) ** (word & ((1 << mode) - 1)).bit_count()
        result.append(operator)
    return result


def second_quantize(matrix, annihilators):
    result = np.zeros_like(annihilators[0])
    for a, ca in enumerate(annihilators):
        for b, cb in enumerate(annihilators):
            result += matrix[a, b] * ca.conj().T @ cb
    return result


class InstrumentTests(unittest.TestCase):
    def test_literal_128_dimension_pointer_dilation_matches_all_records(self):
        graph = c.encode_graph(5, tuple(combinations(range(5), 2)))
        graph.update(nsource=3, F=[3], R=[4])
        source = np.zeros((5, 5), complex)
        source[:3, :3] = [[1, .4j, .2], [-.4j, -1, .3], [.2, .3, -.7]]
        target = np.array([1, 1j, -.4], complex)
        target /= np.linalg.norm(target)
        frame, _ = np.linalg.qr(np.array([[1, .3j], [.2j, 1], [.5, -.7j]], complex))
        covariance = frame @ frame.conj().T
        rows, checks = i.native_family(graph, source, covariance, target, queries=(0, 4, 8), strength=4.0)
        operators = fermions(5)
        pulse = source.copy()
        pulse[:3, 3], pulse[3, :3] = 4 * target, 4 * target.conj()
        U = expm(-1j * math.pi / 16 * second_quantize(pulse, operators))
        H = second_quantize(source, operators)
        V = second_quantize(pulse - source, operators)
        # A coherent two-particle Slater state, including the two separately
        # retained fuel/reference choices inside the same even parity block.
        source_state = np.zeros(32, complex)
        for a, b in combinations(range(3), 2):
            source_state[(1 << a) | (1 << b)] = np.linalg.det(frame[[a, b], :])
        pointer_unitary = np.zeros((128, 128))
        for pointers, word in product(range(4), range(32)):
            flag = (word >> 3) & 1
            output = pointers ^ (3 * flag)
            pointer_unitary[32 * output + word, 32 * pointers + word] = 1
        np.testing.assert_array_equal(pointer_unitary.T @ pointer_unitary, np.eye(128))
        largest_pulse_work, largest_second_switch_on = 0.0, 0.0
        for row in rows:
            state = source_state.copy()
            if row['initial_fuel']:
                state = operators[3].conj().T @ operators[4].conj().T @ state
            phase = np.diag([np.exp(-1j * math.pi * row['phase'] / 2 * ((word >> 3) & 1)) for word in range(32)])
            first = U @ state
            free = expm(-1j * row['query'] / 8 * H) @ first
            phased = phase @ free
            final = U @ phased
            literal = pointer_unitary @ np.r_[final, np.zeros(96)]
            probabilities = (abs(literal.reshape(4, 32)) ** 2).sum(axis=1)
            np.testing.assert_allclose(probabilities, row['pointer_probabilities'], atol=3e-14, rtol=0)
            ref = sum(abs(final[word]) ** 2 for word in range(32) if word & 16)
            self.assertAlmostEqual(float(ref), row['initial_reference'], places=13)
            energies = [float(np.vdot(vector, H @ vector).real)
                        for vector in (state, first, free, phased, final)]
            energies.append(float(np.vdot(literal, np.kron(np.eye(4), H) @ literal).real))
            ledger = row['source_energy_ledger']
            self.assertEqual([stage['stage'] for stage in ledger['stages']],
                             ['pulse_1', 'free_interval', 'fuel_phase', 'pulse_2', 'terminal_pointer'])
            for index, stage in enumerate(ledger['stages']):
                self.assertAlmostEqual(stage['source_energy_before'], energies[index], places=13)
                self.assertAlmostEqual(stage['source_energy_after'], energies[index + 1], places=13)
                self.assertAlmostEqual(stage['external_driving_work'], energies[index + 1] - energies[index], places=13)
            for stage_index, before, after in ((0, state, first), (3, phased, final)):
                stage = ledger['stages'][stage_index]
                on = float(np.vdot(before, V @ before).real)
                off = -float(np.vdot(after, V @ after).real)
                self.assertAlmostEqual(stage['switch_on_work'], on, places=13)
                self.assertAlmostEqual(stage['switch_off_work'], off, places=13)
                self.assertAlmostEqual(stage['external_driving_work'], on + off, places=13)
                largest_pulse_work = max(largest_pulse_work, abs(on + off))
                if stage_index == 3:
                    largest_second_switch_on = max(largest_second_switch_on, abs(on))
            self.assertAlmostEqual(ledger['total_external_driving_work'], energies[-1] - energies[0], places=13)
        self.assertLess(checks['decoded_matrix_error'], 1e-13)
        self.assertGreater(largest_pulse_work, .01)
        self.assertGreater(largest_second_switch_on, .01)

    def test_two_source_even_state_ledger_keeps_source_on(self):
        graph = c.encode_graph(4, tuple(combinations(range(4), 2)))
        graph.update(nsource=2, F=[2], R=[3])
        source = np.zeros((4, 4), complex)
        source[:2, :2] = [[.8, .35 + .2j], [.35 - .2j, -.6]]
        target = np.array([1, 1j]) / math.sqrt(2)
        # Both source modes are occupied, a pure even Slater preparation.
        occupied = 1.0
        rows, _ = i.native_family(graph, source, occupied * np.eye(2), target,
                                  queries=(0, 7), strength=2.0)
        annihilators = fermions(4)
        H = second_quantize(source, annihilators)
        pulse = source.copy()
        pulse[:2, 2], pulse[2, :2] = 2 * target, 2 * target.conj()
        V = second_quantize(pulse - source, annihilators)
        U = expm(-1j * math.pi / 8 * (H + V))
        converter_only = expm(-1j * math.pi / 8 * V)
        source_on_effect, free_state_change = 0.0, 0.0
        for row in rows:
            density = np.zeros((16, 16), complex)
            fuel_reference = 12 * row['initial_fuel']
            density[fuel_reference, fuel_reference] = 1 - occupied
            density[fuel_reference + 3, fuel_reference + 3] = occupied
            first = U @ density @ U.conj().T
            free_gate = expm(-1j * row['query'] / 8 * H)
            free = free_gate @ first @ free_gate.conj().T
            phase = np.diag([np.exp(-1j * math.pi * row['phase'] / 2 * ((word >> 2) & 1))
                             for word in range(16)])
            phased = phase @ free @ phase.conj().T
            final = U @ phased @ U.conj().T
            stages = row['source_energy_ledger']['stages']
            energies = [float(np.trace(H @ state).real) for state in (density, first, free, phased, final)]
            for index, stage in enumerate(stages[:4]):
                self.assertAlmostEqual(stage['source_energy_before'], energies[index], places=13)
                self.assertAlmostEqual(stage['source_energy_after'], energies[index + 1], places=13)
            self.assertAlmostEqual(energies[2], energies[1], places=13)
            self.assertEqual(stages[1]['external_driving_work'], 0.0)
            self.assertEqual(stages[1]['duration'], row['query'] / 8)
            wrong = converter_only @ density @ converter_only.conj().T
            source_on_effect = max(source_on_effect, abs(float(np.trace(H @ wrong).real) - energies[1]))
            free_state_change = max(free_state_change, float(np.linalg.norm(free - first)))
        self.assertGreater(source_on_effect, 1e-5)
        self.assertGreater(free_state_change, .01)

    def test_finite_source_on_pulses_have_detuned_transfer_probabilities(self):
        graph = c.encode_graph(3, tuple(combinations(range(3), 2)))
        graph.update(nsource=1, F=[1], R=[2])
        source = np.diag([1, 0, 0]).astype(complex)
        strength = 64.0
        rows, _ = i.native_family(graph, source, np.zeros((1, 1)),
                                  np.array([1]), queries=(0,), strength=strength)
        row = next(row for row in rows if row['initial_fuel'] == 1 and row['phase'] == 0)
        stages = row['source_energy_ledger']['stages']
        # Independently solve the two-level source/fuel problem. A source
        # energy of one makes its energy equal its excitation probability.
        frequency = math.sqrt(strength ** 2 + .25)
        def transfer(area):
            return strength ** 2 / frequency ** 2 * math.sin(area / strength * frequency) ** 2
        half = transfer(math.pi / 4)
        full = transfer(math.pi / 2)
        self.assertAlmostEqual(stages[0]['source_energy_after'], half, places=13)
        self.assertAlmostEqual(stages[3]['source_energy_after'], full, places=13)
        self.assertGreater(abs(half - .5), 1e-6)
        self.assertGreater(abs(full - 1), 1e-5)
        np.testing.assert_allclose(row['pointer_probabilities'],
                                   [full, 0, 0, 1 - full], atol=3e-14, rtol=0)

    def test_source_preparation_rejects_odd_or_fractional_covariances(self):
        graph = c.encode_graph(5, tuple(combinations(range(5), 2)))
        graph.update(nsource=3, F=[3], R=[4])
        source = np.diag([1, -1, -1, 0, 0]).astype(complex)
        target = np.array([1, 0, 0])
        with self.assertRaisesRegex(ValueError, 'even particle number'):
            i.native_family(graph, source, np.diag([1, 0, 0]), target)
        with self.assertRaisesRegex(ValueError, 'pure Slater projector'):
            i.native_family(graph, source, np.diag([.5, .5, 0]), target)
        # Even a physically admissible even-parity mixture needs a different
        # preparation witness; this helper implements pure Slater sources.
        with self.assertRaisesRegex(ValueError, 'pure Slater projector'):
            i.native_family(graph, source, np.diag([.37, .37, 0]), target)
        # Entrywise roundoff accumulates in a trace. The projector bound
        # controls that accumulation without accepting an odd-rank source.
        _, checks = i.native_family(graph, source, np.diag([1+1.1e-11, 1+1.1e-11, 0]),
                                    target, queries=(0,))
        self.assertEqual(checks['source_particle_number'], 2)
        self.assertGreater(checks['source_particle_number_error'], i.TOL)
        self.assertLessEqual(checks['source_particle_number_error'],
                             checks['source_particle_number_tolerance'])

    def test_work_ledger_rejects_unmodeled_apparatus_energy(self):
        graph = c.encode_graph(4, tuple(combinations(range(4), 2)))
        graph.update(nsource=2, F=[2], R=[3])
        for a, b in ((2, 2), (3, 3), (0, 2), (1, 3)):
            source = np.zeros((4, 4))
            source[a, b] = source[b, a] = 1e-30
            with self.subTest(a=a, b=b), self.assertRaisesRegex(ValueError, 'only on source modes'):
                i.native_family(graph, source, np.eye(2), np.array([1, 0]))

    def test_pointer_gate_is_reversible_on_dirty_memories_and_retains_failures(self):
        for flags in product((False, True), repeat=2):
            for fuel, word in product(range(2), range(4)):
                output = i.pointer_permutation(fuel, word, *flags)
                self.assertEqual(i.pointer_permutation(fuel, output, *flags), word)
        self.assertEqual(i.pointer_law(.4), [.6, 0.0, 0.0, .4])
        self.assertEqual(i.pointer_law(.4, read_b=False), [.6, .4, 0.0, 0.0])
        self.assertEqual(i.pointer_law(.4, word=1), [0.0, .6, .4, 0.0])
        # Copying orthogonal fuel records dephases that fuel if pointers are
        # ignored; the full gate still preserves its external reference.
        rng = np.random.default_rng(82)
        vector = rng.normal(size=16) + 1j * rng.normal(size=16)
        vector /= np.linalg.norm(vector)
        gate = np.zeros((8, 8))
        for fuel, word in product(range(2), range(4)):
            gate[fuel + 2 * i.pointer_permutation(fuel, word), fuel + 2 * word] = 1
        joint = np.kron(np.eye(2), gate)
        np.testing.assert_allclose(joint.conj().T @ joint @ vector, vector, atol=1e-15, rtol=0)

    def test_compile_decode_detects_unsupported_and_malformed_laws(self):
        graph = c.build()
        source = i.source_matrix(graph)
        decoded, residual = i.decode_paulis(graph, i.compile_matrix(graph, source))
        np.testing.assert_array_equal(decoded, source)
        self.assertEqual(residual, 0.0)
        wrong = source.copy()
        wrong[0, 1] = wrong[1, 0] = 1e-30
        with self.assertRaisesRegex(ValueError, 'outside fixed local graph'):
            i.compile_matrix(graph, wrong)
        wrong = source.copy()
        wrong[0, 30] += .1j
        with self.assertRaisesRegex(ValueError, 'Hermitian'):
            i.compile_matrix(graph, wrong)
        bad = [lambda: i.prepare_source(True), lambda: i.prepare_source(0),
               lambda: i.pointer_law(-1e-20), lambda: i.pointer_law(float('nan')),
               lambda: i.pointer_permutation(True, 0), lambda: i.pointer_permutation(0, 4),
               lambda: i.source_matrix(graph, 'unknown'), lambda: i.run(1),
               lambda: i.decode_paulis(graph, {(0, 0, 0): float('nan')}),
               lambda: i.decode_paulis(graph, {(1, 1, 0): 1.0})]
        for index, call in enumerate(bad):
            with self.subTest(index=index), self.assertRaises(ValueError):
                call()

    def test_default_replay_keeps_complete_outcomes_and_same_apparatus_controls(self):
        result = i.run()
        self.assertEqual(result['preparation']['steps_per_star'], 16)
        self.assertFalse(result['preparation']['manual_preparation_reproduced'])
        self.assertFalse(result['preparation']['exact_sea_eigenvectors_used_to_initialize'])
        self.assertEqual(result['schema'], 'oph.scalar_channel.finite_instruments.v2')
        accounting = result['source_work_accounting']
        self.assertFalse(accounting['controller_energy_modeled'])
        self.assertFalse(accounting['autonomous_energy_conservation_claimed'])
        self.assertFalse(accounting['phase_and_pointer_durations_declared'])
        self.assertFalse(accounting['includes_preparation_ramp_work'])
        self.assertEqual(len(result['cases']), 12)
        for case in result['cases']:
            self.assertEqual(len(case['records']), 24)
            self.assertEqual({(row['query'], row['phase'], row['initial_fuel']) for row in case['records']},
                             set(product((0, 4, 8), range(4), (0, 1))))
            for row in case['records']:
                probabilities = row['pointer_probabilities']
                self.assertEqual(len(probabilities), 4)
                self.assertGreaterEqual(min(probabilities), 0)
                self.assertAlmostEqual(sum(probabilities), 1.0, places=14)
                self.assertEqual(probabilities[1:3], [0.0, 0.0])
                ledger = row['source_energy_ledger']
                stages = ledger['stages']
                self.assertLessEqual(ledger['work_balance_residual'], i.TOL)
                self.assertAlmostEqual(ledger['total_external_driving_work'],
                                       ledger['net_source_energy_change'], delta=i.TOL)
                for index, stage in enumerate(stages):
                    self.assertLessEqual(stage['work_balance_residual'], i.TOL)
                    if index:
                        self.assertEqual(stage['source_energy_before'], stages[index - 1]['source_energy_after'])
                    if index in (0, 3):
                        self.assertEqual(stage['source_strength'], 1.0)
                        self.assertEqual(stage['converter_strength'], 64.0)
                        self.assertEqual(stage['duration'], math.pi / 256)
                        self.assertEqual(stage['converter_action'], math.pi / 4)
                    else:
                        self.assertEqual(stage['external_driving_work'], 0.0)
                        self.assertEqual(stage['converter_action'], 0.0)
                    if index in (2, 4):
                        self.assertIsNone(stage['duration'])
                        self.assertIsNone(stage['source_action'])
            self.assertLessEqual(case['checks']['parity_reference_error'], i.TOL)
            self.assertEqual(case['checks']['source_particle_number'], 120)
            self.assertLessEqual(case['checks']['source_projector_error'], i.TOL)
            self.assertLessEqual(case['checks']['largest_source_work_balance_residual'], i.TOL)
            self.assertTrue(case['checks']['source_commutes_with_fuel_phase_and_pointer'])
        for control in result['controls'].values():
            self.assertGreater(control['largest_pointer_probability_difference'], .01)
            self.assertTrue(control['same_preparation_targets_pulses_and_readouts'])


if __name__ == '__main__':
    unittest.main()
