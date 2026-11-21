"""Literal clock, counter and reader evolution with retained quantum references."""
from fractions import Fraction
from itertools import product
import math
import unittest

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1 import reader

from . import controller as c


def amplitudes(length, time):
    """Signed spin-rotation formula, independent of the numerical clock solver."""
    sine, cosine = math.sin(time / length), math.cos(time / length)
    return np.array([(-1j) ** k * math.sqrt(math.comb(length, k)) *
                     sine ** k * cosine ** (length-k) for k in range(length+1)])


def literal_clock(gates):
    """Assemble forward edges and their adjoints by clock matrix units."""
    length, dimension = len(gates), gates[0].shape[0]
    result = sparse.csr_matrix(((length+1)*dimension, (length+1)*dimension), dtype=complex)
    for k, gate in enumerate(gates):
        edge = sparse.csr_matrix(([math.sqrt((k+1)*(length-k))/length],
                                  ([k+1], [k])), shape=(length+1, length+1))
        forward = sparse.kron(edge, gate, format='csr')
        result += forward + forward.conj().T
    return result


def literal_reader_stage(flag, word, stage):
    """The four reversible word updates without the primary permutation arrays."""
    if stage == 0:
        return word ^ (2 | flag)
    if stage == 1:
        return word ^ ((2 | flag) << 2)
    if stage == 2:
        return word ^ ((word & 3) << 4)
    other, courier = (word >> 2) & 3, (word >> 4) & 3
    valid = bool(other & 2) and bool(courier & 2)
    agree = (other & 1) == (courier & 1)
    decision = 1 if valid and agree else (2 if valid else 3)
    public = 7 * (other & 1) if valid and agree else 0
    return word ^ (decision << 6) ^ (public << 8)


class ControllerTests(unittest.TestCase):
    def test_clock_spectrum_and_complete_signed_amplitudes(self):
        for length in (1, 2, 4, 7, 19, 64):
            with self.subTest(length=length):
                weights = np.array([math.sqrt((k+1)*(length-k))/length
                                    for k in range(length)])
                clock = np.diag(weights, 1) + np.diag(weights, -1)
                eigenvalues, eigenvectors = np.linalg.eigh(clock)
                np.testing.assert_allclose(c.couplings(length), weights, atol=1e-15, rtol=0)
                np.testing.assert_allclose(eigenvalues, np.linspace(-1, 1, length+1),
                                           atol=2e-15, rtol=0)
                # Include negative times and signs beyond transfer and return.
                for fraction in (-.23, 0., 1/6, .5, .73, 1., 1.25):
                    time = fraction * math.pi * length
                    spectral = eigenvectors @ (np.exp(-1j*time*eigenvalues) * eigenvectors[0])
                    exact = amplitudes(length, time)
                    np.testing.assert_allclose(c.clock_amplitudes(length, time), spectral,
                                               atol=2e-13, rtol=0)
                    np.testing.assert_allclose(spectral, exact, atol=2e-13, rtol=0)
                    np.testing.assert_allclose(c.clock_law(length, time), abs(exact)**2,
                                               atol=2e-13, rtol=0)
                    self.assertAlmostEqual(float(sum(c.clock_law(length, time))), 1., places=12)
                # The return is an operator identity on every initial clock position.
                returned = (eigenvectors * np.exp(-1j*math.pi*length*eigenvalues)) @ eigenvectors.T
                np.testing.assert_allclose(returned, (-1)**length*np.eye(length+1),
                                           atol=2e-13, rtol=0)

    def test_exact_observation_window_tail_and_full_length_numerical_tail(self):
        for length in (2, 4, 16, 64, 2672):
            with self.subTest(length=length):
                # At either endpoint p=3/4. Integer arithmetic retains every
                # unfinished position, with the strict completion cutoff L/2.
                numerator = sum(math.comb(length, k)*3**k for k in range(length//2))
                tail = Fraction(numerator, 4**length)
                markov_bound = Fraction(3, 4)**(length//2)
                self.assertLessEqual(tail, markov_bound)
                for endpoint in (1, 2):
                    time = endpoint*math.pi*length/3
                    numerical = float(c.clock_law(length, time)[:length//2].sum())
                    self.assertGreater(numerical, 0.)
                    # This checks the represented tiny probability, allowing
                    # rounded endpoint time and p. It does not certify the
                    # implementation to the exact Hamiltonian's error ceiling.
                    np.testing.assert_allclose(numerical, float(tail),
                                               rtol=5e-12, atol=0.)

    def test_numeric_domains_reject_nonreal_nonfinite_boolean_and_noninteger_inputs(self):
        for bad in (True, np.bool_(False), float('nan'), float('inf'), -float('inf'),
                    1+2j, np.complex128(1+2j), np.complex128(1), 10**400, '1', None):
            with self.subTest(real=repr(bad)):
                with self.assertRaises(ValueError):
                    c.finite(bad, 'value')
        for bad in (False, 0, -1, 1., np.int64(2), float('nan'), '2'):
            with self.subTest(integer=repr(bad)):
                for operation in (lambda: c.positive_integer(bad, 'length'),
                                  lambda: c.couplings(bad), lambda: c.clock_law(bad, 1),
                                  lambda: c.program(rounds=bad)):
                    with self.assertRaises(ValueError):
                        operation()
        for operation in (lambda: c.native_gate(float('nan'), .1),
                          lambda: c.native_gate(8, float('inf')),
                          lambda: c.clock_amplitudes(4, False),
                          lambda: c.clock_law(4, float('nan')),
                          lambda: c.program(angle=True)):
            with self.assertRaises(ValueError):
                operation()

    def test_competing_clock_laws_match_literal_matrices_and_change_transfer(self):
        length, time = 8, 4*math.pi
        for kind in ('uniform', 'broken'):
            weights = (np.full(length, .5) if kind == 'uniform' else
                       np.array([math.sqrt((k+1)*(length-k))/length for k in range(length)]))
            if kind == 'broken':
                weights[length//2-1] = 0
            matrix = np.diag(weights, 1) + np.diag(weights, -1)
            actual = c.clock_amplitudes(length, time, kind)
            expected = expm(-1j*time*matrix)[:, 0]
            np.testing.assert_allclose(actual, expected, atol=3e-14, rtol=0)
            self.assertLess(abs(actual[-1])**2, .9)
            if kind == 'broken':
                np.testing.assert_allclose(actual[length//2:], 0, atol=1e-15, rtol=0)

    def test_literal_hamiltonian_rejects_invalid_gates_and_edge_weights(self):
        for gates in ([], [np.diag([1., .9])], [np.ones((2, 2))],
                      [np.eye(2), np.eye(3)], [np.full((2, 2), np.nan)],
                      [np.array([[1., 0.], [0., np.inf]])]):
            with self.subTest(gates=str(gates)):
                with self.assertRaises(ValueError):
                    c.controller_matrix(gates)
        for weights in ([1j], [1+0j], [True], [np.nan], [np.inf], [1., 1.]):
            with self.subTest(weights=weights):
                with self.assertRaises(ValueError):
                    c.controller_matrix([np.eye(2)], weights)
        # Negative real couplings remain Hermitian; they are not a norm defect.
        matrix = c.controller_matrix([np.array([[0, 1j], [1, 0]])], [-.5])
        np.testing.assert_array_equal(matrix.toarray(), matrix.conj().T.toarray())

    def test_full_source_prefix_operator_check_retains_all_reference_columns(self):
        result = c.direct_clock_check()
        self.assertEqual(result['dimension'], 5*78)
        self.assertEqual(result['input_columns'], 78)
        self.assertEqual(len(result['operator_errors']), 3)
        self.assertLess(max(result['operator_errors']), 3e-12)

    def test_modular_counters_are_local_permutations_on_dirty_words(self):
        capacities = tuple(2 + k % 3 for k in range(12))
        for step in (c.Step('seam', (2, 9)), c.Step('reference', (7,)),
                     c.Step('reader', (0, 1), reader=0, stage=3), c.Step('padding')):
            ports = set(step.ports) if step.kind in ('seam', 'reference') else set()
            # Every touched-register input and every value of an untouched
            # register occurs, with nonzero values in the other registers.
            for word in (tuple((j+k) % cap for k, cap in enumerate(capacities)) for j in range(12)):
                expected = tuple((v+1) % capacities[k] if k in ports else v
                                 for k, v in enumerate(word))
                actual = c.counter_gate(word, step, capacities)
                self.assertEqual(actual, expected)
                self.assertEqual(c.counter_gate(actual, step, capacities, inverse=True), word)
            local_words = list(product(*(range(capacities[k]) for k in sorted(ports))))
            images = set()
            for local in local_words:
                word = [cap-1 for cap in capacities]
                for k, value in zip(sorted(ports), local):
                    word[k] = value
                output = c.counter_gate(tuple(word), step, capacities)
                images.add(output)
                self.assertEqual(c.counter_gate(output, step, capacities, inverse=True), tuple(word))
            self.assertEqual(len(images), len(local_words))
        for index, value in ((0, -1), (1, capacities[1]), (2, True)):
            word = [0]*12
            word[index] = value
            with self.assertRaises(ValueError):
                c.counter_gate(tuple(word), c.Step('seam', (0, 1)), capacities)
        for invalid in (0, False, 1.5):
            with self.assertRaises(ValueError):
                c.counter_gate((0,)*12, c.Step('reference', (0,)), (invalid,)*12)
        with self.assertRaises(ValueError):
            c.counter_gate((0,)*12, c.Step('reference', (0,)), capacities, inverse=1)

    def test_program_counts_execute_actual_increments_and_uncompute(self):
        tape = c.program()
        counts = c.count_table(tape)
        expected = np.zeros_like(counts)
        for k, step in enumerate(tape):
            expected[k+1] = expected[k]
            if step.kind in ('seam', 'reference'):
                for port in step.ports:
                    expected[k+1, port] += 1
        np.testing.assert_array_equal(counts, expected)
        np.testing.assert_array_equal(counts[-1], [192]*11 + [176])
        capacities = tuple(int(v)+1 for v in counts[-1])
        word = tuple(int(v) for v in counts[-1])
        for step in reversed(tape):
            word = c.counter_gate(word, step, capacities, inverse=True)
        self.assertEqual(word, (0,)*12)

    def test_literal_counter_tensor_source_hamiltonian(self):
        # Six actual counter words, all nine seam inputs and four clock edges.
        # Counters are tensor factors, including dirty inputs and wrapping;
        # they are not inferred from an observed clock position.
        capacities = (3, 2) + (1,)*10
        words = [(a, b) + (0,)*10 for a, b in product(range(3), range(2))]
        steps = (c.Step('seam', (0, 1)), c.Step('reference', (0,)),
                 c.Step('reader', (0, 1), reader=0, stage=0), c.Step('seam', (0, 1)))
        loads = list(product(range(3), repeat=2))
        generator = np.zeros((9, 9))
        for i, (a, b) in enumerate(loads):
            generator[i, i] += .5 + 8*(a-b)**2
            generator[loads.index((b, a)), i] -= .5
        source = (expm(-.17j*generator),
                  np.diag([np.exp(-.11j*8*a*a) for a, _ in loads]),
                  np.eye(9), expm(.09j*generator))
        gates, expected_gates = [], []
        for step, gate in zip(steps, source):
            actual = np.zeros((6, 6))
            expected = np.zeros((6, 6))
            for col, word in enumerate(words):
                actual[words.index(c.counter_gate(word, step, capacities)), col] = 1
                updated = list(word)
                if step.kind in ('seam', 'reference'):
                    for port in step.ports:
                        updated[port] = (updated[port]+1) % capacities[port]
                expected[words.index(tuple(updated)), col] = 1
            gates.append(sparse.kron(sparse.csr_matrix(actual), gate, format='csr'))
            expected_gates.append(np.kron(expected, gate))
        hamiltonian = c.controller_matrix(gates)
        np.testing.assert_array_equal(hamiltonian.toarray(), literal_clock(gates).toarray())
        # All 54 data columns establish the full operator, so their arbitrary
        # superpositions and untouched references are covered as well.
        initial = np.zeros((5*54, 54), complex)
        initial[:54] = np.eye(54)
        prefixes = [np.eye(54, dtype=complex)]
        for gate in expected_gates:
            prefixes.append(gate @ prefixes[-1])
        for time in (.61*math.pi, 2*math.pi, 4*math.pi):
            actual = expm_multiply(-1j*time*hamiltonian, initial, traceA=0)
            expected = np.vstack([a*p for a, p in zip(amplitudes(4, time), prefixes)])
            self.assertLess(np.linalg.norm(actual-expected, 2), 2e-12)
        no_counters = np.eye(6)
        wrong_final = np.kron(no_counters, source[3] @ source[2] @ source[1] @ source[0])
        # Some prefixes differ even when the chosen finite registers wrap at
        # the final edge. The second prefix retains a nontrivial counter shift.
        self.assertGreater(np.linalg.norm(prefixes[2]-np.kron(no_counters, source[1] @ source[0]), 2), 1.)
        np.testing.assert_allclose(prefixes[-1], wrong_final, atol=2e-15, rtol=0)

    def test_literal_joint_reader_clock_retains_dirty_words_and_reference_coherence(self):
        words = np.arange(reader.APP)
        dimension = 9*reader.APP
        loads = list(product(range(3), repeat=2))
        primary = [reader.stages(flag) for flag in (0, 1)]
        gates, independent_maps = [], []
        columns = np.arange(dimension)
        for stage in range(4):
            permutation, independent = [], []
            for index, (a, b) in enumerate(loads):
                flag = int(a == b == 1)
                permutation.extend(index*reader.APP + primary[flag][stage])
                independent.extend(index*reader.APP + literal_reader_stage(flag, int(w), stage)
                                   for w in words)
            independent = np.asarray(independent)
            np.testing.assert_array_equal(permutation, independent)
            gates.append(sparse.csr_matrix((np.ones(dimension), (permutation, columns)),
                                           shape=(dimension, dimension)))
            independent_maps.append(independent)
        hamiltonian = c.controller_matrix(gates)
        difference = hamiltonian-literal_clock(gates)
        self.assertEqual(difference.nnz, 0)
        self.assertEqual(hamiltonian.shape, (92160, 92160))
        apparatus = np.exp(2j*math.pi*words**2/(2*reader.APP+1))/math.sqrt(reader.APP)
        source = np.column_stack([np.exp(2j*math.pi*np.arange(9)*frequency/9)/3
                                  for frequency in (1, 2)])
        data = np.einsum('sr,w->swr', source, apparatus).reshape(dimension, 2)
        np.testing.assert_allclose(data.conj().T @ data, np.eye(2), atol=2e-14, rtol=0)
        np.testing.assert_allclose(abs(data[:, 0])**2, abs(data[:, 1])**2, atol=1e-18, rtol=0)
        self.assertEqual(np.count_nonzero(data), 2*dimension)
        initial = np.zeros((5*dimension, 2), complex)
        initial[:dimension] = data
        prefixes = [data]
        for permutation in independent_maps:
            output = np.empty_like(data)
            output[permutation] = prefixes[-1]
            prefixes.append(output)
        for time in (.71*math.pi, 2*math.pi, 4*math.pi):
            with self.subTest(time=time):
                actual = expm_multiply(-1j*time*hamiltonian, initial, traceA=0)
                expected = np.vstack([a*p for a, p in zip(amplitudes(4, time), prefixes)])
                self.assertLess(np.linalg.norm(actual-expected, 2), 2e-12)
                joint = actual.reshape(5, 9, reader.APP, 2)/math.sqrt(2)
                ideal = expected.reshape(5, 9, reader.APP, 2)/math.sqrt(2)
                # Trace only clock and apparatus; retain source plus reference.
                rho = np.einsum('kswr,ktwq->srtq', joint, joint.conj()).reshape(18, 18)
                wanted = np.einsum('kswr,ktwq->srtq', ideal, ideal.conj()).reshape(18, 18)
                np.testing.assert_allclose(rho, wanted, atol=2e-13, rtol=0)
                self.assertAlmostEqual(float(np.trace(rho).real), 1., places=12)
                # A reference-dephased substitute has the same populations
                # and fails this retained off-diagonal block.
                self.assertGreater(np.linalg.norm(rho[0::2, 1::2]), .3)

    def test_hermitian_controller_cannot_have_one_way_absorbing_halt(self):
        gate = np.array([[0, 1j], [1, 0]])
        hamiltonian = c.controller_matrix([gate]).toarray()
        run_to_halt = hamiltonian[2:, :2]
        halt_to_run = hamiltonian[:2, 2:]
        np.testing.assert_array_equal(halt_to_run, run_to_halt.conj().T)
        self.assertGreater(np.linalg.norm(run_to_halt), 1.)
        one_way = hamiltonian.copy()
        one_way[:2, 2:] = 0
        self.assertGreater(np.linalg.norm(one_way-one_way.conj().T), 1.)
        # If a Hermitian halt subspace is invariant, its reverse block also
        # vanishes, eliminating arrivals from the orthogonal running space.
        absorbing = one_way.copy()
        absorbing[2:, :2] = absorbing[:2, 2:].conj().T
        np.testing.assert_array_equal(expm(-1j*absorbing)[2:, :2], np.zeros((2, 2)))


if __name__ == '__main__':
    unittest.main()
