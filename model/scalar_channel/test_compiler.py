"""Independent dense small-domain oracles for the imported finite compiler."""

from itertools import product
import math
import unittest
from unittest.mock import patch

import numpy as np
from scipy.linalg import expm

from . import carrier
from . import compiler as c
from . import singlet_fusion as f


def dense_pauli(word, qubits):
    phase, x, z = word
    out = np.zeros((1 << qubits, 1 << qubits), complex)
    for column in range(1 << qubits):
        out[column ^ x, column] = 1j ** phase * (-1) ** (column & z).bit_count()
    return out


def annihilators(modes):
    out = []
    for mode in range(modes):
        matrix = np.zeros((1 << modes, 1 << modes), complex)
        for column in range(1 << modes):
            if column & (1 << mode):
                matrix[column ^ (1 << mode), column] = (-1) ** (column & ((1 << mode) - 1)).bit_count()
        out.append(matrix)
    return out


def dense_mps(state):
    value = np.ones((1, 1), complex)
    for tensor in state:
        value = np.einsum('ai,pij->apj', value, tensor).reshape(-1, tensor.shape[2])
    return value[:, 0]


def dense_sum(terms):
    return sum(amplitude * dense_mps(state) for amplitude, state in terms)


class CompilerTests(unittest.TestCase):
    def test_exact_carrier_incidence_circulation_and_cell_spectrum(self):
        data = carrier.load_source()
        self.assertEqual(data['B'].shape, (30, 120))
        self.assertEqual(len(set(map(tuple, data['wedges']))), 120)
        self.assertEqual(np.bincount(data['edges'].ravel()).tolist(), [5] * 12)
        adjacency = np.zeros((30, 30))
        for i, edge in enumerate(data['edges']):
            for j, other in enumerate(data['edges']):
                adjacency[i, j] = int(i != j and bool(set(edge) & set(other)))
        # Every seam has eight adjacent seams. Each wedge contributes half
        # the rank-one difference form, independently of the stored B.
        expected = 4 * np.eye(30) - adjacency / 2
        np.testing.assert_allclose(data['B'] @ data['B'].T, expected, atol=2e-15, rtol=0)
        for center in range(12):
            chi = carrier.chi(data, center)
            self.assertAlmostEqual(float(np.vdot(chi, chi).real), 1.0)
            np.testing.assert_allclose(data['D'] @ chi, -chi, atol=1e-15, rtol=0)
        values = np.linalg.eigvalsh(data['D'])
        self.assertEqual(int(np.sum(values < 0)), 120)
        self.assertGreaterEqual(float(np.min(abs(values))), 1 - 1e-13)
        self.assertAlmostEqual(float(values[values < 0].sum()), -156.16778455062658, places=10)
        # Mutating one generated carrier cannot alter another preparation.
        data['B'][0, 0] = 99
        self.assertNotEqual(carrier.load_source()['B'][0, 0], 99)

    def test_small_even_parity_intertwiner_against_dense_car(self):
        for links in (((0, 1), (1, 2)), ((0, 1), (0, 2), (1, 2))):
            graph = c.encode_graph(3, links)
            qubits = graph['M']
            identity = np.eye(1 << qubits)
            nums = [dense_pauli(word, qubits) for word in graph['B']]
            edge_ops = [dense_pauli(word, qubits) for word in graph['A']]
            cycles, _ = c.fundamental_cycles(graph)
            projector = identity.copy()
            for number in nums[:-1]:
                projector = projector @ ((identity + number) / 2)
            for _, _, word in cycles:
                projector = projector @ ((identity + dense_pauli(word, qubits)) / 2)
            vacuum = projector[:, np.argmax(np.linalg.norm(projector, axis=0))]
            vacuum /= np.linalg.norm(vacuum)
            xs = [dense_pauli(word, qubits) for word in c.logical_pauli_generators(graph)]
            fermions = annihilators(3)
            gammas = [a + a.conj().T for a in fermions]
            fnums = [np.eye(8) - 2 * a.conj().T @ a for a in fermions]
            fxs = []
            for i in range(2):
                operator = -1j * gammas[i] @ gammas[2]
                for j in range(i + 1, 2):
                    operator = operator @ fnums[j]
                fxs.append(operator)
            code, fock = [], []
            for word in range(4):
                left, right = vacuum.copy(), np.eye(8, dtype=complex)[:, 0]
                for bit in range(2):
                    if word & (1 << bit):
                        left, right = xs[bit] @ left, fxs[bit] @ right
                code.append(left)
                fock.append(right)
            V, W = np.column_stack(code), np.column_stack(fock)
            np.testing.assert_allclose(V.conj().T @ V, np.eye(4), atol=2e-15, rtol=0)
            np.testing.assert_allclose(W.conj().T @ W, np.eye(4), atol=2e-15, rtol=0)
            for i in range(3):
                np.testing.assert_allclose(V.conj().T @ nums[i] @ V, W.conj().T @ fnums[i] @ W, atol=2e-15)
            for k, (i, j) in enumerate(graph['links']):
                target_a = -1j * gammas[i] @ gammas[j]
                np.testing.assert_allclose(V.conj().T @ edge_ops[k] @ V, W.conj().T @ target_a @ W, atol=2e-15)
                z = .3 + .7j
                native = sum(a * dense_pauli(p, qubits) for a, p in c.hopping_words(graph, k, z))
                target = z * fermions[i].conj().T @ fermions[j] + z.conjugate() * fermions[j].conj().T @ fermions[i]
                np.testing.assert_allclose(V.conj().T @ native @ V, W.conj().T @ target @ W, atol=2e-15)
                # A reference entangled with all four code states is retained.
                reference = np.eye(4, dtype=complex).reshape(-1) / 2
                left = np.kron(np.eye(4), V.conj().T @ expm(-.37j * native) @ V) @ reference
                right = np.kron(np.eye(4), W.conj().T @ expm(-.37j * target) @ W) @ reference
                np.testing.assert_allclose(left, right, atol=3e-15, rtol=0)

    def test_tableau_measurement_branches_match_dense_projectors(self):
        graph = c.encode_graph(3, ((0, 1), (0, 2), (1, 2)))
        chord, _, cycle = c.fundamental_cycles(graph)[0][0]
        matrix = dense_pauli(cycle, 3)
        for sign in (-1, 1):
            tableau = c.Tableau(3)
            vector = (np.eye(8) + sign * matrix)[:, 0] / 2
            probability = float(np.vdot(vector, vector).real)
            self.assertEqual(tableau.measure(cycle, sign), probability)
            vector /= math.sqrt(probability)
            if sign < 0:
                correction = (0, 0, 1 << chord)
                tableau.pauli(correction)
                vector = dense_pauli(correction, 3) @ vector
            self.assertEqual(tableau.expectation(cycle), 1)
            for phase, x, z in product(range(4), range(8), range(8)):
                if (phase - (x & z).bit_count()) % 2 == 0:
                    word = (phase, x, z)
                    expected = np.vdot(vector, dense_pauli(word, 3) @ vector).real
                    self.assertAlmostEqual(tableau.expectation(word), expected)

    def test_all_routed_branches_compose_with_gaussian_pointer_instrument(self):
        from . import instruments, routed

        # Three source modes, fuel and parity reference. Seven graph bits
        # have three independent loops, so all eight sign branches fit a
        # literal calculation. Both stars participate in the routed reads.
        links = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3), (3, 4))
        graph = c.encode_graph(5, links, {edge: sum(edge) % 2 for edge in links})
        graph.update(N=2, data={'edges': [(0, 1)]})
        qubits = graph['M']
        identity = np.eye(1 << qubits)
        loops, _ = c.fundamental_cycles(graph)
        vacua, prepared, probabilities = [], [], []

        def returned_data_density(tableau):
            for star in range(2):
                self.assertEqual(tableau.expectation((0, 0, 1 << (qubits + star))), 1)
            density = identity.copy()
            for phase, x, z in tableau.rows:
                # Returned scratch is |00>; restricting its Z factors to
                # +1 gives the complete data projector, without a trace cut.
                self.assertEqual(x >> qubits, 0)
                word = (phase, x, z & ((1 << qubits) - 1))
                density = density @ (identity + dense_pauli(word, qubits)) / 2
            self.assertAlmostEqual(float(np.trace(density).real), 1)
            return density

        for signs in product((-1, 1), repeat=len(loops)):
            tableau = c.Tableau(qubits + 2)
            probability = 1
            for (chord, _, cycle), sign in zip(loops, signs):
                record = routed.measure_loop(tableau, graph, chord, cycle, sign)
                self.assertEqual(record['route'], [graph['owner'][chord],
                                                   1-graph['owner'][chord],
                                                   graph['owner'][chord]])
                probability *= record['probability']
            vacua.append(returned_data_density(tableau))
            probabilities.append(probability)
            # Fill the two source transit modes through the actual graph A.
            tableau.pauli(graph['oriented'](1, 2))
            empty_fuel_rows = list(tableau.rows)
            for fuel in (0, 1):
                tableau.rows = list(empty_fuel_rows)
                if fuel:
                    tableau.pauli(graph['oriented'](3, 4))
                prepared.append((fuel, returned_data_density(tableau)))
        self.assertEqual(probabilities, [.125] * 8)
        self.assertEqual(sum(probabilities), 1)
        for density in vacua:
            np.testing.assert_array_equal(density, vacua[0])

        # Construct the intertwiner independently in graph and literal Fock
        # coordinates. Its entire even-sector operator equality retains any
        # untouched reference, including coherence across number sectors.
        vacuum = vacua[0][:, np.argmax(np.linalg.norm(vacua[0], axis=0))].copy()
        vacuum /= np.linalg.norm(vacuum)
        fermions = annihilators(5)
        gammas = [a + a.conj().T for a in fermions]
        fock_numbers = [np.eye(32) - 2*a.conj().T @ a for a in fermions]
        graph_x = [dense_pauli(word, qubits) for word in c.logical_pauli_generators(graph)]
        fock_x = []
        for mode in range(4):
            operator = -1j * gammas[mode] @ gammas[4]
            for other in range(mode + 1, 4):
                operator = operator @ fock_numbers[other]
            fock_x.append(operator)
        native_columns, fock_columns = [], []
        for word in range(16):
            native, literal = vacuum.copy(), np.eye(32, dtype=complex)[:, 0]
            for mode in range(4):
                if word & (1 << mode):
                    native, literal = graph_x[mode] @ native, fock_x[mode] @ literal
            native_columns.append(native)
            fock_columns.append(literal)
        V, W = np.column_stack(native_columns), np.column_stack(fock_columns)
        transfer = V @ W.conj().T

        source = np.zeros((5, 5), complex)
        source[:3, :3] = [[1, .4j, .2], [-.4j, -1, .3], [.2, .3, -.7]]
        ramp = .4 * source + .6 * np.diag(np.diag(source))
        target = np.array([1, 1j, -.4]); target /= np.linalg.norm(target)
        pulse = source.copy()
        pulse[:3, 3], pulse[3, :3] = 4*target, 4*target.conj()
        phase = np.diag([0, 0, 0, 1, 0])

        def native_generator(matrix):
            return sum(value * dense_pauli(word, qubits)
                       for word, value in instruments.compile_matrix(graph, matrix).items())

        def fock_generator(matrix):
            return sum(matrix[a, b] * fermions[a].conj().T @ fermions[b]
                       for a in range(5) for b in range(5))

        program = [(ramp, .3), (pulse, math.pi/16), (source, .5),
                   (phase, math.pi/2), (pulse, math.pi/16)]
        native_U, literal_U, one_body_U = identity.copy(), np.eye(32), np.eye(5)
        for matrix, action in program:
            native_U = expm(-1j*action*native_generator(matrix)) @ native_U
            literal_U = expm(-1j*action*fock_generator(matrix)) @ literal_U
            one_body_U = expm(-1j*action*matrix) @ one_body_U
        np.testing.assert_allclose(native_U @ V,
                                   V @ (W.conj().T @ literal_U @ W), atol=2e-14, rtol=0)
        native_number = (identity - dense_pauli(graph['B'][3], qubits)) / 2
        literal_number = fermions[3].conj().T @ fermions[3]
        completeness = np.zeros((16, 16), complex)
        for native_projector, literal_projector in (
                (identity-native_number, np.eye(32)-literal_number),
                (native_number, literal_number)):
            native_K = V.conj().T @ native_projector @ native_U @ V
            literal_K = W.conj().T @ literal_projector @ literal_U @ W
            np.testing.assert_allclose(native_K, literal_K, atol=2e-14, rtol=0)
            completeness += native_K.conj().T @ native_K
        np.testing.assert_allclose(completeness, np.eye(16), atol=2e-14, rtol=0)
        for fuel, density in prepared:
            occupied_word = 6 + (24 if fuel else 0)
            expected = transfer[:, occupied_word]
            np.testing.assert_allclose(density, np.outer(expected, expected.conj()),
                                       atol=2e-14, rtol=0)
            covariance = np.diag([0, 1, 1, fuel, fuel])
            gaussian_probability = (one_body_U @ covariance @ one_body_U.conj().T)[3, 3].real
            native_probability = np.trace(native_number @ native_U @ density @ native_U.conj().T).real
            self.assertAlmostEqual(float(native_probability), float(gaussian_probability), places=13)

    def test_fusion_tensors_have_the_declared_physical_singlets(self):
        vacuum = dense_mps(f.mps_for_rank(2, 0))
        expected = np.zeros(49)
        expected[0] = 1
        np.testing.assert_allclose(vacuum, expected, atol=2e-15, rtol=0)
        singlet = dense_mps(f.mps_for_rank(2, 1))
        expected = np.zeros(49)
        expected[[4 * 7 + 4, 5 * 7 + 5, 6 * 7 + 6]] = 1 / math.sqrt(3)
        self.assertAlmostEqual(abs(np.vdot(expected, singlet)), 1.0)
        J = np.zeros((7, 7), complex)
        J[4, 5], J[5, 4] = -1j, 1j
        np.testing.assert_allclose((np.kron(J, np.eye(7)) + np.kron(np.eye(7), J)) @ singlet,
                                   np.zeros(49), atol=1e-15, rtol=0)
        self.assertAlmostEqual(float(np.linalg.norm(np.kron(J, np.eye(7)) @ singlet) ** 2), 2 / 3)
        for r in range(1, 7):
            self.assertEqual(f.completions(2 * r, 0), math.comb(2 * r, r) // (r + 1))
            for rank in range(f.completions(2 * r, 0)):
                self.assertEqual(f.rank_path(f.unrank(r, rank)), rank)
        self.assertLess(f.schur_orthogonality_local(), 2e-15)

    def test_raw_gate_retains_unused_singlet_color_and_quantum_reference(self):
        r = 5
        code = [f.mps_for_rank(r, k) for k in range(32)]
        V = np.column_stack([dense_mps(state) for state in code])
        np.testing.assert_allclose(V.conj().T @ V, np.eye(32), atol=3e-15, rtol=0)
        colored = f.product_mps([1, 0, 0, 0, 0])
        unused = f.mps_for_rank(r, 32)
        weak = f.product_mps([4, 0, 0, 0, 0])
        branches = [[(.5, code[3]), (.5j, colored)], [(.5, unused), (.5j, weak)]]
        permutation = np.eye(32)[:, np.arange(32) ^ 1]
        actual, expected = [], []
        for terms in branches:
            raw = dense_sum(terms)
            transformed, _, _ = f.apply_complete_gate_to_mps_sum(r, terms, lambda v: f.logical_pauli(v, xmask=1))
            actual.append(dense_sum(transformed))
            expected.append(raw + V @ (permutation - np.eye(32)) @ (V.conj().T @ raw))
        actual, expected = np.concatenate(actual), np.concatenate(expected)
        np.testing.assert_allclose(actual, expected, atol=3e-15, rtol=0)
        self.assertAlmostEqual(float(np.vdot(actual, actual).real), 1.0)
        self.assertAlmostEqual(abs(np.vdot(actual, expected)), 1.0)
        # A raw owner permutation does not preserve this truncated fusion code.
        permuted = V[:, 31].reshape((7,) * r).swapaxes(0, 4).reshape(-1)
        leakage = np.linalg.norm(permuted - V @ (V.conj().T @ permuted)) ** 2
        self.assertGreater(leakage, .01)

    def test_invalid_domains_and_checks_survive_optimization(self):
        bad_calls = [
            lambda: carrier.load_source('sphere'), lambda: carrier.chi(carrier.load_source(), True),
            lambda: f.unrank(0, 0), lambda: f.unrank(5, True), lambda: f.unrank(5, 42),
            lambda: f.rank_path((1, 0, 1)), lambda: f.rank_path((0, -1, 0)),
            lambda: f.rank_path((0, 1, 2)), lambda: f.cg_tensor(0, 0),
            lambda: f.mps_for_rank(5, 1.5), lambda: f.product_mps([7]),
            lambda: f.product_amplitude(5, 0, [0]),
            lambda: f.logical_pauli([1, 0, 0]), lambda: f.logical_pauli([1, 0], xmask=2),
            lambda: f.logical_pauli([float('nan'), 0]),
            lambda: f.overlap(f.mps_for_rank(5, 0), f.mps_for_rank(6, 0)),
            lambda: f.charted_product_amplitude(5, 0, [0] * 5, [0] * 5),
            lambda: f.apply_complete_gate_to_mps_sum(2, [(1, f.product_mps([0, 0]))], lambda v: v),
            lambda: f.apply_complete_gate_to_mps_sum(5, [(1, f.mps_for_rank(5, 0))], lambda v: 2 * v),
            lambda: c.encode_graph(3, ((0, 1),)), lambda: c.encode_graph(2, ((0, 1), (0, 1))),
            lambda: c.encode_graph(2, ((1, 0),)), lambda: c.Tableau(0),
            lambda: c.Tableau(2).measure((0, 1, 0), 0),
            lambda: c.Tableau(2).measure((1, 1, 0), 1),
            lambda: c.Tableau(2).expectation((0, 4, 0)),
            lambda: c.Tableau(2).measure((0, 0, 1), -1),
            lambda: list(c.bits(-1)), lambda: c.gf2_rank([-1]),
            lambda: c.hopping_words(c.encode_graph(2, ((0, 1),)), 0, float('nan')),
        ]
        for index, action in enumerate(bad_calls):
            with self.subTest(index=index), self.assertRaises(ValueError):
                action()
        with patch.object(c, 'gf2_rank', return_value=0), self.assertRaisesRegex(ValueError, 'Compiler invariant'):
            c.run()

    def test_full_cell_replay_is_explicit_about_execution_scope(self):
        result = c.run(include_artifacts=True)
        self.assertEqual((result['total_modes'], result['total_BKSF_qubits'], result['independent_loop_constraints']),
                         (174, 432, 259))
        self.assertEqual(len(result['preparation_trace']['loop_measurements']), 259)
        self.assertEqual(len(result['operator_manifest']['complete_logical_X_basis']), 173)
        self.assertTrue(result['theorem_scope']['even_parity_algebra_and_local_instruments_intertwine'])
        self.assertFalse(result['execution_boundary']['local_ancilla_and_ack_service_executed'])
        self.assertFalse(result['execution_boundary']['physical_owner_reset_executed'])
        self.assertGreater(result['full_raw_domain_energy']['conservative_invalid_energy_gap_bound'], 1)


if __name__ == '__main__':
    unittest.main()
