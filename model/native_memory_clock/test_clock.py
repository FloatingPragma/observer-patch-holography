"""Hostile finite checks of the literal rail and complete native instruments.

These tests are exploratory validation, not confirmatory run receipts. The
independent oracle below uses incidence-times-flux and never calls forward,
inverse, native_event, step or column when deriving the event history.
"""

from collections import defaultdict
from dataclasses import replace
from itertools import combinations, product
from math import comb, cos, pi, sin, sqrt
from pathlib import Path
import subprocess
import sys
import unittest

import numpy as np

from . import clock


_EDGES = tuple(combinations(range(4), 2))
_ADDRESSES = tuple(pair for edge in _EDGES for pair in (edge, edge[::-1]))


def _event(state, address, background=(0, 0, 0, 0), law="native"):
    """Separate primary-rule calculation, with both neutral states retained."""
    q = list(background)
    for e, (a, b) in enumerate(_EDGES):
        flux = ((state >> (2 * e + 1)) & 1) - ((state >> (2 * e)) & 1)
        q[a] += flux
        q[b] -= flux
    tail, head = _ADDRESSES[address]
    old = (state >> address) & 1
    difference = q[tail] - q[head]
    if law == "uphill":
        difference = -difference
    changes = law != "hold" and ((old == 0 and difference > 0)
                                 or (old == 1 and difference < 0))
    return state ^ ((1 << address) if changes else 0), old


def _reference_preparation(c):
    # Nonorthogonal reference vectors leave detectable source coherences.
    coefficients = {0: np.array([0.5, 0.5]),
                    1: np.array([0.5, 0.5j])}
    states, h, _ = clock.prepared_span(c, {0: 1 / sqrt(2), 1: 1 / sqrt(2)})
    initial = np.zeros((len(states), 2), complex)
    for i, b in enumerate(states):
        if b.clock == 0:
            initial[i] = coefficients[b.buffers[0]] / sqrt(12 ** c.horizon)
    return states, h, initial, coefficients


def _expected_vector(c, states, coefficients, tau):
    index = {b: i for i, b in enumerate(states)}
    result = np.zeros((len(states), 2), complex)
    for source, reference in coefficients.items():
        for addresses in product(range(12), repeat=c.horizon):
            current, records, sign = source, [0] * c.horizon, 1
            for n in range(c.horizon + 1):
                buffers = [0] * (c.horizon + 1)
                buffers[n] = current
                b = clock.Basis(n, tuple(buffers), addresses, tuple(records))
                coefficient = (sqrt(comb(c.horizon, n))
                               * (-1j * sin(tau / c.horizon)) ** n
                               * cos(tau / c.horizon) ** (c.horizon - n))
                result[index[b]] = sign * coefficient * reference / sqrt(12 ** c.horizon)
                if n < c.horizon:
                    current, old = _event(current, addresses[n], c.background, c.law)
                    records[n] = old
                    if c.phase and old:
                        sign *= -1
    return result


def _conditional_density(states, amplitudes):
    # Trace unread future addresses; preserve source and reference coherences.
    fibers = defaultdict(dict)
    for b, reference in zip(states, amplitudes):
        history = tuple(zip(b.addresses[:b.clock], b.records[:b.clock]))
        key = (b.clock, history, b.addresses[b.clock:])
        source = b.buffers[b.clock]
        for r, amplitude in enumerate(reference):
            fibers[key][(source, r)] = amplitude
    result = defaultdict(complex)
    for (n, history, _), entries in fibers.items():
        for left, a in entries.items():
            for right, b in entries.items():
                result[(n, history, left, right)] += a * b.conjugate()
    return result


def _native_density(c, coefficients, tau):
    result = defaultdict(complex)
    for n in range(c.horizon + 1):
        branches = defaultdict(lambda: defaultdict(complex))
        for addresses in product(range(12), repeat=n):
            for source, reference in coefficients.items():
                current, history, sign = source, [], 1
                for address in addresses:
                    current, old = _event(current, address, c.background, c.law)
                    history.append((address, old))
                    if c.phase and old:
                        sign *= -1
                for r, amplitude in enumerate(reference):
                    branches[tuple(history)][(current, r)] += sign * amplitude / sqrt(12 ** n)
        probability = (comb(c.horizon, n) * sin(tau / c.horizon) ** (2 * n)
                       * cos(tau / c.horizon) ** (2 * (c.horizon - n)))
        for history, entries in branches.items():
            for left, a in entries.items():
                for right, b in entries.items():
                    result[(n, history, left, right)] += probability * a * b.conjugate()
    return result


class ClockTests(unittest.TestCase):
    def test_dirty_columns_inverse_phase_hermiticity_and_patch_support(self):
        rng = np.random.default_rng(84201)
        for law, phase, background in product(("native", "hold", "uphill"), (False, True),
                                              ((0, 0, 0, 0), (1, 0, 0, 0))):
            c = clock.Construction(4, background=background, law=law, phase=phase)
            for _ in range(20):
                b = clock.Basis(int(rng.integers(5)), tuple(map(int, rng.integers(4096, size=5))),
                                tuple(map(int, rng.integers(12, size=4))),
                                tuple(map(int, rng.integers(2, size=4))))
                self.assertEqual(len(clock.column(c, b)), 1 if b.clock in (0, 4) else 2)
                for destination, weight in clock.column(c, b).items():
                    direction = destination.clock - b.clock
                    recovered, back_phase = clock.step(c, destination, -direction)
                    again, phase_value = clock.step(c, b, direction)
                    self.assertEqual(again, destination)
                    self.assertEqual(recovered, b)
                    self.assertEqual(back_phase * phase_value, 1)
                    self.assertEqual(clock.column(c, destination)[b], weight)
                    j = min(b.clock, destination.clock)
                    self.assertEqual(destination.addresses, b.addresses)
                    for k in range(5):
                        if k not in (j, j + 1):
                            self.assertEqual(destination.buffers[k], b.buffers[k])
                    for k in range(4):
                        if k != j:
                            self.assertEqual(destination.records[k], b.records[k])
                    # Vary the entire complement, checking dependence as well
                    # as the absence of writes outside the two-cell support.
                    complement = replace(
                        b,
                        buffers=tuple(value if k in (j, j + 1) else value ^ 4095
                                      for k, value in enumerate(b.buffers)),
                        addresses=tuple(value if k == j else (value + 5) % 12
                                        for k, value in enumerate(b.addresses)),
                        records=tuple(value if k == j else value ^ 1
                                      for k, value in enumerate(b.records)),
                    )
                    varied, varied_phase = clock.step(c, complement, direction)
                    self.assertEqual(varied_phase, phase_value)
                    self.assertEqual(varied.buffers[j:j + 2], destination.buffers[j:j + 2])
                    self.assertEqual(varied.records[j], destination.records[j])

    def test_disjoint_actual_data_gates_commute_including_phases(self):
        c = clock.Construction(4, phase=True)
        rng = np.random.default_rng(95301)

        def gate(b, j):
            output, phase = clock.step(c, replace(b, clock=j), 1)
            return replace(output, clock=0), phase

        for _ in range(64):
            b = clock.Basis(0, tuple(map(int, rng.integers(4096, size=5))),
                            tuple(map(int, rng.integers(12, size=4))),
                            tuple(map(int, rng.integers(2, size=4))))
            for j, k in ((0, 2), (0, 3), (1, 3)):
                left, p = gate(b, j)
                left, q = gate(left, k)
                right, r = gate(b, k)
                right, s = gate(right, j)
                self.assertEqual(left, right)
                self.assertEqual(p * q, r * s)

    def test_literal_evolution_reference_and_full_native_density(self):
        for horizon in (1, 2, 3):
            for phase in (False, True):
                c = clock.Construction(horizon, phase=phase)
                states, h, initial, coefficients = _reference_preparation(c)
                state_set = set(states)
                self.assertEqual(len(states), 2 * (horizon + 1) * 12 ** horizon)
                for b in states:
                    self.assertTrue(set(clock.column(c, b)).issubset(state_set))
                    self.assertTrue(all(x == 0 for k, x in enumerate(b.buffers) if k != b.clock))
                    self.assertTrue(all(x == 0 for x in b.records[b.clock:]))
                tau = 0.731 * horizon
                actual = clock.evolve(h, initial, tau)
                expected = _expected_vector(c, states, coefficients, tau)
                np.testing.assert_allclose(actual, expected, atol=2e-13, rtol=2e-12)
                marginal = np.array([sum(np.vdot(row, row).real
                                         for b, row in zip(states, actual) if b.clock == n)
                                     for n in range(horizon + 1)])
                np.testing.assert_allclose(marginal, clock.clock_probabilities(horizon, tau), atol=2e-12)
                observed = _conditional_density(states, actual)
                oracle = _native_density(c, coefficients, tau)
                keys = observed.keys() | oracle.keys()
                self.assertLess(max(abs(observed[k] - oracle[k]) for k in keys), 3e-13)
                self.assertGreater(sum(abs(value) for (_, _, left, right), value in oracle.items()
                                       if left != right), 0.01)
                np.testing.assert_allclose(clock.evolve(h, actual, -tau), initial, atol=2e-12)
                returned = clock.evolve(h, initial, pi * horizon)
                np.testing.assert_allclose(returned, (-1) ** horizon * initial, atol=3e-12)

    def test_full_carrier_gauge_on_dirty_orbits_and_spectral_multiplicity(self):
        # Full dimension and multiplicities are integer theorem consequences.
        # The enormous full product matrix is not allocated or claimed executed.
        for horizon in (1, 2, 3, 7):
            c = clock.Construction(horizon, phase=True)
            resource = c.resources()
            data_dimension = 2 ** (12 * (horizon + 1)) * 24 ** horizon
            self.assertEqual(int(resource["full_data_dimension"]), data_dimension)
            self.assertEqual(int(resource["full_hilbert_dimension"]), (horizon + 1) * data_dimension)
            values, vectors = np.linalg.eigh(clock.clock_matrix(c))
            expected = np.array([-1 + 2 * k / horizon for k in range(horizon + 1)])
            np.testing.assert_allclose(values, expected, atol=6e-15)
            multiplicities = [data_dimension] * len(expected)
            self.assertEqual(sum(multiplicities), int(resource["full_hilbert_dimension"]))
            self.assertEqual(expected[0], -1)
            self.assertEqual(expected[-1], 1)
            if horizon % 2 == 0:
                self.assertEqual(expected[horizon // 2], 0)
            b = clock.Basis(0, tuple((103 * j + 4095) % 4096 for j in range(horizon + 1)),
                            tuple((5 * j + 3) % 12 for j in range(horizon)),
                            tuple(j % 2 for j in range(horizon)))
            orbit, signs = [b], [1]
            for _ in range(horizon):
                b, phase = clock.step(c, b, 1)
                orbit.append(b)
                signs.append(signs[-1] * phase)
            for eigenvalue, vector in zip(values, vectors.T):
                eigenstate = {b: signs[n] * vector[n] for n, b in enumerate(orbit)}
                action = defaultdict(complex)
                for b, amplitude in eigenstate.items():
                    for output, value in clock.column(c, b).items():
                        action[output] += value * amplitude
                self.assertEqual(set(action), set(eigenstate))
                self.assertLess(max(abs(action[b] - eigenvalue * amplitude)
                                    for b, amplitude in eigenstate.items()), 3e-15)

    def test_interior_clock_measurement_retains_failed_return_outcomes(self):
        for horizon in (1, 2, 3):
            c = clock.Construction(horizon)
            states, h, initial = clock.prepared_span(c, {1: 1.0})
            read_time, total_time = pi * horizon / 4, pi * horizon
            outcomes = clock.measured_return(c, states, h, initial, read_time, total_time)
            self.assertAlmostEqual(float(outcomes.sum()), 1.0, places=12)
            self.assertAlmostEqual(float(outcomes[0]), comb(2 * horizon, horizon) / 4 ** horizon, places=12)
            self.assertGreater(float(outcomes[1:].sum()), 0.4)
            # An endpoint read is an uninformative dephasing control.
            endpoint = clock.measured_return(c, states, h, initial, pi * horizon / 2, total_time)
            self.assertAlmostEqual(float(endpoint[0]), 1.0, places=12)

    def test_phase_completion_preserves_basis_records_and_changes_interference(self):
        cases = []
        for phase in (False, True):
            c = clock.Construction(1, phase=phase)
            states, h, initial = clock.prepared_span(c, {0: 1 / sqrt(2), 1: 1 / sqrt(2)})
            final = clock.evolve(h, initial, pi / 2)
            cases.append((states, final, clock.readout(states, final)))
        self.assertEqual(cases[0][0], cases[1][0])
        keys = cases[0][2].keys() | cases[1][2].keys()
        self.assertLess(max(abs(cases[0][2].get(k, 0) - cases[1][2].get(k, 0)) for k in keys), 1e-14)
        # Fixed source+archive interference effect: baseline endpoint vector.
        projector = cases[0][1]
        baseline_probability = abs(np.vdot(projector, cases[0][1])) ** 2
        alternative_probability = abs(np.vdot(projector, cases[1][1])) ** 2
        self.assertAlmostEqual(baseline_probability, 1.0, places=12)
        self.assertGreater(baseline_probability - alternative_probability, 0.1)

    def test_uniform_clock_competing_profile_has_own_complete_spectrum(self):
        c = clock.Construction(3, profile="uniform")
        values, vectors = np.linalg.eigh(clock.clock_matrix(c))
        np.testing.assert_allclose(values, sorted(cos(k * pi / 5) for k in range(1, 5)), atol=2e-15)
        states, h, initial = clock.prepared_span(c, {1: 1.0})
        tau = 3 * pi / 2
        actual = clock.clock_marginal(states, clock.evolve(h, initial, tau), 3)
        expected_amplitudes = vectors @ (np.exp(-1j * tau * values) * vectors[0].conjugate())
        np.testing.assert_allclose(actual, abs(expected_amplitudes) ** 2, atol=3e-13)
        self.assertLess(actual[-1], 0.999)
        self.assertGreater(float(actual[:-1].sum()), 0.001)

    def test_inputs_fail_without_implicit_normalization(self):
        for value in (0, -1, True, 1.5, "1", None):
            with self.assertRaises(ValueError):
                clock.Construction(value)
        for background in ((2, 0, 0, 0), (False, 0, 0, 0), [0, 0, 0, 0], (0, 0, 0)):
            with self.assertRaises(ValueError):
                clock.Construction(1, background=background)
        c = clock.Construction(1)
        for amplitudes in ({}, {0: 0}, {0: 2}, {0: np.nan}, {0: np.inf}, {4096: 1}, {True: 1}, {0.0: 1}):
            with self.subTest(amplitudes=amplitudes), self.assertRaises(ValueError):
                clock.prepared_span(c, amplitudes)
        valid = clock.Basis(0, (1, 4095), (0,), (1,))
        for malformed in (replace(valid, clock=True), replace(valid, clock=2),
                          replace(valid, buffers=(1, -1)), replace(valid, buffers=[1, 0]),
                          replace(valid, records=(2,)), replace(valid, addresses=(12,))):
            with self.assertRaises(ValueError):
                clock.column(c, malformed)
        for direction in (True, 0, 2, -2, 1.0):
            with self.assertRaises(ValueError):
                clock.step(c, valid, direction)
        with self.assertRaises(ValueError):
            clock.step(c, valid, -1)
        states, h, initial = clock.prepared_span(c, {1: 1.0})
        for tau in (np.nan, np.inf, 1j):
            with self.assertRaises(ValueError):
                clock.evolve(h, initial, tau)
        small = initial / 3
        output = clock.evolve(h, small, 0.3)
        self.assertAlmostEqual(float(np.vdot(output, output).real), 1 / 9, places=14)

    def test_optimized_python_real_scientific_guards_reject_mutations(self):
        root = Path(__file__).resolve().parents[2]
        programs = (
            ("from model.native_memory_clock.clock import *\n"
             "prepared_span(Construction(1), {0: 2.0})\n", "unit norm"),
            ("from model.native_memory_clock import clock\n"
             "original = clock.column\n"
             "clock.column = lambda c,b: {} if b.clock == 1 else original(c,b)\n"
             "clock.prepared_span(clock.Construction(1), {0: 1.0})\n", "non-Hermitian literal rail"),
            ("from model.native_memory_clock import clock\n"
             "original = clock.column\n"
             "def corrupted(c,b):\n"
             "    out = original(c,b)\n"
             "    if b.clock == 0:\n"
             "        out[clock.Basis(1,(9,9),(0,),(0,))] = 0.25\n"
             "    return out\n"
             "clock.column = corrupted\n"
             "clock.prepared_span(clock.Construction(1), {0: 1.0})\n", "nonzero complementary output"),
        )
        for program, message in programs:
            result = subprocess.run([sys.executable, "-O", "-c", program], cwd=root,
                                    text=True, capture_output=True, check=False)
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn(message, result.stderr)


if __name__ == "__main__":
    unittest.main()
