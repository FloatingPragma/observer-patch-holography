"""Independent finite-clock checks and pre-freeze numerical controls."""

from dataclasses import replace
from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import product
from math import comb, sqrt
import random
import unittest
from unittest.mock import patch

import numpy as np
from scipy.sparse.linalg import expm_multiply

from . import packet as clock
from .program import (Data, Program, edge, gate_support, prefixes, preparation,
                      resources)


def _fraction(record):
    return Fraction(int(record["numerator"]), int(record["denominator"]))


def _exact_moments(order):
    # Integer autocorrelations of the literal binomial packet, extended by zero.
    weights = [comb(order, j) for j in range(order + 1)]
    z = sum(value * value for value in weights)
    mean = Fraction(sum(j * value * value for j, value in enumerate(weights)), z)
    second = Fraction(sum(j * j * value * value for j, value in enumerate(weights)), z)
    overlap_one = sum(weights[j] * weights[j + 1] for j in range(order))
    overlap_two = sum(weights[j] * weights[j + 2] for j in range(order - 1))
    velocity = Fraction(overlap_one, z)
    velocity_second = Fraction(1, 2) + Fraction(overlap_two, 2 * z)
    mixed = sum((Fraction(2 * j + 1, 2) * weights[j] * weights[j + 1]
                 for j in range(order)), Fraction()) / z
    return (mean, second - mean * mean, velocity,
            velocity_second - velocity * velocity, mixed - mean * velocity)


def _dirty(rng):
    return Data(tuple(rng.randrange(4096) for _ in range(5)),
                tuple(rng.randrange(12) for _ in range(4)),
                tuple(rng.randrange(4) for _ in range(3)),
                rng.randrange(4), rng.randrange(4),
                tuple(rng.randrange(2) for _ in range(2)), rng.randrange(12))


def _registers(data):
    values = {f"S{k}": v for k, v in enumerate(data.buffers)}
    values.update({f"A{k}": v for k, v in enumerate(data.addresses)})
    values.update({f"C{k + 1}": v for k, v in enumerate(data.couriers)})
    values.update(PA=data.pointer_a, PB=data.pointer_b, b1=data.old_bits[0],
                  b2=data.old_bits[1], D=data.second_address)
    return values


def _pointer_energy(states, h, vector):
    labels = {(data.pointer_a, data.pointer_b) for _, data in states}
    return sum(float(np.vdot(branch, h @ branch).real)
               for label in labels
               for branch in [vector * np.array([(data.pointer_a, data.pointer_b) == label
                                                 for _, data in states])])


class PacketTests(unittest.TestCase):
    def test_integer_binomial_moments_and_zero_symmetric_covariance(self):
        for order in (1, 2, 3, 4, 8, 16, 31, 64, 127):
            direct = _exact_moments(order)
            self.assertEqual(direct[:4], clock.moments(order))
            self.assertEqual(direct[4], 0)
            self.assertEqual(sum(comb(order, j) ** 2 for j in range(order + 1)),
                             comb(2 * order, order))

    def test_dst_evolution_against_literal_open_rail_and_complete_spectrum(self):
        for rail in (clock.Rail(48, 8, 4), clock.Rail(96, 16, 8)):
            h = clock.clock_matrix(rail)
            expected_spectrum = np.linalg.eigvalsh(h.toarray())
            actual_spectrum = clock.spectrum(rail)
            np.testing.assert_allclose(sorted(actual_spectrum), expected_spectrum, atol=1e-14, rtol=0)
            self.assertEqual(len(set(actual_spectrum)), rail.sites)
            full_dimension = rail.sites * int(resources()["full_data_dimension"])
            self.assertEqual(sum(int(resources()["full_data_dimension"]) for _ in actual_spectrum),
                             full_dimension)
            for direction, tau in product((1, -1), (0, 2.375, 17.31)):
                initial = clock.packet(rail, direction)
                expected = expm_multiply(-1j * tau * h, initial, traceA=0)
                actual = clock.evolve_clock(rail, tau, direction)
                np.testing.assert_allclose(actual, expected, atol=6e-14, rtol=2e-12)
                self.assertAlmostEqual(float(np.vdot(actual, actual).real), 1, places=13)
        full_rail = clock.Rail()
        all_values = clock.spectrum(full_rail)
        self.assertEqual(len(all_values), 4096)
        self.assertTrue(np.all(np.diff(all_values) < 0))
        self.assertLess(all_values[-1], 0)
        self.assertGreater(all_values[0], 0)
        self.assertEqual(len(all_values) * int(resources()["full_data_dimension"]),
                         full_rail.sites * 2 ** 72 * 12 ** 5)

    def test_literal_joint_columns_match_gauge_for_complex_sources_and_controls(self):
        sources = {0: 1 / 3, 1: 1j / 3, 17: (1 + 1j) / 3, 4095: (2 - 1j) / 3}
        configs = (Program(), Program(background=(1, 0, 0, 0), law="uphill"),
                   Program(communication=False), Program(write_a=False, feedback=False))
        for rail, config in product((clock.Rail(48, 8, 4), clock.Rail(96, 16, 8)), configs):
            states, h, initial = clock.prepared_span(rail, config, sources)
            self.assertEqual(len(states), rail.sites * 12 * len(sources))
            self.assertEqual(h.nnz, 2 * (rail.sites - 1) * 12 * len(sources))
            self.assertEqual((h - h.T.conjugate()).nnz, 0)
            self.assertAlmostEqual(float(np.vdot(initial, initial).real), 1, places=13)
            actual = clock.literal_evolve(h, initial, 11.375)
            expected = clock.factored_vector(rail, config, sources, states, 11.375)
            np.testing.assert_allclose(actual, expected, atol=1e-13, rtol=2e-12)
            np.testing.assert_allclose(clock.literal_evolve(h, actual, -11.375), initial,
                                       atol=2e-13, rtol=2e-12)

    def test_dirty_joint_hermiticity_and_actual_local_support(self):
        rng = random.Random(30491)
        rail = clock.Rail(48, 8, 4)
        configs = (Program(), Program(law="uphill"), Program(communication=False),
                   Program(feedback=False), Program(write_a=False))
        positions = (0, rail.start - 1, rail.start, rail.start + 1,
                     rail.start + 2, rail.start + 3, rail.completed, rail.sites - 1)
        for config, position in product(configs, positions):
            for _ in range(6):
                data = _dirty(rng)
                basis = (position, data)
                outputs = clock.column(rail, config, basis)
                self.assertEqual(len(outputs), 1 if position in (0, rail.sites - 1) else 2)
                for destination, value in outputs.items():
                    self.assertEqual(value, 0.5)
                    self.assertEqual(clock.column(rail, config, destination)[basis], 0.5)
                    target, changed = destination
                    gate = min(position, target) - rail.start
                    before, after = _registers(data), _registers(changed)
                    support = set(gate_support(gate)) if 0 <= gate < 4 else set()
                    for name in before:
                        if name not in support:
                            self.assertEqual(before[name], after[name])
                # A forward program followed by its actual adjoints restores
                # all populated memories and inactive source buffers.
                orbit = prefixes(config, data)
                restored = orbit[-1]
                for gate in reversed(range(4)):
                    restored = edge(config, restored, gate, -1)
                self.assertEqual(restored, data)

    def test_reducing_span_rejects_complement_and_missing_reverse_column(self):
        rail, config = clock.Rail(48, 8, 4), Program()
        original = clock.column

        def complement(r, c, basis):
            out = original(r, c, basis)
            if basis[0] == 0:
                data = replace(basis[1], pointer_a=3)
                out[(1, data)] = 0.25
            return out

        with patch.object(clock, "column", complement), self.assertRaisesRegex(ValueError, "outside declared reducing span"):
            clock.prepared_span(rail, config, {1: 1.0})

        def no_reverse(r, c, basis):
            out = original(r, c, basis)
            return {destination: value for destination, value in out.items()
                    if not (basis[0] == 1 and destination[0] == 0)}

        with patch.object(clock, "column", no_reverse), self.assertRaisesRegex(ValueError, "missing Hermitian column"):
            clock.prepared_span(rail, config, {1: 1.0})

    def test_entire_interval_certificate_uses_exact_moments_and_boundary_tail(self):
        rail = clock.Rail()
        # Certificate construction must not depend on sampled-time propagation.
        with patch.object(clock, "evolve_clock", side_effect=RuntimeError("time-grid substitution")):
            result = clock.interval_certificate(rail, 256, 512)
        mean, variance, velocity, velocity_variance, covariance = _exact_moments(64)
        self.assertEqual(covariance, 0)
        self.assertEqual(result["window"], [256, 512])
        self.assertEqual(result["completion_position"], 2116)
        self.assertEqual(_fraction(result["cantelli_failure"]), Fraction(226536448, 278318778205))
        self.assertEqual(_fraction(result["position_variance"]), variance)
        self.assertEqual(_fraction(result["mean_velocity"]), velocity)
        self.assertEqual(_fraction(result["velocity_variance"]), velocity_variance)
        # Algebraic derivative numerator of V(t)/(V(t)+(v*t-c)^2).
        # For all t>=256, v*t-c>0 and v*A+B*c*t>0, so it decreases.
        c = Fraction(64) - mean + 3
        self.assertGreater(velocity * 256 - c, 0)
        self.assertGreater(velocity * variance, 0)
        self.assertGreater(velocity_variance * c, 0)
        for t in (Fraction(256), Fraction(513, 2), Fraction(512)):
            v = variance + velocity_variance * t * t
            separation = velocity * t - c
            derivative_numerator = 2 * velocity_variance * t * separation ** 2 - 2 * velocity * separation * v
            self.assertEqual(derivative_numerator,
                             -2 * separation * (velocity * variance + velocity_variance * c * t))
        distance = min(rail.offset + 1, rail.sites - rail.offset - rail.order)
        boundary = 2 * Fraction(3 * 512, distance) ** distance / (1 - Fraction(512, distance + 1))
        self.assertEqual(result["boundary_distance"], distance)
        self.assertEqual(_fraction(result["boundary_ratio"]), Fraction(24, 31))
        self.assertEqual(result["boundary_power"], distance)
        self.assertEqual(_fraction(result["boundary_prefactor"]), Fraction(3970, 1473))
        self.assertLessEqual(boundary, _fraction(result["boundary_upper"]))
        self.assertLessEqual(_fraction(result["boundary_upper"]), Fraction(1, 2 ** 659))
        self.assertEqual(_fraction(result["failure_upper"]),
                         Fraction(226536448, 278318778205) + _fraction(result["boundary_upper"]))
        self.assertLess(_fraction(result["failure_upper"]), Fraction(1, 1000))
        with self.assertRaises(ValueError):
            clock.interval_certificate(rail, 1, 512)
        with self.assertRaises(ValueError):
            clock.interval_certificate(rail, 256, 1000)

    def test_packet_representation_bound_covers_actual_binary_amplitudes(self):
        for order in (1, 4, 8, 64):
            rail = clock.Rail(256, 32, order)
            values = clock.packet(rail)[rail.offset:rail.start + 1]
            z = comb(2 * order, order)
            represented_root = Fraction.from_float(sqrt(z))
            actual = [Fraction.from_float(float(abs(value.real or value.imag))) for value in values]
            division = sum((abs(value - Fraction(comb(order, j), 1) / represented_root)
                            for j, value in enumerate(actual)), Fraction())
            expected = abs(represented_root ** 2 - z) / represented_root ** 2 + division
            self.assertEqual(clock.preparation_error_bound(rail), expected)
            with localcontext() as ctx:
                ctx.prec = 110
                root = Decimal(z).sqrt()
                norm_squared = sum((Decimal(value.numerator) / Decimal(value.denominator)
                                    - Decimal(comb(order, j)) / root) ** 2
                                   for j, value in enumerate(actual))
                decimal_bound = Decimal(expected.numerator) / Decimal(expected.denominator)
                self.assertLessEqual(norm_squared.sqrt(), decimal_bound)
            self.assertLess(float(expected), 1e-13)

    def test_direction_and_late_return_are_separate_exploratory_controls(self):
        rail = clock.Rail()
        for tau in (256, 512):
            forward = clock.evolve_clock(rail, tau)
            backward = clock.evolve_clock(rail, tau, -1)
            self.assertGreater(float(np.sum(abs(forward[rail.completed:]) ** 2)), 0.999)
            self.assertLess(float(np.sum(abs(backward[rail.completed:]) ** 2)), 1e-12)
        late = clock.evolve_clock(rail, 5120)
        self.assertLess(float(np.sum(abs(late[rail.completed:]) ** 2)), 0.001)
        self.assertGreater(float(np.sum(abs(late[:rail.completed]) ** 2)), 0.999)

    def test_pointer_dephasing_energy_for_preparation_and_complex_clock_control(self):
        rail, config = clock.Rail(48, 8, 4), Program()
        sources = {0: 1 / sqrt(2), 1: 1j / sqrt(2)}
        states, h, initial = clock.prepared_span(rail, config, sources)
        for tau in (0, 7.25, 19.5):
            vector = clock.literal_evolve(h, initial, tau)
            energy = float(np.vdot(vector, h @ vector).real)
            self.assertAlmostEqual(energy, 0, places=13)
            self.assertAlmostEqual(_pointer_energy(states, h, vector), energy, places=13)
        # A different coherent clock preparation has nonzero edge energy.
        # Reading B destroys that edge coherence, so zero work is not universal.
        history = prefixes(config, preparation(1, 1))
        before, after = (rail.start + 3, history[3]), (rail.completed, history[4])
        self.assertNotEqual(before[1].pointer_b, after[1].pointer_b)
        index = {state: i for i, state in enumerate(states)}
        vector = np.zeros(len(states), complex)
        vector[index[before]] = 1 / sqrt(2)
        vector[index[after]] = (1 + 1j) / 2
        energy = float(np.vdot(vector, h @ vector).real)
        self.assertAlmostEqual(energy, 1 / (2 * sqrt(2)), places=14)
        self.assertAlmostEqual(_pointer_energy(states, h, vector), 0, places=14)


if __name__ == "__main__":
    unittest.main()
