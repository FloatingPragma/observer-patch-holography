"""Composition attacks for the autonomous two-observer program.

These are exploratory implementation tests, not confirmatory receipts. The
source-85 discrimination witness was chosen by native-event enumeration.
The independent oracle calls only the inherited native event component.
"""

from collections import Counter, defaultdict
from dataclasses import replace
from fractions import Fraction
from itertools import product
from math import sqrt
import unittest

import numpy as np

from model.native_memory_clock.native import BACKGROUNDS, inverse, native_event

from . import oracle, program


SOURCES = (0, 1, 2, 3, 17, 85, 1365, 2730)


def _literal_columns(config, sources, stage):
    """Read the executed register outputs, preserving their source factor."""
    columns = {}
    for source in sources:
        for address in range(12):
            data = program.preparation(source, address)
            for j in range(stage):
                data = program.edge(config, data, j)
            # Deliberately do not call the oracle's archive-label function.
            label = (data.pointer_a, data.pointer_b, data.second_address,
                     data.addresses, data.couriers, data.old_bits)
            key = (data.buffers, source)
            block = columns.setdefault(label, {})
            block[key] = block.get(key, 0.0) + 1 / sqrt(12)
    return columns


def _matrices(columns, sources):
    source_index = {source: index for index, source in enumerate(sources)}
    result = {}
    for outcome, entries in columns.items():
        outputs = sorted({buffers for buffers, _ in entries})
        output_index = {buffers: index for index, buffers in enumerate(outputs)}
        matrix = np.zeros((len(outputs), len(sources)), complex)
        for (buffers, source), value in entries.items():
            matrix[output_index[buffers], source_index[source]] += value
        result[outcome] = (tuple(outputs), matrix)
    return result


def _source_density(amplitudes, cell):
    """Partial trace of every other actual register, without projection."""
    groups = defaultdict(dict)
    for data, amplitude in amplitudes.items():
        source = data.buffers[cell]
        buffers = list(data.buffers)
        buffers[cell] = 0
        environment = replace(data, buffers=tuple(buffers))
        group = groups[environment]
        group[source] = group.get(source, 0j) + amplitude
    density = defaultdict(complex)
    for vector in groups.values():
        for left, a in vector.items():
            for right, b in vector.items():
                density[left, right] += a * b.conjugate()
    return density


class ChannelTests(unittest.TestCase):
    def test_every_prefix_matches_independent_all_register_semantics(self):
        for background, law, communication, feedback, write_a in product(
                BACKGROUNDS, ("native", "uphill"), (False, True),
                (False, True), (False, True)):
            config = program.Program(background, law, communication, feedback, write_a)
            for source, address in product(SOURCES, range(12)):
                expected = oracle.expected_prefixes(config, source, address)
                actual = program.prefixes(config, program.preparation(source, address))
                self.assertEqual(actual, expected)
                for output in actual:
                    self.assertEqual(program.public_record(output), oracle.record(output))

    def test_complete_kraus_maps_keep_reference_and_incomplete_observer_branches(self):
        # Nonorthogonal complex reference vectors expose coherences within
        # each archive label, rather than only classical source frequencies.
        coefficients = np.array([[1 + i / 9, (i % 3 - 1) * 1j, (i + 1j) / 7]
                                 for i in range(len(SOURCES))], complex)
        coefficients /= np.linalg.norm(coefficients)
        input_reference = coefficients.T @ coefficients.conjugate()
        retained_coherence = 0.0
        for background, law, communication, feedback, write_a in product(
                BACKGROUNDS, ("native", "uphill"), (False, True),
                (False, True), (False, True)):
            config = program.Program(background, law, communication, feedback, write_a)
            for stage in range(5):
                actual = _matrices(_literal_columns(config, SOURCES, stage), SOURCES)
                expected = _matrices(oracle.branch_columns(config, SOURCES, stage), SOURCES)
                self.assertEqual(actual.keys(), expected.keys())
                gram = np.zeros((len(SOURCES), len(SOURCES)), complex)
                reference = np.zeros((3, 3), complex)
                total = 0.0
                for label, (outputs, matrix) in actual.items():
                    expected_outputs, wanted = expected[label]
                    self.assertEqual(outputs, expected_outputs)
                    np.testing.assert_array_equal(matrix, wanted)
                    gram += matrix.conjugate().T @ matrix
                    vector = matrix @ coefficients
                    expected_vector = wanted @ coefficients
                    density = np.outer(vector.ravel(), vector.ravel().conjugate())
                    expected_density = np.outer(expected_vector.ravel(), expected_vector.ravel().conjugate())
                    np.testing.assert_allclose(density, expected_density, atol=2e-15, rtol=0)
                    total += float(np.vdot(vector, vector).real)
                    reference += vector.T @ vector.conjugate()
                    if stage == 4 and len(outputs) > 1:
                        shaped = density.reshape(len(outputs), 3, len(outputs), 3)
                        retained_coherence += float(sum(np.abs(shaped[i, :, j, :]).sum()
                                                        for i in range(len(outputs))
                                                        for j in range(len(outputs)) if i != j))
                np.testing.assert_allclose(gram, np.eye(len(SOURCES)), atol=2e-15, rtol=0)
                np.testing.assert_allclose(reference, input_reference, atol=2e-15, rtol=0)
                self.assertAlmostEqual(total, 1.0, places=14)
                if stage < 4 or not communication or not write_a:
                    self.assertTrue(all(not (label[1] & 2) for label in actual))
        self.assertGreater(retained_coherence, 0.1)

    def test_selected_bit_feedback_and_same_apparatus_law_controls(self):
        source = 85  # Exploratory discriminator, not a held-out physical input.
        native = program.Program()
        variants = {"native": native, "uphill": replace(native, law="uphill"),
                    "open_loop": replace(native, feedback=False),
                    "no_message": replace(native, communication=False),
                    "no_write": replace(native, write_a=False)}
        laws = {}
        for name, config in variants.items():
            outcomes = [program.prefixes(config, program.preparation(source, a))[-1]
                        for a in range(12)]
            self.assertEqual(outcomes,
                             [oracle.expected_prefixes(config, source, a)[-1] for a in range(12)])
            laws[name] = Counter(data.buffers[4] for data in outcomes)
            self.assertEqual(sum(laws[name].values()), 12)
            if name == "native":
                self.assertEqual(Counter(data.pointer_a & 1 for data in outcomes), {0: 6, 1: 6})
                self.assertEqual(Counter(data.second_address for data in outcomes), {0: 6, 6: 6})
                self.assertTrue(all(data.pointer_a == data.pointer_b and data.pointer_a & 2
                                    for data in outcomes))
            elif name == "no_message":
                self.assertTrue(all(data.pointer_a & 2 and data.pointer_b == 0 for data in outcomes))
            elif name == "no_write":
                self.assertTrue(all(data.pointer_a == data.pointer_b == 0 for data in outcomes))
        def tv(left, right):
            return Fraction(sum(abs(left[x] - right[x]) for x in left.keys() | right.keys()), 24)
        self.assertEqual(tv(laws["native"], laws["open_loop"]), Fraction(1, 2))
        self.assertEqual(tv(laws["native"], laws["uphill"]), 1)
        self.assertEqual(laws["open_loop"], laws["no_message"])
        self.assertEqual(laws["open_loop"], laws["no_write"])

    def test_pointer_only_instrument_preserves_coherence_between_private_histories(self):
        coefficients = np.array([[1 + i / 5, (i + 1) * 1j / 7, (i % 3 - 1 + 1j) / 4]
                                 for i in range(len(SOURCES))], complex)
        coefficients /= np.linalg.norm(coefficients)
        input_reference = coefficients.T @ coefficients.conjugate()
        controls = ({}, {"law": "uphill"}, {"communication": False},
                    {"feedback": False}, {"write_a": False})
        cross_archive_coherence = 0.0
        dephasing_discrepancy = 0.0
        for background, options in product(BACKGROUNDS, controls):
            config = program.Program(background=background, **options)
            for stage in range(5):
                actual, expected = {}, {}
                for index, source in enumerate(SOURCES):
                    for address in range(12):
                        data = program.preparation(source, address)
                        for j in range(stage):
                            data = program.edge(config, data, j)
                        wanted = oracle.expected_prefixes(config, source, address)[stage]
                        for output, blocks in ((data, actual), (wanted, expected)):
                            pointer = (output.pointer_a, output.pointer_b)
                            block = blocks.setdefault(pointer, {})
                            block[output] = block.get(output, np.zeros(3, complex)) + coefficients[index] / sqrt(12)
                self.assertEqual(actual.keys(), expected.keys())
                reference = np.zeros((3, 3), complex)
                trace = 0.0
                for pointer, block in actual.items():
                    # Every private register, including distinct first-address
                    # histories, is retained inside this one public outcome.
                    outputs = sorted(block)
                    self.assertEqual(set(outputs), set(expected[pointer]))
                    vector = np.array([block[data] for data in outputs])
                    wanted = np.array([expected[pointer][data] for data in outputs])
                    density = np.outer(vector.ravel(), vector.ravel().conjugate())
                    target = np.outer(wanted.ravel(), wanted.ravel().conjugate())
                    np.testing.assert_array_equal(density, target)
                    trace += float(np.trace(density).real)
                    reference += vector.T @ vector.conjugate()
                    if stage == 4 and not options and background == BACKGROUNDS[0]:
                        labels = [oracle.archive_label(data) for data in outputs]
                        same_archive = np.array([[left == right for right in labels] for left in labels])
                        mask = np.repeat(np.repeat(same_archive, 3, axis=0), 3, axis=1)
                        dephased = density * mask
                        difference = density - dephased
                        cross_archive_coherence += float(np.abs(difference).sum())
                        dephasing_discrepancy += float(np.linalg.norm(difference))
                    if stage < 4 or not config.communication or not config.write_a:
                        self.assertFalse(pointer[1] & 2)
                self.assertAlmostEqual(trace, 1.0, places=14)
                np.testing.assert_allclose(reference, input_reference, atol=2e-15, rtol=0)
        self.assertGreater(cross_archive_coherence, 1.0)
        self.assertGreater(dephasing_discrepancy, 0.1)

    def test_feedback_follows_received_pointer_under_source_counterfactuals(self):
        config = program.Program()
        template = oracle.expected_prefixes(config, 85, 0)[3]
        for source, receiver, courier in product((1, 17, 85, 2730), range(4), range(4)):
            received = receiver ^ courier
            data = replace(template, buffers=(0, 0, 0, source, 0),
                           couriers=(0, 0, courier), pointer_b=receiver,
                           second_address=0, old_bits=(0, 0))
            actual = program.edge(config, data, 3)
            chosen = 6 if received == 3 else 0
            self.assertEqual(actual.pointer_b, received)
            self.assertEqual(actual.second_address, chosen)
            self.assertEqual(actual.buffers[4], native_event(source, chosen).target)
            self.assertEqual(program.edge(replace(config, feedback=False), data, 3).second_address, 0)
        # Identical source, distinct received record, distinct source action.
        endpoints = []
        for message in (2, 3):
            data = replace(template, buffers=(0, 0, 0, 85, 0), pointer_b=0,
                           couriers=(0, 0, message), old_bits=(0, 0), second_address=0)
            endpoints.append(program.edge(config, data, 3).buffers[4])
        self.assertEqual(endpoints, [84, 21])

    def test_retained_orthogonal_copies_block_native_only_coherence_recovery(self):
        config = program.Program()
        outputs = {program.prefixes(config, program.preparation(source, 0))[-1]: 1 / sqrt(2)
                   for source in (0, 3)}
        self.assertEqual({data.pointer_a for data in outputs}, {2, 3})
        self.assertTrue(all(data.pointer_a == data.pointer_b == data.couriers[2] for data in outputs))
        # Undo the two native events, preserving every observer/message copy.
        native_only = {}
        for data, amplitude in outputs.items():
            first, blank2 = inverse(data.buffers[4], data.second_address,
                                    data.old_bits[1], config.background, config.law)
            source, blank1 = inverse(first, data.addresses[3], data.old_bits[0],
                                    config.background, config.law)
            self.assertEqual((blank1, blank2), (0, 0))
            restored = replace(data, buffers=(0, 0, 0, 0, source), old_bits=(0, 0))
            native_only[restored] = amplitude
        reduced = _source_density(native_only, 4)
        self.assertAlmostEqual(reduced[0, 0].real, 0.5)
        self.assertAlmostEqual(reduced[3, 3].real, 0.5)
        self.assertEqual(reduced[0, 3], 0)
        # Reversibly uncopy B, the courier, and the controlled address while
        # retaining A. That one public copy alone obstructs source coherence.
        only_a = {}
        for data, amplitude in native_only.items():
            second_address = (data.second_address - 6 * (data.pointer_b & 1)) % 12
            cleaned = replace(data, second_address=second_address,
                              pointer_b=data.pointer_b ^ data.couriers[2],
                              couriers=(0, 0, data.couriers[2] ^ data.pointer_a))
            self.assertEqual(cleaned.pointer_b, 0)
            self.assertEqual(cleaned.couriers, (0, 0, 0))
            self.assertEqual(cleaned.second_address, 0)
            only_a[cleaned] = amplitude
        self.assertEqual(_source_density(only_a, 4)[0, 3], 0)
        self.assertEqual({data.pointer_a for data in only_a}, {2, 3})
        # The full inverse restores coherence by uncomputing every public copy.
        fully_reversed = {}
        for data, amplitude in outputs.items():
            for j in reversed(range(4)):
                data = program.edge(config, data, j, -1)
            self.assertEqual(data.pointer_a, 0)
            self.assertEqual(data.pointer_b, 0)
            self.assertEqual(data.couriers, (0, 0, 0))
            fully_reversed[data] = amplitude
        recovered = _source_density(fully_reversed, 0)
        self.assertAlmostEqual(recovered[0, 3].real, 0.5)
        self.assertAlmostEqual(sum(recovered[x, y].real for x, y in product((0, 3), repeat=2)) / 2, 1)


if __name__ == "__main__":
    unittest.main()
