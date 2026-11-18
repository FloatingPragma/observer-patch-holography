"""Exact gate, support and source tests; no confirmatory receipt is written."""

from dataclasses import replace
from itertools import combinations, product
import random
import unittest

from .program import (Data, Program, REGISTER_CELLS, edge, gate_column,
                      gate_support, prefixes, preparation, public_record,
                      resources, validate_data)


_EDGES = tuple(combinations(range(4), 2))
_ADDRESSES = tuple(pair for pair in _EDGES for pair in (pair, pair[::-1]))


def _event(source, address, background, law):
    """Incidence/flux and strict crossing, without the production event code."""
    q = list(background)
    for e, (left, right) in enumerate(_EDGES):
        flux = ((source >> (2 * e + 1)) & 1) - ((source >> (2 * e)) & 1)
        q[left] += flux
        q[right] -= flux
    tail, head = _ADDRESSES[address]
    if law == "uphill":
        tail, head = head, tail
    bit = (source >> address) & 1
    current = 1 if not bit and q[tail] > q[head] else -1 if bit and q[head] > q[tail] else 0
    return source ^ ((1 << address) if current else 0), bit


def _registers(data):
    out = {f"S{k}": value for k, value in enumerate(data.buffers)}
    out.update({f"A{k}": value for k, value in enumerate(data.addresses)})
    out.update({f"C{k + 1}": value for k, value in enumerate(data.couriers)})
    out.update(PA=data.pointer_a, PB=data.pointer_b, b1=data.old_bits[0],
               b2=data.old_bits[1], D=data.second_address)
    return out


def _vary_complement(data, support):
    return Data(
        tuple(v if f"S{k}" in support else v ^ 4095 for k, v in enumerate(data.buffers)),
        tuple(v if f"A{k}" in support else (v + 5) % 12 for k, v in enumerate(data.addresses)),
        tuple(v if f"C{k + 1}" in support else v ^ 3 for k, v in enumerate(data.couriers)),
        data.pointer_a if "PA" in support else data.pointer_a ^ 3,
        data.pointer_b if "PB" in support else data.pointer_b ^ 3,
        tuple(v if f"b{k + 1}" in support else v ^ 1 for k, v in enumerate(data.old_bits)),
        data.second_address if "D" in support else (data.second_address + 5) % 12,
    )


def _random_data(rng):
    return Data(tuple(rng.randrange(4096) for _ in range(5)),
                tuple(rng.randrange(12) for _ in range(4)),
                tuple(rng.randrange(4) for _ in range(3)),
                rng.randrange(4), rng.randrange(4),
                tuple(rng.randrange(2) for _ in range(2)), rng.randrange(12))


class ProgramTests(unittest.TestCase):
    def test_first_event_matches_primary_rule_on_every_source_and_address(self):
        checked = 0
        for background, law in product(((0, 0, 0, 0), (1, 0, 0, 0)), ("native", "uphill")):
            config = Program(background=background, law=law)
            for source in range(4096):
                for address in range(12):
                    expected, old = _event(source, address, background, law)
                    output = edge(config, preparation(source, address), 0)
                    self.assertEqual(output.buffers, (0, expected, 0, 0, 0))
                    self.assertEqual(output.addresses, (0, address, 0, 0))
                    self.assertEqual(output.old_bits, (old, 0))
                    self.assertEqual((output.pointer_a, output.pointer_b, output.couriers), (0, 0, (0, 0, 0)))
                    checked += 1
        self.assertEqual(checked, 196608)

    def test_pointer_writer_and_sender_complete_small_domain(self):
        # All selected-bit, address, pointer and two courier combinations.
        for enabled, communication in product((False, True), repeat=2):
            config = Program(write_a=enabled, communication=communication)
            for bit, address, pa, c1, c2 in product(range(2), range(12), range(4), range(4), range(4)):
                data = Data((91, bit << address, 87, 23, 41), (8, address, 7, 6),
                            (c1, c2, 3), pa, 1, (1, 1), 9)
                output = edge(config, data, 1)
                expected_pa = pa ^ (2 | bit) if enabled else pa
                sent = c1 ^ expected_pa if communication else c1
                self.assertEqual(output.pointer_a, expected_pa)
                self.assertEqual(output.couriers, (c2, sent, 3))
                self.assertEqual(output.buffers, (91, 87, bit << address, 23, 41))
                self.assertEqual(output.addresses, (8, 7, address, 6))
                self.assertEqual(edge(config, output, 1, -1), data)

    def test_receiver_feedback_complete_pointer_address_domain(self):
        for law, background, feedback in product(("native", "uphill"),
                                                  ((0, 0, 0, 0), (1, 0, 0, 0)),
                                                  (False, True)):
            config = Program(background=background, law=law, feedback=feedback)
            for source, pb, courier, address, old in product((0, 1, 17, 1365, 2730, 4095),
                                                            range(4), range(4), range(12), range(2)):
                data = Data((9, 13, 11, source, 123), (1, 2, 3, 4), (1, 2, courier),
                            2, pb, (1, old), address)
                received = pb ^ courier
                actual_address = (address + 6 * (received & 1)) % 12 if feedback and received & 2 else address
                target, previous = _event(source, actual_address, background, law)
                target ^= old << actual_address
                output = edge(config, data, 3)
                self.assertEqual(output.pointer_b, received)
                self.assertEqual(output.second_address, actual_address)
                self.assertEqual(output.buffers, (9, 13, 11, 123, target))
                self.assertEqual(output.old_bits, (1, previous))
                self.assertEqual(edge(config, output, 3, -1), data)

    def test_dirty_inverses_columns_and_local_support_all_controls(self):
        rng = random.Random(1204)
        for law, background, communication, feedback, write_a in product(
            ("native", "uphill"), ((0, 0, 0, 0), (1, 0, 0, 0)),
            (False, True), (False, True), (False, True),
        ):
            config = Program(background, law, communication, feedback, write_a)
            for _ in range(12):
                data = _random_data(rng)
                trajectory = prefixes(config, data)
                recovered = trajectory[-1]
                for j in reversed(range(4)):
                    recovered = edge(config, recovered, j, -1)
                self.assertEqual(recovered, data)
                for j, direction in product(range(4), (1, -1)):
                    output = edge(config, data, j, direction)
                    self.assertEqual(edge(config, output, j, -direction), data)
                    self.assertEqual(gate_column(config, data, j, direction), {output: 1})
                    support = set(gate_support(j))
                    self.assertEqual({REGISTER_CELLS[name] for name in support}, {j, j + 1})
                    altered = _vary_complement(data, support)
                    varied = edge(config, altered, j, direction)
                    before, after, alternate = _registers(data), _registers(output), _registers(varied)
                    for name in before:
                        if name in support:
                            self.assertEqual(alternate[name], after[name])
                        else:
                            self.assertEqual(after[name], before[name])

    def test_actual_received_memory_controls_source_with_fixed_source_input(self):
        base = Data((0, 0, 0, 1, 0), (0, 0, 0, 0), (0, 0, 0), 0, 0, (0, 0), 0)
        outcomes = [edge(Program(), replace(base, couriers=(0, 0, value)), 3)
                    for value in range(4)]
        self.assertEqual([data.second_address for data in outcomes], [0, 0, 0, 6])
        self.assertEqual([data.buffers[4] for data in outcomes], [0, 0, 0, 65])
        disabled = edge(Program(feedback=False), replace(base, couriers=(0, 0, 3)), 3)
        self.assertEqual((disabled.pointer_b, disabled.second_address, disabled.buffers[4]), (3, 0, 0))

    def test_fresh_records_are_actual_register_values_and_controls_retain_invalid(self):
        for source, address in product((0, 1, 17, 1365, 2730, 4095), range(12)):
            trajectory = prefixes(Program(), preparation(source, address))
            after_first = trajectory[1].buffers[1]
            observed = (after_first >> address) & 1
            final = trajectory[-1]
            self.assertEqual((final.pointer_a, final.pointer_b, final.couriers),
                             (2 | observed, 2 | observed, (0, 0, 2 | observed)))
            self.assertEqual(final.second_address, 6 * observed)
            for config in (Program(communication=False), Program(write_a=False)):
                blocked = prefixes(config, preparation(source, address))[-1]
                self.assertEqual((blocked.pointer_b, blocked.second_address), (0, 0))
            self.assertEqual(public_record(trajectory[0])[:4], (0, 0, 0, 0))
        dirty = replace(preparation(1, 0), pointer_a=1, pointer_b=3)
        self.assertEqual(public_record(dirty)[:4], (0, 1, 1, 1))
        self.assertNotEqual(public_record(dirty), public_record(replace(dirty, buffers=(1, 0, 0, 0, 1))))
        self.assertEqual(int(resources()["full_data_dimension"]), 4096 ** 5 * 12 ** 5 * 4 ** 5 * 2 ** 2)

    def test_invalid_inputs_fail_under_optimized_python(self):
        for source in (-1, 4096, True, 1.0, "1", None):
            with self.assertRaises(ValueError):
                preparation(source, 0)
        for address in (-1, 12, True, 1.0, "1", None):
            with self.assertRaises(ValueError):
                preparation(1, address)
        for arguments in ({"background": [0, 0, 0, 0]}, {"background": (False, 0, 0, 0)},
                          {"background": (2, 0, 0, 0)}, {"law": "hold"}, {"law": None},
                          {"communication": 1}, {"feedback": 0}, {"write_a": None}):
            with self.assertRaises(ValueError):
                Program(**arguments)
        data = preparation(1, 0)
        for bad in (replace(data, buffers=[1, 0, 0, 0, 0]), replace(data, buffers=(4096,) * 5),
                    replace(data, addresses=(True, 0, 0, 0)), replace(data, couriers=(4, 0, 0)),
                    replace(data, old_bits=(0, 2)), replace(data, pointer_a=True),
                    replace(data, pointer_b=4), replace(data, second_address=12)):
            with self.assertRaises(ValueError):
                validate_data(bad)
            with self.assertRaises(ValueError):
                edge(Program(), bad, 0)
        for j in (-1, 4, True, 1.0):
            with self.assertRaises(ValueError):
                edge(Program(), data, j)
        for direction in (0, 2, True, 1.0):
            with self.assertRaises(ValueError):
                edge(Program(), data, 0, direction)


if __name__ == "__main__":
    unittest.main()
