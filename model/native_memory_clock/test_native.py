"""Exact finite unit tests; scientific receipts are produced by the package replay."""

from fractions import Fraction
import unittest

from .native import (
    ADDRESSES,
    ADDRESS_COUNT,
    BACKGROUNDS,
    EDGES,
    LAWS,
    STATE_COUNT,
    decode_history,
    dependency_addresses,
    forward,
    inverse,
    loads,
    native_event,
    no_old_bit_collision,
    omitted_record_overlap,
    reused_record_counterexample,
)


class _PrimaryCrossing:
    """Primary method body from the revision and path declared in native.py."""

    kind = "rule"

    def crossing(self, bit, tail, head):
        if self.kind == "ties only" and abs(tail - head) != 1:
            return 0
        return 1 if not bit and tail > head else -1 if bit and head > tail else 0


def _independent_loads(state, background):
    # Incidence-times-flux formulation, independent of address accumulation.
    result = list(background)
    for edge_index, (left, right) in enumerate(EDGES):
        flux = ((state >> (2 * edge_index + 1)) & 1) - ((state >> (2 * edge_index)) & 1)
        result[left] += flux
        result[right] -= flux
    return tuple(result)


class NativeTests(unittest.TestCase):
    def test_primary_crossing_and_native_ledger_exhaustive(self):
        primary = _PrimaryCrossing()
        checked = 0
        for background in BACKGROUNDS:
            for state in range(STATE_COUNT):
                before = _independent_loads(state, background)
                self.assertEqual(loads(state, background), before)
                for address, (tail, head) in enumerate(ADDRESSES):
                    bit = (state >> address) & 1
                    current = primary.crossing(bit, before[tail], before[head])
                    expected = state ^ (1 << address) if current else state
                    event = native_event(state, address, background)
                    self.assertEqual(event.target, expected)
                    self.assertEqual(event.current, current)
                    self.assertEqual(event.new_bit, (expected >> address) & 1)
                    self.assertEqual(event.after_loads, _independent_loads(expected, background))
                    self.assertEqual(sum(event.before_loads), sum(event.after_loads))
                    self.assertGreaterEqual(event.energy_fall, 0)
                    self.assertEqual(
                        event.energy_fall,
                        2 * current * (before[tail] - before[head]) - 2 * current * current,
                    )
                    self.assertEqual(forward(state, address, 0, background), (expected, bit))
                    seam, coin = divmod(address, 2)
                    decoded = decode_history(expected, (address,), (bit,), background)
                    self.assertEqual(decoded, (state, ((seam, coin, -current if coin else current),)))
                    checked += 1
        self.assertEqual(checked, 98_304)

    def test_terminal_decoder_all_two_address_pairs_and_controls(self):
        primary = _PrimaryCrossing()
        sources = (0, 1, 2, 3, 17, 0x555, 0xAAA, 0xFFF)
        for law in LAWS:
            for background in BACKGROUNDS:
                for source in sources:
                    for first in range(ADDRESS_COUNT):
                        for second in range(ADDRESS_COUNT):
                            addresses = (first, second)
                            endpoint, old_bits, expected = source, [], []
                            for address in addresses:
                                before = _independent_loads(endpoint, background)
                                tail, head = ADDRESSES[address]
                                bit = (endpoint >> address) & 1
                                old_bits.append(bit)
                                if law == "hold":
                                    current = 0
                                elif law == "uphill":
                                    current = primary.crossing(bit, before[head], before[tail])
                                else:
                                    current = primary.crossing(bit, before[tail], before[head])
                                if current:
                                    endpoint ^= 1 << address
                                seam, coin = divmod(address, 2)
                                expected.append((seam, coin, -current if coin else current))
                            decoded = decode_history(endpoint, addresses, tuple(old_bits), background, law)
                            self.assertEqual(decoded, (source, tuple(expected)))
        self.assertEqual(decode_history(17, (), ()), (17, ()))

    def test_terminal_decoder_rejects_impossible_or_malformed_histories(self):
        for endpoint, addresses, old_bits in ((1, (0,), (0,)), (1, (0,), (1,)),
                                              (0, (0, 0), (0, 1))):
            with self.subTest(history=(endpoint, addresses, old_bits)), self.assertRaises(ValueError):
                decode_history(endpoint, addresses, old_bits)
        # These same terminal labels are possible under a different declared source.
        self.assertEqual(decode_history(1, (0,), (1,), law="hold"), (1, ((0, 0, 0),)))
        self.assertEqual(decode_history(1, (0,), (0,), BACKGROUNDS[1]), (0, ((0, 0, 1),)))
        for arguments in ((True, (), ()), (STATE_COUNT, (), ()), (0, [], ()),
                          (0, (), []), (0, (0,), ()), (0, (), (0,)),
                          (0, (-1,), (0,)), (0, (12,), (0,)), (0, (False,), (0,)),
                          (0, (0,), (True,)), (0, (0,), (2,)), (0, (0,), (0.0,))):
            with self.subTest(arguments=arguments), self.assertRaises(ValueError):
                decode_history(*arguments)
        with self.assertRaises(ValueError):
            decode_history(0, (), (), law="unknown")
        with self.assertRaises(ValueError):
            decode_history(0, (), (), (0.0, 0, 0, 0))

    def test_total_permutation_including_dirty_records_exhaustive(self):
        # Each control uses the same complete source/address/record domain.
        for law in LAWS:
            checked = 0
            for background in BACKGROUNDS:
                for address in range(ADDRESS_COUNT):
                    outputs = set()
                    for state in range(STATE_COUNT):
                        for record in (0, 1):
                            output = forward(state, address, record, background, law)
                            self.assertEqual(inverse(output[0], address, output[1], background, law), (state, record))
                            back = inverse(state, address, record, background, law)
                            self.assertEqual(forward(back[0], address, back[1], background, law), (state, record))
                            outputs.add(output)
                            checked += 1
                    self.assertEqual(len(outputs), STATE_COUNT * 2)
            self.assertEqual(checked, 196_608)

    def test_endpoint_star_support_exhaustive(self):
        for address in range(ADDRESS_COUNT):
            support = dependency_addresses(address)
            outside = tuple(index for index in range(ADDRESS_COUNT) if index not in support)
            self.assertEqual(len(support), 10)
            self.assertEqual(len(outside), 2)
            # Every support configuration and all complement assignments.
            for background in BACKGROUNDS:
                for encoded in range(1 << len(support)):
                    source = sum(((encoded >> i) & 1) << bit for i, bit in enumerate(support))
                    for law in LAWS:
                        expected = native_event(source, address, background, law).new_bit
                        for complement in range(1 << len(outside)):
                            state = source | sum(
                                ((complement >> i) & 1) << bit for i, bit in enumerate(outside)
                            )
                            self.assertEqual(native_event(state, address, background, law).new_bit, expected)

    def test_explicit_information_loss_and_reuse_witnesses(self):
        witness = no_old_bit_collision()
        self.assertNotEqual(witness["left"], witness["right"])
        self.assertEqual(witness["overlap"], Fraction(1, 12))
        self.assertEqual(omitted_record_overlap(0, 0), 1)
        self.assertEqual(forward(0, 0, 0), (0, 0))
        self.assertEqual(forward(1, 0, 0), (0, 1))
        reuse = reused_record_counterexample()
        self.assertEqual(reuse["first"], (0, 1))
        self.assertEqual(reuse["reused"], (1, 0))
        self.assertEqual(reuse["fresh"], (0, 0))
        self.assertNotEqual(reuse["reused"][0], reuse["fresh"][0])

    def test_controls_are_distinct_with_same_records(self):
        native = native_event(1, 0)
        hold = native_event(1, 0, law="hold")
        uphill = native_event(0, 0, BACKGROUNDS[1], "uphill")
        self.assertEqual((native.target, hold.target), (0, 1))
        self.assertEqual(uphill.target, 0)
        # A nonzero uphill event raises the same squared-load diagnostic.
        uphill = native_event(1, 2, law="uphill")
        self.assertNotEqual(uphill.current, 0)
        self.assertLess(uphill.energy_fall, 0)
        for state in (0, 1, 2, 10, STATE_COUNT - 1):
            for address in range(ADDRESS_COUNT):
                self.assertEqual(forward(state, address, 0, law="hold"), (state, (state >> address) & 1))

    def test_invalid_inputs_rejected_under_optimized_python(self):
        for value in (-1, STATE_COUNT, 0.0, True, "0", None):
            with self.subTest(state=value), self.assertRaises(ValueError):
                loads(value)
            with self.subTest(target=value), self.assertRaises(ValueError):
                inverse(value, 0, 0)
        for value in (-1, ADDRESS_COUNT, 0.0, True, "0", None):
            for operation in (native_event, dependency_addresses):
                arguments = (0, value) if operation is native_event else (value,)
                with self.subTest(address=value, operation=operation), self.assertRaises(ValueError):
                    operation(*arguments)
        for value in (-1, 2, 0.0, True, "0", None):
            for operation in (forward, inverse):
                with self.subTest(record=value, operation=operation), self.assertRaises(ValueError):
                    operation(0, 0, value)
        for background in ([0, 0, 0, 0], (0, 0, 0), (False, 0, 0, 0), (0.0, 0, 0, 0), (2, 0, 0, 0)):
            with self.subTest(background=background), self.assertRaises(ValueError):
                native_event(0, 0, background)
        for law in (None, 0, True, [], "unknown"):
            with self.subTest(law=law), self.assertRaises(ValueError):
                forward(0, 0, 0, law=law)


if __name__ == "__main__":
    unittest.main()
