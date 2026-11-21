from fractions import Fraction

import numpy as np
import pytest

from model.executable_model import cell, federation

from .clock import Registers, closing_tick, interval_response, material_replay


def reference_suppression(q):
    """Independent killed-walk propagation, counting all probability mass at first return."""
    if q == 1:
        return Fraction(1)
    choices = 1 << (q - 1).bit_length()
    success = (1 - Fraction(choices - q, choices)**12) / q
    states = {-1: Fraction(1, 2), 1: Fraction(1, 2)}
    suppressed = Fraction(0)
    for _ in range(11):
        after = {}
        for x, mass in states.items():
            for y in (x - 1, x + 1):
                step = success * mass / 2
                if y == 0:
                    suppressed += step
                else:
                    after[y] = after.get(y, 0) + step
        states = after
    return suppressed


@pytest.mark.parametrize('bits', (3, 4))
def test_complete_interval_range_and_boundary(bits):
    for q in range(1, 2**bits):
        row = interval_response(q, age_bits=bits)
        assert row['qualified']
        assert row['closing_candidate_suppression'] == pytest.approx(float(reference_suppression(q)), abs=1e-15)
    overflow = interval_response(2**bits, age_bits=bits)
    assert overflow['finite_rate'] == 1
    assert not overflow['qualified']


def test_clock_reads_only_bounded_draws():
    def never():
        pytest.fail('deterministic branch drew a bit')

    assert not closing_tick(1, 1, never)
    assert closing_tick(8, 1, never)
    assert not closing_tick(2, 0, never)
    assert not closing_tick(2, 1, iter([1, 1, 0]).__next__)
    # At q=3, eleven is an invalid two-bit address. Twelve such rejections release the primary.
    assert closing_tick(3, 1, iter([1] * 25).__next__)
    # First signed step +1, q=8 closing address 111, second signed step -1.
    assert not closing_tick(8, 1, iter([1, 1, 1, 1, 0]).__next__, age_bits=4)
    assert closing_tick(16, 1, never, age_bits=4)


@pytest.mark.parametrize('level,bits', [(0, 3), (0, 4), (1, 4), (2, 3)])
def test_native_material_replay_and_overlap(level, bits):
    schedule = (8, 1, 16, 2)
    row = material_replay(federation(level), schedule, age_bits=bits, seed=31)
    events = row['local_events']
    handoffs = [k + 1 for k, e in enumerate(events) if e['handoff']]
    assert handoffs == np.cumsum(schedule).tolist()
    assert len(events) == sum(schedule)
    assert row['reunion']['home'] == sum(e['home_primary'] for e in events)
    assert row['reunion']['traveller'] == row['home_retained_traveller'] == sum(e['tick'] for e in events)
    assert row['moving_cells'] == (12 if bits == 3 else 13)
    assert row['transfer_cells'] == row['moving_cells'] + 1
    assert row['gates']['intercell_pi_contacts'] > 0
    assert row['max_native_column_error'] < 1e-11
    assert events[7]['primary'] == 1
    if bits == 3:
        assert events[7]['tick'] == 1 and events[7]['extra_bits'] == []
    else:
        assert events[7]['tick'] == 0 and events[7]['extra_bits'] == [1, 1, 1, 1, 0]
    for i, event in enumerate(events):
        if event['handoff']:
            assert event['age'] == 0
        else:
            assert event['tick'] == event['primary']
            assert event['age'] == min(events[i - 1]['age'] + 1 if i else 1, 2**bits - 1)
    assert events[8]['tick'] == 0  # The unit handoff resets after either overflow or a valid q=8.


def test_native_custody_and_guards():
    car = federation(0)
    reg = Registers(car, 31)
    a, b = car.external[0]
    owner = (a // 12, b // 12)
    with pytest.raises(ValueError, match='outside'):
        reg.swap(a, b, owner[:1])
    with pytest.raises(ValueError, match='outside'):
        reg.flag(0, 2, (0,))
    before = np.random.default_rng(101).choice([0, 2], len(reg.loads))
    reg.loads[:] = before
    reg.shift(owner, True)
    for src, dst in (owner, owner[::-1]):
        assert np.array_equal(reg.loads[12*dst:12*dst+12], before[12*src:12*src+12])
    reg.shift(owner, False)
    assert np.array_equal(reg.loads, before)


def test_explicit_home_reunion_required():
    with pytest.raises(ValueError, match='pairs of positive'):
        material_replay(federation(0), [3])
    with pytest.raises(ValueError, match='at least'):
        material_replay(cell(), [2, 3])


def test_two_clock_inputs_obey_the_same_physical_publication_even_when_inverted():
    car = federation(0)
    reg = Registers(car, seed=51)
    c, d = (int(p) // 12 for p in car.external[0])
    private = next(q // 12 for a, b in car.external for p, q in ((a, b), (b, a))
                   if p // 12 == d and q // 12 != c)
    owner = {c, d, private}
    public = (12*c, 12*c + 1, 12*c + 10)
    elapsed = [(12*private, 12*private + 1, 12*private + 10)]
    probe = reg.prepare_oscillator(c, owner)
    reg.prepare_front(d)
    reg.initialize_digits([public, *elapsed])
    reg.publish_tick(public, 1, owner)
    # A real rail exchange corrupts the committed one before the listeners are serviced.
    # Both listeners must read zero even though the producer admitted a one.
    reg.swap(*public[:2], owner)
    assert reg.post_tick(public, d, elapsed, owner) == probe.on_tick(reg, public, owner) == 0
    assert reg.reading(d, elapsed, owner) == 0
    assert np.array_equal(probe.state, [1, 0])
    reg.publish_tick(public, 1, owner)
    assert reg.post_tick(public, d, elapsed, owner) == probe.on_tick(reg, public, owner) == 1
    assert reg.reading(d, elapsed, owner) == 1
    assert abs(probe.state[1]) > .09
    before = probe.state.copy()
    reg.publish_tick(public, 0, owner)
    assert reg.post_tick(public, d, elapsed, owner) == probe.on_tick(reg, public, owner) == 0
    assert np.array_equal(probe.state, before)
    assert reg.reading(d, elapsed, owner) == 1
    assert reg.read_tick(public, {c}) == 0
    with pytest.raises(ValueError, match='outside'):
        reg.read_tick(public, {d})
    reg.loads[public[0]], reg.loads[public[1]] = 1, 1
    with pytest.raises(ValueError, match='orthogonal'):
        reg.publish_tick(public, 1, owner)


def test_public_rail_survives_transfer_and_rejects_a_positive_after_overflow():
    car = federation(0)
    rows = [material_replay(car, (8, 1, 16, 2), age_bits=w, seed=31) for w in (3, 4)]
    for row in rows:
        events = row['local_events']
        assert [e['oscillator_tick'] for e in events] == [e['tick'] for e in events]
        assert [e['home_tick'] for e in events] == [e['home_oscillator_tick'] for e in events]
        assert [e['home_tick'] for e in events] == [e['home_primary'] for e in events]
        assert events[8]['handoff'] and events[8]['age'] == events[8]['tick'] == 0
        assert row['oscillators']['traveller']['owners'] == [row['oscillators']['traveller']['bit']]*2
    assert rows[0]['local_events'][7]['age'] == 0
    assert rows[0]['local_events'][7]['oscillator_tick'] == 1
    assert rows[1]['local_events'][7]['oscillator_tick'] == 0


@pytest.mark.parametrize('width', [3, 4, 5])
def test_saturated_age_releases_one_primary_then_fresh_unit_handoff_clears_it(width):
    from .clock import _fields, _digits, _material_update

    car = federation(0)
    bank = tuple(range(width + 8))
    owner = set(bank)
    reg = Registers(car, 0)
    _, fields = _fields(bank, width)
    reg.initialize_digits([digit for c in bank for digit in _digits(c)])
    digit = fields['accepted'][0]
    reg.put(fields['age'], (1 << width) - 1, owner)
    # Both flags are measured on local orthogonal digits. No extra fair draw is legal
    # on either the saturated branch or the immediately following one-service close.
    for name in ('closing', 'address', 'depth', 'attempt'):
        reg.put(fields[name], 1, owner)
    primary_digit = _digits(bank[0])[0]
    reg.assign(primary_digit, 1, owner)
    for expected in (1, 0):
        primary = reg.read(primary_digit, owner)
        result = _material_update(reg, owner, fields, True, primary, width,
                                  lambda: pytest.fail('saturated or unit closing borrowed a fair draw'))
        reg.publish_tick(digit, result, owner)
        assert reg.read_tick(digit, owner) == expected
        assert reg.get(fields['age'], owner) == 0
        for name in ('closing', 'address', 'depth', 'attempt'):
            assert reg.get(fields[name], owner) == 1
