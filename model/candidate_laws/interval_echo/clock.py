"""Interval Echo's bounded clock and an owned, sharp-register replay on the common federation.

The fixture selects when to attempt an adjacent handoff. Only a marker read on the departed and arrived
patches reaches the clock. Born draws select sharp outcomes of a checked native seam instrument.
"""
from __future__ import annotations

from collections import Counter, deque
from functools import lru_cache
from math import comb, sqrt
from operator import index

import numpy as np
from scipy.linalg import expm

from .candidate import KAPPA, native_law as law


STEPS = 12
ADDRESS_TRIES = 12
# (first rail, second rail, adjacent vacancy used to read the first rail)
DIGITS = ((0, 1, 10), (2, 3, 10), (5, 11, 4))
EXTRA_DIGITS = ((6, 7, 10), (9, 8, 4))
RING = (1, 5, 11, 10, 7)


def _capacity(age_bits):
    if age_bits not in (3, 4, 5):
        raise ValueError('physical clock layouts have three, four or five interval bits')
    return (1 << age_bits) - 1


def closing_tick(q, primary, fair, *, age_bits=3):
    """Fast scalar clock decision conditional on the measured length, primary mark and fair draws."""
    capacity = _capacity(age_bits)
    if q < 1:
        raise ValueError('a handoff must close a positive interval')
    if not primary or q > capacity:
        return bool(primary)
    if q == 1:
        return False
    position = 2 * fair() - 1
    width = (q - 1).bit_length()
    for _ in range(1, STEPS):
        for _ in range(ADDRESS_TRIES):
            address = 0
            for _ in range(width):
                address = 2 * address + fair()
            if address < q:
                break
        if address != q - 1:
            break
        position += 2 * fair() - 1
        if position == 0:
            return False
    return True


def interval_response(q, *, age_bits=3):
    """Exact finite-walk expectation up to floating-point rounding. This comparator is outside the clock."""
    capacity = _capacity(age_bits)
    if q < 1:
        raise ValueError('a handoff must close a positive interval')
    suppression = 0.0
    if q == 1:
        suppression = 1.0
    elif q <= capacity:
        width = (q - 1).bit_length()
        continuation = (1 - (1 - q / (1 << width)) ** ADDRESS_TRIES) / q
        suppression = sum(
            comb(2 * n, n) * continuation ** (2 * n - 1) / ((2 * n - 1) * 4**n)
            for n in range(1, STEPS // 2 + 1)
        )
    rate = 1 - suppression / q
    chord = sqrt(1 - 1 / q**2)
    return dict(interval=q, qualified=q <= capacity, closing_candidate_suppression=suppression,
                finite_rate=rate, resolved_chord_rate=chord, finite_minus_chord=rate - chord)


@lru_cache(maxsize=1)
def _primitives():
    """Check the reduced gates against all relevant columns of the 99-state seam, including pointer undo."""
    echo = law.echo()
    h0, h1 = law.local_term_parts()
    h = h0 + KAPPA * h1
    swap = expm(-1j * np.pi * h)
    half = expm(-0.5j * np.pi * h)
    pointer = np.arange(law.DIM) % len(law.SRC) == law.SRC.index((1, 1))
    error = 0.
    rates = {}
    for a in (0, 2):
        for b in (0, 2):
            ket = np.eye(law.DIM)[:, law.SRC.index((a, b))]
            swapped = np.eye(law.DIM)[:, law.SRC.index((b, a))]
            error = max(error, np.max(abs(swap @ ket - swapped)))
            written = echo @ ket
            # On sharp rails the native agreement is certain or impossible and is nondestructively reversed.
            error = max(error, np.linalg.norm(written[pointer != (a + b == 2)]),
                        np.linalg.norm(echo @ written - ket))
            if a + b == 2:
                output = half @ ket
                at_first, at_second = (law.SRC.index(x) for x in ((2, 0), (0, 2)))
                rates[a] = float(abs(output[at_second])**2)
                error = max(error, abs(abs(output[at_first])**2 - .5),
                            abs(rates[a] - .5),
                            np.linalg.norm(np.delete(output, [at_first, at_second])))
    if error > 1.e-11:
        raise ArithmeticError(f'native sharp-register reduction failed at {error}')
    return echo, rates, float(error)


def _guide(car, needed):
    neighbours = [[] for _ in range(car.cells)]
    for a, b in car.external:
        neighbours[a // 12].append(b // 12)
        neighbours[b // 12].append(a // 12)
    if needed > car.cells:
        raise ValueError(f'this clock fixture requires a federation with at least {needed} cells')

    def extend(path, used):
        if len(path) == needed:
            return path
        for nxt in sorted(neighbours[path[-1]]):
            if nxt not in used:
                found = extend(path + [nxt], used | {nxt})
                if found:
                    return found
        return None

    for first in range(car.cells):
        path = extend([first], {first})
        if path:
            return path
    raise ValueError('no native cell path can hold the enrolled registers')


def _digits(cell, layout=DIGITS):
    return [tuple(12 * cell + x for x in d) for d in layout]


@lru_cache(maxsize=2)
def _dimer_pulse(duration):
    """Restrict the actual 99-state native operator to its invariant sharp double rails."""
    h0, h1 = law.local_term_parts()
    generator = h0 + KAPPA * h1
    basis = np.eye(law.DIM)[:, [law.SRC.index((2, 0)), law.SRC.index((0, 2))]]
    reduced = basis.T @ generator @ basis
    if np.linalg.norm(generator @ basis - basis @ reduced) > 1e-11:
        raise ArithmeticError('the native generator does not preserve the oscillator rails')
    return expm(-1j * duration * reduced)


class SeamOscillator:
    """A single coherent double on two physical rails, kept in its exact reducing subspace."""

    def __init__(self, cell):
        self.rails = [12 * cell + 6, 12 * cell + 7]
        self.state = np.array([1., 0.], dtype=complex)

    def on_tick(self, registers, public_digit, owner):
        a, b = self.rails
        registers.native(a, b, owner)
        cell = a // 12
        if b // 12 != cell or public_digit[0] // 12 != cell:
            raise ValueError('the oscillator requires its local event register')
        tick = registers.read_tick(public_digit, owner)
        if tick:
            self.state = _dimer_pulse(np.pi / 16) @ self.state
            registers.gates['dimer_exposures'] += 1
        return tick

    def terminal(self, registers, owner, analyzer):
        """Test the coherent pointer at one selected setting and leave the two projected rails sharp."""
        if analyzer not in (0, 1):
            raise ValueError('select the in-phase (0) or quarter-period (1) terminal analyzer')
        a, b = self.rails
        registers.native(a, b, owner)
        rotated = _dimer_pulse(np.pi / 2) @ self.state
        probabilities = [float(abs(v[0])**2) for v in (self.state, rotated)]
        norm_error = abs(float(np.vdot(self.state, self.state).real) - 1.)
        if norm_error > 1e-10:
            raise ArithmeticError('oscillator norm lost during transport')
        outcome = int(registers.rng.random() < probabilities[analyzer])
        for port in self.rails:
            del registers.coherent[port]
        registers.loads[a], registers.loads[b] = 2 * outcome, 2 * (1 - outcome)
        bit = registers.flag(a, (a // 12) * 12 + 10, owner)
        if bit != outcome:
            raise ArithmeticError('terminal pointer disagrees with the projected coherent rails')
        return {'rails': list(self.rails), 'analyzer': analyzer, 'bit': int(bit),
                'first_rail_probabilities': [round(value, 12) for value in probabilities],
                'norm_error_below_1e-10': True}


class Registers:
    """Only adjacent contact permutations and certified pointer reads manipulate prepared sharp rails."""

    def __init__(self, car, seed):
        self.car = car
        self.loads = np.zeros(car.n, dtype=np.int8)
        self.rng = np.random.default_rng(seed)
        self.gates = Counter()
        self.routes = {}
        self.coherent = {}
        self.echo, self.rates, self.error = _primitives()

    def native(self, a, b, owner):
        if a // 12 not in owner or b // 12 not in owner or b not in self.car.neighbours[a]:
            raise ValueError('local contact lies outside the patch or off the native incidence')

    def swap(self, a, b, owner):
        self.native(a, b, owner)
        self.loads[a], self.loads[b] = self.loads[b], self.loads[a]
        if self.coherent:
            left = self.coherent.pop(a, None)
            right = self.coherent.pop(b, None)
            for port, item in ((a, right), (b, left)):
                if item is not None:
                    probe, index = item
                    probe.rails[index] = port
                    self.coherent[port] = item
        self.gates['pi_contacts'] += 1
        self.gates['intercell_pi_contacts'] += a // 12 != b // 12

    def flag(self, a, b, owner):
        self.native(a, b, owner)
        if a in self.coherent or b in self.coherent:
            raise ValueError('sharp read encountered an unprojected oscillator rail')
        self.gates['echo_read_undo'] += 1
        return int(self.loads[a] + self.loads[b]) == 2

    def fair(self, cell, owner):
        a, b, vacancy = (12 * cell + i for i in (2, 3, 4))
        self.native(a, b, owner)
        if self.loads[a] + self.loads[b] != 2 or self.loads[vacancy] != 0:
            raise ValueError('clock source or vacancy is unprepared')
        # Half contact entangles the two rails. Resolving rail b against its vacancy projects them.
        outcome = self.rng.random() < self.rates[int(self.loads[a])]
        self.loads[a], self.loads[b] = (0, 2) if outcome else (2, 0)
        self.gates['half_contacts'] += 1
        return int(self.flag(b, vacancy, owner))

    def read(self, digit, owner):
        a, _, vacancy = digit
        if self.loads[vacancy]:
            raise ValueError('occupied readout vacancy')
        return int(self.flag(a, vacancy, owner))

    def assign(self, digit, value, owner):
        a, b, _ = digit
        if self.read(digit, owner) != value:
            self.swap(a, b, owner)

    def prepare_oscillator(self, cell, owner):
        probe = SeamOscillator(cell)
        self.native(*probe.rails, owner)
        if any(self.loads[p] or p in self.coherent for p in probe.rails):
            raise ValueError('oscillator needs two unused rails')
        for i, p in enumerate(probe.rails):
            self.coherent[p] = probe, i
        return probe

    def _check_tick_digit(self, digit, owner):
        a, b, vacancy = digit
        self.native(a, b, owner)
        self.native(a, vacancy, owner)
        if (self.loads[vacancy] or self.loads[a] not in (0, 2) or self.loads[a] + self.loads[b] != 2
                or any(p in self.coherent for p in digit)):
            raise ValueError('public tick requires orthogonal double rails and a sharp vacant reader')

    def publish_tick(self, digit, admitted, owner):
        """Commit both outcomes by native read and conditional exchange, overwriting any preceding tick."""
        self._check_tick_digit(digit, owner)
        if admitted not in (0, 1):
            raise ValueError('a clock publication is one recorded bit')
        self.assign(digit, admitted, owner)

    def read_tick(self, digit, owner, *, at_hand=None):
        """Read the committed bit. A remote hand temporarily brings its first rail to a local vacancy."""
        self._check_tick_digit(digit, owner)
        if at_hand is None or at_hand == digit[0] // 12:
            return self.read(digit, owner)
        tip, vacancy = 12 * at_hand + 9, 12 * at_hand + 4
        self.native(tip, vacancy, owner)
        if self.loads[tip] or self.loads[vacancy] or tip in self.coherent or vacancy in self.coherent:
            raise ValueError('hand receiver requires a vacant local pair')
        swap_word = self.route(digit[0], tip, owner=owner)
        for a, b in swap_word:
            self.swap(a, b, owner)
        try:
            return int(self.flag(tip, vacancy, owner))
        finally:
            for a, b in reversed(swap_word):
                self.swap(a, b, owner)

    def get(self, digits, owner):
        return sum(self.read(d, owner) << i for i, d in enumerate(digits))

    def put(self, digits, value, owner):
        if value < 0 or value >= (1 << len(digits)):
            raise OverflowError('finite material register exceeded')
        for i, d in enumerate(digits):
            self.assign(d, (value >> i) & 1, owner)

    def initialize_digits(self, digits):
        for _, b, _ in digits:
            self.loads[b] = 2

    def prepare_front(self, cell):
        for label in (2, RING[0], 8):
            self.loads[12 * cell + label] = 2

    def post(self, cell, digits, owner):
        count = self.get(digits, owner)
        if count == (1 << len(digits)) - 1:
            raise OverflowError('elapsed word exhausted')
        ring = [12 * cell + i for i in RING]
        for i in range(len(ring) - 1, 0, -1):
            self.swap(ring[i], ring[i - 1], owner)
        self.put(digits, count + 1, owner)

    def post_tick(self, public_digit, cell, digits, owner):
        """Drive the five-position hand from a transferred physical event pointer."""
        tick = self.read_tick(public_digit, owner, at_hand=cell)
        if tick:
            self.post(cell, digits, owner)
        return tick

    def reading(self, cell, digits, owner):
        count = self.get(digits, owner)
        ring = [12 * cell + i for i in RING]
        marks = [self.flag(ring[i-1], p, owner) for i, p in enumerate(ring)]
        if marks != [count % 5 in ((i - 1) % 5, i) for i in range(5)]:
            raise AssertionError('phase hand and elapsed register disagree')
        return count

    def route(self, a, b, owner=None):
        cells = frozenset((a // 12, b // 12) if owner is None else owner)
        key = (a, b, cells)
        if key not in self.routes:
            via = {a: None}
            queue = deque([a])
            while b not in via:
                if not queue:
                    raise ValueError('no native route within the acting patch')
                x = queue.popleft()
                for y in self.car.neighbours[x]:
                    if y // 12 in cells and y not in via:
                        queue.append(y)
                        via[y] = x
            path = [b]
            while path[-1] != a:
                path.append(via[path[-1]])
            path.reverse()
            hops = list(zip(path, path[1:]))
            self.routes[key] = hops + hops[-2::-1]
        return self.routes[key]

    def shift(self, cells, right):
        owner = set(cells)
        pairs = list(zip(cells, cells[1:]))
        for src, dst in reversed(pairs) if right else pairs:
            for i in range(12):
                for a, b in self.route(12 * src + i, 12 * dst + i):
                    self.swap(a, b, owner)

    def signed_record(self, digit, value, writer, readers):
        """Prepare a strict native write. Each owner unwrites, reads a rail and restores the same sign."""
        self.assign(digit, value, writer)
        a, b, _ = digit
        self.native(a, b, writer)
        column = law.SRC.index((int(self.loads[a]), int(self.loads[b])))
        written = self.echo[:, column]
        if np.linalg.norm(written[np.arange(law.DIM) % 9 != law.SRC.index((1, 1))]) > 1.e-12:
            raise AssertionError('strict write lost its public agreement')
        # Keep the two open loads in the written state while this seam carries its collar vector.
        self.loads[a] = self.loads[b] = 1
        self.gates['strict_writes'] += 1
        outcomes = []
        for reader in readers:
            self.native(a, b, reader)
            recovered = self.echo @ written
            self.loads[a], self.loads[b] = law.SRC[int(np.argmax(abs(recovered)))]
            outcomes.append(self.read(digit, reader))
            written = self.echo[:, law.SRC.index((int(self.loads[a]), int(self.loads[b])))]
            self.loads[a] = self.loads[b] = 1
            self.gates['strict_rereads'] += 1
        self.native(a, b, writer)
        restored = self.echo @ written
        self.loads[a], self.loads[b] = law.SRC[int(np.argmax(abs(restored)))]
        return outcomes


def _fields(banks, age_bits):
    digits = [d for c in banks for d in _digits(c)]
    count, work = digits[:9], digits[9:]
    fields = {}
    for name, width in (('age', age_bits), ('closing', age_bits), ('address', age_bits), ('position', 5),
                        ('depth', 4), ('attempt', 4), ('accepted', 1), ('suppression', 1)):
        fields[name], work = work[:width], work[width:]
    if work or len(fields['suppression']) != 1:
        raise ValueError('work bank has the wrong size')
    return count, fields


def _material_update(registers, owner, fields, moved, primary, age_bits, draw):
    def get(name):
        return registers.get(fields[name], owner)

    def put(name, value):
        registers.put(fields[name], value, owner)

    age = get('age')
    capacity = _capacity(age_bits)
    put('age', 0 if moved else min(age + 1, capacity))
    if not moved or not primary or age == capacity:
        return bool(primary)
    if age == 0:
        return False
    put('closing', age + 1)
    put('suppression', 0)
    put('position', (2 * draw() - 1) % 32)
    put('depth', 1)
    while get('depth') < STEPS:
        q = get('closing')
        put('attempt', 0)
        while get('attempt') < ADDRESS_TRIES:
            put('attempt', get('attempt') + 1)
            put('address', 0)
            for _ in range((q - 1).bit_length()):
                put('address', 2 * get('address') + draw())
            if get('address') < q:
                break
        if get('address') != q - 1:
            break
        put('position', (get('position') + 2 * draw() - 1) % 32)
        put('depth', get('depth') + 1)
        if get('position') == 0:
            put('suppression', 1)
            break
    return not get('suppression')


def material_replay(car, intervals=(2, 6), *, age_bits=3, seed=20261121, analyzer=0, rest_after=0):
    """Carry the entire clock one adjacent step and back per interval pair. Read two hands at the shared seam."""
    _capacity(age_bits)
    intervals = tuple(index(q) for q in intervals)
    if not intervals or len(intervals) % 2 or any(q < 1 for q in intervals):
        raise ValueError('supply one or more pairs of positive out-and-return intervals')
    rest_after = index(rest_after)
    if rest_after < 0:
        raise ValueError('rest after return must be nonnegative')
    work_bits = 15 + 3 * age_bits
    n_bank = (9 + work_bits) // 3
    size = 1 + n_bank
    # The guide and quiet windows are supplied. Home can read the overlap and has five private cells.
    path = _guide(car, size + 1 + 5)
    rail = path[:size + 1]
    private = path[size + 1:]
    shared = rail[size - 1]
    home = set(private) | set(rail[-2:])
    home_front = private[0]
    home_count = [d for c in private[1:4] for d in _digits(c)]
    receipt = [d for c in private[1:4] for d in _digits(c, EXTRA_DIGITS)] + _digits(private[4])
    reg = Registers(car, seed)
    reg.prepare_front(shared)
    reg.prepare_front(home_front)
    home_digit = _digits(private[-1], EXTRA_DIGITS[1:])[0]
    reg.initialize_digits(home_count + receipt + [home_digit])
    home_oscillator = reg.prepare_oscillator(private[-1], home)
    moving_oscillator = reg.prepare_oscillator(rail[-3], set(rail[:-1]))
    for bank in rail[:n_bank]:
        reg.initialize_digits(_digits(bank))
    # Strict movement record uses the spare (6,7) pair on the trailing bank. Its address moves with the patch.
    reg.initialize_digits(_digits(rail[0], EXTRA_DIGITS[:1]))
    position = 0
    events = []
    for segment, q in enumerate((*intervals, rest_after)):
        for opportunity in range(1, q + 1):
            transfer_requested = segment < len(intervals) and opportunity == q
            cells = rail[position:position + size]
            prior = cells[-1]
            owner = set(cells)
            if not reg.flag(12 * prior + 8, 12 * prior + 9, owner):
                raise AssertionError('departure marker absent before service')
            if transfer_requested:
                reg.shift(rail, right=(position == 0))
            moved = not reg.flag(12 * prior + 8, 12 * prior + 9, owner)
            if moved:
                position = 1 - position
            cells = rail[position:position + size]
            owner = set(cells)
            front = cells[-1]
            if not reg.flag(12 * front + 8, 12 * front + 9, owner):
                raise AssertionError('destination failed its native marker read')
            observed, = reg.signed_record(_digits(cells[0], EXTRA_DIGITS[:1])[0], int(moved), owner, (owner,))
            count, fields = _fields(cells[:-1], age_bits)
            primary = reg.fair(front, owner)
            extra = []

            def draw():
                extra.append(reg.fair(front, owner))
                return extra[-1]

            reg.publish_tick(fields['accepted'][0],
                             _material_update(reg, owner, fields, observed, primary, age_bits, draw), owner)
            tick = reg.post_tick(fields['accepted'][0], front, count, owner)
            oscillator_tick = moving_oscillator.on_tick(reg, fields['accepted'][0], owner)
            if oscillator_tick != tick:
                raise ArithmeticError('the two traveller instruments read different elapsed events')
            tape = iter(extra)
            expected = closing_tick(opportunity, primary, lambda: next(tape), age_bits=age_bits) if moved else primary
            if tick != expected or next(tape, None) is not None or moved != transfer_requested:
                raise AssertionError('independent finite-walk replay or marker disagrees')
            home_primary = reg.fair(home_front, home)
            reg.publish_tick(home_digit, home_primary, home)
            home_tick = reg.post_tick(home_digit, home_front, home_count, home)
            home_oscillator_tick = home_oscillator.on_tick(reg, home_digit, home)
            if home_oscillator_tick != home_tick:
                raise ArithmeticError('the two home instruments read different elapsed events')
            events.append(dict(cell=front, handoff=bool(observed), primary=primary, extra_bits=extra,
                               tick=int(tick), oscillator_tick=oscillator_tick,
                               home_tick=home_tick, home_oscillator_tick=home_oscillator_tick, age=reg.get(fields['age'], owner),
                               traveller=reg.reading(front, count, owner), home_primary=home_primary,
                               home=reg.reading(home_front, home_count, home)))
    if position != 0 or front != shared:
        raise AssertionError('the travellers did not reunite at the original seam')
    seam = (12 * shared + 8, 12 * shared + 9, 12 * shared + 6)
    publication = {}
    for name, digits, writer in (('traveller', count, owner), ('home', home_count, home)):
        word = []
        for i, digit in enumerate(digits):
            value = reg.read(digit, writer)
            one, two = reg.signed_record(seam, value, writer, (owner, home))
            if one != two:
                raise AssertionError('overlap observers read different signs')
            word.append(one)
            if name == 'traveller':
                reg.assign(receipt[i], two, home)
        publication[name] = sum(bit << i for i, bit in enumerate(word))
    oscillator_readings = {}
    for name, probe, writer in (('traveller', moving_oscillator, owner), ('home', home_oscillator, home)):
        terminal = probe.terminal(reg, writer, analyzer)
        # Deliver the projected sharp rail to the common seam, where the local
        # publication controller receives it through a native pointer.
        reg.assign(seam, 1, writer)  # Leave port 9 vacant for its travelling rail.
        published = reg.read_tick((*probe.rails, probe.rails[0] // 12 * 12 + 10), writer, at_hand=shared)
        terminal['owners'] = reg.signed_record(seam, published, writer, (owner, home))
        if terminal['owners'] != [terminal['bit']] * 2:
            raise ArithmeticError('owners disagree on the terminal dimer pointer')
        oscillator_readings[name] = terminal
    reg.assign(seam, 1, owner)  # Restore the prepared departure marker after the signed publications.
    if reg.get(receipt, home) != publication['traveller']:
        raise AssertionError('home lost the copied traveller word')
    return dict(seed=seed, carrier=car.name, intervals=list(intervals), age_bits=age_bits,
                capacity=_capacity(age_bits), work_bits=work_bits, moving_cells=size, transfer_cells=size + 1,
                guide=rail, home_private=private, shared_cell=shared, seam=list(seam[:2]),
                gates=dict(reg.gates), max_native_column_error=reg.error, local_events=events,
                reunion=publication, home_retained_traveller=reg.get(receipt, home),
                elapsed_event_registers={'home': list(home_digit),
                                         'traveller_before_transfer': list(_fields(rail[:n_bank], age_bits)[1]['accepted'][0]),
                                         'traveller_after_transfer': list(_fields(rail[1:1+n_bank], age_bits)[1]['accepted'][0])},
                hand_period=5, counter_bits=9, oscillators=oscillator_readings,
                **({'rest_after': rest_after} if rest_after else {}))


def trip_cli(argv=None, *, default_bits=3):
    import argparse
    import json
    from model.executable_model import federation

    parser = argparse.ArgumentParser(description=material_replay.__doc__)
    parser.add_argument('--level', type=int, choices=(0, 1, 2), default=0)
    parser.add_argument('--age-bits', type=int, choices=(3, 4, 5), default=default_bits)
    parser.add_argument('--seed', type=int, default=20261121)
    parser.add_argument('--rest-after', type=int, default=0, help='continue local repairs at home before final readout')
    parser.add_argument('intervals', nargs='*', type=int, default=[2, 6])
    args = parser.parse_args(argv)
    if args.rest_after < 0:
        parser.error('rest after return must be nonnegative')
    if args.seed < 0 or len(args.intervals) % 2 or any(q <= 0 for q in args.intervals):
        parser.error('choose a nonnegative seed and pairs of positive out-and-return intervals')
    replay = material_replay(federation(args.level), args.intervals, age_bits=args.age_bits, seed=args.seed, rest_after=args.rest_after)
    print(json.dumps(replay, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(trip_cli())
