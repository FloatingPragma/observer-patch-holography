"""Independent occupation, exterior-order and complete-history controls."""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import comb

import numpy as np
import pytest

from . import many_body as m
from . import statistics as two


def _ordered_oracle(mask, seam, coin, law):
    """Move a named occupied token, then sort its whole exterior basis."""
    occupied = [p for p in range(12) if mask & (1 << p)]
    a, b = seam
    left, right = int(a in occupied), int(b in occupied)
    target_left = (left+right)//2+coin*((left+right) % 2)
    current = left-target_left
    moved = occupied.copy()
    if current:
        origin, target = (a, b) if current == 1 else (b, a)
        moved[moved.index(origin)] = target
    inversions = sum(x > y for i, x in enumerate(moved) for y in moved[i+1:])
    raw = {0: 0, 1: 1, -1: 2}[current]
    if law == 'E':
        raw += 4*((left+right) % 2)
    return raw, sum(1 << p for p in moved), (-1)**inversions, current


@pytest.fixture(scope='module')
def receipt():
    return m.certificate()


@pytest.mark.parametrize('law', ('C', 'E'))
def test_all_4096_states_every_seam_coin_and_both_phase_lifts(law):
    for seam, coin in product(m.SEAMS, (0, 1)):
        incoming, outgoing = Counter(), Counter()
        for mask in range(4096):
            expected = _ordered_oracle(mask, seam, coin, law)
            assert m.column(mask, law, seam, coin) == expected
            assert m.column(mask, law, seam, coin, 1) == (*expected[:2], 1, expected[3])
            incoming[mask] += 1
            outgoing[expected[1]] += 1
        assert set(incoming.values()) == {1}
        # Each coin separately is a recorded isometry, not a bijection after
        # discarding the record. Two coins together give the unital channel.
        assert sum(outgoing.values()) == 4096


def test_complete_raw_menu_all_number_sectors_and_original_h66(receipt):
    assert receipt['carrier']['sector_dimensions'] == [comb(12, n) for n in range(13)]
    assert sum(receipt['carrier']['sector_dimensions']) == 4096
    for law, count in (('C', 240), ('E', 480)):
        rows = m.branches(law)
        assert len(rows) == count
        assert len({(r.seam, r.coin, r.raw) for r in rows}) == count
        assert receipt['native'][law]['columns'] == 30*2*4096
        assert receipt['native'][law]['positive_columns_checked'] == 30*2*4096
        assert receipt['native'][law]['columns_by_number'] == [60*comb(12, n) for n in range(13)]
        assert all(not row.inputs for row in rows if row.abort)
        for old in two.branches(law):
            row = m.branch(law, old.seam, old.coin, old.raw, number=2)
            expected = sorted((sum(1 << p for p in two.PAIRS[a]),
                               sum(1 << p for p in two.PAIRS[b]), sign)
                              for a, b, sign in zip(old.inputs, old.outputs, old.signs))
            assert list(zip(row.inputs, row.outputs, row.signs)) == expected


def _apply(row, state):
    output = np.zeros_like(state)
    output[list(row.outputs)] = np.asarray(row.signs)[:, None]*state[list(row.inputs)]
    return output


def test_complete_adaptive_histories_arbitrary_reference_and_cross_number_coherence():
    rng = np.random.default_rng(102910)
    initial = rng.normal(size=(4096, 3))+1j*rng.normal(size=(4096, 3))
    initial /= np.linalg.norm(initial)
    reference = np.zeros((3, 3), complex)
    total = 0.
    differing = 0
    for coin, raw in product((0, 1), range(4)):
        second_seam = (0, 11) if raw == 1 else (0, 1)
        for second_coin, second_raw in product((0, 1), range(8)):
            labels = (('C', (0, 5), coin, raw),
                      ('E', second_seam, second_coin, second_raw))
            positive, negative = (m.history_map(labels, character) for character in (1, -1))
            assert positive.inputs == negative.inputs and positive.outputs == negative.outputs
            assert positive.squared_amplitude == negative.squared_amplitude == Fraction(1, 4)
            left, right = _apply(positive, initial)/2, _apply(negative, initial)/2
            assert np.allclose(np.abs(left)**2, np.abs(right)**2, atol=1e-18)
            assert np.allclose(left.conj().T@left, right.conj().T@right, atol=1e-16)
            reference += right.conj().T@right
            total += float(np.linalg.norm(right)**2)
            differing += np.linalg.norm(left-right) > .01
    assert total == pytest.approx(1, abs=2e-15)
    assert np.allclose(reference, initial.conj().T@initial, atol=2e-15)
    assert differing > 0


def test_full_car_algebra_and_sector_raising_are_checked(receipt):
    assert receipt['CAR']['anticommutator_columns'] == 2*12**2*4096
    for mask in range(4096):
        for port in range(12):
            for creation in (False, True):
                value = m.ladder(mask, port, creation=creation)
                if value is None:
                    assert bool(mask & (1 << port)) == creation
                else:
                    occupied = [p for p in range(12) if mask & (1 << p)]
                    count = sum(p < port for p in occupied)
                    assert value == (mask ^ (1 << port), (-1)**count)
                    assert value[0].bit_count() == mask.bit_count()+(1 if creation else -1)
    assert receipt['creation_continuation']['spectator_ladder_columns'] == 30*2*4096*10*2
    control = receipt['creation_continuation']['positive_control_failure']
    assert control['left_sign'] == -control['right_sign']


def test_particle_hole_full_operator_and_nonbipartite_triangle(receipt):
    value = receipt['particle_hole']
    assert value['all_native_column_checks'] == 2*30*2*4096
    assert value['uniform_moving_minus_removable_by_site_phases'] is False
    loops = value['odd_triangle_loop_controls']
    assert [(row['number'], row['phase']) for row in loops] == [(1, 1), (11, -1)]
    for row in loops:
        labels = tuple((law, tuple(seam), coin, raw) for law, seam, coin, raw in row['history'])
        history = m.history_map(labels, number=row['number'])
        i = history.inputs.index(row['initial_mask'])
        assert history.outputs[i] == row['initial_mask']
        assert history.signs[i] == row['phase']
        assert history.squared_amplitude == Fraction(1, 8)
    # Any site-sign gauge has trivial product around the actual triangle.
    for signs in product((-1, 1), repeat=3):
        assert product_sign(signs) == 1
    for law, seam, coin in product(('C', 'E'), m.SEAMS, (0, 1)):
        raw_slots = range(4 if law == 'C' else 8)
        for raw in raw_slots:
            other = m.particle_hole_label(law, seam, coin, raw)
            assert m.particle_hole_label(*other) == (law, seam, coin, raw)
            assert ((raw in m.circuit.record_values(law, coin)) ==
                    (other[3] in m.circuit.record_values(law, 1-coin)))


def product_sign(signs):
    return (signs[0]*signs[1])*(signs[1]*signs[2])*(signs[2]*signs[0])


@pytest.mark.parametrize('kwargs', (
    {'mask': True}, {'mask': -1}, {'mask': 4096}, {'mask': 1.},
    {'law': None}, {'law': 'c'}, {'seam': [0, 5]}, {'seam': (0., 5)},
    {'seam': (0, 2)}, {'coin': True}, {'coin': 2}, {'coin': .0},
    {'character': True}, {'character': 0}, {'character': -1.},
))
def test_column_rejects_malformed_public_arguments(kwargs):
    arguments = {'mask': 1, 'law': 'C', 'seam': (0, 5), 'coin': 0, 'character': -1}
    arguments.update(kwargs)
    with pytest.raises(ValueError):
        m.column(**arguments)


@pytest.mark.parametrize('raw', (None, True, 1., -1, 4))
def test_branch_rejects_malformed_raw_addresses(raw):
    with pytest.raises(ValueError):
        m.branch('C', (0, 5), 0, raw)


@pytest.mark.parametrize('number', (True, -1, 13, 2.))
def test_number_sector_rejects_coercion(number):
    with pytest.raises(ValueError):
        m.states(number)


@pytest.mark.parametrize('creation', (None, 0, 1, 'yes'))
def test_ladder_requires_explicit_boolean(creation):
    with pytest.raises(ValueError):
        m.ladder(0, 0, creation=creation)


@pytest.mark.parametrize('entry', (0, 1, 2, 3))
def test_consumed_native_column_rejects_noninteger_payload(monkeypatch, entry):
    original = m.full_domain.branch_action
    def poisoned(*args, **kwargs):
        row = list(original(*args, **kwargs))
        row[entry] = float(row[entry])
        return tuple(row)
    monkeypatch.setattr(m.full_domain, 'branch_action', poisoned)
    with pytest.raises(ArithmeticError, match='exact native output'):
        m.column(3, 'C', (0, 5), 0)


def test_native_phase_mutation_is_caught_independently(monkeypatch):
    original = m.full_domain.branch_action
    def positive(*args, **kwargs):
        raw, output, _, current = original(*args, **kwargs)
        return raw, output, 1, current
    monkeypatch.setattr(m.full_domain, 'branch_action', positive)
    with pytest.raises(ArithmeticError, match='CAR bilinear'):
        m.complete_certificate()


@pytest.mark.parametrize('law', ('C', 'E'))
@pytest.mark.parametrize('entry', range(4))
def test_positive_only_column_mutation_is_caught_by_certificate(monkeypatch, law, entry):
    original = m.column
    def poisoned(mask, actual_law, seam, coin, character=-1):
        row = list(original(mask, actual_law, seam, coin, character))
        if character == 1 and actual_law == law and mask == 1 and seam == (0, 1) and coin == 0:
            row[entry] = (0, 3, -1, -1)[entry]
        return tuple(row)
    monkeypatch.setattr(m, 'column', poisoned)
    # The negative lift and the existing positive countercontrol on (0,5)
    # are unchanged. The replay must inspect this other positive column.
    with pytest.raises(ArithmeticError, match='positive.*native'):
        m.complete_certificate()


def test_census_and_ladder_mutations_fail_not_only_probability_tests(monkeypatch):
    original = m._ladder
    def unsigned(mask, port, creation):
        value = original(mask, port, creation)
        return None if value is None else (value[0], 1)
    monkeypatch.setattr(m, '_ladder', unsigned)
    with pytest.raises(ArithmeticError, match='anticommutator'):
        m.car_certificate()


def test_empty_and_zero_histories_are_retained():
    identity = m.history_map(())
    assert identity.inputs == identity.outputs == tuple(range(4096))
    assert set(identity.signs) == {1} and identity.squared_amplitude == 1
    abort = m.history_map((('C', (0, 5), 0, 3),))
    assert abort.inputs == abort.outputs == abort.signs == ()
    assert abort.squared_amplitude == Fraction(1, 2)
    with pytest.raises(ValueError):
        m.history_map([('C', (0, 5), 0, 0)])
