"""Exact controls for complete local instruments and their capacity bounds."""
from fractions import Fraction as F
from itertools import product
from math import comb
import json

import pytest

from . import limit as l, packet


@pytest.mark.parametrize('seams', [630, 2520, 10080])
def test_one_common_clock_and_factor_two(seams):
    bins = (F(0), F(1, 126), F(1, 63))
    ticks = l.aligned_clock(seams, bins)
    assert ticks == (0, seams//126, seams//63)
    assert l.tick_bin(ticks[1], seams, bins) == 0
    assert l.tick_bin(ticks[1]+1, seams, bins) == 1
    assert F(ticks[-1], 2*seams) == bins[-1]/2


@pytest.mark.parametrize('bins', [(F(0), F(1, 64)), (0, F(1)),
                                 (F(0), F(0)), (F(1), F(2))])
def test_unproved_clock_changes_rejected(bins):
    with pytest.raises(ValueError):
        l.aligned_clock(630, bins)


@pytest.mark.parametrize('seams,inside', [(1, 1), (30, 3), (630, 33), (10080, 33)])
def test_marked_selector_exact_tv_and_multiple_outcomes(seams, inside):
    result = l.selector_law_control(seams, inside, (F(0), F(1, seams)))
    rate = F(inside, seams)
    lo, hi = l.exp_minus_interval(rate)
    assert result['one_tick_TV_interval'] == (rate*(1-hi), rate*(1-lo))
    # Positive Poisson multiple-event mass must not be silently normalized away.
    assert 1-(1+rate)*hi > 0
    assert result['multiple_marks_retained'] is True


@pytest.mark.parametrize('population', range(2, 12))
def test_state_bound_for_every_local_sector(population):
    for q in range(population+1):
        for retained in range(population+1):
            for density in (F(0), F(1, 3), F(1)):
                exact = sum(comb(retained, k)*abs(
                    l.restriction.pattern_weight(population, q, retained, k)
                    - density**k*(1-density)**(retained-k))
                    for k in range(retained+1))/2
                assert exact <= l.canonical_state_bound(population, q, retained, density)


def family():
    return l.cone.native.Family((1,)*4, ((0, 1), (1, 2), (2, 3)), (1,)*3)


def trace(block):
    return sum(v for (i, a, j, b), v in block.items() if i == j and a == b)


def test_complete_raw_observation_preserves_reference_and_zero_slots():
    graph = family()
    # A source/reference entangled state; remote port 3 has a definite zero.
    rho = {(1 << i, i, 1 << j, j): F(1, 2) for i in (0, 1) for j in (0, 1)}
    observed = ((0, 1),)
    schedules = {word: F(1, 9) for word in product(graph.seams, repeat=2)}
    result = l.observe_schedules(graph, rho, schedules, (0, 1), observed,
                                 3, (F(0), F(2, 3)))
    assert sum(map(trace, result.values())) == 1
    assert any(not block for block in result.values())
    assert {len(key) for key in result} == {0, 1, 2}
    # Complete records include the native tag/coin/current and raw slot.
    assert all(len(event) == 3 and event[0] == 0 for key in result for event in key)
    # Coherences entangled with the external reference are actually present.
    assert any(a != b and value for block in result.values()
               for (i, a, j, b), value in block.items())


def test_within_bin_chronology_cannot_be_discarded():
    graph = family(); rho = {(1, 0, 1, 0): F(1)}
    first = ((0, 1), (1, 2)); second = tuple(reversed(first))
    options = (graph, rho)
    a = l.observe_schedules(*options, {first: F(1)}, (0, 1, 2), (), 3, (F(0), F(2, 3)))
    b = l.observe_schedules(*options, {second: F(1)}, (0, 1, 2), (), 3, (F(0), F(2, 3)))
    assert a != b
    assert sum(map(trace, a.values())) == sum(map(trace, b.values())) == 1


def test_disjoint_unread_repair_is_exact_identity_on_local_observation():
    graph = family()
    rho = {(1 << i, i, 1 << j, j): F(1, 2) for i in (0, 1) for j in (0, 1)}
    remote = l.observe_schedules(graph, rho, {((2, 3),): F(1)}, (0, 1), (),
                                 3, (F(0), F(1, 3)))
    assert remote == {(): rho}


@pytest.mark.parametrize('radius', [0, 1, 2, 3, 4, 10])
def test_radius_factorials_and_native_attempt_count(radius):
    errors = l.radius_errors(2, 3, radius, 30, (F(0), F(1, 3)))
    assert 0 <= errors['native_cone'] <= errors['Poisson_cone'] <= 1
    if radius >= 10:
        assert errors['native_cone'] == 0


@pytest.mark.parametrize('edge', [(False, 1), (0, True), (0.0, 1), [0, 1]])
def test_observed_source_addresses_reject_scalar_aliases(edge):
    with pytest.raises(ValueError):
        l.observe_schedules(family(), {(1, 0, 1, 0): F(1)},
                            {((0, 1),): F(1)}, (0, 1), (edge,),
                            3, (F(0), F(1, 3)))


def test_actual_source_halos_are_orientation_preserving():
    rows = [l.halo_signature(level) for level in (0, 1, 2)]
    assert all(row['internal_edges'] == rows[0]['internal_edges'] for row in rows)
    assert all(row['attachment_ports'] == (1, 3, 10) for row in rows)
    assert all(row['endpoint_orientation_preserved'] for row in rows)
    assert all(len(row['window']) == 15 and len(row['internal_edges']) == 33 for row in rows)


def test_large_capacity_bound_does_not_erase_geometric_error():
    data = packet.geometry.census(200)
    row = l.local_budget(data['ports'], data['ports']//3, data['seams'],
                        15, 33, 12, 6, 1, (F(0), F(1, 63)), F(1, 3))
    assert row['initial_half_trace_upper'] < F(1, 10**120)
    assert row['complete_window_instrument_halfdiamond_upper'] < F(1, 10**120)
    assert row['native_to_common_finite_window_upper'] > F(1, 100)
    assert data['seams']//63 > l.cone.source.MAX_ATTEMPTS


def test_poisson_truncation_keeps_a_separate_overflow_budget():
    assert l.poisson_tail(F(0), 0) == 0
    assert l.poisson_tail(F(1), 0) == 1
    assert l.poisson_tail(F(1), 3) == F(1, 24)


@pytest.mark.parametrize('value', [True, 1.0, '1', -1])
def test_strict_integer_domain(value):
    with pytest.raises(ValueError):
        l.aligned_clock(value, (F(0), F(1)))


def test_receipt_exact_replay_and_false_premise_guard(tmp_path):
    value = packet.produce()
    assert value['scope']['complete_local_raw_records_and_quantum_disturbance'] is True
    assert value['scope']['physical_capacity_density_clock_or_dimension_selected'] is False
    path = tmp_path/'mutated.json'
    value['scope']['visible_native_coins_retained'] = 1
    path.write_text(json.dumps(value))
    with pytest.raises(ValueError):
        packet.verify(path)
