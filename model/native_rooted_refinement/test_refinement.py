from copy import deepcopy
from fractions import Fraction as F
import json
import struct

import numpy as np
import pytest

from . import refinement as r, routing as q, attachment as a, packet


@pytest.mark.parametrize('radius', (0, 1, 2, 4, 8, 16, 32))
def test_exact_corner_addresses_stabilize_with_proved_collar(radius):
    level = r.minimum_level(radius)
    assert 2**level > radius+2
    assert r.signature(level, radius) == r.signature(level+2, radius)
    assert r.infinite_ball(radius) == r.signature(level, radius)
    topology = r.Topology()
    for c in topology.ball(level, radius):
        for side in range(3):
            other, remote_side = topology.neighbour(level, c, side)
            assert topology.neighbour(level, other, remote_side) == (c, side)


def test_infinite_topology_restrictions_and_original_five_corner_faces():
    small = dict(r.infinite_ball(4)); large = dict(r.infinite_ball(8))
    assert all(large[c] == row for c, row in small.items())
    assert {c[0] for c in large} == set(range(5))
    with pytest.raises(ValueError, match='collar'):
        r.signature(1, 4)


@pytest.mark.parametrize('level', (0, 1, 2, 3))
def test_literal_global_prefix_matches_local_rows_and_coordinates(level):
    result = a.source_prefix(level)
    assert result['midpoint_coordinates_bit_identical'] is True
    assert result['complete_routed_neighbour_rows_identical'] is True


@pytest.mark.parametrize('radius', (1, 2, 4, 8))
def test_finite_routed_tables_match_and_keep_source_orientation(radius):
    first = q.routed_window(r.minimum_level(radius), radius)
    last = q.routed_window(14, radius)
    assert q.same_window(first, last)
    assert all(a < b for a, b in last['internal_seams'])
    assert len(last['source_ports']) == 12*len(r.infinite_ball(radius))
    assert all(type(v) is int for v in last['source_ports'])


@pytest.mark.parametrize('mutation', ('endpoint', 'source_address', 'slot', 'bool_level',
                                    'bool_orientation', 'list_tuple', 'extra'))
def test_complete_window_authentication_rejects_typed_mutations(mutation):
    value = deepcopy(q.routed_window(3, 2))
    if mutation == 'endpoint':
        edge = value['internal_seams'][0]
        value['internal_seams'] = ((edge[1], edge[0]),)+value['internal_seams'][1:]
    elif mutation == 'source_address':
        value['source_ports'] = (999,)+value['source_ports'][1:]
    elif mutation == 'slot':
        value['routing_certificates'][0]['ports'] = (0, 1, 2)
    elif mutation == 'bool_level':
        value['level'] = True
    elif mutation == 'bool_orientation':
        value['orientation_preserved'] = 1
    elif mutation == 'list_tuple':
        value['root_ports'] = list(value['root_ports'])
    else:
        value['extra'] = False
    with pytest.raises(ValueError):
        q.authenticate_window(value)


def test_exact_binary32_rounding_includes_nearest_even_and_sign():
    assert q.contraction_error_bound() < q.SCORE_ERROR
    one = F(1); next_value = F(struct.unpack('<f', struct.pack('<I', 0x3f800001))[0])
    midpoint = (one+next_value)/2
    assert q.round32(midpoint) == one
    assert q.round32(midpoint+F(1, 2**60)) == next_value
    assert q.round32(-midpoint) == -one
    assert q.round32(F()) == 0


def test_assignment_certifies_stable_ties_and_rejects_ambiguous_scores():
    scores = np.full((3, 12), -0.5, dtype=np.float64)
    scores[0, 0] = scores[0, 1] = 0.75
    scores[1, 2] = 0.75; scores[2, 3] = 0.75
    assert q.certified_assignment(scores)['ports'] == (0, 2, 3)
    midpoint = 1+2**-24
    scores[0, 0] = float(np.float32(midpoint)); scores[0, 1] = midpoint
    with pytest.raises(ValueError, match='argmax'):
        q.certified_assignment(scores)
    with pytest.raises(ValueError, match='unique'):
        q.certified_assignment(np.full((3, 12), 0.75, dtype=np.float64))


def test_collision_repair_has_a_unique_certified_optimum():
    scores = np.full((3, 12), -0.5, dtype=np.float64)
    scores[:, 0] = (0.9, 0.8, 0.7)
    scores[1, 1] = 0.79; scores[2, 2] = 0.69
    result = q.certified_assignment(scores)
    assert result['ports'] == (0, 1, 2)
    assert result['collision_repaired'] is True and result['certified_score_gap'] > 0


def test_finite_state_and_record_bound_uses_all_numbers_and_counted_clock():
    window = q.routed_window(14, 4)
    row = a.history_budget(window, F(1, 3), (F(), F(1, 1260), F(1, 630)))
    assert len(window['canonical_ports']) == 348
    assert row['population'] == 64424509440
    assert row['attempts'] == 268435456
    assert len(row['observed_root_internal_seams']) == 30
    assert all(set(e) <= set(row['final_observed_root_ports'])
               for e in row['observed_root_internal_seams'])
    assert row['other_window_record_values_traced'] is True
    assert row['budget']['native_to_common_finite_window_upper'] < F(1, 10**6)
    assert 'native_to_matching_infinite_local_process_upper' not in row['budget']
    assert row['arbitrary_correlated_exterior_reference_asserted'] is False
    assert a.occupation(4, F(2, 7))[1] == 240*4**4*2//7
    with pytest.raises(ValueError, match='align'):
        a.history_budget(window, F(1, 3), (F(), F(1, 11)))


def test_host_stop_and_tiny_arithmetic_controls_are_distinct():
    row = a.host_and_arithmetic_controls()
    assert row['admitted_level200_host_status'] == 'RESOURCE_STOP_UNEVALUATED_CARRIER'
    assert row['attempts_before_stop'] == 0
    assert row['four_byte_level_tag_overflow'] is True
    assert row['corner_zero_spherical_area_within64_steps'] is True
    assert row['astronomic_global_graph_execution_claimed'] is False
    with pytest.raises(ValueError):
        q.LocalSource().cell(15, 0)


@pytest.mark.parametrize('value', (True, -1, 1.0, '1'))
def test_integer_domains_reject_aliases(value):
    with pytest.raises(ValueError):
        r.minimum_level(value)


@pytest.fixture(scope='module')
def receipt():
    return packet.produce()


def test_retained_receipt_and_strict_numeric_types(receipt):
    assert packet.canonical(json.loads(packet.OUTPUT.read_text())) == packet.canonical(receipt)
    changed = deepcopy(receipt)
    changed['scope']['finite_local_source_routing_certificates'] = 1
    assert packet.canonical(changed) != packet.canonical(receipt)


def test_source_and_semantic_receipt_mutations_fail(tmp_path, receipt, monkeypatch):
    for path in ('schema', 'source_bindings'):
        changed = deepcopy(receipt)
        if path == 'schema':
            changed[path] = 'other'
        else:
            key = next(iter(changed[path])); changed[path][key] = '0'*64
        output = tmp_path/'mutated.json'; output.write_text(packet.canonical(changed))
        with pytest.raises(ValueError, match='differs'):
            packet.verify(output)
    original = packet.source_bindings
    calls = [0]
    def changing():
        result = original(); calls[0] += 1
        if calls[0] > 1:
            result[next(iter(result))] = '0'*64
        return result
    monkeypatch.setattr(packet, 'source_bindings', changing)
    with pytest.raises(ArithmeticError, match='changed'):
        packet.produce()
