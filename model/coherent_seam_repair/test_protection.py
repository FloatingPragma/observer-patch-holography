"""Protection is relative to the executable boundary and move grammar."""
import json

import pytest

from .protection import BoundaryFlux, certificate, OUTPUT


def test_actual_flux_moves_and_cannot_leave_fixed_boundary():
    law = BoundaryFlux()
    state = 1 << law.edges.index((0, 1))
    edge = law.edges.index((0, 5))
    output, record = law.transport(state, edge)
    assert record['moved'] and output != state
    assert law.region_charge(state) == law.region_charge(output) == 1
    assert law.transport(output, edge)[0] == state
    with pytest.raises(ValueError, match='boundary'):
        law.transport(state, law.edges.index((0, 1)))


def test_empty_adjacent_pair_has_retained_no_operation():
    law = BoundaryFlux()
    state = 1 << law.edges.index((0, 1))
    output, record = law.transport(state, law.edges.index((4, 5)))
    assert output == state
    assert not record['moved']
    assert record['before_face_syndrome'] == record['after_face_syndrome']


@pytest.mark.parametrize('edge', [-1, 30, True, 0.5, '0'])
def test_invalid_edge_addresses_rejected(edge):
    law = BoundaryFlux()
    for operation in (law.flip, law.transport):
        with pytest.raises(ValueError):
            operation(0, edge)


@pytest.mark.parametrize('state', [-1, 1 << 30, 0.0, True])
def test_out_of_domain_link_words_rejected(state):
    with pytest.raises(ValueError, match='thirty bits'):
        BoundaryFlux().syndrome(state)


def test_replay_retains_protection_scope_and_counterexample():
    result = certificate()
    assert json.loads(OUTPUT.read_text()) == result
    assert result['protection']['fixed_boundary_odd_flux_cannot_reach_vacuum']
    assert not result['protection']['arbitrary_link_updates_protect_initial_pair']
    assert result['filled_sphere_first_homology_dimension'] == 0
    assert not any(result['scope'].values())
    assert result['single_defect_positions'] == 5
    assert result['involutory_transport_move_checks'] + result['retained_no_operation_checks'] == 50
