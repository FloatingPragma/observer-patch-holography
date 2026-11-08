"""Exact complete preparation/readout composition and resource boundaries."""
from fractions import Fraction as F
from math import comb

import pytest

from . import clock as c, composition as j


def test_complete_endpoint_unitality_and_number_conservation():
    result = j.unitality_certificate()
    assert result['number_preserving_column_checks'] == 100
    assert result['sum_K_K_adjoint_exact_identity'] is True
    assert result['external_reference_independence_required'] is False


def test_unread_local_preparation_preserves_full_canonical_marginal():
    # Every global q sector of a four-port graph; no fixed local number.
    graph = j.native_graph.Family((1,)*4, ((0, 1), (1, 2), (2, 3)), (1,)*3)
    for q in range(5):
        density = {(x, 0, x, 0): F(1, comb(4, q)) for x in range(16) if x.bit_count() == q}
        for edge in graph.seams:
            rows = j.prep.readout.cone.fixed_instrument(graph, density, edge, law='E')
            output = {}
            for _, block in rows:
                for key, value in block.items(): output[key] = output.get(key, F())+value
            assert output == density


def test_constant_detector_weight_preserved_on_arbitrary_input():
    graph = j.native_graph.Family((1,)*3, ((0, 1), (1, 2)), (1,)*2)
    rho = {(0, 0, 0, 0): F(1, 4), (7, 1, 7, 1): F(1, 4),
           (0, 0, 7, 1): F(1, 4), (7, 1, 0, 0): F(1, 4),
           (1, 2, 1, 2): F(1, 2)}
    for edge in graph.seams:
        output = {}
        for _, block in j.prep.readout.cone.fixed_instrument(graph, rho, edge, law='E'):
            for key, value in block.items(): output[key] = output.get(key, F())+value
        eta = sum(value for (x, a, y, b), value in output.items() if x == y and a == b and x in (0, 7))
        assert eta == F(1, 2)
        assert output[0, 0, 7, 1] == F(1, 4)


def test_actual_long_canonical_programme_has_no_extra_preparation_tick():
    result = j.canonical_protocol()
    plan = result['plan']; M = 2580480
    assert plan['total_attempts'] == 165154571
    assert plan['fixed_attempts'] == 3851
    assert plan['uniform_attempts'] == 165150720
    assert plan['program']['declared_resources']['AST_nodes'] == 3852
    assert plan['partition'] == (F(), F(3851, M), F(3851, M)+64)
    assert result['detector_ports'] == 348
    assert result['detector_incident_seams'] == 922
    assert result['prep_mixing_record_values_marginalized'] == 3840
    assert result['external_reference_may_be_correlated'] is True
    bound = result['error']['actual_protocol_to_reference_preparation_oracle_upper']
    assert bound <= F(348, 2**64)+result['initial_and_postprep_ambiguous_weight_upper']
    assert bound < F(1, 10**16)
    assert result['short_time_Poisson_budget_used'] is False


@pytest.mark.parametrize('level', [0, 2, 6, 14, 200])
def test_sparse_two_source_protocol_is_capacity_independent_in_attempts(level):
    result = j.sparse_protocol(level)
    plan = result['plan']
    assert plan['total_attempts'] == plan['fixed_attempts'] == 3862
    assert plan['program']['declared_resources']['AST_nodes'] == 3862
    assert result['feasible_local_numbers'] == (0, 1, 2)
    assert result['readout_oracle_error'] == 0
    assert result['actual_protocol_to_reference_preparation_oracle_upper'] == F(1, 2**64)
    assert plan['partition'][-1] == F(3862, 630*4**level)
    assert result['arbitrary_source_reference_input_in_declared_global_sector'] is True


@pytest.mark.parametrize('total', [0, 1, 2, 11, 229, 238, 239, 240])
def test_public_global_sector_resolves_constant_corner(total):
    result = j.sparse_protocol(0, total, blocks=0)
    assert result['constant_feasibility']['ambiguous_constant'] is False
    assert result['readout_oracle_error'] == 0


@pytest.mark.parametrize('total', [12, 80, 120, 228])
def test_ambiguous_global_sector_cannot_claim_exact_root_readout(total):
    with pytest.raises(ValueError, match='ambiguity'):
        j.sparse_protocol(0, total, blocks=0)


def test_two_tree_node_limit_is_not_preparation_alone_limit():
    assert j.sparse_protocol(0, blocks=67)['plan']['total_attempts'] == 4042
    with pytest.raises(ValueError, match='bounded literal command'):
        j.sparse_protocol(0, blocks=68)


def test_high_capacity_sparse_protocol_keeps_actual_host_stop():
    plan = j.sparse_protocol(200)['plan']
    marker = {'undeclared_density_marker': True}
    stopped = c.native.execute(plan['program'], marker)
    assert stopped.status == 'RESOURCE_STOP_UNEVALUATED_CARRIER'
    assert stopped.attempted_steps == 0
    assert stopped.branches[0].density == marker


def test_arbitrary_input_composition_does_not_infer_a_small_erasure():
    value = j.combined_error(F(1, 2**64), F(1, 2**64), F(1))
    assert value['actual_preparation_to_its_own_oracle_upper'] == 1
    assert value['actual_protocol_to_reference_preparation_oracle_upper'] == 1
    value = j.combined_error(F(1, 2**64), F(0), F(0))
    assert value['actual_protocol_to_reference_preparation_oracle_upper'] == F(1, 2**64)


@pytest.mark.parametrize('bad', [True, 0, 0.5, '0'])
def test_composition_rejects_untyped_error_premise(bad):
    with pytest.raises(ValueError): j.combined_error(bad, F(), F())
