"""Selected-state persistence controls in the local-coarse modular instrument."""
from fractions import Fraction
import pytest
from oph_openwave.compact_native import run


@pytest.fixture(scope='module')
def result():
    return run()


def test_canonical_closed_cycle_and_complete_outcomes(result):
    for control in result['canonical_flat_controls']:
        assert not any(control['integer_cycle_divergence'])
        assert len(control['spatial_support'])==3
        assert control['energy_exact']=='-23/100'
        assert control['every_reduced_kraus_selected_state_max_error']<1e-12
        outcomes=control['complete_native_outcomes_conditional_on_commit']
        assert len(outcomes)==10*control['seams']
        assert sum(Fraction(row['probability_conditional_on_commit']) for row in outcomes)==1
        assert sum(row['record_value']!=0 for row in outcomes)==6*control['seams']


def test_same_actual_instrument_and_record_preparation_controls(result):
    triangle=result['explicit_actual_modular_triangle']
    assert triangle['intertwiner_max_error']==0
    assert triangle['lifted_eigenstate_residual']==0
    assert triangle['local_coarse_step_trace_distance']<1e-12
    assert triangle['definite_per_position_record_lift_port_leak_t037']>.08
    assert triangle['full_monitor_then_coherent_port_leak_t037']>.23
    for control in result['canonical_flat_controls']:
        assert control['broken_cycle_port_leakage_t037']>.1


def test_enormous_record_support_and_scope_not_hidden(result):
    for control in result['canonical_flat_controls']:
        fiber=control['record_fiber']
        assert int(fiber['prepared_nonzero_joint_basis_count_exact'])==3*6**fiber['cycle_rank']
        assert int(fiber['fiber_cardinality_exact'])>10**30
        assert fiber['enumerated'] is False
    assert all(result['checks'].values())
    assert all(result['mutations'].values())
    assert result['full_criterion_met'] is False
    assert result['official_rows_changed']==[]
