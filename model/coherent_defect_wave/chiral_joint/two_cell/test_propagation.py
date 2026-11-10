import sympy as s
from .propagation import transport_control,pair_preparation,arrival


def test_complete_open_source_has_positive_controlled_transport():
    result=transport_control(False)
    lo,hi=map(s.Rational,result['complete_arrival_yes_probability_interval'])
    assert 0<lo<hi<1
    assert result['exact_first_derivative_arrival_norm_squared']=='3/4'
    assert result['no_hopping_control_arrival_probability_at_all_times']=='0'
    assert not result['prepared_state_is_physical_ground']


def test_periodic_quantum_source_transports_in_both_directions():
    result=transport_control(True)
    lo,hi=map(s.Rational,result['complete_arrival_yes_probability_interval'])
    assert 0<lo<hi<1
    assert result['exact_first_derivative_arrival_norm_squared']=='3/16'
    assert result['opposite_direction_complete_transport_laws_equal_by_exact_translation']
    assert not result['periodic_reference_is_open_source_ground']


def test_preparation_is_local_and_the_receiver_starts_even():
    for periodic in (False,True):
        assert not arrival(pair_preparation(1,periodic),0)
        assert not arrival(pair_preparation(0,periodic),1)
