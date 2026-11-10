"""Full source encoding and complete finite joined-clock controls."""
import pytest
import sympy as s

from .certificate import experiment
from .model import (State,seed,route,native_generator,joined_taylor,
                    weak_taylor,complete_binary,changed_face,taylor_error)


def test_actual_common_clock_receipt_and_complete_error_bounds():
    result=experiment()
    assert result['all_columns_checked_by_full_source_reencoding']
    assert result['complete_local_routing_columns_with_native_spectators']==121
    assert result['coincident_object_states_retained']==100
    no,yes=result['exact_unitary_probability_enclosures']
    assert s.sympify(yes['lower'])>0
    assert s.sympify(no['lower'])+s.sympify(yes['upper'])==1
    assert s.sympify(no['upper'])+s.sympify(yes['lower'])==1
    assert s.sympify(result['Taylor_output_norm_squared'])>1
    assert result['Taylor_output_renormalized'] is False


def test_empty_native_exact_route_and_native_number_control():
    state=seed(with_native=False);routed=route(state,3,s.pi/12)
    assert native_generator(state).norm_squared()==0
    assert native_generator(routed).norm_squared()==0
    assert weak_taylor(routed,s.Rational(1,7),3)==routed
    assert dict(complete_binary(routed,changed_face))[True].norm_squared()==0
    assert taylor_error(state,s.Rational(1,7),3)==0


def test_error_bound_input_contract_and_zero_coupling():
    state=seed();out,error=joined_taylor(state,3,s.pi/12,0)
    assert out==route(state,3,s.pi/12) and error==0
    doubled=State(dict(state.amplitudes),4)
    with pytest.raises(ValueError,match='normalized'):
        joined_taylor(doubled,3,s.pi/12,1)
    with pytest.raises(ValueError,match='order'):
        weak_taylor(state,1,-1)
