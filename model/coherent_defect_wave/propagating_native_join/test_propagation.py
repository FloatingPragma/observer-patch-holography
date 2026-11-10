"""Source intertwiners, retained collisions and complete probability boundaries."""
import sympy as s
import pytest

from .certificate import matching_experiment,spinful_experiment,exact_scaling_checks
from .model import matching_edges


def test_complete_matching_and_same_link_response():
    result=matching_experiment()
    assert result['literal_Gauss_intertwiner_checked']
    records=result['exact_joint_records']
    lower=s.sympify(records[-1]['lower']);upper=s.sympify(records[-1]['upper'])
    assert 0<lower<upper<s.Rational(1,2)
    assert s.sympify(records[2]['lower'])+upper==s.Rational(1,2)
    assert s.sympify(records[2]['upper'])+lower==s.Rational(1,2)
    for phase in range(6):assert len(matching_edges(4,phase))==64
    with pytest.raises(ValueError):matching_edges(3,0)


def test_spinful_common_clock_and_full_native_occupation():
    result=spinful_experiment()
    assert result['native_one_body_modes_per_lane']==28800
    assert result['exact_two_native_exterior_inverse_checked']
    assert result['initial_native_receiver_probability']=='0'
    assert result['propagated_native_receiver_probability']=='1'
    assert result['native_propagation_off_changed_face_probability']=='0'
    records=result['exact_V_clock_joint_records'];yes=records[-1];no=records[-2]
    assert s.sympify(yes['lower'])>0
    assert s.sympify(no['lower'])+s.sympify(yes['upper'])==1
    assert s.sympify(no['upper'])+s.sympify(yes['lower'])==1
    assert result['Taylor_output_renormalized'] is False


def test_dispersion_and_vertex_exclusion_controls():
    result=exact_scaling_checks()
    assert result['standard_spinor_density_vertex_boost_residual_squared']==['16']*3
    assert result['scalar_beta_tensor_beta_residual']=='0'
    assert not result['interacting_state_replaced_by_free_probability']
    assert not result['dense_sea_or_strong_scaling_excluded_by_this_bound']
