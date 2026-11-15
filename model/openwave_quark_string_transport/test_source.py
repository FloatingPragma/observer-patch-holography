from fractions import Fraction as F
import numpy as np
from scipy.linalg import expm
from .source import jacobi,evolve,full_car_control,records,carrier_comparison

def test_complete_dynamics_and_no_hopping_control():
    for q in(-1,2):
        H=jacobi(q,4);v,eps,_=evolve(H,F(1));actual=expm(-1j*np.array(H,float))[:,0]
        assert max(abs(complex(float(a),float(b))-actual[i])for i,(a,b)in enumerate(v))<2e-14
        zero,_,_=evolve(jacobi(q,4,F(0)),F(1))
        assert all(a==b==0 for a,b in zero[1:])
        assert float(records(q,4,F(1))['conditional_endpoint_probability_intervals'][1][0])>0

def test_full_carrier_hole_sign_and_gauge_unitarity():
    assert full_car_control(2)['full_active_CAR_dimension']==64
    assert full_car_control(-1)['all_hole_columns']==6

def test_enlarged_carrier_comparison_keeps_tail():
    c=carrier_comparison(2,F(1))
    assert F(c['complete_embedded_vector_difference_interval'][0])>0
    assert F(c['complete_embedded_vector_difference_interval'][1])<F(1,100)
