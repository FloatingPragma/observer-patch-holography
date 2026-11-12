from fractions import Fraction as F
from .certificate import calculate

def test_exact_current_and_certified_held_out_failure():
    r=calculate();assert r['exact_local_current_algebra']['gA']=='5/3'
    lo,hi=map(F,r['finite_hopping_certificate']['axial_ratio_interval_rational'])
    assert lo>F(8,5)and lo<F(5,3)<hi
    assert lo<r['floating_evaluation']['gA']<hi
    assert r['floating_evaluation']['current_leakage_operator_norm'][3]>.9
    assert not r['observed_axial_ratio_reproduced']
