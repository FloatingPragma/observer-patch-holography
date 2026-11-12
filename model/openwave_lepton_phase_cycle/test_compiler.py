from fractions import Fraction as F
import numpy as np
import pytest
from .compiler import formal_order_certificate,trace_error_bound,rounds_for_error,active_matrices,high_precision_cycle,mass_intervals


def test_exact_noncommuting_order_and_coefficient_mutation():
    report=formal_order_certificate()
    assert report['exact_noncommuting_words']==341
    assert report['wrong_p_two_fifths_nonzero_word_errors']>0


def test_full_fock_twirl_and_negative_duration_control():
    for e in (F(1,20000),F(1,7),F(1,3)):
        v,hs,j,target,bad=active_matrices(e)
        sign=np.diag([(-1)**((mask&3).bit_count()) for mask in range(16)])
        assert np.linalg.norm(sign@v@sign+v)<1e-14
        assert bad>.1
        assert np.linalg.norm(sum(hs)/4@j-j@target)<1e-14


def test_exact_budget_and_denominator_guards():
    n=rounds_for_error()
    assert n==15229
    fractional=rounds_for_error(time=F(1,3))
    assert trace_error_bound(fractional,time=F(1,3))<=F(1,10**12)
    assert trace_error_bound(n)<=F(1,10**12)<trace_error_bound(n-1)
    with pytest.raises(ValueError):trace_error_bound(13)
    with pytest.raises(ValueError):rounds_for_error(F(0))
    assert F(mass_intervals(n)['source_energy_error_exact'])<F(1,20000)


def test_complete_high_precision_rest_readout():
    row=high_precision_cycle(F(1,7),100,precision=65)
    assert float(row['nonvacuum_probability'])<1e-12
    assert abs(sum(map(float,row['complete_final_ancillary_probabilities']))-1)<1e-14
    assert abs(sum(map(float,row['complete_Ramsey_Xplus_Xminus_Yplus_Yminus']))-1)<1e-14
    assert abs(float(row['executed_Ramsey_energy'])-1/7)<1e-10
