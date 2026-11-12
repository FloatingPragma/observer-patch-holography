from fractions import Fraction as F
import math
import numpy as np
import mpmath as mp
from . import relaxation as r


def test_boltzmann_bounds_extreme_values():
    with mp.workdps(100):
        for x in (F(0),F(1,10),F(10),F(1000)):
            lo,hi=r.negative_exp_endpoint(x);value=mp.exp(-mp.mpf(x.numerator)/x.denominator)
            assert mp.mpf(lo.numerator)/lo.denominator<=value<=mp.mpf(hi.numerator)/hi.denominator


def test_complete_boson_instrument_and_actual_stationarity():
    b=r.boson_bath(F(1));_,p,_,_=r.thermal_state(F(1))
    assert len(b['branches'])==339 and sum(x['exact_zero_operator'] for x in b['branches'])==112
    assert np.max(abs(np.asarray(p@b['matrix']).ravel()-p))<2e-14
    assert np.max(abs(np.asarray(b['matrix'].sum(axis=1)).ravel()-1))<2e-14
    assert b['delta']>F(1,1683)


def test_complete_colour_reference_not_classically_erased():
    control=r.reference_controls(F(1))
    assert control['quantum_bath_reference_entanglement_transferred']
    assert abs(control['fine_bath_colour_measurement_negative_control_unhalved_distance']-1)<1e-13


def test_cooling_produces_full_nonbinomial_group_laws():
    result=r.execute(1)
    first=result['stages'][0]['complete_joint_count_readout']
    assert first['largest_fermion_non_binomial_second_moment']>1e-4
    assert len(first['all435_boson_group_count_laws'])==435
    assert all(len(stage['complete1653_boson_state_probabilities'])==1653 for stage in result['stages'])
    assert all(stage['full_source_state_L1_error_bound_numeric']<1e-10 for stage in result['stages'])
