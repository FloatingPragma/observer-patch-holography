from fractions import Fraction as F
import numpy as np
from scipy.linalg import expm
import pytest
from . import bounds as b,algebra as a,precision


def test_constructive_resource_budgets_and_precision_refinement():
    first=b.resource_plan();second=b.resource_plan(error=F(1,100))
    for row in (first,second):
        assert row['total_trace_error']<=row['requested_trace_error']
        assert row['infinite_ground_prepared'] is False
        assert row['native_pulses_upper']>100000
    assert second['radius']>first['radius']
    assert second['port_fraction_bits']>first['port_fraction_bits']
    assert b.cutoff(100,0,1,F(1),F())==0
    assert b.cutoff(300,0,1,F(1),F(2))<2  # beyond old exponent-64 helper cap


def test_overlap_product_formula_bound_against_noncommuting_chain():
    h0=(sum(a.tensor(x+x+'I') for x in 'XYZ')-np.eye(8))/4
    h1=(sum(a.tensor('I'+x+x) for x in 'XYZ')-np.eye(8))/4
    exact=expm(-.7j*(h0+h1))
    for n in (1,3,10):
        approx=np.linalg.matrix_power(expm(-.7j*h0/n)@expm(-.7j*h1/n),n)
        assert np.linalg.norm(exact-approx,2)<=.7**2/n
    assert np.linalg.norm(exact-expm(-.7j*h0)@expm(-.7j*h1),2)>.01


def test_invalid_resource_domains_and_reference_uniform_factor():
    with pytest.raises(ValueError):b.cutoff(1,1,1,F(1),F(1))
    with pytest.raises(ValueError):b.resource_plan(interventions=1,repairs=1,probes=1)
    with pytest.raises(ValueError):b.resource_plan(error=.1)
    tail=b.cutoff(30,0,1,F(1),F(1,32))
    assert b.history_cutoff(30,0,1,F(1),F(1,32),2)==5*tail
