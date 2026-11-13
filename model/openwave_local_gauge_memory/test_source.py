from fractions import Fraction as F
import pytest
import sympy as s
from .source import graph,local_bound,zero_time_records,complete_noncommuting_control,cosine_interval

def test_third_patch_incidence_is_local_and_exact():
    assert [(len(graph(L)['faces']),len(graph(L)['incident_faces'])) for L in (1,2,3,4)]==[(6,5),(36,12),(108,20),(240,28)]
    assert len(graph(3)['links'])==144

def test_fixed_physical_frequencies_and_nonvacuous_interacting_bounds():
    for representation,frequency in (('fundamental','96'),('spin_one','256')):
        rows=[local_bound(L,F(1,4096),representation) for L in (1,2,3)]
        assert {r['free_frequency_native_units_exact'] for r in rows}=={frequency}
        assert all(F(r['normalized_interacting_real_correlation_interval_exact'][0])>F(97,100) for r in rows)
        assert F(rows[2]['heisenberg_operator_error_upper_exact'])>F(rows[0]['heisenberg_operator_error_upper_exact'])

def test_all_records_reference_and_noncommuting_control():
    for r in ('fundamental','spin_one'):
        records=zero_time_records(r)
        assert len(records)==9 and sum(F(x['probability_exact']) for x in records)==1
    assert complete_noncommuting_control()['reference_entanglement_retained']

def test_exact_cosine_bracket_and_invalid_domain():
    lo,hi=cosine_interval(F(1,16))
    assert lo<F(1) and hi<F(1) and hi-lo>0
    with pytest.raises(ValueError):cosine_interval(2)
    with pytest.raises(ValueError):graph(0)

def test_work_derivative_identity_and_missing_zero_outcome():
    x=s.Symbol('x',real=True);eta=s.Rational(1,2)
    k=((1+eta*x)/2,(1-eta*x)/2,s.sqrt((1-eta**2*x*x)/2))
    assert s.simplify(sum(s.diff(a,x)**2 for a in k)-eta**2/(2*(1-eta**2*x*x)))==0
    assert s.simplify(sum(a*a for a in k))==1
    assert s.simplify(sum(a*a for a in k[:2])-1)!=0
