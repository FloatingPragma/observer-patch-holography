import pytest
import sympy as s
from . import measure as m


def test_actual_casimir_coherence_and_diagonal_nonclosure():
    c=m.closure_certificate()
    assert c['scalar_sector_adjoint_closure']is False
    assert c['coherent_witness']['initial_half_probabilities']==['1/2','1/2']
    assert len(set(c['coherent_witness']['next_half_probabilities']))==2
    assert c['within_j_zero_diagonal_exit_probabilities'][1]=='0'
    assert sum(row['weak_C2']!=0 for row in m.source())==2


def test_stationary_family_keeps_cross_family_quantum_coherence():
    result=m.stationary_certificate();assert result['complex_fixed_operator_dimension']==7
    indices=[i for i,row in enumerate(m.source())if row['d']==3]
    rho=s.zeros(10)
    for i in indices:
        for j in indices:
            for z in (0,1):rho[2*i+z,2*j+z]=s.Rational(1,4)
    assert m.channel(rho)==rho and s.trace(rho)==1
    assert result['tracial_ten_dimensional_preparation']['centre_balance']=='3/5'
    assert result['native_default_d6_preparation']['limiting_centre_balance']=='0'


def test_complete_records_and_parent_gram_replay():
    for row in m.histories(3):
        assert len(row['complete_records'])==8
        assert sum(map(s.Rational,row['complete_records'].values()))==1
    for row in m.source():
        p=row['p'];d=row['d']
        assert p==s.Rational(1024*d,(d+256)**2+256)
        assert 0<p<1


def test_domain_and_complete_effect_controls():
    with pytest.raises(ValueError):m.histories(0)
    with pytest.raises(ValueError):m.pulse(2,s.Float(.125))
    assert m.adjoint(s.eye(10))==s.eye(10)
    # Removing the return outcome loses probability on a legitimate input.
    _,u=m.pulse(2);bad=u.H*s.diag(0,1)*u
    assert not m.zero(bad-s.eye(2))
