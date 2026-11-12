from fractions import Fraction as F
import numpy as np
import pytest
from . import bath as b


def test_minimal_reversible_support_and_exact_all_record_laws():
    d=b.source()
    assert d['native_directed_edges']==4104 and d['reverse_edges_added']==3924
    assert len(d['components'])==26
    for q in (F(1),F(1,2)):
        rows,records,_=b.instrument(q)
        total=[F(0)]*1296
        for x,proposal,target,label,p in records:total[x]+=p
        assert total==[F(1)]*1296 and len(records)==17352
        assert all(row[x]>=F(1,2) for x,row in enumerate(rows))


def test_actual_equilibrium_is_not_native_stationarity():
    pi,_=b.equilibrium(F(1,2));d=b.source();kernel=b.native.kernel();out=[F(0)]*1296
    for x,p in enumerate(pi):
        for state,q in kernel[d['states'][x]].items():out[d['index'][state]]+=p*q
    assert out!=pi
    # A pure repeated-load nonzero-flux singleton has independently free weight.
    assert len([x for x in d['components'] if len(x)==1])==12


def test_complete_cooling_preserves_mass_symmetry_and_entropy_law():
    r=b.cooling(4)
    assert r['total_entropy_production']>0 and not r['state_probabilities_renormalized']
    for row in r['stages']:
        p=np.array(row['complete_state_probabilities'])
        assert np.min(p)>=0 and abs(p.sum()-1)<2e-13
        assert row['state_probability_L1_roundoff_bound_numeric']<1e-10
        assert abs(row['entropy_production']-(row['entropy_change']-row['beta_in_declared_cost_units']*row['cost_heat']))<1e-11


def test_invalid_thermal_parameters_rejected():
    for q in (0,-1,F(3,2)):
        with pytest.raises(ValueError):b.instrument(q)
