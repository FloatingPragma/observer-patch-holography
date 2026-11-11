"""Cutoff and finite-size diagnostics use matched ensembles and expose exact-fit designs."""
import copy
import numpy as np
import pytest
from .finite_size import fit_cutoff,effective_exponent,analyse
from .certificate import load_production


def test_declared_cutoff_forms_recover_known_polynomial_without_volume_mixing():
    points=[dict(size=8,beta_G=b,value=160+3/b+7/b**2,err=.1) for b in (6.,8.,12.)]
    fit=fit_cutoff(points,(1,2))
    assert np.allclose(fit['coefficients'],[160,3,7]) and fit['dof']==0
    points[0]['size']=12
    with pytest.raises(ValueError):fit_cutoff(points)


def test_rank_and_error_mutations_are_rejected():
    points=[dict(size=8,beta_G=b,value=160,err=.1) for b in (6.,8.,12.)]
    bad=copy.deepcopy(points);bad[0]['err']=0
    with pytest.raises(ValueError):fit_cutoff(bad)
    bad=copy.deepcopy(points);bad[0]['beta_G']=8
    with pytest.raises(ValueError):fit_cutoff(bad)
    with pytest.raises(ValueError):fit_cutoff(points[:2],(1,2))


def test_volume_exponent_is_independent_of_susceptibility_units():
    a=effective_exponent((8,2,.1),(16,16,.2))
    b=effective_exponent((8,20,1),(16,160,2))
    assert a['value']==3 and a==b


def test_all_retained_cutoff_alternatives_are_exposed_not_selected_by_agreement():
    out=analyse(load_production())
    for quantity in out['matched_size_cutoff'].values():
        assert len(quantity['variants'])==6
        assert len(quantity['zero_dof_variants'])==4
        assert quantity['sensitivity_envelope']>=max(abs(q['value']-quantity['variants']['linear']['value'])
                                                   for q in quantity['variants'].values())
    assert all(p['value']+3*p['stat']<1 for p in out['crossover_volume_exponents'])
    assert all(p['value']-3*p['stat']>1 for p in out['first_order_volume_exponents'])
