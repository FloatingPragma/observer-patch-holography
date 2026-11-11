"""Independent algebra, normalization, mutation and numerical-flow controls."""
from pathlib import Path
import ast

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from . import running,threshold


def test_literal_complex_car_and_full_source_columns():
    assert threshold.exact_local_car()['identity_exact']
    assert threshold.exact_source_commutators()['complete_exact_operator_identities']==72


def test_threshold_positive_rationalized_and_complete_records():
    y=np.diag([1e-5,.05,.9])
    report,a=threshold.analytic_records(y)
    low=np.asarray([r['low_frequency'] for r in report['rows']])
    assert np.all(low>0)
    assert np.all(low[1:]/low[:-1] >= (np.diag(y)[1:]/np.diag(y)[:-1])**2)
    assert report['rows'][0]['upper_weight_per_spin']>0
    assert report['complete_success_plus_blocked_max']<1e-14
    assert np.linalg.eigvalsh(a['blocked_effect']).min()>=-1e-14


def test_inverse_uses_failures_and_resolves_reference_offset():
    y=np.diag([.01,.1,.8])
    _,a=threshold.analytic_records(y,coupling=.5,scalar=2.,reference=.125)
    result=threshold.reconstruct_coupling(a['frequencies'],a['blocked_effect'])
    np.testing.assert_allclose(result['coupling_times_singular_values'],.5*np.diag(y),atol=1e-12)
    assert abs(result['scalar_energy']-2)<1e-12
    assert abs(result['reference_energy']-.125)<1e-12
    wrong=a['frequencies'].copy();wrong[1]+=.1
    with pytest.raises(ValueError):
        threshold.reconstruct_coupling(wrong,a['blocked_effect'])
    changed=threshold.reconstruct_coupling(a['frequencies'],2*a['blocked_effect'])
    assert abs(changed['scalar_energy']-2)>.01


def test_inverse_rejects_invalid_complete_effects_and_nonfinite_parameters():
    _,a=threshold.analytic_records(np.diag([.01,.1,.8]))
    blocked=a['blocked_effect']
    nonhermitian=blocked.copy();nonhermitian[0,1]=.1j
    indefinite=blocked.copy();indefinite[0,0]-=.2;indefinite[1,1]+=.2
    nan=blocked.copy();nan[0,0]=np.nan
    for bad in (nonhermitian,indefinite,nan):
        with pytest.raises(ValueError):threshold.reconstruct_coupling(a['frequencies'],bad)
    bad_lines=a['frequencies'].copy();bad_lines[0]=np.nan
    with pytest.raises(ValueError):threshold.reconstruct_coupling(bad_lines,blocked)
    for option in ({'coupling':np.nan},{'coupling':1j},{'scalar':np.inf},{'reference':np.nan}):
        with pytest.raises(ValueError):threshold.analytic_records(np.eye(3),**option)


def test_flavour_rank_and_zero_controls():
    control=running.algebra_controls()
    assert control['flavour_covariance_max']<1e-14
    assert control['zero_Yukawa_beta_exact']
    # A rank-deficient input is checked through the integrated matrix law.
    y=np.zeros((3,3,3),complex);y[:,0,0]=[.5,.3,.2]
    solution=solve_ivp(running.beta,(0,-10),running.pack([.5]*3,y),rtol=1e-10,atol=1e-12)
    _,out=running.unpack(solution.y[:,-1]);out[:,0,0]=0
    assert np.count_nonzero(out)==0


def test_exact_lepton_ratio_derivative_independent_of_other_channels():
    y=np.zeros((3,3,3),complex)
    y[0]=np.diag([.3,.4,.5]);y[1]=np.diag([.7,.2,.1]);y[2]=np.diag([.01,.1,.8])
    _,b=running.unpack(running.beta(0,running.pack([.4,.6,.8],y)))
    for low,high in ((0,1),(1,2)):
        direct=b[2,high,high]/y[2,high,high]-b[2,low,low]/y[2,low,low]
        exact=1.5*(y[2,high,high]**2-y[2,low,low]**2)/running.LOOP
        assert abs(direct-exact)<1e-16


def test_complete_complex_source_record_match():
    u=np.array([[1,1j,0],[1j,1,0],[0,0,np.sqrt(2)]])/np.sqrt(2)
    y=u@np.diag([.125,.25,.5])@u.conj().T
    result,_=threshold.complete_control(y,np.eye(3))
    assert max(result['analytic_fixed_probe_correlation_checks'].values())<2e-10


def test_empirical_data_is_confined_to_separate_comparator():
    root=Path(__file__).parent
    for name in ('certificate.py','threshold.py','running.py'):
        tree=ast.parse((root/name).read_text())
        values=[x.value for x in ast.walk(tree) if isinstance(x,ast.Constant) and isinstance(x.value,str)]
        assert not any('calibrated_matter_parameters.json' in value or 'pdg_202' in value for value in values)
