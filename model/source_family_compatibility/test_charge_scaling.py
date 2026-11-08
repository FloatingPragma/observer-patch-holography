"""Complete source checks for volume-uniform strong-contact response bounds."""
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm

from . import charge_scaling as s
from . import composite


@pytest.mark.parametrize('n',[3,4,8])
@pytest.mark.parametrize('kind',['ballistic','atomic'])
def test_exact_full_source_dimer_and_all_original_probabilities(n,kind):
    r=F(1,2) if kind=='ballistic' else F(1,n*n-4)
    row=s.dimer(n,r)
    assert row['positive_shift_exact_rank']==15
    assert row['ground_degeneracy']==1
    for branch in row['complete_local_charge_records'].values():
        assert sum(branch['original_probabilities'])==F(1,2)
        assert branch['preparation_failure_probability']==F(1,2)
        assert branch['other_energy_outcomes_probability']==0
        assert branch['total_probability']==1
    assert row['exact_residual_Gram_norm_upper']==F(5,2)*row['c']**2


def _actual(n,scale):
    row=s.dimer(n,scale)
    data=composite.source((2,1,1),row['c'],row['U'])
    num,den,_=composite.ground(data,row['ground_energy'])
    rho=num/den;k=data['H']/data['scale']-float(row['ground_energy'])*np.eye(len(rho))
    size=len(rho)
    def embed(ops):
        return np.stack([np.concatenate([(a.T@rho).ravel(),(a@rho).conj().ravel()])for a in ops],axis=1)
    v,w=embed(data['a']),embed(data['b'])
    def time(columns,t):
        blocks=columns.reshape(2,size,size*columns.shape[1])
        return np.concatenate([expm(-1j*t*k)@blocks[0],expm(1j*t*k.conj())@blocks[1]]).reshape(columns.shape)
    return row,data,rho,k,v,w,time


@pytest.mark.parametrize('n',[3,8])
def test_independent_full_diagonalization_and_signed_response(n):
    row,data,rho,k,v,w,time=_actual(n,F(1,2))
    energies,vectors=np.linalg.eigh(k)
    assert abs(energies[0])<1e-12
    assert np.count_nonzero(abs(energies)<1e-12)==1
    np.testing.assert_allclose(vectors[:,0,None]@vectors[:,0,None].T,rho,atol=4e-15,rtol=0)
    for operator in (data['a'][0],data['a'][0].T):
        coefficients=vectors.T@operator@rho
        probabilities=np.einsum('ij,ij->i',coefficients,coefficients)
        branch=row['complete_local_charge_records']['addition']
        for energy,weight in zip(branch['energy_atoms'],branch['original_probabilities'],strict=True):
            assert abs(probabilities[abs(energies-float(energy))<1e-10].sum()-float(weight))<1e-14
        assert abs(probabilities.sum()-.5)<1e-14
    np.testing.assert_allclose(v.T@v,np.eye(4),atol=4e-15,rtol=0)
    np.testing.assert_allclose(v.T@w,0,atol=4e-15,rtol=0)
    assert np.linalg.norm(time(v,.17)-v,2)>0


@pytest.mark.parametrize('n',[3,8])
@pytest.mark.parametrize('t',[F(1,16),F(1,2)])
def test_actual_scaled_response_obeys_atomic_comparison(n,t):
    row,_,_,_,v,w,time=_actual(n,F(1,2))
    exact=time(v,float(t/row['alpha']))
    ideal=v*np.cos(float(t))-1j*w*np.sin(float(t))
    error=np.linalg.norm(exact-ideal,2)
    # The exact dimer has C=-3h and R=0, so gamma=3c/2.
    bound=float(abs(t)*3*row['c']/(2*row['alpha']))
    assert 0<error<bound
    correlation=v.T@exact
    assert np.linalg.norm(correlation-np.cos(float(t))*np.eye(4),2)<=error+1e-14


def test_nonlocal_smearing_does_not_have_exact_separate_half_weights():
    _,data,rho,_,_,_,_=_actual(3,F(1,2))
    a=(data['a'][0]+data['a'][2])/np.sqrt(2)
    removal=np.trace(rho@a.T@a)
    addition=np.trace(rho@a@a.T)
    assert abs(removal+addition-1)<1e-14
    assert abs(removal-.5)>.1
    assert abs(addition-.5)>.1


def test_ballistic_strong_contact_vanishes_in_fixed_physical_window():
    values=[s.response_budget(n,F(n*n-4,2),6,1)['original_total_fixed_window_weight_upper']
            for n in(1024,4096,16384)]
    assert F(0)<values[2]<values[1]<values[0]<1
    # The exact positive coefficient gives the n^-2 rate without a fitted exponent.
    n=16384
    assert values[-1] < F(600000,n*n)


def test_ratio_only_atomic_limit_keeps_finite_original_energy():
    n=32;row=s.dimer(n,F(1,n*n-4))
    assert row['U']==1 and row['c']<F(1,10)
    for branch in row['complete_local_charge_records'].values():
        assert all(F(0)<e<1 for e in branch['energy_atoms'])
        assert sum(branch['original_probabilities'])==F(1,2)
    assert row['fixed_energy_window_1_bound']==1


def test_cell_contact_identity_and_uniform_a_fourth_rate():
    spacing=F(1,64);g=F(3,2);u=g/spacing**3
    # Exact onsite occupation eigenvalues prove the whole diagonal identity.
    for up in(0,1):
        for down in(0,1):
            integrated_quartic=spacing**3*(up/spacing**3)*(down/spacing**3)
            assert g*integrated_quartic==u*up*down
            integrated_number=spacing**3*((up+down)/spacing**3)
            assert -(g/(2*spacing**3))*integrated_number==-(u/2)*(up+down)
    values=[s.cell_contact_budget(F(1,n),1,1,6,1) for n in(32,64,128)]
    assert all(row['U']*row['spacing']**3==1 for row in values)
    assert values[2]['original_total_fixed_window_weight_upper']<values[1]['original_total_fixed_window_weight_upper']/15
    assert values[1]['original_total_fixed_window_weight_upper']<values[0]['original_total_fixed_window_weight_upper']/15


@pytest.mark.parametrize('args',[(0,1,6,1),(1,0,6,1),(1,1,0,1),(1,1,True,1),(1,1,6,-1)])
def test_invalid_bound_domains_rejected(args):
    with pytest.raises(ValueError):s.response_budget(*args)


def test_band_and_scaled_time_bounds_preserve_complements():
    row=s.response_budget(1,1000,6,1,band_width=10000,scaled_time=F(1,100))
    assert 0<row['original_total_outside_two_bands_weight_upper']<F(1,2)
    assert 0<row['original_local_each_charge_inside_band_weight_lower']<F(1,2)
    assert row['interaction_clock_vector_error_squared_upper']<4
    assert not row['hard_many_body_gap_derived'] and not row['arbitrary_input_channel_bound']
    with pytest.raises(ValueError):s.response_budget(1,1000,6,1,band_width=0)
    with pytest.raises(ValueError):s.cell_contact_budget(0,1,1,6,1)
    with pytest.raises(ValueError):s.dimer(2)


def test_exact_source_gram_rejects_wrong_interaction_center(monkeypatch):
    s.dimer.cache_clear()
    original=composite.source
    def changed(*args,**kwargs):
        result=original(*args,**kwargs)
        result['H']=result['H']+np.eye(len(result['H']),dtype=np.int64)
        return result
    monkeypatch.setattr(composite,'source',changed)
    with pytest.raises(ArithmeticError):s.dimer(4)
    s.dimer.cache_clear()


def test_certificate_retains_the_two_distinct_scalings_and_no_particle_promotion():
    row=s.certificate()
    assert row['volume_uniform'] and not row['requires_total_edge_norm_projection']
    assert row['physical_species_assigned']==0 and not row['physical_masses_calculated']
    assert not row['all_smearings_have_exact_separate_half_weights']
    assert not row['hard_many_body_gap_derived'] and not row['arbitrary_input_channel_bound']


@pytest.mark.parametrize('mu',[F(-2),F(-1,2),F(0),F(1,2),F(2)])
def test_charge_clock_shift_preserves_at_least_one_escaping_branch(mu):
    row=s.charge_clock_shift_budget(1,1000,6,1,mu)
    assert row['increased_branch']==('removal' if mu>=0 else 'addition')
    assert row['original_window_upper']==1-abs(mu)
    assert not row['shifted_reference_is_assumed_ground']
    unshifted=s.response_budget(1,1000,6,1)['original_local_each_charge_window_weight_upper']
    assert row['original_local_selected_branch_window_weight_upper']<=unshifted
    if abs(mu)>1:assert row['original_local_selected_branch_window_weight_upper']==0
    # Apply the actual local dimer measure, without relabeling its reference.
    d=s.dimer(4,F(1,2));branch=d['complete_local_charge_records'][row['increased_branch']]
    increased=[energy+abs(mu)for energy in branch['energy_atoms']]
    retained=sum((weight for energy,weight in zip(increased,branch['original_probabilities'],strict=True)
                  if F(0)<=energy<=4),F())
    original=sum((weight for energy,weight in zip(branch['energy_atoms'],branch['original_probabilities'],strict=True)
                  if F(0)<=energy<=4),F())
    assert original>0
    assert retained<=original


@pytest.mark.parametrize('n',[3,4,8])
def test_homogeneous_full_ground_and_original_measure_pushforward(n):
    row=s.homogeneous_dimer(n)
    assert row['positive_clock_factor']==3
    assert row['lambda']==F(n*n-4,2*n)
    assert row['same_complete_ground_projector'] and row['same_original_charge_weights']
    assert row['charge_energy_pushforward_checked'] and row['residual_Gram_scales_quadratically']
    assert not row['graph_or_ground_subsequence_changed']
    assert not row['fixed_ratio_volume_uniform_gap_inferred']
    first=s.dimer(n,F(1,2));second=s.dimer(n,F(3,2))
    a=first['complete_local_charge_records']['addition']
    b=second['complete_local_charge_records']['addition']
    endpoint=(a['energy_atoms'][0]+a['energy_atoms'][1])/2
    probability=lambda branch,cutoff:sum((weight for energy,weight in
       zip(branch['energy_atoms'],branch['original_probabilities'],strict=True)if energy<=cutoff),F())
    assert 0<probability(a,endpoint)<F(1,2)
    assert probability(a,endpoint)==probability(b,3*endpoint)


def test_changing_dimensionless_ratio_does_not_preserve_ground_weights():
    first=s.dimer(3);second=s.dimer(4)
    assert first['U']/first['c']!=second['U']/second['c']
    assert first['complete_local_charge_records']['addition']['original_probabilities']!=second['complete_local_charge_records']['addition']['original_probabilities']


def test_homogeneous_comparison_uses_unbounded_exact_cross_products():
    row=s.homogeneous_dimer(3,F(1,10**20),F(3,2))
    assert row['positive_clock_factor']==3*10**20//2
    assert row['same_complete_ground_projector']
    assert row['charge_energy_pushforward_checked']
