"""Native source, complete preparation and thermodynamic response checks."""
from fractions import Fraction as F
import numpy as np
import pytest
import scipy.linalg as la
import sympy as s

from . import state_response as r


def array(matrix):
    return np.array(matrix.tolist(),dtype=complex)


def test_unchanged_full_native_action_components_and_current():
    c=r.exact_control()
    assert c['dimension']==208 and c['charge']==6
    assert c['H']==c['HE']+c['HW']+c['HD']
    assert c['H']==c['H'].H
    assert c['H']*c['Gauss']==c['Gauss']*c['H']
    assert s.I*(c['H']*c['E']-c['E']*c['H'])==c['current']
    assert c['current']!=s.zeros(208)


def test_source_coordinate_bound_includes_every_dirac_branch():
    c=r.field_route.rotor_coefficients()
    z=s.Rational(3,5)+s.I*s.Rational(4,5)
    h=c['beta']+c['one_hop']*z+c['one_hop'].H*s.conjugate(z)
    assert (h*h-s.Rational(4,5)*s.eye(4)).applyfunc(s.simplify)==s.zeros(4)
    energies=la.eigvalsh(array(r.exact_control()['H']))
    assert energies.min()>=-4 and energies.max()<=22


def test_reducing_blocks_keep_every_source_basis_state():
    blocks=r.connected_blocks()
    assert sorted(i for block in blocks for i in block)==list(range(208))
    assert dict(sorted(r.Counter(map(len,blocks)).items()))=={1:52,4:20,6:4,8:5,12:1}
    h=r.exact_control()['H']
    for block in blocks:
        outside=tuple(i for i in range(208) if i not in block)
        assert h.extract(block,outside)==s.zeros(len(block),len(outside))


def test_real_arithmetic_chart_restores_native_action_gauge_and_flux():
    phases,hr=r.real_computation_chart(); h=r.exact_control()['H']
    assert all(s.conjugate(p)*p==1 for p in phases)
    assert all(v.is_Rational for v in hr.todok().values())
    for ij,value in h.todok().items():
        assert s.expand(phases[ij[0]]*hr[ij]*s.conjugate(phases[ij[1]])-value)==0
    phase=s.diag(*phases)
    assert phase*r.exact_control()['Gauss']==r.exact_control()['Gauss']*phase
    assert phase*r.exact_control()['E']==r.exact_control()['E']*phase


def test_positive_preparation_is_exactly_normalized_stationary_and_reference_safe():
    prepared=r.prepared_density(); rho=prepared['rho']; c=r.exact_control()
    assert s.trace(rho)==1 and rho==rho.H
    assert c['H']*rho==rho*c['H']
    assert c['Gauss']*rho==rho*c['Gauss']
    assert prepared['preparation_trace_error_upper']<F(1,10**18)
    assert all(q==q.H for _,q in prepared['amplitude_blocks'])
    assert all(all(part.is_Rational for part in v.as_real_imag()) for v in rho.todok().values())
    # Every occupation0..4 has positive mass. No Gauss/number sector filter.
    n=13
    for occupation in range(5):
        weight=sum(rho[i,i] for i in range(208) if (i//n).bit_count()==occupation)
        assert weight>0


def test_canonical_vector_tail_bounds_a_nontrivial_low_order_approximation():
    # Low order makes the analytic bound larger than floating roundoff.
    p=r.prepared_density(order=1)
    h=array(r.exact_control()['H']); e=la.expm(-h/128)
    q=np.zeros_like(h)
    phases=np.array(p['computational_basis_phases'],dtype=complex)
    for ix,block in p['amplitude_blocks']:q[np.ix_(ix,ix)]=array(block)
    q=phases[:,None]*q*phases.conjugate()[None,:]
    target=e.ravel()/la.norm(e); actual=q.ravel()/la.norm(q)
    distance=np.sqrt(max(0.,1-abs(np.vdot(target,actual))**2))
    assert distance<float(p['preparation_trace_error_upper'])
    assert distance>1e-5


def test_complete_cayley_unitarity_and_exact_clock_bound_on_all_blocks():
    for ix in r.connected_blocks():
        re,im=r.cayley_parts(ix)
        assert re.T*re+im.T*im==s.eye(len(ix))
        assert re.T*im-im.T*re==s.zeros(len(ix))
        u=array(r.cayley_block(ix))
        exact=la.expm(-1j*array(r.exact_control()['H'].extract(ix,ix))/64)
        assert la.norm(u-exact,2)<=float(r.cayley_bound())+1e-13
    assert r.cayley_bound()==F(1331,6291456)


def test_all_original_flux_pairs_and_zeros_are_retained():
    law=r.record_law()
    assert law['measurement_times']==[F(0),F(1,64)]
    assert len(law['records'])==169
    assert {(x['first_flux'],x['second_flux']) for x in law['records']}=={
        (f,g) for f in range(-6,7) for g in range(-6,7)}
    assert sum(x['original_probability'] for x in law['records'])==1
    assert all(x['original_probability']>=0 for x in law['records'])
    assert law['zero_probability_labels']>0
    assert any(x['first_flux']!=x['second_flux'] and x['original_probability']>0 for x in law['records'])
    assert law['Gibbs_exact_evolution_output_half_trace_error_upper']<F(1,4000)


def test_flux_records_against_independent_full_matrix_reference():
    p=r.prepared_density(); rho=array(p['rho']); h=array(r.exact_control()['H'])
    eye=np.eye(208); step=(eye-1j*h/512)@la.inv(eye+1j*h/512)
    u=np.linalg.matrix_power(step,4)
    law=r.record_law(); flux=np.tile(np.arange(-6,7),16)
    for row in law['records']:
        first=flux==row['first_flux']; second=flux==row['second_flux']
        reduced=u[np.ix_(second,first)]
        probability=np.trace(reduced@rho[np.ix_(first,first)]@reduced.conj().T).real
        assert abs(probability-float(row['original_probability']))<3e-14


def test_gibbs_kms_heat_and_strain_identities_use_the_same_H():
    c=r.exact_control(); h=array(c['H']); beta=1/64
    en,v=la.eigh(h); p=np.exp(-beta*en); p/=p.sum(); rho=(v*p)@v.conj().T
    a=array(c['E']); b=array(c['current']); t=1/32
    def alpha(z,op):
        return (v*np.exp(1j*z*en))@v.conj().T@op@(v*np.exp(-1j*z*en))@v.conj().T
    assert abs(np.trace(rho@a@alpha(t+1j*beta,b))-np.trace(rho@alpha(t,b)@a))<2e-12
    k=array(c['HE']-c['HD']-2*c['HW'])
    expectation=np.trace(rho@k).real
    covariance=np.trace(rho@h@k).real-(p@en)*expectation
    ds_target=-beta**2*covariance
    def entropy(tau):
        ev=la.eigvalsh(h+tau*k); q=np.exp(-beta*ev);q/=q.sum()
        return -q@np.log(q),-np.log(np.exp(-beta*ev).sum())/beta
    ds=(entropy(1e-4)[0]-entropy(-1e-4)[0])/2e-4
    df=(entropy(1e-4)[1]-entropy(-1e-4)[1])/2e-4
    assert abs(ds-ds_target)<2e-8
    assert abs(df-expectation)<2e-7


def test_mixed_spatial_energy_witness_and_source_scope():
    witness=r.mixed_ground_obstruction(side=2,spacing='1/2')
    assert witness['same_H_energy_decrease']==F(12000,17)
    beta=s.kronecker_product(r.field_route.rotor_coefficients()['beta'],s.eye(15))
    filled=(s.eye(60)-beta)/2
    assert filled*filled==filled and s.trace(filled)==30
    assert s.trace(beta*filled)==-30
    assert s.trace(beta/2)==0
    assert sum(r.source_data.load()['charges'])==0
    assert not witness['physical_vacuum_continuum_impossibility_claimed']


@pytest.mark.parametrize('call',[
    lambda:r.preparation_bound(order=True),
    lambda:r.preparation_bound(beta='0'),
    lambda:r.preparation_bound(order=0,beta='1'),
    lambda:r.cayley_bound(steps=0),
    lambda:r.mixed_ground_obstruction(side=1),
    lambda:r.mixed_ground_obstruction(spacing='0'),
])
def test_invalid_source_and_error_data_are_rejected(call):
    with pytest.raises(ValueError):call()


def test_certificate_separates_gibbs_target_finite_preparation_and_physical_claims():
    cert=r.certificate()
    assert cert['source_thermal_identities']['polynomial_preparation_exactly_Gibbs'] is False
    assert cert['source_thermal_identities']['arbitrary_reference_operations_passive'] is False
    assert cert['finite_state_preparation']['joint_purification_stationary'] is False
    assert cert['physical_vacuum_or_continuum_derived'] is False
    assert cert['physical_temperature_energy_or_mass_calibrations']==[]
