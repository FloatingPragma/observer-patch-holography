"""Source restrictions, support, apparatus and independent entropy checks."""
from fractions import Fraction as F
import math
import copy

import numpy as np
import pytest
import scipy.linalg as la
import sympy as s

from . import regional as r


def entropy(matrix):
    values=la.eigvalsh(np.asarray(matrix,dtype=complex))
    assert np.min(values)>0
    return -np.sum(values*np.log(values))


def test_full_spatial_car_chart_and_outgoing_link_partition():
    masks=(0,(1<<240)-1,(1<<3)|(1<<71)|(1<<125)|(1<<239))
    for mask in masks:
        data=r.spatial_chart(mask,tuple(range(12)),4,(1,3))
        restored=sum(((data['inside_CAR_mask']>>j)&1)<<i
                     for j,i in enumerate(data['inside_mode_order']))
        restored+=sum(((data['outside_CAR_mask']>>j)&1)<<i
                      for j,i in enumerate(data['outside_mode_order']))
        assert restored==mask
        assert data['inside_link_words']==(3,4,5,9,10,11)
        assert data['outside_link_words']==(0,1,2,6,7,8)
    swapped=r.spatial_chart((1<<0)|(1<<60),tuple(range(6)),2,(1,))
    assert swapped['CAR_sign']==-1


def test_faithfulness_floor_uses_regional_dimension_and_actual_source_width():
    data=r.spatial_faithfulness(side=2,region_cells=1,beta='1/64')
    assert data['regional_dimension_bits_upper']==87
    exponent=data['common_Boltzmann_floor_factored']['negative_exponent']
    assert exponent>=2*F(1,64)*data['source']['spectral_width_upper']
    assert data['ambient_marginal_minimum_factored']['negative_exponent']==exponent+87
    assert not data['physical_Gauss_projection_preserves_ambient_full_rank']
    assert not data['floor_uniform_in_volume_or_zero_temperature_claimed']


def test_selected_gauss_contract_bound_keeps_actual_state_and_support():
    from model.particle_common_contract import contract
    member=contract.member(0);receipt=r.common_contract_bound(member)
    contract.verify_result(member,receipt)
    data=receipt['result']
    assert F(data['beta'])==F(member['state']['beta'])
    assert data['Taylor_order']==member['state']['Taylor_order']
    assert data['dyadic_scale']==member['state']['dyadic_scale']
    assert data['regional_density_floor_factored']['operator']=='Tr_out(P_G)'
    assert not data['regional_full_ambient_faithfulness_claimed']
    assert not data['warm_native_control_used_as_selected_state']
    changed=copy.deepcopy(member);changed['state']['beta']='1/64'
    with pytest.raises(ValueError):r.common_contract_bound(changed)


def test_native_boundary_projector_and_distinct_entropy_conventions():
    data=r.native_boundary_pair()
    assert data['physical_singlet_rank']==5 and data['full_pair_dimension']==225
    assert sum(data['partial_projector_diagonal'])==5
    rows=data['block_data']['rows']
    # The tracial state on the5-dimensional singlet carrier has weights1/5.
    extended=np.asarray([float(x)/5 for x in data['partial_projector_diagonal']])
    edge=sum(math.log(row['boundary_irrep_dimension'])/5 for row in rows)
    assert abs(-np.sum(extended*np.log(extended))-math.log(5)-edge)<1e-14
    assert not data['spatial_Wilson_state_substituted']


def test_boundary_algebra_keeps_supported_center_and_unavailable_sectors():
    data=r.boundary_blocks((('a',2,3,4),('b',1,0,2)),beta_width='1/8')
    assert data['physical_rank']==6
    a,b=data['rows']
    assert a['partial_trace_physical_projector_eigenvalue']==F(3,4)
    assert a['extended_marginal_floor_factor']==F(1,8)
    assert a['observable_block_floor_factor']==F(1,2)
    assert not b['supported'] and b['observable_block_floor_factor']==0
    assert not data['source_region_assumed_to_be_separately_Gauss_singlet']


@pytest.mark.parametrize('physical',[False,True])
def test_native_state_has_same_H_and_correct_gauss_support(physical):
    data=r.native_state(physical)
    assert data['dimension']==(78 if physical else 208)
    assert s.trace(data['rho'])==1
    assert r.clean(data['rho']*data['H']-data['H']*data['rho'])==s.zeros(data['dimension'])
    assert data['minimum_source_weight']>0
    assert not data['postselection']
    if physical:assert {m.bit_count() for m in data['masks']}=={2}


@pytest.mark.parametrize('physical,keep',[(False,'spinor'),(False,'flux'),(True,'spinor'),(True,'flux')])
def test_rational_modular_bounds_against_independent_matrix_log(physical,keep):
    data=r.native_region(physical,keep);sigma=np.asarray(data['sigma'],dtype=complex)
    values,vecs=la.eigh(sigma);actual=(vecs*(-np.log(values)))@vecs.conj().T
    polynomial=np.asarray(data['log']['matrix'],dtype=complex)+math.log(sigma.shape[0])*np.eye(sigma.shape[0])
    error=la.norm(actual-polynomial,2)
    assert np.min(values)>=float(data['regional_minimum_weight'])-2e-15
    assert error<=float(data['log']['operator_error_upper'])+4e-15
    slope=np.trace(np.asarray(data['tangent'],dtype=complex)@actual).real
    lo,hi=map(float,data['modular_slope_interval'])
    assert lo-2e-16<=slope<=hi+2e-16
    assert data['modular_correction']['modular_equals_bare_local_energy_plus_scalar_excluded']
    assert not data['spatial_region']


@pytest.mark.parametrize('physical,keep',[(False,'spinor'),(False,'flux'),(True,'spinor'),(True,'flux')])
def test_complete_source_apparatus_and_finite_relative_entropy(physical,keep):
    data=r.native_region(physical,keep)
    assert len(data['complete_pointer_records'])==4
    assert sum(x['original_probability'] for x in data['complete_pointer_records'])==1
    assert all(x['original_probability']>0 for x in data['complete_pointer_records'])
    assert data['unitary'].H*data['unitary']==s.eye(data['regional_dimension'])
    sigma=np.asarray(data['sigma'],dtype=complex);tangent=np.asarray(data['tangent'],dtype=complex)
    t=float(data['amplitude']);after=sigma+t*tangent
    vals,vecs=la.eigh(sigma);K=(vecs*(-np.log(vals)))@vecs.conj().T
    change=entropy(after)-entropy(sigma)
    remainder=t*np.trace(tangent@K).real-change
    assert -3e-15<=remainder<=float(data['finite_relative_entropy_upper'])+3e-15
    lo,hi=map(float,data['entropy_change_interval'])
    assert lo-3e-15<=change<=hi+3e-15
    assert data['finite_relative_entropy_upper']<=data['inherited_general_tangent_relative_entropy_upper']
    if keep=='spinor':
        N=s.diag(*(m.bit_count() for m in data['state']['masks']))
        assert data['unitary']*N==N*data['unitary']


def test_partial_trace_retains_an_untouched_quantum_reference():
    psi=s.Matrix([1,0,0,1])/s.sqrt(2);rho=psi*psi.H
    assert r.partial_trace(rho,2,2)==s.eye(2)/2
    # A local matrix remains a multiplicative injection on the full factor.
    a=s.Matrix([[0,1],[1,0]]);b=s.diag(1,-1)
    assert s.kronecker_product(a*b,s.eye(2))==s.kronecker_product(a,s.eye(2))*s.kronecker_product(b,s.eye(2))


def test_source_clock_is_distinct_from_polynomial_modular_time():
    source=r.native_state();rho=source['rho'];H=source['H']
    assert r.clean(rho*H-H*rho)==s.zeros(208)
    # Three actual vacuum/flux eigenvalues, with exact rational amplitudes.
    energies=(F(0),F(1,2),F(2))
    for flux,e in zip((0,1,2),energies,strict=True):
        index=flux+6
        assert H[index,index]==s.Rational(str(e))
        assert sum(bool(H[index,j]) for j in range(208))==(0 if e==0 else 1)
    Q=[sum((F(22)-e)**j/F(128)**j/math.factorial(j) for j in range(13)) for e in energies]
    # An affine K(E) would require4K(1/2)=3K(0)+K(2).
    assert Q[1]**4 != Q[0]**3*Q[2]


@pytest.mark.parametrize('bad',[(),(('a',1,0,2),),(('a',1,1,0),),(('a',1,1,2),('a',1,1,2))])
def test_invalid_boundary_data_is_rejected(bad):
    with pytest.raises(ValueError):r.boundary_blocks(bad)


def test_certificate_keeps_geometric_obligations_explicit():
    data=r.certificate()
    assert data['scope']['finite_regional_entropy_first_law_on_actual_source']
    assert not data['scope']['fixed_proper_volume_area_variation_supplied']
    assert not data['scope']['Einstein_dynamics_derived']
    assert len(data['native_internal_controls'])==4
